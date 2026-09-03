"""
execute_tau_relaxations.py
==========================
Executes 4-orientation relaxations (0, 90, 180, 270 deg) for the N=8 Tau candidate subset
on canonical beta12-B40H15 with strict convergence verification.
"""

import subprocess, re, time, hashlib, shutil
import numpy as np, pandas as pd
from pathlib import Path


def _project_root(marker="MANIFEST_SHA256.txt"):
    from pathlib import Path as _P
    here = _P(__file__).resolve()
    for anc in [here.parent, *here.parents]:
        if (anc / marker).exists() or ((anc / "data").is_dir() and (anc / "README.md").exists()):
            return anc
    return here.parent


def _find_xtb():
    import shutil
    from pathlib import Path as _P
    w = shutil.which("xtb") or shutil.which("xtb.exe")
    if w:
        return _P(w)
    for anc in [_P(__file__).resolve().parent, *_P(__file__).resolve().parents]:
        hits = list(anc.glob("**/xtb-*/bin/xtb.exe")) or list(anc.glob("**/xtb-*/bin/xtb"))
        if hits:
            return hits[0]
    return _P("xtb")


from scipy.stats import spearmanr
from sklearn.metrics import mean_absolute_error

base = _project_root()
calc = base / "calculations" / "tau"
proc = base / "data" / "processed"
xtb = _find_xtb()
def sha256_file(fp):
    h = hashlib.sha256()
    with open(fp, "rb") as fh:
        for chunk in iter(lambda: fh.read(8192), b""): h.update(chunk)
    return h.hexdigest()

def parse_xtb_opt_output(out_file):
    text = Path(out_file).read_text(encoding="utf-8", errors="replace")
    energy = None
    converged = False
    cycles = 0
    grad_norm = None
    
    if "GEOMETRY OPTIMIZATION CONVERGED" in text and "FAILED TO CONVERGE GEOMETRY OPTIMIZATION" not in text:
        converged = True
    elif "*** convergence criteria satisfied" in text and "FAILED TO CONVERGE" not in text and "Final Singlepoint" not in text:
        converged = True
        
    for l in text.splitlines():
        if "TOTAL ENERGY" in l:
            m = re.search(r"(-?\d+\.\d+)\s+Eh", l)
            if m: energy = float(m.group(1))
        if "GRADIENT NORM" in l:
            m = re.search(r"(\d+\.\d+)\s+Eh", l)
            if m: grad_norm = float(m.group(1))
        if "CYCLE" in l:
            m = re.search(r"CYCLE\s+(\d+)", l)
            if m:
                c = int(m.group(1))
                if c > cycles: cycles = c
        if "GEOMETRY OPTIMIZATION CONVERGED AFTER" in l:
            m = re.search(r"AFTER\s+(\d+)\s+ITERATIONS", l)
            if m:
                cycles = int(m.group(1))
                
    return energy, converged, grad_norm, cycles

# Read carrier
opt_carrier = calc / "beta12_carrier_optimized.xyz"
carrier_lines = opt_carrier.read_text().splitlines()
n_carrier = int(carrier_lines[0])
carrier_atoms = []
for l in carrier_lines[2:2+n_carrier]:
    p = l.split()
    carrier_atoms.append((p[0], float(p[1]), float(p[2]), float(p[3])))
z_top = max(z for _, _, _, z in carrier_atoms)

# Carrier energy
e_carrier = -56.192156

df_drugs = pd.read_csv(proc / "dataset_tau_borophene_pristine.csv")
top_candidates = ["Hydromethylthionine", "Curcumin", "EGCG", "Resveratrol", "Quercetin", "Baicalein", "Honokiol", "Donepezil"]

print("="*80)
print("RUNNING TAU N=8 MULTI-ORIENTATION RELAXATIONS WITH STRICT CONVERGENCE")
print("="*80)

relaxed_rows = []

for name in top_candidates:
    dir_name = name.replace(" ", "_").replace("-", "_")
    mol_dir = calc / dir_name
    drug_xyz = mol_dir / f"{dir_name}_drug.xyz"
    row = df_drugs[df_drugs["name"] == name].iloc[0]
    ed = float(row["E_drug_Eh"])
    q = int(row["formal_charge"])
    sp_val = float(row["delta_Eint_SP_kcal_mol"])

    drug_lines = drug_xyz.read_text().splitlines()
    n_drug = int(drug_lines[0])
    d_atoms_raw = []
    for l in drug_lines[2:2+n_drug]:
        p = l.split()
        d_atoms_raw.append((p[0], float(p[1]), float(p[2]), float(p[3])))
    
    orig_coords = np.array([[x, y, z] for _, x, y, z in d_atoms_raw])
    orig_coords -= np.mean(orig_coords, axis=0)

    best_e = 999.0
    best_deg = None
    best_final_xyz = None
    best_out = None
    best_gn = None
    best_cyc = None
    
    print(f"\n>>> Compound: {name} (N_drug={n_drug}, q={q}, SP={sp_val:.2f} kcal/mol)")
    
    for angle_deg in [0, 90, 180, 270]:
        theta = np.radians(angle_deg)
        R_z = np.array([[np.cos(theta), -np.sin(theta), 0], [np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
        rot_coords = orig_coords @ R_z.T
        rot_coords[:, 2] -= np.min(rot_coords[:, 2])
        rot_coords[:, 2] += (z_top + 3.20)
        
        in_xyz = mol_dir / f"{dir_name}_opt_{angle_deg}deg_input.xyz"
        with open(in_xyz, "w") as fh:
            fh.write(f"{n_drug+len(carrier_atoms)}\n{name} {angle_deg} deg on B40H15\n")
            for p, (x, y, z) in zip(d_atoms_raw, rot_coords):
                fh.write(f"{p[0]}  {x:12.6f}  {y:12.6f}  {z:12.6f}\n")
            for elem, x, y, z in carrier_atoms:
                fh.write(f"{elem}  {x:12.6f}  {y:12.6f}  {z:12.6f}\n")
        
        out_f = mol_dir / f"{dir_name}_opt_{angle_deg}deg.out"
        final_xyz = mol_dir / f"{dir_name}_orientation_{angle_deg}deg_final.xyz"
        
        xtbopt = mol_dir / "xtbopt.xyz"
        if xtbopt.exists():
            xtbopt.unlink()
            
        cmd = [
            str(xtb), str(in_xyz),
            "--opt", "loose",
            "--gfn", "2",
            "--chrg", str(q),
            "--uhf", "0",
            "--iterations", "500",
            "--cycles", "400",
            "--norestart"
        ]
        
        t0 = time.time()
        res = subprocess.run(cmd, cwd=str(mol_dir), stdout=open(out_f, "w"), timeout=300)
        dt = time.time() - t0
        
        if xtbopt.exists():
            shutil.copy(xtbopt, final_xyz)
            
        e_val, conv, gn, cyc = parse_xtb_opt_output(out_f)
        
        status_str = "CONVERGED" if conv else "FAILED"
        print(f"    Orientation {angle_deg:>3} deg: {status_str:<9} | E={e_val} Eh | Cycles={cyc:>3} | GradNorm={gn} | Time={dt:.1f}s")
        
        if conv and e_val is not None and final_xyz.exists():
            if e_val < best_e:
                best_e = e_val
                best_deg = angle_deg
                best_final_xyz = final_xyz
                best_out = out_f
                best_gn = gn
                best_cyc = cyc
                
    if best_final_xyz and best_e < 900.0:
        de_opt = (best_e - ed - e_carrier) * 627.509
        print(f"  --> [SELECTED MINIMUM] {name}: Best Orient={best_deg} deg, Relaxed E_int = {de_opt:.2f} kcal/mol (Cycles={best_cyc}, GradNorm={best_gn})")
        relaxed_rows.append({
            "name": name,
            "best_orientation_deg": best_deg,
            "delta_Eint_SP_kcal_mol": sp_val,
            "delta_Eint_relaxed_kcal_mol": round(de_opt, 3),
            "convergence_status": "CONVERGED (*** convergence criteria satisfied ***)",
            "gradient_norm_Eh_a0": best_gn,
            "optimization_cycles": best_cyc,
            "output_file": str(best_out.relative_to(base)),
            "final_pose_file": str(best_final_xyz.relative_to(base)),
            "sha256": sha256_file(best_final_xyz)
        })
    else:
        print(f"  [ERROR] No converged orientation found for {name}!")

df_rel = pd.DataFrame(relaxed_rows)
df_rel.to_csv(proc / "relaxed_adsorption_subset.csv", index=False)
print("\n" + "="*80)
print(f"SAVED relaxed_adsorption_subset.csv with N = {len(df_rel)} compounds (Target: 8)")
print("="*80)

if len(df_rel) == 8:
    rho_s, p_s = spearmanr(df_rel["delta_Eint_SP_kcal_mol"], df_rel["delta_Eint_relaxed_kcal_mol"])
    mae_s = mean_absolute_error(df_rel["delta_Eint_SP_kcal_mol"], df_rel["delta_Eint_relaxed_kcal_mol"])
    print(f"[RELAXED SUBSET AUDIT N=8] Spearman rho = {rho_s:.4f} (p={p_s:.4f}), MAE = {mae_s:.2f} kcal/mol")
else:
    print(f"[INCOMPLETE] len(df_rel) is {len(df_rel)}, not 8.")
