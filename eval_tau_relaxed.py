import subprocess, re, numpy as np, pandas as pd
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
e_bor = -67.658968

df = pd.read_csv(proc / "dataset_tau_borophene_pristine.csv")
top_candidates = ["Hydromethylthionine", "Curcumin", "EGCG", "Resveratrol", "Quercetin", "Baicalein", "Honokiol", "Donepezil"]

b_lines = (calc / "B48H12_optimized.xyz").read_text().splitlines()
n_b = int(b_lines[0])
b_atoms = []
for l in b_lines[2:2+n_b]:
    p = l.split()
    b_atoms.append((p[0], float(p[1]), float(p[2]), float(p[3])))
b_coords = np.array([[x, y, z] for _, x, y, z in b_atoms])
z_top = np.max(b_coords[:, 2])

relaxed_rows = []
for name in top_candidates:
    dir_name = name.replace(" ", "_").replace("-", "_")
    mol_dir = calc / dir_name
    d_xyz = mol_dir / f"{dir_name}_drug.xyz"
    row = df[df["name"] == name].iloc[0]
    ed = row["E_drug_Eh"]
    q = row["formal_charge"]

    drug_lines = d_xyz.read_text().splitlines()
    n_drug = int(drug_lines[0])
    d_atoms_raw = []
    for l in drug_lines[2:2+n_drug]:
        p = l.split()
        d_atoms_raw.append((p[0], float(p[1]), float(p[2]), float(p[3])))
    
    orig_coords = np.array([[x, y, z] for _, x, y, z in d_atoms_raw])
    orig_coords -= np.mean(orig_coords, axis=0)

    best_e = 999.0
    for angle_deg in [0, 90, 180, 270]:
        theta = np.radians(angle_deg)
        R_z = np.array([[np.cos(theta), -np.sin(theta), 0], [np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
        rot_coords = orig_coords @ R_z.T
        rot_coords[:, 2] -= np.min(rot_coords[:, 2])
        rot_coords[:, 2] += (z_top + 3.20)
        
        c_xyz = mol_dir / f"{dir_name}_opt_{angle_deg}deg.xyz"
        if not c_xyz.exists():
            with open(c_xyz, "w") as fh:
                fh.write(f"{n_drug+len(b_atoms)}\n{name} {angle_deg} deg\n")
                for p, (x, y, z) in zip(d_atoms_raw, rot_coords):
                    fh.write(f"{p[0]}  {x:12.6f}  {y:12.6f}  {z:12.6f}\n")
                for elem, x, y, z in b_atoms:
                    fh.write(f"{elem}  {x:12.6f}  {y:12.6f}  {z:12.6f}\n")
        
        opt_out = mol_dir / f"{dir_name}_opt_{angle_deg}deg.out"
        if not opt_out.exists() or "TOTAL ENERGY" not in opt_out.read_text(encoding="utf-8", errors="replace"):
            cmd = [str(xtb), str(c_xyz), "--opt", "vloose", "--gfn", "2", "--chrg", str(q), "--uhf", "0", "--etemp", "300", "--iterations", "500", "--cycles", "15", "--norestart"]
            subprocess.run(cmd, cwd=str(mol_dir), stdout=open(opt_out, "w"), timeout=60)
        
        for l in opt_out.read_text(encoding="utf-8", errors="replace").splitlines():
            if "TOTAL ENERGY" in l:
                m = re.search(r"(-?\d+\.\d+)\s+Eh", l)
                if m:
                    val = float(m.group(1))
                    if val < best_e: best_e = val
    
    de_opt = (best_e - ed - e_bor) * 627.509 if best_e < 900.0 else None
    sp_val = row["delta_Eint_SP_kcal_mol"]
    print(f"{name:<22} SP = {sp_val:>7.2f} kcal/mol | Relaxed Min = {de_opt:>7.2f} kcal/mol")
    relaxed_rows.append({"name": name, "delta_Eint_SP_kcal_mol": sp_val, "delta_Eint_relaxed_kcal_mol": de_opt})

df_rel = pd.DataFrame(relaxed_rows).dropna()
df_rel.to_csv(proc / "relaxed_adsorption_subset.csv", index=False)
rho_s, p_s = spearmanr(df_rel["delta_Eint_SP_kcal_mol"], df_rel["delta_Eint_relaxed_kcal_mol"])
mae_s = mean_absolute_error(df_rel["delta_Eint_SP_kcal_mol"], df_rel["delta_Eint_relaxed_kcal_mol"])
print(f"\n[TAU RELAXED SUBSET VALIDATION] Spearman rho = {rho_s:.4f} (p={p_s:.4f}), MAE = {mae_s:.2f} kcal/mol")
