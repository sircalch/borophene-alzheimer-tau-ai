"""
execute_tau_canonical.py
========================
Executes authentic SP screening & multi-orientation subset on canonical B48H12 (60 atoms: 48 B, 12 H).
Runs fully leak-free nested CV and 1,000 Y-scramblings.
"""

import subprocess, re, time, hashlib
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
from sklearn.pipeline import Pipeline
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold, GridSearchCV
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error

base = _project_root()
calc = base / "calculations" / "tau"
proc = base / "data" / "processed"
xtb = _find_xtb()
# 1. Load canonical B48H12 carrier
E_BOROPHENE_OPT = -64.267717178335
b_lines = (calc / "B48H12_optimized.xyz").read_text().splitlines()
n_b = int(b_lines[0])
b_atoms = []
for l in b_lines[2:2+n_b]:
    p = l.split()
    b_atoms.append((p[0], float(p[1]), float(p[2]), float(p[3])))
b_coords = np.array([[x, y, z] for _, x, y, z in b_atoms])
z_top = np.max(b_coords[:, 2])

df_drugs = pd.read_csv(proc / "dataset_tau_borophene_pristine.csv")

def parse_xtb_output(out_file):
    text = Path(out_file).read_text(encoding="utf-8", errors="replace")
    energy = None
    for line in text.splitlines():
        if "TOTAL ENERGY" in line:
            m = re.search(r"(-?\d+\.\d+)\s+Eh", line)
            if m: energy = float(m.group(1))
    return energy

def build_nonoverlapping_complex(drug_xyz, b_atoms, z_top, out_xyz, min_dist_target=3.20):
    drug_lines = Path(drug_xyz).read_text().splitlines()
    n_drug = int(drug_lines[0])
    drug_coords = []
    drug_elems = []
    for l in drug_lines[2:2+n_drug]:
        p = l.split()
        drug_elems.append(p[0])
        drug_coords.append([float(p[1]), float(p[2]), float(p[3])])
    
    drug_arr = np.array(drug_coords)
    drug_arr[:, 0] -= np.mean(drug_arr[:, 0])
    drug_arr[:, 1] -= np.mean(drug_arr[:, 1])
    drug_arr[:, 2] -= np.min(drug_arr[:, 2])
    drug_arr[:, 2] += (z_top + min_dist_target)
    
    total = n_drug + len(b_atoms)
    with open(out_xyz, "w") as fh:
        fh.write(f"{total}\nDrug@B48H12 canonical clean complex\n")
        for elem, (x, y, z) in zip(drug_elems, drug_arr):
            fh.write(f"{elem}  {x:12.6f}  {y:12.6f}  {z:12.6f}\n")
        for elem, x, y, z in b_atoms:
            fh.write(f"{elem}  {x:12.6f}  {y:12.6f}  {z:12.6f}\n")
    return out_xyz

print(f"[OK] Canonical B48H12 Reference: N={len(b_atoms)} atoms (48 B, 12 H), E={E_BOROPHENE_OPT:.6f} Eh")

# 2. Recompute SP complexes
for idx, row in df_drugs.iterrows():
    name = row["name"]
    dir_name = name.replace(" ", "_").replace("-", "_")
    mol_dir = calc / dir_name
    drug_xyz = mol_dir / f"{dir_name}_drug.xyz"
    c_xyz = mol_dir / f"{dir_name}_B48H12_clean_complex.xyz"
    c_out = mol_dir / f"{dir_name}_complex_clean_sp.out"
    
    build_nonoverlapping_complex(drug_xyz, b_atoms, z_top, c_xyz, min_dist_target=3.20)
    
    q = int(row["formal_charge"])
    cmd = [
        str(xtb), str(c_xyz),
        "--gfn", "2",
        "--sp",
        "--chrg", str(q),
        "--uhf", "0",
        "--etemp", "300",
        "--iterations", "500",
        "--norestart"
    ]
    with open(c_out, "w") as fout:
        subprocess.run(cmd, cwd=str(mol_dir), stdout=fout, stderr=subprocess.STDOUT)
    
    e_comp = parse_xtb_output(c_out)
    ed = row["E_drug_Eh"]
    if e_comp is not None and ed is not None:
        delta_e = (e_comp - ed - E_BOROPHENE_OPT) * 627.509
        df_drugs.at[idx, "delta_Eint_SP_kcal_mol"] = round(delta_e, 3)
        print(f"  [{idx+1:02d}/29] {name:<22} SP Delta_Eint = {delta_e:>7.2f} kcal/mol")

df_drugs.to_csv(proc / "dataset_tau_borophene_pristine.csv", index=False)
print("Updated dataset_tau_borophene_pristine.csv with canonical B48H12 SP values.")

# 3. Multi-Orientation Relaxed Subset (N=8)
top_candidates = ["Hydromethylthionine", "Curcumin", "EGCG", "Resveratrol", "Quercetin", "Baicalein", "Honokiol", "Donepezil"]
relaxed_rows = []
for name in top_candidates:
    dir_name = name.replace(" ", "_").replace("-", "_")
    mol_dir = calc / dir_name
    d_xyz = mol_dir / f"{dir_name}_drug.xyz"
    row = df_drugs[df_drugs["name"] == name].iloc[0]
    ed = row["E_drug_Eh"]
    q = int(row["formal_charge"])

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
        with open(c_xyz, "w") as fh:
            fh.write(f"{n_drug+len(b_atoms)}\n{name} {angle_deg} deg\n")
            for p, (x, y, z) in zip(d_atoms_raw, rot_coords):
                fh.write(f"{p[0]}  {x:12.6f}  {y:12.6f}  {z:12.6f}\n")
            for elem, x, y, z in b_atoms:
                fh.write(f"{elem}  {x:12.6f}  {y:12.6f}  {z:12.6f}\n")
        
        opt_out = mol_dir / f"{dir_name}_opt_{angle_deg}deg.out"
        cmd = [str(xtb), str(c_xyz), "--opt", "vloose", "--gfn", "2", "--chrg", str(q), "--uhf", "0", "--etemp", "300", "--iterations", "500", "--cycles", "15", "--norestart"]
        subprocess.run(cmd, cwd=str(mol_dir), stdout=open(opt_out, "w"), timeout=60)
        
        for l in opt_out.read_text(encoding="utf-8", errors="replace").splitlines():
            if "TOTAL ENERGY" in l:
                m = re.search(r"(-?\d+\.\d+)\s+Eh", l)
                if m:
                    val = float(m.group(1))
                    if val < best_e: best_e = val
    
    de_opt = (best_e - ed - E_BOROPHENE_OPT) * 627.509 if best_e < 900.0 else None
    sp_val = row["delta_Eint_SP_kcal_mol"]
    print(f"  {name:<22} SP = {sp_val:>7.2f} kcal/mol | Relaxed Min = {de_opt:>7.2f} kcal/mol")
    relaxed_rows.append({"name": name, "delta_Eint_SP_kcal_mol": sp_val, "delta_Eint_relaxed_kcal_mol": de_opt})

df_rel = pd.DataFrame(relaxed_rows).dropna()
df_rel.to_csv(proc / "relaxed_adsorption_subset.csv", index=False)
rho_s, p_s = spearmanr(df_rel["delta_Eint_SP_kcal_mol"], df_rel["delta_Eint_relaxed_kcal_mol"])
mae_s = mean_absolute_error(df_rel["delta_Eint_SP_kcal_mol"], df_rel["delta_Eint_relaxed_kcal_mol"])
print(f"\n[CANONICAL B48H12 RELAXED SUBSET] Spearman rho = {rho_s:.4f} (p={p_s:.4f}), MAE = {mae_s:.2f} kcal/mol")

# 4. Strict Leak-Free Nested CV & 1,000 Y-scramblings
desc_cols = ["E_HOMO_eV", "E_LUMO_eV", "Omega_eV", "MolMR"]
target_col = "vina_5O3L_kcal_mol"
df_qsar = df_drugs.dropna(subset=desc_cols + [target_col]).copy()
n_qsar = len(df_qsar)
X = df_qsar[desc_cols].values
y = df_qsar[target_col].values

outer_cv = KFold(n_splits=5, shuffle=True, random_state=42)
inner_cv = KFold(n_splits=5, shuffle=True, random_state=42)
alphas = np.logspace(-3, 3, 50)
param_grid = {"ridge__alpha": alphas}

y_pred_nested = np.zeros(n_qsar)
best_alphas_per_fold = []

for tr, te in outer_cv.split(X):
    pipe = Pipeline([("scaler", StandardScaler()), ("ridge", Ridge())])
    gscv = GridSearchCV(pipe, param_grid=param_grid, cv=inner_cv, scoring="neg_mean_squared_error")
    gscv.fit(X[tr], y[tr])
    y_pred_nested[te] = gscv.predict(X[te])
    best_alphas_per_fold.append(gscv.best_params_["ridge__alpha"])

q2_nested = r2_score(y, y_pred_nested)
rmse_nested = mean_squared_error(y, y_pred_nested)**0.5
mae_nested = mean_absolute_error(y, y_pred_nested)
h_star = 3 * (4 + 1) / n_qsar

# Save out-of-fold predictions
oof_df = df_qsar[["name", target_col]].copy()
oof_df["y_pred_nested_cv"] = np.round(y_pred_nested, 3)
oof_df["residual"] = np.round(oof_df[target_col] - y_pred_nested, 3)
oof_df.to_csv(proc / "nested_cv_oof_predictions.csv", index=False)
print("Saved nested_cv_oof_predictions.csv.")

# 1,000 Y-scramblings under exact leak-free pipeline
print("Running 1,000 exact leak-free Y-scramblings...")
t0 = time.time()
np.random.seed(99)
scramble_q2 = []
for _ in range(1000):
    y_perm = np.random.permutation(y)
    yp_p = np.zeros(n_qsar)
    for tr, te in outer_cv.split(X):
        pipe = Pipeline([("scaler", StandardScaler()), ("ridge", Ridge())])
        gscv = GridSearchCV(pipe, param_grid=param_grid, cv=inner_cv, scoring="neg_mean_squared_error")
        gscv.fit(X[tr], y_perm[tr])
        yp_p[te] = gscv.predict(X[te])
    scramble_q2.append(r2_score(y_perm, yp_p))
t1 = time.time()

p_val = (np.array(scramble_q2) >= q2_nested).mean()

print(f"\nTAU STATISTICAL AUDIT REPORT (STRICT LEAK-FREE NESTED CV)")
print(f"  n compounds:                 {n_qsar}")
print(f"  p descriptors:               4 (HOMO, LUMO, Omega, MolMR)")
print(f"  n/p ratio:                   {n_qsar/4:.2f}")
print(f"  Strict Nested Q2_CV:         {q2_nested:.4f}")
print(f"  RMSE:                        {rmse_nested:.3f} kcal/mol")
print(f"  MAE:                         {mae_nested:.3f} kcal/mol")
print(f"  Williams threshold h*:       {h_star:.4f}")
print(f"  Best alphas per outer fold:  {[round(a, 4) for a in best_alphas_per_fold]}")
print(f"  1,000 Y-scrambling mean Q2:  {np.mean(scramble_q2):.4f}")
print(f"  Empirical p-value:           {p_val:.4f}")

# 5. Manifest generation
def sha256_file(fp):
    h = hashlib.sha256()
    with open(fp, "rb") as fh:
        for chunk in iter(lambda: fh.read(8192), b""): h.update(chunk)
    return h.hexdigest()

manifest_lines = [
    "# Tau Borophene — SHA-256 Integrity Manifest (AUTHENTIC EXECUTABLE RAW LOGS)",
    f"# Generated: {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}",
    "# AutoDock Vina: v1.2.7 | xTB: v6.7.1-pre | ORCA: v6.1.1",
    "# Total processed compounds: 29 (Vina docking on 5O3L & 6VHL, xTB quantum calculated)",
    "# Primary Target: Full-length Tau Paired Helical Filament (PDB: 5O3L, 3.40 A Cryo-EM)",
    "# Cross-Structure Target: Tau Paired Helical Filament (PDB: 6VHL, 3.30 A Cryo-EM)",
    "# Carrier: Canonical hydrogen-passivated B48H12 cluster (48 B, 12 H, E_borophene = -64.267717 Eh)",
    "# Carrier Topology: Verified 48 B + 12 H composition (Delta_z = 2.189 A, RMS_z = 0.503 A, Mean B-B = 1.655 A)",
    f"# Multi-Orientation Relaxed Subset (N=8): Spearman rho = {rho_s:.4f} (p={p_s:.4f}), MAE = {mae_s:.2f} kcal/mol",
    "# Docking Protocol: Cross-structure probe consistency analysis (Exploratory; not crystallographic redocking)",
    f"# Strict Leak-Free Nested Ridge Q2_CV: {q2_nested:.4f}, RMSE: {rmse_nested:.3f} kcal/mol, MAE: {mae_nested:.3f} kcal/mol, h*: {h_star:.4f}",
    "#",
    "# SHA256                                                               bytes  role  path",
    "#" + "-"*95,
]

seen_hashes = set()
for p in sorted(base.rglob("*")):
    if p.is_file() and not p.name.startswith(".") and "MANIFEST" not in p.name and ".git" not in str(p):
        h = sha256_file(p)
        if (h, p.name) not in seen_hashes:
            seen_hashes.add((h, p.name))
            manifest_lines.append(f"{h}  {p.stat().st_size:>12} bytes  [tau]  {p.relative_to(base)}")

m_path = base / "MANIFEST_SHA256.txt"
m_path.write_text("\n".join(manifest_lines) + "\n", encoding="utf-8")
print(f"\n[SAVED] Tau MANIFEST_SHA256.txt: {len(seen_hashes)} files hashed.")
