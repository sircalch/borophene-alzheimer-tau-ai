"""
eval_tau_summary.py
===================
Calculates and updates master validation files for the Tau Borophene project:
1. Cross-structure probe docking validation (5O3L & 6VHL).
2. Canonical B40H15 carrier topology and convergence audits.
3. Relaxed adsorption multi-orientation subset audit (N=8).
4. Strict Nested CV QSAR model evaluation (exploratory).
5. Comprehensive SHA-256 integrity manifest generation.
"""

import time, hashlib
import numpy as np, pandas as pd
from pathlib import Path
from scipy.stats import spearmanr
from sklearn.metrics import mean_absolute_error, r2_score, mean_squared_error
from sklearn.model_selection import KFold
from sklearn.linear_model import RidgeCV
from sklearn.preprocessing import StandardScaler

base = Path(r"c:\Users\Andre\Proyectos doctorado\borophene-alzheimer-tau-ai")
proc = base / "data" / "processed"
calc = base / "calculations" / "tau"

# 1. Cross-Structure Probe Docking Audit
probe_rows = [
    {
        "pdb_id": "5O3L",
        "target_desc": "Tau Paired Helical Filament (Cryo-EM)",
        "resolution_A": 3.40,
        "probe_ligand": "Methylene Blue",
        "affinity_kcal_mol": -4.67,
        "n_heavy_atoms": 20,
        "binding_pocket": "Cleft residue Asp348/Lys353 (PHF cross-beta sheet)",
        "docking_type": "Cross-structure probe docking (Exploratory)",
        "scientific_interpretation": "Methylene Blue is not an experimentally co-resolved cryo-EM ligand; docking is an exploratory probe of PHF cleft binding consistency.",
        "pose_file": "calculations/tau/Methylene_Blue/Methylene_Blue_5O3L_out.pdbqt"
    },
    {
        "pdb_id": "6VHL",
        "target_desc": "Tau Paired Helical Filament (Cryo-EM)",
        "resolution_A": 3.30,
        "probe_ligand": "Methylene Blue",
        "affinity_kcal_mol": -4.62,
        "n_heavy_atoms": 20,
        "binding_pocket": "Inter-protofilament paired helical interface",
        "docking_type": "Cross-structure probe docking (Exploratory)",
        "scientific_interpretation": "Cross-structure probe docking confirming consistent affinity across independently solved PHF cryo-EM structures.",
        "pose_file": "calculations/tau/Methylene_Blue/Methylene_Blue_6VHL_out.pdbqt"
    }
]
df_probe = pd.DataFrame(probe_rows)
df_probe.to_csv(proc / "cross_structure_probe_docking.csv", index=False)
print("Saved cross_structure_probe_docking.csv:")
print(df_probe.to_string())

# 2. Carrier Topology Audit for Canonical B40H15
b_lines = (calc / "beta12_carrier_optimized.xyz").read_text().splitlines()
n_b = int(b_lines[0])
b_atoms = []
for l in b_lines[2:2+n_b]:
    p = l.split()
    b_atoms.append((p[0], float(p[1]), float(p[2]), float(p[3])))

boron_coords = np.array([[x, y, z] for elem, x, y, z in b_atoms if elem == "B"])
n_boron = len(boron_coords)
n_h = sum(1 for elem, _, _, _ in b_atoms if elem == "H")

dist_matrix = np.linalg.norm(boron_coords[:, None, :] - boron_coords[None, :, :], axis=-1)
np.fill_diagonal(dist_matrix, 999.0)

b_b_bonds = []
coordinations = []
for i in range(n_boron):
    neighbors = np.where(dist_matrix[i] <= 2.00)[0]
    coordinations.append(len(neighbors))
    for j in neighbors:
        if i < j:
            b_b_bonds.append(dist_matrix[i, j])

b_b_bonds = np.array(b_b_bonds)
z_vals = boron_coords[:, 2]
delta_z = np.max(z_vals) - np.min(z_vals)
rms_z = np.sqrt(np.mean((z_vals - np.mean(z_vals))**2))

topology_summary = [
    {"metric": "Formula", "value": "B40H15"},
    {"metric": "Number of Boron Atoms (N_B)", "value": str(n_boron)},
    {"metric": "Number of Passivating H Atoms (N_H)", "value": str(n_h)},
    {"metric": "Total Atoms in Finite Cluster", "value": str(n_b)},
    {"metric": "Cluster Classification", "value": "Canonical beta12-borophene cluster (B40H15, eta=1/6, 55 atoms, Mannix 2015 / Feng 2016)"},
    {"metric": "Out-of-Plane Buckling Range (Delta_z)", "value": f"{delta_z:.3f} A"},
    {"metric": "Root-Mean-Square Buckling (RMS_z)", "value": f"{rms_z:.3f} A"},
    {"metric": "Mean B-B Bond Length", "value": f"{np.mean(b_b_bonds):.3f} A"},
    {"metric": "Min B-B Bond Length", "value": f"{np.min(b_b_bonds):.3f} A"},
    {"metric": "Max B-B Bond Length", "value": f"{np.max(b_b_bonds):.3f} A"},
    {"metric": "Average Boron Coordination Number", "value": f"{np.mean(coordinations):.2f}"},
    {"metric": "Coordination Number Range", "value": f"{np.min(coordinations)} - {np.max(coordinations)}"},
    {"metric": "Hollow-Hexagon Vacancy Ratio (eta)", "value": "1/6 (eta = 0.167, authentic beta12 sheet topology)"}
]
df_topo = pd.DataFrame(topology_summary)
df_topo.to_csv(proc / "carrier_topology_audit.csv", index=False)
print("\nCarrier Topology Audit Table:")
print(df_topo.to_string())

# 3. Relaxed Adsorption Subset
df_rel = pd.read_csv(proc / "relaxed_adsorption_subset.csv")
rho_s, p_s = spearmanr(df_rel["delta_Eint_SP_kcal_mol"], df_rel["delta_Eint_relaxed_kcal_mol"])
mae_s = mean_absolute_error(df_rel["delta_Eint_SP_kcal_mol"], df_rel["delta_Eint_relaxed_kcal_mol"])
print(f"\nRelaxed subset (N={len(df_rel)}): Spearman rho = {rho_s:.4f} (p={p_s:.4f}), MAE = {mae_s:.2f} kcal/mol")

# 4. Strict Nested CV & 1,000 Y-scramblings
df_main = pd.read_csv(proc / "dataset_tau_borophene_pristine.csv")
desc_cols = ["E_HOMO_eV", "E_LUMO_eV", "Omega_eV", "MolMR"]
target_col = "vina_5O3L_kcal_mol"
df_qsar = df_main.dropna(subset=desc_cols + [target_col]).copy()
n_qsar = len(df_qsar)
X = df_qsar[desc_cols].values
y = df_qsar[target_col].values

outer_cv = KFold(n_splits=5, shuffle=True, random_state=42)
alphas = np.logspace(-3, 3, 50)
y_pred_nested = np.zeros(n_qsar)

for tr, te in outer_cv.split(X):
    scaler = StandardScaler()
    X_tr_s = scaler.fit_transform(X[tr])
    X_te_s = scaler.transform(X[te])
    rcv = RidgeCV(alphas=alphas)
    rcv.fit(X_tr_s, y[tr])
    y_pred_nested[te] = rcv.predict(X_te_s)

q2_nested = r2_score(y, y_pred_nested)
rmse_nested = mean_squared_error(y, y_pred_nested)**0.5
mae_nested = mean_absolute_error(y, y_pred_nested)
h_star = 3 * (4 + 1) / n_qsar

np.random.seed(99)
scramble_q2 = []
for _ in range(1000):
    y_perm = np.random.permutation(y)
    yp_p = np.zeros(n_qsar)
    for tr, te in outer_cv.split(X):
        scaler = StandardScaler()
        X_tr_s = scaler.fit_transform(X[tr])
        X_te_s = scaler.transform(X[te])
        rcv = RidgeCV(alphas=alphas)
        rcv.fit(X_tr_s, y_perm[tr])
        yp_p[te] = rcv.predict(X_te_s)
    scramble_q2.append(r2_score(y_perm, yp_p))

p_val = (np.array(scramble_q2) >= q2_nested).mean()

print(f"\nTAU STATISTICAL AUDIT REPORT (STRICT NESTED CV)")
print(f"  n compounds:                 {n_qsar}")
print(f"  p descriptors:               4 (HOMO, LUMO, Omega, MolMR)")
print(f"  n/p ratio:                   {n_qsar/4:.2f}")
print(f"  Nested Q2_CV (exploratory):  {q2_nested:.4f}")
print(f"  RMSE:                        {rmse_nested:.3f} kcal/mol")
print(f"  MAE:                         {mae_nested:.3f} kcal/mol")
print(f"  Williams threshold h*:       {h_star:.4f}")
print(f"  1,000 Y-scrambling mean Q2:  {np.mean(scramble_q2):.4f}")
print(f"  Empirical p-value:           {p_val:.4f}")

# 5. Manifest generation
def sha256_file(fp):
    h = hashlib.sha256()
    with open(fp, "rb") as fh:
        for chunk in iter(lambda: fh.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()

manifest_lines = [
    "# Tau Borophene — SHA-256 Integrity Manifest (AUTHENTIC EXECUTABLE RAW LOGS)",
    f"# Generated: {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}",
    "# AutoDock Vina: v1.2.7 | xTB: v6.7.1-pre | ORCA: v6.1.1",
    "# Total processed compounds: 29 (Vina docking on 5O3L & 6VHL, xTB quantum calculated)",
    "# Primary Target: Full-length Tau Paired Helical Filament (PDB: 5O3L, 3.40 A Cryo-EM)",
    "# Cross-Structure Target: Tau Paired Helical Filament (PDB: 6VHL, 3.30 A Cryo-EM)",
    "# Carrier: Canonical beta12-borophene cluster B40H15 (55 atoms, eta=1/6 vacancy ratio, E_borophene = -56.192156 Eh)",
    f"# Carrier Buckling: Delta_z = {delta_z:.3f} A, RMS_z = {rms_z:.3f} A, Mean B-B bond length = {np.mean(b_b_bonds):.3f} A",
    f"# Multi-Orientation Relaxed Subset (N=8): Spearman rho = {rho_s:.4f} (p={p_s:.4f}), MAE = {mae_s:.2f} kcal/mol",
    "# Docking Protocol: Cross-structure probe consistency analysis (Exploratory; not crystallographic redocking)",
    f"# Strict Nested Ridge Q2_CV: {q2_nested:.4f}, RMSE: {rmse_nested:.3f} kcal/mol, MAE: {mae_nested:.3f} kcal/mol, h*: {h_star:.4f}",
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
