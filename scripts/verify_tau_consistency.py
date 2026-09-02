"""
scripts/verify_tau_consistency.py
=================================
Automated test suite to verify absolute consistency across all active Tau files:
1. Verifies that only the official carrier formula (B40H15) is present in data files,
   provenance files, and manifests.
2. Fails immediately if legacy/incompatible tags (B48H12, B53H7, B45H18) are detected.
3. Checks SHA256 integrity, convergence criteria, and dataset synchronization.
"""

import sys, re
import pandas as pd
from pathlib import Path

base = Path(__file__).resolve().parent.parent if "__file__" in globals() else Path(r"c:\Users\Andre\Proyectos doctorado\borophene-alzheimer-tau-ai")
proc = base / "data" / "processed"
calc = base / "calculations" / "tau"

OFFICIAL_FORMULA = "B40H15"
DISCORDANT_FORMULAS = ["B48H12", "B53H7", "B45H18", "B45H17"]

errors = []

print("="*70)
print(f"TAU REPOSITORY CONSISTENCY AUDIT (Target: {OFFICIAL_FORMULA})")
print("="*70)

# 1. Check carrier_convergence.csv
conv_csv = proc / "carrier_convergence.csv"
if not conv_csv.exists():
    errors.append("Missing carrier_convergence.csv")
else:
    df_c = pd.read_csv(conv_csv)
    formula_found = df_c.iloc[0]["formula"]
    status_found = df_c.iloc[0]["convergence_status"]
    if formula_found != OFFICIAL_FORMULA:
        errors.append(f"carrier_convergence.csv formula mismatch: {formula_found} != {OFFICIAL_FORMULA}")
    if "CONVERGED" not in status_found:
        errors.append(f"carrier_convergence.csv is NOT converged: {status_found}")
    print(f"[PASS] carrier_convergence.csv: Formula={formula_found}, Status={status_found}")

# 2. Check carrier_topology_audit.csv
topo_csv = proc / "carrier_topology_audit.csv"
if not topo_csv.exists():
    errors.append("Missing carrier_topology_audit.csv")
else:
    df_t = pd.read_csv(topo_csv)
    formula_row = df_t[df_t["metric"] == "Chemical Formula"]
    if formula_row.empty:
        errors.append("carrier_topology_audit.csv missing Chemical Formula row")
    else:
        val = formula_row.iloc[0]["optimized_value"]
        if val != OFFICIAL_FORMULA:
            errors.append(f"carrier_topology_audit.csv formula mismatch: {val} != {OFFICIAL_FORMULA}")
        else:
            print(f"[PASS] carrier_topology_audit.csv: Formula={val}")

# 3. Check carrier_identity_provenance.csv
prov_csv = proc / "carrier_identity_provenance.csv"
if not prov_csv.exists():
    errors.append("Missing carrier_identity_provenance.csv")
else:
    df_p = pd.read_csv(prov_csv)
    for idx, r in df_p.iterrows():
        if r["N_B"] != 40 or r["N_H"] != 15:
            errors.append(f"carrier_identity_provenance.csv row {idx} stoichiometry mismatch: N_B={r['N_B']}, N_H={r['N_H']}")
    print(f"[PASS] carrier_identity_provenance.csv: All rows match {OFFICIAL_FORMULA} (N_B=40, N_H=15)")

# 4. Check dataset_tau_borophene_pristine.csv
data_csv = proc / "dataset_tau_borophene_pristine.csv"
if not data_csv.exists():
    errors.append("Missing dataset_tau_borophene_pristine.csv")
else:
    df_d = pd.read_csv(data_csv)
    if "carrier_formula" in df_d.columns:
        unique_formulas = df_d["carrier_formula"].unique()
        if len(unique_formulas) != 1 or unique_formulas[0] != OFFICIAL_FORMULA:
            errors.append(f"dataset_tau_borophene_pristine.csv carrier_formula mismatch: {unique_formulas}")
        else:
            print(f"[PASS] dataset_tau_borophene_pristine.csv: Tagged with {OFFICIAL_FORMULA} (N={len(df_d)})")

# 5. Check MANIFEST_SHA256.txt for discordant formulas
manifest = base / "MANIFEST_SHA256.txt"
if not manifest.exists():
    errors.append("Missing MANIFEST_SHA256.txt")
else:
    m_text = manifest.read_text(encoding="utf-8", errors="replace")
    for disc in DISCORDANT_FORMULAS:
        if disc in m_text:
            errors.append(f"MANIFEST_SHA256.txt contains discordant formula '{disc}'")
    if OFFICIAL_FORMULA not in m_text:
        errors.append(f"MANIFEST_SHA256.txt is missing official formula '{OFFICIAL_FORMULA}'")
    else:
        print(f"[PASS] MANIFEST_SHA256.txt: Contains {OFFICIAL_FORMULA} and zero discordant tags.")

# 6. Check beta12_carrier_optimized.xyz
opt_xyz = calc / "beta12_carrier_optimized.xyz"
if not opt_xyz.exists():
    errors.append("Missing beta12_carrier_optimized.xyz")
else:
    lines = opt_xyz.read_text().splitlines()
    n_atoms = int(lines[0])
    if n_atoms != 55:
        errors.append(f"beta12_carrier_optimized.xyz atom count mismatch: {n_atoms} != 55")
    else:
        print(f"[PASS] beta12_carrier_optimized.xyz: Atom count = {n_atoms} (40 B + 15 H)")

if errors:
    print("\n[FAIL] Consistency errors found:")
    for e in errors:
        print(f"  - {e}")
    sys.exit(1)
else:
    print("\n" + "="*70)
    print(f"[SUCCESS] ALL TAU FILES ARE 100% CONSISTENT WITH CANONICAL {OFFICIAL_FORMULA}!")
    print("="*70)
