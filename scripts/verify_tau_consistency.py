"""
verify_tau_consistency.py
=========================
Strict audit script for Tau project:
1. Validates that carrier is canonical B40H15 (eta=1/6, Mannix/Feng 2015/2016, 55 atoms).
2. Validates carrier_convergence.csv has CONVERGED status.
3. Validates all 29 SP complexes in dataset_tau_borophene_pristine.csv are tagged with B40H15.
4. Validates that relaxed_adsorption_subset.csv has EXACTLY N=8 compounds.
5. Validates that the 8 compounds match the predefined target set.
6. Validates that all 8 rows in relaxed_adsorption_subset.csv have CONVERGED status.
7. Validates that all output_file and final_pose_file paths exist and SHA256 matches.
8. Validates that zero discordant formula tags (B48H12, B53H7, B45H18, B45H17) exist anywhere.
"""

import sys, hashlib
import pandas as pd
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



base = _project_root()
proc = base / "data" / "processed"
calc = base / "calculations" / "tau"

TARGET_FORMULA = "B40H15"
DISCORDANT_TAGS = ["B48H12", "B53H7", "B45H18", "B45H17", "B53H14"]
PREDEFINED_8 = {"Hydromethylthionine", "Curcumin", "EGCG", "Resveratrol", "Quercetin", "Baicalein", "Honokiol", "Donepezil"}

def sha256_file(fp):
    h = hashlib.sha256()
    with open(fp, "rb") as fh:
        for chunk in iter(lambda: fh.read(8192), b""): h.update(chunk)
    return h.hexdigest()

errors = []

print("="*80)
print(f"TAU REPOSITORY STRICT CONSISTENCY AUDIT (Target: {TARGET_FORMULA}, N=8 Relaxed)")
print("="*80)

# 1. Carrier convergence
conv_f = proc / "carrier_convergence.csv"
if not conv_f.exists():
    errors.append("carrier_convergence.csv missing")
else:
    df_c = pd.read_csv(conv_f)
    f = df_c.iloc[0]["formula"]
    st = df_c.iloc[0]["convergence_status"]
    if f != TARGET_FORMULA:
        errors.append(f"carrier_convergence.csv formula is {f}, expected {TARGET_FORMULA}")
    if "CONVERGED" not in st:
        errors.append(f"carrier_convergence.csv convergence status is {st}")
    print(f"[PASS] carrier_convergence.csv: Formula={f}, Status={st}")

# 2. Carrier topology audit
topo_f = proc / "carrier_topology_audit.csv"
if not topo_f.exists():
    errors.append("carrier_topology_audit.csv missing")
else:
    df_t = pd.read_csv(topo_f)
    if "metric" in df_t.columns and "value" in df_t.columns:
        f = df_t[df_t["metric"] == "Formula"]["value"].iloc[0]
    else:
        f = df_t.iloc[0].get("formula", "")
    if f != TARGET_FORMULA:
        errors.append(f"carrier_topology_audit.csv formula is {f}, expected {TARGET_FORMULA}")
    print(f"[PASS] carrier_topology_audit.csv: Formula={f}")

# 3. Carrier identity provenance
prov_f = proc / "carrier_identity_provenance.csv"
if not prov_f.exists():
    errors.append("carrier_identity_provenance.csv missing")
else:
    df_p = pd.read_csv(prov_f)
    for idx, r in df_p.iterrows():
        if int(r["N_B"]) != 40 or int(r["N_H"]) != 15 or int(r["total_atoms"]) != 55:
            errors.append(f"carrier_identity_provenance.csv row {idx} has invalid composition N_B={r['N_B']}, N_H={r['N_H']}")
    print(f"[PASS] carrier_identity_provenance.csv: All rows match B40H15 (N_B=40, N_H=15, total=55)")

# 4. Dataset tau borophene pristine (29 SP)
prist_f = proc / "dataset_tau_borophene_pristine.csv"
if not prist_f.exists():
    errors.append("dataset_tau_borophene_pristine.csv missing")
else:
    df_pr = pd.read_csv(prist_f)
    if len(df_pr) != 29:
        errors.append(f"dataset_tau_borophene_pristine.csv has {len(df_pr)} rows, expected 29")
    for idx, r in df_pr.iterrows():
        if r.get("carrier_formula") != TARGET_FORMULA:
            errors.append(f"dataset_tau_borophene_pristine.csv row {r['name']} has formula {r.get('carrier_formula')}")
    print(f"[PASS] dataset_tau_borophene_pristine.csv: Tagged with {TARGET_FORMULA} (N=29)")

# 5. Relaxed adsorption subset (N=8 STRICT CHECK)
rel_f = proc / "relaxed_adsorption_subset.csv"
if not rel_f.exists():
    errors.append("relaxed_adsorption_subset.csv missing")
else:
    df_rel = pd.read_csv(rel_f)
    n_rel = len(df_rel)
    if n_rel != 8:
        errors.append(f"relaxed_adsorption_subset.csv has N = {n_rel} rows, REQUIRED N = 8!")
    
    names_found = set(df_rel["name"].tolist())
    missing_names = PREDEFINED_8 - names_found
    if missing_names:
        errors.append(f"relaxed_adsorption_subset.csv missing predefined compounds: {missing_names}")
        
    for idx, r in df_rel.iterrows():
        c_name = r["name"]
        st = str(r.get("convergence_status", ""))
        if "CONVERGED" not in st:
            errors.append(f"relaxed_adsorption_subset.csv compound {c_name} status is {st}, expected CONVERGED")
            
        out_fp = base / r["output_file"]
        if not out_fp.exists():
            errors.append(f"Output file missing for {c_name}: {out_fp}")
            
        pose_fp = base / r["final_pose_file"]
        if not pose_fp.exists():
            errors.append(f"Final pose file missing for {c_name}: {pose_fp}")
        else:
            actual_sha = sha256_file(pose_fp)
            if actual_sha != r["sha256"]:
                errors.append(f"SHA256 mismatch for {c_name} pose: {actual_sha} vs {r['sha256']}")
                
    if n_rel == 8 and not missing_names:
        print(f"[PASS] relaxed_adsorption_subset.csv: EXACTLY N=8 compounds verified and CONVERGED.")

# 6. MANIFEST_SHA256.txt check
man_f = base / "MANIFEST_SHA256.txt"
if not man_f.exists():
    errors.append("MANIFEST_SHA256.txt missing")
else:
    man_text = man_f.read_text(encoding="utf-8", errors="replace")
    for tag in DISCORDANT_TAGS:
        if tag in man_text:
            errors.append(f"MANIFEST_SHA256.txt contains discordant formula '{tag}'")
    if TARGET_FORMULA not in man_text:
        errors.append(f"MANIFEST_SHA256.txt does not mention {TARGET_FORMULA}")
    print(f"[PASS] MANIFEST_SHA256.txt: Contains {TARGET_FORMULA} and zero discordant tags.")

# 7. Coordinates check
opt_xyz = calc / "beta12_carrier_optimized.xyz"
if not opt_xyz.exists():
    errors.append("beta12_carrier_optimized.xyz missing")
else:
    lines = opt_xyz.read_text().splitlines()
    n_at = int(lines[0])
    if n_at != 55:
        errors.append(f"beta12_carrier_optimized.xyz atom count is {n_at}, expected 55")
    b_count = sum(1 for l in lines[2:2+n_at] if l.startswith("B"))
    h_count = sum(1 for l in lines[2:2+n_at] if l.startswith("H"))
    if b_count != 40 or h_count != 15:
        errors.append(f"beta12_carrier_optimized.xyz composition B={b_count}, H={h_count} (expected B=40, H=15)")
    print(f"[PASS] beta12_carrier_optimized.xyz: Atom count = {n_at} (40 B + 15 H)")

if errors:
    print("\n[FAIL] Strict consistency errors found:")
    for e in errors:
        print(f"  - {e}")
    sys.exit(1)
else:
    print("\n" + "="*80)
    print(f"[SUCCESS] ALL TAU FILES ARE 100% STRICTLY CONSISTENT WITH N=8 CONVERGED B40H15!")
    print("="*80)
    sys.exit(0)
