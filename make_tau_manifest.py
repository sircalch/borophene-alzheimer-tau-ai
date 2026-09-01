import hashlib, time, pandas as pd
from pathlib import Path

base = Path(r"c:\Users\Andre\Proyectos doctorado\borophene-alzheimer-tau-ai")
calc = base / "calculations" / "tau"
proc = base / "data" / "processed"

# 1. Calculation provenance
df = pd.read_csv(proc / "dataset_tau_borophene_pristine.csv")
prov_rows = []
for name in df["name"]:
    dir_name = name.replace(" ", "_").replace("-", "_")
    mol_dir = calc / dir_name
    out_f = mol_dir / f"{dir_name}_drug_sp.out"
    log_f = mol_dir / f"{dir_name}_5O3L_vina.log"
    c_out = mol_dir / f"{dir_name}_complex_clean_sp.out"
    prov_rows.append({
        "compound": name,
        "drug_xtb_log": str(out_f.relative_to(base)) if out_f.exists() else "N/A",
        "vina_5o3l_log": str(log_f.relative_to(base)) if log_f.exists() else "N/A",
        "complex_sp_log": str(c_out.relative_to(base)) if c_out.exists() else "N/A"
    })

pd.DataFrame(prov_rows).to_csv(proc / "calculation_provenance.csv", index=False)
print("Provenance CSV written.")

# 2. SHA-256 Manifest
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
    "# Total processed compounds: 29 (Vina docking on 5O3L cryo-EM 2.90 A, xTB quantum calculated)",
    "# Primary Target: Full-length Tau Paired Helical Filament (PDB: 5O3L, 2.90 A Cryo-EM)",
    "# Secondary Target: Tau Straight Filament (PDB: 6VHL, 3.40 A Cryo-EM)",
    "# Carrier: Finite hydrogen-passivated B48H12 boron cluster (60 atoms, beta12 motif, E_borophene = -67.658968 Eh)",
    "# Out-of-plane buckling: Delta_z = 5.128 A, RMS buckling = 1.233 A, Mean B-B bond length = 1.670 A",
    "# Multi-Orientation Relaxed Subset (N=8): Spearman rho = 0.9524 (p=0.0003), MAE = 14.13 kcal/mol",
    "# Nested Ridge Q2_CV (exploratory): -0.1282, RMSE: 0.667 kcal/mol, MAE: 0.510 kcal/mol, h*: 0.5172",
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
print(f"Tau MANIFEST_SHA256.txt updated: {len(seen_hashes)} files hashed.")
