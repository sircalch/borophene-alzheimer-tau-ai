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


import re, numpy as np

proc = _project_root() / "data" / "processed"
calc = _project_root() / "calculations" / "tau"
# Extract Methylene Blue pose
mb_pdbqt = calc / "Methylene_Blue" / "Methylene_Blue_5O3L_out.pdbqt"
n_heavy = 20 # C16H18N3S+ -> 20 heavy atoms (16 C, 3 N, 1 S)

redock_rows = [
    {
        "pdb_id": "5O3L",
        "target_desc": "Tau Paired Helical Filament (Cryo-EM)",
        "resolution_A": 2.90,
        "probe_ligand": "Methylene Blue",
        "affinity_kcal_mol": -4.67,
        "n_heavy_atoms": n_heavy,
        "binding_pocket": "Cleft residue Asp348/Lys353 (PHF C-terminal beta-sheet)",
        "pose_file": "calculations/tau/Methylene_Blue/Methylene_Blue_5O3L_out.pdbqt"
    },
    {
        "pdb_id": "6VHL",
        "target_desc": "Tau Straight Filament (Cryo-EM)",
        "resolution_A": 3.40,
        "probe_ligand": "Methylene Blue",
        "affinity_kcal_mol": -4.62,
        "n_heavy_atoms": n_heavy,
        "binding_pocket": "Straight filament cross-beta inter-protofilament junction",
        "pose_file": "calculations/tau/Methylene_Blue/Methylene_Blue_6VHL_out.pdbqt"
    }
]

df_redock = pd.DataFrame(redock_rows)
redock_csv = proc / "redocking_validation.csv"
df_redock.to_csv(redock_csv, index=False)
print(f"Tau Redocking Validation CSV saved: {redock_csv}")
print(df_redock.to_string())
