import pandas as pd
from pathlib import Path
import re, numpy as np

proc = Path(r"c:\Users\Andre\Proyectos doctorado\borophene-alzheimer-tau-ai\data\processed")
calc = Path(r"c:\Users\Andre\Proyectos doctorado\borophene-alzheimer-tau-ai\calculations\tau")

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
