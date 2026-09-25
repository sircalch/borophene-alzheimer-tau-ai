"""
run_entire_tau_study.py
Reproduces every number, table and figure of the Tau / beta12-borophene article
from the raw inputs (PubChem structures, PDB 8FUG).

Steps (each resumable; the docking and quantum steps take hours on a desktop):
  1. docking at the GTP-1 site of PDB 8FUG, four redocking controls
  2. residue contacts of the docked poses
  3. planar beta12 B44H16 carrier (boron restrained, H relaxed) and lattice scan
  4. GFN2-xTB adsorption of 28 drugs on the restrained sheet
  5. reference energy of the restrained carrier (for dE_ads)
  6. QSPR: nested CV, Y-scrambling, applicability domain
  7. figures, manuscript, supporting information, cover letter
The compound library (data/processed/compound_library_pubchem.csv) was built by
src/descriptors/build_library_from_pubchem.py.
"""
import os
import subprocess
import sys
import time

BASE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable

STEPS = [
    ("AutoDock Vina docking (PDB 8FUG)", ["src/docking/run_vina_docking.py"]),
    ("Residue contacts", ["src/docking/contacts.py", "data/raw/8FUG_receptor_H.pdb", "results/docking/real_poses",
                          "data/processed/compound_library_pubchem.csv", "results/docking/residue_contacts.csv"]),
    ("Supported beta12 carrier", ["src/quantum/build_beta12_flake.py"]),
    ("beta12 lattice scan", ["src/quantum/scan_beta12_lattice.py"]),
    ("GFN2-xTB adsorption", ["src/quantum/run_adsorption.py", "--workers", "4"]),
    ("Carrier reference energy", ["src/quantum/carrier_reference.py"]),
    ("QSPR", ["src/ml_models/qspr_nested_cv.py"]),
    ("Figures", ["src/figures/make_figures.py"]),
    ("Manuscript", ["src/manuscript/build_manuscript.py"]),
    ("Supporting information", ["src/manuscript/build_si.py"]),
    ("Cover letter", ["src/manuscript/build_cover_letter.py"]),
]


def main():
    for i, (title, cmd) in enumerate(STEPS, 1):
        print(f"\n{'=' * 70}\n  [{i}/{len(STEPS)}] {title}\n{'=' * 70}", flush=True)
        t0 = time.time()
        if subprocess.run([PY, *cmd], cwd=BASE).returncode:
            sys.exit(f"step failed: {title}")
        print(f"  done in {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
