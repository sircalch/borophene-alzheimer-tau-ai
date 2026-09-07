"""
run_entire_tau_study.py
Master end-to-end pipeline for Article 4 (Alzheimer Tau / 2D beta-12 borophene).

Reproduces every real number and figure in the manuscript. Steps 1-4 curate the
cohort, compute real descriptors and run the real Tau-filament docking; step 5
is the corrected adsorption screening (relaxed complexes - see the module
docstring for why the earlier single-point version was wrong); steps 6-10 build
the leak-free QSPR, applicability domain, figures, manuscript and SI.
"""
import os
import sys
import time

BASE = os.path.dirname(os.path.abspath(__file__))


def run_step(n, total, title, rel_path, args=""):
    script = os.path.join(BASE, rel_path)
    print(f"\n{'='*70}\n  [Step {n}/{total}] {title}\n{'='*70}")
    t0 = time.time()
    ret = os.system(f'python "{script}" {args}')
    if ret != 0:
        print(f"[ERROR] Step {n}: {title} (exit {ret})")
        return False
    print(f"[OK] Step {n} in {time.time()-t0:.1f}s")
    return True


def main():
    print("=" * 70)
    print("  BOROPHENE-ALZHEIMER-TAU-AI : MASTER REPRODUCIBILITY PIPELINE")
    print("=" * 70)
    steps = [
        ("Tau drug-library curation", "src/descriptors/curate_tau_dataset.py", ""),
        ("RDKit + GFN2-xTB descriptors", "src/descriptors/compute_tau_descriptors.py", ""),
        ("Real AutoDock Vina docking (Tau filament, PDB 5O3L)", "src/docking/run_tau_real_docking.py", ""),
        ("Residue-level contact analysis", "src/docking/analyze_tau_interactions.py", ""),
        ("Adsorption screening on B40H15 (relaxed complexes)", "recompute_tau_adsorption.py", "--commit-datasets"),
        ("Isolated-drug ML table", "src/ml_models/build_tau_ml_datasets.py", ""),
        ("Leak-free nested 5x5 CV + Y-scrambling", "scripts/run_nested_cv_leakfree.py", ""),
        ("OECD applicability domain (Williams)", "src/ml_models/compute_tau_oecd_applicability_domain.py", ""),
        ("Master figure suite", "src/visualization/generate_tau_master_figures.py", ""),
        ("Word manuscript", "src/visualization/generate_tau_word_manuscript.py", ""),
        ("Supporting information", "src/visualization/generate_supporting_information.py", ""),
    ]
    for i, (title, path, args) in enumerate(steps, 1):
        if not run_step(i, len(steps), title, path, args):
            sys.exit(1)
    print("\n" + "=" * 70)
    print(">>> PIPELINE COMPLETE <<<")
    print("  manuscript/Beilstein_Manuscript_Tau_Borophene_Monreal_Hernandez_et_al.docx")
    print("=" * 70)


if __name__ == "__main__":
    main()
