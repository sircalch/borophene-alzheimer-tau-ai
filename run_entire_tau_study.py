"""
run_entire_tau_study.py
Master End-to-End Pipeline Runner for 100% Reproducibility of Article 4:
Alzheimer's Disease Tau Fibril Disaggregation & 2D Borophene Nanosheets.
"""

import os
import sys
import time

def run_step(step_num, title, script_rel_path):
    base_dir = os.path.dirname(os.path.abspath(__file__))
    script_path = os.path.join(base_dir, script_rel_path)
    print(f"\n=======================================================")
    print(f"  [Step {step_num}/8] {title}")
    print(f"=======================================================")
    t0 = time.time()
    ret = os.system(f'python "{script_path}"')
    t_elapsed = time.time() - t0
    if ret != 0:
        print(f"[ERROR] Step {step_num}: {title} (Exit Code: {ret})")
        return False
    print(f"[OK] Step {step_num} completed in {t_elapsed:.2f} seconds.")
    return True

def main():
    print("=" * 65)
    print("  BOROPHENE-ALZHEIMER-TAU-AI: MASTER REPRODUCIBILITY PIPELINE")
    print("  Authors: Andrés Monreal Hernández et al.")
    print("=" * 65)
    
    steps = [
        (1, "Tau Drug Library Curation", "src/descriptors/curate_tau_dataset.py"),
        (2, "20-Descriptor RDKit & Quantum Calculation", "src/descriptors/compute_tau_descriptors.py"),
        (3, "Parallel Real AutoDock Vina Docking (PDB 6VHL)", "src/docking/run_tau_real_docking.py"),
        (4, "Residue-Level Contact Analysis", "src/docking/analyze_tau_interactions.py"),
        (5, "Machine Learning Training & SHAP XAI", "src/ml_models/train_tau_qsar_models.py"),
        (6, "OECD Applicability Domain (Williams Plot)", "src/ml_models/compute_tau_oecd_applicability_domain.py"),
        (7, "Publication-Grade Figures Suite (300+ DPI)", "src/visualization/generate_tau_q1_figures.py"),
        (8, "Word Manuscript Compilation & Submission Packaging", "src/visualization/generate_tau_word_manuscript.py")
    ]
    
    for s_num, title, path in steps:
        success = run_step(s_num, title, path)
        if not success:
            sys.exit(1)
            
    print("\n" + "=" * 65)
    print(">>> FULL REPRODUCIBILITY PIPELINE EXECUTED SUCCESSFULLY! <<<")
    print("  Manuscript Word File: manuscript/Beilstein_Manuscript_Tau_Borophene_Monreal_Hernandez_et_al.docx")
    print("  Submission ZIP File:  borophene-alzheimer-tau-ai-FINAL-SUBMISSION-READY.zip")
    print("=" * 65)

if __name__ == "__main__":
    main()
