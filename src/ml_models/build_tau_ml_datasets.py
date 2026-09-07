"""
build_tau_ml_datasets.py
========================
Builds the isolated-drug ML table for the Tau study by merging the real RDKit +
GFN2-xTB descriptors with the real AutoDock Vina docking scores. Nothing is
estimated from an empirical formula.

Replaces train_tau_qsar_models.py, which additionally wrote a FABRICATED
`Target_DeltaG_bind = Vina - 4.20 + 0.045*Delta_E_ads` (with Delta_E_ads itself
an empirical function of RDKit descriptors) into
dataset_drug_borophene_pristine.csv / _functionalized.csv. Those two files are no
longer used by any figure, model or the manuscript - the borophene endpoint is
the real GFN2-xTB `delta_Eint_SP_kcal_mol` in dataset_tau_borophene_pristine.csv
(see recompute_tau_adsorption.py). This script deletes the stale fabricated files
if present.
"""
import os
import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PROC = os.path.join(BASE, "data", "processed")


def main():
    # Use the 5O3L docking run (vina_5O3L_kcal_mol) that the rest of the paper
    # uses - results/docking/real_vina_docking_summary.csv is a stale run on a
    # different (BACE / gamma-secretase) cohort, mean -4.36.
    boro = pd.read_csv(os.path.join(PROC, "dataset_tau_borophene_pristine.csv"))
    desc = pd.read_csv(os.path.join(PROC, "tau_isolated_descriptors.csv"))
    vcol = "Real_Vina_Docking_Score_kcal_mol"
    v = boro[["name", "vina_5O3L_kcal_mol"]].rename(columns={"vina_5O3L_kcal_mol": vcol})
    merged = pd.merge(desc, v, on="name")
    merged.to_csv(os.path.join(PROC, "dataset_isolated_tau_drugs.csv"), index=False)
    print(f"[OK] dataset_isolated_tau_drugs.csv  n={len(merged)}  "
          f"Vina(5O3L) {merged[vcol].min():.2f} to {merged[vcol].max():.2f} (mean {merged[vcol].mean():.2f})")

    for stale in ("dataset_drug_borophene_pristine.csv", "dataset_drug_borophene_functionalized.csv"):
        p = os.path.join(PROC, stale)
        if os.path.exists(p):
            os.remove(p)
            print(f"[removed stale fabricated file] {stale}")


if __name__ == "__main__":
    main()
