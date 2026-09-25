"""
qspr_nested_cv.py - QSPR model of the GFN2-xTB drug / beta12-borophene
interaction energy (protocol: qspr_core.py - ridge, nested 5x5 CV, 1,000
Y-permutations, leverage applicability domain).

  target      : dE_int on the B40H15 flake, drugs whose own bonding was
                unchanged by adsorption
  descriptors : fixed before any fit - MW and TPSA (RDKit, PubChem structure),
                static polarizability alpha(0) (GFN2-xTB, relaxed isolated
                drug) and formal charge. Conceptual-DFT indices were not used:
                several drugs are cations, for which omega is not comparable.

Outputs in results/qspr/: dEint_{summary.json, oof.csv, y_scrambling.csv}
"""
import json
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

sys.path.insert(0, str(Path(__file__).resolve().parent))
import qspr_core  # noqa: E402

BASE = Path(__file__).resolve().parents[2]
FEATURES = ["MW", "TPSA", "alpha", "charge"]


def main():
    lib = pd.read_csv(BASE / "data" / "processed" / "compound_library_pubchem.csv")
    iso = pd.read_csv(BASE / "results" / "quantum" / "isolated_drugs_qm_results.csv")
    ads = pd.read_csv(BASE / "results" / "quantum" / "adsorption_results.csv")
    ads = ads[ads.status == "OK"]
    df = lib[["name", "smiles"]].merge(iso[["name", "charge", "alpha_au"]], on="name").merge(
        ads[["name", "delta_Eint_kcal_mol", "drug_intact"]], on="name")
    m = df.smiles.map(Chem.MolFromSmiles)
    df["MW"], df["TPSA"] = m.map(Descriptors.MolWt), m.map(Descriptors.TPSA)
    df = df.rename(columns={"alpha_au": "alpha", "delta_Eint_kcal_mol": "dE_int"})
    excluded = df.loc[~df.drug_intact.astype(bool), "name"].tolist()
    train = df[df.drug_intact.astype(bool)].reset_index(drop=True)
    s = qspr_core.run(train, FEATURES, "dE_int", BASE / "results" / "qspr", "dEint")
    s["excluded_not_intact"] = excluded
    (BASE / "results" / "qspr" / "dEint_summary.json").write_text(json.dumps(s, indent=2))
    print(json.dumps(s, indent=2))


if __name__ == "__main__":
    main()
