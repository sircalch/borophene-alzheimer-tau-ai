"""
build_si.py - Supporting Information (Word) for the Tau / beta12-borophene study,
generated from the pipeline outputs.

writes manuscript/submission/Supporting_Information_Tau_Borophene.docx
"""
import re
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import docx_kit as k  # noqa: E402
from build_manuscript import AUTHOR, BASE, TITLE, load  # noqa: E402

OUT = BASE / "manuscript" / "submission"


def num(x, nd=1):
    return "–" if pd.isna(x) else f"{x:.{nd}f}".replace("-", "−")


def bond_types(x):
    import ast
    if not isinstance(x, str) or x in ("[]", ""):
        return "–"
    out = []
    for t in ast.literal_eval(x):
        a, b = (re.sub(r"\d", "", u) for u in t.split("-"))
        if f"{a}–{b}" not in out:
            out.append(f"{a}–{b}")
    return ", ".join(out)


def main():
    d = load()
    m = d["m"]
    doc = k.new_document()
    k.para(doc, "**Supporting Information**", align="left", size=15, space_after=4)
    k.para(doc, TITLE, align="left", size=11, space_after=4)
    k.para(doc, AUTHOR, align="left", size=10, space_after=14)

    lib = d["lib"]
    fam = m.set_index("name").family
    rows = [[r.name, fam.get(r.name, "excluded (mixture)"), str(r.pubchem_cid), r.formula, str(r.formal_charge)]
            for r in lib.itertuples()]
    k.table(doc, ("S1", "Compounds. Structures from PubChem, identity checked by InChIKey; counter-ions removed, "
                        "acids and bases neutralised where a neutral form exists."),
            ["Compound", "Family", "PubChem CID", "Formula", "Charge"], rows, align="llllc", font=7.5)

    rows = []
    for r in m.sort_values("delta_Eint_kcal_mol").itertuples():
        rows.append([r.name, num(r.vina_8FUG_kcal_mol, 2), num(r.delta_Eint_kcal_mol), num(r.delta_Eads_kcal_mol),
                     num(r.min_contact_A, 2), r.adsorption_mode[:4] + ".", bond_types(r.drug_carrier_bonds),
                     "yes" if r.drug_intact else "no"])
    k.table(doc, ("S2", "Per-drug results: Vina score at the GTP-1 site of PDB 8FUG; GFN2-xTB Δ*E*_{int} and "
                        "Δ*E*_{ads} on the supported β_{12} B_{44}H_{16} sheet (kcal mol^{−1}); closest "
                        "drug–carrier heavy-atom contact *d*_{min} (Å); regime; drug–carrier bond types; whether "
                        "the drug kept its own bonding."),
            ["Compound", "Vina", "Δ*E*_{int}", "Δ*E*_{ads}", "*d*_{min}", "Regime", "Bonds", "Intact"],
            rows, align="lccccccc", font=7.5)

    sc = pd.read_csv(BASE / "results" / "quantum" / "beta12_lattice_scan.csv")
    rows = [[f"{r.scale:.2f}", f"{r.d_BB_A:.3f}", "singlet" if r.uhf == 0 else "triplet", num(r.dE_kcal_mol)]
            for r in sc.itertuples()]
    k.table(doc, ("S3", "GFN2-xTB single-point energies of the planar B_{44}H_{16} flake with the boron lattice "
                        "scaled uniformly (hydrogens moved with their boron), relative to the lowest value."),
            ["Scale", "*d*_{BB} (Å)", "State", "Δ*E* (kcal mol^{−1})"], rows, align="cccc")

    rd = d["redock"]
    rows = [[r.control, num(r.affinity_kcal_mol, 2), num(r.rmsd_heavy_atom_A, 2), r.docking_status.split(" ")[0]]
            for r in rd.itertuples()]
    k.table(doc, ("S4", "Redocking controls for GTP-1 (Y9H) at chain K of PDB 8FUG."),
            ["Control", "Vina (kcal mol^{−1})", "RMSD (Å)", "Status"], rows, align="lccc", font=8)

    ref = d["ref_scan"].copy()
    e0 = ref.E_Eh.min()
    rows = [[r.from_complex, num((r.E_Eh - e0) * 627.509, 2)] for r in ref.itertuples()]
    k.table(doc, ("S5", "Re-relaxation of the restrained carrier taken from each physisorbed complex (drug "
                        "removed); energy relative to the lowest, which is the reference for Δ*E*_{ads}."),
            ["Carrier from complex with", "Δ*E* (kcal mol^{−1})"], rows, align="lc", font=8)

    old = d["old"]
    import json
    E0 = json.loads((BASE / "data" / "processed" / "carrier_B40H15_collapsing_2026-09-24.json").read_text())["E_Eh"]
    old = old.assign(dcar=(old.E_carrier_frozen_Eh - E0) * 627.509).sort_values("dcar")
    rows = [[r.name, num(r.dcar)] for r in old.itertuples()]
    k.table(doc, ("S6", "First carrier model (unconstrained B_{40}H_{15} flake): energy of the carrier fragment in "
                        "each adsorption complex relative to the isolated relaxed flake. Negative values mean the "
                        "flake reconstructed into a lower-energy structure during adsorption; this model was "
                        "therefore replaced by the supported sheet."),
            ["Drug", "Δ*E*_{carrier} (kcal mol^{−1})"], rows, align="lc", font=8)

    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / "Supporting_Information_Tau_Borophene.docx"
    doc.save(out)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
