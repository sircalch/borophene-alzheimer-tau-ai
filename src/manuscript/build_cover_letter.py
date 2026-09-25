"""build_cover_letter.py - cover letter for Journal of Molecular Modeling (Word)."""
import sys
from datetime import date
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import docx_kit as k  # noqa: E402
from build_manuscript import AFFIL, AUTHOR, EMAIL, ORCID, TITLE, f2, load, stats_  # noqa: E402

OUT = HERE.parents[1] / "manuscript" / "submission"


def main():
    d = load()
    s, q = stats_(d), d["q"]
    doc = k.new_document()
    for t in (AUTHOR, AFFIL, EMAIL, "", date.today().strftime("%d %B %Y"), "",
              "The Editor-in-Chief", "Journal of Molecular Modeling", ""):
        k.para(doc, t, align="left", space_after=0)
    k.para(doc, "Dear Editor,", align="left")
    k.para(doc, f"I submit the manuscript \"{TITLE}\" for consideration as an Original Paper in the Journal of "
                "Molecular Modeling.")
    k.para(doc,
           f"The study examines {s['n']} tau-directed and Alzheimer drugs at the GTP-1 site of an Alzheimer "
           "paired-helical-filament cryo-EM structure and on β_{12} borophene as a candidate carrier. Two "
           "methodological findings are reported openly: docking into the ligand-free fibril does not "
           "reproduce the stacked tracer pose, so the scores are presented as exploratory, and an "
           "unconstrained borophene flake collapses during adsorption at the GFN2-xTB level, so the carrier is "
           "modelled as a planar sheet held at the lattice of Ag-supported β_{12} borophene. On this sheet "
           f"{len(s['phys'])} drugs physisorb and {len(s['chem'])} chemisorb, with every complex checked for "
           "changes in the drug's bonding. A descriptor-based QSPR model, evaluated with nested "
           f"cross-validation and Y-scrambling, explains part of the variance (*Q*^{{2}}_{{CV}} = {f2(q['Q2_CV'])}).")
    k.para(doc,
           "All structures, relaxed geometries, docking poses and the complete pipeline that regenerates every "
           "number and figure are openly available. The manuscript is original, has not been published and is "
           "not under consideration elsewhere. The author declares no competing interests.")
    k.para(doc, "Sincerely,", align="left", space_after=0)
    k.para(doc, f"{AUTHOR} (ORCID {ORCID})", align="left")
    out = OUT / "Cover_Letter_JMM.docx"
    doc.save(out)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
