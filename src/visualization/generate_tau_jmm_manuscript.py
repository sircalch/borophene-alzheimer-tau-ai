"""
generate_tau_jmm_manuscript.py
================================
Builds a Journal of Molecular Modeling (Springer) submission variant of the
Tau/borophene manuscript. Does NOT touch the canonical Beilstein/Molecular
Diversity body produced by generate_tau_word_manuscript.py -- it
post-processes a freshly regenerated copy of that docx, same recipe as the
GBM JMM variant (generate_gbm_jmm_manuscript.py):

  1. Structured Context/Methods abstract (JMM requirement, 150-250 words).
     Context is a condensed re-derivation of the original abstract with the
     dynamically-computed docking numbers (Vina range/mean, Q2_CV) extracted
     via regex from the freshly generated text, never retyped by hand.
  2. Section reorder: JMM wants Methods to follow the Introduction. The
     Beilstein body instead puts "4. Experimental" after Conclusions
     (Beilstein house style). Moved + renumbered:
         1. Introduction            (unchanged)
         2. Experimental            (was "4.", subsections 4.1-4.3 -> 2.1-2.3)
         3. Results and Discussion  (was "2.", subsections 2.1-2.5 -> 3.1-3.5)
         4. Summary                 (was "3. Conclusions")
     No in-text cross-references to fix here (unlike GBM's "Section 2.1").
"""

import os
import re
import sys

base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from docx import Document

import generate_tau_word_manuscript as full_gen

JMM_METHODS = (
    "All calculations use GFN2-xTB (xtb v6.7.1, D4 dispersion). Drugs and the hydrogen-terminated β12 borophene "
    "cluster (B40H15, C1) were geometry-optimized; each complex was built by centring the drug 3.2 Å above the "
    "sheet in four orientations, retaining the lowest-energy pose. Contacts below 1.9 Å are classified as "
    "chemisorption, the rest as physisorption. AutoDock Vina v1.2.7 docked against the cryo-EM Alzheimer Tau "
    "filament core (PDB ID: 5O3L); as a cross-β assembly, scores are a relative ranking. A StandardScaler+RidgeCV "
    "surrogate was trained inside a leak-free nested 5×5 cross-validation, feature importance via SHAP, "
    "applicability domain via OECD Principle 3 (Williams leverage)."
)


def _build_condensed_context(original_abstract):
    """Re-derive a <=140-word Context paragraph, extracting the dynamically
    computed docking numbers (Vina range/mean, Q2_CV) with regex so they can
    never drift stale, and dropping sentences that duplicate the Methods
    paragraph or are out of scope (BBB-transcytosis future-work aside).
    """
    m_vina = re.search(r"Vina scores of ([\-0-9.]+ to [\-0-9.]+) kcal/mol \(mean ([\-0-9.]+)\)", original_abstract)
    m_q2 = re.search(r"Q2_CV = ([0-9.n/a]+) for the Tau-filament Vina docking score", original_abstract)
    if not (m_vina and m_q2):
        raise RuntimeError("Could not extract one or more dynamic values from the original Tau abstract "
                            "-- source text likely changed; update the regexes in _build_condensed_context.")
    vina_range, vina_mean = m_vina.group(1), m_vina.group(2)
    q2_vina = m_q2.group(1)

    return (
        "Tau aggregation into paired helical filaments is a defining neuropathological hallmark of Alzheimer's "
        "disease. We screen 29 Tau-directed therapeutics against pristine β12 borophene (B40H15) with GFN2-xTB, "
        "relaxing every drug-surface complex, and dock them against the cryo-EM Alzheimer Tau filament core "
        "(PDB ID: 5O3L). The complexes fall into two separated regimes: 12 of 29 ligands chemisorb, forming a "
        "covalent B-C or B-O bond (1.36-1.69 Å, ΔE_int,SP = -81 to -233 kcal/mol), while the remaining 17 "
        "physisorb (2.6-3.7 Å, -8 to -44 kcal/mol) -- pristine borophene is a chemically reactive surface, not "
        f"a reversible carrier, for most of these ligands. Docking gave Vina scores of {vina_range} kcal/mol "
        f"(mean {vina_mean}); a leak-free nested surrogate reached Q2_CV = {q2_vina} for the docking score "
        "(Y-scrambled Q2≈0.15, p=0.009) and 0.06 for the physisorption energy -- the chemisorption/physisorption "
        "dichotomy, not descriptor-based prediction, is the robust result."
    )


def _find_heading(doc, text_exact=None):
    for i, p in enumerate(doc.paragraphs):
        if not p.style.name.startswith("Heading"):
            continue
        if text_exact is not None and p.text.strip() == text_exact:
            return i, p
    raise RuntimeError(f"Heading not found: {text_exact}")


def generate_tau_jmm_manuscript():
    full_gen.generate_tau_word_manuscript()

    src_docx = os.path.join(base_dir, "manuscript", "Beilstein_Manuscript_Tau_Borophene_Monreal_Hernandez_et_al.docx")
    doc = Document(src_docx)

    # ---------------------------------------------------------------
    # 1) Structured Context/Methods abstract
    # ---------------------------------------------------------------
    abs_head_idx, _ = _find_heading(doc, text_exact="Abstract")
    abstract_para = doc.paragraphs[abs_head_idx + 1]
    keywords_para = doc.paragraphs[abs_head_idx + 2]
    if "Keywords" not in keywords_para.text:
        raise RuntimeError("Unexpected structure: paragraph after Abstract body is not Keywords.")

    original_abstract = abstract_para.text
    context_text = _build_condensed_context(original_abstract)

    p_context = abstract_para.insert_paragraph_before()
    p_context.paragraph_format.space_after = abstract_para.paragraph_format.space_after
    r1 = p_context.add_run("Context ")
    r1.font.bold = True
    p_context.add_run(context_text)

    p_methods = abstract_para.insert_paragraph_before()
    p_methods.paragraph_format.space_after = abstract_para.paragraph_format.space_after
    r2 = p_methods.add_run("Methods ")
    r2.font.bold = True
    p_methods.add_run(JMM_METHODS)

    abstract_para._element.getparent().remove(abstract_para._element)

    # ---------------------------------------------------------------
    # 2) Move "4. Experimental" (+ its 3 paragraphs) to right after Introduction
    # ---------------------------------------------------------------
    exp_idx, _ = _find_heading(doc, text_exact="4. Experimental")
    data_avail_idx, _ = _find_heading(doc, text_exact="Data Availability")
    results_idx, results_head = _find_heading(doc, text_exact="2. Results and Discussion")

    elements_to_move = [doc.paragraphs[i]._p for i in range(exp_idx, data_avail_idx)]
    anchor_element = results_head._p
    for el in elements_to_move:
        anchor_element.addprevious(el)

    # ---------------------------------------------------------------
    # 3) Renumber headings and the 4.x -> 2.x inline paragraph labels
    # ---------------------------------------------------------------
    _, exp_head = _find_heading(doc, text_exact="4. Experimental")
    exp_head.runs[0].text = "2. Experimental"

    _, results_head2 = _find_heading(doc, text_exact="2. Results and Discussion")
    results_head2.runs[0].text = "3. Results and Discussion"

    subsection_renumber = {
        "2.1 Two adsorption regimes on pristine beta-12 borophene": "3.1 Two adsorption regimes on pristine beta-12 borophene",
        "2.2 Docking against the cryo-EM Tau filament core": "3.2 Docking against the cryo-EM Tau filament core",
        "2.3 Nano-QSAR surrogate model": "3.3 Nano-QSAR surrogate model",
        "2.4 Applicability domain (OECD Principle 3)": "3.4 Applicability domain (OECD Principle 3)",
        "2.5 Interfacial charge redistribution": "3.5 Interfacial charge redistribution",
    }
    for old, new in subsection_renumber.items():
        _, h = _find_heading(doc, text_exact=old)
        h.runs[0].text = new

    _, concl_head = _find_heading(doc, text_exact="3. Conclusions")
    concl_head.runs[0].text = "4. Summary"

    label_renumber = {
        "4.1 Quantum-chemical framework": "2.1 Quantum-chemical framework",
        "4.2 Molecular docking": "2.2 Molecular docking",
        "4.3 Surrogate model and applicability domain": "2.3 Surrogate model and applicability domain",
    }
    for p in doc.paragraphs:
        for old, new in label_renumber.items():
            if p.text.strip().startswith(old) and p.runs:
                p.runs[0].text = p.runs[0].text.replace(old, new, 1)

    out_dir = os.path.join(base_dir, "manuscript", "submission_ready")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "02_Manuscript_Tau_Borophene_JMM_Submission.docx")
    doc.save(out_path)
    print(f"[SUCCESS] Generated JMM submission manuscript: {out_path}")
    return out_path


if __name__ == "__main__":
    generate_tau_jmm_manuscript()
