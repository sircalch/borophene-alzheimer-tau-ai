"""
generate_tau_word_manuscript.py
Builds the complete, publication-grade Microsoft Word (.docx) manuscript
with all 9 figures embedded, formatted tables, and 45 verified citations for Article 4 (Tau & Borophene).
"""

import os
import json
import pandas as pd
from docx import Document
from docx.shared import Inches, Pt, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml import OxmlElement, parse_xml
from docx.oxml.ns import nsdecls, qn

def set_cell_background(cell, fill_color):
    tcPr = cell._element.get_or_add_tcPr()
    shd = parse_xml(f'<w:shd {nsdecls("w")} w:fill="{fill_color}"/>')
    tcPr.append(shd)

def set_cell_margins(cell, top=100, bottom=100, left=150, right=150):
    tcPr = cell._element.get_or_add_tcPr()
    tcMar = OxmlElement('w:tcMar')
    for m, val in [('top', top), ('bottom', bottom), ('left', left), ('right', right)]:
        node = OxmlElement(f'w:{m}')
        node.set(qn('w:w'), str(val))
        node.set(qn('w:type'), 'dxa')
        tcMar.append(node)
    tcPr.append(tcMar)

def add_heading_styled(doc, text, level=1):
    h = doc.add_heading(text, level=level)
    h.paragraph_format.space_before = Pt(14)
    h.paragraph_format.space_after = Pt(6)
    h.paragraph_format.keep_with_next = True
    for r in h.runs:
        r.font.name = 'Times New Roman'
        r.font.bold = True
        if level == 1:
            r.font.size = Pt(14)
            r.font.color.rgb = RGBColor(74, 20, 140)
        elif level == 2:
            r.font.size = Pt(12)
            r.font.color.rgb = RGBColor(49, 27, 146)
        else:
            r.font.size = Pt(11)
            r.font.color.rgb = RGBColor(33, 33, 33)
    return h

def add_image_if_exists(doc, img_path, caption_text, width=Inches(6.2)):
    if os.path.exists(img_path):
        p_img = doc.add_paragraph()
        p_img.alignment = WD_ALIGN_PARAGRAPH.CENTER
        p_img.paragraph_format.space_before = Pt(10)
        p_img.paragraph_format.space_after = Pt(4)
        run = p_img.add_run()
        run.add_picture(img_path, width=width)
        
        p_cap = doc.add_paragraph()
        p_cap.paragraph_format.space_after = Pt(12)
        p_cap.paragraph_format.line_spacing = 1.15
        r_num = p_cap.add_run(caption_text.split(':')[0] + ": ")
        r_num.font.bold = True
        r_num.font.size = Pt(9.5)
        r_num.font.color.rgb = RGBColor(74, 20, 140)
        
        r_desc = p_cap.add_run(':'.join(caption_text.split(':')[1:]))
        r_desc.font.size = Pt(9.5)
        r_desc.font.italic = True
    else:
        print(f"Warning: image {img_path} not found.")

def generate_tau_word_manuscript():
    base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    fig_dir = os.path.join(base_dir, "figures")
    doc = Document()
    
    for s in doc.sections:
        s.top_margin = Inches(1.0)
        s.bottom_margin = Inches(1.0)
        s.left_margin = Inches(1.0)
        s.right_margin = Inches(1.0)
        
    style = doc.styles['Normal']
    font = style.font
    font.name = 'Times New Roman'
    font.size = Pt(11)
    font.color.rgb = RGBColor(33, 33, 33)
    
    # Title & Authors
    p_title = doc.add_paragraph()
    p_title.paragraph_format.space_after = Pt(12)
    p_title.paragraph_format.line_spacing = 1.15
    r_title = p_title.add_run("Machine Learning-Driven Nano-QSAR and Quantum Chemical Design of Functionalized 2D Borophene Nanosheets for Targeted Disaggregation of Pathological Tau Fibrils in Alzheimer's Disease")
    r_title.font.size = Pt(16)
    r_title.font.bold = True
    r_title.font.color.rgb = RGBColor(74, 20, 140)
    
    p_auth = doc.add_paragraph()
    p_auth.paragraph_format.space_after = Pt(4)
    r_a1 = p_auth.add_run("Andrés Monreal Hernández")
    r_a1.font.bold = True
    p_auth.add_run("1,*, ")
    r_a2 = p_auth.add_run("Sara Lizbeth Franco Amaya")
    r_a2.font.bold = True
    p_auth.add_run("2, and ")
    r_a3 = p_auth.add_run("Carlos Ivanhoe Martínez Osorio")
    r_a3.font.bold = True
    p_auth.add_run("3")
    
    p_aff = doc.add_paragraph()
    p_aff.paragraph_format.space_after = Pt(14)
    p_aff.add_run(
        "1 Universidad Estatal de Sonora, Hermosillo, Sonora, Mexico. ORCID: 0009-0009-1207-8597\n"
        "2 Doctorado en Nanotecnología, Universidad de Sonora, Hermosillo, Sonora, Mexico. ORCID: 0009-0005-0272-0241\n"
        "3 Doctorado en Ciencia de Materiales, Universidad de Sonora, Hermosillo, Sonora, Mexico. ORCID: 0009-0003-7872-4965\n"
        "* Corresponding author: andres.monreal@ues.mx"
    )
    p_aff.runs[0].font.size = Pt(9.5)
    p_aff.runs[0].font.italic = True
    
    # Graphical Abstract
    add_image_if_exists(doc, os.path.join(fig_dir, "fig1_graphical_abstract.png"),
                        "Graphical Abstract: Multi-Scale Quantum, Docking, and Machine Learning Framework for 2D Borophene Transcytosis and Targeted Disaggregation of Alzheimer's Tau Fibrils.")
    
    # Abstract
    add_heading_styled(doc, "Abstract", level=1)
    p_abs = doc.add_paragraph()
    p_abs.paragraph_format.space_after = Pt(8)
    p_abs.paragraph_format.line_spacing = 1.15
    p_abs.add_run(
        "Hyperphosphorylation and aggregation of microtubule-associated protein Tau into paired helical filaments (PHFs) is a defining neuropathological "
        "hallmark of Alzheimer's disease (AD) and correlates closely with cognitive decline [1,3,6]. Here we present a computational framework combining "
        "GFN2-xTB tight-binding quantum chemistry (with D4 dispersion) [24,26], physical molecular docking (AutoDock Vina v1.2.7 [31,32] against the "
        "cryo-EM Alzheimer Tau filament core, PDB ID: 5O3L [7]), and a leak-free cross-validated explainable Nano-QSAR surrogate, for a curated cohort "
        "of 29 clinical-stage and experimental Tau-directed therapeutics (including hydromethylthionine/LMTX, EGCG, curcumin, tideglusib and AZD1080). "
        "Real GFN2-xTB single-point interaction energies of the 26 non-anomalous compounds on the pristine beta-12 borophene cluster (B40H15) span "
        "-13.8 to -0.9 kcal/mol; three cationic phenothiazine dyes (methylene blue, azure A, toluidine blue O) give large positive single-point energies "
        "from steric clash and are reported as such and flagged out-of-domain. An Angiopep-2 or chi3-PEG-Tf functionalized borophene for BBB transcytosis "
        "is discussed only as future work, since no real structural or quantum data for it exist in this study. Docking against the Tau filament core gave "
        "Vina scores of -3.79 to -6.83 kcal/mol (mean -5.17), with recurrent contacts at the cross-beta residues Gly335, Leu357, Gln336, Val337 and Pro332. "
        "A leak-free nested 5x5 cross-validated RidgeCV surrogate on the real data reached Q2_CV = 0.46 (isolated descriptors) and 0.07 (pristine-borophene "
        "interaction energy); the feature-importance analysis is reported as exploratory. Every value is computed from the deposited pipeline; no descriptor "
        "or energy is estimated from an empirical formula."
    )
    
    p_kw = doc.add_paragraph()
    p_kw.paragraph_format.space_after = Pt(14)
    r_kwt = p_kw.add_run("Keywords: ")
    r_kwt.font.bold = True
    p_kw.add_run("2D Borophene; Alzheimer's Disease; Tau Paired Helical Filaments; LMTX; Blood-Brain Barrier; AutoDock Vina; Explainable AI (SHAP); OECD Validation.")
    
    # Sections
    add_heading_styled(doc, "1. Introduction", level=1)
    doc.add_paragraph(
        "Alzheimer's disease (AD) is the most common neurodegenerative disorder worldwide [1]. Neurofibrillary tangles built from hyperphosphorylated Tau "
        "filaments follow a stereotyped spatiotemporal progression (Braak staging) that tracks clinical severity more tightly than amyloid burden [2,3,6]. "
        "Cryo-electron microscopy of patient-derived filament cores [7-10] has provided an atomic template for structure-based design of Tau ligands and "
        "aggregation modulators, of which the phenothiazine leuco-methylthioninium (hydromethylthionine / LMTX) is the most clinically advanced [11,12]; "
        "natural polyphenols such as EGCG and curcumin also remodel or inhibit Tau assembly in vitro [14,15]."
    )
    doc.add_paragraph(
        "Two-dimensional boron (borophene), synthesized as several polymorphs including the beta-12 lattice [21-23], is metallic, strongly polarizable and "
        "forms delocalized multicentre B-B bonds, motivating its evaluation as a drug-loading surface. In this study the carrier is a hydrogen-terminated "
        "beta-12 borophene cluster (B40H15); a peptide-functionalized allotrope for receptor-mediated BBB transcytosis is considered only as a prospective "
        "extension (Conclusions), and all quantum results refer to the pristine surface."
    )

    add_image_if_exists(doc, os.path.join(fig_dir, "fig1_tau_workflow_methodology.png"),
                        "Figure 1: Multi-scale computational workflow: GFN2-xTB quantum-chemical adsorption on pristine beta-12 borophene (B40H15), real AutoDock Vina docking against the cryo-EM Tau filament core (PDB 5O3L), and a leak-free cross-validated explainable Nano-QSAR surrogate.")

    add_heading_styled(doc, "2. Computational and Experimental Section", level=1)
    doc.add_paragraph(
        "2.1 Quantum-chemical framework: Each isolated drug, the pristine beta-12 borophene cluster (B40H15), and every drug-borophene complex were "
        "geometry-optimized and evaluated at single point with GFN2-xTB (xtb v6.7.1) including the D4 dispersion correction [24,26]. The standardized "
        "single-point interaction energy is Delta_E_int,SP = E(complex) - E(borophene) - E(drug), both fragments taken at the complex geometry. "
        "Frontier-orbital energies and conceptual-DFT global reactivity indices (hardness eta = gap/2, softness, electronegativity, electrophilicity "
        "omega = mu^2/2eta) [39,40] were read directly from the xtb output; no descriptor is estimated from an empirical formula."
    )
    doc.add_paragraph(
        "2.2 Molecular docking: Docking used AutoDock Vina v1.2.7 [31,32] against the cryo-EM structure of the Alzheimer Tau filament core (PDB ID: 5O3L [7]), "
        "with ligands protonated at pH 7.4 and prepared with Meeko. Because the paired-helical-filament core is a cross-beta assembly rather than a globular "
        "pocket, the Vina scores are reported as a relative ranking of surface / cleft affinity rather than an absolute binding free energy."
    )
    doc.add_paragraph(
        "2.3 Surrogate model and applicability domain: A StandardScaler + RidgeCV model was trained inside a leak-free nested 5x5 cross-validation on the "
        "real observed data (four descriptors: MolWt, MolMR, E_HOMO, omega). Feature importance was inspected with an ExtraTrees estimator and SHAP [38] "
        "and is reported as exploratory only. The applicability domain follows OECD Principle 3 [34-36] via Williams hat-matrix leverage."
    )
    
    add_image_if_exists(doc, os.path.join(fig_dir, "fig2_tau_quantum_cdft_architecture.png"),
                        "Figure 2: Real quantum conceptual-DFT electronic reactivity of the isolated Tau therapeutics (real GFN2-xTB single points, n=29): (a) HOMO/LUMO frontier-orbital distribution; (b) chemical hardness vs. electrophilicity index. No real complex-level frontier-orbital calculation exists for the borophene surface.")
    
    add_heading_styled(doc, "3. Results and Discussion", level=1)

    add_heading_styled(doc, "3.1 Quantum interaction energies on pristine beta-12 borophene", level=2)
    doc.add_paragraph(
        "Real GFN2-xTB single-point interaction energies for the 26 non-anomalous therapeutics on the B40H15 cluster range from -0.9 to -13.8 kcal/mol "
        "(curcumin -7.8; EGCG -5.3 kcal/mol), consistent with dispersion-assisted physisorption of the drug pi-systems on the polarizable metallic boron "
        "lattice. Three cationic phenothiazine dyes (methylene blue, azure A, toluidine blue O) return large positive single-point energies (up to "
        "+190 kcal/mol): at the fixed complex geometry the rigid planar cation is forced into steric overlap with the lattice, so these points are "
        "physically meaningless and are excluded from the model and flagged out-of-domain rather than removed silently."
    )

    add_heading_styled(doc, "3.2 Docking against the cryo-EM Tau filament core", level=2)
    doc.add_paragraph(
        "AutoDock Vina scores against the Tau filament core (PDB 5O3L) span -3.79 to -6.83 kcal/mol (mean -5.17). The highest-ranked ligands are "
        "chrysamine G (-6.83), EGCG (-5.88), luteolin (-5.87), fisetin (-5.75) and donepezil (-5.75 kcal/mol); hydromethylthionine/LMTX scores -4.77 kcal/mol. "
        "The cross-beta assembly offers no deep pocket, so these values rank relative surface / cleft affinity rather than absolute binding free energy."
    )

    add_image_if_exists(doc, os.path.join(fig_dir, "fig3_tau_docking_vina_statistical_profiles.png"),
                        "Figure 3: Molecular docking statistical profiles against the cryo-EM Tau filament core (PDB 5O3L): (a) distribution of real Vina scores; (b) ranking of the top-10 compounds (chrysamine G -6.83, EGCG -5.88 kcal/mol; hydromethylthionine/LMTX -4.77 kcal/mol).")

    add_image_if_exists(doc, os.path.join(fig_dir, "fig4_tau_residue_contact_frequency.png"),
                        "Figure 4: Residue-level contact frequencies on the Tau filament core (real Vina poses, contact distance <= 3.8 A): most frequent contacts are the cross-beta residues Gly335, Leu357, Gln336, Val337 and Pro332.")
    
    # Table 1: Descriptors. MW/LogP/PSA are real RDKit descriptors (always
    # computed from SMILES). E_HOMO/omega previously came from
    # tau_isolated_descriptors.csv, whose E_HOMO was an empirical-formula
    # placeholder ("-5.20 - 0.18*LogP - ...") never overwritten with real
    # data; merged here with real GFN2-xTB frontier orbitals parsed from
    # calculations/tau/*/*_drug_sp.out.
    desc_csv = os.path.join(base_dir, "data", "processed", "tau_isolated_descriptors.csv")
    homo_lumo_csv = os.path.join(base_dir, "data", "processed", "tau_isolated_real_homo_lumo.csv")
    if os.path.exists(desc_csv):
        df_desc = pd.read_csv(desc_csv)
        if os.path.exists(homo_lumo_csv):
            df_real = pd.read_csv(homo_lumo_csv)
            df_desc = df_desc.merge(df_real, on="name", how="inner")
            df_desc["E_HOMO"] = df_desc["E_HOMO_real_eV"]
            gap = df_desc["E_LUMO_real_eV"] - df_desc["E_HOMO_real_eV"]
            mu = -(df_desc["E_HOMO_real_eV"] + df_desc["E_LUMO_real_eV"]) / 2.0
            df_desc["Electrophilicity_omega"] = mu ** 2 / (2.0 * (gap / 2.0))
        doc.add_paragraph()
        p_t1 = doc.add_paragraph()
        r_t1 = p_t1.add_run("Table 1: Physicochemical, Topological, and Quantum CDFT Descriptors for Representative Alzheimer/Tau Therapeutics.")
        r_t1.font.bold = True
        r_t1.font.size = Pt(10)
        
        table1 = doc.add_table(rows=1, cols=7)
        table1.alignment = WD_TABLE_ALIGNMENT.CENTER
        hdr_cells = table1.rows[0].cells
        hdr_titles = ["Compound", "Class", "MW (g/mol)", "LogP", "PSA (Å²)", "E_HOMO (eV)", "omega (eV)"]
        for idx, title in enumerate(hdr_titles):
            hdr_cells[idx].text = title
            set_cell_background(hdr_cells[idx], "4A148C")
            set_cell_margins(hdr_cells[idx], 80, 80, 100, 100)
            for r in hdr_cells[idx].paragraphs[0].runs:
                r.font.bold = True
                r.font.color.rgb = RGBColor(255, 255, 255)
                r.font.size = Pt(9)
                
        for _, row in df_desc.head(10).iterrows():
            row_cells = table1.add_row().cells
            row_vals = [
                str(row['name']), str(row['drug_class'])[:22], f"{row['MW']:.1f}",
                f"{row['LogP']:.2f}", f"{row['PSA']:.1f}", f"{row['E_HOMO']:.2f}", f"{row['Electrophilicity_omega']:.2f}"
            ]
            for c_idx, val in enumerate(row_vals):
                row_cells[c_idx].text = val
                set_cell_margins(row_cells[c_idx], 60, 60, 80, 80)
                for r in row_cells[c_idx].paragraphs[0].runs:
                    r.font.size = Pt(8.5)
                    
    add_heading_styled(doc, "3.3 Nano-QSAR surrogate model", level=2)
    doc.add_paragraph(
        "A StandardScaler + RidgeCV surrogate evaluated by leak-free nested 5x5 cross-validation reached Q2_CV = 0.46 for the isolated-descriptor model and "
        "0.07 for the pristine-borophene interaction-energy model (n = 29, four descriptors: MolWt, MolMR, E_HOMO, omega). The borophene model is therefore "
        "essentially non-predictive; the exploratory ExtraTrees / SHAP ranking (Figure 6), led by E_HOMO, is reported only as a qualitative indication and "
        "not as a validated structure-property relationship [42,43]."
    )

    add_image_if_exists(doc, os.path.join(fig_dir, "fig5_tau_parity_models_evaluation.png"),
                        "Figure 5: Leak-free nested 5x5 CV parity plots (real observed vs out-of-fold predicted) for the isolated and pristine-borophene systems. The chi3-PEG-Tf functionalized system has no real data and is not shown.")

    add_image_if_exists(doc, os.path.join(fig_dir, "fig6_tau_shap_xai_importance_rankings.png"),
                        "Figure 6: Exploratory feature-importance ranking on the real GFN2-xTB pristine-borophene interaction energy.")

    add_image_if_exists(doc, os.path.join(fig_dir, "fig7_tau_descriptor_correlation_matrix.png"),
                        "Figure 7: Pearson inter-descriptor correlation heatmap (real descriptor matrix, 29 Alzheimer/Tau therapeutics).")

    add_heading_styled(doc, "3.4 Applicability domain (OECD Principle 3)", level=2)
    doc.add_paragraph(
        "Williams hat-matrix leverage on the real 8-descriptor matrix gives a warning leverage h* = 0.93; 28 of the 29 compounds fall inside the domain "
        "(leverage below h* and standardized residual within +/-3sigma) for both real-data systems [34-36]."
    )

    add_image_if_exists(doc, os.path.join(fig_dir, "fig8_tau_williams_applicability_domain.png"),
                        "Figure 8: OECD Principle 3 Williams plots defining the applicability domain for the Tau therapeutics on beta-12 borophene (real data only).")

    add_image_if_exists(doc, os.path.join(fig_dir, "fig9_tau_3d_spatial_binding_modes.png"),
                        "Figure 9: Representative binding modes (schematic): (a) EGCG at the Tau filament cleft (PDB 5O3L); (b) hydromethylthionine/LMTX pose; (c) a drug on the pristine beta-12 borophene surface with its real GFN2-xTB Delta_E_int,SP.")

    add_heading_styled(doc, "4. Conclusions", level=1)
    doc.add_paragraph(
        "We report a quantum-informed, explainable Nano-QSAR analysis of pristine beta-12 borophene (B40H15) as a candidate loading surface for "
        "Tau-directed therapeutics. Real GFN2-xTB single-point interaction energies show dispersion-assisted physisorption of the 26 non-anomalous drugs "
        "(Delta_E_int,SP = -0.9 to -13.8 kcal/mol), while three rigid cationic phenothiazine dyes are sterically incompatible at the fixed geometry and are "
        "flagged rather than hidden. Docking against the cryo-EM Tau filament core ranks polyphenols (chrysamine G, EGCG, luteolin, fisetin) highest. "
        "The surrogate model is not predictive for the borophene interaction energy, so the descriptor rankings are exploratory. A peptide-functionalized "
        "borophene for LRP-1-mediated BBB transcytosis is a natural extension but has no real structural or quantum data here and would require dedicated "
        "complex-geometry modeling."
    )
    
    add_heading_styled(doc, "Acknowledgements & Data Availability", level=1)
    doc.add_paragraph("Supported by Universidad Estatal de Sonora and Universidad de Sonora. Full code and docking PDBQT files are available in the repository.")
    
    add_heading_styled(doc, "References", level=1)
    from build_tau_verified_references import TAU_VERIFIED_REFERENCES as VERIFIED_REFERENCES
    for idx, ref in enumerate(VERIFIED_REFERENCES, 1):
        p_ref = doc.add_paragraph()
        p_ref.paragraph_format.left_indent = Inches(0.4)
        p_ref.paragraph_format.space_after = Pt(3)
        r_num = p_ref.add_run(f"{idx}. ")
        r_num.font.bold = True
        p_ref.add_run(ref['citation'] + " ")
        if ref.get('doi'):
            r_doi = p_ref.add_run(f"doi:{ref['doi']}")
            r_doi.font.italic = True
            r_doi.font.size = Pt(9.0)
            r_doi.font.color.rgb = RGBColor(74, 20, 140)
        
    out_docx = os.path.join(base_dir, "manuscript", "Beilstein_Manuscript_Tau_Borophene_Monreal_Hernandez_et_al.docx")
    doc.save(out_docx)
    print(f"Generated Comprehensive Tau Word Manuscript: {out_docx}")
    return out_docx

if __name__ == "__main__":
    generate_tau_word_manuscript()
