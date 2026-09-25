"""
build_manuscript.py - Journal of Molecular Modeling submission (Word) for the
Tau / beta12-borophene study. Every number in the text, tables and captions is
read from the result files written by the pipeline.

usage: python src/manuscript/build_manuscript.py
writes manuscript/submission/Manuscript_Tau_Borophene_JMM.docx
"""
import json
import re
import sys
from pathlib import Path

import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "figures"))
import docx_kit as k  # noqa: E402
from references import REFS  # noqa: E402

BASE = HERE.parents[1]
FIG = BASE / "figures"
OUT = BASE / "manuscript" / "submission"

AUTHOR = "Andrés Monreal Hernández"
AFFIL = "Universidad Estatal de Sonora, Ley Federal del Trabajo S/N, Col. Apolo, 83100 Hermosillo, Sonora, Mexico"
EMAIL = "andres.monreal@ues.mx"
ORCID = "0009-0009-1207-8597"
REPO = "https://github.com/sircalch/borophene-alzheimer-tau-ai"
TITLE = ("Tau-directed drugs on β_{12} borophene: docking at the GTP-1 site of Alzheimer paired helical "
         "filaments and GFN2-xTB adsorption on a substrate-supported sheet")


class Cites:
    def __init__(self):
        self.order = []

    def __call__(self, *keys):
        nums = sorted({self._n(k_) for k_ in keys})
        spans, start = [], nums[0]
        for a, b in zip(nums, nums[1:] + [None]):
            if b != a + 1:
                spans.append(f"{start}" if start == a else f"{start}–{a}" if a - start > 1 else f"{start}, {a}")
                start = b
        return "[" + ", ".join(spans) + "]"

    def _n(self, key):
        if key not in REFS:
            raise KeyError(key)
        if key not in self.order:
            self.order.append(key)
        return self.order.index(key) + 1

    def list(self):
        return [REFS[k_] for k_ in self.order]


def f1(x):
    return f"{x:.1f}".replace("-", "−")


def f2(x):
    return "0.00" if abs(x) < 0.005 else f"{x:.2f}".replace("-", "−")


def sci(p):
    if p >= 0.001:
        return f"{p:.3f}"
    m_, e = f"{p:.1e}".split("e")
    return f"{m_} × 10^{{−{int(e[1:])}}}"


def dn(name):
    return name[0].lower() + name[1:] if name[0].isupper() and name[1:2].islower() and " " not in name else name


def load():
    import make_figures as F
    d = F.data()
    d["carrier"] = json.loads((BASE / "data" / "processed" / "carrier.json").read_text())
    d["ref_scan"] = pd.read_csv(BASE / "calculations" / "tau" / "carrier_reference" / "carrier_reference_scan.csv")
    d["old"] = pd.read_csv(BASE / "results" / "quantum" / "adsorption_results_B40H15_collapsing_2026-09-24.csv")
    d["lib"] = pd.read_csv(BASE / "data" / "processed" / "compound_library_pubchem.csv")
    return d


def stats_(d):
    m = d["m"]
    rd = d["redock"].set_index("control")
    s = {"n": len(m), "n_lib": len(d["lib"])}
    s["r_free_x"] = rd.loc["ligand-free fibril: self-redock, crystal conformation", "rmsd_heavy_atom_A"]
    s["r_free_s"] = rd.loc["ligand-free fibril: production protocol, from SMILES", "rmsd_heavy_atom_A"]
    s["r_stack_x"] = rd.loc["tracer stack present: self-redock, crystal conformation", "rmsd_heavy_atom_A"]
    s["r_stack_s"] = rd.loc["tracer stack present: production protocol, from SMILES", "rmsd_heavy_atom_A"]
    ct = d["contacts"]
    n = ct.name.nunique()
    s["freq"] = (ct.groupby("residue").name.nunique() / n * 100).sort_values(ascending=False)
    s["polar"] = (ct[ct.polar].groupby("residue").name.nunique() / n * 100).sort_values(ascending=False)
    s["chains"] = sorted(ct.chain.unique())
    v = m.set_index("name").vina_8FUG_kcal_mol.sort_values()
    s["vina"] = v
    s["fam_v"] = m.groupby("family").vina_8FUG_kcal_mol.median()
    s["fam_e"] = m.groupby("family").delta_Eint_kcal_mol.median()
    s["kw_e"] = stats.kruskal(*[g.delta_Eint_kcal_mol for _, g in m.groupby("family")])
    s["kw_v"] = stats.kruskal(*[g.vina_8FUG_kcal_mol for _, g in m.groupby("family")])
    s["chem"] = m[m.adsorption_mode == "chemisorption"].set_index("name")
    s["phys"] = m[m.adsorption_mode == "physisorption"].set_index("name")
    s["react"] = m[~m.drug_intact].set_index("name")
    ph = s["phys"]
    s["cat"] = ph[ph.formal_charge > 0]
    s["neu"] = ph[ph.formal_charge == 0]
    s["rho_charge"] = stats.spearmanr(ph.formal_charge, ph.delta_Eint_kcal_mol)
    s["rho_size"] = stats.spearmanr(ph.n_heavy_atoms, ph.delta_Eint_kcal_mol)
    s["rho_vina"] = stats.spearmanr(m.vina_8FUG_kcal_mol, m.delta_Eint_kcal_mol)
    s["rho_vina_size"] = stats.spearmanr(m.vina_8FUG_kcal_mol, m.n_heavy_atoms)
    c = d["carrier"]
    s["E_shift"] = (c["E_build_Eh"] - c["E_Eh"]) * 627.509
    o = d["old"]
    E0 = json.loads((BASE / "data" / "processed" / "carrier_B40H15_collapsing_2026-09-24.json").read_text())["E_Eh"]
    s["old_drop"] = ((o.E_carrier_frozen_Eh - E0) * 627.509).min()
    return s


# ------------------------------------------------------------------ sections
def front(doc):
    k.para(doc, f"**{TITLE}**", align="left", size=15, space_after=12)
    k.para(doc, f"{AUTHOR}^{{*}}", align="left", space_after=2)
    k.para(doc, AFFIL, align="left", size=10, space_after=2)
    k.para(doc, f"^{{*}}Corresponding author: {EMAIL}; ORCID {ORCID}", align="left", size=10, space_after=14)


def abstract(doc, d, c):
    s, q = stats_(d), d["q"]
    ph, ch = s["phys"], s["chem"]
    k.heading(doc, "Abstract")
    k.labelled(doc, "Context",
               "Tau aggregation inhibitors and imaging probes act on paired helical filaments (PHF) in the "
               f"Alzheimer brain. For {s['n']} tau-directed and Alzheimer drugs we docked each compound at the "
               "GTP-1 tracer site of a PHF cryo-EM structure and computed its adsorption on β_{12} borophene, a "
               "candidate carrier. Docking into the "
               "ligand-free fibril does not reproduce the tracer pose (RMSD "
               f"{f1(s['r_free_x'])} and {f1(s['r_free_s'])} Å); the docking scores are therefore exploratory. "
               "An unconstrained finite borophene flake collapses into a compact boron cluster during adsorption, "
               "so the carrier was modelled as a planar B_{44}H_{16} sheet held at the β_{12} lattice, as it is "
               f"when grown on Ag(111). On this sheet {len(ph)} drugs physisorb and {len(ch)} chemisorb; Congo Red "
               "and tideglusib react. Physisorption is strongest for the cationic phenothiazinium and "
               f"benzothiazolium dyes (Δ*E*_{{int}} {f1(s['cat'].delta_Eint_kcal_mol.max())} to "
               f"{f1(s['cat'].delta_Eint_kcal_mol.min())} kcal mol^{{−1}}) and weakest for the symptomatic "
               "Alzheimer drugs. A ridge model on four pre-selected descriptors explains part of the variance "
               f"(*Q*^{{2}}_{{CV}} = {f2(q['Q2_CV'])}, Y-scrambling *p* = {q['Y_scrambling']['p']:.3f}), mainly "
               "through formal charge.")
    k.labelled(doc, "Methods",
               "Structures were taken from PubChem. Drugs were docked with AutoDock Vina 1.2.7 into the PHF "
               "segment of PDB 8FUG at the GTP-1 site, with four redocking controls. The carrier and all "
               "complexes were computed with GFN2-xTB (xtb 6.7.1), with the boron atoms of the sheet restrained "
               "to the β_{12} lattice; each drug was adsorbed from four relaxed orientations and the regime was "
               "assigned from drug–carrier bond formation. Ridge QSPR models were assessed by nested 5×5 "
               "cross-validation, 1,000-fold Y-scrambling and a leverage applicability domain.")
    k.para(doc, "**Keywords** Tau · Paired helical filament · Borophene · GFN2-xTB · Molecular docking · "
                "Drug delivery", align="left")


def introduction(doc, c):
    k.heading(doc, "Introduction")
    k.para(doc,
           "In Alzheimer's disease the microtubule-associated protein tau assembles into paired helical "
           "filaments (PHF), whose amount and spread track neurodegeneration " + c("wang2016", "congdon2018") +
           ". Cryo-EM structures of PHF isolated from patient brain show a C-shaped core of residues 306–378, "
           "identical across patients " + c("fitzpatrick2017", "shi2021") + ". Small molecules that bind these "
           "filaments include imaging probes and aggregation inhibitors: phenothiazines such as methylene blue "
           "and its reduced form hydromethylthionine " + c("wischik1996", "baddeley2015", "gauthier2016") +
           ", polyphenols such as curcumin and EGCG " + c("rane2017", "seidler2018") + ", and oligomer "
           "modulators such as anle138b " + c("wagner2013") + ". The PET tracer GTP-1 has been resolved bound "
           "to PHF at 2.7 Å, stacked along the fibril axis in a cleft of each protofilament (PDB 8FUG) " +
           c("merz2023") + ", which provides a structurally defined site for docking.", indent=True)
    k.para(doc,
           "Many of these compounds reach the brain poorly, and two-dimensional materials have been proposed as "
           "carriers. Borophene, a single atomic layer of boron, has been grown on Ag(111) in several "
           "polymorphs, of which β_{12} is a triangular lattice with one sixth of the sites vacant "
           + c("mannix2015", "feng2016", "kong2017") + ". It is metallic and, unlike graphene or boron nitride, "
           "its boron atoms are electron-deficient, so it may bind drugs through lone-pair donation as well as "
           "by dispersion. Freestanding borophene is not stable: it exists on its growth substrate, which keeps "
           "it planar " + c("kong2017") + ".", indent=True)
    k.para(doc,
           "This study characterises a cohort of tau-directed and Alzheimer drugs on two fronts: docking at the "
           "GTP-1 site of PHF, with explicit redocking controls, and GFN2-xTB adsorption on a β_{12} borophene "
           "model whose geometry is held at the supported lattice. It reports which drugs physisorb or "
           "chemisorb, which react with the sheet, and whether the interaction energy can be anticipated from "
           "four molecular descriptors under leak-free validation.", indent=True)


def methods(doc, d, c):
    s = stats_(d)
    car = d["carrier"]
    k.heading(doc, "Methods")
    k.heading(doc, "Compound set", 2)
    fam = d["m"].family.value_counts()
    k.para(doc,
           f"The cohort comprises {s['n_lib']} compounds. Every structure was retrieved from PubChem by name "
           + c("kim2021_pubchem") + " and its identity checked by InChIKey (Table S1). Counter-ions were removed "
           "and acids and bases neutralised where a neutral form exists; the phenothiazinium dyes and "
           "thioflavin T keep their permanent positive charge. Thioflavin S, a mixture of sulfonated oligomers "
           f"rather than a single compound, was excluded, leaving {s['n']} drugs in four families: "
           f"{fam['Dyes and imaging probes']} dyes and imaging probes, {fam['Polyphenols']} polyphenols, "
           f"{fam['Aggregation / kinase modulators']} aggregation or kinase modulators and "
           f"{fam['Symptomatic AD drugs']} symptomatic Alzheimer drugs.", indent=True)

    k.heading(doc, "Docking at the GTP-1 site of PHF", 2)
    k.para(doc,
           "The receptor is the 23-chain PHF segment of PDB 8FUG " + c("merz2023", "berman2000") + ", "
           "prepared with PDBFixer " + c("eastman2017") + " (all GTP-1 copies removed, missing atoms added, "
           "protonation at pH 7.4) and typed with Meeko. The docking box, a 22 Å cube, was centred on the GTP-1 "
           "copy bound to chain K, in the middle of the stack of protofilament 1 and away from the segment ends. "
           "Ligands were built from SMILES with RDKit " + c("rdkit") + " (ETKDGv3 " + c("wang2020_etkdg") +
           ", MMFF94 " + c("halgren1996") + "); up to five ring conformers per drug were docked and the best "
           "score was kept. Docking used AutoDock Vina 1.2.7 " + c("trott2010", "eberhardt2021") +
           " (exhaustiveness 16, nine modes, fixed seed). Four controls were run: GTP-1 redocked from its "
           "cryo-EM conformation and from SMILES, each into the ligand-free fibril and into the fibril with "
           "the neighbouring GTP-1 copies of the stack kept as receptor atoms. RMSDs were computed on heavy "
           "atoms with symmetry correction and without re-alignment. The cohort was docked into the "
           "ligand-free fibril. Residues within 4.0 Å of the top pose were counted as contacts, and N/O pairs "
           "within 3.5 Å as polar contacts.", indent=True)

    k.heading(doc, "β_{12} borophene model", 2)
    k.para(doc,
           "A first model, an H-terminated B_{40}H_{15} flake relaxed without constraints, was a local minimum "
           "(no imaginary frequencies) but not a stable carrier: in every adsorption complex it contracted into "
           f"a more compact boron cluster, lowering its own energy by up to {f1(-s['old_drop'])} kcal mol^{{−1}}, "
           "so that the computed binding energies were dominated by the reconstruction of the carrier. Because "
           "borophene exists only on a supporting metal, which holds it planar " + c("feng2016", "kong2017") +
           ", the carrier was rebuilt as a substrate-supported sheet. A planar β_{12} lattice (triangular "
           f"lattice, B–B {car['lattice_A']['d_BB']:.2f} Å, rectangular cell *a* = {car['lattice_A']['a']:.2f} Å, "
           f"*b* = {car['lattice_A']['b']:.2f} Å, one vacancy per six sites " + c("feng2016", "kong2017") + ") was cut to "
           f"a {car['formula'].replace('B44', 'B_{44}').replace('H16', 'H_{16}')} flake, and every edge boron "
           "that had lost a neighbour was capped with hydrogen. When the whole lattice is scaled uniformly, the "
           "GFN2-xTB energy of this planar flake is lowest at the reference B–B distance (Table S3), so the "
           "model is not strained against the method. Throughout, the boron atoms were held at the lattice by "
           "harmonic restraints on their interatomic distances (xtb $constrain, force constant 5 E_{h} "
           "bohr^{−2}), and hydrogens and drugs relaxed freely. The singlet lies "
           f"{f1((car['E_triplet_Eh'] - car['E_singlet_Eh']) * 627.509)} kcal mol^{{−1}} below the triplet, and "
           "the vanishing HOMO–LUMO gap reflects the metallic character of the sheet; all calculations on the "
           "carrier used Fermi smearing (1500 K).", indent=True)

    k.heading(doc, "Adsorption calculations", 2)
    k.para(doc,
           "All quantum-chemical calculations used GFN2-xTB " + c("bannwarth2019") + " (xtb 6.7.1, D4 "
           "dispersion " + c("caldeweyher2019") + ") in the gas phase, with the drug in its PubChem charge "
           "state. Each drug was relaxed in isolation from its lowest MMFF94 conformer, placed with its mean "
           "plane parallel to the sheet 3.2 Å above it, and four rotations about the surface normal (0, 90, 180 "
           "and 270°) were relaxed; the lowest-energy converged pose with a closest drug–carrier heavy-atom "
           "contact of 1.25–4.0 Å was kept. On that pose the interaction energy Δ*E*_{int} = *E*_{complex} − "
           "*E*_{carrier} − *E*_{drug} was computed with both fragments frozen at the complex geometry. The "
           "adsorption energy Δ*E*_{ads} is referenced to the relaxed drug and to the relaxed restrained sheet. "
           "Within the restraints the sheet can buckle by about 0.1 Å, and it has several such minima; the "
           "reference energy is the lowest one, obtained by re-relaxing the carrier from each of the "
           f"{len(d['ref_scan'])} physisorbed complexes ({f1(s['E_shift'])} kcal mol^{{−1}} below the "
           "initial planar geometry). A complex was classed as chemisorbed when at least one drug–carrier pair "
           "was closer than 1.15 times the sum of the covalent radii. The bonding of every drug was compared "
           "before and after adsorption to detect proton transfer or bond breaking.", indent=True)

    k.heading(doc, "Descriptors and QSPR model", 2)
    k.para(doc,
           "Four descriptors were fixed before any model was fitted: molecular weight and topological polar "
           "surface area (RDKit), the static polarizability α(0) of the relaxed isolated drug (GFN2-xTB) and "
           "the formal charge. Conceptual-DFT indices such as electrophilicity " + c("parr1999") + " were not "
           "used, because they are not comparable between neutral and cationic molecules. The target is "
           "Δ*E*_{int} for drugs whose own bonding was unchanged on adsorption. A ridge-regression model "
           "(scikit-learn " + c("pedregosa2011") + ", standardisation and penalty inside one pipeline) was "
           "evaluated by nested cross-validation, with five outer folds for performance and a five-fold inner "
           "loop for the penalty " + c("cawley2010") + ". Chance correlation was tested by repeating the whole "
           "procedure on 1,000 permutations of the target " + c("rucker2007") + ", and the applicability domain "
           "follows OECD principle 3 " + c("oecd2007", "gramatica2007", "tropsha2010") + " (leverage warning "
           "*h*^{*} = 3(*p* + 1)/*n*, standardised residuals within ±3).", indent=True)


def results(doc, d, c):
    s, q, m = stats_(d), d["q"], d["m"]
    ph, ch, rx = s["phys"], s["chem"], s["react"]
    k.heading(doc, "Results and discussion")
    k.para(doc, "The workflow is summarised in Fig. 1. All quantities are computed; none is fitted to "
                "experimental data.", indent=True)
    k.figure(doc, FIG / "Fig1.png", 1,
             "Workflow of the study: PubChem compound set, docking at the GTP-1 site of Alzheimer PHF (PDB "
             "8FUG), the supported β_{12} borophene model, GFN2-xTB adsorption, and the QSPR model with its "
             "validation. The inset of the last stage is the out-of-fold parity plot of Fig. 7a.")

    k.heading(doc, "Docking at the GTP-1 site", 2)
    fr, po, v, fv = s["freq"], s["polar"], s["vina"], s["fam_v"]
    k.para(doc,
           "Redocking shows how far docking into a fibril can be trusted. Into the ligand-free fibril, GTP-1 "
           f"docked {f1(s['r_free_x'])} Å (cryo-EM conformation) and {f1(s['r_free_s'])} Å (from SMILES) away "
           "from its observed pose, both outside the usual 2 Å criterion. With the neighbouring GTP-1 copies of "
           f"the stack kept, the cryo-EM conformation was recovered at {f2(s['r_stack_x'])} Å, whereas the "
           f"ligand rebuilt from SMILES reached {f2(s['r_stack_s'])} Å (Fig. 2a). The observed pose therefore "
           "depends on stacking against the ligands above and below it, which a single-ligand docking into the "
           "empty cleft does not reproduce " + c("merz2023") + ". The cohort scores reported below were "
           "obtained on the ligand-free fibril and are an exploratory ranking of how well each drug fits the "
           "cleft, not predicted binding modes.", indent=True)
    k.para(doc,
           f"Nearly all top poses sit on the same three residues: {fr.index[0]} ({fr.iloc[0]:.0f}% of drugs), "
           f"{fr.index[1]} ({fr.iloc[1]:.0f}%) and {fr.index[2]} ({fr.iloc[2]:.0f}%), with polar contacts to "
           f"Ser352 in {po['Ser352']:.0f}% and to Gln351 in {po['Gln351']:.0f}% of drugs (Fig. 2b). These "
           "residues line the GTP-1 cleft on consecutive rungs of the stack (chains "
           f"{', '.join(s['chains'])}). Scores ranged from {f1(v.iloc[-1])} to {f1(v.iloc[0])} kcal mol^{{−1}} "
           f"(Fig. 3a); the extended azo dyes {v.index[0].replace('Chrysamine', 'chrysamine')} and "
           f"{v.index[1]} scored best, and the small "
           "symptomatic drugs worst (family median "
           f"{f1(fv['Symptomatic AD drugs'])} against {f1(fv['Dyes and imaging probes'])} for dyes and probes). "
           "Over the whole cohort the score follows molecular size (Spearman "
           f"ρ = {f2(s['rho_vina_size'].statistic)} with heavy-atom count), as expected for an empirical score "
           + c("kitchen2004") + ".", indent=True)
    k.figure(doc, FIG / "Fig2.png", 2,
             "Docking at the GTP-1 site of PHF (PDB 8FUG). **a** Cryo-EM GTP-1 bound to chain K (grey), with the "
             "stacked copies above and below (green), the pose redocked from its cryo-EM conformation with the "
             f"stack present (orange, {f2(s['r_stack_x'])} Å) and the pose from SMILES (teal, "
             f"{f2(s['r_stack_s'])} Å); contacting side chains as thin sticks. **b** Percentage of the "
             f"{s['n']} drugs whose top pose contacts each residue (light bars, any heavy atom within 4.0 Å; "
             "dark bars, N/O pairs within 3.5 Å).")
    k.figure(doc, FIG / "Fig3.png", 3,
             "Docking scores and their relation to adsorption. **a** Vina score of every drug in the ligand-free "
             "fibril, coloured by family. **b** −Δ*E*_{int} on the β_{12} sheet against the Vina score; open "
             "symbols, chemisorbed complexes.")

    k.heading(doc, "The supported β_{12} sheet", 2)
    car = d["carrier"]
    n_above = int((m.E_LUMO_eV < (car["HOMO_eV"] + car["LUMO_eV"]) / 2).sum())
    k.para(doc,
           "The planar B_{44}H_{16} flake (Fig. 4a) is metallic at the GFN2-xTB level, with a Fermi level near "
           f"{f2((car['HOMO_eV'] + car['LUMO_eV']) / 2)} eV. For {n_above} of the {s['n']} drugs this lies above "
           "the drug LUMO, and for the cationic dyes it lies several eV above both frontier orbitals (Fig. 4b); "
           "at this level of theory the sheet is therefore an electron donor toward most drugs. Absolute "
           "GFN2-xTB orbital energies are not quantitative, and the comparison is qualitative.", indent=True)
    k.figure(doc, FIG / "Fig4.png", 4,
             "Carrier and frontier orbitals. **a** Planar β_{12} B_{44}H_{16} flake (B pink, H white). **b** "
             "GFN2-xTB HOMO (circles) and LUMO (squares) of the relaxed drugs, sorted by gap and coloured by "
             "family; dashed line, Fermi level of the sheet.")

    k.heading(doc, "Adsorption", 2)
    cat, neu = s["cat"], s["neu"]
    k.para(doc,
           f"On the supported sheet {len(ph)} of the {s['n']} drugs physisorbed, with closest contacts of "
           f"{f2(ph.min_contact_A.min())}–{f2(ph.min_contact_A.max())} Å and Δ*E*_{{int}} from "
           f"{f1(ph.delta_Eint_kcal_mol.max())} to {f1(ph.delta_Eint_kcal_mol.min())} kcal mol^{{−1}} (Fig. 5, "
           "Table 1). The strongest physisorbed drugs are the four cations (methylene blue, azure A, toluidine "
           f"blue O and thioflavin T, {f1(cat.delta_Eint_kcal_mol.max())} to {f1(cat.delta_Eint_kcal_mol.min())} "
           "kcal mol^{−1}), which bind a metallic sheet through electrostatic polarisation in addition to "
           f"dispersion; among neutral physisorbed drugs Δ*E*_{{int}} ranged from "
           f"{f1(neu.delta_Eint_kcal_mol.max())} to {f1(neu.delta_Eint_kcal_mol.min())} kcal mol^{{−1}}. Across "
           f"physisorbed drugs the energy follows the formal charge (ρ = {f2(s['rho_charge'].statistic)}) but "
           f"not size (ρ = {f2(s['rho_size'].statistic)} with heavy-atom count). The families differ "
           f"(Kruskal–Wallis *p* = {sci(s['kw_e'].pvalue)}; Fig. 5b): the median Δ*E*_{{int}} is "
           f"{f1(s['fam_e']['Dyes and imaging probes'])} kcal mol^{{−1}} for dyes and probes, "
           f"{f1(s['fam_e']['Polyphenols'])} for polyphenols, {f1(s['fam_e']['Aggregation / kinase modulators'])} "
           f"for modulators and {f1(s['fam_e']['Symptomatic AD drugs'])} for symptomatic drugs. In the gas "
           "phase the cation–sheet attraction is not screened; in water it would be much weaker, so the gap "
           "between cationic and neutral drugs is an upper bound.", indent=True)
    bonds = {n_: ", ".join(re.sub(r"\d", "", b) for b in eval(r.drug_carrier_bonds))
             for n_, r in ch.iterrows()}
    k.para(doc,
           f"{['No', 'One', 'Two', 'Three', 'Four', 'Five'][len(ch)] if len(ch) < 6 else len(ch)} drugs chemisorbed. Chrysamine G bonded to boron through three of its four azo nitrogens "
           f"(Δ*E*_{{int}} {f1(ch.loc['Chrysamine G', 'delta_Eint_kcal_mol'])} kcal "
           "mol^{−1}) without any change in its own bonding (Fig. 6a). Congo Red transferred the proton of a "
           "sulfonic acid group to a boron atom and bonded to boron through a sulfonate oxygen "
           f"({f1(ch.loc['Congo Red', 'delta_Eint_kcal_mol'])}; Fig. 6c); as a disodium salt in practice, it "
           "would not carry that proton. Tideglusib opened its thiadiazolidinedione ring at the S–N bond and "
           f"bonded to boron through sulfur ({f1(ch.loc['Tideglusib', 'delta_Eint_kcal_mol'])}). For these two "
           "drugs Δ*E*_{int} contains a reaction energy. The electron-deficient boron sheet thus forms B–N bonds "
           "with azo nitrogens, takes up an acidic proton and opens a strained S–N heterocycle, but binds the "
           "phenothiazines, flavonoids and cholinesterase inhibitors by physisorption alone.", indent=True)
    k.figure(doc, FIG / "Fig5.png", 5,
             "Adsorption on the supported β_{12} sheet. **a** −Δ*E*_{int} against the closest drug–carrier "
             "heavy-atom contact, coloured by regime (chemisorption: at least one drug–carrier pair closer than "
             "1.15 times the sum of covalent radii; crosses, drugs whose own bonding changed). **b** −Δ*E*_{int} "
             "by family; boxes show the median and interquartile range.")
    k.figure(doc, FIG / "Fig6.png", 6,
             "Relaxed GFN2-xTB complexes. **a** Chrysamine G (chemisorbed, strongest intact binder). **b** "
             "Azure A (strongest physisorbed). **c** Congo Red, which transferred a proton to the sheet. Colours: "
             "B pink, C grey, N blue, O red, S yellow, H white.")
    rows = []
    for lab, g in (("Dyes and imaging probes", "Dyes and imaging probes"), ("Polyphenols", "Polyphenols"),
                   ("Aggregation / kinase modulators", "Aggregation / kinase modulators"),
                   ("Symptomatic AD drugs", "Symptomatic AD drugs")):
        x = m[m.family == g]
        rows.append([lab, str(len(x)), f"{int((x.adsorption_mode == 'physisorption').sum())} / "
                                       f"{int((x.adsorption_mode == 'chemisorption').sum())}",
                     f1(x.vina_8FUG_kcal_mol.median()),
                     f"{f1(x.delta_Eint_kcal_mol.median())} ({f1(x.delta_Eint_kcal_mol.max())} to "
                     f"{f1(x.delta_Eint_kcal_mol.min())})", f1(x.delta_Eads_kcal_mol.median())])
    k.table(doc, (1, "Docking and adsorption by drug family."),
            ["Family", "*n*", "Phys. / chem.", "Vina", "Δ*E*_{int}", "Δ*E*_{ads}"], rows, align="lccccc",
            font=8.5, note="Energies in kcal mol^{−1}; Vina, median score in PDB 8FUG; Δ*E*_{int}, median (range); "
                           "Δ*E*_{ads}, median. Per-drug values are in Table S2.")

    k.heading(doc, "Docking score and adsorption", 2)
    k.para(doc,
           "Drugs that score well at the GTP-1 site also bind the sheet more strongly (Spearman "
           f"ρ = {f2(s['rho_vina'].statistic)}, *p* = {sci(s['rho_vina'].pvalue)}; Fig. 3b). The trend is carried "
           "by the extended, planar dyes, which fit the elongated cleft and stack on the sheet, and by the "
           "small aliphatic symptomatic drugs, which do neither well; it should not be read as a mechanistic "
           "link between the two endpoints.", indent=True)

    k.heading(doc, "QSPR model of the interaction energy", 2)
    k.para(doc,
           f"For the {q['n']} drugs that kept their bonding, the ridge model reaches an out-of-fold "
           f"*Q*^{{2}}_{{CV}} of {f2(q['Q2_CV'])} (RMSE {f1(q['RMSE'])}, MAE {f1(q['MAE'])} kcal mol^{{−1}}; "
           f"Table 2, Fig. 7). Only {q['Y_scrambling']['p'] * 100:.1f}% of 1,000 permuted targets did as well, "
           "so the signal is real, but it is weak and uneven: the outer-fold *Q*^{2} values range from "
           f"{f2(min(q['Q2_folds']))} to {f2(max(q['Q2_folds']))}. The largest standardised coefficient is that "
           f"of formal charge ({f2(q['coef_std']['charge'])}), which separates the four cations from the rest; "
           f"chrysamine G, bound through three N–B bonds, is under-predicted by about "
           f"{f1(abs(d['oof'].set_index('name').loc['Chrysamine G', 'dE_int'] - d['oof'].set_index('name').loc['Chrysamine G', 'oof_pred']))} "
           "kcal mol^{−1} and is a residual outlier, and anle138b lies beyond the leverage threshold. The model "
           "captures the cation effect but not the chemistry of individual contacts, and it should not be used "
           "to rank new compounds.", indent=True)
    k.figure(doc, FIG / "Fig7.png", 7,
             "Validation of the QSPR model for Δ*E*_{int}. **a** Out-of-fold predictions from nested 5×5 "
             "cross-validation against the GFN2-xTB values; the shaded band is ±10% of the range. **b** "
             "Williams plot (dotted line, *h*^{*}; open symbols, outside the applicability domain). **c** "
             "Distribution of *Q*^{2}_{CV} over 1,000 permuted targets, with the value of the model.")
    qn = {"alpha": "α(0)", "charge": "charge"}
    k.table(doc, (2, "Ridge QSPR model of Δ*E*_{int} (MW, TPSA, α(0), formal charge) under nested 5×5 "
                     "cross-validation."),
            ["Quantity", "Value"],
            [["*n*, *p*", f"{q['n']}, {q['p']}"],
             ["*Q*^{2}_{CV} (out-of-fold)", f2(q["Q2_CV"])],
             ["*Q*^{2} of the outer folds", ", ".join(f2(x) for x in q["Q2_folds"])],
             ["RMSE / MAE (kcal mol^{−1})", f"{f1(q['RMSE'])} / {f1(q['MAE'])}"],
             ["Ridge penalty λ (inner CV)", f"{q['alpha_final']:g}"],
             ["Standardised coefficients", ", ".join(f"{qn.get(a, a)} {f2(b)}" for a, b in q["coef_std"].items())],
             ["Y-scrambling, 1,000 runs: mean *Q*^{2}; *p*",
              f"{f2(q['Y_scrambling']['mean_Q2'])}; {q['Y_scrambling']['p']:.3f}"],
             ["Applicability domain: *h*^{*}; outside", f"{f2(q['AD']['h_star'])}; {', '.join(q['AD']['outside'])}"],
             ["Excluded (bonding changed)", ", ".join(q["excluded_not_intact"])]],
            align="ll", font=8.5)

    k.heading(doc, "Limitations", 2)
    k.para(doc,
           "The docking controls fail on the ligand-free fibril, so the cohort scores are exploratory. The "
           "carrier is a finite flake whose boron atoms are held at the lattice of Ag-supported β_{12} "
           "borophene; the metal substrate itself, which interacts electronically with the sheet " + c("feng2016") +
           ", is not included. All energies are gas-phase GFN2-xTB values " + c("bannwarth2019") + "; the "
           "interactions of the cationic dyes and the bond-forming events in particular should be revisited "
           "with DFT and implicit or explicit solvent. None of the results has been compared with experiment.",
           indent=True)


def conclusions(doc, d, c):
    s, q = stats_(d), d["q"]
    k.heading(doc, "Conclusions")
    k.para(doc,
           "Docking into the ligand-free PHF of PDB 8FUG places nearly every drug on Gln351, Ser352 and Lys353, "
           "but it does not reproduce the stacked GTP-1 pose, and its scores are exploratory. A finite, "
           "unconstrained borophene flake is not a usable carrier model at the GFN2-xTB level, because it "
           "collapses during adsorption; a planar sheet held at the β_{12} lattice, as on Ag(111), is. On this "
           f"sheet {len(s['phys'])} of {s['n']} drugs physisorb, the cationic dyes most strongly, and "
           f"{len(s['chem'])} chemisorb, of which Congo Red and tideglusib react. The interaction energy is "
           f"partly predictable from four descriptors (*Q*^{{2}}_{{CV}} = {f2(q['Q2_CV'])}), mainly through "
           "formal charge. Calculations with the metal substrate and with solvent are the next step before "
           "any experimental test.", indent=True)


def declarations(doc):
    k.heading(doc, "Declarations")
    for label, text in (
        ("Author contribution", f"{AUTHOR} conceived the study, performed all calculations and analyses, "
                                "and wrote the manuscript."),
        ("Funding", "No funding was received for this work."),
        ("Data availability", "All input structures, relaxed geometries, docking poses, xtb outputs and "
                              f"result tables are available at {REPO}."),
        ("Code availability", f"The complete pipeline, which regenerates every number, table and figure "
                              f"of this article from the raw inputs, is available at {REPO} under the MIT "
                              "licence. Software: xtb 6.7.1, AutoDock Vina 1.2.7, Meeko, RDKit, PDBFixer/"
                              "OpenMM, scikit-learn, PyMOL (open source)."),
        ("Ethics approval", "Not applicable."),
        ("Consent to participate", "Not applicable."),
        ("Consent for publication", "Not applicable."),
        ("Competing interests", "The author declares no competing interests."),
    ):
        k.labelled(doc, label, text)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    d = load()
    c = Cites()
    doc = k.new_document()
    front(doc)
    abstract(doc, d, c)
    introduction(doc, c)
    methods(doc, d, c)
    results(doc, d, c)
    conclusions(doc, d, c)
    declarations(doc)
    k.references(doc, c.list())
    out = OUT / "Manuscript_Tau_Borophene_JMM.docx"
    doc.save(out)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
