"""
curate_tau_dataset.py
Curates a library of 35 clinical-stage and experimental Tau aggregation inhibitors,
neuroprotective modulators, and FDA-approved Alzheimer's disease therapeutics.
"""

import os
import pandas as pd

TAU_DRUG_LIBRARY = [
    # Direct Tau Disaggregators & Phenothiazines
    {"name": "Hydromethylthionine", "class": "Direct Tau Aggregation Inhibitor (LMTX)", "smiles": "CN(C)C1=CC2=C(C=C1)SC3=C2C=CC(=C3)N(C)C", "drugbank_id": "DB12411"},
    {"name": "Methylene Blue", "class": "Phenothiazine Tau Disaggregator", "smiles": "CN(C)c1ccc2nc3ccc(=[N+](C)C)cc3sc2c1.[Cl-]", "drugbank_id": "DB09241"},
    {"name": "Curcumin", "class": "Beta-Sheet Intercalator & Disaggregator", "smiles": "COc1cc(/C=C/C(=O)CC(=O)/C=C/c2ccc(O)c(OC)c2)ccc1O", "drugbank_id": "DB01654"},
    {"name": "EGCG", "class": "Polyphenolic Fibril Remodeler", "smiles": "O=C(Oc1cc(O)cc(O)c1)[C@@H]2[C@H](Oc3cc(O)cc(O)c3C2)c4cc(O)c(O)c(O)c4", "drugbank_id": "DB03561"},
    {"name": "Resveratrol", "class": "SIRT1 Activator & Tau Phosphorylation Reducer", "smiles": "Oc1ccc(/C=C/c2cc(O)cc(O)c2)cc1", "drugbank_id": "DB02709"},
    {"name": "Rosmarinic Acid", "class": "Polyphenolic Tau Destabilizer", "smiles": "O=C(O[C@H](C(=O)O)Cc1ccc(O)c(O)c1)/C=C/c2ccc(O)c(O)c2", "drugbank_id": "DB07971"},
    {"name": "Baicalein", "class": "Flavonoid Tau Oligomerization Inhibitor", "smiles": "O=C1C=C(c2ccccc2)Oc3c1c(O)c(O)c(O)c3", "drugbank_id": "DB04285"},
    {"name": "Myricetin", "class": "Flavonoid Tau Assembly Inhibitor", "smiles": "O=C1C(O)=C(c2cc(O)c(O)c(O)c2)Oc3cc(O)cc(O)c13", "drugbank_id": "DB03975"},
    {"name": "Quercetin", "class": "Flavonoid Neuroprotective Agent", "smiles": "O=C1C(O)=C(c2ccc(O)c(O)c2)Oc3cc(O)cc(O)c13", "drugbank_id": "DB04216"},
    {"name": "Fisetin", "class": "Senolytic & Tau Disaggregator", "smiles": "O=C1C(O)=C(c2ccc(O)c(O)c2)Oc3cc(O)ccc13", "drugbank_id": "DB07795"},
    
    # Kinase (GSK-3beta) Inhibitors targeting Tau Hyperphosphorylation
    {"name": "Tideglusib", "class": "Irreversible GSK-3beta Inhibitor", "smiles": "O=C1N(Cc2ccccc2)C(=O)S/C1=N\\c3cccc4ccccc34", "drugbank_id": "DB12129"},
    {"name": "AZD1080", "class": "GSK-3alpha/beta Inhibitor", "smiles": "CC1=NC(=C(C(=N1)C2=CC=NC=C2)N3CCOCC3)C4=CC=CC=C4F", "drugbank_id": "DB14981"},
    
    # FDA-Approved Symptomatic AD Therapeutics
    {"name": "Memantine", "class": "NMDA Receptor Antagonist", "smiles": "CC12CC3CC(C1)(CC(C3)(C2)N)C", "drugbank_id": "DB00994"},
    {"name": "Donepezil", "class": "AChE Inhibitor", "smiles": "COc1cc2C(=O)C(Cc3ccn(Cc4ccccc4)cc3)Cc2cc1OC", "drugbank_id": "DB00843"},
    {"name": "Rivastigmine", "class": "Dual AChE/BuChE Inhibitor", "smiles": "CCN(C)C(=O)Oc1cccc(c1)[C@@H](C)N(C)C", "drugbank_id": "DB00989"},
    {"name": "Galantamine", "class": "Allosteric nAChR Modulator & AChE Inhibitor", "smiles": "CN1CC[C@]23c4c5ccc(OC)c4O[C@H]2[C@@H](O)C=C[C@H]3CC1C5", "drugbank_id": "DB00674"},
    {"name": "Tacrine", "class": "Cholinesterase Inhibitor", "smiles": "Nc1c2ccccc2nc3CCCCc13", "drugbank_id": "DB00393"},
    
    # BACE1 & Gamma-Secretase Regulators
    {"name": "Verubecestat", "class": "BACE1 Inhibitor", "smiles": "CN1C(=N)N(C)C(=O)C1(C)c2cc(NC(=O)c3ncc(F)cn3)ccc2F", "drugbank_id": "DB12459"},
    {"name": "Atabecestat", "class": "BACE1 Inhibitor", "smiles": "CC1(N=C(N)SC12c3cc(NC(=O)c4ccc(OC(F)F)nc4)ccc3F)C(F)(F)F", "drugbank_id": "DB15021"},
    {"name": "Lanabecestat", "class": "BACE1 Inhibitor", "smiles": "CN1C=C(C=N1)c2cc(NC(=O)c3ncc(OCC(F)(F)F)cn3)ccc2F", "drugbank_id": "DB12458"},
    {"name": "Umibecestat", "class": "BACE1 Modulator", "smiles": "CC1(NC(=N)OC12c3cc(NC(=O)c4ncc(C#N)cn4)ccc3F)C(F)(F)F", "drugbank_id": "DB15233"},
    {"name": "Semagacestat", "class": "Gamma-Secretase Inhibitor", "smiles": "CCC(C)(C)C(=O)N[C@@H](C)C(=O)N[C@H]1CCc2ccccc2C1", "drugbank_id": "DB05459"},
    {"name": "Avagacestat", "class": "Gamma-Secretase Inhibitor", "smiles": "NS(=O)(=O)c1ccc(cc1)C(F)(F)c2ccc(Cl)c(NC(=O)C(F)(F)C(F)(F)F)c2", "drugbank_id": "DB12128"},
    
    # Novel Mechanisms & Neuroprotectors
    {"name": "Bexarotene", "class": "RXR Agonist Clearing Aggregates", "smiles": "CC1(C)CCC(C)(C)c2cc(c(C)cc12)/C(=C)c3ccc(C(=O)O)cc3", "drugbank_id": "DB00611"},
    {"name": "Sodium Phenylbutyrate", "class": "Chemical Chaperone", "smiles": "O=C(O)CCCCc1ccccc1", "drugbank_id": "DB06819"},
    {"name": "TUDCA", "class": "Neuroprotective Bile Acid", "smiles": "C[C@H](CCC(=O)NCCS(=O)(=O)O)[C@H]1CC[C@H]2[C@@H]3[C@H](O)C[C@@H]4C[C@H](O)CC[C@]4(C)[C@H]3CC[C@@]12C", "drugbank_id": "DB09568"},
    {"name": "Anavex 2-73", "class": "Sigma-1 Receptor Agonist", "smiles": "CN(C)CCC(c1ccccc1)C2CCCO2", "drugbank_id": "DB12247"},
    {"name": "GV-971", "class": "Marine Oligosaccharide Surrogate", "smiles": "O=C(O)[C@@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O", "drugbank_id": "DB15520"},
    {"name": "Bryostatin-1", "class": "PKC Epsilon Activator", "smiles": "CCC(C)C(=O)OC1CC(=O)OCC2(O)CC(OC(=O)CC3(O)CC(O)CC(=O)O3)C(O)C(C)O2", "drugbank_id": "DB06161"},
    {"name": "Rifampicin", "class": "Oligomerization Blocker", "smiles": "COc1c2C(=O)c3c(O)c4c(O)c(NC(=O)/C(C)=C/C=C/[C@H](C)[C@H](O)[C@@H](C)[C@@H](O)[C@H](C)[C@H](OC(=O)C)[C@@H](C)/C=C/CO)c(O)c4c(O)c3C(=O)c2c(O)c(CN5CCN(C)CC5)c1O", "drugbank_id": "DB01045"},
    {"name": "Trehalose", "class": "Autophagy Inducer", "smiles": "OC[C@H]1O[C@H](O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O", "drugbank_id": "DB03327"},
    {"name": "Melatonin", "class": "Antioxidant & Tau Kinase Inhibitor", "smiles": "CC(=O)NCCc1c[nH]c2ccc(OC)cc12", "drugbank_id": "DB01065"},
    {"name": "Apigenin", "class": "Flavonoid Anti-inflammatory", "smiles": "O=C1C=C(c2ccc(O)cc2)Oc3cc(O)cc(O)c13", "drugbank_id": "DB07352"},
    {"name": "Luteolin", "class": "Flavonoid Kinase Inhibitor", "smiles": "O=C1C=C(c2ccc(O)c(O)c2)Oc3cc(O)cc(O)c13", "drugbank_id": "DB06841"},
    {"name": "Honokiol", "class": "Neolignan BBB-Penetrating Neuroprotector", "smiles": "C=CCc1ccc(O)c(c1)c2cc(CC=C)ccc2O", "drugbank_id": "DB08149"}
]

def curate():
    base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    out_csv = os.path.join(base_dir, "data", "raw", "tau_drug_library.csv")
    df = pd.DataFrame(TAU_DRUG_LIBRARY)
    df.to_csv(out_csv, index=False)
    print(f"Successfully curated {len(df)} Alzheimer/Tau therapeutics to: {out_csv}")

if __name__ == "__main__":
    curate()
