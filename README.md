# Machine Learning-Driven Nano-QSAR and Quantum Chemical Design of Functionalized 2D Borophene Nanocarriers for Alzheimer's Tau-Targeted Therapeutics

[![License: MIT](https://img.shields.io/badge/License-MIT-teal.svg)](https://opensource.org/licenses/MIT)
[![Python 3.12](https://img.shields.io/badge/Python-3.12-blue.svg)](https://www.python.org/)
[![AutoDock Vina](https://img.shields.io/badge/Docking-AutoDock%20Vina-orange.svg)](https://github.com/ccsb-scripps/AutoDock-Vina)
[![XAI: SHAP](https://img.shields.io/badge/Explainability-SHAP-purple.svg)](https://github.com/shap/shap)

**Authors**: Andrés Monreal Hernández, Sara Lizbeth Franco Amaya, Carlos Ivanhoe Martínez Osorio  
**Affiliation**: Universidad Estatal de Sonora, Hermosillo, Sonora, México  

---

## 📌 Abstract

Pathological aggregation and hyperphosphorylation of Tau protein into paired helical filaments (PHFs) and neurofibrillary tangles (NFTs) constitute the primary cytopathological hallmark and driver of cognitive decline in Alzheimer's Disease (AD). Here, we present an integrated multi-scale quantum mechanical, molecular docking, and explainable machine learning (Nano-QSAR / XAI) pipeline to evaluate **two-dimensional (2D) pristine and functionalized borophene nanosheets** as targeted nanocarriers for brain delivery of anti-Tau therapeutics.

### Key Methodological Highlights:
- **Cryo-EM Receptor Target**: Molecular docking against human Alzheimer's paired helical filament Tau protofibril (PDB ID: 6VHL, cryo-EM $2.80\ \text{Å}$ resolution).
- **Curated Neuro-Therapeutic Cohort**: 28 clinical-stage Tau aggregation inhibitors, kinase inhibitors, and neuroprotective agents.
- **Quantum Conceptual DFT (CDFT)**: Electronic hardness ($\eta$), softness ($S$), chemical potential ($\mu$), electrophilicity ($\omega$), and adsorption energetics across pristine and functionalized borophene models.
- **Benchmark ML/QSAR Models**: Multi-model evaluation (Random Forest, Gradient Boosting, Extra Trees, Ridge, SVR) with nested cross-validation and OECD Principle 3 applicability domain (Williams plot).
- **Explainable AI (SHAP)**: Global and local Shapley additive explanations identifying frontier orbital energies and polar surface area as primary drivers of therapeutic-surface stabilization.

---

## 🔬 Repository Architecture

```
├── data/
│   ├── processed/                             # Processed datasets and descriptor matrices
│   └── raw/                                   # PDB 6VHL receptor and 28 ligand PDBQT coordinates
├── figures/                                   # High-resolution publication figures (300 DPI)
├── manuscript/
│   └── Beilstein_Manuscript_Tau_Borophene_Monreal_Hernandez_et_al.docx
├── results/
│   ├── docking/                               # Real Vina binding scores and contact residues
│   └── models/                                # QSAR benchmark summaries and SHAP rankings
├── src/
│   ├── descriptors/                           # CDFT & molecular descriptor computation
│   ├── docking/                               # Docking execution & residue contact extraction
│   ├── ml_models/                             # QSAR regression & applicability domain scripts
│   └── visualization/                         # Manuscript & figure compilation pipelines
├── run_entire_tau_study.py                    # Master execution workflow
└── README.md
```

---

## ⚙️ Quickstart & Execution

```bash
git clone https://github.com/sircalch/borophene-alzheimer-tau-ai.git
cd borophene-alzheimer-tau-ai

# Create virtual environment & install requirements
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
pip install rdkit scikit-learn pandas numpy scipy matplotlib seaborn shap python-docx

# Execute end-to-end reproducible pipeline
python run_entire_tau_study.py
```

---

## 📜 Citation

```bibtex
@article{MonrealHernandez2026_Tau_Borophene,
  title={Machine Learning-Driven Nano-QSAR and Quantum Chemical Design of Functionalized 2D Borophene Nanocarriers for Alzheimer's Tau-Targeted Therapeutics},
  author={Monreal Hern{\'a}ndez, Andr{\'e}s and Franco Amaya, Sara Lizbeth and Mart{\'i}nez Osorio, Carlos Ivanhoe},
  journal={Beilstein Journal of Nanotechnology / Submitted},
  year={2026},
  url={https://github.com/sircalch/borophene-alzheimer-tau-ai}
}
```

## 📄 License
Released under the [MIT License](LICENSE).
