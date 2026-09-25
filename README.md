# Tau-directed drugs on β12 borophene: docking at the GTP-1 site of Alzheimer paired helical filaments and GFN2-xTB adsorption on a substrate-supported sheet

Code and data for the article above (single author: Andrés Monreal Hernández, Universidad Estatal de Sonora,
ORCID 0009-0009-1207-8597). Target journal: *Journal of Molecular Modeling*.

## Main results (all computed by this pipeline)

| Quantity | Value |
|---|---|
| Drugs | 28 (Thioflavin S excluded: mixture) |
| GTP-1 redocking, ligand-free fibril (cryo-EM conf. / SMILES) | 6.73 / 6.28 Å (fail; cohort docking exploratory) |
| GTP-1 redocking, tracer stack present (cryo-EM conf. / SMILES) | 0.71 / 2.81 Å |
| Carrier | planar β12 B44H16, boron restrained to the lattice (Ag-supported model) |
| Physisorbed / chemisorbed / drug reacted | 25 / 3 / 2 |
| Physisorbed ΔE_int range, kcal/mol | -11.1 to -99.1 |
| QSPR (ridge; MW, TPSA, α(0), charge): Q²_CV / Y-scrambling p | 0.23 / 0.004 |

An unconstrained B40H15 flake (earlier model) collapses into a compact cluster during adsorption; its results
are kept only as evidence (`results/quantum/adsorption_results_B40H15_collapsing_2026-09-24.csv`, Table S6).

## Reproduce

```
pip install -r requirements.txt      # plus xtb 6.7.1 on PATH (or XTB_EXE), PyMOL for renders
python run_entire_tau_study.py       # every step is resumable
```

## Layout

| Path | Content |
|---|---|
| `data/processed/compound_library_pubchem.csv` | 29 compounds: PubChem CID, SMILES, InChIKey, charge |
| `data/processed/carrier.json` | supported β12 carrier (geometry, energies, restraint settings) |
| `data/raw/` | PDB 8FUG, prepared receptors, ligand PDBQT ring-conformer ensembles |
| `calculations/` | xtb inputs/outputs: carrier, lattice scan, carrier reference, every adsorption complex |
| `results/docking/` | redocking controls, poses, scores, residue contacts |
| `results/quantum/` | adsorption results, isolated-drug descriptors, lattice scan |
| `results/qspr/` | QSPR summary, out-of-fold predictions, Y-scrambling |
| `figures/` | Fig1–Fig7 (PDF, 600 dpi PNG and TIFF) and cached 3D renders |
| `manuscript/submission/` | manuscript, supporting information, cover letter (Word) |
| `src/` | pipeline code |

## Licence

MIT (see `LICENSE`).
