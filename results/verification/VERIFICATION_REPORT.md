# Independent verification of the Tau / β12-borophene manuscript (2026-09-28)

Each check recomputes the reported values from the raw outputs (`src/verification/` and
the inline checks recorded here). The only reused function is the extraction of the
crystal ligand from 8FUG.

| What | How | Result |
|---|---|---|
| Adsorption energies (28 complexes) | xtb 6.7.1 GFN2 single points re-run (complex and carrier with Fermi smearing at 1500 K, as in the pipeline) | identical (max 9e-12 Eh); ΔE_int identical in 28/28 |
| Pose choice, regime, drug integrity, closest contact | from the optimisation logs; bond analysis re-implemented | 28/28 each |
| Carrier reference for ΔE_ads | lowest of the 25 restrained re-relaxations | FDDNP, -61.510526 Eh = value used; ΔE_ads reproduced (max 0.0005 kcal/mol) |
| Singlet-triplet order of the sheet | carrier.json energies | singlet 2.1 kcal/mol below the triplet (text identical) |
| Docking scores | best REMARK VINA RESULT of each pose file | 28/28 |
| The four redocking controls | RMSD recomputed with RDKit | 6.731, 6.282, 0.713 and 2.814 Å, identical |
| Identity of the 29 structures | InChIKey vs PubChem record of the CID | 24 identical; the other 5 are dyes stored by PubChem as salts (Cl⁻, Na⁺), and they are the same molecule once the counter-ion is removed, as the Methods state |
| QSPR | full nested CV and 1,000 permutations re-run (outputs rewritten 2026-09-28) | identical to the reported results (no change in any file) |
| References | DOI/title against Crossref; relevance read sentence by sentence | 0 problems after one fix. "Polyphenols such as curcumin and EGCG" cited Seidler et al. 2018, which is about peptide inhibitors, not EGCG; it was replaced by Wobst et al. 2015 (EGCG inhibits tau aggregation). The 2.7 Å resolution of 8FUG was confirmed in the PDB. |

## Open before submission
- The rebuilt code and data are on the local branch `rebuild-2026-09` only. The
  availability statement is true only after that branch is pushed (with the author's approval).
