"""
run_vina_docking.py
===================
AutoDock Vina docking of the Tau cohort into the paired-helical-filament (PHF)
core from Alzheimer's disease brain, PDB 8FUG (cryo-EM, 2.7 A), following
src/docking/vina_protocol.py.

8FUG carries the PET tracer GTP-1 (Y9H) stacked along the fibril axis in a
cleft of each protofilament, which allows the protocol to be validated by
redocking - the earlier receptor (5O3L, 3.4 A) has no ligand, and its box had
been centred on previously docked poses.

Receptor: the whole 23-chain filament segment, all GTP-1 copies removed,
protonated at pH 7.4. Redocking controls are run on this ligand-free fibril
and, because GTP-1 binds as a stack along the fibril axis, also with the
neighbouring GTP-1 copies kept (only the docked copy removed). The cohort is
docked into the ligand-free fibril; its scores are an exploratory ranking. Box: 22 A cube centred on the GTP-1 copy in the middle
of the stack of protofilament 1 (chain K), away from the segment ends.

Usage:  python src/docking/run_vina_docking.py [--controls-only]
Outputs
  data/raw/8FUG_receptor.pdbqt (+ _H.pdb)
  results/docking/real_poses/<drug>_out.pdbqt, <drug>_vina.log
  results/docking/real_vina_docking_summary.csv, redocking_validation.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import vina_protocol as vp  # noqa: E402

BASE = Path(__file__).resolve().parents[2]
RAW = BASE / "data" / "raw"
LIB = BASE / "data" / "processed" / "compound_library_pubchem.csv"
POSES = BASE / "results" / "docking" / "real_poses"
LIGS = RAW / "ligands_pdbqt"
VINA = (BASE / "src" / "docking" / "vina.exe").resolve()

PDB = RAW / "8FUG.pdb"
REF, REF_CHAIN = "Y9H", "K"             # GTP-1, middle of the protofilament-1 stack
REF_SMILES = "FCCC1CCN(CC1)c1ccn2c(n1)nc1ccccc12"
BOX = 22.0
EXCLUDE = {"Thioflavin-S"}              # heterogeneous dye mixture, not a single compound


def main(controls_only=False):
    POSES.mkdir(parents=True, exist_ok=True)
    LIGS.mkdir(parents=True, exist_ok=True)
    center = vp.box_center(PDB, REF, REF_CHAIN)
    _, receptor = vp.prepare_receptor(PDB, None, RAW / "8FUG_receptor")
    print(f"receptor {receptor.name}; box centre {np.round(center, 3)} ({REF}, chain {REF_CHAIN})", flush=True)

    if controls_only or not (BASE / "results" / "docking" / "redocking_validation.csv").exists():
        controls, nha = vp.redock(VINA, PDB, receptor, center, BOX, REF, REF_CHAIN, REF_SMILES, POSES)
        controls = [(f"ligand-free fibril: {c}", s, r, p) for c, s, r, p in controls]
        # GTP-1 binds as a stack along the fibril axis: repeat both controls with the neighbouring
        # GTP-1 copies kept as rigid receptor atoms (only the chain-K copy removed)
        stack = RAW / "8FUG_receptor_tracer_stack.pdbqt"
        chains = sorted({l[21] for l in open(PDB) if l.startswith("HETATM") and l[17:20] == REF})
        extra = [ln for c in chains if c != REF_CHAIN for ln in vp._cofactor_pdbqt_lines(PDB, REF, c, REF_SMILES)]
        stack.write_text(receptor.read_text() + "".join(extra))
        stack_dir = POSES / "tracer_stack_control"
        stack_dir.mkdir(exist_ok=True)
        c2, _ = vp.redock(VINA, PDB, stack, center, BOX, REF, REF_CHAIN, REF_SMILES, stack_dir)
        controls += [(f"tracer stack present: {c}", s, r, f"tracer_stack_control/{p}") for c, s, r, p in c2]
        pd.DataFrame([{
            "pdb_id": "8FUG", "chain": REF_CHAIN, "target_desc": "AD paired helical filament (cryo-EM, 2.7 A)",
            "probe_ligand": "GTP-1 (Y9H)", "control": c, "affinity_kcal_mol": s, "n_heavy_atoms": nha,
            "rmsd_heavy_atom_A": round(r, 3),
            "docking_status": "PASSED (RMSD <= 2.0 A)" if r <= 2.0 else "FAILED (RMSD > 2.0 A)",
            "mapping_method": "RDKit CalcRMS (symmetry-aware, no re-alignment)",
            "pose_file": f"results/docking/real_poses/{p}"} for c, s, r, p in controls]).to_csv(
            BASE / "results" / "docking" / "redocking_validation.csv", index=False)
        for c, s, r, _ in controls:
            print(f"redocking {REF} [{c}]: {s:.2f} kcal/mol, RMSD {r:.2f} A", flush=True)
    if controls_only:
        return

    lib = pd.read_csv(LIB)
    lib = lib[~lib.name.isin(EXCLUDE)]
    rows = []
    for i, r in enumerate(lib.itertuples(), 1):
        s, n_conf = vp.dock_ensemble(VINA, receptor, vp.ligand_pdbqts(r.name, r.smiles, LIGS),
                                     center, BOX, POSES, vp.slug(r.name))
        rows.append({"name": r.name, "class": r[2], "vina_8FUG_kcal_mol": s,
                     "ligand_efficiency": round(-s / r.n_heavy_atoms, 4), "n_ring_conformers": n_conf,
                     "pose_file": f"results/docking/real_poses/{vp.slug(r.name)}_out.pdbqt"})
        print(f"[{i:02d}/{len(lib)}] {r.name:<20} {s:7.3f} kcal/mol  ({n_conf} conformer(s))", flush=True)
    pd.DataFrame(rows).to_csv(BASE / "results" / "docking" / "real_vina_docking_summary.csv", index=False)


if __name__ == "__main__":
    main(controls_only="--controls-only" in sys.argv)
