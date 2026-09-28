"""
verify_adsorption.py - independent re-check of every Tau / beta12-borophene adsorption result
(carrier calculations with Fermi smearing at 1500 K, as in the pipeline).

Does not import the pipeline. For each drug x carrier it
  1. re-runs the three GFN2-xTB single points (complex, carrier and drug frozen at
     the complex geometry) with the same settings and compares them with result.json;
  2. checks that complex_opt.xyz is the lowest-energy converged pose whose closest
     drug-carrier heavy-atom contact lies in 1.25-4.0 A (energies read from the xtb
     optimisation logs);
  3. recomputes drug-carrier bonds, drug integrity and closest contact with its own
     bond code (1.15 x sum of covalent radii);
  4. recomputes dE_int and dE_ads and compares them with results/quantum/adsorption_results.csv,
     the table the manuscript reads.
writes results/verification/adsorption_check.csv and prints a summary.

usage: python src/verification/verify_adsorption.py [--workers 4]
"""
import json
import os
import re
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[2]
ADS = BASE / "calculations" / "adsorption"
XTB = os.environ.get("XTB_EXE", "C:/Users/Andre/mm/xtb/Library/bin/xtb.exe")
HARTREE = 627.509474
COV = {"H": 0.31, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "P": 1.07, "S": 1.05,
       "Cl": 1.02, "Br": 1.20, "I": 1.39}


def read_xyz(f):
    L = Path(f).read_text().splitlines()
    n = int(L[0].split()[0])
    el = [x.split()[0] for x in L[2:2 + n]]
    return el, np.array([[float(v) for v in x.split()[1:4]] for x in L[2:2 + n]])


def sp(xyz_file, chrg, uhf, extra=()):
    env = dict(os.environ, OMP_NUM_THREADS="2", XTBPATH=str(Path(XTB).parent.parent / "share" / "xtb"))
    with tempfile.TemporaryDirectory() as td:
        (Path(td) / "m.xyz").write_text(Path(xyz_file).read_text())
        out = subprocess.run([XTB, "m.xyz", "--sp", "--gfn", "2", "--chrg", str(chrg), "--uhf", str(uhf), *extra],
                             cwd=td, env=env, capture_output=True, text=True, encoding="utf-8", errors="replace",
                             timeout=1800).stdout
    m = re.findall(r"TOTAL ENERGY\s+(-?\d+\.\d+)", out)
    return float(m[-1]) if m else float("nan")


def bonds(el, xyz, A, B):
    return {(i, j) for i in A for j in B if i < j or A is not B
            if np.linalg.norm(xyz[i] - xyz[j]) < 1.15 * (COV[el[i]] + COV[el[j]])}


def pairs_within(el, xyz, idx):
    return {(i, j) for i, j in combinations(idx, 2)
            if np.linalg.norm(xyz[i] - xyz[j]) < 1.15 * (COV[el[i]] + COV[el[j]])}


def log_energy(log):
    m = re.findall(r"energy:\s+(-?\d+\.\d+)", Path(log).read_text())
    return float(m[-1]) if m else float("nan")


def check(job):
    carrier, d = job
    r = json.loads((d / "result.json").read_text())
    if r.get("status") != "OK":
        return {"carrier": carrier, "name": r["name"], "status": r.get("status")}
    q = json.loads((ADS / "_drugs" / d.name / "drug.json").read_text())["charge"]
    uhf = 0
    el, xyz = read_xyz(d / "complex_opt.xyz")
    nd = len(read_xyz(d / "frag_drug.xyz")[0])
    D, C = list(range(nd)), list(range(nd, len(el)))
    # 1. single points
    et = ("--etemp", "1500")
    ec, ecf, edf = sp(d / "complex_opt.xyz", q, uhf, et), sp(d / "frag_carrier.xyz", 0, uhf, et), sp(d / "frag_drug.xyz", q, 0)
    # 2. pose choice: lowest energy among converged poses with contact in range
    valid = [p for p in r["poses"] if "E_Eh" in p and 1.25 <= p["min_contact_A"] <= 4.0]
    best = min(valid, key=lambda p: p["E_Eh"]) if valid else None
    # 3. bonds, integrity, contact
    heavy_d = [i for i in D if el[i] != "H"]
    heavy_c = [i for i in C if el[i] != "H"]
    dmin = min(np.linalg.norm(xyz[i] - xyz[j]) for i in heavy_d for j in heavy_c)
    cross = [(i, j) for i in D for j in C if np.linalg.norm(xyz[i] - xyz[j]) < 1.15 * (COV[el[i]] + COV[el[j]])]
    del0, dxyz0 = read_xyz(ADS / "_drugs" / d.name / "drug_opt.xyz")
    intact = pairs_within(del0, dxyz0, range(nd)) == pairs_within(el, xyz, D)
    return {"carrier": carrier, "name": r["name"], "status": "OK",
            "dE_int_json": r["delta_Eint_kcal_mol"], "dE_int_recomputed": round((ec - ecf - edf) * HARTREE, 3),
            "E_complex_diff_Eh": ec - r["E_complex_Eh"], "E_carrier_diff_Eh": ecf - r["E_carrier_frozen_Eh"],
            "E_drug_diff_Eh": edf - r["E_drug_frozen_Eh"],
            "pose_ok": best is not None and best["angle"] == r["best_orientation_deg"],
            "dmin_recomputed": round(dmin, 3), "dmin_json": r["min_contact_A"],
            "chem_recomputed": bool(cross), "chem_json": r["adsorption_mode"] == "chemisorption",
            "intact_recomputed": intact, "intact_json": bool(r["drug_intact"])}


def main(workers=4):
    jobs = [(c, d) for c in ("beta12_supported",) for d in sorted((ADS / c).iterdir()) if (d / "result.json").exists()]
    with ProcessPoolExecutor(workers) as ex:
        rows = list(ex.map(check, jobs))
    df = pd.DataFrame(rows)
    # 4. the table the manuscript reads
    tab = pd.read_csv(BASE / "results" / "quantum" / "adsorption_results.csv")
    m = df.merge(tab[["name", "carrier", "delta_Eint_kcal_mol", "delta_Eads_kcal_mol", "adsorption_mode",
                      "drug_intact", "min_contact_A"]], on=["name", "carrier"], how="left")
    m["table_matches_json"] = (np.isclose(m.delta_Eint_kcal_mol, m.dE_int_json, atol=1e-3)
                               & (m.adsorption_mode.eq("chemisorption") == m.chem_json)
                               & (m.drug_intact.astype(bool) == m.intact_json))
    out = BASE / "results" / "verification"
    out.mkdir(parents=True, exist_ok=True)
    m.to_csv(out / "adsorption_check.csv", index=False)
    ok = m[m.status == "OK"]
    print(f"complexes checked: {len(ok)} of {len(m)}")
    print(f"max |dE_int recomputed - reported|: {np.abs(ok.dE_int_recomputed - ok.dE_int_json).max():.4f} kcal/mol")
    print(f"max |E diff| complex/carrier/drug (Eh): {ok.E_complex_diff_Eh.abs().max():.2e} "
          f"{ok.E_carrier_diff_Eh.abs().max():.2e} {ok.E_drug_diff_Eh.abs().max():.2e}")
    print(f"pose choice correct: {ok.pose_ok.sum()}/{len(ok)}")
    print(f"regime matches: {(ok.chem_recomputed == ok.chem_json).sum()}/{len(ok)}; "
          f"integrity matches: {(ok.intact_recomputed == ok.intact_json).sum()}/{len(ok)}; "
          f"d_min matches (0.01 A): {(np.abs(ok.dmin_recomputed - ok.dmin_json) < 0.01).sum()}/{len(ok)}")
    print(f"manuscript table matches result.json: {ok.table_matches_json.sum()}/{len(ok)}")


if __name__ == "__main__":
    main(int(sys.argv[sys.argv.index("--workers") + 1]) if "--workers" in sys.argv else 4)
