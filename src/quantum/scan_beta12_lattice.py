"""
scan_beta12_lattice.py - GFN2-xTB single points of the planar beta12 flake with
the boron lattice scaled uniformly (0.90-1.10 of d = 1.69 A), singlet and
triplet, hydrogens moved with their boron. Shows that the planar model is not
strained against the method (minimum at the reference lattice).

writes results/quantum/beta12_lattice_scan.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import adsorption_engine as ae  # noqa: E402

BASE = Path(__file__).resolve().parents[2]
W = BASE / "calculations" / "tau" / "beta12_lattice_scan"


def main():
    el, x = ae.read_xyz(BASE / "calculations" / "tau" / "beta12_flat" / "start.xyz")
    nb = el.count("B")
    c = x[:nb].mean(0)
    W.mkdir(parents=True, exist_ok=True)
    rows = []
    for s in (0.90, 0.95, 1.00, 1.05, 1.10):
        y = x.copy()
        y[:nb] = (x[:nb] - c) * s + c
        for i in range(nb, len(el)):
            j = np.argmin(np.linalg.norm(x[:nb] - x[i], axis=1))
            y[i] = y[j] + (x[i] - x[j])
        ae.write_xyz(W / f"s{s:.2f}.xyz", el, y)
        for u in (0, 2):
            out, _ = ae.run_xtb([f"s{s:.2f}.xyz", "--sp", "--gfn", "2", "--uhf", str(u), "--etemp", "1500",
                                 "--namespace", f"s{s:.2f}_u{u}"], W)
            rows.append({"scale": s, "d_BB_A": round(1.69 * s, 3), "uhf": u, "E_Eh": ae.parse(out)["E"]})
    df = pd.DataFrame(rows)
    df["dE_kcal_mol"] = (df.E_Eh - df.E_Eh.min()) * ae.HARTREE
    df.to_csv(BASE / "results" / "quantum" / "beta12_lattice_scan.csv", index=False)
    print(df.to_string())


if __name__ == "__main__":
    main()
