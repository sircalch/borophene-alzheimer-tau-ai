"""
build_beta12_flake.py - flat, H-terminated beta12 borophene flake used as a
substrate-supported carrier model.

beta12 borophene is a triangular boron lattice with 1/6 of the sites vacant,
the vacancies forming rows along a (rectangular cell a = sqrt(3) d,
b = 3 d, 5 B per cell). Freestanding borophene is not stable; beta12 is grown on
Ag(111), which keeps it planar. A finite, unconstrained flake at the GFN2-xTB
level collapses into a compact boron cluster (the earlier B40H15 model did so
in every adsorption complex), so here the boron atoms are held at the ideal
planar geometry (harmonic position restraints, xtb $constrain, force
constant 5 Eh/bohr^2; xtb's exact $fix makes the optimiser diverge for this
system) and only the edge hydrogens are relaxed. The same boron
positions are kept fixed in every adsorption complex.

  d = b / 3 with the DFT lattice constants a = 2.93 A, b = 5.07 A
  flake: sites inside a 4 x 2 cell rectangle centred on the sheet, B atoms with
         fewer than 3 in-flake neighbours removed; every B that lost a
         neighbour relative to the infinite sheet capped with one H (1.19 A,
         in-plane, pointing away from its neighbours; caps closer than 1.6 A
         to another H dropped)

writes calculations/tau/beta12_flat/ (xtb run) and data/processed/carrier.json
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import adsorption_engine as ae  # noqa: E402

BASE = Path(__file__).resolve().parents[2]
WORK = BASE / "calculations" / "tau" / "beta12_flat"
A, B_ = 2.93, 5.07
D = B_ / 3
NB = 1.15 * 2 * 0.84             # B-B bond criterion (A)
NCELL = (4, 2)


def sheet(nx, ny):
    v1, v2 = np.array([0, D]), np.array([np.sqrt(3) * D / 2, D / 2])
    basis = [v1, 2 * v1, v2, v2 + v1, v2 + 2 * v1]          # (0,0) is the vacancy
    pts = []
    for i in range(-nx, nx + 1):
        for j in range(-ny, ny + 1):
            o = np.array([i * A, j * B_])
            pts += [o + b for b in basis]
    return np.array(pts)


def neighbours(p):
    dd = np.linalg.norm(p[:, None] - p[None], axis=2)
    return (dd > 0.1) & (dd < NB)


def main():
    full = sheet(8, 6)
    cn_full = neighbours(full).sum(1)
    c = np.array([NCELL[0] * A / 2, NCELL[1] * B_ / 2]) * 0
    hx, hy = NCELL[0] * A / 2 + 0.01, NCELL[1] * B_ / 2 + 0.01
    keep = (np.abs(full[:, 0] - c[0]) <= hx) & (np.abs(full[:, 1] - c[1]) <= hy)
    idx = np.where(keep)[0]
    while True:                                                # prune weakly bound edge atoms
        nb = neighbours(full[idx]).sum(1)
        if (nb >= 3).all():
            break
        idx = idx[nb >= 3]
    P = full[idx]
    nbm = neighbours(P)
    lost = cn_full[idx] - nbm.sum(1)
    H = []
    for k in np.where(lost >= 1)[0]:
        v = -(P[nbm[k]] - P[k]).sum(0)
        h = P[k] + 1.19 * v / np.linalg.norm(v)
        if all(np.linalg.norm(h - x) > 1.6 for x in H) and np.linalg.norm(P - h, axis=1).min() > 1.1:
            H.append(h)
    el = ["B"] * len(P) + ["H"] * len(H)
    xyz = np.hstack([np.vstack([P] + H if H else [P]), np.zeros((len(el), 1))])
    if (3 * len(P) + len(H)) % 2:                              # keep an even electron count
        el, xyz = el[:-1], xyz[:-1]
    xyz -= xyz.mean(0)
    WORK.mkdir(parents=True, exist_ok=True)
    for f in WORK.glob("*"):                                   # no stale xtb restart files
        if f.is_file():
            f.unlink()
    ae.write_xyz(WORK / "start.xyz", el, xyz, "flat beta12 flake")
    nbor = el.count("B")
    (WORK / "fix.inp").write_text("$constrain\n   force constant=5.0\n   atoms: 1-" + str(nbor) + "\n$end\n")
    res = {}
    for uhf in (0, 2):
        ns = f"u{uhf}"
        ae.run_xtb(["start.xyz", "--opt", "--gfn", "2", "--uhf", str(uhf), "--etemp", "1500",
                    "--input", "fix.inp", "--namespace", ns, "--iterations", "500"], WORK)
        out, _ = ae.run_xtb([f"{ns}.xtbopt.xyz", "--sp", "--gfn", "2", "--uhf", str(uhf), "--etemp", "1500",
                             "--namespace", f"{ns}sp"], WORK)
        (WORK / f"{ns}_sp.out").write_text(out, encoding="utf-8")
        res[uhf] = ae.parse(out)
    uhf = min(res, key=lambda u: res[u]["E"])
    best = res[uhf]
    final = WORK / f"u{uhf}.xtbopt.xyz"
    formula = f"B{nbor}H{el.count('H')}"
    rec = {"name": "beta12_supported", "formula": formula, "xyz": str(final.relative_to(BASE)).replace("\\", "/"),
           "E_Eh": best["E"], "uhf": uhf, "xtb_args": ["--etemp", "1500"], "fix_elements": ["B"],
           "HOMO_eV": best.get("HOMO"), "LUMO_eV": best.get("LUMO"),
           "E_singlet_Eh": res[0]["E"], "E_triplet_Eh": res[2]["E"],
           "lattice_A": {"a": A, "b": B_, "d_BB": round(D, 4)},
           "note": "planar beta12 flake, B fixed (substrate-supported model), H relaxed",
           "exclude_drugs": ["Thioflavin-S"]}
    (BASE / "data" / "processed" / "carrier.json").write_text(json.dumps(rec, indent=2))
    print(json.dumps(rec, indent=2))


if __name__ == "__main__":
    main()
