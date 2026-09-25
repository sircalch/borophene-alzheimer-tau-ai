"""
carrier_reference.py - reference energy of the restrained beta12 carrier for dE_ads.

The flat, restrained B44H16 flake has several nearby minima (the sheet can
buckle by ~0.1 A within the restraints). The build geometry is one of them;
carriers extracted from physisorbed complexes relax, under the same restraints,
into lower ones. dE_ads must be referenced to the lowest restrained carrier, so
the carrier of every physisorbed complex is re-relaxed (drug removed, same
$constrain, GFN2-xTB, etemp 1500) and the lowest energy is kept.
dE_int (frozen fragments) does not depend on this reference.

updates data/processed/carrier.json (E_Eh; build energy kept as E_build_Eh) and
rewrites delta_Eads_kcal_mol in results/quantum/adsorption_results.csv
"""
import json
import shutil
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import adsorption_engine as ae  # noqa: E402

BASE = Path(__file__).resolve().parents[2]
CJ = BASE / "data" / "processed" / "carrier.json"
WORK = BASE / "calculations" / "tau" / "carrier_reference"


def main():
    c = json.loads(CJ.read_text())
    ads = pd.read_csv(BASE / "results" / "quantum" / "adsorption_results.csv")
    phys = ads[(ads.adsorption_mode == "physisorption") & ads.drug_intact.astype(bool)].name
    WORK.mkdir(parents=True, exist_ok=True)
    nb = sum(1 for l in (BASE / c["xyz"]).read_text().splitlines()[2:] if l.split() and l.split()[0] == "B")
    rows = []
    for n in phys:
        wd = WORK / n.replace(" ", "_")
        wd.mkdir(exist_ok=True)
        if not (wd / "sp.out").exists():
            shutil.copy(BASE / "calculations" / "adsorption" / c["name"] / n.replace(" ", "_") / "frag_carrier.xyz",
                        wd / "carrier.xyz")
            (wd / "fix.inp").write_text("$constrain\n   force constant=5.0\n   atoms: 1-" + str(nb) + "\n$end\n")
            ae.run_xtb(["carrier.xyz", "--opt", "--gfn", "2", "--uhf", str(c["uhf"]), "--etemp", "1500",
                        "--input", "fix.inp", "--iterations", "500", "--cycles", "800"], wd)
            out, _ = ae.run_xtb(["xtbopt.xyz", "--sp", "--gfn", "2", "--uhf", str(c["uhf"]), "--etemp", "1500",
                                 "--namespace", "sp"], wd)
            (wd / "sp.out").write_text(out, encoding="utf-8")
        e = ae.parse((wd / "sp.out").read_text(encoding="utf-8"))["E"]
        rows.append({"from_complex": n, "E_Eh": e})
        print(f"{n:<22} {e:.6f}", flush=True)
    r = pd.DataFrame(rows).sort_values("E_Eh")
    r.to_csv(WORK / "carrier_reference_scan.csv", index=False)
    c.setdefault("E_build_Eh", c["E_Eh"])
    c["E_Eh"] = float(r.E_Eh.iloc[0])
    c["E_ref_source"] = f"lowest restrained re-relaxation of the carrier from {len(r)} physisorbed complexes " \
                        f"(best: {r.from_complex.iloc[0]})"
    CJ.write_text(json.dumps(c, indent=2))
    ads = pd.read_csv(BASE / "results" / "quantum" / "adsorption_results.csv")
    iso = pd.read_csv(BASE / "results" / "quantum" / "isolated_drugs_qm_results.csv").set_index("name")
    ok = ads.status == "OK"
    ads.loc[ok, "delta_Eads_kcal_mol"] = ((ads.loc[ok, "E_complex_Eh"] - c["E_Eh"]
                                           - ads.loc[ok, "name"].map(iso.E_drug_Eh)) * ae.HARTREE).round(3)
    ads.to_csv(BASE / "results" / "quantum" / "adsorption_results.csv", index=False)
    for r in ads[ok].itertuples():             # keep the per-complex cache consistent
        f = BASE / "calculations" / "adsorption" / c["name"] / r.name.replace(" ", "_") / "result.json"
        j = json.loads(f.read_text())
        j["delta_Eads_kcal_mol"] = r.delta_Eads_kcal_mol
        f.write_text(json.dumps(j, indent=2))
    print("E_ref", c["E_Eh"], "build", c["E_build_Eh"], "diff kcal/mol", (c["E_build_Eh"] - c["E_Eh"]) * ae.HARTREE)


if __name__ == "__main__":
    main()
