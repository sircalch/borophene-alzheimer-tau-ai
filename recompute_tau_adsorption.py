"""
recompute_tau_adsorption.py
===========================
Rebuild the Tau / beta-12 borophene adsorption dataset from scratch, correctly.

Why: the shipped `delta_Eint_SP_kcal_mol` column of
`data/processed/dataset_tau_borophene_pristine.csv` was a *single point* on an
UNRELAXED geometry where the drug was dropped ~3.2 A above the carrier's highest
edge atom - which for the buckled B40H15 flake leaves the drug 5-6 A off the
sheet (verified). No optimisation was ever run. Two drug structures were also
wrong (Tideglusib and Congo Red SMILES were truncated).

This script, for each of the 29 drugs:
  1. builds the drug from (corrected) SMILES, RDKit ETKDG + MMFF, then
     GFN2-xTB --opt  ->  E_drug (relaxed, isolated), drug_opt.xyz
  2. places it flat 3.2 A above the local sheet surface, tries z-rotations
     0/90/180/270 deg, GFN2-xTB --opt on each (chrg = formal charge)
  3. keeps converged poses whose closest drug-carrier contact is 1.3-4.0 A
     (bound, not clashing); picks the lowest-energy one
  4. on that pose: SP for E_complex, then SP on the carrier and drug fragments
     frozen at the complex geometry  ->  E_cf, E_df
  5. Delta_Eint_SP  = (E_complex - E_cf - E_df) * 627.509      (interaction)
     Delta_Eads     = (E_complex - E_carrier_iso - E_drug_iso) * 627.509  (adsorption)

Outputs (idempotent / resumable - re-run to fill gaps):
  calculations/tau_recompute/<Drug>/...            per-drug workspace
  calculations/tau_recompute/results.csv           machine-readable summary
  and, when --commit-datasets is passed and every drug is done:
  data/processed/dataset_tau_borophene_pristine.csv   (delta_Eint_SP + E_drug)
  data/processed/relaxed_adsorption_subset.csv        (all 29)
"""
import os, sys, json, time, shutil, subprocess, hashlib
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parent
PROC = BASE / "data" / "processed"
WORK = BASE / "calculations" / "tau_recompute"
WORK.mkdir(parents=True, exist_ok=True)

XTB = os.environ.get("XTB_EXE", "")
if not XTB:
    XTB = shutil.which("xtb") or shutil.which("xtb.exe") or ""
if not XTB:
    hits = list(BASE.glob("**/xtb-*/bin/xtb.exe")) or [Path("C:/Users/Andre/mm/xtb/Library/bin/xtb.exe")]
    XTB = str(hits[0])
XTB_ENV = dict(os.environ)
_share = Path(XTB).parent.parent / "share" / "xtb"
if _share.is_dir():
    XTB_ENV["XTBPATH"] = str(_share)
XTB_ENV.setdefault("OMP_NUM_THREADS", "4")

HARTREE = 627.509474

# --- corrected structures (the two that shipped truncated) -------------------
SMILES_FIX = {
    "Tideglusib": "O=C1SC(=O)N(Cc2ccccc2)N1c1cccc2ccccc12",   # NP-12, C19H14N2O2S
    "Congo Red":  "OS(=O)(=O)c1cc2ccc(N=Nc3ccc(-c4ccc(N=Nc5ccc6cc(S(=O)(=O)O)c(N)cc6c5)cc4)cc3)cc2c(N)c1",  # free acid, C32H24N6O6S2
}


# --------------------------------------------------------------------------- io
def read_xyz(p):
    L = Path(p).read_text().splitlines()
    n = int(L[0].split()[0])
    el = [x.split()[0] for x in L[2:2 + n]]
    xyz = np.array([[float(v) for v in x.split()[1:4]] for x in L[2:2 + n]])
    return el, xyz


def write_xyz(p, el, xyz, comment=""):
    with open(p, "w") as f:
        f.write(f"{len(el)}\n{comment}\n")
        for e, (x, y, z) in zip(el, xyz):
            f.write(f"{e:2s} {x:15.8f} {y:15.8f} {z:15.8f}\n")


def xtb_energy(out_text):
    e = None
    for l in out_text.splitlines():
        if "TOTAL ENERGY" in l:
            for tok in l.split():
                try:
                    e = float(tok); break
                except ValueError:
                    continue
    return e


def xtb_converged(out_text):
    return ("GEOMETRY OPTIMIZATION CONVERGED" in out_text
            and "FAILED TO CONVERGE GEOMETRY OPTIMIZATION" not in out_text)


def run_xtb(args, cwd, timeout=900):
    p = subprocess.run([XTB, *args], cwd=str(cwd), env=XTB_ENV,
                       capture_output=True, text=True, errors="replace",
                       timeout=timeout)
    return p.stdout + "\n" + p.stderr


# ------------------------------------------------------------------ chemistry
def drug_from_smiles(smiles, out_xyz):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise ValueError(f"bad SMILES: {smiles}")
    m = Chem.AddHs(m)
    cids = AllChem.EmbedMultipleConfs(m, numConfs=8, randomSeed=0xC0FFEE)
    if not cids:
        AllChem.EmbedMolecule(m, randomSeed=1)
        cids = [0]
    best, best_e = None, 1e9
    for c in cids:
        try:
            ff = AllChem.MMFFGetMoleculeForceField(
                m, AllChem.MMFFGetMoleculeProperties(m), confId=c)
            ff.Minimize(maxIts=2000)
            e = ff.CalcEnergy()
        except Exception:
            e = 0.0
        if e < best_e:
            best_e, best = e, c
    conf = m.GetConformer(best)
    el = [a.GetSymbol() for a in m.GetAtoms()]
    xyz = np.array([list(conf.GetAtomPosition(i)) for i in range(m.GetNumAtoms())])
    write_xyz(out_xyz, el, xyz, "rdkit ETKDG+MMFF")
    return el, xyz


def place(drug_el, drug_xyz, carr_el, carr_xyz, angle_deg, gap=3.2):
    d = drug_xyz.copy()
    d[:, :2] -= d[:, :2].mean(axis=0)
    th = np.radians(angle_deg)
    R = np.array([[np.cos(th), -np.sin(th), 0],
                  [np.sin(th), np.cos(th), 0], [0, 0, 1]])
    d = d @ R.T
    # local sheet height under the drug footprint
    cb = carr_xyz[[i for i, e in enumerate(carr_el) if e == "B"]]
    cb_xy_c = cb[:, :2].mean(axis=0)
    d[:, :2] += cb_xy_c
    near = carr_xyz[np.linalg.norm(carr_xyz[:, :2] - cb_xy_c, axis=1) < 7.0]
    z_local = near[:, 2].max() if len(near) else carr_xyz[:, 2].max()
    d[:, 2] += (z_local + gap) - d[:, 2].min()
    el = list(drug_el) + list(carr_el)
    xyz = np.vstack([d, carr_xyz])
    return el, xyz, len(drug_el)


def min_contact(el, xyz, n_drug):
    dd = xyz[:n_drug]
    cc = xyz[n_drug:]
    hv_d = [i for i, e in enumerate(el[:n_drug]) if e != "H"]
    hv_c = [i for i, e in enumerate(el[n_drug:]) if e != "H"]
    if not hv_d or not hv_c:
        return 99.0
    D = np.linalg.norm(dd[hv_d][:, None, :] - cc[hv_c][None, :, :], axis=2)
    return float(D.min())


# --------------------------------------------------------------------- worker
def process(args):
    name, smiles, q = args
    dslug = name.replace(" ", "_").replace("-", "_")
    wd = WORK / dslug
    wd.mkdir(exist_ok=True)
    res_p = wd / "result.json"
    if res_p.exists():
        try:
            return json.loads(res_p.read_text())
        except Exception:
            pass
    t0 = time.time()
    rec = {"name": name, "formal_charge": q, "smiles": smiles}
    try:
        # 1. drug
        draw = wd / "drug_raw.xyz"
        drug_from_smiles(smiles, draw)
        run_xtb([str(draw), "--opt", "--gfn", "2", "--chrg", str(q),
                 "--uhf", "0", "--iterations", "500", "--namespace", "dopt"],
                wd, timeout=900)
        dopt = wd / "dopt.xtbopt.xyz"
        if not dopt.exists():
            dopt = wd / "xtbopt.xyz"
        del_el, del_xyz = read_xyz(dopt)
        e_drug = xtb_energy(run_xtb([str(dopt), "--sp", "--gfn", "2", "--chrg",
                            str(q), "--uhf", "0", "--namespace", "dsp"], wd, 300))
        write_xyz(wd / "drug_opt.xyz", del_el, del_xyz, f"{name} GFN2 opt")
        rec["E_drug_Eh"] = e_drug
        rec["drug_formula"] = "".join(sorted(set(del_el)))
        rec["n_drug"] = len(del_el)

        carr_el, carr_xyz = read_xyz(BASE / "calculations" / "tau" / "beta12_carrier_optimized.xyz")
        e_carr = -56.192156  # GFN2 opt, matches beta12_carrier_optimized.xyz

        # 2-3. orientations
        best = None
        for ang in (0, 90, 180, 270):
            tag = f"o{ang}"
            el, xyz, nd = place(del_el, del_xyz, carr_el, carr_xyz, ang)
            cin = wd / f"cin_{ang}.xyz"
            write_xyz(cin, el, xyz, f"{name} {ang}deg")
            for f in wd.glob(f"{tag}.xtbopt.xyz"):
                f.unlink()
            txt = run_xtb([str(cin), "--opt", "--gfn", "2", "--chrg", str(q),
                           "--uhf", "0", "--iterations", "500", "--cycles", "500",
                           "--namespace", tag], wd, timeout=1200)
            cx = wd / f"{tag}.xtbopt.xyz"
            if not cx.exists() or not xtb_converged(txt):
                continue
            fel, fxyz = read_xyz(cx)
            e_c = xtb_energy(txt)
            mc = min_contact(fel, fxyz, nd)
            entry = {"ang": ang, "E": e_c, "contact": mc, "path": str(cx)}
            if 1.25 <= mc <= 4.0 and e_c is not None:
                if best is None or e_c < best["E"]:
                    best = entry
            elif best is None:
                best = {**entry, "unbound": True}

        if best is None or best.get("unbound") or best.get("E") is None:
            rec.update(status="NO_STABLE_ADSORPTION",
                       best=best, seconds=round(time.time() - t0, 1))
            res_p.write_text(json.dumps(rec, indent=2))
            return rec

        # 4. SP decomposition at complex geometry
        fel, fxyz = read_xyz(best["path"])
        nd = rec["n_drug"]
        write_xyz(wd / "complex_opt.xyz", fel, fxyz, f"{name}/B40H15 GFN2 opt")
        e_complex = xtb_energy(run_xtb([str(wd / "complex_opt.xyz"), "--sp",
                    "--gfn", "2", "--chrg", str(q), "--uhf", "0",
                    "--namespace", "csp"], wd, 300))
        write_xyz(wd / "frag_carrier.xyz", fel[nd:], fxyz[nd:], "carrier @complex")
        write_xyz(wd / "frag_drug.xyz", fel[:nd], fxyz[:nd], "drug @complex")
        e_cf = xtb_energy(run_xtb([str(wd / "frag_carrier.xyz"), "--sp", "--gfn",
                          "2", "--chrg", "0", "--uhf", "0", "--namespace", "cf"], wd, 300))
        e_df = xtb_energy(run_xtb([str(wd / "frag_drug.xyz"), "--sp", "--gfn", "2",
                          "--chrg", str(q), "--uhf", "0", "--namespace", "df"], wd, 300))

        d_int = (e_complex - e_cf - e_df) * HARTREE
        d_ads = (e_complex - e_carr - e_drug) * HARTREE
        rec.update(status="OK", best_orientation_deg=best["ang"],
                   min_contact_A=round(best["contact"], 3),
                   E_complex_Eh=e_complex, E_carrier_frozen_Eh=e_cf,
                   E_drug_frozen_Eh=e_df,
                   delta_Eint_SP_kcal_mol=round(d_int, 3),
                   delta_Eads_kcal_mol=round(d_ads, 3),
                   final_pose_file=str((wd / "complex_opt.xyz").relative_to(BASE)),
                   sha256=hashlib.sha256((wd / "complex_opt.xyz").read_bytes()).hexdigest(),
                   seconds=round(time.time() - t0, 1))
    except Exception as e:
        import traceback
        rec.update(status="ERROR", error=repr(e), tb=traceback.format_exc()[-1500:],
                   seconds=round(time.time() - t0, 1))
    res_p.write_text(json.dumps(rec, indent=2))
    return rec


def main():
    df = pd.read_csv(PROC / "dataset_tau_borophene_pristine.csv")
    jobs = []
    for _, r in df.iterrows():
        smi = SMILES_FIX.get(r["name"], r["smiles"])
        jobs.append((r["name"], smi, int(r["formal_charge"])))

    only = sys.argv[1:] and sys.argv[1] != "--commit-datasets"
    if only:
        want = set(sys.argv[1:])
        jobs = [j for j in jobs if j[0] in want]

    nproc = int(os.environ.get("TAU_NPROC", "5"))
    print(f"xtb: {XTB}\n{len(jobs)} drugs, {nproc} workers", flush=True)
    done = []
    with ProcessPoolExecutor(max_workers=nproc) as ex:
        futs = {ex.submit(process, j): j[0] for j in jobs}
        for fut in as_completed(futs):
            r = fut.result()
            done.append(r)
            print(f"  [{len(done):2d}/{len(jobs)}] {r['name']:<22} {r.get('status'):<20} "
                  f"dEint={r.get('delta_Eint_SP_kcal_mol','-')}  "
                  f"contact={r.get('min_contact_A','-')}  {r.get('seconds','?')}s", flush=True)

    allres = [json.loads((WORK / n.replace(' ', '_').replace('-', '_') / "result.json").read_text())
              for n in df["name"]
              if (WORK / n.replace(' ', '_').replace('-', '_') / "result.json").exists()]
    pd.DataFrame(allres).to_csv(WORK / "results.csv", index=False)
    print(f"\nwrote {WORK / 'results.csv'}  ({len(allres)}/{len(df)} drugs)", flush=True)

    if "--commit-datasets" in sys.argv:
        ok = {r["name"]: r for r in allres if r.get("status") == "OK"}
        if len(ok) < len(df):
            print(f"NOT committing: only {len(ok)}/{len(df)} drugs OK", flush=True)
            return
        df2 = df.copy()
        for i, r in df2.iterrows():
            o = ok[r["name"]]
            df2.at[i, "E_drug_Eh"] = o["E_drug_Eh"]
            df2.at[i, "delta_Eint_SP_kcal_mol"] = o["delta_Eint_SP_kcal_mol"]
            df2.at[i, "carrier_formula"] = "B40H15"
        df2.to_csv(PROC / "dataset_tau_borophene_pristine.csv", index=False)
        sub = pd.DataFrame([{
            "name": r["name"], "best_orientation_deg": ok[r["name"]]["best_orientation_deg"],
            "delta_Eint_SP_kcal_mol": ok[r["name"]]["delta_Eint_SP_kcal_mol"],
            "delta_Eads_kcal_mol": ok[r["name"]]["delta_Eads_kcal_mol"],
            "min_contact_A": ok[r["name"]]["min_contact_A"],
            "convergence_status": "CONVERGED",
            "final_pose_file": ok[r["name"]]["final_pose_file"],
            "sha256": ok[r["name"]]["sha256"],
        } for _, r in df2.iterrows()])
        sub.to_csv(PROC / "relaxed_adsorption_subset.csv", index=False)
        print("committed dataset_tau_borophene_pristine.csv + relaxed_adsorption_subset.csv", flush=True)


if __name__ == "__main__":
    main()
