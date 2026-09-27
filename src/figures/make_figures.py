"""
make_figures.py - every figure of the Tau / beta12-borophene manuscript, drawn
from the pipeline outputs at Springer print size (see style.py).

usage: python src/figures/make_figures.py [fig ...]      (default: all)
writes figures/FigN.{pdf,png,tif}; 3D renders are cached in figures/_renders
"""
import json
import sys
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import figkit as K  # noqa: E402
import render3d as R  # noqa: E402
import style as S  # noqa: E402

BASE = HERE.parents[1]
FIG = BASE / "figures"
REN = FIG / "_renders"
POSES = BASE / "results" / "docking" / "real_poses"
CALC = BASE / "calculations"
CARRIER = json.loads((BASE / "data" / "processed" / "carrier.json").read_text())
CARRIER_XYZ = BASE / CARRIER["xyz"]
STACK = "G+I+K+M+O"          # protofilament-1 chains that line the GTP-1 site (chain K)

FAMILY = {"Dyes and imaging probes": S.GROUPS[0], "Polyphenols": S.GROUPS[1],
          "Aggregation / kinase modulators": S.GROUPS[2], "Symptomatic AD drugs": S.GROUPS[3]}
SHORT = {"Dyes and imaging probes": "Dyes /\nprobes", "Polyphenols": "Poly-\nphenols",
         "Aggregation / kinase modulators": "Modu-\nlators", "Symptomatic AD drugs": "Sympto-\nmatic"}
DYES = {"Methylene Blue", "Hydromethylthionine", "Azure A", "Toluidine Blue O", "Thioflavin-T", "Congo Red",
        "Chrysamine G", "FDDNP"}
MODS = {"Anle138b", "Tideglusib", "AZD1080", "Bexarotene"}
SYMPT = {"Donepezil", "Rivastigmine", "Galantamine", "Memantine", "Tacrine"}


def family(n):
    return ("Dyes and imaging probes" if n in DYES else "Aggregation / kinase modulators" if n in MODS
            else "Symptomatic AD drugs" if n in SYMPT else "Polyphenols")


def carrier_levels():
    """HOMO/LUMO of the restrained beta12 carrier (build_beta12_flake.py)."""
    return {"HOMO_eV": CARRIER["HOMO_eV"], "LUMO_eV": CARRIER["LUMO_eV"]}


def data():
    lib = pd.read_csv(BASE / "data" / "processed" / "compound_library_pubchem.csv")
    dock = pd.read_csv(BASE / "results" / "docking" / "real_vina_docking_summary.csv")
    ads = pd.read_csv(BASE / "results" / "quantum" / "adsorption_results.csv")
    iso = pd.read_csv(BASE / "results" / "quantum" / "isolated_drugs_qm_results.csv")
    m = lib[["name", "class", "smiles", "formal_charge", "n_heavy_atoms"]].merge(
        dock[["name", "vina_8FUG_kcal_mol", "ligand_efficiency"]], on="name").merge(
        ads[ads.status == "OK"], on="name").merge(iso.drop(columns=["charge"]), on="name")
    m["drug_intact"] = m.drug_intact.astype(bool)
    m["family"] = m.name.map(family)
    m["color"] = m.family.map(FAMILY)
    d = {"m": m, "carrier": carrier_levels()}
    d["redock"] = pd.read_csv(BASE / "results" / "docking" / "redocking_validation.csv")
    d["contacts"] = pd.read_csv(BASE / "results" / "docking" / "residue_contacts.csv")
    q = BASE / "results" / "qspr"
    d["q"] = json.loads((q / "dEint_summary.json").read_text())
    d["oof"] = pd.read_csv(q / "dEint_oof.csv")
    d["perm"] = pd.read_csv(q / "dEint_y_scrambling.csv")
    return d


def render(name, fn, *a, **kw):
    REN.mkdir(parents=True, exist_ok=True)
    png, meta = REN / f"{name}.png", REN / f"{name}.json"
    if not png.exists():
        out = fn(*a, out_png=str(png), **kw)
        meta.write_text(json.dumps(out if isinstance(out, dict) else {}))
        S.autocrop(png)
    return png, json.loads(meta.read_text()) if meta.exists() else {}


def family_legend(fig, y=0.0):
    K.legend_row(fig, [("dot", c, f) for f, c in FAMILY.items()], y=y, fontsize=6.2)


def complex_path(name):
    return CALC / "adsorption" / CARRIER["name"] / name.replace(" ", "_") / "complex_opt.xyz"


def carrier_render():
    return render("carrier", R.molecule, str(CARRIER_XYZ), tilt=0, size=(1300, 1100))


def pocket_render():
    x, s = REN / "redock_xtal.pdb", REN / "redock_smiles.pdb"
    REN.mkdir(parents=True, exist_ok=True)
    ctl = POSES / "tracer_stack_control"
    R.first_model_pdb(str(ctl / "Y9H_redock_out.pdbqt"), str(x))
    R.first_model_pdb(str(ctl / "Y9H_smiles_out.pdbqt"), str(s))
    return render("pocket", R.pocket_closeup, str(BASE / "data" / "raw" / "8FUG.pdb"),
                  [{"sel": "resn Y9H and chain K", "color": (0.70, 0.72, 0.75), "radius": 0.30, "name": "xtal"},
                   {"path": str(x), "color": (0.91, 0.55, 0.16), "radius": 0.16, "name": "rdx"},
                   {"path": str(s), "color": (0.20, 0.62, 0.60), "radius": 0.16, "name": "rds"}],
                  chain=STACK, cofactors="resn Y9H", size=(1600, 1300), slab=26, zoom=2.2)


# ------------------------------------------------------------------ Fig. 1
def fig1(d):
    """Workflow."""
    m, q = d["m"], d["q"]
    car, _ = carrier_render()
    pocket, _ = pocket_render()
    best = m.nsmallest(1, "vina_8FUG_kcal_mol").name.iloc[0]
    cplx, _ = render(f"cplx_{best}", R.molecule, str(complex_path(best)), tilt=70, size=(1300, 1100))
    mol2d = REN / "mol2d.png"
    if not mol2d.exists():
        from rdkit import Chem
        from rdkit.Chem import Draw
        Draw.MolToFile(Chem.MolFromSmiles(m.set_index("name").loc["Methylene Blue", "smiles"]), str(mol2d),
                       size=(700, 480))
        S.autocrop(mol2d)
    mini = REN / "qspr_mini.png"
    o = d["oof"]
    fm, am = plt.subplots(figsize=(1.2, 1.2))
    lo, hi = min(o.dE_int.min(), o.oof_pred.min()), max(o.dE_int.max(), o.oof_pred.max())
    am.plot([lo, hi], [lo, hi], color=S.MUTED, lw=0.8, ls=(0, (3, 2)))
    am.scatter(o.dE_int, o.oof_pred, s=9, color=S.ADS, edgecolor="white", lw=0.3)
    am.set_xticks([]); am.set_yticks([]); am.set_xlabel("observed", fontsize=6); am.set_ylabel("predicted", fontsize=6)
    fm.savefig(mini, dpi=400, bbox_inches="tight"); plt.close(fm)
    nchem = int((m.adsorption_mode == "chemisorption").sum())
    stages = [
        ("Drug set", [f"{len(m)} tau-directed and", "AD drugs, 4 families", "from PubChem"]),
        ("Fibril docking", ["PDB 8FUG", "PHF, GTP-1 site", "Vina 1.2.7", "2 redock controls"]),
        ("β$_{12}$ borophene", ["B$_{44}$H$_{16}$ planar flake", "B held at β$_{12}$ lattice", "(Ag-supported)"]),
        ("Adsorption", ["4 poses per drug", "Δ$E_{int}$, Δ$E_{ads}$", "bond check"]),
        ("QSPR", ["4 descriptors, ridge", "nested 5×5 CV", "Y-scrambling, AD"]),
    ]
    fig = plt.figure(figsize=(S.DOUBLE, 58 * S.MM))
    K.workflow(fig, stages, ("Outcome", [f"{len(m) - nchem} physisorbed", f"{nchem} chemisorbed",
                                         f"Δ$E_{{int}}$ $Q^2_{{CV}}$ = {q['Q2_CV']:.2f}"]),
               images={0: mol2d, 1: pocket, 2: car, 3: cplx, 4: mini}, accent=S.DOCK)
    S.save(fig, FIG, "Fig1")


# ------------------------------------------------------------------ Fig. 2
def fig2(d):
    """Docking controls in the GTP-1 site and contact residues."""
    png, _ = pocket_render()
    rd = d["redock"].set_index("control")
    ct = d["contacts"]
    n = ct.name.nunique()
    freq = ct.groupby("residue").name.nunique().sort_values(ascending=False).head(8)
    polar = ct[ct.polar].groupby("residue").name.nunique().reindex(freq.index).fillna(0)
    fig = plt.figure(figsize=(S.DOUBLE, 78 * S.MM))
    gs = GridSpec(1, 2, width_ratios=[1.35, 1], wspace=0.28, left=0.01, right=0.99, top=0.95, bottom=0.2)
    ax = fig.add_subplot(gs[0])
    ax.imshow(mpimg.imread(png))
    ax.set_axis_off()
    S.panel(ax, "a", x=0.02, y=0.97)
    kx = "tracer stack present: self-redock, crystal conformation"
    ks = "tracer stack present: production protocol, from SMILES"
    K.legend_row(fig, [("line", "#b3b8bf", "cryo-EM GTP-1 (Y9H), chain K"),
                       ("line", "#e88c29", f"redock, cryo-EM conf. ({rd.loc[kx, 'rmsd_heavy_atom_A']:.2f} Å)"),
                       ("line", "#339e99", f"from SMILES ({rd.loc[ks, 'rmsd_heavy_atom_A']:.2f} Å)")],
                   y=0.0, fontsize=6)
    ax2 = fig.add_subplot(gs[1])
    y = np.arange(len(freq))[::-1]
    ax2.barh(y, freq.values / n * 100, color=S.DOCK, alpha=0.25, height=0.7, label="any contact")
    ax2.barh(y, polar.values / n * 100, color=S.DOCK, height=0.7, label="polar contact")
    ax2.set_yticks(y)
    ax2.set_yticklabels(freq.index, fontsize=6.2)
    ax2.set_xlabel(f"Drugs in contact (% of {n})")
    ax2.set_xlim(0, 100)
    ax2.tick_params(axis="y", length=0)
    K.light_grid(ax2, "x")
    ax2.legend(loc="lower right", frameon=False)
    S.panel(ax2, "b", x=-0.22)
    S.save(fig, FIG, "Fig2")


# ------------------------------------------------------------------ Fig. 3
def fig3(d):
    """Docking scores by drug and docking vs adsorption."""
    m = d["m"]
    fig, axs = plt.subplots(1, 2, figsize=(S.DOUBLE, 92 * S.MM),
                            gridspec_kw=dict(width_ratios=[1, 1.15], wspace=0.35))
    fig.subplots_adjust(bottom=0.16)
    K.ranked_dots(axs[0], m.name.values, m.vina_8FUG_kcal_mol.values, m.color.values,
                  "Vina score (kcal mol$^{-1}$)")
    S.panel(axs[0], "a", x=-0.45)
    ax = axs[1]
    chem = (m.adsorption_mode == "chemisorption").values
    ax.scatter(m.vina_8FUG_kcal_mol[~chem], -m.delta_Eint_kcal_mol[~chem], c=m.color[~chem], s=22,
               edgecolor="white", lw=0.4, zorder=3)
    ax.scatter(m.vina_8FUG_kcal_mol[chem], -m.delta_Eint_kcal_mol[chem], facecolor="white",
               edgecolor=m.color[chem], s=22, lw=1.1, zorder=3)
    rho, p = spearmanr(m.vina_8FUG_kcal_mol, m.delta_Eint_kcal_mol)
    K.stat_box(ax, [f"Spearman ρ = {rho:.2f} (p = {p:.2f})", "open: chemisorbed"], loc="upper right")
    ax.set_xlabel("Vina score, 8FUG (kcal mol$^{-1}$)")
    ax.set_ylabel("−Δ$E_{int}$, β$_{12}$ flake (kcal mol$^{-1}$)")
    K.light_grid(ax)
    S.panel(ax, "b", x=-0.2)
    family_legend(fig, y=0.0)
    S.save(fig, FIG, "Fig3")


# ------------------------------------------------------------------ Fig. 4
def fig4(d):
    """Carrier model and frontier levels of drugs vs carrier."""
    m = d["m"].sort_values("Gap_eV").reset_index(drop=True)
    c = d["carrier"]
    car, _ = carrier_render()
    fig = plt.figure(figsize=(S.DOUBLE, 80 * S.MM))
    gs = GridSpec(1, 2, width_ratios=[0.8, 1.7], wspace=0.12, left=0.01, right=0.9, top=0.9, bottom=0.3)
    ax = fig.add_subplot(gs[0])
    K.render_panel(ax, car, sub="β$_{12}$ B$_{44}$H$_{16}$ (metallic, $E_F$ ≈ "
                   f"{(c['HOMO_eV'] + c['LUMO_eV']) / 2:.2f} eV)")
    S.panel(ax, "a", x=0.02, y=1.02)
    ax = fig.add_subplot(gs[1])
    x = np.arange(len(m))
    ax.vlines(x, m.E_HOMO_eV, m.E_LUMO_eV, color=m.color, lw=2.2, alpha=0.9)
    ax.scatter(x, m.E_HOMO_eV, s=7, color=m.color, zorder=3)
    ax.scatter(x, m.E_LUMO_eV, s=7, color=m.color, zorder=3, marker="s")
    ax.axhspan(c["HOMO_eV"], c["LUMO_eV"], color=S.PANEL_BG, zorder=0)
    for k, lab in (("HOMO_eV", "carrier\nHOMO"), ("LUMO_eV", "carrier\nLUMO")):
        ax.axhline(c[k], color=S.MUTED, lw=0.7, ls=(0, (4, 2)))
    ax.text(len(m) - 0.3, (c["HOMO_eV"] + c["LUMO_eV"]) / 2, " carrier\n HOMO/LUMO", ha="left", va="center",
            fontsize=6, color=S.MUTED, clip_on=False)
    ax.set_xlim(-0.7, len(m) - 0.3)
    ax.set_xticks(x)
    ax.set_xticklabels(m.name, rotation=90, fontsize=5.4)
    ax.set_ylabel("Orbital energy (eV)")
    ax.tick_params(axis="x", length=0)
    K.light_grid(ax, "y")
    S.panel(ax, "b", x=-0.09)
    family_legend(fig, y=0.93)
    S.save(fig, FIG, "Fig4")


# ------------------------------------------------------------------ Fig. 5
def fig5(d):
    """Adsorption landscape and energies by family."""
    m = d["m"]
    fig, axs = plt.subplots(1, 2, figsize=(S.DOUBLE, 72 * S.MM),
                            gridspec_kw=dict(width_ratios=[1.35, 1], wspace=0.3))
    chem = (m.adsorption_mode == "chemisorption").values
    K.landscape(axs[0], m.min_contact_A.values, m.delta_Eint_kcal_mol.values, chem, names=m.name.values,
                max_labels=5)
    bad = m[~m.drug_intact]
    if len(bad):
        axs[0].scatter(bad.min_contact_A, -bad.delta_Eint_kcal_mol, marker="x", color=S.INK, s=18, lw=0.8,
                       zorder=5)
    S.panel(axs[0], "a", x=-0.14)
    fams = list(FAMILY)
    K.group_strip(axs[1], m.family.values, -m.delta_Eint_kcal_mol.values, [FAMILY[f] for f in fams],
                  "−Δ$E_{int}$ (kcal mol$^{-1}$)", order=fams)
    axs[1].set_yscale("log")
    K.log_ticks(axs[1])
    axs[1].set_xticklabels([SHORT[f] for f in fams], fontsize=6)
    S.panel(axs[1], "b", x=-0.22)
    S.save(fig, FIG, "Fig5")


# ------------------------------------------------------------------ Fig. 6
def fig6(d):
    """Representative complexes: strongest chemisorbed, strongest physisorbed, reacted."""
    m = d["m"].set_index("name")
    ok = m[m.drug_intact]
    picks = [ok[ok.adsorption_mode == "chemisorption"].delta_Eint_kcal_mol.idxmin(),
             ok[ok.adsorption_mode == "physisorption"].delta_Eint_kcal_mol.idxmin()]
    picks += list(m[~m.drug_intact].index[:1])
    items = []
    for nm in picks:
        png, _ = render(f"cplx_{nm}", R.molecule, str(complex_path(nm)), tilt=70, size=(1300, 1100))
        r = m.loc[nm]
        mode = "reacted" if not r.drug_intact else r.adsorption_mode
        items.append((png, f"{nm} ({mode})", f"Δ$E_{{int}}$ = {r.delta_Eint_kcal_mol:.1f} kcal mol$^{{-1}}$"))
    fig, axs = plt.subplots(1, len(items), figsize=(S.DOUBLE, 62 * S.MM), gridspec_kw=dict(wspace=0.05))
    for ax, (png, t, sub), l in zip(np.atleast_1d(axs), items, "abcd"):
        K.render_panel(ax, png, title=t, sub=sub)
        S.panel(ax, l, x=0.02, y=1.02)
    S.save(fig, FIG, "Fig6")


# ------------------------------------------------------------------ Fig. 7
def fig7(d):
    """QSPR of dE_int."""
    q, oof, perm = d["q"], d["oof"], d["perm"]
    fig, axs = plt.subplots(1, 3, figsize=(S.DOUBLE, 62 * S.MM),
                            gridspec_kw=dict(wspace=0.55, width_ratios=[1, 1, 0.9]))
    K.parity(axs[0], oof.dE_int.values, oof.oof_pred.values, S.ADS, stats_lines=None,
             names=oof.name.values, label_extremes=2)
    K.stat_box(axs[0], [f"$n$ = {q['n']}, $p$ = {q['p']}", f"$Q^2_{{CV}}$ = {q['Q2_CV']:.2f}",
                        f"RMSE = {q['RMSE']:.1f}"], loc="lower right")
    K.williams(axs[1], oof.leverage.values, oof.std_residual.values, q["AD"]["h_star"], S.ADS,
               names=oof.name.values)
    K.scrambling(axs[2], perm.Q2_perm.values, q["Q2_CV"], S.ADS)
    fig.canvas.draw()
    top = max(ax.get_position().y1 for ax in axs)
    for ax, l in zip(axs, "abc"):
        fig.text(ax.get_position().x0 - 0.055, top + 0.02, l, fontsize=9, fontweight="bold", va="bottom")
    S.save(fig, FIG, "Fig7")


FIGS = {"1": fig1, "2": fig2, "3": fig3, "4": fig4, "5": fig5, "6": fig6, "7": fig7}

if __name__ == "__main__":
    S.apply()
    d = data()
    for key in (sys.argv[1:] or FIGS):
        FIGS[key](d)
        print(f"Fig{key} done")
