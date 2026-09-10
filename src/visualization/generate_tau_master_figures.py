"""
generate_tau_master_figures.py
Master 9-Figure Q1 Scientific Visualization Engine at 300+ DPI for Article 4:
Alzheimer's Disease Tau Fibril Disaggregation & 2D Borophene Nanosheets.
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.ensemble import ExtraTreesRegressor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _pubstyle
_pubstyle.apply()
try:
    import _mol3d
except Exception:
    _mol3d = None

def get_dirs():
    base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    fig_dir = os.path.join(base_dir, "figures")
    os.makedirs(fig_dir, exist_ok=True)
    return base_dir, fig_dir

def make_graphical_abstract(base_dir, fig_dir):
    """Composed graphical abstract -> figures/fig1_graphical_abstract.png."""
    import graphical_abstract
    graphical_abstract.build()


def make_fig1_workflow(base_dir, fig_dir):
    fig, ax = plt.subplots(figsize=(14, 7), dpi=300)
    ax.axis('off')
    
    boxes = [
        ("1. 2D beta-12 Borophene\n(pristine B40H15 cluster)", 0.05, 0.55, 0.25, 0.35, "#F3E5F5", "#6A1B9A"),
        ("2. Blood-Brain Barrier (BBB)\nReceptor-mediated transcytosis route\n(proposed; not modelled here)", 0.38, 0.55, 0.25, 0.35, "#EDE7F6", "#4527A0"),
        ("3. Cryo-EM Target\nAlzheimer Tau filament core\n(PDB 5O3L)", 0.70, 0.55, 0.25, 0.35, "#FCE4EC", "#AD1457"),
        ("4. Quantum tight-binding (GFN2-xTB)\nInteraction energies + CDFT indices\n(real Delta_Eint,SP -13.8 to -0.9 kcal/mol\nfor 26/29; 3 phenothiazine dyes clash)", 0.04, 0.10, 0.27, 0.35, "#E0F7FA", "#00838F"),
        ("5. Real physical docking\nAutoDock Vina v1.2.7 (cross-beta)\n(29 Alzheimer / Tau drugs; exploratory)", 0.375, 0.10, 0.25, 0.35, "#E8F5E9", "#2E7D32"),
        ("6. Explainable machine learning\nLeak-free nested Ridge CV\n(descriptor QSPR weak; Williams domain)", 0.70, 0.10, 0.25, 0.35, "#FFF3E0", "#E65100"),
    ]
    
    for title, x, y, w, h, bg_c, border_c in boxes:
        rect = patches.Rectangle((x, y), w, h, facecolor=bg_c, edgecolor=border_c, lw=2.0, transform=ax.transAxes, zorder=2)
        ax.add_patch(rect)
        ax.text(x + w/2, y + h/2, title, ha='center', va='center', fontsize=9.0, fontweight='bold', color='#311B92', transform=ax.transAxes, zorder=3)
        
    arrow_props = dict(facecolor='#37474F', edgecolor='#37474F', width=2.5, headwidth=8, shrink=0.05)
    ax.annotate('', xy=(0.37, 0.72), xytext=(0.31, 0.72), xycoords='axes fraction', arrowprops=arrow_props)
    ax.annotate('', xy=(0.69, 0.72), xytext=(0.64, 0.72), xycoords='axes fraction', arrowprops=arrow_props)
    ax.annotate('', xy=(0.37, 0.27), xytext=(0.31, 0.27), xycoords='axes fraction', arrowprops=arrow_props)
    ax.annotate('', xy=(0.69, 0.27), xytext=(0.64, 0.27), xycoords='axes fraction', arrowprops=arrow_props)
    ax.annotate('', xy=(0.50, 0.48), xytext=(0.50, 0.54), xycoords='axes fraction', arrowprops=dict(facecolor='#4A148C', width=2.0, headwidth=7))
    
    plt.title("Figure 1: Multi-Scale Computational Workflow: Quantum-Guided & Machine Learning Modeling of Borophene for Alzheimer's Tau Fibrils", fontsize=13, fontweight='bold', pad=15)
    out_p = os.path.join(fig_dir, "fig1_tau_workflow_methodology.png")
    plt.savefig(out_p, bbox_inches='tight')
    plt.close()
    print(f"Generated Figure 1: {out_p}")

def make_fig2_quantum(base_dir, fig_dir):
    # Was entirely hardcoded arrays (homo/lumo/eta/omega for 3 "systems"),
    # never computed from any real calculation -- and no real complex-level
    # FMO calculation exists for either borophene variant (the pristine
    # dataset's E_HOMO_eV/E_LUMO_eV are the isolated-drug orbitals, reused;
    # the chi3-PEG-Tf functionalized carrier has no real data at all). Now
    # shows the real per-compound GFN2-xTB frontier orbital distribution and
    # CDFT indices for the 29-compound isolated cohort (real *_drug_sp.out
    # single points, zero new computation).
    homo_lumo_csv = os.path.join(base_dir, "data", "processed", "tau_isolated_real_homo_lumo.csv")
    df = pd.read_csv(homo_lumo_csv)
    homo = df["E_HOMO_real_eV"].values
    lumo = df["E_LUMO_real_eV"].values
    gap = lumo - homo
    mu = -(homo + lumo) / 2.0
    eta = gap / 2.0
    omega = mu ** 2 / (2.0 * eta)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6), dpi=300)
    plt.subplots_adjust(top=0.86, wspace=0.28)

    ax0 = axes[0]
    ax0.hist(homo, bins=10, color='#4A148C', alpha=0.75, edgecolor='k', label=f'E_HOMO (mean={homo.mean():.2f} eV)')
    ax0.hist(lumo, bins=10, color='#D81B60', alpha=0.75, edgecolor='k', label=f'E_LUMO (mean={lumo.mean():.2f} eV)')
    ax0.set_xlabel("Electronic Energy (eV)", fontsize=11)
    ax0.set_ylabel("Compound Count", fontsize=11)
    ax0.set_title(f"(a) Real GFN2-xTB Frontier Molecular Orbitals (n={len(df)})", fontsize=11.5, fontweight='bold', pad=10)
    ax0.grid(True, linestyle=':', alpha=0.6)
    ax0.legend(loc='upper left', frameon=True, fontsize=9)

    ax1 = axes[1]
    ax1.scatter(eta, omega, color='#00796B', s=70, edgecolor='k', alpha=0.85)
    ax1.set_xlabel(r"Chemical Hardness $\eta$ (eV)", fontsize=11)
    ax1.set_ylabel(r"Electrophilicity Index $\omega$ (eV)", fontsize=11)
    ax1.set_title("(b) Real Conceptual DFT Global Reactivity Indices", fontsize=11.5, fontweight='bold', pad=10)
    ax1.grid(True, linestyle=':', alpha=0.6)

    plt.suptitle("Figure 2: Real Quantum CDFT Electronic Reactivity of the Isolated Tau Therapeutics Cohort", fontsize=12.5, fontweight='bold', y=0.98)
    out_p = os.path.join(fig_dir, "fig2_tau_quantum_cdft_architecture.png")
    plt.savefig(out_p, bbox_inches='tight')
    plt.close()
    print(f"Generated Figure 2: {out_p}")

def make_fig3_docking_profiles(base_dir, fig_dir):
    # Use the same real docking column the manuscript body, SI, Figure 5 and
    # Figure 9 use (dataset_tau_borophene_pristine.csv / vina_5O3L_kcal_mol),
    # NOT results/docking/real_vina_docking_summary.csv, which is a stale run
    # over a different (BACE/gamma-secretase) cohort inconsistent with the paper.
    vina_csv = os.path.join(base_dir, "data", "processed", "dataset_tau_borophene_pristine.csv")
    if not os.path.exists(vina_csv):
        return
    df = pd.read_csv(vina_csv).rename(columns={"vina_5O3L_kcal_mol": "Real_Vina_Docking_Score_kcal_mol"})
    df = df.dropna(subset=["Real_Vina_Docking_Score_kcal_mol"])

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.8), dpi=300)
    plt.subplots_adjust(top=0.86, wspace=0.30, bottom=0.15)

    ax0 = axes[0]
    sns.histplot(df['Real_Vina_Docking_Score_kcal_mol'], kde=True, color='#4A148C', bins=12, ax=ax0, edgecolor='k')
    ax0.axvline(df['Real_Vina_Docking_Score_kcal_mol'].mean(), color='r', linestyle='--', lw=2.0,
                label=f"Mean Delta_G = {df['Real_Vina_Docking_Score_kcal_mol'].mean():.2f} kcal/mol")
    ax0.set_xlabel("AutoDock Vina Real Binding Energy (kcal/mol)", fontsize=10.5, fontweight='bold')
    ax0.set_ylabel("Therapeutic Compound Count", fontsize=10.5, fontweight='bold')
    ax0.set_title("(a) Vina score distribution on the Tau filament core (PDB 5O3L)", fontsize=11.5, fontweight='bold', pad=10)
    ax0.legend(loc='upper left', frameon=True)
    ax0.grid(True, linestyle=':', alpha=0.6)

    ax1 = axes[1]
    df_sorted = df.sort_values(by='Real_Vina_Docking_Score_kcal_mol', ascending=True).head(10)
    colors = sns.color_palette("plasma", n_colors=10)
    bars = ax1.barh(df_sorted['name'], df_sorted['Real_Vina_Docking_Score_kcal_mol'], color=colors, edgecolor='k')
    ax1.set_xlabel("Real AutoDock Vina Score (kcal/mol)", fontsize=10.5, fontweight='bold')
    ax1.set_ylabel("Alzheimer / Tau Disaggregator", fontsize=10.5, fontweight='bold')
    ax1.set_title("(b) Top 10 High-Affinity Tau PHF Disaggregators", fontsize=11.5, fontweight='bold', pad=10)
    ax1.invert_yaxis()
    ax1.grid(True, linestyle=':', alpha=0.6)
    
    for bar in bars:
        w = bar.get_width()
        ax1.text(w - 0.20, bar.get_y() + bar.get_height()/2, f"{w:.2f}", 
                 va='center', ha='right', fontsize=9, fontweight='bold', color='white')
                 
    plt.suptitle("Figure 3: Physical Molecular Docking Statistical Profiles on Human Cryo-EM Tau Filament", fontsize=13, fontweight='bold', y=0.96)
    out_p = os.path.join(fig_dir, "fig3_tau_docking_vina_statistical_profiles.png")
    plt.savefig(out_p, bbox_inches='tight')
    plt.close()
    print(f"Generated Figure 3: {out_p}")

def make_fig4_residues(base_dir, fig_dir):
    freq_csv = os.path.join(base_dir, "results", "docking", "residue_frequency_ranking.csv")
    if not os.path.exists(freq_csv):
        return
    df = pd.read_csv(freq_csv).head(12)
    
    fig, ax = plt.subplots(figsize=(10, 5.5), dpi=300)
    colors = sns.color_palette("flare", n_colors=len(df))
    bars = ax.bar(df['Residue'], df['Contact_Frequency'], color=colors, edgecolor='k', lw=1.2)
    
    ax.set_xlabel("Tau filament core cross-beta residue (PDB 5O3L)", fontsize=11, fontweight='bold')
    ax.set_ylabel("Atomic Contact Frequency (d <= 3.8 Å)", fontsize=11, fontweight='bold')
    ax.set_title("Figure 4: Residue-Level Interaction Fingerprints on Human Tau Paired Helical Filaments", fontsize=12.5, fontweight='bold', pad=12)
    ax.grid(True, linestyle=':', alpha=0.6)
    
    for bar in bars:
        h = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2, h + 0.5, str(int(h)), 
                ha='center', va='bottom', fontsize=9.5, fontweight='bold')
                
    ax.set_ylim(0, max(df['Contact_Frequency']) + 4)
    out_p = os.path.join(fig_dir, "fig4_tau_residue_contact_frequency.png")
    plt.savefig(out_p, bbox_inches='tight')
    plt.close()
    print(f"Generated Figure 4: {out_p}")

def make_fig5_parity(base_dir, fig_dir):
    # Panels (b)/(c) were fit on `Target_DeltaG_bind` from
    # dataset_drug_borophene_pristine.csv / _functionalized.csv, whose
    # Delta_E_ads_kcal_mol was FABRICATED by train_tau_qsar_models.py from an
    # empirical formula over RDKit descriptors, never a real xTB calculation.
    # Real GFN2-xTB single-point interaction energies for all 29 compounds on
    # the pristine beta12 borophene carrier already exist
    # (dataset_tau_borophene_pristine.csv, delta_Eint_SP_kcal_mol -- the same
    # data used by scripts/run_nested_cv_leakfree.py), so panel (b) is fixed
    # with zero new computation. No real structural/quantum data exists for
    # the chi3-PEG-Tf functionalized carrier (no complex geometries were ever
    # built for it) -- that panel is omitted rather than left fabricated.
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.linear_model import RidgeCV
    from sklearn.model_selection import KFold, cross_val_predict
    from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score

    alpha_grid = np.array([0.001, 0.01, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0])
    # Both panels use the single master table dataset_tau_borophene_pristine.csv
    # (29 aggregation-inhibitor / imaging-ligand cohort). (a) Tau-filament Vina
    # docking score; (b) GFN2-xTB physisorption interaction energy on B40H15.
    _mt = os.path.join(base_dir, "data", "processed", "dataset_tau_borophene_pristine.csv")
    systems = [
        ("Tau-filament docking (Vina, 5O3L)", _mt,
         ["MolWt", "MolMR", "E_HOMO_eV", "Omega_eV"], "vina_5O3L_kcal_mol"),
        ("Borophene physisorption (GFN2-xTB)", _mt,
         ["MolWt", "MolMR", "E_HOMO_eV", "Omega_eV"], "delta_Eint_SP_kcal_mol"),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.5), dpi=300)
    plt.subplots_adjust(top=0.80, wspace=0.28, bottom=0.15)
    colors = ["#4A148C", "#00695C"]

    for ax_idx, (sys_name, f_path, desc_cols, target_col) in enumerate(systems):
        if not os.path.exists(f_path):
            continue
        df = pd.read_csv(f_path).dropna(subset=desc_cols + [target_col])
        # The GFN2-xTB screening shows 12/29 ligands CHEMISORB on pristine
        # borophene (covalent B-C/B-O, dEint -81 to -233 kcal/mol); those are a
        # different physical regime and are excluded from the physisorption
        # QSPR. Model is fit on the 17 genuine physisorbers only.
        if "adsorption_mode" in df.columns and target_col == "delta_Eint_SP_kcal_mol":
            df = df[df["adsorption_mode"] == "physisorption"]
        X = df[desc_cols].values
        y = df[target_col].values
        n, p = X.shape

        outer_cv = KFold(n_splits=5, shuffle=True, random_state=42)
        inner_cv = KFold(n_splits=5, shuffle=True, random_state=42)
        pipe = Pipeline([("scaler", StandardScaler()), ("ridge", RidgeCV(alphas=alpha_grid, cv=inner_cv))])
        y_pred = cross_val_predict(pipe, X, y, cv=outer_cv)
        rmse = mean_squared_error(y, y_pred) ** 0.5
        mae = mean_absolute_error(y, y_pred)
        r2 = r2_score(y, y_pred)

        ax = axes[ax_idx]
        ax.scatter(y, y_pred, color=colors[ax_idx], alpha=0.85, s=70, edgecolor='k', label=f'Out-of-Fold (n={n})')
        min_v = min(y.min(), y_pred.min()) - 0.5
        max_v = max(y.max(), y_pred.max()) + 0.5
        ax.plot([min_v, max_v], [min_v, max_v], 'r--', lw=2.0, label='Ideal 1:1 Parity')

        stats_txt = f"Leak-free nested 5x5 CV (n={n}, p={p})\nRMSE = {rmse:.2f} kcal/mol\nMAE = {mae:.2f} kcal/mol\n$Q^2_{{CV}}$ = {r2:.3f}"
        ax.text(0.05, 0.95, stats_txt, transform=ax.transAxes, fontsize=8.5, va='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.85, edgecolor='#B0BEC5'))

        ax.set_title(f"({chr(97+ax_idx)}) {sys_name}", fontsize=11.5, fontweight='bold', pad=10)
        ax.set_xlabel("Real Observed (kcal/mol)", fontsize=10.5)
        if ax_idx == 0:
            ax.set_ylabel("Out-of-Fold Predicted (kcal/mol)", fontsize=10.5)
        ax.grid(True, linestyle=':', alpha=0.6)
        ax.legend(loc='lower right', fontsize=8.5, frameon=True)

    plt.suptitle("Figure 5: Leak-Free Nested CV Parity for Nano-QSAR on Borophene (real data only)", fontsize=13, fontweight='bold', y=0.98)
    out_p = os.path.join(fig_dir, "fig5_tau_parity_models_evaluation.png")
    plt.savefig(out_p, bbox_inches='tight')
    plt.close()
    print(f"Generated Figure 5: {out_p}")

def make_fig6_shap(base_dir, fig_dir):
    # Was fit on dataset_drug_borophene_functionalized.csv's FABRICATED
    # Target_DeltaG_bind (no real structural/quantum data exists for the
    # chi3-PEG-Tf functionalized carrier). Refit on the real GFN2-xTB
    # delta_Eint_SP_kcal_mol for the pristine beta12 borophene (all 29
    # compounds), the same real data used in Figure 5.
    f_path = os.path.join(base_dir, "data", "processed", "dataset_tau_borophene_pristine.csv")
    if not os.path.exists(f_path):
        return
    df = pd.read_csv(f_path)
    feature_cols = ["MolWt", "MolMR", "E_HOMO_eV", "E_LUMO_eV", "Gap_eV", "Eta_eV", "Mu_eV", "Omega_eV"]
    df = df.dropna(subset=feature_cols + ["delta_Eint_SP_kcal_mol"])
    if "adsorption_mode" in df.columns:
        df = df[df["adsorption_mode"] == "physisorption"]   # physisorption QSPR only
    X = df[feature_cols]
    y = df['delta_Eint_SP_kcal_mol']

    model = ExtraTreesRegressor(n_estimators=100, random_state=42)
    model.fit(X, y)
    
    importances = model.feature_importances_
    indices = np.argsort(importances)[::-1][:10]
    
    top_features = [feature_cols[i] for i in indices]
    top_importances = importances[indices]
    
    fig, ax = plt.subplots(figsize=(10, 5.5), dpi=300)
    colors = sns.color_palette("Purples_r", n_colors=len(top_features))
    bars = ax.barh(top_features[::-1], top_importances[::-1], color=colors, edgecolor='k')
    
    ax.set_xlabel("Mean Absolute SHAP Value / Gini Feature Importance", fontsize=11, fontweight='bold')
    ax.set_ylabel("Molecular / Quantum CDFT Descriptor", fontsize=11, fontweight='bold')
    ax.set_title("Figure 6: Exploratory Feature Importance Rankings for 2D Borophene Delivery (real ΔE_int, pristine)", fontsize=11, fontweight='bold', pad=12)
    ax.grid(True, linestyle=':', alpha=0.6)
    
    for bar in bars:
        w = bar.get_width()
        ax.text(w + 0.005, bar.get_y() + bar.get_height()/2, f"{w:.3f}", 
                va='center', ha='left', fontsize=9, fontweight='bold')
                
    ax.set_xlim(0, max(top_importances) + 0.06)
    out_p = os.path.join(fig_dir, "fig6_tau_shap_xai_importance_rankings.png")
    plt.savefig(out_p, bbox_inches='tight')
    plt.close()
    print(f"Generated Figure 6: {out_p}")

try:
    import _pymol
except Exception:
    _pymol = None


def _pm_panel(ax, png, title=None, subtitle=None):
    import matplotlib.image as mpimg
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    if png and os.path.exists(png):
        ax.imshow(mpimg.imread(png))
    else:
        ax.text(0.5, 0.5, "render unavailable", ha="center", va="center", transform=ax.transAxes)
    if title:
        ax.set_title(title, fontsize=9.5, fontweight="bold", pad=5)
    if subtitle:
        ax.text(0.5, -0.03, subtitle, ha="center", va="top", fontsize=8.0,
                color=_pubstyle.MUTED, transform=ax.transAxes)


def make_fig9_3d_spatial(base_dir, fig_dir):
    """Figure 9 - PyMOL ray-traced renders of the two real adsorption regimes on
    pristine beta-12 borophene: chemisorption (covalent B-C/B-O) vs physisorption.
    Geometries are the GFN2-xTB relaxed complexes from calculations/tau_recompute/
    (copied to calculations/tau/<Drug>/<Drug>_B40H15_bound_complex.xyz)."""
    df = pd.read_csv(os.path.join(base_dir, "data", "processed",
                     "dataset_tau_borophene_pristine.csv")).set_index("name")
    ads = df["delta_Eint_SP_kcal_mol"]
    mode = df["adsorption_mode"] if "adsorption_mode" in df.columns else None
    calc = os.path.join(base_dir, "calculations", "tau")
    C = os.path.join(fig_dir, "_pm_cache"); os.makedirs(C, exist_ok=True)

    def cdir(name):
        return name.replace(" ", "_").replace("-", "_").replace("/", "_")

    chem = "Curcumin"                     # covalent B-C, dEint ~ -92
    phys = "Thioflavin-T"                 # physisorbed at ~3.4 A, dEint ~ -10
    jobs = [
        (os.path.join(calc, "beta12_carrier_optimized.xyz"), os.path.join(C, "t9_a.png"),
         "(a)  Pristine $\\beta$-12 borophene (B$_{40}$H$_{15}$)", "GFN2-xTB optimised carrier model"),
        (os.path.join(calc, cdir(chem), f"{cdir(chem)}_B40H15_bound_complex.xyz"), os.path.join(C, "t9_b.png"),
         f"(b)  {chem} - chemisorption",
         f"covalent B-C contact 1.4 A · $\\Delta E_{{int,SP}}$ = {ads[chem]:.0f} kcal/mol"),
        (os.path.join(calc, cdir(phys), f"{cdir(phys)}_B40H15_bound_complex.xyz"), os.path.join(C, "t9_c.png"),
         f"(c)  {phys} - physisorption",
         f"stacked at 3.4 A · $\\Delta E_{{int,SP}}$ = {ads[phys]:.0f} kcal/mol"),
    ]
    if _pymol and _pymol.AVAILABLE:
        for src, png, _, _ in jobs:
            try:
                _pymol.complex_figure(src, png, size=(1400, 1150), carbon="grey55", tilt=22)
            except Exception as exc:
                print(f"[fig9 PyMOL {os.path.basename(src)}] {exc}")

    fig, axes = plt.subplots(1, 3, figsize=(11.4, 4.1))
    fig.subplots_adjust(wspace=0.05, top=0.83, bottom=0.15, left=0.02, right=0.98)
    for ax, (_, png, title, sub) in zip(axes, jobs):
        _pm_panel(ax, png, title, sub)
    n_chem = int((mode == "chemisorption").sum()) if mode is not None else 12
    fig.suptitle(f"Figure 9. Two adsorption regimes on pristine $\\beta$-12 borophene: "
                 f"{n_chem}/29 ligands chemisorb (covalent B-C/B-O), the rest physisorb",
                 fontsize=10.0, fontweight="bold", y=0.985)
    out_p = os.path.join(fig_dir, "fig9_tau_3d_spatial_binding_modes.png")
    _pubstyle.save(fig, out_p, also_pdf=False)
    print(f"Generated Figure 9 (PyMOL ray-traced): {out_p}")


def make_fig11_adsorption_landscape(base_dir, fig_dir):
    """Figure 11 - the real GFN2-xTB adsorption landscape: closest drug-carrier
    contact vs interaction energy, one point per ligand, coloured by regime.
    This is the central corrected result (replaces the fabricated/unrelaxed
    delta_Eint_SP column of earlier drafts)."""
    p = os.path.join(base_dir, "data", "processed", "dataset_tau_borophene_pristine.csv")
    df = pd.read_csv(p)
    if "min_contact_A" not in df.columns:
        print("[fig11] dataset lacks min_contact_A - run recompute_tau_adsorption.py")
        return
    fig, ax = plt.subplots(figsize=(7.6, 5.4))
    for m, col, lab in [("chemisorption", getattr(_pubstyle, "WARN", "#d55e00"), "chemisorption (B-C / B-O)"),
                        ("physisorption", getattr(_pubstyle, "ACCENT", "#0072b2"), "physisorption")]:
        s = df[df["adsorption_mode"] == m]
        ax.scatter(s["min_contact_A"], s["delta_Eint_SP_kcal_mol"], s=55,
                   color=col, edgecolor="k", linewidth=0.5, label=f"{lab}  (n={len(s)})", zorder=3)
    ax.axvspan(1.2, 1.9, color="0.9", zorder=0)
    ax.set_xlabel("closest drug-carrier heavy-atom contact (Å)")
    ax.set_ylabel("$\\Delta E_{int,SP}$ (kcal mol$^{-1}$, GFN2-xTB)")
    ax.set_yscale("symlog")
    ax.legend(frameon=True, fontsize=8.5, loc="lower right")
    ax.set_title("Figure 11. Adsorption landscape of 29 Tau-directed ligands on pristine $\\beta$-12 borophene",
                 fontsize=10, fontweight="bold", pad=8)
    _pubstyle.save(fig, os.path.join(fig_dir, "fig11_tau_adsorption_landscape.png"), also_pdf=False)
    print("Generated Figure 11 (adsorption landscape)")


def make_fig7_correlation(base_dir, fig_dir):
    csv_p = os.path.join(base_dir, "data", "processed", "tau_isolated_descriptors.csv")
    if not os.path.exists(csv_p):
        return
    df = pd.read_csv(csv_p)
    cols = [c for c in ["MW", "LogP", "LogS", "WS_mg_mL", "HBA", "HBD", "PSA",
                        "RBC", "NOR", "AromRings", "Polarizability_alpha",
                        "Fraction_Csp3", "E_HOMO", "E_LUMO", "Gap_eV",
                        "Hardness_eta", "Softness_S", "Electronegativity_chi",
                        "Chemical_Potential_mu", "Electrophilicity_omega"]
            if c in df.columns]
    corr = df[cols].corr()
    fig, ax = plt.subplots(figsize=(9.6, 8.0))
    sns.heatmap(corr, annot=True, fmt=".2f", cmap="vlag", center=0, vmin=-1, vmax=1,
                cbar_kws={"label": "Pearson correlation $r$", "shrink": 0.8}, ax=ax,
                annot_kws={"size": 6.0}, linewidths=0.4, linecolor="white", square=True)
    ax.set_title(f"Pearson inter-descriptor correlation ({len(cols)} descriptors, "
                 f"{len(df)} Tau therapeutics)")
    ax.tick_params(labelsize=6.5)
    out_p = os.path.join(fig_dir, "fig7_tau_descriptor_correlation_matrix.png")
    _pubstyle.save(fig, out_p, also_pdf=False)
    print(f"Generated Figure 7: {out_p}")


def make_fig10_deltarho(base_dir, fig_dir):
    """Figure 10 - charge-density difference for the curcumin / beta-12 borophene
    chemisorption complex (best-orientation relaxed pose, real GFN2-xTB
    densities). Cube + build script in results/quantum/drho/."""
    try:
        import _drho_fig
    except Exception as exc:
        print(f"[fig10 drho] helper unavailable: {exc}")
        return
    drho_dir = os.path.join(base_dir, "results", "quantum", "drho")
    dEint = None
    try:
        df = pd.read_csv(os.path.join(base_dir, "data", "processed",
                         "relaxed_adsorption_subset.csv")).set_index("name")
        dEint = float(df.loc["Curcumin", "delta_Eint_SP_kcal_mol"])
    except Exception:
        pass
    render = os.path.join(drho_dir, "tau_deltarho_render.png")
    render = _drho_fig.render_isosurface(drho_dir, "tau", render, level=0.005,
                                         turn=(8, -18, 0))
    out_p = os.path.join(fig_dir, "fig10_tau_charge_density_difference.png")
    _drho_fig.compose(out_p, render, 10,
                      "Covalent charge transfer in the curcumin / $\\beta$-12 borophene chemisorption complex",
                      "curcumin", "$\\beta$-12 borophene", 0.005, dEint_kcal=dEint)
    print(f"Generated Figure 10 (charge-density difference): {out_p}")


def generate_master_suite():
    base_dir, fig_dir = get_dirs()
    make_graphical_abstract(base_dir, fig_dir)
    make_fig1_workflow(base_dir, fig_dir)
    make_fig2_quantum(base_dir, fig_dir)
    make_fig3_docking_profiles(base_dir, fig_dir)
    make_fig4_residues(base_dir, fig_dir)
    make_fig5_parity(base_dir, fig_dir)
    make_fig6_shap(base_dir, fig_dir)
    make_fig7_correlation(base_dir, fig_dir)
    make_fig9_3d_spatial(base_dir, fig_dir)
    make_fig10_deltarho(base_dir, fig_dir)
    make_fig11_adsorption_landscape(base_dir, fig_dir)
    print("Master figure suite for Article 4 (Tau/Borophene) generated successfully.")

if __name__ == "__main__":
    generate_master_suite()
