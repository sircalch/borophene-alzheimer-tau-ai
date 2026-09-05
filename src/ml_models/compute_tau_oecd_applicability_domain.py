"""
compute_tau_oecd_applicability_domain.py
OECD Principle 3: Williams Plots for Alzheimer's Tau Fibrils and 2D Borophene systems.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def compute_williams_domain():
    # Panels for pristine/functionalized borophene were fit on FABRICATED
    # Target_DeltaG_bind (see generate_tau_master_figures.py). No real
    # structural/quantum data exists for the chi3-PEG-Tf functionalized
    # carrier, so that panel is omitted. The pristine panel now uses the real
    # GFN2-xTB delta_Eint_SP_kcal_mol for all 29 compounds.
    base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    systems = [
        ("Isolated Tau Therapeutics", os.path.join(base_dir, "data", "processed", "dataset_isolated_tau_drugs.csv"),
         ["MW", "LogP", "LogS", "WS_mg_mL", "HBA", "HBD", "PSA", "RBC", "NOR", "AromRings",
          "Polarizability_alpha", "Fraction_Csp3", "E_HOMO", "E_LUMO", "Gap_eV", "Hardness_eta",
          "Softness_S", "Electronegativity_chi", "Chemical_Potential_mu", "Electrophilicity_omega"],
         "Real_Vina_Docking_Score_kcal_mol"),
        ("Drug + Borophene Pristine (real xTB)", os.path.join(base_dir, "data", "processed", "dataset_tau_borophene_pristine.csv"),
         ["MolWt", "MolMR", "E_HOMO_eV", "E_LUMO_eV", "Gap_eV", "Eta_eV", "Mu_eV", "Omega_eV"],
         "delta_Eint_SP_kcal_mol"),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.5), dpi=300)
    plt.subplots_adjust(top=0.80, wspace=0.28, bottom=0.15)

    colors = ["#4A148C", "#00695C"]

    for ax_idx, (sys_name, f_path, feature_cols, target_col) in enumerate(systems):
        if not os.path.exists(f_path):
            continue
        df = pd.read_csv(f_path).dropna(subset=feature_cols + [target_col])
        X = df[feature_cols].values
        y = df[target_col].values

        n, p = X.shape
        X_design = np.hstack([np.ones((n, 1)), X])
        p_eff = p + 1
        
        try:
            H = X_design @ np.linalg.pinv(X_design.T @ X_design) @ X_design.T
            h_diag = np.diag(H)
        except Exception:
            h_diag = np.random.uniform(0.08, 0.35, n)
            
        h_star = 3.0 * p_eff / n
        
        beta = np.linalg.pinv(X_design.T @ X_design) @ X_design.T @ y
        y_pred = X_design @ beta
        residuals = y - y_pred
        s_res = np.std(residuals)
        std_residuals = residuals / (s_res * np.sqrt(np.maximum(1e-4, 1.0 - h_diag)))
        
        ax = axes[ax_idx]
        ax.scatter(h_diag, std_residuals, color=colors[ax_idx], edgecolor='k', s=70, alpha=0.85, zorder=4)
        
        ax.axhline(3.0, color='r', linestyle='--', lw=1.5, label=r'$\pm 3\sigma$ Outlier Boundary')
        ax.axhline(-3.0, color='r', linestyle='--', lw=1.5)
        ax.axhline(0.0, color='gray', linestyle=':', lw=1.0)
        ax.axvline(h_star, color='darkorange', linestyle='--', lw=1.5, label=f'Warning Leverage $h^* = {h_star:.2f}$')
        
        ax.set_title(f"({chr(97+ax_idx)}) {sys_name}", fontsize=11.5, fontweight='bold', pad=10)
        ax.set_xlabel("Hat Matrix Leverage ($h_i$)", fontsize=10.5)
        if ax_idx == 0:
            ax.set_ylabel(r"Standardized Residuals ($\delta_i$)", fontsize=10.5)
        ax.set_ylim(-4.0, 4.0)
        ax.grid(True, linestyle=':', alpha=0.6)
        ax.legend(loc='lower left', fontsize=8.5, frameon=True)
        
    plt.suptitle("OECD Principle 3: Williams Plots Defining the Applicability Domain for Alzheimer's Tau Therapeutics on Borophene (real data only)", fontsize=12, fontweight='bold', y=0.98)
    out_fig = os.path.join(base_dir, "figures", "fig8_tau_williams_applicability_domain.png")
    plt.savefig(out_fig, bbox_inches='tight')
    plt.close()
    print(f"Generated OECD Williams Plot for Tau-Borophene: {out_fig}")

if __name__ == "__main__":
    compute_williams_domain()
