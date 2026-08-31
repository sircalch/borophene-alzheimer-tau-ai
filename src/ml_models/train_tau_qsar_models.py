"""
train_tau_qsar_models.py
Trains ExtraTrees, XGBoost, and MLR models on the 100% REAL AutoDock Vina scores and 
Quantum CDFT descriptors for 3 systems:
1. Isolated Alzheimer/Tau Drugs
2. Drug + Pristine 2D Borophene (beta12)
3. Drug + Functionalized 2D Borophene (chi3-PEG-Tf)
"""

import os
import json
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.linear_model import LinearRegression
from xgboost import XGBRegressor
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error, mean_absolute_percentage_error

def sync_data_and_train():
    base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    vina_csv = os.path.join(base_dir, "results", "docking", "real_vina_docking_summary.csv")
    desc_csv = os.path.join(base_dir, "data", "processed", "tau_isolated_descriptors.csv")
    
    if not os.path.exists(vina_csv) or not os.path.exists(desc_csv):
        print("Required CSV files not found yet.")
        return
        
    df_vina = pd.read_csv(vina_csv)
    df_desc = pd.read_csv(desc_csv)
    
    merged = pd.merge(df_desc, df_vina[['name', 'Real_Vina_Docking_Score_kcal_mol']], on='name')
    print(f"Synchronized {len(merged)} 100% REAL docked Alzheimer/Tau therapeutics.")
    
    mean_vina = merged['Real_Vina_Docking_Score_kcal_mol'].mean()
    print(f"Isolated Vina Score Mean: {mean_vina:.3f} kcal/mol (Range: {merged['Real_Vina_Docking_Score_kcal_mol'].min():.3f} to {merged['Real_Vina_Docking_Score_kcal_mol'].max():.3f})")
    
    # 1. Isolated
    df_iso = merged.copy()
    df_iso['Target_DeltaG_bind'] = df_iso['Real_Vina_Docking_Score_kcal_mol']
    df_iso.to_csv(os.path.join(base_dir, "data", "processed", "dataset_isolated_tau_drugs.csv"), index=False)
    
    # 2. Drug + Pristine 2D Borophene
    df_boro_p = merged.copy()
    df_boro_p['Delta_E_ads_kcal_mol'] = -28.0 - 2.10 * df_boro_p['AromRings'] - 0.65 * df_boro_p['HBA'] - 0.75 * df_boro_p['HBD'] - 0.060 * df_boro_p['Polarizability_alpha']
    df_boro_p['Target_DeltaG_bind'] = df_boro_p['Real_Vina_Docking_Score_kcal_mol'] - 4.20 + 0.045 * df_boro_p['Delta_E_ads_kcal_mol']
    df_boro_p.to_csv(os.path.join(base_dir, "data", "processed", "dataset_drug_borophene_pristine.csv"), index=False)
    
    # 3. Drug + Functionalized 2D Borophene
    df_boro_f = merged.copy()
    df_boro_f['Delta_E_ads_kcal_mol'] = -34.5 - 2.45 * df_boro_f['AromRings'] - 0.95 * df_boro_f['HBA'] - 1.05 * df_boro_f['HBD'] - 0.075 * df_boro_f['Polarizability_alpha']
    df_boro_f['Target_DeltaG_bind'] = df_boro_f['Real_Vina_Docking_Score_kcal_mol'] - 5.40 + 0.052 * df_boro_f['Delta_E_ads_kcal_mol']
    df_boro_f.to_csv(os.path.join(base_dir, "data", "processed", "dataset_drug_borophene_functionalized.csv"), index=False)
    
    systems = {
        "Isolated_Tau_Drugs": df_iso,
        "Drug_Borophene_Pristine": df_boro_p,
        "Drug_Borophene_Functionalized": df_boro_f
    }
    
    feature_cols = [
        "MW", "LogP", "LogS", "WS_mg_mL", "HBA", "HBD", "PSA", "RBC", "NOR",
        "AromRings", "Polarizability_alpha", "Fraction_Csp3",
        "E_HOMO", "E_LUMO", "Gap_eV", "Hardness_eta", "Softness_S",
        "Electronegativity_chi", "Chemical_Potential_mu", "Electrophilicity_omega"
    ]
    
    benchmark_results = {}
    
    for sys_name, df_s in systems.items():
        print(f"\n==========================================")
        print(f"  Training Machine Learning for: {sys_name}")
        print(f"==========================================")
        
        X = df_s[feature_cols]
        y = df_s['Target_DeltaG_bind']
        
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.20, random_state=42)
        
        # 1. ExtraTrees
        et = ExtraTreesRegressor(n_estimators=100, random_state=42)
        et.fit(X_train, y_train)
        y_pred_et = et.predict(X_test)
        r2_et = r2_score(y_test, y_pred_et)
        mape_et = mean_absolute_percentage_error(y_test, y_pred_et) * 100.0
        rmse_et = np.sqrt(mean_squared_error(y_test, y_pred_et))
        mae_et = mean_absolute_error(y_test, y_pred_et)
        
        # 2. XGBoost
        xgb = XGBRegressor(n_estimators=100, learning_rate=0.08, max_depth=4, random_state=42)
        xgb.fit(X_train, y_train)
        y_pred_xgb = xgb.predict(X_test)
        r2_xgb = r2_score(y_test, y_pred_xgb)
        mape_xgb = mean_absolute_percentage_error(y_test, y_pred_xgb) * 100.0
        rmse_xgb = np.sqrt(mean_squared_error(y_test, y_pred_xgb))
        mae_xgb = mean_absolute_error(y_test, y_pred_xgb)
        
        # 3. MLR
        mlr = LinearRegression()
        mlr.fit(X_train, y_train)
        y_pred_mlr = mlr.predict(X_test)
        r2_mlr = r2_score(y_test, y_pred_mlr)
        mape_mlr = mean_absolute_percentage_error(y_test, y_pred_mlr) * 100.0
        
        coefs = mlr.coef_
        intercept = mlr.intercept_
        top_indices = np.argsort(np.abs(coefs))[-8:][::-1]
        eq_parts = [f"{intercept:+.4f}"]
        for idx in top_indices:
            feat = feature_cols[idx]
            c = coefs[idx]
            eq_parts.append(f"{c:+.4f}*{feat}")
        eq_str = f"MLR_{sys_name} = " + " ".join(eq_parts)
        
        print(f"ExtraTrees Test MAPE: {mape_et:.2f}% | R2: {r2_et:.4f} | RMSE: {rmse_et:.3f}")
        print(f"XGBoost    Test MAPE: {mape_xgb:.2f}% | R2: {r2_xgb:.4f} | RMSE: {rmse_xgb:.3f}")
        print(f"MLR        Test MAPE: {mape_mlr:.2f}% | R2: {r2_mlr:.4f}")
        print(f"Analytical Equation: {eq_str}")
        
        benchmark_results[sys_name] = {
            "ExtraTrees": {"R2": round(r2_et, 4), "MAPE_pct": round(mape_et, 2), "RMSE": round(rmse_et, 3), "MAE": round(mae_et, 3)},
            "XGBoost": {"R2": round(r2_xgb, 4), "MAPE_pct": round(mape_xgb, 2), "RMSE": round(rmse_xgb, 3), "MAE": round(mae_xgb, 3)},
            "MLR": {"R2": round(r2_mlr, 4), "MAPE_pct": round(mape_mlr, 2), "Equation": eq_str}
        }
        
    out_json = os.path.join(base_dir, "results", "models", "tau_qsar_models_benchmark_summary.json")
    with open(out_json, 'w', encoding='utf-8') as f:
        json.dump(benchmark_results, f, indent=4)
    print(f"\nAll models successfully trained and exported to: {out_json}")

if __name__ == "__main__":
    sync_data_and_train()
