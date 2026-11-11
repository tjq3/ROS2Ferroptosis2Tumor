# Spatial tumour CTMC

The event-driven direct stochastic simulation algorithm is in `macro_sim_core_hazard_SSA_FINAL.m`. The drivers use the original locked values of proliferation, clearance and domain settings.

## Main results

1. After generating the fate-route and repair-time scan outputs in [LPO-FPT](../LPO-FPT/), run `Meso_to_Macro_Risk_Interface_FINAL_v2` in this directory. It computes the RMST-equivalent commitment hazards used by the CTMC.
2. Run `Macro_Q1_TargetedRefinement_FINAL` to simulate the commitment–clearance landscape. It writes `Macro_Q1_TargetedRefinement_FINAL_OUT/Macro_Q1_TargetedRefinement_FINAL.mat`. Run `Fig5B_Q1_Main_REVISED_WITH_TRANSLATION_FINAL` to plot extinction, growth phenotypes and transition boundaries. The [external consistency check](../External-Validation/) also uses the Q1 MAT file.
3. Run `Macro_Q2_SpatialArrangement_FINAL` after Q1. Then run `Macro_Q2_NeighbourhoodAccess_Mechanism_FINAL` and `Macro_Q2_TargetedConditionMap_FINAL` to compute neighbourhood-access and cross-condition results. Run `Fig5C_Q2_REVISED_3D_MECHANISM_FINAL_v2` after the spatial-arrangement and neighbourhood-access runs.

## Robustness and sensitivity

| Entry points | Analysis and prerequisite |
| --- | --- |
| `SI_M1_MesoMacro_Interface_Robustness_FINAL.m` | Interface window and time-scale sensitivity using the two LPO-FPT outputs |
| `SI_M2_Macro_Parameter_Benchmarks_FINAL.m` | Proliferation and clearance benchmarks using the frozen locks |
| `SI_M3_Macro_SSA_Numerical_Validation_FINAL.m` | SSA regression, Monte Carlo precision, domain sensitivity and supportive legacy time-step comparison |
| `SI_M5_Q1_StatisticalValidation_FINAL.m` | Statistical checks using the Q1 output |
| `SI_Q2_1_SpatialArrangement_Robustness_FINAL.m`, `SI_Q2_2_NeighbourhoodAccess_Mechanism_FINAL.m` | Spatial-pattern robustness and neighbourhood access using the Q2 arrangement and mechanism outputs |
| `SI_Q2_3_ConditionGenerality_FINAL.m` | Sensitivity across commitment conditions using the targeted condition map |

## Frozen inputs

| File | Use |
| --- | --- |
| `Macro_r_Calibration_OUT/Macro_r_LOCKED.mat` | Locked proliferation benchmark |
| `Macro_mu_Calibration_OUT/Macro_mu_LOCKED.mat` | Locked clearance benchmark |
| `Macro_Hazard_SSA_Domain_Confirm_OUT/Macro_SSA_LOCKED.mat` | Locked SSA domain settings |
| `Macro_Hazard_PreQ1_DT_Refinement_OUT/Macro_Hazard_PreQ1_DT_Refinement_FINAL.mat` | Legacy fixed-step comparator for the SSA numerical audit |
