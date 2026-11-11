# Spatial tumour CTMC

The event-driven direct stochastic simulation algorithm is in `CTMC_EventDrivenSimulation.m`. The drivers use the original locked values of proliferation, clearance and domain settings.

## Main results

1. After generating the fate-route and repair-time scan outputs in [LPO-FPT](../LPO-FPT/), run `CTMC_MesoscopicRiskInterface` in this directory. It computes the RMST-equivalent commitment hazards used by the CTMC.
2. Run `CTMC_CommitmentClearanceScan` to simulate the commitment–clearance landscape. It writes `Macro_Q1_TargetedRefinement_FINAL_OUT/Macro_Q1_TargetedRefinement_FINAL.mat`. Run `CTMC_CommitmentClearanceFigure` to plot extinction, growth phenotypes and transition boundaries. The [external consistency check](../External-Validation/) also uses this MAT file.
3. Run `CTMC_SpatialArrangementSimulation` after the commitment and clearance scan. Then run `CTMC_NeighbourhoodAccessMechanism` and `CTMC_SpatialArrangementConditionScan` to compute neighbourhood-access and cross-condition results. Run `CTMC_SpatialArrangementFigure` after the spatial-arrangement and neighbourhood-access runs.

## Robustness and sensitivity

| Entry points | Analysis and prerequisite |
| --- | --- |
| `CTMC_RiskInterfaceRobustness.m` | Interface window and time-scale sensitivity using the two LPO-FPT outputs |
| `CTMC_ParameterBenchmarks.m` | Proliferation and clearance benchmarks using the frozen locks |
| `CTMC_NumericalValidation.m` | SSA regression, Monte Carlo precision, domain sensitivity and supportive legacy time-step comparison |
| `CTMC_CommitmentClearanceValidation.m` | Statistical checks using the commitment and clearance scan output |
| `CTMC_SpatialArrangementRobustness.m`, `CTMC_NeighbourhoodAccessFigure.m` | Spatial-pattern robustness and neighbourhood access using the spatial arrangement and neighbourhood access outputs |
| `CTMC_SpatialArrangementGenerality.m` | Sensitivity across commitment conditions using the targeted condition map |

## Frozen inputs

| File | Use |
| --- | --- |
| `Macro_r_Calibration_OUT/Macro_r_LOCKED.mat` | Locked proliferation benchmark |
| `Macro_mu_Calibration_OUT/Macro_mu_LOCKED.mat` | Locked clearance benchmark |
| `Macro_Hazard_SSA_Domain_Confirm_OUT/Macro_SSA_LOCKED.mat` | Locked SSA domain settings |
| `Macro_Hazard_PreQ1_DT_Refinement_OUT/Macro_Hazard_PreQ1_DT_Refinement_FINAL.mat` | Legacy fixed-step comparator for the SSA numerical audit |
