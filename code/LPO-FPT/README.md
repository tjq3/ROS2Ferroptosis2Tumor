# Mesoscopic LPO–FPT

1. Generate `TDCME_Fig2ABC_results_final.mat` in [TD-CME](../TD-CME/) and `Fig3BC_predictorIndependent_data.mat` in [SD-RDME](../SD-RDME/).
2. Run `LPO_FateRoutes` in this directory. It generates the microscopic input histories used for the LPO decoder and writes `Fig4B_fate_routes_final/Fig4B_fateRoutes_data.mat`.
3. Run `Fig4D_CellStateDependentFateRewiring_FINAL2` for the repair-time scan. It writes the MAT file in `Fig4D_FINAL_outputs/` used by the [CTMC interface](../CTMC/).

| Entry point | Analysis and prerequisite |
| --- | --- |
| `LPO_FateRoutes.m` | LPO dynamics, first-passage times and fate-route classification |
| `LPO_JointMechanism.m` | Joint-required mechanism using fate-route output |
| `Fig4D_CellStateDependentFateRewiring_FINAL2.m` | Repair-time scan using fate-route output |
| `LPO_JointMechanismRobustness.m` | Timing and threshold diagnostics using joint-mechanism and fate-route outputs |
| `LPO_ThresholdAndWeightRobustness.m` | LPO threshold and temporal/spatial weight sensitivity using fate-route output |
| `LPO_SpatialCoarseGraining.m` | Spatial-window and pooling sensitivity using fate-route and SD-RDME outputs |
| `LPO_TemporalMemoryAndHorizon.m` | Temporal-window, lipid-memory and observation-horizon sensitivity using fate-route output |
