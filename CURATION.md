# Selection and verification

## Selection

The release contains 36 MATLAB files across the five model directories. The original TD-CME archive contributes five files, SD-RDME contributes eight, LPO-FPT contributes seven, and the macro archive contributes ten. The external validation script contributes one. Five CTMC dependencies were recovered from the accompanying production directory.

The excluded archive entries are `LQSA6.m`, a standalone stationary-approximation illustration, `pathvalid5.m`, a standalone reduced path-integral illustration, and `SDRDME_cristalumen5chinese.m`, a schematic-only geometry drawing. Generated images, exploratory versions and duplicate archived results are outside the code selection.

Additional theory-only or secondary display scripts were excluded to focus the release.

- `TDCME_Stage6_ModelHierarchyAblation_FINAL.m`. Stationary, reduced-dynamic and full-state hierarchy comparison supporting the theory.
- `SDRDME_meancov2.m`. Spatial-average and covariance illustration supporting the theory.
- `cellstate_SI.m`. Companion diagnostic figures from the already retained cell-state scan.
- `Fig5A_MesoMacro_Interface_FINAL.m`. Interface visualization from the already retained interface generator.
- `SI_M4_Q1_MechanisticSupport_FINAL.m`. Additional Q1 mechanistic plots from the already retained Q1 output.

## Source mapping

Paths in the source column are relative to the supplied archive's `finalfinal/` directory. A `macro/v2/` prefix identifies recovered production dependencies. The validation source is `validation/finalfinal/validation.m`.

| Released file | Source |
| --- | --- |
| [`TDCME_fixedC2_ablation_FINAL.m`](code/TD-CME/TDCME_fixedC2_ablation_FINAL.m) | `TDCME_fixedC2_ablation_FINAL.m` |
| [`TDCME_Fig2ABC.m`](code/TD-CME/TDCME_Fig2ABC.m) | `TDCME_matchtohighorderandSI2.m` |
| [`TDCME_RegimeMap.m`](code/TD-CME/TDCME_RegimeMap.m) | `regimemap.m` |
| [`TDCME_CoreRobustness.m`](code/TD-CME/TDCME_CoreRobustness.m) | `robustness.m` |
| [`TDCME_TimeStepConvergence.m`](code/TD-CME/TDCME_TimeStepConvergence.m) | `robustness3.m` |
| [`SDRDME_RedoxBufferScan.m`](code/SD-RDME/SDRDME_RedoxBufferScan.m) | `SDRDME_Fe3Fe2redox.m` |
| [`SDRDME_RedoxBufferFigure.m`](code/SD-RDME/SDRDME_RedoxBufferFigure.m) | `SDRDME_Fe3Fe2redoxhighlight.m` |
| [`SDRDME_Counterfactuals.m`](code/SD-RDME/SDRDME_Counterfactuals.m) | `SDRDME_UHBFandSI.m` |
| [`SDRDME_LocalBurdenAndFPT.m`](code/SD-RDME/SDRDME_LocalBurdenAndFPT.m) | `SDRDME_burdenandFPT.m` |
| [`SDRDME_HeldOutValidation.m`](code/SD-RDME/SDRDME_HeldOutValidation.m) | `SDRDME_c1.m` |
| [`SDRDME_DescriptorScreening.m`](code/SD-RDME/SDRDME_DescriptorScreening.m) | `SDRDME_c1BDXHBandSI.m` |
| [`SDRDME_GeometryRobustness.m`](code/SD-RDME/SDRDME_GeometryRobustness.m) | `SDRDME_geometry.m` |
| [`SDRDME_CoreRobustness.m`](code/SD-RDME/SDRDME_CoreRobustness.m) | `SDRDME_robustness.m` |
| [`LPO_ThresholdAndWeightRobustness.m`](code/LPO-FPT/LPO_ThresholdAndWeightRobustness.m) | `MLWTWS_SIfinal.m` |
| [`Fig4D_CellStateDependentFateRewiring_FINAL2.m`](code/LPO-FPT/Fig4D_CellStateDependentFateRewiring_FINAL2.m) | `cellstate.m` |
| [`LPO_FateRoutes.m`](code/LPO-FPT/LPO_FateRoutes.m) | `fateroute.m` |
| [`LPO_JointMechanism.m`](code/LPO-FPT/LPO_JointMechanism.m) | `jointmechanism_final.m` |
| [`LPO_JointMechanismRobustness.m`](code/LPO-FPT/LPO_JointMechanismRobustness.m) | `jointmechanism_finalSI.m` |
| [`LPO_SpatialCoarseGraining.m`](code/LPO-FPT/LPO_SpatialCoarseGraining.m) | `spatialcoarseSI.m` |
| [`LPO_TemporalMemoryAndHorizon.m`](code/LPO-FPT/LPO_TemporalMemoryAndHorizon.m) | `temporalmemorySI.m` |
| [`Macro_Q1_TargetedRefinement_FINAL.m`](code/CTMC/Macro_Q1_TargetedRefinement_FINAL.m) | `Fig5B_Q1_mainphenotype_calu.m` |
| [`Fig5B_Q1_Main_REVISED_WITH_TRANSLATION_FINAL.m`](code/CTMC/Fig5B_Q1_Main_REVISED_WITH_TRANSLATION_FINAL.m) | `Fig5B_Q1_mainphenotype_plot.m` |
| [`Fig5C_Q2_REVISED_3D_MECHANISM_FINAL_v2.m`](code/CTMC/Fig5C_Q2_REVISED_3D_MECHANISM_FINAL_v2.m) | `Fig5C_Q2_REVISED_3D_MECHANISM_FINAL_v2.m` |
| [`SI_M1_MesoMacro_Interface_Robustness_FINAL.m`](code/CTMC/SI_M1_MesoMacro_Interface_Robustness_FINAL.m) | `SI_M1_MesoMacro_Interface_Robustness_FINAL.m` |
| [`SI_M2_Macro_Parameter_Benchmarks_FINAL.m`](code/CTMC/SI_M2_Macro_Parameter_Benchmarks_FINAL.m) | `SI_M2_Macro_Parameter_Benchmarks_FINAL.m` |
| [`SI_M3_Macro_SSA_Numerical_Validation_FINAL.m`](code/CTMC/SI_M3_Macro_SSA_Numerical_Validation_FINAL.m) | `SI_M3_Macro_SSA_Numerical_Validation_FINAL.m` |
| [`SI_M5_Q1_StatisticalValidation_FINAL.m`](code/CTMC/SI_M5_Q1_StatisticalValidation_FINAL.m) | `SI_M5_Q1_StatisticalValidation_FINAL.m` |
| [`SI_Q2_1_SpatialArrangement_Robustness_FINAL.m`](code/CTMC/SI_Q2_1_SpatialArrangement_Robustness_FINAL.m) | `SI_Q2_1_SpatialArrangement_Robustness_FINAL.m` |
| [`SI_Q2_2_NeighbourhoodAccess_Mechanism_FINAL.m`](code/CTMC/SI_Q2_2_NeighbourhoodAccess_Mechanism_FINAL.m) | `SI_Q2_2_NeighbourhoodAccess_Mechanism_FINAL.m` |
| [`SI_Q2_3_ConditionGenerality_FINAL.m`](code/CTMC/SI_Q2_3_ConditionGenerality_FINAL.m) | `SI_Q2_3_ConditionGenerality_FINAL.m` |
| [`Fig6_ExternalConsistency_FSP1_LIGHT_FINAL_TOPJOURNAL.m`](code/External-Validation/Fig6_ExternalConsistency_FSP1_LIGHT_FINAL_TOPJOURNAL.m) | `validation/finalfinal/validation.m` |
| [`Meso_to_Macro_Risk_Interface_FINAL_v2.m`](code/CTMC/Meso_to_Macro_Risk_Interface_FINAL_v2.m) | `macro/v2/Meso_to_Macro_Risk_Interface_FINAL_v2.m` |
| [`Macro_Q2_SpatialArrangement_FINAL.m`](code/CTMC/Macro_Q2_SpatialArrangement_FINAL.m) | `macro/v2/Macro_Q2_SpatialArrangement_FINAL.m` |
| [`Macro_Q2_NeighbourhoodAccess_Mechanism_FINAL.m`](code/CTMC/Macro_Q2_NeighbourhoodAccess_Mechanism_FINAL.m) | `macro/v2/Macro_Q2_NeighbourhoodAccess_Mechanism_FINAL.m` |
| [`Macro_Q2_TargetedConditionMap_FINAL.m`](code/CTMC/Macro_Q2_TargetedConditionMap_FINAL.m) | `macro/v2/Macro_Q2_TargetedConditionMap_FINAL.m` |
| [`macro_sim_core_hazard_SSA_FINAL.m`](code/CTMC/macro_sim_core_hazard_SSA_FINAL.m) | `macro/v2/macro_sim_core_hazard_SSA_FINAL.m` |

The selected v2 meso–macro interface is the author's orientation-compatible implementation. It handles the condition-by-repair-time layout saved by the finalized cell-state analysis.

## Edits

All released MATLAB files use UTF-8 and LF line endings. Function-file names match the functions they expose. Script names provide descriptive model and analysis roles.

- `regimemap.m` — Use the operating-system temporary directory for default MATLAB worker storage.
- `MLWTWS_SIfinal.m` — Resolve the fate-route MAT file relative to this model directory.
- `cellstate.m` — Resolve the fate-route MAT file relative to this model directory.
- `cellstate_SI.m` — Read the LAYOUTFIX output actually written by the selected cell-state driver.
- `fateroute.m` — Add peer TD-CME and SD-RDME directories to the existing input search.
- `jointmechanism_final.m` — Add the final fate-route output directory to the existing input search.
- `jointmechanism_finalSI.m` — Add the final fate-route output directory to the existing input search.
- `spatialcoarseSI.m` — Add the final fate-route and peer SD-RDME directories to the existing input search.
- `temporalmemorySI.m` — Add the final fate-route output directory to the existing input search.
- `SI_M1_MesoMacro_Interface_Robustness_FINAL.m` — Read the two inputs directly from the peer LPO-FPT output directories.
- `Meso_to_Macro_Risk_Interface_FINAL_v2.m` — Read the two inputs directly from the peer LPO-FPT output directories.
- `validation.m` — Add the peer CTMC Q1 output and original publisher workbook names and sheet names.

Every released file was compared with its decoded source and the explicit path/input changes above. Parameters, reaction propensities, random seeds, numerical solvers, thresholds, bootstrap calculations and statistical decision rules retain their source implementations.

## Verification

MATLAB R2022a with Statistics and Machine Learning Toolbox and Parallel Computing Toolbox was used for verification.

- Code Analyzer checked all 41 MATLAB files for syntax errors.
- The TD-CME and SD-RDME baseline entry points completed using their existing smoke settings in a disposable working copy.
- The external consistency analysis completed using the original publisher workbooks and archived Q1 output. Its three CSV outputs matched the archived study outputs byte for byte.
- All 320 cells read from the study-curated external workbooks matched the corresponding publisher workbook cells.
- The orientation-compatible meso–macro interface completed its raw-FPT reconstruction, reference-state and time-scale checks. The 48-state interface CSV matched the archived production CSV byte for byte.
- The threshold/weight robustness, interface robustness, Q1 figure, Q2 mechanism figure and three Q2 SI plotting entry points completed using archived study intermediates in a disposable working copy.
- The LPO fate-route driver completed with its original smoke settings, loading archived final microscopic inputs from the peer directories and using its existing serial fallback after MATLAB worker startup failed.
- The CTMC core completed a small scalar-hazard simulation with exact seeded repeatability. Site-map and inherited-cell-mark hazard modes also completed with non-negative populations and consistent total counts.

The verification used reduced simulations and archived production intermediates. Full production ensembles and every sensitivity scan have additional runtime requirements. Some original graphics labels produced TeX interpreter warnings under R2022a while exporting their output files.
