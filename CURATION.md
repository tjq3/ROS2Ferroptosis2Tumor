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

| Released file | Source | Comparison before the naming update |
| --- | --- | --- |
| [`TDCME_IronRecyclingAblation.m`](code/TD-CME/TDCME_IronRecyclingAblation.m) | `TDCME_fixedC2_ablation_FINAL.m` | Byte-identical |
| [`TDCME_ROSBurdenAndFPT.m`](code/TD-CME/TDCME_ROSBurdenAndFPT.m) | `TDCME_matchtohighorderandSI2.m` | Byte-identical |
| [`TDCME_RegimeMap.m`](code/TD-CME/TDCME_RegimeMap.m) | `regimemap.m` | Path, workbook or worker-storage edits |
| [`TDCME_CoreRobustness.m`](code/TD-CME/TDCME_CoreRobustness.m) | `robustness.m` | Byte-identical |
| [`TDCME_TimeStepConvergence.m`](code/TD-CME/TDCME_TimeStepConvergence.m) | `robustness3.m` | Byte-identical |
| [`SDRDME_RedoxBufferScan.m`](code/SD-RDME/SDRDME_RedoxBufferScan.m) | `SDRDME_Fe3Fe2redox.m` | Byte-identical |
| [`SDRDME_RedoxBufferFigure.m`](code/SD-RDME/SDRDME_RedoxBufferFigure.m) | `SDRDME_Fe3Fe2redoxhighlight.m` | Byte-identical |
| [`SDRDME_Counterfactuals.m`](code/SD-RDME/SDRDME_Counterfactuals.m) | `SDRDME_UHBFandSI.m` | Byte-identical |
| [`SDRDME_LocalBurdenAndFPT.m`](code/SD-RDME/SDRDME_LocalBurdenAndFPT.m) | `SDRDME_burdenandFPT.m` | Byte-identical |
| [`SDRDME_HeldOutValidation.m`](code/SD-RDME/SDRDME_HeldOutValidation.m) | `SDRDME_c1.m` | Byte-identical |
| [`SDRDME_DescriptorScreening.m`](code/SD-RDME/SDRDME_DescriptorScreening.m) | `SDRDME_c1BDXHBandSI.m` | Byte-identical |
| [`SDRDME_GeometryRobustness.m`](code/SD-RDME/SDRDME_GeometryRobustness.m) | `SDRDME_geometry.m` | Byte-identical |
| [`SDRDME_CoreRobustness.m`](code/SD-RDME/SDRDME_CoreRobustness.m) | `SDRDME_robustness.m` | Byte-identical |
| [`LPO_ThresholdAndWeightRobustness.m`](code/LPO-FPT/LPO_ThresholdAndWeightRobustness.m) | `MLWTWS_SIfinal.m` | Path, workbook or worker-storage edits |
| [`LPO_RepairTimeAndFateRoutes.m`](code/LPO-FPT/LPO_RepairTimeAndFateRoutes.m) | `cellstate.m` | Path, workbook or worker-storage edits |
| [`LPO_FateRoutes.m`](code/LPO-FPT/LPO_FateRoutes.m) | `fateroute.m` | Path, workbook or worker-storage edits |
| [`LPO_JointMechanism.m`](code/LPO-FPT/LPO_JointMechanism.m) | `jointmechanism_final.m` | Path, workbook or worker-storage edits |
| [`LPO_JointMechanismRobustness.m`](code/LPO-FPT/LPO_JointMechanismRobustness.m) | `jointmechanism_finalSI.m` | Path, workbook or worker-storage edits |
| [`LPO_SpatialCoarseGraining.m`](code/LPO-FPT/LPO_SpatialCoarseGraining.m) | `spatialcoarseSI.m` | Path, workbook or worker-storage edits |
| [`LPO_TemporalMemoryAndHorizon.m`](code/LPO-FPT/LPO_TemporalMemoryAndHorizon.m) | `temporalmemorySI.m` | Path, workbook or worker-storage edits |
| [`CTMC_CommitmentClearanceScan.m`](code/CTMC/CTMC_CommitmentClearanceScan.m) | `Fig5B_Q1_mainphenotype_calu.m` | Byte-identical |
| [`CTMC_CommitmentClearanceFigure.m`](code/CTMC/CTMC_CommitmentClearanceFigure.m) | `Fig5B_Q1_mainphenotype_plot.m` | Byte-identical |
| [`CTMC_SpatialArrangementFigure.m`](code/CTMC/CTMC_SpatialArrangementFigure.m) | `Fig5C_Q2_REVISED_3D_MECHANISM_FINAL_v2.m` | Byte-identical |
| [`CTMC_RiskInterfaceRobustness.m`](code/CTMC/CTMC_RiskInterfaceRobustness.m) | `SI_M1_MesoMacro_Interface_Robustness_FINAL.m` | Path, workbook or worker-storage edits |
| [`CTMC_ParameterBenchmarks.m`](code/CTMC/CTMC_ParameterBenchmarks.m) | `SI_M2_Macro_Parameter_Benchmarks_FINAL.m` | Byte-identical |
| [`CTMC_NumericalValidation.m`](code/CTMC/CTMC_NumericalValidation.m) | `SI_M3_Macro_SSA_Numerical_Validation_FINAL.m` | Byte-identical |
| [`CTMC_CommitmentClearanceValidation.m`](code/CTMC/CTMC_CommitmentClearanceValidation.m) | `SI_M5_Q1_StatisticalValidation_FINAL.m` | Byte-identical |
| [`CTMC_SpatialArrangementRobustness.m`](code/CTMC/CTMC_SpatialArrangementRobustness.m) | `SI_Q2_1_SpatialArrangement_Robustness_FINAL.m` | Byte-identical |
| [`CTMC_NeighbourhoodAccessFigure.m`](code/CTMC/CTMC_NeighbourhoodAccessFigure.m) | `SI_Q2_2_NeighbourhoodAccess_Mechanism_FINAL.m` | Byte-identical |
| [`CTMC_SpatialArrangementGenerality.m`](code/CTMC/CTMC_SpatialArrangementGenerality.m) | `SI_Q2_3_ConditionGenerality_FINAL.m` | Byte-identical |
| [`ExternalValidation_FSP1Consistency.m`](code/External-Validation/ExternalValidation_FSP1Consistency.m) | `validation/finalfinal/validation.m` | Path, workbook or worker-storage edits |
| [`CTMC_MesoscopicRiskInterface.m`](code/CTMC/CTMC_MesoscopicRiskInterface.m) | `macro/v2/Meso_to_Macro_Risk_Interface_FINAL_v2.m` | Path, workbook or worker-storage edits |
| [`CTMC_SpatialArrangementSimulation.m`](code/CTMC/CTMC_SpatialArrangementSimulation.m) | `macro/v2/Macro_Q2_SpatialArrangement_FINAL.m` | Byte-identical |
| [`CTMC_NeighbourhoodAccessMechanism.m`](code/CTMC/CTMC_NeighbourhoodAccessMechanism.m) | `macro/v2/Macro_Q2_NeighbourhoodAccess_Mechanism_FINAL.m` | Byte-identical |
| [`CTMC_SpatialArrangementConditionScan.m`](code/CTMC/CTMC_SpatialArrangementConditionScan.m) | `macro/v2/Macro_Q2_TargetedConditionMap_FINAL.m` | Byte-identical |
| [`CTMC_EventDrivenSimulation.m`](code/CTMC/CTMC_EventDrivenSimulation.m) | `macro/v2/macro_sim_core_hazard_SSA_FINAL.m` | Byte-identical |

The selected v2 meso–macro interface is the author's orientation-compatible implementation. It handles the condition-by-repair-time layout saved by the finalized cell-state analysis.

## Edits

Public MATLAB filenames identify the model and analysis. Function declarations, calls, help comments, run instructions, and documentation use the corresponding entry-point names. The original MAT filenames and directory names allow the same archived results to be loaded.

- `TDCME_RegimeMap.m` — Use the operating-system temporary directory for default MATLAB worker storage.
- `LPO_ThresholdAndWeightRobustness.m` — Resolve the fate-route MAT file relative to this model directory.
- `LPO_RepairTimeAndFateRoutes.m` — Resolve the fate-route MAT file relative to this model directory.
- `LPO_FateRoutes.m` — Add peer TD-CME and SD-RDME directories to the existing input search.
- `LPO_JointMechanism.m` — Add the final fate-route output directory to the existing input search.
- `LPO_JointMechanismRobustness.m` — Add the final fate-route output directory to the existing input search.
- `LPO_SpatialCoarseGraining.m` — Add the final fate-route and peer SD-RDME directories to the existing input search.
- `LPO_TemporalMemoryAndHorizon.m` — Add the final fate-route output directory to the existing input search.
- `CTMC_RiskInterfaceRobustness.m` — Read the two inputs directly from the peer LPO-FPT output directories.
- `CTMC_MesoscopicRiskInterface.m` — Read the two inputs directly from the peer LPO-FPT output directories.
- `ExternalValidation_FSP1Consistency.m` — Add the peer CTMC Q1 output and original publisher workbook names and sheet names.

Every released file was compared with its decoded source and the explicit path/input changes above. Parameters, reaction propensities, random seeds, numerical solvers, thresholds, bootstrap calculations and statistical decision rules retain their source implementations.

## Direct comparison with the supplied sources

Before the naming update, 25 of the 36 MATLAB files matched their source files byte for byte. The other 11 differed in the path, workbook-selection, and worker-storage changes documented above. The naming update changes script filenames and references, including MATLAB function declarations, calls, and file checks.

The four frozen MAT files were compared as binary files against their originals in `macro/v2/`. All four matched byte for byte.

| Frozen input | SHA-256 |
| --- | --- |
| `Macro_Hazard_PreQ1_DT_Refinement_FINAL.mat` | `d63360742ca3adcba27b1e89a1172149e50efcd8fb434f9633b16c20f2a1ffcf` |
| `Macro_SSA_LOCKED.mat` | `ce766b0dbaa4bde73bef1347a3f33f3b47ff2e730030dd0209da56f58982a19e` |
| `Macro_mu_LOCKED.mat` | `bf8c2c8d1582d981cbd93cb07c7e16ed6594c8db43dc5c9acdda6b27f1acff2f` |
| `Macro_r_LOCKED.mat` | `ea55353b61cbb09f41abba9a04125d61598044c85011a94e17eab4dd1cd4462c` |

The selected interface's source is `macro/v2/Meso_to_Macro_Risk_Interface_FINAL_v2.m`. Relative to the author's v1 interface, the author's v2 normalizes 4-by-12 or 12-by-4 state arrays, accepts both Fig4B structure layouts, and checks the sizes of first-passage samples before comparing them. Its empirical RMST and RMST-equivalent hazard formulas are the same as v1. The release adds the peer-directory paths shown above to v2.

Input resolution is part of the edit scope. The CTMC interface and its window-robustness script first check the peer LPO-FPT outputs. The external check first checks the peer CTMC Q1 output. Other added locations act as fallbacks after the original current-directory and script-directory checks. The external check retains the curated workbook names as its first choice and adds the publisher workbook names as fallbacks.

The eight excluded MATLAB entry points remain available in the author's original archives. Their standalone outputs are outside this release. The archived figure files, exploratory versions, intermediate results, checkpoints and experimental workbooks are also outside the release, apart from the four explicitly listed frozen inputs. The external-data README provides publisher download links and exact worksheets and ranges.

## Verification

MATLAB R2022a with Statistics and Machine Learning Toolbox and Parallel Computing Toolbox was used for verification.

- Code Analyzer checked the final 36 MATLAB files for syntax errors.
- The TD-CME and SD-RDME baseline entry points completed using their existing smoke settings in a disposable working copy.
- The external consistency analysis completed using the original publisher workbooks and archived Q1 output. Its three CSV outputs matched the archived study outputs byte for byte.
- All 320 cells read from the study-curated external workbooks matched the corresponding publisher workbook cells.
- The orientation-compatible meso–macro interface completed its raw-FPT reconstruction, reference-state and time-scale checks. The 48-state interface CSV matched the archived production CSV byte for byte.
- The threshold/weight robustness, interface robustness, Q1 figure, Q2 mechanism figure and three Q2 SI plotting entry points completed using archived study intermediates in a disposable working copy.
- The LPO fate-route driver completed with its original smoke settings, loading archived final microscopic inputs from the peer directories and using its existing serial fallback after MATLAB worker startup failed.
- The CTMC core completed a small scalar-hazard simulation with exact seeded repeatability. Site-map and inherited-cell-mark hazard modes also completed with non-negative populations and consistent total counts.

The verification used reduced simulations and archived production intermediates. Full production ensembles and every sensitivity scan have additional runtime requirements. Some original graphics labels produced TeX interpreter warnings under R2022a while exporting their output files.
