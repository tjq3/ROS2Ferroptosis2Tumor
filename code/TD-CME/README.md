# Cytosolic TD-CME

Run `TDCME_Fig2ABC` in this directory with `run_level = 'final'` to compare matched CIR and Gamma-OU environments, low-order ROS statistics, and full-state ROS tail and first-passage risk. The script writes `TDCME_Fig2ABC_results_final.mat` here for the mesoscopic workflow.

| Entry point | Analysis |
| --- | --- |
| `TDCME_Fig2ABC.m` | Core cytosolic model and baseline comparison |
| `TDCME_RegimeMap.m` | Timescale–dispersion sensitivity of ROS first-passage risk |
| `TDCME_CoreRobustness.m` | Threshold quantile, observation horizon, time step and seed robustness |
| `TDCME_TimeStepConvergence.m` | Fixed-threshold time-step convergence |
| `TDCME_fixedC2_ablation_FINAL.m` | Constant versus stochastic Fe³⁺ recycling coefficient |

Each analysis can be run independently. Each writes its own figures, tables or MAT files according to its export settings.
