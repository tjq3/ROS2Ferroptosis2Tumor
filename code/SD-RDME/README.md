# Mitochondrial SD-RDME

Run `SDRDME_LocalBurdenAndFPT` in this directory with `RUN_MODE = 'final'`. The script runs the next-subvolume stochastic simulation and writes `Fig3BC_predictorIndependent_data.mat` here. Descriptor screening and held-out validation read that file. The mesoscopic model also uses it for the spatial input.

| Entry point | Analysis |
| --- | --- |
| `SDRDME_LocalBurdenAndFPT.m` | Global burden matching, local ROS burden and first passage |
| `SDRDME_DescriptorScreening.m` | Candidate-field screening, quantile and scale sensitivity |
| `SDRDME_HeldOutValidation.m` | Independent trajectories for held-out reaction-field validation |
| `SDRDME_Counterfactuals.m` | Uniform, reaction-only, diffusion-only and full heterogeneity |
| `SDRDME_RedoxBufferScan.m` | Recycling-capacity and diffusion-heterogeneity scan |
| `SDRDME_RedoxBufferFigure.m` | Mechanism figure from the redox scan CSVs under `condition_redox_buffer_results/final/` |
| `SDRDME_CoreRobustness.m` | Four robustness analyses of the spatial model |
| `SDRDME_GeometryRobustness.m` | Sensitivity to mitochondrial outline and crista geometry |

The counterfactual, scan, geometry and robustness scripts are standalone.
