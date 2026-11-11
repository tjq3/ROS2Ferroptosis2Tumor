%% TDCME_IronRecyclingAblation
% =========================================================================
% Targeted fixed-c2 ablation for the cytosolic strict TD-CME.
%
% Scientific question
% -------------------
% Does the main CIR-versus-GammaOU path-risk separation persist when the
% Fe3+ -> Fe2+ recycling coefficient is no longer stochastic, but is fixed
% at its stationary mean mu_c2?
%
% Two model arms are compared:
%
%   FULL-c2      : c1(t) stochastic and c2(t) stochastic
%                  (identical to the finalized strict TD-CME kernel)
%
%   FIXED-c2     : c1(t) stochastic and c2(t) == mu_c2
%
% Within each arm, CIR and GammaOU are matched in the stationary Gamma
% marginal and exponential correlation time of c1(t); they differ only in
% higher-order temporal path architecture. The biochemical reaction network,
% iron pool, observation horizon, time step, integrated-hazard reaction
% update, and c1 parameters are unchanged from the finalized TD-CME analysis.
%
% Threshold strategy
% ------------------
% Two threshold schemes are reported:
%
%   1) Arm-specific pooled Q0.90 threshold (PRIMARY ABLATION COMPARISON)
%      Each arm receives its own independently calibrated pooled Q0.90
%      endpoint threshold. This compares CIR-vs-GammaOU separation at the
%      same nominal operating point within each model arm.
%
%   2) Full-model absolute threshold (SECONDARY FIXED-ESTIMAND AUDIT)
%      The absolute threshold calibrated in FULL-c2 is also applied unchanged
%      to FIXED-c2. This checks whether the conclusion depends on threshold
%      re-calibration after the ablation.
%
% Primary effect signs
% --------------------
%   DeltaF_60    = P(T_GammaOU <= 60 min) - P(T_CIR <= 60 min)
%   DeltaRMST_60 = RMST_CIR(60 min) - RMST_GammaOU(60 min)
%
% Positive values indicate earlier/higher GammaOU ROS-burden crossing.
%
% Important interpretation
% ------------------------
% These are ROS-burden first-passage times, NOT ferroptosis FPTs.
%
% Numerical kernel
% ----------------
%   - N = 1200 total iron molecules
%   - dt = 2 s, Tobs = 3600 s
%   - matched stationary-binomial Fe2+/Fe3+ initialization
%   - matched stationary-Gamma c1 initialization
%   - FULL-c2 additionally uses matched stationary-Gamma c2 initialization
%   - CIR: exact finite-step noncentral-chi-square transition; trapezoidal
%          approximation of the within-step c1/c2 integral
%   - GammaOU: exact exponential decay plus compound-Poisson jumps; exact
%              within-step integrated rate conditional on jump times/sizes
%   - FIXED-c2: exact within-step integral I2 = mu_c2 * dt
%   - chemical reactions: integrated hazards -> binomial finite-pool updates
%
% Outputs
% -------
%   TDCME_fixedC2_ablation_FINAL_figure.{pdf,png,tif,fig}
%   TDCME_fixedC2_ablation_summary.csv
%   TDCME_fixedC2_ablation_thresholds.csv
%   TDCME_fixedC2_ablation_effectComparison.csv
%   TDCME_fixedC2_ablation_results.mat
%   TDCME_fixedC2_ablation_run.log
%
% Recommended workflow
% --------------------
%   run_level = 'smoke' -> 'pilot' -> 'final'
%
% Required toolbox
% ----------------
%   Statistics and Machine Learning Toolbox
% =========================================================================

clear; clc; close all;

%% =========================================================================
% 0. Run configuration
% =========================================================================
run_level = 'final';      % 'smoke' | 'pilot' | 'final'

cfg = struct();
cfg.export_figures = true;
cfg.export_tables  = true;
cfg.save_results   = true;
cfg.export_dpi     = 600;

switch lower(run_level)
    case 'smoke'
        cfg.n_cal  = 1200;
        cfg.n_eval = 2500;
        cfg.n_boot = 200;
    case 'pilot'
        cfg.n_cal  = 5000;
        cfg.n_eval = 10000;
        cfg.n_boot = 500;
    case 'final'
        % Matched to the scale used in the finalized strict TD-CME analyses.
        cfg.n_cal  = 12000;
        cfg.n_eval = 30000;
        cfg.n_boot = 1000;
    otherwise
        error('run_level must be smoke, pilot, or final.');
end

cfg.rec_dt = 10;          % s; output recording only, not model dt
cfg.q_Mb   = 0.90;
cfg.smooth_fano = 5;      % 5 x 10 s = 50 s display smoothing

% Fixed seeds. FULL and FIXED arms deliberately reuse the same calibration
% and evaluation seeds so that their x1(0) and c1(0) draws are aligned.
% The resulting arm-to-arm comparison is still treated descriptively rather
% than as a paired-random-number estimator.
seed_cal  = 20261021;
seed_eval = 20261022;
seed_boot_full  = 20261023;
seed_boot_fixed = 20261024;

script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir), script_dir = pwd; end
out_dir = fullfile(script_dir, ['results_TDCME_fixedC2_' lower(run_level)]);
if ~exist(out_dir,'dir'), mkdir(out_dir); end

out_stem = fullfile(out_dir, 'TDCME_fixedC2_ablation_FINAL_figure');
out_summary = fullfile(out_dir, 'TDCME_fixedC2_ablation_summary.csv');
out_thresholds = fullfile(out_dir, 'TDCME_fixedC2_ablation_thresholds.csv');
out_compare = fullfile(out_dir, 'TDCME_fixedC2_ablation_effectComparison.csv');
out_mat = fullfile(out_dir, 'TDCME_fixedC2_ablation_results.mat');
log_file = fullfile(out_dir, 'TDCME_fixedC2_ablation_run.log');
if exist(log_file,'file'), delete(log_file); end
diary(log_file);
cleanup_diary = onCleanup(@() diary('off')); %#ok<NASGU>

%% =========================================================================
% 1. Baseline parameters: identical to finalized strict TD-CME
% =========================================================================
p = struct();
p.N = 1200;

% Oxidation / Fenton channel c1
p.kappa1 = 3.0e-4;        % s^-1
p.shape1 = 0.20;
p.mu_c1  = 1.5e-5;       % s^-1
p.theta1 = p.mu_c1 / p.shape1;
p.a1     = p.shape1 * p.kappa1;  % GammaOU jump-event intensity

% Reduction / recycling channel c2
p.kappa2 = 4.0e-4;        % s^-1
p.shape2 = 3.50;
p.mu_c2  = 1.7e-5;       % s^-1
p.theta2 = p.mu_c2 / p.shape2;
p.a2     = p.shape2 * p.kappa2;

% Mean-rate iron partition used in the finalized TD-CME initialization.
p.lambda  = p.mu_c1 + p.mu_c2;
p.p_x1_ss = p.mu_c2 / p.lambda;
p.x1_ss   = p.N * p.p_x1_ss;
p.x2_ss   = p.N - p.x1_ss;

% Low-order reference quantities.
p.var_c1 = p.shape1 * p.theta1^2;
p.var_c2 = p.shape2 * p.theta2^2;
p.v_eff  = p.mu_c1 * p.x1_ss;
p.F_int  = 1 - 2 * p.mu_c1 * p.mu_c2 / p.lambda^2;

dh_dc1 = p.N * p.mu_c2^2 / p.lambda^2;
dh_dc2 = p.N * p.mu_c1^2 / p.lambda^2;
p.A1 = dh_dc1^2 * p.var_c1;
p.A2 = dh_dc2^2 * p.var_c2;

%% =========================================================================
% 2. Simulation settings
% =========================================================================
sim = struct();
sim.dt      = 2.0;        % s
sim.T_total = 3600;       % s = 60 min

validate_configuration(sim,cfg);

fprintf('============================================================\n');
fprintf('TD-CME fixed-c2 ablation | %s\n', upper(run_level));
fprintf('============================================================\n');
fprintf('N               = %d\n',p.N);
fprintf('c1 mean/shape   = %.6e / %.3f\n',p.mu_c1,p.shape1);
fprintf('c1 kappa        = %.6e s^-1\n',p.kappa1);
fprintf('c2 mean/shape   = %.6e / %.3f\n',p.mu_c2,p.shape2);
fprintf('c2 kappa        = %.6e s^-1\n',p.kappa2);
fprintf('CV1^2 / CV2^2  = %.4f / %.4f\n',1/p.shape1,1/p.shape2);
fprintf('dt / Tobs       = %.1f s / %.1f min\n',sim.dt,sim.T_total/60);
fprintf('n_cal/n_eval    = %d / %d\n',cfg.n_cal,cfg.n_eval);
fprintf('n_boot          = %d\n',cfg.n_boot);
fprintf('============================================================\n\n');

%% =========================================================================
% 3. Independent threshold calibration
% =========================================================================
fprintf('Calibration: FULL-c2 arm...\n');
rng(seed_cal,'twister');
CAL_full = simulate_calibration_arm(p,sim,cfg.n_cal,'stochastic');

fprintf('Calibration: FIXED-c2 arm...\n');
rng(seed_cal,'twister');
CAL_fixed = simulate_calibration_arm(p,sim,cfg.n_cal,'fixed');

M_full = max(1, round(quantile([CAL_full.CIR.Y60; CAL_full.GOU.Y60], cfg.q_Mb)));
M_fixed = max(1, round(quantile([CAL_fixed.CIR.Y60; CAL_fixed.GOU.Y60], cfg.q_Mb)));
M_common = M_full;

fprintf('\nThreshold calibration:\n');
fprintf('  FULL-c2  own pooled Q%.2f threshold = %d\n',cfg.q_Mb,M_full);
fprintf('  FIXED-c2 own pooled Q%.2f threshold = %d\n',cfg.q_Mb,M_fixed);
fprintf('  Common absolute audit threshold     = %d (FULL-c2 threshold)\n\n',M_common);

thresholdTable = table( ...
    ["FULL-c2";"FIXED-c2"], ...
    [M_full;M_fixed], ...
    repmat(M_common,2,1), ...
    [mean(CAL_full.CIR.Y60);mean(CAL_fixed.CIR.Y60)], ...
    [mean(CAL_full.GOU.Y60);mean(CAL_fixed.GOU.Y60)], ...
    'VariableNames',{'ModelArm','OwnQ90Threshold','CommonFullThreshold', ...
    'CalibrationMeanY_CIR','CalibrationMeanY_GammaOU'});

%% =========================================================================
% 4. Independent evaluation runs
% =========================================================================
% Threshold vector columns:
%   col 1 = arm-specific threshold
%   col 2 = common FULL-c2 absolute threshold

fprintf('Evaluation: FULL-c2 arm...\n');
rng(seed_eval,'twister');
R_full = simulate_evaluation_arm( ...
    p,sim,cfg,cfg.n_eval,'stochastic',[M_full M_common]);

fprintf('Evaluation: FIXED-c2 arm...\n');
rng(seed_eval,'twister');
R_fixed = simulate_evaluation_arm( ...
    p,sim,cfg,cfg.n_eval,'fixed',[M_fixed M_common]);

%% =========================================================================
% 5. Summaries at arm-specific and common thresholds
% =========================================================================
S_full_own   = summarize_arm(R_full,1,sim.T_total,M_full);
S_full_common = summarize_arm(R_full,2,sim.T_total,M_common);
S_fixed_own  = summarize_arm(R_fixed,1,sim.T_total,M_fixed);
S_fixed_common = summarize_arm(R_fixed,2,sim.T_total,M_common);

% FULL own and FULL common are identical by construction, but both labels are
% retained so the exported table is structurally explicit.

fprintf('Bootstrap: FULL-c2...\n');
rng(seed_boot_full,'twister');
B_full_own = bootstrap_effects(R_full,1,sim.T_total,cfg.n_boot);
rng(seed_boot_full+100,'twister');
B_full_common = bootstrap_effects(R_full,2,sim.T_total,cfg.n_boot);

fprintf('Bootstrap: FIXED-c2...\n');
rng(seed_boot_fixed,'twister');
B_fixed_own = bootstrap_effects(R_fixed,1,sim.T_total,cfg.n_boot);
rng(seed_boot_fixed+100,'twister');
B_fixed_common = bootstrap_effects(R_fixed,2,sim.T_total,cfg.n_boot);

CI_full_own = bootstrap_ci(B_full_own);
CI_full_common = bootstrap_ci(B_full_common);
CI_fixed_own = bootstrap_ci(B_fixed_own);
CI_fixed_common = bootstrap_ci(B_fixed_common);

%% =========================================================================
% 6. Low-order diagnostics
% =========================================================================
t = R_full.t(:);
mean_reference = p.v_eff .* t;

% Full stochastic-c2 LNA includes c1 and c2 environmental terms.
fano_LNA_full = finite_time_fano_general( ...
    t,p.v_eff,p.F_int, ...
    [p.A1 p.A2],[p.kappa1 p.kappa2]);

% Fixed-c2 LNA keeps only the c1 environmental term.
fano_LNA_fixed = finite_time_fano_general( ...
    t,p.v_eff,p.F_int,p.A1,p.kappa1);

fano_full_c_s  = movmean(R_full.CIR.fano,cfg.smooth_fano,'omitnan');
fano_full_g_s  = movmean(R_full.GOU.fano,cfg.smooth_fano,'omitnan');
fano_fixed_c_s = movmean(R_fixed.CIR.fano,cfg.smooth_fano,'omitnan');
fano_fixed_g_s = movmean(R_fixed.GOU.fano,cfg.smooth_fano,'omitnan');

% Endpoint low-order comparison.
low_full = low_order_summary(R_full);
low_fixed = low_order_summary(R_fixed);

%% =========================================================================
% 7. Effect-retention comparison: FIXED-c2 relative to FULL-c2
% =========================================================================
% Independent-bootstrap approximation to the difference between model-arm
% effect estimates. The two arms are not claimed to be paired estimators.
DeltaDeltaF_own_boot = B_fixed_own.dF - B_full_own.dF;
DeltaDeltaR_own_boot = B_fixed_own.dRMST_min - B_full_own.dRMST_min;

DeltaDeltaF_common_boot = B_fixed_common.dF - B_full_common.dF;
DeltaDeltaR_common_boot = B_fixed_common.dRMST_min - B_full_common.dRMST_min;

compare = struct();
compare.own.DeltaDeltaF = S_fixed_own.DeltaF - S_full_own.DeltaF;
compare.own.DeltaDeltaF_CI = prctile(DeltaDeltaF_own_boot,[2.5 97.5]);
compare.own.DeltaDeltaRMST_min = ...
    S_fixed_own.DeltaRMST_min - S_full_own.DeltaRMST_min;
compare.own.DeltaDeltaRMST_CI = prctile(DeltaDeltaR_own_boot,[2.5 97.5]);

compare.common.DeltaDeltaF = S_fixed_common.DeltaF - S_full_common.DeltaF;
compare.common.DeltaDeltaF_CI = prctile(DeltaDeltaF_common_boot,[2.5 97.5]);
compare.common.DeltaDeltaRMST_min = ...
    S_fixed_common.DeltaRMST_min - S_full_common.DeltaRMST_min;
compare.common.DeltaDeltaRMST_CI = prctile(DeltaDeltaR_common_boot,[2.5 97.5]);

compare.own.RetentionF = safe_ratio(S_fixed_own.DeltaF,S_full_own.DeltaF);
compare.own.RetentionRMST = safe_ratio( ...
    S_fixed_own.DeltaRMST_min,S_full_own.DeltaRMST_min);
compare.common.RetentionF = safe_ratio( ...
    S_fixed_common.DeltaF,S_full_common.DeltaF);
compare.common.RetentionRMST = safe_ratio( ...
    S_fixed_common.DeltaRMST_min,S_full_common.DeltaRMST_min);

fprintf('\n==================== PRIMARY RESULTS ====================\n');
print_summary('FULL-c2 | own Q90',S_full_own,CI_full_own);
print_summary('FIXED-c2 | own Q90',S_fixed_own,CI_fixed_own);
fprintf('Effect retention (FIXED/FULL, own thresholds):\n');
fprintf('  DeltaF retention       = %.4f\n',compare.own.RetentionF);
fprintf('  DeltaRMST retention    = %.4f\n',compare.own.RetentionRMST);
fprintf('  FIXED-FULL DeltaF      = %.5f, 95%% CI [%.5f, %.5f]\n', ...
    compare.own.DeltaDeltaF,compare.own.DeltaDeltaF_CI(1),compare.own.DeltaDeltaF_CI(2));
fprintf('  FIXED-FULL DeltaRMST   = %.4f min, 95%% CI [%.4f, %.4f]\n', ...
    compare.own.DeltaDeltaRMST_min, ...
    compare.own.DeltaDeltaRMST_CI(1),compare.own.DeltaDeltaRMST_CI(2));

fprintf('\nFixed absolute-threshold audit (M=%d):\n',M_common);
print_summary('FULL-c2 | common M',S_full_common,CI_full_common);
print_summary('FIXED-c2 | common M',S_fixed_common,CI_fixed_common);
fprintf('=========================================================\n\n');

%% =========================================================================
% 8. Export summary tables
% =========================================================================
summaryTable = build_summary_table( ...
    S_full_own,S_full_common,S_fixed_own,S_fixed_common, ...
    CI_full_own,CI_full_common,CI_fixed_own,CI_fixed_common, ...
    low_full,low_fixed,M_full,M_fixed,M_common);

comparisonTable = table( ...
    ["Arm-specific Q90";"Full absolute M"], ...
    [compare.own.RetentionF;compare.common.RetentionF], ...
    [compare.own.RetentionRMST;compare.common.RetentionRMST], ...
    [compare.own.DeltaDeltaF;compare.common.DeltaDeltaF], ...
    [compare.own.DeltaDeltaF_CI(1);compare.common.DeltaDeltaF_CI(1)], ...
    [compare.own.DeltaDeltaF_CI(2);compare.common.DeltaDeltaF_CI(2)], ...
    [compare.own.DeltaDeltaRMST_min;compare.common.DeltaDeltaRMST_min], ...
    [compare.own.DeltaDeltaRMST_CI(1);compare.common.DeltaDeltaRMST_CI(1)], ...
    [compare.own.DeltaDeltaRMST_CI(2);compare.common.DeltaDeltaRMST_CI(2)], ...
    'VariableNames',{'ThresholdScheme','Retention_DeltaF', ...
    'Retention_DeltaRMST','FixedMinusFull_DeltaF', ...
    'FixedMinusFull_DeltaF_CI_L','FixedMinusFull_DeltaF_CI_U', ...
    'FixedMinusFull_DeltaRMST_min', ...
    'FixedMinusFull_DeltaRMST_CI_L','FixedMinusFull_DeltaRMST_CI_U'});

if cfg.export_tables
    writetable(summaryTable,out_summary);
    writetable(thresholdTable,out_thresholds);
    writetable(comparisonTable,out_compare);
end

%% =========================================================================
% 9. Publication-oriented SI figure
% =========================================================================
% Panel logic:
%   A  Mean cumulative ROS Y(t)
%   B  Fano factor
%   C  FPT CDF at arm-specific Q0.90 thresholds
%   D  Endpoint ROS-load survival at 60 min
%   E  DeltaF_60 with 95% bootstrap CI
%   F  DeltaRMST_60 with 95% bootstrap CI

colCIR = [0.12 0.34 0.62];
colGOU = [0.74 0.18 0.12];
colFull = [0.20 0.21 0.23];
colFixed = [0.48 0.24 0.62];
colTheory = [0.48 0.48 0.48];

fig = figure('Color','w','Units','centimeters','Position',[2 2 18.3 10.8]);
tl = tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');

% A. Mean cumulative ROS
ax = nexttile(tl,1); hold(ax,'on');
plot(ax,t/60,R_full.CIR.meanY,'-','Color',colCIR,'LineWidth',1.15);
plot(ax,t/60,R_full.GOU.meanY,'-','Color',colGOU,'LineWidth',1.15);
plot(ax,t/60,R_fixed.CIR.meanY,'--','Color',colCIR,'LineWidth',1.15);
plot(ax,t/60,R_fixed.GOU.meanY,'--','Color',colGOU,'LineWidth',1.15);
plot(ax,t/60,mean_reference,':','Color',colTheory,'LineWidth',1.0);
xlabel(ax,'Time (min)'); ylabel(ax,'Mean cumulative ROS events');
title(ax,'Mean burden');
legend(ax,{'CIR, full c_2','Gamma-OU, full c_2', ...
    'CIR, fixed c_2','Gamma-OU, fixed c_2','Mean-rate reference'}, ...
    'Location','northwest','Box','off','FontSize',5.5);
format_axes(ax); panel_label(ax,'A');

% B. Fano factor
ax = nexttile(tl,2); hold(ax,'on');
plot(ax,t/60,fano_full_c_s,'-','Color',colCIR,'LineWidth',1.10);
plot(ax,t/60,fano_full_g_s,'-','Color',colGOU,'LineWidth',1.10);
plot(ax,t/60,fano_fixed_c_s,'--','Color',colCIR,'LineWidth',1.10);
plot(ax,t/60,fano_fixed_g_s,'--','Color',colGOU,'LineWidth',1.10);
plot(ax,t/60,fano_LNA_full,':','Color',colFull,'LineWidth',1.0);
plot(ax,t/60,fano_LNA_fixed,':','Color',colFixed,'LineWidth',1.0);
xlabel(ax,'Time (min)'); ylabel(ax,'Fano factor');
title(ax,'Low-order dispersion');
format_axes(ax); panel_label(ax,'B');

% C. FPT CDF at arm-specific thresholds
ax = nexttile(tl,3); hold(ax,'on');
tgrid = linspace(0,sim.T_total,241)';
plot(ax,tgrid/60,empirical_fpt_cdf(R_full.CIR.FPT(:,1),tgrid), ...
    '-','Color',colCIR,'LineWidth',1.15);
plot(ax,tgrid/60,empirical_fpt_cdf(R_full.GOU.FPT(:,1),tgrid), ...
    '-','Color',colGOU,'LineWidth',1.15);
plot(ax,tgrid/60,empirical_fpt_cdf(R_fixed.CIR.FPT(:,1),tgrid), ...
    '--','Color',colCIR,'LineWidth',1.15);
plot(ax,tgrid/60,empirical_fpt_cdf(R_fixed.GOU.FPT(:,1),tgrid), ...
    '--','Color',colGOU,'LineWidth',1.15);
xlabel(ax,'Time (min)'); ylabel(ax,'P(T_{M_b} \leq t)');
title(ax,sprintf('ROS-burden FPT (M_{full}=%d, M_{fixed}=%d)',M_full,M_fixed));
format_axes(ax); panel_label(ax,'C');

% D. Endpoint ROS-load survival at 60 min
ax = nexttile(tl,4); hold(ax,'on');
Yall = [R_full.CIR.Y60;R_full.GOU.Y60;R_fixed.CIR.Y60;R_fixed.GOU.Y60];
mmax = max([M_full M_fixed ceil(quantile(Yall,0.997))]);
mgrid = 0:mmax;
plot(ax,mgrid,tail_from_samples(R_full.CIR.Y60,mgrid),'-','Color',colCIR,'LineWidth',1.10);
plot(ax,mgrid,tail_from_samples(R_full.GOU.Y60,mgrid),'-','Color',colGOU,'LineWidth',1.10);
plot(ax,mgrid,tail_from_samples(R_fixed.CIR.Y60,mgrid),'--','Color',colCIR,'LineWidth',1.10);
plot(ax,mgrid,tail_from_samples(R_fixed.GOU.Y60,mgrid),'--','Color',colGOU,'LineWidth',1.10);
xline(ax,M_common,':','Color',colTheory,'LineWidth',0.9);
xlabel(ax,'Cumulative ROS events at 60 min'); ylabel(ax,'P(Y_{60} \geq m)');
title(ax,'Endpoint burden tail');
format_axes(ax); panel_label(ax,'D');

% E. DeltaF_60
ax = nexttile(tl,5); hold(ax,'on');
xpos = [1 2];
yOwn = [S_full_own.DeltaF S_fixed_own.DeltaF];
ciOwn = [CI_full_own.dF;CI_fixed_own.dF];
yCom = [S_full_common.DeltaF S_fixed_common.DeltaF];
ciCom = [CI_full_common.dF;CI_fixed_common.dF];
errorbar(ax,xpos-0.08,yOwn,yOwn-ciOwn(:,1)',ciOwn(:,2)'-yOwn, ...
    'o','Color',colFull,'MarkerFaceColor','w','MarkerSize',4.2,'LineWidth',1.0,'CapSize',4);
errorbar(ax,xpos+0.08,yCom,yCom-ciCom(:,1)',ciCom(:,2)'-yCom, ...
    's','Color',colFixed,'MarkerFaceColor','w','MarkerSize',4.0,'LineWidth',1.0,'CapSize',4);
yline(ax,0,':','Color',[0.65 0.65 0.65],'LineWidth',0.8);
set(ax,'XTick',xpos,'XTickLabel',{'Full c_2','Fixed c_2'});
ylabel(ax,'\DeltaF_{60}');
title(ax,'Crossing-probability separation');
legend(ax,{'Arm-specific Q_{0.90}','Common full-model M_b'}, ...
    'Location','best','Box','off','FontSize',5.5);
format_axes(ax); panel_label(ax,'E');

% F. DeltaRMST_60
ax = nexttile(tl,6); hold(ax,'on');
yOwn = [S_full_own.DeltaRMST_min S_fixed_own.DeltaRMST_min];
ciOwn = [CI_full_own.dRMST;CI_fixed_own.dRMST];
yCom = [S_full_common.DeltaRMST_min S_fixed_common.DeltaRMST_min];
ciCom = [CI_full_common.dRMST;CI_fixed_common.dRMST];
errorbar(ax,xpos-0.08,yOwn,yOwn-ciOwn(:,1)',ciOwn(:,2)'-yOwn, ...
    'o','Color',colFull,'MarkerFaceColor','w','MarkerSize',4.2,'LineWidth',1.0,'CapSize',4);
errorbar(ax,xpos+0.08,yCom,yCom-ciCom(:,1)',ciCom(:,2)'-yCom, ...
    's','Color',colFixed,'MarkerFaceColor','w','MarkerSize',4.0,'LineWidth',1.0,'CapSize',4);
yline(ax,0,':','Color',[0.65 0.65 0.65],'LineWidth',0.8);
set(ax,'XTick',xpos,'XTickLabel',{'Full c_2','Fixed c_2'});
ylabel(ax,'\DeltaRMST_{60} (min)');
title(ax,'First-passage-time separation');
format_axes(ax); panel_label(ax,'F');

if cfg.export_figures
    export_figure_all(fig,out_stem,cfg.export_dpi);
end

%% =========================================================================
% 10. Save full results
% =========================================================================
metadata = struct();
metadata.run_level = run_level;
metadata.seed_cal = seed_cal;
metadata.seed_eval = seed_eval;
metadata.seed_boot_full = seed_boot_full;
metadata.seed_boot_fixed = seed_boot_fixed;
metadata.description = [ ...
    'Targeted TD-CME ablation: stochastic c1 with either stochastic c2 ', ...
    'or c2 fixed at mu_c2.'];
metadata.primary_threshold_scheme = 'arm-specific pooled Q0.90';
metadata.secondary_threshold_scheme = 'FULL-c2 absolute threshold';

if cfg.save_results
    save(out_mat,'metadata','p','sim','cfg', ...
        'CAL_full','CAL_fixed','R_full','R_fixed', ...
        'M_full','M_fixed','M_common', ...
        'S_full_own','S_full_common','S_fixed_own','S_fixed_common', ...
        'B_full_own','B_full_common','B_fixed_own','B_fixed_common', ...
        'CI_full_own','CI_full_common','CI_fixed_own','CI_fixed_common', ...
        'low_full','low_fixed','compare', ...
        'summaryTable','thresholdTable','comparisonTable', ...
        'mean_reference','fano_LNA_full','fano_LNA_fixed','-v7.3');
end

fprintf('Completed. Outputs are in:\n  %s\n',out_dir);
fprintf('  Figure: %s.*\n',out_stem);
fprintf('  Summary: %s\n',out_summary);
fprintf('  MAT: %s\n',out_mat);

%% =========================================================================
% Local functions
% =========================================================================

function validate_configuration(sim,cfg)
    assert(sim.dt>0 && sim.T_total>0,'dt and T_total must be positive.');
    assert(abs(sim.T_total/sim.dt-round(sim.T_total/sim.dt))<1e-10, ...
        'T_total must be divisible by dt.');
    assert(abs(cfg.rec_dt/sim.dt-round(cfg.rec_dt/sim.dt))<1e-10, ...
        'rec_dt must be divisible by dt.');
    assert(cfg.q_Mb>0 && cfg.q_Mb<1,'q_Mb must lie in (0,1).');
end

function OUT = simulate_calibration_arm(p,sim,n,c2_mode)
    n_steps = round(sim.T_total/sim.dt);

    % Same initialization convention as the finalized strict TD-CME.
    x1_0 = binornd(p.N,p.p_x1_ss,n,1);
    c1_0 = gamrnd(p.shape1,p.theta1,n,1);

    if strcmpi(c2_mode,'stochastic')
        c2_0 = gamrnd(p.shape2,p.theta2,n,1);
    elseif strcmpi(c2_mode,'fixed')
        c2_0 = p.mu_c2*ones(n,1);
    else
        error('Unknown c2_mode: %s',c2_mode);
    end

    % CIR label.
    x1_c = x1_0;
    x2_c = p.N-x1_c;
    y_c  = zeros(n,1);
    c1_c = c1_0;
    c2_c = c2_0;

    % GammaOU label.
    x1_g = x1_0;
    x2_g = p.N-x1_g;
    y_g  = zeros(n,1);
    c1_g = c1_0;
    c2_g = c2_0;

    for it=1:n_steps
        % ------------------------------- CIR ------------------------------
        c1_old = c1_c;
        c1_c = cir_exact_transition( ...
            c1_c,p.kappa1,p.shape1,p.theta1,sim.dt);
        I1_c = 0.5*(c1_old+c1_c)*sim.dt;

        if strcmpi(c2_mode,'stochastic')
            c2_old = c2_c;
            c2_c = cir_exact_transition( ...
                c2_c,p.kappa2,p.shape2,p.theta2,sim.dt);
            I2_c = 0.5*(c2_old+c2_c)*sim.dt;
        else
            I2_c = p.mu_c2*sim.dt*ones(n,1);
        end

        pr1_c = min(max(1-exp(-I1_c),0),1);
        pr2_c = min(max(1-exp(-I2_c),0),1);
        r1_c = binornd(x1_c,pr1_c);
        r2_c = binornd(x2_c,pr2_c);
        x1_c = x1_c-r1_c+r2_c;
        x2_c = p.N-x1_c;
        y_c = y_c+r1_c;

        % ----------------------------- GammaOU ---------------------------
        [c1_g,I1_g] = gamma_ou_step_with_integral( ...
            c1_g,p.a1,p.kappa1,p.theta1,sim.dt);

        if strcmpi(c2_mode,'stochastic')
            [c2_g,I2_g] = gamma_ou_step_with_integral( ...
                c2_g,p.a2,p.kappa2,p.theta2,sim.dt);
        else
            I2_g = p.mu_c2*sim.dt*ones(n,1);
        end

        pr1_g = min(max(1-exp(-I1_g),0),1);
        pr2_g = min(max(1-exp(-I2_g),0),1);
        r1_g = binornd(x1_g,pr1_g);
        r2_g = binornd(x2_g,pr2_g);
        x1_g = x1_g-r1_g+r2_g;
        x2_g = p.N-x1_g;
        y_g = y_g+r1_g;
    end

    OUT.CIR.Y60 = y_c;
    OUT.GOU.Y60 = y_g;
end

function OUT = simulate_evaluation_arm(p,sim,cfg,n,c2_mode,thresholds)
    n_steps = round(sim.T_total/sim.dt);
    rec_stride = round(cfg.rec_dt/sim.dt);
    rec_idx = 0:rec_stride:n_steps;
    t_rec = rec_idx(:)*sim.dt;
    n_thr = numel(thresholds);

    x1_0 = binornd(p.N,p.p_x1_ss,n,1);
    c1_0 = gamrnd(p.shape1,p.theta1,n,1);

    if strcmpi(c2_mode,'stochastic')
        c2_0 = gamrnd(p.shape2,p.theta2,n,1);
    elseif strcmpi(c2_mode,'fixed')
        c2_0 = p.mu_c2*ones(n,1);
    else
        error('Unknown c2_mode: %s',c2_mode);
    end

    x1_c = x1_0; x2_c = p.N-x1_c; y_c = zeros(n,1);
    c1_c = c1_0; c2_c = c2_0;

    x1_g = x1_0; x2_g = p.N-x1_g; y_g = zeros(n,1);
    c1_g = c1_0; c2_g = c2_0;

    FPT_c = nan(n,n_thr);
    FPT_g = nan(n,n_thr);

    mean_c = zeros(numel(rec_idx),1);
    mean_g = zeros(numel(rec_idx),1);
    var_c  = zeros(numel(rec_idx),1);
    var_g  = zeros(numel(rec_idx),1);
    ir=2;

    for it=1:n_steps
        t_now = it*sim.dt;

        % ------------------------------- CIR ------------------------------
        c1_old = c1_c;
        c1_c = cir_exact_transition( ...
            c1_c,p.kappa1,p.shape1,p.theta1,sim.dt);
        I1_c = 0.5*(c1_old+c1_c)*sim.dt;

        if strcmpi(c2_mode,'stochastic')
            c2_old = c2_c;
            c2_c = cir_exact_transition( ...
                c2_c,p.kappa2,p.shape2,p.theta2,sim.dt);
            I2_c = 0.5*(c2_old+c2_c)*sim.dt;
        else
            I2_c = p.mu_c2*sim.dt*ones(n,1);
        end

        pr1_c = min(max(1-exp(-I1_c),0),1);
        pr2_c = min(max(1-exp(-I2_c),0),1);
        r1_c = binornd(x1_c,pr1_c);
        r2_c = binornd(x2_c,pr2_c);
        x1_c = x1_c-r1_c+r2_c;
        x2_c = p.N-x1_c;
        y_c = y_c+r1_c;

        for j=1:n_thr
            hit = isnan(FPT_c(:,j)) & y_c>=thresholds(j);
            FPT_c(hit,j)=t_now;
        end

        % ----------------------------- GammaOU ---------------------------
        [c1_g,I1_g] = gamma_ou_step_with_integral( ...
            c1_g,p.a1,p.kappa1,p.theta1,sim.dt);

        if strcmpi(c2_mode,'stochastic')
            [c2_g,I2_g] = gamma_ou_step_with_integral( ...
                c2_g,p.a2,p.kappa2,p.theta2,sim.dt);
        else
            I2_g = p.mu_c2*sim.dt*ones(n,1);
        end

        pr1_g = min(max(1-exp(-I1_g),0),1);
        pr2_g = min(max(1-exp(-I2_g),0),1);
        r1_g = binornd(x1_g,pr1_g);
        r2_g = binornd(x2_g,pr2_g);
        x1_g = x1_g-r1_g+r2_g;
        x2_g = p.N-x1_g;
        y_g = y_g+r1_g;

        for j=1:n_thr
            hit = isnan(FPT_g(:,j)) & y_g>=thresholds(j);
            FPT_g(hit,j)=t_now;
        end

        if ir<=numel(rec_idx) && it==rec_idx(ir)
            mean_c(ir)=mean(y_c);
            mean_g(ir)=mean(y_g);
            var_c(ir)=var(y_c,0);
            var_g(ir)=var(y_g,0);
            ir=ir+1;
        end

        if mod(it,round(600/sim.dt))==0
            fprintf('  %s: t = %.1f min\n',upper(c2_mode),t_now/60);
        end
    end

    fano_c = var_c./max(mean_c,eps);
    fano_g = var_g./max(mean_g,eps);
    fano_c(1)=nan; fano_g(1)=nan;

    OUT.mode = c2_mode;
    OUT.thresholds = thresholds(:)';
    OUT.t = t_rec;
    OUT.CIR.meanY = mean_c;
    OUT.CIR.varY = var_c;
    OUT.CIR.fano = fano_c;
    OUT.CIR.Y60 = y_c;
    OUT.CIR.FPT = FPT_c;
    OUT.GOU.meanY = mean_g;
    OUT.GOU.varY = var_g;
    OUT.GOU.fano = fano_g;
    OUT.GOU.Y60 = y_g;
    OUT.GOU.FPT = FPT_g;
end

function S = summarize_arm(R,thr_col,tau,threshold)
    fpt_c = R.CIR.FPT(:,thr_col);
    fpt_g = R.GOU.FPT(:,thr_col);

    event_c = ~isnan(fpt_c) & fpt_c<=tau;
    event_g = ~isnan(fpt_g) & fpt_g<=tau;

    S.threshold = threshold;
    S.meanY_CIR = mean(R.CIR.Y60);
    S.meanY_GOU = mean(R.GOU.Y60);
    S.fanoY_CIR = var(R.CIR.Y60,0)/max(mean(R.CIR.Y60),eps);
    S.fanoY_GOU = var(R.GOU.Y60,0)/max(mean(R.GOU.Y60),eps);
    S.P60_CIR = mean(event_c);
    S.P60_GOU = mean(event_g);
    S.DeltaF = S.P60_GOU-S.P60_CIR;
    S.RMST_CIR_min = rmst_from_fpt(fpt_c,tau)/60;
    S.RMST_GOU_min = rmst_from_fpt(fpt_g,tau)/60;
    S.DeltaRMST_min = S.RMST_CIR_min-S.RMST_GOU_min;
end

function B = bootstrap_effects(R,thr_col,tau,n_boot)
    fpt_c = R.CIR.FPT(:,thr_col);
    fpt_g = R.GOU.FPT(:,thr_col);
    n_c = numel(fpt_c);
    n_g = numel(fpt_g);

    event_c = ~isnan(fpt_c) & fpt_c<=tau;
    event_g = ~isnan(fpt_g) & fpt_g<=tau;

    tc = fpt_c; tg = fpt_g;
    tc(isnan(tc))=tau; tg(isnan(tg))=tau;
    tc=min(tc,tau); tg=min(tg,tau);

    B.dF = zeros(n_boot,1);
    B.dRMST_min = zeros(n_boot,1);

    for b=1:n_boot
        ic = randi(n_c,n_c,1);
        ig = randi(n_g,n_g,1);
        B.dF(b) = mean(event_g(ig))-mean(event_c(ic));
        B.dRMST_min(b) = (mean(tc(ic))-mean(tg(ig)))/60;
    end
end

function CI = bootstrap_ci(B)
    CI.dF = prctile(B.dF,[2.5 97.5]);
    CI.dRMST = prctile(B.dRMST_min,[2.5 97.5]);
end

function L = low_order_summary(R)
    L.meanY_CIR = mean(R.CIR.Y60);
    L.meanY_GOU = mean(R.GOU.Y60);
    L.relMeanDiff = abs(L.meanY_GOU-L.meanY_CIR) / ...
        max(0.5*(L.meanY_GOU+L.meanY_CIR),eps);

    L.fano_CIR = var(R.CIR.Y60,0)/max(L.meanY_CIR,eps);
    L.fano_GOU = var(R.GOU.Y60,0)/max(L.meanY_GOU,eps);
    L.relFanoDiff = abs(L.fano_GOU-L.fano_CIR) / ...
        max(0.5*(L.fano_GOU+L.fano_CIR),eps);
end

function T = build_summary_table( ...
    Sfo,Sfc,Sxo,Sxc,CIfo,CIfc,CIxo,CIxc, ...
    lowFull,lowFixed,Mfull,Mfixed,Mcommon)

    ModelArm = ["FULL-c2";"FULL-c2";"FIXED-c2";"FIXED-c2"];
    ThresholdScheme = ["Arm-specific Q0.90";"Full absolute M"; ...
        "Arm-specific Q0.90";"Full absolute M"];
    Threshold = [Mfull;Mcommon;Mfixed;Mcommon];
    S = {Sfo;Sfc;Sxo;Sxc};
    C = {CIfo;CIfc;CIxo;CIxc};

    n=4;
    MeanY_CIR=zeros(n,1); MeanY_GammaOU=zeros(n,1);
    Fano_CIR=zeros(n,1); Fano_GammaOU=zeros(n,1);
    P60_CIR=zeros(n,1); P60_GammaOU=zeros(n,1); DeltaF=zeros(n,1);
    DeltaF_CI_L=zeros(n,1); DeltaF_CI_U=zeros(n,1);
    RMST_CIR_min=zeros(n,1); RMST_GammaOU_min=zeros(n,1);
    DeltaRMST_min=zeros(n,1); DeltaRMST_CI_L=zeros(n,1); DeltaRMST_CI_U=zeros(n,1);

    for i=1:n
        s=S{i}; c=C{i};
        MeanY_CIR(i)=s.meanY_CIR; MeanY_GammaOU(i)=s.meanY_GOU;
        Fano_CIR(i)=s.fanoY_CIR; Fano_GammaOU(i)=s.fanoY_GOU;
        P60_CIR(i)=s.P60_CIR; P60_GammaOU(i)=s.P60_GOU;
        DeltaF(i)=s.DeltaF;
        DeltaF_CI_L(i)=c.dF(1); DeltaF_CI_U(i)=c.dF(2);
        RMST_CIR_min(i)=s.RMST_CIR_min; RMST_GammaOU_min(i)=s.RMST_GOU_min;
        DeltaRMST_min(i)=s.DeltaRMST_min;
        DeltaRMST_CI_L(i)=c.dRMST(1); DeltaRMST_CI_U(i)=c.dRMST(2);
    end

    RelMeanDiff = [lowFull.relMeanDiff;lowFull.relMeanDiff; ...
        lowFixed.relMeanDiff;lowFixed.relMeanDiff];
    RelFanoDiff = [lowFull.relFanoDiff;lowFull.relFanoDiff; ...
        lowFixed.relFanoDiff;lowFixed.relFanoDiff];

    T = table(ModelArm,ThresholdScheme,Threshold, ...
        MeanY_CIR,MeanY_GammaOU,RelMeanDiff,Fano_CIR,Fano_GammaOU,RelFanoDiff, ...
        P60_CIR,P60_GammaOU,DeltaF,DeltaF_CI_L,DeltaF_CI_U, ...
        RMST_CIR_min,RMST_GammaOU_min,DeltaRMST_min, ...
        DeltaRMST_CI_L,DeltaRMST_CI_U);
end

function F = finite_time_fano_general(t,V,F_int,A,kappa)
    t=t(:);
    A=A(:)'; kappa=kappa(:)';
    F=nan(size(t));
    idx=t>0;
    tt=t(idx);

    var_ext=zeros(size(tt));
    for j=1:numel(A)
        G = tt./kappa(j) - (1-exp(-kappa(j).*tt))./(kappa(j)^2);
        var_ext = var_ext + 2*A(j).*G;
    end
    F(idx)=F_int + var_ext./(V.*tt);
    F(~idx)=0;
end

function F = empirical_fpt_cdf(fpt,tgrid)
    F=arrayfun(@(t) mean(~isnan(fpt) & fpt<=t),tgrid);
    F=F(:);
end

function tail = tail_from_samples(y,mgrid)
    tail=arrayfun(@(m) mean(y>=m),mgrid);
end

function value = rmst_from_fpt(fpt,tau)
    z=fpt;
    z(isnan(z))=tau;
    z=min(z,tau);
    value=mean(z);
end

function c_next = cir_exact_transition(c,kappa,shape,theta,dt)
    if dt==0
        c_next=c;
        return;
    end

    e=exp(-kappa*dt);
    one_minus_e=1-e;
    scale=0.5*theta*one_minus_e;
    df=2*shape;
    nc=2*e.*max(c,0)./(theta*one_minus_e);

    if exist('ncx2rnd','file')==2
        c_next=scale.*ncx2rnd(df,nc);
    else
        K=poissrnd(nc./2);
        c_next=scale.*gamrnd(df/2+K,2);
    end
    c_next=max(c_next,0);
end

function [c_new,I] = gamma_ou_step_with_integral(c_old,a,kappa,theta,dt)
    if dt==0
        c_new=c_old;
        I=zeros(size(c_old));
        return;
    end

    n=numel(c_old);
    decay=exp(-kappa*dt);
    c_new=c_old.*decay;
    I=c_old.*(1-decay)./kappa;

    Nj=poissrnd(a*dt,n,1);
    idx=find(Nj>0);

    for ii=1:numel(idx)
        j=idx(ii);
        m=Nj(j);
        U=dt.*rand(m,1);
        Y=gamrnd(1,theta,m,1);

        c_new(j)=c_new(j)+sum(Y.*exp(-kappa.*(dt-U)));
        I(j)=I(j)+sum(Y.*(1-exp(-kappa.*(dt-U)))./kappa);
    end

    c_new=max(c_new,0);
    I=max(I,0);
end

function r = safe_ratio(num,den)
    if abs(den)<1e-12
        r=nan;
    else
        r=num/den;
    end
end

function print_summary(label,S,CI)
    fprintf('%s\n',label);
    fprintf('  threshold       = %d\n',S.threshold);
    fprintf('  mean Y60 C/G    = %.4f / %.4f\n',S.meanY_CIR,S.meanY_GOU);
    fprintf('  Fano Y60 C/G    = %.4f / %.4f\n',S.fanoY_CIR,S.fanoY_GOU);
    fprintf('  Pcross60 C/G    = %.4f / %.4f\n',S.P60_CIR,S.P60_GOU);
    fprintf('  DeltaF60        = %.5f, 95%% CI [%.5f, %.5f]\n', ...
        S.DeltaF,CI.dF(1),CI.dF(2));
    fprintf('  RMST60 C/G      = %.4f / %.4f min\n',S.RMST_CIR_min,S.RMST_GOU_min);
    fprintf('  DeltaRMST60     = %.4f min, 95%% CI [%.4f, %.4f]\n', ...
        S.DeltaRMST_min,CI.dRMST(1),CI.dRMST(2));
end

function panel_label(ax,label)
    text(ax,0.01,0.99,label,'Units','normalized','FontWeight','bold', ...
        'FontSize',10,'VerticalAlignment','top','HorizontalAlignment','left');
end

function format_axes(ax)
    set(ax,'Box','off','TickDir','out','LineWidth',0.75, ...
        'FontName','Arial','FontSize',7);
end

function export_figure_all(fig,stem,dpi)
    if exist('exportgraphics','file')==2
        exportgraphics(fig,[stem '.pdf'],'ContentType','vector');
        exportgraphics(fig,[stem '.png'],'Resolution',dpi);
        exportgraphics(fig,[stem '.tif'],'Resolution',dpi);
    else
        print(fig,[stem '.pdf'],'-dpdf','-bestfit');
        print(fig,[stem '.png'],'-dpng',sprintf('-r%d',dpi));
        print(fig,[stem '.tif'],'-dtiff',sprintf('-r%d',dpi));
    end
    savefig(fig,[stem '.fig']);
end
