% =========================================================================
% TDCME_Fig2ABC_FullHorizonPanelC_FINAL.m
%
% Purpose
% -------
% Generate the first three main-text panels for the cytosolic strict TD-CME:
%
%   Panel A. Matched baseline
%       Matched stationary Gamma marginals and matched exponential memory
%       for CIR and GammaOU rate environments, while their sample paths
%       retain continuous-diffusive versus jump-intermittent architectures.
%
%   Panel B. Low-order approximate consistency
%       Full-state TD-CME mean cumulative ROS and Fano factor remain closely
%       aligned between CIR and GammaOU and are compared with the shared
%       first-moment closure / finite-time LNA reference.
%
%   Panel C. High-order risk separation over the complete 0-60 min horizon
%       The panel displays the cumulative ROS-load tail Y(60 min) and the
%       ROS-burden FPT CDF over the same 0-60 min observation interval.
%       A common independently calibrated threshold M_b is used for both.
%
% The same script also creates a companion SI figure:
%
%   SI-A. Shared low-order closure/LNA prediction
%   SI-B. Full-state CIR-GammaOU low-order deviations in 0-60 min
%   SI-C. Long-time Fano behavior relative to the LNA reference
%   SI-D. Shared LNA Fano decomposition
%
% Key consistency with previous strict TD-CME code
% ------------------------------------------------
%   N = 1200
%   c1: kappa1=3.0e-4 s^-1, shape1=0.20, mu_c1=1.5e-5 s^-1
%   c2: kappa2=4.0e-4 s^-1, shape2=3.50, mu_c2=1.7e-5 s^-1
%   dt=2 s, T=3600 s
%   matched stationary-Gamma rate initialization
%   matched stationary-binomial Fe2+/Fe3+ initialization
%   CIR exact transition
%   GammaOU exact exponential decay + jump-time integral
%   binomial reaction update from integrated hazards
%
% Important interpretation
% ------------------------
% Panel C reports cumulative ROS burden and ROS-burden first-passage time
% over 0-60 min, NOT ferroptosis FPT. Ferroptosis FPT is defined later at
% the mesoscopic cell-fate layer.
%
% Required toolbox
% ----------------
% Statistics and Machine Learning Toolbox:
%   binornd, gamrnd, poissrnd, ncx2rnd (fallback provided if unavailable)
%
% Recommended workflow
% --------------------
% Set run_level = 'smoke', then 'pilot', then 'final'.
% =========================================================================

clear; clc; close all;

%% =========================================================================
% 0. Run configuration
% =========================================================================
run_level = 'final';       % 'smoke' | 'pilot' | 'final'

make_main_figure = true;
make_SI_figure   = true;
export_figures   = true;
export_tables    = true;
save_results_mat = true;

switch lower(run_level)
    case 'smoke'
        cfg.n_cal   = 3000;
        cfg.n_eval  = 5000;
        cfg.n_env   = 5000;
        cfg.n_long  = 1000;
        cfg.n_boot  = 200;
    case 'pilot'
        cfg.n_cal   = 8000;
        cfg.n_eval  = 12000;
        cfg.n_env   = 10000;
        cfg.n_long  = 3000;
        cfg.n_boot  = 500;
    case 'final'
        cfg.n_cal   = 12000;
        cfg.n_eval  = 30000;
        cfg.n_env   = 20000;
        cfg.n_long  = 8000;
        cfg.n_boot  = 1000;
    otherwise
        error('run_level must be ''smoke'', ''pilot'', or ''final''.');
end

cfg.n_show_paths = 16;
cfg.rec_dt_main  = 10;      % output recording only; model dt remains 2 s
cfg.rec_dt_long  = 60;
cfg.T_long       = 28800;   % 480 min, SI diagnostic
cfg.smooth_main  = 5;       % 5 x 10 s = 50 s smoothing for displayed Fano
cfg.smooth_long  = 5;       % 5 x 60 s = 5 min smoothing for displayed Fano

seed_env  = 20260801;
seed_cal  = 20260802;
seed_eval = 20260803;
seed_long = 20260804;
seed_boot = 20260805;

%% =========================================================================
% 1. Parameters: identical to the previous strict TD-CME main setting
% =========================================================================
p = struct();

p.N = 1200;

% Oxidation / Fenton channel c1
p.kappa1 = 3.0e-4;
p.shape1 = 0.20;
p.mu_c1  = 1.5e-5;
p.theta1 = p.mu_c1 / p.shape1;
p.a1     = p.shape1 * p.kappa1;

% Reduction / recycling channel c2
p.kappa2 = 4.0e-4;
p.shape2 = 3.50;
p.mu_c2  = 1.7e-5;
p.theta2 = p.mu_c2 / p.shape2;
p.a2     = p.shape2 * p.kappa2;

% Mean iron partition
p.lambda  = p.mu_c1 + p.mu_c2;
p.p_x1_ss = p.mu_c2 / p.lambda;
p.x1_ss   = p.N * p.p_x1_ss;
p.x2_ss   = p.N - p.x1_ss;

% Shared low-order closure / LNA quantities
p.var_c1 = p.shape1 * p.theta1^2;
p.var_c2 = p.shape2 * p.theta2^2;

p.v_eff = p.mu_c1 * p.x1_ss;  % N*mu1*mu2/(mu1+mu2)

p.F_int = 1 - 2 * p.mu_c1 * p.mu_c2 / p.lambda^2;

dh_dc1 = p.N * p.mu_c2^2 / p.lambda^2;
dh_dc2 = p.N * p.mu_c1^2 / p.lambda^2;

p.A1 = dh_dc1^2 * p.var_c1;
p.A2 = dh_dc2^2 * p.var_c2;

p.F_ext1 = 2 * p.A1 / (p.v_eff * p.kappa1);
p.F_ext2 = 2 * p.A2 / (p.v_eff * p.kappa2);
p.F_ext  = p.F_ext1 + p.F_ext2;
p.F_LNA_long = p.F_int + p.F_ext;

%% =========================================================================
% 2. Simulation and threshold settings
% =========================================================================
sim = struct();

sim.dt      = 2.0;       % s
sim.T_total = 3600;      % s = 60 min

% One common full-horizon burden threshold, independently calibrated from
% pooled Y(60 min), is used for both the endpoint tail and 0-60 min FPT.
sim.q_Mb = 0.90;

validate_time_grid(sim, cfg);

fprintf('============================================================\n');
fprintf('TD-CME Fig.2A-C: matched baseline -> low-order -> high-order\n');
fprintf('============================================================\n');
fprintf('Run level: %s\n', run_level);
fprintf('N = %d\n', p.N);
fprintf('c1: kappa=%.6e, shape=%.3f, mean=%.6e, theta=%.6e\n', ...
    p.kappa1, p.shape1, p.mu_c1, p.theta1);
fprintf('c2: kappa=%.6e, shape=%.3f, mean=%.6e, theta=%.6e\n', ...
    p.kappa2, p.shape2, p.mu_c2, p.theta2);
fprintf('mean-field Fe2+ = %.3f, Fe3+ = %.3f\n', p.x1_ss, p.x2_ss);
fprintf('dt=%.1f s, full observation horizon=0-%.1f min\n', ...
    sim.dt, sim.T_total/60);
fprintf('n_cal=%d, n_eval=%d, n_env=%d, n_long=%d, n_boot=%d\n', ...
    cfg.n_cal, cfg.n_eval, cfg.n_env, cfg.n_long, cfg.n_boot);
fprintf('============================================================\n\n');

%% =========================================================================
% 3. Panel A input audit: matched marginals and memory, distinct paths
% =========================================================================
rng(seed_env, 'twister');
fprintf('Environmental matching audit starts...\n');
ENV = simulate_environment_audit(p, sim, cfg);
fprintf('Environmental matching audit finished.\n\n');

%% =========================================================================
% 4. Independent threshold calibration
% =========================================================================
rng(seed_cal, 'twister');
fprintf('Independent calibration run starts...\n');
CAL = simulate_calibration(p, sim, cfg.n_cal);
fprintf('Independent calibration run finished.\n\n');

Mb = max(1, round(quantile([CAL.CIR.Y60; CAL.GOU.Y60], sim.q_Mb)));

fprintf('Independent pooled full-horizon threshold:\n');
fprintf('  M_b = round(Q%.2f pooled Y(60 min)) = %d\n\n', sim.q_Mb, Mb);

%% =========================================================================
% 5. Full-state evaluation run
% =========================================================================
rng(seed_eval, 'twister');
fprintf('Full-state evaluation run starts...\n');
R = simulate_full_evaluation(p, sim, cfg, cfg.n_eval, Mb);
fprintf('Full-state evaluation run finished.\n\n');

%% =========================================================================
% 6. Shared low-order analytic predictions
% =========================================================================
t_main = R.t(:);

mean_LNA = p.v_eff .* t_main;
fano_LNA = finite_time_fano( ...
    t_main, p.v_eff, p.F_int, p.A1, p.kappa1, p.A2, p.kappa2);

fano_c_s = movmean(R.CIR.fano, cfg.smooth_main, 'omitnan');
fano_g_s = movmean(R.GOU.fano, cfg.smooth_main, 'omitnan');

den_mean = 0.5 .* (R.CIR.meanY + R.GOU.meanY);
den_fano = 0.5 .* (fano_c_s + fano_g_s);

rel_mean_CG = abs(R.GOU.meanY - R.CIR.meanY) ./ max(den_mean, eps);
rel_fano_CG = abs(fano_g_s - fano_c_s) ./ max(den_fano, eps);

rel_mean_CG(~isfinite(rel_mean_CG)) = nan;
rel_fano_CG(~isfinite(rel_fano_CG)) = nan;

valid_low = t_main >= 120 & den_mean > 1 & ...
    isfinite(rel_mean_CG) & isfinite(rel_fano_CG);

max_rel_mean_CG = max(rel_mean_CG(valid_low), [], 'omitnan');
max_rel_fano_CG = max(rel_fano_CG(valid_low), [], 'omitnan');

% Relative deviations of each full-state model from the shared closure/LNA.
rel_mean_c_theory = abs(R.CIR.meanY - mean_LNA) ./ max(mean_LNA, eps);
rel_mean_g_theory = abs(R.GOU.meanY - mean_LNA) ./ max(mean_LNA, eps);
rel_fano_c_theory = abs(fano_c_s - fano_LNA) ./ max(fano_LNA, eps);
rel_fano_g_theory = abs(fano_g_s - fano_LNA) ./ max(fano_LNA, eps);

rel_mean_c_theory(~isfinite(rel_mean_c_theory)) = nan;
rel_mean_g_theory(~isfinite(rel_mean_g_theory)) = nan;
rel_fano_c_theory(~isfinite(rel_fano_c_theory)) = nan;
rel_fano_g_theory(~isfinite(rel_fano_g_theory)) = nan;

%% =========================================================================
% 7. Panel C high-order risk separation over the full 0-60 min horizon
% =========================================================================
Y60_c = R.CIR.Y60;
Y60_g = R.GOU.Y60;

% Endpoint cumulative-ROS survival curves at the common 60-min horizon.
m_max = max([Mb + 3, ceil(quantile([Y60_c; Y60_g], 0.997))]);
m_grid = 0:m_max;

tail_c = arrayfun(@(m) mean(Y60_c >= m), m_grid);
tail_g = arrayfun(@(m) mean(Y60_g >= m), m_grid);

se_tail_c = sqrt(tail_c .* (1 - tail_c) / cfg.n_eval);
se_tail_g = sqrt(tail_g .* (1 - tail_g) / cfg.n_eval);

tail_c_lo = max(tail_c - 1.96 .* se_tail_c, 0.5 / cfg.n_eval);
tail_c_hi = min(tail_c + 1.96 .* se_tail_c, 1);
tail_g_lo = max(tail_g - 1.96 .* se_tail_g, 0.5 / cfg.n_eval);
tail_g_hi = min(tail_g + 1.96 .* se_tail_g, 1);

Pfull_c = mean(Y60_c >= Mb);
Pfull_g = mean(Y60_g >= Mb);
dP_full = Pfull_g - Pfull_c;

% FPT distributions and RMST over the same complete 0-60 min horizon.
t_fpt = linspace(0, sim.T_total, 241)';
F_c = empirical_fpt_cdf(R.CIR.FPT, t_fpt);
F_g = empirical_fpt_cdf(R.GOU.FPT, t_fpt);

se_F_c = sqrt(F_c .* (1 - F_c) / cfg.n_eval);
se_F_g = sqrt(F_g .* (1 - F_g) / cfg.n_eval);

F_c_lo = max(F_c - 1.96 .* se_F_c, 0);
F_c_hi = min(F_c + 1.96 .* se_F_c, 1);
F_g_lo = max(F_g - 1.96 .* se_F_g, 0);
F_g_hi = min(F_g + 1.96 .* se_F_g, 1);

F60_c = mean(~isnan(R.CIR.FPT) & R.CIR.FPT <= sim.T_total);
F60_g = mean(~isnan(R.GOU.FPT) & R.GOU.FPT <= sim.T_total);
dF60  = F60_g - F60_c;

RMST_c = rmst_from_fpt(R.CIR.FPT, sim.T_total);
RMST_g = rmst_from_fpt(R.GOU.FPT, sim.T_total);
dRMST  = RMST_c - RMST_g;  % positive means GammaOU crosses earlier

% For monotone cumulative Y(t), endpoint exceedance and terminal FPT are
% the same event when the threshold and observation horizon are identical.
identity_error_c = abs(Pfull_c - F60_c);
identity_error_g = abs(Pfull_g - F60_g);

rng(seed_boot, 'twister');
BOOT = bootstrap_high_order_effects( ...
    Y60_c, Y60_g, Mb, R.CIR.FPT, R.GOU.FPT, sim.T_total, cfg.n_boot);

CI_dP_full = prctile(BOOT.dP_full, [2.5 97.5]);
CI_dF60    = prctile(BOOT.dF60,   [2.5 97.5]);
CI_dRMST   = prctile(BOOT.dRMST,  [2.5 97.5]);

fprintf('Main-window low-order diagnostics:\n');
fprintf('  max symmetric relative mean difference = %.3f %%\n', ...
    100 * max_rel_mean_CG);
fprintf('  max symmetric relative Fano difference = %.3f %%\n\n', ...
    100 * max_rel_fano_CG);

fprintf('Full-horizon high-order risk summary:\n');
fprintf('  accumulation/FPT horizon = [0, %.1f] min\n', sim.T_total/60);
fprintf('  P_CIR(Y60 >= M_b)     = %.5f\n', Pfull_c);
fprintf('  P_GammaOU(Y60 >= M_b) = %.5f\n', Pfull_g);
fprintf('  DeltaP_60             = %.5f [%.5f, %.5f]\n', ...
    dP_full, CI_dP_full(1), CI_dP_full(2));
fprintf('  DeltaF_60             = %.5f [%.5f, %.5f]\n', ...
    dF60, CI_dF60(1), CI_dF60(2));
fprintf('  DeltaRMST_60          = %.3f min [%.3f, %.3f]\n', ...
    dRMST/60, CI_dRMST(1)/60, CI_dRMST(2)/60);
fprintf('  endpoint-FPT identity error CIR/GOU = %.3e / %.3e\n\n', ...
    identity_error_c, identity_error_g);

%% =========================================================================
% 8. Long-time full-state simulation for SI
% =========================================================================
rng(seed_long, 'twister');
fprintf('Long-time SI diagnostic run starts...\n');
LONG = simulate_moments_only(p, sim.dt, cfg.T_long, cfg.rec_dt_long, cfg.n_long);
fprintf('Long-time SI diagnostic run finished.\n\n');

fano_c_long_s = movmean(LONG.CIR.fano, cfg.smooth_long, 'omitnan');
fano_g_long_s = movmean(LONG.GOU.fano, cfg.smooth_long, 'omitnan');

fano_LNA_long = finite_time_fano( ...
    LONG.t, p.v_eff, p.F_int, p.A1, p.kappa1, p.A2, p.kappa2);

tail_idx_long = LONG.t >= 0.80 * cfg.T_long;
plateau_values = [fano_c_long_s(tail_idx_long); fano_g_long_s(tail_idx_long)];
plateau_values = plateau_values(isfinite(plateau_values));
F_emp_plateau = mean(plateau_values);

%% =========================================================================
% 9. Plot main Fig.2A-C
% =========================================================================
c_cir = hex2rgb('#4DBBD5');
c_gou = hex2rgb('#E64B35');
c_blk = [0.10 0.10 0.10];
c_gray = [0.45 0.45 0.45];
c_light = [0.92 0.92 0.92];

if make_main_figure
    figMain = figure('Color', 'w', 'Position', [50 50 1900 900]);

    % ---------------- Panel A: path architecture ----------------
    axA1 = axes(figMain, 'Position', [0.055 0.675 0.250 0.205]);
    hold(axA1, 'on');
    for j = 1:size(ENV.path_CIR, 2)
        plot(axA1, ENV.path_t/60, ENV.path_CIR(:,j)/p.mu_c1, ...
            'Color', blend_with_white(c_cir, 0.70), 'LineWidth', 0.75);
    end
    set(axA1, 'YScale', 'log');
    xlim(axA1, [0 sim.T_total/60]);
    ylim(axA1, ENV.path_ylim);
    ylabel(axA1, 'CIR: c_1(t)/\mu_{c1}');
    set(axA1, 'XTickLabel', []);
    title(axA1, 'Continuous diffusive paths');
    format_axes(axA1);

    axA2 = axes(figMain, 'Position', [0.055 0.485 0.250 0.160]);
    hold(axA2, 'on');
    for j = 1:size(ENV.path_GOU, 2)
        plot(axA2, ENV.path_t/60, ENV.path_GOU(:,j)/p.mu_c1, ...
            'Color', blend_with_white(c_gou, 0.68), 'LineWidth', 0.75);
    end
    set(axA2, 'YScale', 'log');
    xlim(axA2, [0 sim.T_total/60]);
    ylim(axA2, ENV.path_ylim);
    xlabel(axA2, 'Time (min)');
    ylabel(axA2, 'GammaOU: c_1(t)/\mu_{c1}');
    title(axA2, 'Jump-intermittent paths');
    format_axes(axA2);

    % ---------------- Panel A: stationary marginal ----------------
    axA3 = axes(figMain, 'Position', [0.055 0.105 0.115 0.275]);
    hold(axA3, 'on');
    plot(axA3, ENV.z_centers, ENV.pdf_CIR, 'o', ...
        'Color', c_cir, 'MarkerSize', 3.5, 'LineWidth', 0.9);
    plot(axA3, ENV.z_centers, ENV.pdf_GOU, 's', ...
        'Color', c_gou, 'MarkerSize', 3.2, 'LineWidth', 0.9);
    plot(axA3, ENV.z_grid, ENV.pdf_theory, '-', ...
        'Color', c_blk, 'LineWidth', 2.0);
    xlabel(axA3, 'Normalized rate c_1/\mu_{c1}');
    ylabel(axA3, 'Density');
    title(axA3, 'Matched stationary marginal');
    xlim(axA3, [0 ENV.z_plot_max]);
    legend(axA3, {'CIR','GammaOU','Gamma theory'}, ...
        'Box','off','Location','northeast','FontSize',8);
    format_axes(axA3);

    % ---------------- Panel A: memory ----------------
    axA4 = axes(figMain, 'Position', [0.190 0.105 0.115 0.275]);
    hold(axA4, 'on');
    plot(axA4, ENV.lags/60, ENV.acf_CIR, 'o-', ...
        'Color', c_cir, 'LineWidth', 1.5, 'MarkerSize', 4);
    plot(axA4, ENV.lags/60, ENV.acf_GOU, 's--', ...
        'Color', c_gou, 'LineWidth', 1.5, 'MarkerSize', 4);
    plot(axA4, ENV.lags/60, ENV.acf_theory, ':', ...
        'Color', c_blk, 'LineWidth', 2.0);
    xlabel(axA4, 'Lag (min)');
    ylabel(axA4, 'Corr[c_1(t),c_1(t+\tau)]');
    title(axA4, 'Matched memory');
    xlim(axA4, [0 max(ENV.lags)/60]);
    ylim(axA4, [-0.05 1.05]);
    legend(axA4, {'CIR','GammaOU','e^{-\kappa_1\tau}'}, ...
        'Box','off','Location','northeast','FontSize',8);
    format_axes(axA4);

    % ---------------- Panel B: mean ----------------
    axB1 = axes(figMain, 'Position', [0.365 0.545 0.270 0.335]);
    hold(axB1, 'on');
    plot(axB1, t_main/60, R.CIR.meanY, ...
        'Color', c_cir, 'LineWidth', 2.5);
    plot(axB1, t_main/60, R.GOU.meanY, '--', ...
        'Color', c_gou, 'LineWidth', 2.3);
    plot(axB1, t_main/60, mean_LNA, ':', ...
        'Color', c_blk, 'LineWidth', 2.2);
    xlabel(axB1, 'Time (min)');
    ylabel(axB1, 'Mean cumulative ROS, E[Y_t]');
    title(axB1, 'Full-state mean versus shared closure');
    legend(axB1, {'CIR full TD-CME','GammaOU full TD-CME', ...
        'Shared first-moment closure'}, ...
        'Box','off','Location','northwest');
    xlim(axB1, [0 sim.T_total/60]);
    format_axes(axB1);

    % ---------------- Panel B: Fano ----------------
    axB2 = axes(figMain, 'Position', [0.365 0.105 0.270 0.335]);
    hold(axB2, 'on');
    plot(axB2, t_main/60, fano_c_s, ...
        'Color', c_cir, 'LineWidth', 2.3);
    plot(axB2, t_main/60, fano_g_s, '--', ...
        'Color', c_gou, 'LineWidth', 2.2);
    plot(axB2, t_main/60, fano_LNA, ':', ...
        'Color', c_blk, 'LineWidth', 2.2);
    xlabel(axB2, 'Time (min)');
    ylabel(axB2, 'Fano factor Var(Y_t)/E[Y_t]');
    title(axB2, 'Low-order fluctuations remain closely aligned');
    legend(axB2, {'CIR full TD-CME','GammaOU full TD-CME', ...
        'Shared finite-time LNA'}, ...
        'Box','off','Location','southeast');
    xlim(axB2, [0 sim.T_total/60]);
    format_axes(axB2);

    ymaxB2 = max([fano_c_s(:); fano_g_s(:); fano_LNA(:)], [], 'omitnan');
    text(axB2, 0.04*sim.T_total/60, 0.91*ymaxB2, ...
        sprintf('max CIR-GammaOU difference\nmean %.2f%%; Fano %.2f%%', ...
        100*max_rel_mean_CG, 100*max_rel_fano_CG), ...
        'FontName','Arial','FontSize',10,'Color',c_blk, ...
        'VerticalAlignment','top');

    % ---------------- Panel C: full-horizon cumulative-ROS tail ----------------
    axC1 = axes(figMain, 'Position', [0.690 0.545 0.270 0.335]);
    hold(axC1, 'on');
    fill(axC1, [m_grid fliplr(m_grid)], ...
        [tail_c_lo fliplr(tail_c_hi)], ...
        c_cir, 'FaceAlpha', 0.10, 'EdgeColor', 'none', ...
        'HandleVisibility','off');
    fill(axC1, [m_grid fliplr(m_grid)], ...
        [tail_g_lo fliplr(tail_g_hi)], ...
        c_gou, 'FaceAlpha', 0.10, 'EdgeColor', 'none', ...
        'HandleVisibility','off');
    plot(axC1, m_grid, tail_c, ...
        'Color', c_cir, 'LineWidth', 2.5);
    plot(axC1, m_grid, tail_g, '--', ...
        'Color', c_gou, 'LineWidth', 2.4);
    xline(axC1, Mb, ':', 'Color', c_blk, 'LineWidth', 1.5);
    set(axC1, 'YScale', 'log');
    ylim(axC1, [max(0.5/cfg.n_eval, 1e-4), 1]);
    xlabel(axC1, 'Cumulative ROS at 60 min, Y(60)');
    ylabel(axC1, 'Tail probability P(Y(60) \geq m)');
    title(axC1, 'Full-horizon ROS-load tail separates');
    legend(axC1, {'CIR','GammaOU'}, ...
        'Box','off','Location','southwest');
    format_axes(axC1);

    txtC1 = sprintf('M_b=%d\n\\DeltaP_{60}=%.4f\n95%% CI [%.4f, %.4f]', ...
        Mb, dP_full, CI_dP_full(1), CI_dP_full(2));
    text(axC1, 0.97*m_max, 0.18, txtC1, ...
        'HorizontalAlignment','right','VerticalAlignment','top', ...
        'FontName','Arial','FontSize',10,'Color',c_blk);

    % ---------------- Panel C: FPT ----------------
    axC2 = axes(figMain, 'Position', [0.690 0.105 0.270 0.335]);
    hold(axC2, 'on');
    fill(axC2, [t_fpt; flipud(t_fpt)]/60, ...
        [F_c_lo; flipud(F_c_hi)], ...
        c_cir, 'FaceAlpha', 0.10, 'EdgeColor', 'none', ...
        'HandleVisibility','off');
    fill(axC2, [t_fpt; flipud(t_fpt)]/60, ...
        [F_g_lo; flipud(F_g_hi)], ...
        c_gou, 'FaceAlpha', 0.10, 'EdgeColor', 'none', ...
        'HandleVisibility','off');
    plot(axC2, t_fpt/60, F_c, ...
        'Color', c_cir, 'LineWidth', 2.5);
    plot(axC2, t_fpt/60, F_g, '--', ...
        'Color', c_gou, 'LineWidth', 2.4);
    xlabel(axC2, 'Time to first ROS-burden crossing (min)');
    ylabel(axC2, sprintf('F_{T_{M_b}}(t), M_b=%d', Mb));
    title(axC2, 'ROS-burden first-passage risk separates');
    legend(axC2, {'CIR','GammaOU'}, ...
        'Box','off','Location','southeast');
    xlim(axC2, [0 sim.T_total/60]);
    ymaxC2 = max([F_c; F_g], [], 'omitnan');
    ylim(axC2, [0 max(0.15, 1.15*ymaxC2)]);
    format_axes(axC2);

    txtC2 = sprintf(['\\DeltaF_{60}=%.4f [%.4f, %.4f]\n' ...
        '\\DeltaRMST_{60}=%.2f min [%.2f, %.2f]'], ...
        dF60, CI_dF60(1), CI_dF60(2), ...
        dRMST/60, CI_dRMST(1)/60, CI_dRMST(2)/60);
    text(axC2, 0.04*sim.T_total/60, 0.84*max(0.15, 1.15*ymaxC2), ...
        txtC2, 'FontName','Arial','FontSize',10, ...
        'Color',c_blk,'VerticalAlignment','top');

    % Column labels
    annotation(figMain, 'textbox', [0.055 0.915 0.250 0.045], ...
        'String', 'A  Matched stochastic baseline', ...
        'FontName','Arial','FontSize',15,'FontWeight','bold', ...
        'EdgeColor','none','HorizontalAlignment','left');
    annotation(figMain, 'textbox', [0.365 0.915 0.270 0.045], ...
        'String', 'B  Low-order approximate consistency', ...
        'FontName','Arial','FontSize',15,'FontWeight','bold', ...
        'EdgeColor','none','HorizontalAlignment','left');
    annotation(figMain, 'textbox', [0.690 0.915 0.270 0.045], ...
        'String', 'C  High-order ROS-risk separation', ...
        'FontName','Arial','FontSize',15,'FontWeight','bold', ...
        'EdgeColor','none','HorizontalAlignment','left');

    if export_figures
        export_figure_safe(figMain, ...
            sprintf('Fig2ABC_TDCME_%s', lower(run_level)));
    end
end

%% =========================================================================
% 10. Plot SI low-order closure/LNA companion
% =========================================================================
if make_SI_figure
    figSI = figure('Color','w','Position',[80 80 1500 900]);
    tlSI = tiledlayout(figSI, 2, 2, ...
        'TileSpacing','compact','Padding','compact');

    % SI-A: shared analytic prediction
    axS1 = nexttile(tlSI, 1);
    yyaxis(axS1, 'left');
    plot(axS1, t_main/60, mean_LNA, ...
        'Color', c_blk, 'LineWidth', 2.4);
    ylabel(axS1, 'Shared mean closure, E[Y_t]');
    yyaxis(axS1, 'right');
    plot(axS1, t_main/60, fano_LNA, '--', ...
        'Color', c_gray, 'LineWidth', 2.2);
    ylabel(axS1, 'Shared finite-time LNA Fano');
    xlabel(axS1, 'Time (min)');
    title(axS1, 'A. Matched marginals and memory give one low-order prediction');
    xlim(axS1, [0 sim.T_total/60]);
    format_axes(axS1);

    % SI-B: empirical model-model deviations
    axS2 = nexttile(tlSI, 2);
    hold(axS2, 'on');
    plot(axS2, t_main(valid_low)/60, ...
        100*rel_mean_CG(valid_low), ...
        'Color', c_cir, 'LineWidth', 2.3);
    plot(axS2, t_main(valid_low)/60, ...
        100*rel_fano_CG(valid_low), ...
        'Color', c_gou, 'LineWidth', 2.3);
    yline(axS2, 5, ':', 'Color', c_gray, 'LineWidth', 1.4);
    xlabel(axS2, 'Time (min)');
    ylabel(axS2, 'Symmetric relative difference (%)');
    title(axS2, 'B. Full-state low-order deviations in the main window');
    legend(axS2, {'Mean cumulative ROS','Fano factor','5% guide'}, ...
        'Box','off','Location','northwest');
    xlim(axS2, [0 sim.T_total/60]);
    ymaxS2 = max(6, 1.18*100*max([max_rel_mean_CG, max_rel_fano_CG]));
    ylim(axS2, [0 ymaxS2]);
    text(axS2, 0.97*sim.T_total/60, 0.92*ymaxS2, ...
        sprintf('max mean = %.2f%%\nmax Fano = %.2f%%', ...
        100*max_rel_mean_CG, 100*max_rel_fano_CG), ...
        'HorizontalAlignment','right','VerticalAlignment','top', ...
        'FontName','Arial','FontSize',10);
    format_axes(axS2);

    % SI-C: long-time Fano diagnostic
    axS3 = nexttile(tlSI, 3);
    hold(axS3, 'on');
    ymaxS3 = 1.12 * max([fano_c_long_s(:); fano_g_long_s(:); ...
        fano_LNA_long(:); F_emp_plateau], [], 'omitnan');
    patch(axS3, [0 60 60 0], [0 0 ymaxS3 ymaxS3], ...
        c_light, 'EdgeColor','none','FaceAlpha',0.70, ...
        'HandleVisibility','off');
    plot(axS3, LONG.t/60, fano_c_long_s, ...
        'Color', c_cir, 'LineWidth', 2.1);
    plot(axS3, LONG.t/60, fano_g_long_s, '--', ...
        'Color', c_gou, 'LineWidth', 2.1);
    plot(axS3, LONG.t/60, fano_LNA_long, ':', ...
        'Color', c_blk, 'LineWidth', 2.2);
    yline(axS3, F_emp_plateau, '-.', ...
        'Color', c_gray, 'LineWidth', 1.6);
    xlabel(axS3, 'Time (min)');
    ylabel(axS3, 'Fano factor');
    title(axS3, 'C. Long-time full-state Fano versus the LNA reference');
    legend(axS3, {'CIR full TD-CME','GammaOU full TD-CME', ...
        'Finite-time LNA','Empirical late-time plateau'}, ...
        'Box','off','Location','southeast');
    xlim(axS3, [0 cfg.T_long/60]);
    ylim(axS3, [0 ymaxS3]);
    format_axes(axS3);

    % SI-D: shared decomposition
    axS4 = nexttile(tlSI, 4);
    hold(axS4, 'on');
    comp = [p.F_int, p.F_ext1, p.F_ext2];
    bh = barh(axS4, 1, comp, 0.35, 'stacked');
    bh(1).FaceColor = [0.62 0.62 0.62];
    bh(2).FaceColor = [0.93 0.69 0.13];
    bh(3).FaceColor = [0.49 0.18 0.56];
    for j = 1:numel(bh)
        bh(j).EdgeColor = 'none';
    end
    xline(axS4, p.F_LNA_long, ':', ...
        'Color', c_blk, 'LineWidth', 2.2);
    set(axS4, 'YTick', 1, ...
        'YTickLabel', {'Shared CIR/GammaOU LNA'});
    xlabel(axS4, 'Long-time Fano contribution');
    title(axS4, 'D. Shared closure/LNA noise decomposition');
    legend(axS4, {'Intrinsic reaction noise', ...
        'Extrinsic c_1 contribution', ...
        'Extrinsic c_2 contribution', ...
        'Total LNA Fano'}, ...
        'Box','off','Location','southoutside', ...
        'Orientation','horizontal');
    xlim(axS4, [0 1.12*p.F_LNA_long]);
    ylim(axS4, [0.55 1.45]);
    text(axS4, 0.98*p.F_LNA_long, 1.22, ...
        sprintf('F_{total}=%.2f', p.F_LNA_long), ...
        'HorizontalAlignment','right','FontSize',10);
    format_axes(axS4);

    if export_figures
        export_figure_safe(figSI, ...
            sprintf('SI_TDCME_LowOrder_Closure_LNA_%s', lower(run_level)));
    end
end

%% =========================================================================
% 11. Export summaries and compact result object
% =========================================================================
summaryMain = table( ...
    string(run_level), p.N, p.kappa1, p.shape1, p.mu_c1, p.theta1, ...
    p.kappa2, p.shape2, p.mu_c2, p.theta2, ...
    sim.dt, sim.T_total, ...
    cfg.n_cal, cfg.n_eval, Mb, ...
    max_rel_mean_CG, max_rel_fano_CG, ...
    Pfull_c, Pfull_g, dP_full, CI_dP_full(1), CI_dP_full(2), ...
    F60_c, F60_g, dF60, CI_dF60(1), CI_dF60(2), ...
    RMST_c/60, RMST_g/60, dRMST/60, CI_dRMST(1)/60, CI_dRMST(2)/60, ...
    'VariableNames', { ...
    'run_level','N','kappa1','shape1','mu_c1','theta1', ...
    'kappa2','shape2','mu_c2','theta2', ...
    'dt_s','T_total_s', ...
    'n_cal','n_eval','Mb', ...
    'max_rel_mean_CIR_GOU','max_rel_fano_CIR_GOU', ...
    'Pfull60_CIR','Pfull60_GOU','DeltaP60','DeltaP60_CIlo','DeltaP60_CIhi', ...
    'F60_CIR','F60_GOU','DeltaF60','DeltaF60_CIlo','DeltaF60_CIhi', ...
    'RMST60_CIR_min','RMST60_GOU_min','DeltaRMST60_min', ...
    'DeltaRMST60_CIlo_min','DeltaRMST60_CIhi_min'});

summarySI = table( ...
    p.v_eff, p.F_int, p.F_ext1, p.F_ext2, p.F_ext, p.F_LNA_long, ...
    F_emp_plateau, ...
    max(rel_mean_c_theory(valid_low), [], 'omitnan'), ...
    max(rel_mean_g_theory(valid_low), [], 'omitnan'), ...
    max(rel_fano_c_theory(valid_low), [], 'omitnan'), ...
    max(rel_fano_g_theory(valid_low), [], 'omitnan'), ...
    'VariableNames', { ...
    'mean_flux_closure','F_intrinsic','F_ext_c1','F_ext_c2', ...
    'F_ext_total','F_LNA_long','F_empirical_plateau', ...
    'max_rel_mean_CIR_vs_closure','max_rel_mean_GOU_vs_closure', ...
    'max_rel_fano_CIR_vs_LNA','max_rel_fano_GOU_vs_LNA'});

inputAudit = table( ...
    ["c1";"c2"], ...
    [p.mu_c1;p.mu_c2], ...
    [p.var_c1;p.var_c2], ...
    [p.kappa1;p.kappa2], ...
    [1/p.kappa1;1/p.kappa2], ...
    [p.shape1;p.shape2], ...
    [p.theta1;p.theta2], ...
    'VariableNames', {'channel','stationary_mean','stationary_variance', ...
    'kappa_per_s','memory_time_s','Gamma_shape','Gamma_scale'});

if export_tables
    writetable(summaryMain, ...
        sprintf('Fig2ABC_TDCME_summary_%s.csv', lower(run_level)));
    writetable(summarySI, ...
        sprintf('SI_TDCME_LowOrder_summary_%s.csv', lower(run_level)));
    writetable(inputAudit, ...
        sprintf('Fig2A_input_matching_audit_%s.csv', lower(run_level)));

    curveMain = table( ...
        t_main/60, R.CIR.meanY, R.GOU.meanY, mean_LNA, ...
        fano_c_s, fano_g_s, fano_LNA, ...
        'VariableNames', {'time_min','mean_CIR','mean_GOU','mean_closure', ...
        'Fano_CIR','Fano_GOU','Fano_LNA'});
    writetable(curveMain, ...
        sprintf('Fig2B_loworder_curves_%s.csv', lower(run_level)));

    tailTable = table( ...
        m_grid(:), tail_c(:), tail_g(:), ...
        tail_c_lo(:), tail_c_hi(:), tail_g_lo(:), tail_g_hi(:), ...
        'VariableNames', {'cumulative_ROS_Y60','tail_CIR','tail_GOU', ...
        'tail_CIR_CIlo','tail_CIR_CIhi','tail_GOU_CIlo','tail_GOU_CIhi'});
    writetable(tailTable, ...
        sprintf('Fig2C_fullhorizon_tail_%s.csv', lower(run_level)));

    fptTable = table( ...
        t_fpt/60, F_c, F_g, F_c_lo, F_c_hi, F_g_lo, F_g_hi, ...
        'VariableNames', {'time_min','FPT_CDF_CIR','FPT_CDF_GOU', ...
        'CIR_CIlo','CIR_CIhi','GOU_CIlo','GOU_CIhi'});
    writetable(fptTable, ...
        sprintf('Fig2C_FPT_curves_%s.csv', lower(run_level)));
end

if save_results_mat
    SAVE = struct();
    SAVE.parameters = p;
    SAVE.simulation = sim;
    SAVE.configuration = cfg;
    SAVE.thresholds.Mb = Mb;
    SAVE.ENV = ENV;
    SAVE.R = R;
    SAVE.LONG = LONG;
    SAVE.low_order.mean_LNA = mean_LNA;
    SAVE.low_order.fano_LNA = fano_LNA;
    SAVE.high_order.Y60_grid = m_grid;
    SAVE.high_order.tail_CIR = tail_c;
    SAVE.high_order.tail_GOU = tail_g;
    SAVE.high_order.t_fpt = t_fpt;
    SAVE.high_order.F_CIR = F_c;
    SAVE.high_order.F_GOU = F_g;
    SAVE.bootstrap = BOOT;
    SAVE.summaryMain = summaryMain;
    SAVE.summarySI = summarySI;

    save(sprintf('TDCME_Fig2ABC_results_%s.mat', lower(run_level)), ...
        'SAVE', '-v7.3');
end

fprintf('Finished. Main conclusions must be judged from the final run,\n');
fprintf('not from smoke-level Monte Carlo fluctuations.\n');

%% =========================================================================
% Local functions
% =========================================================================

function validate_time_grid(sim, cfg)
    vals = [sim.T_total, cfg.rec_dt_main, cfg.rec_dt_long, cfg.T_long];
    for k = 1:numel(vals)
        if abs(vals(k)/sim.dt - round(vals(k)/sim.dt)) > 1e-10
            error('All output and analysis times must be divisible by dt.');
        end
    end
end

function ENV = simulate_environment_audit(p, sim, cfg)
    % Representative path ensemble for c1.
    n_show = cfg.n_show_paths;
    n_steps = round(sim.T_total / sim.dt);
    rec_stride = round(20 / sim.dt);
    rec_idx = 0:rec_stride:n_steps;
    t_rec = rec_idx(:) .* sim.dt;

    c0 = gamrnd(p.shape1, p.theta1, n_show, 1);
    c_c = c0;
    c_g = c0;

    path_c = zeros(numel(rec_idx), n_show);
    path_g = zeros(numel(rec_idx), n_show);
    path_c(1,:) = c_c.';
    path_g(1,:) = c_g.';

    ir = 2;
    for it = 1:n_steps
        c_c = cir_exact_transition( ...
            c_c, p.kappa1, p.shape1, p.theta1, sim.dt);
        [c_g, ~] = gamma_ou_step_with_integral( ...
            c_g, p.a1, p.kappa1, p.theta1, sim.dt);

        if ir <= numel(rec_idx) && it == rec_idx(ir)
            path_c(ir,:) = c_c.';
            path_g(ir,:) = c_g.';
            ir = ir + 1;
        end
    end

    all_norm = [path_c(:); path_g(:)] ./ p.mu_c1;
    positive_norm = all_norm(all_norm > 0);
    ylo = max(1e-4, quantile(positive_norm, 0.002));
    yhi = max(10, quantile(positive_norm, 0.995));
    ENV.path_ylim = [ylo yhi];

    ENV.path_t = t_rec;
    ENV.path_CIR = path_c;
    ENV.path_GOU = path_g;

    % Stationary-marginal audit after a 60-min direct transition.
    c0_env = gamrnd(p.shape1, p.theta1, cfg.n_env, 1);
    cT_c = cir_exact_transition( ...
        c0_env, p.kappa1, p.shape1, p.theta1, sim.T_total);
    [cT_g, ~] = gamma_ou_step_with_integral( ...
        c0_env, p.a1, p.kappa1, p.theta1, sim.T_total);

    z_c = cT_c ./ p.mu_c1;
    z_g = cT_g ./ p.mu_c1;

    z_plot_max = max(8, quantile([z_c; z_g], 0.995));
    edges = linspace(0, z_plot_max, 42);
    centers = 0.5 .* (edges(1:end-1) + edges(2:end));

    pdf_c = histcounts(z_c, edges, 'Normalization','pdf');
    pdf_g = histcounts(z_g, edges, 'Normalization','pdf');

    z_grid = linspace(0.01, z_plot_max, 400);
    k = p.shape1;
    pdf_theory = (k^k / gamma(k)) .* z_grid.^(k-1) .* exp(-k .* z_grid);

    ENV.z_centers = centers;
    ENV.pdf_CIR = pdf_c;
    ENV.pdf_GOU = pdf_g;
    ENV.z_grid = z_grid;
    ENV.pdf_theory = pdf_theory;
    ENV.z_plot_max = z_plot_max;

    % Ensemble ACF audit using common stationary starting states.
    lags = [0 300 600 900 1200 1800 2400 3000 3600]';
    c0_acf = gamrnd(p.shape1, p.theta1, cfg.n_env, 1);

    acf_c = zeros(size(lags));
    acf_g = zeros(size(lags));

    acf_c(1) = 1;
    acf_g(1) = 1;

    for j = 2:numel(lags)
        lag = lags(j);

        c_lag_c = cir_exact_transition( ...
            c0_acf, p.kappa1, p.shape1, p.theta1, lag);
        [c_lag_g, ~] = gamma_ou_step_with_integral( ...
            c0_acf, p.a1, p.kappa1, p.theta1, lag);

        cc = corrcoef(c0_acf, c_lag_c);
        gg = corrcoef(c0_acf, c_lag_g);

        acf_c(j) = cc(1,2);
        acf_g(j) = gg(1,2);
    end

    ENV.lags = lags;
    ENV.acf_CIR = acf_c;
    ENV.acf_GOU = acf_g;
    ENV.acf_theory = exp(-p.kappa1 .* lags);

    % Printed empirical audit.
    fprintf('c1 input audit after direct stationary transition:\n');
    fprintf('  mean CIR/GOU/theory = %.6e / %.6e / %.6e\n', ...
        mean(cT_c), mean(cT_g), p.mu_c1);
    fprintf('  var  CIR/GOU/theory = %.6e / %.6e / %.6e\n', ...
        var(cT_c,1), var(cT_g,1), p.var_c1);
    fprintf('  max |empirical ACF - theory| CIR/GOU = %.4f / %.4f\n', ...
        max(abs(acf_c - ENV.acf_theory)), ...
        max(abs(acf_g - ENV.acf_theory)));
end

function OUT = simulate_calibration(p, sim, n)
    % Independent calibration over the complete 0-60 min horizon.
    n_steps = round(sim.T_total / sim.dt);

    x1_0 = binornd(p.N, p.p_x1_ss, n, 1);
    c1_0 = gamrnd(p.shape1, p.theta1, n, 1);
    c2_0 = gamrnd(p.shape2, p.theta2, n, 1);

    x1_c = x1_0;
    x2_c = p.N - x1_c;
    y_c  = zeros(n,1);
    c1_c = c1_0;
    c2_c = c2_0;

    x1_g = x1_0;
    x2_g = p.N - x1_g;
    y_g  = zeros(n,1);
    c1_g = c1_0;
    c2_g = c2_0;

    for it = 1:n_steps
        % CIR
        c1_old = c1_c;
        c2_old = c2_c;

        c1_c = cir_exact_transition( ...
            c1_c, p.kappa1, p.shape1, p.theta1, sim.dt);
        c2_c = cir_exact_transition( ...
            c2_c, p.kappa2, p.shape2, p.theta2, sim.dt);

        I1_c = 0.5 .* (c1_old + c1_c) .* sim.dt;
        I2_c = 0.5 .* (c2_old + c2_c) .* sim.dt;

        pr1_c = min(max(1 - exp(-I1_c), 0), 1);
        pr2_c = min(max(1 - exp(-I2_c), 0), 1);

        r1_c = binornd(x1_c, pr1_c);
        r2_c = binornd(x2_c, pr2_c);

        x1_c = x1_c - r1_c + r2_c;
        x2_c = p.N - x1_c;
        y_c  = y_c + r1_c;

        % GammaOU
        [c1_g, I1_g] = gamma_ou_step_with_integral( ...
            c1_g, p.a1, p.kappa1, p.theta1, sim.dt);
        [c2_g, I2_g] = gamma_ou_step_with_integral( ...
            c2_g, p.a2, p.kappa2, p.theta2, sim.dt);

        pr1_g = min(max(1 - exp(-I1_g), 0), 1);
        pr2_g = min(max(1 - exp(-I2_g), 0), 1);

        r1_g = binornd(x1_g, pr1_g);
        r2_g = binornd(x2_g, pr2_g);

        x1_g = x1_g - r1_g + r2_g;
        x2_g = p.N - x1_g;
        y_g  = y_g + r1_g;
    end

    OUT.CIR.Y60 = y_c;
    OUT.GOU.Y60 = y_g;
end

function OUT = simulate_full_evaluation(p, sim, cfg, n, Mb)
    n_steps = round(sim.T_total / sim.dt);

    rec_stride = round(cfg.rec_dt_main / sim.dt);
    rec_idx = 0:rec_stride:n_steps;
    t_rec = rec_idx(:) .* sim.dt;

    x1_0 = binornd(p.N, p.p_x1_ss, n, 1);
    c1_0 = gamrnd(p.shape1, p.theta1, n, 1);
    c2_0 = gamrnd(p.shape2, p.theta2, n, 1);

    x1_c = x1_0;
    x2_c = p.N - x1_c;
    y_c  = zeros(n,1);
    c1_c = c1_0;
    c2_c = c2_0;

    x1_g = x1_0;
    x2_g = p.N - x1_g;
    y_g  = zeros(n,1);
    c1_g = c1_0;
    c2_g = c2_0;

    FPT_c = nan(n,1);
    FPT_g = nan(n,1);

    mean_c = zeros(numel(rec_idx),1);
    mean_g = zeros(numel(rec_idx),1);
    var_c  = zeros(numel(rec_idx),1);
    var_g  = zeros(numel(rec_idx),1);

    mean_c(1) = 0;
    mean_g(1) = 0;
    var_c(1) = 0;
    var_g(1) = 0;

    ir = 2;

    for it = 1:n_steps
        t_now = it .* sim.dt;

        % CIR
        c1_old = c1_c;
        c2_old = c2_c;

        c1_c = cir_exact_transition( ...
            c1_c, p.kappa1, p.shape1, p.theta1, sim.dt);
        c2_c = cir_exact_transition( ...
            c2_c, p.kappa2, p.shape2, p.theta2, sim.dt);

        I1_c = 0.5 .* (c1_old + c1_c) .* sim.dt;
        I2_c = 0.5 .* (c2_old + c2_c) .* sim.dt;

        pr1_c = min(max(1 - exp(-I1_c), 0), 1);
        pr2_c = min(max(1 - exp(-I2_c), 0), 1);

        r1_c = binornd(x1_c, pr1_c);
        r2_c = binornd(x2_c, pr2_c);

        x1_c = x1_c - r1_c + r2_c;
        x2_c = p.N - x1_c;
        y_c  = y_c + r1_c;

        new_c = isnan(FPT_c) & y_c >= Mb;
        FPT_c(new_c) = t_now;

        % GammaOU
        [c1_g, I1_g] = gamma_ou_step_with_integral( ...
            c1_g, p.a1, p.kappa1, p.theta1, sim.dt);
        [c2_g, I2_g] = gamma_ou_step_with_integral( ...
            c2_g, p.a2, p.kappa2, p.theta2, sim.dt);

        pr1_g = min(max(1 - exp(-I1_g), 0), 1);
        pr2_g = min(max(1 - exp(-I2_g), 0), 1);

        r1_g = binornd(x1_g, pr1_g);
        r2_g = binornd(x2_g, pr2_g);

        x1_g = x1_g - r1_g + r2_g;
        x2_g = p.N - x1_g;
        y_g  = y_g + r1_g;

        new_g = isnan(FPT_g) & y_g >= Mb;
        FPT_g(new_g) = t_now;

        if ir <= numel(rec_idx) && it == rec_idx(ir)
            mean_c(ir) = mean(y_c);
            mean_g(ir) = mean(y_g);
            var_c(ir)  = var(y_c,0);
            var_g(ir)  = var(y_g,0);
            ir = ir + 1;
        end

        if mod(it, round(600/sim.dt)) == 0
            fprintf('  evaluation t = %.1f min\n', t_now/60);
        end
    end

    fano_c = var_c ./ max(mean_c, eps);
    fano_g = var_g ./ max(mean_g, eps);
    fano_c(1) = nan;
    fano_g(1) = nan;

    OUT.t = t_rec;

    OUT.CIR.meanY = mean_c;
    OUT.CIR.varY  = var_c;
    OUT.CIR.fano  = fano_c;
    OUT.CIR.Y60   = y_c;
    OUT.CIR.FPT   = FPT_c;

    OUT.GOU.meanY = mean_g;
    OUT.GOU.varY  = var_g;
    OUT.GOU.fano  = fano_g;
    OUT.GOU.Y60   = y_g;
    OUT.GOU.FPT   = FPT_g;
end

function OUT = simulate_moments_only(p, dt, T_total, rec_dt, n)
    n_steps = round(T_total/dt);
    rec_stride = round(rec_dt/dt);
    rec_idx = 0:rec_stride:n_steps;
    t_rec = rec_idx(:).*dt;

    x1_0 = binornd(p.N, p.p_x1_ss, n, 1);
    c1_0 = gamrnd(p.shape1, p.theta1, n, 1);
    c2_0 = gamrnd(p.shape2, p.theta2, n, 1);

    x1_c = x1_0;
    x2_c = p.N - x1_c;
    y_c  = zeros(n,1);
    c1_c = c1_0;
    c2_c = c2_0;

    x1_g = x1_0;
    x2_g = p.N - x1_g;
    y_g  = zeros(n,1);
    c1_g = c1_0;
    c2_g = c2_0;

    mean_c = zeros(numel(rec_idx),1);
    mean_g = zeros(numel(rec_idx),1);
    var_c  = zeros(numel(rec_idx),1);
    var_g  = zeros(numel(rec_idx),1);

    ir = 2;

    for it = 1:n_steps
        % CIR
        c1_old = c1_c;
        c2_old = c2_c;

        c1_c = cir_exact_transition( ...
            c1_c, p.kappa1, p.shape1, p.theta1, dt);
        c2_c = cir_exact_transition( ...
            c2_c, p.kappa2, p.shape2, p.theta2, dt);

        I1_c = 0.5 .* (c1_old + c1_c) .* dt;
        I2_c = 0.5 .* (c2_old + c2_c) .* dt;

        pr1_c = min(max(1 - exp(-I1_c), 0), 1);
        pr2_c = min(max(1 - exp(-I2_c), 0), 1);

        r1_c = binornd(x1_c, pr1_c);
        r2_c = binornd(x2_c, pr2_c);

        x1_c = x1_c - r1_c + r2_c;
        x2_c = p.N - x1_c;
        y_c  = y_c + r1_c;

        % GammaOU
        [c1_g, I1_g] = gamma_ou_step_with_integral( ...
            c1_g, p.a1, p.kappa1, p.theta1, dt);
        [c2_g, I2_g] = gamma_ou_step_with_integral( ...
            c2_g, p.a2, p.kappa2, p.theta2, dt);

        pr1_g = min(max(1 - exp(-I1_g), 0), 1);
        pr2_g = min(max(1 - exp(-I2_g), 0), 1);

        r1_g = binornd(x1_g, pr1_g);
        r2_g = binornd(x2_g, pr2_g);

        x1_g = x1_g - r1_g + r2_g;
        x2_g = p.N - x1_g;
        y_g  = y_g + r1_g;

        if ir <= numel(rec_idx) && it == rec_idx(ir)
            mean_c(ir) = mean(y_c);
            mean_g(ir) = mean(y_g);
            var_c(ir)  = var(y_c,0);
            var_g(ir)  = var(y_g,0);
            ir = ir + 1;
        end

        if mod(it, round(3600/dt)) == 0
            fprintf('  long-time t = %.1f min\n', it*dt/60);
        end
    end

    fano_c = var_c ./ max(mean_c,eps);
    fano_g = var_g ./ max(mean_g,eps);
    fano_c(1) = nan;
    fano_g(1) = nan;

    OUT.t = t_rec;
    OUT.CIR.meanY = mean_c;
    OUT.CIR.varY  = var_c;
    OUT.CIR.fano  = fano_c;
    OUT.GOU.meanY = mean_g;
    OUT.GOU.varY  = var_g;
    OUT.GOU.fano  = fano_g;
end

function F = finite_time_fano(t, V, F_int, A1, k1, A2, k2)
    t = t(:);
    F = nan(size(t));

    idx = t > 0;
    tt = t(idx);

    G1 = tt ./ k1 - (1 - exp(-k1 .* tt)) ./ (k1^2);
    G2 = tt ./ k2 - (1 - exp(-k2 .* tt)) ./ (k2^2);

    var_ext = 2 .* (A1 .* G1 + A2 .* G2);

    F(idx) = F_int + var_ext ./ (V .* tt);
    F(~idx) = 0;
end

function F = empirical_fpt_cdf(fpt, tgrid)
    F = arrayfun(@(t) mean(~isnan(fpt) & fpt <= t), tgrid);
    F = F(:);
end

function value = rmst_from_fpt(fpt, tau)
    z = fpt;
    z(isnan(z)) = tau;
    z = min(z, tau);
    value = mean(z);
end

function BOOT = bootstrap_high_order_effects( ...
        Y60_c, Y60_g, Mb, fpt_c, fpt_g, tau, n_boot)

    n_c = numel(Y60_c);
    n_g = numel(Y60_g);

    event_c = ~isnan(fpt_c) & fpt_c <= tau;
    event_g = ~isnan(fpt_g) & fpt_g <= tau;

    tc = fpt_c;
    tg = fpt_g;
    tc(isnan(tc)) = tau;
    tg(isnan(tg)) = tau;
    tc = min(tc,tau);
    tg = min(tg,tau);

    BOOT.dP_full = zeros(n_boot,1);
    BOOT.dF60    = zeros(n_boot,1);
    BOOT.dRMST   = zeros(n_boot,1);

    for b = 1:n_boot
        ic = randi(n_c, n_c, 1);
        ig = randi(n_g, n_g, 1);

        BOOT.dP_full(b) = ...
            mean(Y60_g(ig) >= Mb) - mean(Y60_c(ic) >= Mb);

        BOOT.dF60(b) = ...
            mean(event_g(ig)) - mean(event_c(ic));

        BOOT.dRMST(b) = mean(tc(ic)) - mean(tg(ig));
    end
end

function c_next = cir_exact_transition(c, kappa, shape, theta, dt)
    if dt == 0
        c_next = c;
        return;
    end

    e = exp(-kappa .* dt);
    one_minus_e = 1 - e;

    scale = 0.5 .* theta .* one_minus_e;
    df = 2 .* shape;
    nc = 2 .* e .* max(c,0) ./ (theta .* one_minus_e);

    if exist('ncx2rnd','file') == 2
        c_next = scale .* ncx2rnd(df, nc);
    else
        K = poissrnd(nc ./ 2);
        c_next = scale .* gamrnd(df/2 + K, 2);
    end

    c_next = max(c_next,0);
end

function [c_new, I] = gamma_ou_step_with_integral( ...
        c_old, a, kappa, theta, dt)

    if dt == 0
        c_new = c_old;
        I = zeros(size(c_old));
        return;
    end

    n = numel(c_old);

    decay = exp(-kappa .* dt);
    c_new = c_old .* decay;
    I = c_old .* (1 - decay) ./ kappa;

    Nj = poissrnd(a .* dt, n, 1);
    idx = find(Nj > 0);

    for ii = 1:numel(idx)
        j = idx(ii);
        m = Nj(j);

        U = dt .* rand(m,1);
        Y = gamrnd(1, theta, m,1);

        c_new(j) = c_new(j) + ...
            sum(Y .* exp(-kappa .* (dt-U)));

        I(j) = I(j) + ...
            sum(Y .* (1-exp(-kappa .* (dt-U))) ./ kappa);
    end

    c_new = max(c_new,0);
    I = max(I,0);
end

function export_figure_safe(fig, stem)
    if exist('exportgraphics','file') == 2
        exportgraphics(fig, [stem '.png'], 'Resolution', 300);
        exportgraphics(fig, [stem '.pdf'], 'ContentType','vector');
    else
        print(fig, [stem '.png'], '-dpng', '-r300');
        print(fig, [stem '.pdf'], '-dpdf', '-bestfit');
    end
end

function format_axes(ax)
    set(ax, 'Box','off','TickDir','out', ...
        'LineWidth',1.1,'FontSize',10.5,'FontName','Arial');
end

function out = blend_with_white(color, white_fraction)
    out = (1-white_fraction).*color + white_fraction.*[1 1 1];
end

function rgb = hex2rgb(hexStr)
    hexStr = strrep(hexStr,'#','');
    rgb = sscanf(hexStr,'%2x%2x%2x',[1 3])./255;
end
