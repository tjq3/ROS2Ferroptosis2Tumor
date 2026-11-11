% =========================================================================
% TDCME_FixedThreshold_dtConvergence_FINAL.m
%
% Purpose
% -------
% Numerical time-step convergence audit for the cytosolic full-state TD-CME
% solver used in the manuscript.
%
% IMPORTANT NUMERICAL INTERPRETATION
% ----------------------------------
% The continuous-time TD-CME is the model. The numerical implementation used
% here is a time-discretized integrated-hazard splitting solver:
%
%   * CIR environmental rates:
%       exact finite-step endpoint transition;
%       within-step rate integral approximated by the trapezoidal rule.
%
%   * GammaOU environmental rates:
%       exact exponential decay plus compound-Poisson jumps;
%       exact within-step integral conditional on sampled jump times/sizes.
%
%   * Chemical reactions:
%       integrated hazards are converted to one-step reaction probabilities,
%       followed by binomial updates of the finite Fe2+/Fe3+ pools.
%
% Therefore this script does NOT describe the numerical kernel as an exact
% SSA. Its purpose is to demonstrate convergence of the fixed scientific
% estimand as dt is refined.
%
% Key change relative to the previous dt-robustness script
% --------------------------------------------------------
% The ROS threshold is calibrated ONCE at the finest reference time step,
% and the SAME absolute threshold M_ref is then used for every dt.
%
% Thus the quantities
%
%   DeltaF_60    = P(T_GammaOU <= 60 min) - P(T_CIR <= 60 min)
%
%   DeltaRMST_60 = RMST_CIR(60 min) - RMST_GammaOU(60 min)
%
% are fixed estimands across the dt sweep.
%
% Positive values indicate earlier/higher GammaOU ROS-burden crossing.
%
% Baseline parameters are matched to the manuscript TD-CME setting:
%   N = 1200
%   c1: kappa=3.0e-4 s^-1, shape=0.20, mu=1.5e-5 s^-1
%   c2: kappa=4.0e-4 s^-1, shape=3.50, mu=1.7e-5 s^-1
%   Tobs = 3600 s = 60 min
%
% Recommended workflow:
%   run_level = 'smoke' -> 'pilot' -> 'final'
%
% Required:
%   Statistics and Machine Learning Toolbox
%
% =========================================================================

clear;
clc;
close all;

%% =========================================================================
% 0. Run configuration
% =========================================================================
run_level = 'final';       % 'smoke' | 'pilot' | 'final'

cfg = struct();

% Keep the same dt range already used in the study.
% The finest value is the reference discretization.
cfg.dt_list = [4.0 2.0 1.0 0.5];
cfg.dt_ref  = 0.5;
cfg.dt_main = 2.0;

cfg.T_obs = 3600;          % s
cfg.q_ref = 0.90;          % threshold calibrated once at dt_ref
cfg.n_boot = 1000;

% Independent random streams.
cfg.seed_ref_cal  = 20264001;
cfg.seed_eval_base = 20264100;
cfg.seed_boot_base = 20264200;

% Export.
cfg.export_figures = true;
cfg.export_tables  = true;
cfg.save_results   = true;
cfg.export_dpi     = 600;

switch lower(run_level)
    case 'smoke'
        cfg.n_cal_ref = 1200;
        cfg.n_eval    = 2500;
        cfg.n_boot    = 200;

    case 'pilot'
        cfg.n_cal_ref = 5000;
        cfg.n_eval    = 8000;
        cfg.n_boot    = 500;

    case 'final'
        % The reference calibration is independent of all evaluation runs.
        % Increase n_eval if a narrower CI is required.
        cfg.n_cal_ref = 12000;
        cfg.n_eval    = 20000;
        cfg.n_boot    = 1000;

    otherwise
        error('run_level must be smoke, pilot, or final.');
end

% Optional stronger audit:
% Add 0.25 s to cfg.dt_list and set cfg.dt_ref=0.25.
% This is substantially more expensive and is not required for the first
% publication-quality convergence pass if 1.0 and 0.5 s already agree.

validate_dt_configuration(cfg);

%% =========================================================================
% 1. Baseline model parameters
% =========================================================================
p = struct();

p.N = 1200;

% Fenton / oxidation channel c1.
p.kappa1 = 3.0e-4;
p.shape1 = 0.20;
p.mu_c1  = 1.5e-5;
p.theta1 = p.mu_c1 / p.shape1;
p.a1     = p.shape1 * p.kappa1;

% Recycling / reduction channel c2.
p.kappa2 = 4.0e-4;
p.shape2 = 3.50;
p.mu_c2  = 1.7e-5;
p.theta2 = p.mu_c2 / p.shape2;
p.a2     = p.shape2 * p.kappa2;

% Same initialization convention as the manuscript solver.
p.p_x1_ss = p.mu_c2 / (p.mu_c1 + p.mu_c2);

fprintf('============================================================\n');
fprintf('TD-CME fixed-threshold dt convergence audit\n');
fprintf('============================================================\n');
fprintf('Run level      : %s\n', run_level);
fprintf('dt list (s)    : %s\n', mat2str(cfg.dt_list));
fprintf('reference dt   : %.3f s\n', cfg.dt_ref);
fprintf('manuscript dt  : %.3f s\n', cfg.dt_main);
fprintf('T_obs          : %.1f min\n', cfg.T_obs/60);
fprintf('reference q    : %.3f\n', cfg.q_ref);
fprintf('n_cal_ref      : %d\n', cfg.n_cal_ref);
fprintf('n_eval per dt  : %d\n', cfg.n_eval);
fprintf('n_boot         : %d\n', cfg.n_boot);
fprintf('============================================================\n\n');

%% =========================================================================
% 2. Independent threshold calibration at the finest reference dt
% =========================================================================
rng(cfg.seed_ref_cal, 'twister');

fprintf('Calibrating one common absolute threshold at dt_ref=%.3f s...\n', ...
    cfg.dt_ref);

CALref = simulate_tdcme_splitting( ...
    p, cfg.dt_ref, cfg.T_obs, cfg.n_cal_ref, []);

pooled_ref = [CALref.CIR.Yend; CALref.GOU.Yend];

M_ref = max(1, round(quantile(pooled_ref, cfg.q_ref)));

fprintf('Reference threshold:\n');
fprintf('  M_ref = round(Q%.2f pooled Y60 at dt_ref) = %d\n\n', ...
    cfg.q_ref, M_ref);

%% =========================================================================
% 3. Evaluate every dt at the SAME absolute threshold
% =========================================================================
nDt = numel(cfg.dt_list);

S = repmat(empty_summary(), nDt, 1);
CI_F = nan(nDt,2);
CI_R = nan(nDt,2);

for idt = 1:nDt
    dt = cfg.dt_list(idt);

    fprintf('Evaluating dt %d/%d: %.3f s, fixed M_ref=%d\n', ...
        idt, nDt, dt, M_ref);

    rng(cfg.seed_eval_base + 100*idt, 'twister');

    OUT = simulate_tdcme_splitting( ...
        p, dt, cfg.T_obs, cfg.n_eval, M_ref);

    [S(idt), dF_i, dR_i_min] = summarize_fixed_estimand( ...
        OUT.CIR.Yend, OUT.GOU.Yend, ...
        OUT.CIR.FPT, OUT.GOU.FPT, ...
        M_ref, cfg.T_obs);

    rng(cfg.seed_boot_base + 100*idt, 'twister');

    [bootF, bootR] = paired_bootstrap_effect( ...
        dF_i, dR_i_min, cfg.n_boot);

    CI_F(idt,:) = prctile(bootF, [2.5 97.5]);
    CI_R(idt,:) = prctile(bootR, [2.5 97.5]);
end

%% =========================================================================
% 4. Assemble convergence table
% =========================================================================
dt_table = table( ...
    cfg.dt_list(:), ...
    repmat(M_ref,nDt,1), ...
    [S.P_CIR]', ...
    [S.P_GOU]', ...
    [S.DeltaF]', ...
    CI_F(:,1), CI_F(:,2), ...
    [S.RMST_CIR_min]', ...
    [S.RMST_GOU_min]', ...
    [S.DeltaRMST_min]', ...
    CI_R(:,1), CI_R(:,2), ...
    [S.meanY_CIR]', ...
    [S.meanY_GOU]', ...
    [S.fano_CIR]', ...
    [S.fano_GOU]', ...
    [S.endpoint_identity_error_CIR]', ...
    [S.endpoint_identity_error_GOU]', ...
    'VariableNames', { ...
        'dt_s','M_fixed', ...
        'F_CIR','F_GammaOU','DeltaF', ...
        'DeltaF_CI_low','DeltaF_CI_high', ...
        'RMST_CIR_min','RMST_GammaOU_min','DeltaRMST_min', ...
        'DeltaRMST_CI_low_min','DeltaRMST_CI_high_min', ...
        'meanY_CIR','meanY_GammaOU', ...
        'Fano_CIR','Fano_GammaOU', ...
        'endpoint_FPT_identity_error_CIR', ...
        'endpoint_FPT_identity_error_GammaOU'});

% Reference row.
[~, idxRef] = min(abs(dt_table.dt_s - cfg.dt_ref));
assert(abs(dt_table.dt_s(idxRef)-cfg.dt_ref)<1e-12);

% Absolute effect biases relative to the finest reference discretization.
dt_table.DeltaF_bias_vs_ref = ...
    dt_table.DeltaF - dt_table.DeltaF(idxRef);

dt_table.DeltaRMST_bias_vs_ref_min = ...
    dt_table.DeltaRMST_min - dt_table.DeltaRMST_min(idxRef);

% Model-specific low-order numerical deviations from reference.
dt_table.meanY_CIR_relerr_pct = ...
    100*(dt_table.meanY_CIR-dt_table.meanY_CIR(idxRef)) / ...
    max(abs(dt_table.meanY_CIR(idxRef)),eps);

dt_table.meanY_GammaOU_relerr_pct = ...
    100*(dt_table.meanY_GammaOU-dt_table.meanY_GammaOU(idxRef)) / ...
    max(abs(dt_table.meanY_GammaOU(idxRef)),eps);

dt_table.Fano_CIR_relerr_pct = ...
    100*(dt_table.Fano_CIR-dt_table.Fano_CIR(idxRef)) / ...
    max(abs(dt_table.Fano_CIR(idxRef)),eps);

dt_table.Fano_GammaOU_relerr_pct = ...
    100*(dt_table.Fano_GammaOU-dt_table.Fano_GammaOU(idxRef)) / ...
    max(abs(dt_table.Fano_GammaOU(idxRef)),eps);

fprintf('\nFixed-estimand convergence summary:\n');
disp(dt_table(:,{ ...
    'dt_s','M_fixed', ...
    'DeltaF','DeltaF_CI_low','DeltaF_CI_high', ...
    'DeltaF_bias_vs_ref', ...
    'DeltaRMST_min','DeltaRMST_CI_low_min', ...
    'DeltaRMST_CI_high_min','DeltaRMST_bias_vs_ref_min', ...
    'meanY_CIR_relerr_pct','meanY_GammaOU_relerr_pct', ...
    'Fano_CIR_relerr_pct','Fano_GammaOU_relerr_pct'}));

% Manuscript-dt summary.
[~,idxMain] = min(abs(dt_table.dt_s-cfg.dt_main));

fprintf('\nManuscript dt versus reference:\n');
fprintf('  dt_main=%.3f s, dt_ref=%.3f s\n', ...
    dt_table.dt_s(idxMain), dt_table.dt_s(idxRef));
fprintf('  DeltaF bias        = %+.6f\n', ...
    dt_table.DeltaF_bias_vs_ref(idxMain));
fprintf('  DeltaRMST bias     = %+.6f min\n', ...
    dt_table.DeltaRMST_bias_vs_ref_min(idxMain));
fprintf('  meanY CIR error    = %+.3f %%\n', ...
    dt_table.meanY_CIR_relerr_pct(idxMain));
fprintf('  meanY GammaOU err  = %+.3f %%\n', ...
    dt_table.meanY_GammaOU_relerr_pct(idxMain));
fprintf('  Fano CIR error     = %+.3f %%\n', ...
    dt_table.Fano_CIR_relerr_pct(idxMain));
fprintf('  Fano GammaOU error = %+.3f %%\n\n', ...
    dt_table.Fano_GammaOU_relerr_pct(idxMain));

%% =========================================================================
% 5. Plot publication-oriented SI figure
% =========================================================================
style = journal_style();

fig = figure( ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[1 1 18.3 12.2], ...
    'Renderer','painters');

tl = tiledlayout(fig,2,2, ...
    'TileSpacing','compact','Padding','compact');

% -------------------------------------------------------------------------
% a | DeltaF convergence
% -------------------------------------------------------------------------
ax = nexttile(tl,1);
hold(ax,'on');

plot_reference_band(ax, ...
    min(cfg.dt_list), max(cfg.dt_list), ...
    CI_F(idxRef,1), CI_F(idxRef,2), style.reference_fill);

errorbar(ax, dt_table.dt_s, dt_table.DeltaF, ...
    dt_table.DeltaF-dt_table.DeltaF_CI_low, ...
    dt_table.DeltaF_CI_high-dt_table.DeltaF, ...
    'o-', ...
    'Color',style.gou, ...
    'MarkerFaceColor','w', ...
    'MarkerSize',style.marker_size, ...
    'LineWidth',style.data_line, ...
    'CapSize',3);

yline(ax,dt_table.DeltaF(idxRef),'--', ...
    'Color',style.dark_gray,'LineWidth',0.8, ...
    'HandleVisibility','off');

xline(ax,cfg.dt_main,':', ...
    'Color',style.light_gray,'LineWidth',0.8, ...
    'HandleVisibility','off');

xlabel(ax,'Time step, dt (s)');
ylabel(ax,'\DeltaF_{60}');
title(ax,'Fixed-threshold crossing-risk convergence');

set(ax,'XDir','reverse','XTick',sort(cfg.dt_list));
format_axes(ax,style);
panel_label(ax,'a',style);

% -------------------------------------------------------------------------
% b | DeltaRMST convergence
% -------------------------------------------------------------------------
ax = nexttile(tl,2);
hold(ax,'on');

plot_reference_band(ax, ...
    min(cfg.dt_list), max(cfg.dt_list), ...
    CI_R(idxRef,1), CI_R(idxRef,2), style.reference_fill);

errorbar(ax, dt_table.dt_s, dt_table.DeltaRMST_min, ...
    dt_table.DeltaRMST_min-dt_table.DeltaRMST_CI_low_min, ...
    dt_table.DeltaRMST_CI_high_min-dt_table.DeltaRMST_min, ...
    'o-', ...
    'Color',style.cir, ...
    'MarkerFaceColor','w', ...
    'MarkerSize',style.marker_size, ...
    'LineWidth',style.data_line, ...
    'CapSize',3);

yline(ax,dt_table.DeltaRMST_min(idxRef),'--', ...
    'Color',style.dark_gray,'LineWidth',0.8, ...
    'HandleVisibility','off');

xline(ax,cfg.dt_main,':', ...
    'Color',style.light_gray,'LineWidth',0.8, ...
    'HandleVisibility','off');

xlabel(ax,'Time step, dt (s)');
ylabel(ax,'\DeltaRMST_{60} (min)');
title(ax,'Fixed-threshold FPT-timing convergence');

set(ax,'XDir','reverse','XTick',sort(cfg.dt_list));
format_axes(ax,style);
panel_label(ax,'b',style);

% -------------------------------------------------------------------------
% c | mean Y convergence
% -------------------------------------------------------------------------
ax = nexttile(tl,3);
hold(ax,'on');

plot(ax,dt_table.dt_s,dt_table.meanY_CIR,'o-', ...
    'Color',style.cir,'MarkerFaceColor','w', ...
    'MarkerSize',style.marker_size,'LineWidth',style.data_line);

plot(ax,dt_table.dt_s,dt_table.meanY_GammaOU,'s--', ...
    'Color',style.gou,'MarkerFaceColor','w', ...
    'MarkerSize',style.marker_size,'LineWidth',style.data_line);

xline(ax,cfg.dt_main,':', ...
    'Color',style.light_gray,'LineWidth',0.8, ...
    'HandleVisibility','off');

xlabel(ax,'Time step, dt (s)');
ylabel(ax,'Mean cumulative ROS, E[Y_{60}]');
title(ax,'Low-order mean convergence');
legend(ax,{'CIR','GammaOU'},'Box','off','Location','best');

set(ax,'XDir','reverse','XTick',sort(cfg.dt_list));
format_axes(ax,style);
panel_label(ax,'c',style);

% -------------------------------------------------------------------------
% d | Fano convergence
% -------------------------------------------------------------------------
ax = nexttile(tl,4);
hold(ax,'on');

plot(ax,dt_table.dt_s,dt_table.Fano_CIR,'o-', ...
    'Color',style.cir,'MarkerFaceColor','w', ...
    'MarkerSize',style.marker_size,'LineWidth',style.data_line);

plot(ax,dt_table.dt_s,dt_table.Fano_GammaOU,'s--', ...
    'Color',style.gou,'MarkerFaceColor','w', ...
    'MarkerSize',style.marker_size,'LineWidth',style.data_line);

xline(ax,cfg.dt_main,':', ...
    'Color',style.light_gray,'LineWidth',0.8, ...
    'HandleVisibility','off');

xlabel(ax,'Time step, dt (s)');
ylabel(ax,'Fano factor, Var(Y_{60})/E[Y_{60}]');
title(ax,'Low-order fluctuation convergence');
legend(ax,{'CIR','GammaOU'},'Box','off','Location','best');

set(ax,'XDir','reverse','XTick',sort(cfg.dt_list));
format_axes(ax,style);
panel_label(ax,'d',style);

sgtitle(tl, ...
    sprintf(['TD-CME numerical convergence at one fixed threshold ' ...
             '(M_{ref}=%d, calibrated at dt_{ref}=%.1f s)'], ...
             M_ref,cfg.dt_ref), ...
    'FontName','Arial','FontWeight','normal','FontSize',9.5);

%% =========================================================================
% 6. Export
% =========================================================================
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end

out_dir = fullfile(script_dir, ...
    ['results_TDCME_fixed_dt_convergence_' lower(run_level)]);

if ~exist(out_dir,'dir')
    mkdir(out_dir);
end

if cfg.export_tables
    writetable(dt_table,fullfile(out_dir, ...
        ['TDCME_fixed_threshold_dt_convergence_' lower(run_level) '.csv']));
end

if cfg.save_results
    save(fullfile(out_dir, ...
        ['TDCME_fixed_threshold_dt_convergence_' lower(run_level) '.mat']), ...
        'p','cfg','run_level','M_ref','dt_table','CALref','S','CI_F','CI_R', ...
        '-v7.3');
end

if cfg.export_figures
    stem = fullfile(out_dir, ...
        ['TDCME_fixed_threshold_dt_convergence_' lower(run_level)]);

    try
        exportgraphics(fig,[stem '.png'],'Resolution',cfg.export_dpi);
        exportgraphics(fig,[stem '.pdf'],'ContentType','vector');
    catch
        print(fig,[stem '.png'],'-dpng',sprintf('-r%d',cfg.export_dpi));
        print(fig,[stem '.pdf'],'-dpdf','-painters');
    end
end

fprintf('Completed.\nResults directory:\n%s\n',out_dir);

%% =========================================================================
% Local functions
% =========================================================================

function OUT = simulate_tdcme_splitting(p,dt,T,n,threshold)
% Time-discretized integrated-hazard splitting approximation to the
% continuous-time full-state TD-CME.
%
% threshold=[] performs endpoint simulation only.
%
% Shared initial conditions are used between CIR and GammaOU so that the
% within-dt scientific comparison is paired.

    validate_time_grid(dt,T);

    nSteps = round(T/dt);

    % Matched stationary initial conditions.
    x1Init = binornd(p.N,p.p_x1_ss,n,1);
    c1Init = gamrnd(p.shape1,p.theta1,n,1);
    c2Init = gamrnd(p.shape2,p.theta2,n,1);

    % CIR branch.
    x1C = x1Init;
    x2C = p.N-x1C;
    yC  = zeros(n,1);
    c1C = c1Init;
    c2C = c2Init;

    % GammaOU branch.
    x1G = x1Init;
    x2G = p.N-x1G;
    yG  = zeros(n,1);
    c1G = c1Init;
    c2G = c2Init;

    if isempty(threshold)
        fptC = [];
        fptG = [];
    else
        threshold = double(threshold);
        fptC = nan(n,1);
        fptG = nan(n,1);
    end

    for it = 1:nSteps
        tNow = it*dt;

        % ================================================================
        % CIR environment:
        % exact finite-step endpoint; trapezoidal within-step integral.
        % ================================================================
        c1Old = c1C;
        c2Old = c2C;

        c1C = cir_exact_transition( ...
            c1C,p.kappa1,p.shape1,p.theta1,dt);

        c2C = cir_exact_transition( ...
            c2C,p.kappa2,p.shape2,p.theta2,dt);

        I1C = 0.5*(c1Old+c1C)*dt;
        I2C = 0.5*(c2Old+c2C)*dt;

        p1C = min(max(1-exp(-I1C),0),1);
        p2C = min(max(1-exp(-I2C),0),1);

        % Splitting/binomial chemical update.
        r1C = binornd(x1C,p1C);
        r2C = binornd(x2C,p2C);

        x1C = x1C-r1C+r2C;
        x2C = p.N-x1C;
        yC  = yC+r1C;

        % ================================================================
        % GammaOU environment:
        % exact endpoint and exact within-step environmental integral.
        % ================================================================
        [c1G,I1G] = gamma_ou_step_with_integral( ...
            c1G,p.a1,p.kappa1,p.theta1,dt);

        [c2G,I2G] = gamma_ou_step_with_integral( ...
            c2G,p.a2,p.kappa2,p.theta2,dt);

        p1G = min(max(1-exp(-I1G),0),1);
        p2G = min(max(1-exp(-I2G),0),1);

        % Splitting/binomial chemical update.
        r1G = binornd(x1G,p1G);
        r2G = binornd(x2G,p2G);

        x1G = x1G-r1G+r2G;
        x2G = p.N-x1G;
        yG  = yG+r1G;

        % Grid-resolved first passage.
        if ~isempty(threshold)
            hitC = isnan(fptC) & yC>=threshold;
            hitG = isnan(fptG) & yG>=threshold;

            fptC(hitC) = tNow;
            fptG(hitG) = tNow;
        end
    end

    OUT.CIR.Yend = yC;
    OUT.GOU.Yend = yG;
    OUT.CIR.FPT = fptC;
    OUT.GOU.FPT = fptG;
end

function [S,dF_i,dR_i_min] = summarize_fixed_estimand( ...
        YC,YG,fptC,fptG,M,T)

    eventC = ~isnan(fptC) & fptC<=T;
    eventG = ~isnan(fptG) & fptG<=T;

    timeC = restricted_time(fptC,T);
    timeG = restricted_time(fptG,T);

    dF_i = double(eventG)-double(eventC);
    dR_i_min = (timeC-timeG)/60;

    S = empty_summary();

    S.P_CIR = mean(eventC);
    S.P_GOU = mean(eventG);
    S.DeltaF = mean(dF_i);

    S.RMST_CIR_min = mean(timeC)/60;
    S.RMST_GOU_min = mean(timeG)/60;
    S.DeltaRMST_min = mean(dR_i_min);

    S.meanY_CIR = mean(YC);
    S.meanY_GOU = mean(YG);

    S.fano_CIR = var(YC,0)/max(mean(YC),eps);
    S.fano_GOU = var(YG,0)/max(mean(YG),eps);

    % Because Y(t) is monotone, endpoint exceedance and terminal FPT event
    % should coincide up to numerical bookkeeping.
    S.endpoint_identity_error_CIR = ...
        abs(mean(YC>=M)-S.P_CIR);

    S.endpoint_identity_error_GOU = ...
        abs(mean(YG>=M)-S.P_GOU);
end

function S = empty_summary()
    S = struct( ...
        'P_CIR',nan, ...
        'P_GOU',nan, ...
        'DeltaF',nan, ...
        'RMST_CIR_min',nan, ...
        'RMST_GOU_min',nan, ...
        'DeltaRMST_min',nan, ...
        'meanY_CIR',nan, ...
        'meanY_GOU',nan, ...
        'fano_CIR',nan, ...
        'fano_GOU',nan, ...
        'endpoint_identity_error_CIR',nan, ...
        'endpoint_identity_error_GOU',nan);
end

function t = restricted_time(fpt,T)
    t = fpt;
    t(isnan(t)) = T;
    t = min(t,T);
end

function [bootF,bootR] = paired_bootstrap_effect(dF_i,dR_i,nBoot)
    dF_i = dF_i(:);
    dR_i = dR_i(:);

    n = numel(dF_i);

    bootF = zeros(nBoot,1);
    bootR = zeros(nBoot,1);

    for ib = 1:nBoot
        idx = randi(n,n,1);
        bootF(ib) = mean(dF_i(idx));
        bootR(ib) = mean(dR_i(idx));
    end
end

function cNext = cir_exact_transition(c,kappa,shape,theta,dt)
% Exact finite-time endpoint transition for
%   dc = kappa*(mu-c)dt + sqrt(2*kappa*theta*c)dW,
% with mu=shape*theta.

    if dt==0
        cNext = c;
        return;
    end

    decay = exp(-kappa*dt);
    oneMinusDecay = 1-decay;

    scale = 0.5*theta*oneMinusDecay;
    df = 2*shape;

    nc = 2*decay.*max(c,0) ./ ...
        (theta*oneMinusDecay);

    if exist('ncx2rnd','file')==2
        cNext = scale.*ncx2rnd(df,nc);
    else
        poissonMix = poissrnd(nc/2);
        cNext = scale.*gamrnd(df/2+poissonMix,2);
    end

    cNext = max(cNext,0);
end

function [cNew,I] = gamma_ou_step_with_integral( ...
        cOld,a,kappa,theta,dt)
% Exact finite-step endpoint and exact within-step integral for the
% compound-Poisson GammaOU environment conditional on sampled jumps:
%
%   dc = -kappa*c dt + dL(t)
%
% jump arrival rate = a
% jump size ~ Exponential(scale=theta).

    if dt==0
        cNew = cOld;
        I = zeros(size(cOld));
        return;
    end

    n = numel(cOld);
    decay = exp(-kappa*dt);

    % Deterministic decay of the starting value.
    cNew = cOld.*decay;
    I = cOld.*(1-decay)/kappa;

    % Compound-Poisson jumps inside the step.
    nJump = poissrnd(a*dt,n,1);
    totalJump = sum(nJump);

    if totalJump>0
        owner = repelem((1:n)',nJump);

        jumpTime = dt*rand(totalJump,1);
        jumpSize = gamrnd(1,theta,totalJump,1);

        endpointContribution = ...
            jumpSize.*exp(-kappa*(dt-jumpTime));

        integralContribution = ...
            jumpSize.*(1-exp(-kappa*(dt-jumpTime)))/kappa;

        cNew = cNew + accumarray( ...
            owner,endpointContribution,[n 1],@sum,0);

        I = I + accumarray( ...
            owner,integralContribution,[n 1],@sum,0);
    end

    cNew = max(cNew,0);
    I = max(I,0);
end

function validate_dt_configuration(cfg)
    if ~any(abs(cfg.dt_list-cfg.dt_ref)<1e-12)
        error('cfg.dt_ref must be included in cfg.dt_list.');
    end

    if ~any(abs(cfg.dt_list-cfg.dt_main)<1e-12)
        error('cfg.dt_main must be included in cfg.dt_list.');
    end

    for dt = cfg.dt_list
        validate_time_grid(dt,cfg.T_obs);
    end
end

function validate_time_grid(dt,T)
    n = T/dt;
    if abs(n-round(n))>1e-10
        error('T_obs must be an integer multiple of every dt.');
    end
end

function plot_reference_band(ax,xmin,xmax,ylo,yhi,faceColor)
    if ~isfinite(ylo) || ~isfinite(yhi)
        return;
    end

    patch(ax, ...
        [xmin xmax xmax xmin], ...
        [ylo ylo yhi yhi], ...
        faceColor, ...
        'EdgeColor','none', ...
        'FaceAlpha',0.35, ...
        'HandleVisibility','off');
end

function style = journal_style()
    style.font_name = 'Arial';
    style.axis_font = 7.3;
    style.label_font = 7.8;
    style.title_font = 8.2;
    style.panel_font = 10.5;
    style.marker_size = 4.8;
    style.data_line = 1.25;

    style.cir = [0.235 0.329 0.533];
    style.gou = [0.902 0.294 0.208];

    style.dark_gray = [0.25 0.25 0.25];
    style.light_gray = [0.72 0.72 0.72];
    style.reference_fill = [0.88 0.88 0.88];
end

function format_axes(ax,style)
    set(ax, ...
        'FontName',style.font_name, ...
        'FontSize',style.axis_font, ...
        'LineWidth',0.75, ...
        'TickDir','out', ...
        'Box','off', ...
        'Layer','top');
end

function panel_label(ax,label,style)
    text(ax,-0.12,1.06,label, ...
        'Units','normalized', ...
        'FontName',style.font_name, ...
        'FontSize',style.panel_font, ...
        'FontWeight','bold', ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','bottom');
end
