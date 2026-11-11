% =========================================================================
% TDCME_CoreRobustness.m
%
% Core robustness analyses for the cytosolic strict TD-CME study
% -------------------------------------------------------------------------
% Scientific baseline:
%   CIR and GammaOU rate environments are matched in stationary Gamma
%   marginal and exponential autocorrelation, but differ in higher-order
%   temporal path architecture.
%
% This script implements the four robustness analyses that are most directly
% required for the main TD-CME claim:
%
%   1. Threshold-quantile robustness
%      q = 0.85, 0.875, 0.90, 0.925, 0.95 at the fixed 60-min horizon.
%
%   2. Observation-horizon robustness
%      T_obs = 30, 45, 60, 75, 90 min, with an independently calibrated
%      pooled Q0.90 threshold at EACH horizon and FPT/RMST restricted to the
%      SAME horizon.
%
%   3. Numerical time-step convergence
%      dt = 4, 2, 1, 0.5 s at T_obs=60 min and q=0.90.
%
%   4. Independent-seed replication
%      Repeated independent calibration/evaluation batches at the manuscript
%      baseline dt=2 s, T_obs=60 min and q=0.90.
%
% Primary estimands:
%
%   DeltaF_T =
%       P(T_GammaOU <= T) - P(T_CIR <= T)
%
%   DeltaRMST_T =
%       RMST_CIR(T) - RMST_GammaOU(T)
%
% Positive values indicate earlier/higher GammaOU ROS-burden crossing.
%
% Important terminology:
%   These are ROS-burden first-passage times, not ferroptosis FPTs.
%
% Strict kernel retained:
%   N=1200
%   c1: kappa1=3.0e-4 s^-1, shape1=0.20, mu_c1=1.5e-5 s^-1
%   c2: kappa2=4.0e-4 s^-1, shape2=3.50, mu_c2=1.7e-5 s^-1
%   matched stationary-binomial iron initialization
%   matched stationary-Gamma rate initialization
%   exact CIR transition
%   exact GammaOU decay plus within-step jump integral
%   integrated-hazard binomial reaction updates
%   independent calibration and evaluation
%   paired bootstrap using the same trajectory indices for CIR/GammaOU
%
% Recommended workflow:
%   run_level = 'smoke' -> 'pilot' -> 'final'
%
% The script is serial by default to avoid Windows non-ASCII path problems.
% =========================================================================

clear;
clc;
close all;

%% =========================================================================
% 0. Run configuration
% =========================================================================
run_level = 'final';      % 'smoke' | 'pilot' | 'final'

run_threshold_horizon = true;
run_dt_convergence     = true;
run_seed_replication   = true;

export_figures = true;
export_tables  = true;
save_results   = true;
save_raw       = false;

cfg = struct();
cfg.export_dpi = 600;

% Analysis-choice grids.
cfg.q_list = [0.85 0.875 0.90 0.925 0.95];
cfg.q_main = 0.90;

cfg.horizon_min = [30 45 60 75 90];
cfg.horizon_main_min = 60;

% Numerical convergence grid.
cfg.dt_list = [4.0 2.0 1.0 0.5];
cfg.dt_main = 2.0;

% Independent-seed replication.
cfg.n_seed_batches = 5;

% Fixed seeds. Each module uses independent random-number streams.
cfg.seed_choice_cal  = 20263001;
cfg.seed_choice_eval = 20263002;
cfg.seed_choice_boot = 20263003;

cfg.seed_dt_base     = 20263100;
cfg.seed_seed_base   = 20263200;

switch lower(run_level)
    case 'smoke'
        cfg.choice_n_cal  = 1200;
        cfg.choice_n_eval = 2500;
        cfg.choice_n_boot = 200;

        cfg.dt_n_cal  = 700;
        cfg.dt_n_eval = 1400;
        cfg.dt_n_boot = 150;

        cfg.seed_n_cal  = 700;
        cfg.seed_n_eval = 1400;
        cfg.seed_n_boot = 150;

    case 'pilot'
        cfg.choice_n_cal  = 5000;
        cfg.choice_n_eval = 10000;
        cfg.choice_n_boot = 500;

        cfg.dt_n_cal  = 2500;
        cfg.dt_n_eval = 5000;
        cfg.dt_n_boot = 400;

        cfg.seed_n_cal  = 2500;
        cfg.seed_n_eval = 5000;
        cfg.seed_n_boot = 400;

    case 'final'
        % Analysis-choice robustness uses the same evaluation scale as the
        % main full-horizon result.
        cfg.choice_n_cal  = 12000;
        cfg.choice_n_eval = 30000;
        cfg.choice_n_boot = 1000;

        % dt=0.5 s is expensive. These values are large enough for an SI
        % convergence audit without multiplying the main cost excessively.
        cfg.dt_n_cal  = 6000;
        cfg.dt_n_eval = 12000;
        cfg.dt_n_boot = 800;

        % Five independent baseline batches.
        cfg.seed_n_cal  = 5000;
        cfg.seed_n_eval = 10000;
        cfg.seed_n_boot = 600;

    otherwise
        error('run_level must be smoke, pilot, or final.');
end

% Output directory.
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end

out_dir = fullfile(script_dir, ...
    ['results_TDCME_core_robustness_' lower(run_level)]);

if ~exist(out_dir,'dir')
    mkdir(out_dir);
end

log_file = fullfile(out_dir, ...
    ['TDCME_core_robustness_' lower(run_level) '.log']);

if exist(log_file,'file')
    delete(log_file);
end

diary(log_file);
cleanup_diary = onCleanup(@() diary('off')); %#ok<NASGU>

%% =========================================================================
% 1. Baseline parameters
% =========================================================================
p = struct();

p.N = 1200;

% Oxidation/Fenton channel c1.
p.kappa1 = 3.0e-4;
p.shape1 = 0.20;
p.mu_c1  = 1.5e-5;
p.theta1 = p.mu_c1/p.shape1;
p.a1     = p.shape1*p.kappa1;

% Reduction/recycling channel c2.
p.kappa2 = 4.0e-4;
p.shape2 = 3.50;
p.mu_c2  = 1.7e-5;
p.theta2 = p.mu_c2/p.shape2;
p.a2     = p.shape2*p.kappa2;

% Shared stationary iron initialization.
p.p_x1_ss = p.mu_c2/(p.mu_c1+p.mu_c2);

style = journal_style();

fprintf('============================================================\n');
fprintf('TD-CME core SI robustness analyses\n');
fprintf('============================================================\n');
fprintf('Run level: %s\n',run_level);
fprintf('N=%d\n',p.N);
fprintf('c1: kappa=%.6e, shape=%.3f, mean=%.6e\n', ...
    p.kappa1,p.shape1,p.mu_c1);
fprintf('c2: kappa=%.6e, shape=%.3f, mean=%.6e\n', ...
    p.kappa2,p.shape2,p.mu_c2);
fprintf('q values: %s\n',mat2str(cfg.q_list));
fprintf('horizons (min): %s\n',mat2str(cfg.horizon_min));
fprintf('dt values (s): %s\n',mat2str(cfg.dt_list));
fprintf('seed batches: %d\n',cfg.n_seed_batches);
fprintf('============================================================\n\n');

RESULT = struct();
RESULT.parameters = p;
RESULT.configuration = cfg;
RESULT.run_level = run_level;

%% =========================================================================
% 2. Threshold-quantile and horizon robustness
% =========================================================================
if run_threshold_horizon
    fprintf('------------------------------------------------------------\n');
    fprintf('Module 1: threshold-quantile and horizon robustness\n');
    fprintf('------------------------------------------------------------\n');

    dt = cfg.dt_main;
    horizon_s = cfg.horizon_min(:)'*60;
    T_max = max(horizon_s);

    validate_time_points(dt,horizon_s);

    % ---------------- Independent calibration ---------------------------
    rng(cfg.seed_choice_cal,'twister');

    fprintf('Independent multi-horizon calibration starts...\n');

    CAL = simulate_multi_horizon_tdcme( ...
        p,dt,T_max,horizon_s,cfg.choice_n_cal,[]);

    fprintf('Independent multi-horizon calibration finished.\n');

    [~,idx_main_horizon] = min(abs( ...
        cfg.horizon_min-cfg.horizon_main_min));

    pooled_cal_main = [ ...
        CAL.CIR.Y(:,idx_main_horizon); ...
        CAL.GOU.Y(:,idx_main_horizon)];

    Mb_q = zeros(size(cfg.q_list));

    for iq=1:numel(cfg.q_list)
        Mb_q(iq) = max(1,round(quantile( ...
            pooled_cal_main,cfg.q_list(iq))));
    end

    Mb_horizon = zeros(size(horizon_s));

    for ih=1:numel(horizon_s)
        pooled_h = [CAL.CIR.Y(:,ih);CAL.GOU.Y(:,ih)];

        Mb_horizon(ih) = max(1,round(quantile( ...
            pooled_h,cfg.q_main)));
    end

    all_thresholds = unique([Mb_q(:);Mb_horizon(:)])';

    [~,map_q] = ismember(Mb_q,all_thresholds);
    [~,map_h] = ismember(Mb_horizon,all_thresholds);

    fprintf('Calibrated q-sweep thresholds at 60 min:\n');
    disp(table(cfg.q_list(:),Mb_q(:), ...
        'VariableNames',{'q','Mb'}));

    fprintf('Calibrated Q0.90 thresholds across horizons:\n');
    disp(table(cfg.horizon_min(:),Mb_horizon(:), ...
        'VariableNames',{'horizon_min','Mb'}));

    % ---------------- Independent paired evaluation ---------------------
    rng(cfg.seed_choice_eval,'twister');

    fprintf('Independent multi-threshold evaluation starts...\n');

    EVAL = simulate_multi_horizon_tdcme( ...
        p,dt,T_max,horizon_s,cfg.choice_n_eval,all_thresholds);

    fprintf('Independent multi-threshold evaluation finished.\n');

    % Contribution matrices. Columns are all q and horizon conditions.
    nQ = numel(cfg.q_list);
    nH = numel(horizon_s);
    nCond = nQ+nH;

    dF_contrib = zeros(cfg.choice_n_eval,nCond);
    dR_contrib_min = zeros(cfg.choice_n_eval,nCond);

    threshold_summary = repmat(empty_effect_summary(),nQ,1);
    horizon_summary   = repmat(empty_effect_summary(),nH,1);

    % q sweep at the fixed 60-min horizon.
    tau_main = cfg.horizon_main_min*60;

    for iq=1:nQ
        threshold_index = map_q(iq);

        YC = EVAL.CIR.Y(:,idx_main_horizon);
        YG = EVAL.GOU.Y(:,idx_main_horizon);

        fptC = EVAL.CIR.FPT(:,threshold_index);
        fptG = EVAL.GOU.FPT(:,threshold_index);

        [S,dF_i,dR_i] = summarize_effect( ...
            YC,YG,fptC,fptG,Mb_q(iq),tau_main);

        threshold_summary(iq) = S;
        dF_contrib(:,iq) = dF_i;
        dR_contrib_min(:,iq) = dR_i/60;
    end

    % Horizon sweep. Threshold and FPT horizon are time matched.
    for ih=1:nH
        threshold_index = map_h(ih);
        tau = horizon_s(ih);

        YC = EVAL.CIR.Y(:,ih);
        YG = EVAL.GOU.Y(:,ih);

        fptC = EVAL.CIR.FPT(:,threshold_index);
        fptG = EVAL.GOU.FPT(:,threshold_index);

        [S,dF_i,dR_i] = summarize_effect( ...
            YC,YG,fptC,fptG,Mb_horizon(ih),tau);

        horizon_summary(ih) = S;
        dF_contrib(:,nQ+ih) = dF_i;
        dR_contrib_min(:,nQ+ih) = dR_i/60;
    end

    % One paired bootstrap serves all q and horizon conditions.
    rng(cfg.seed_choice_boot,'twister');

    [bootF,bootR] = paired_bootstrap_matrix( ...
        dF_contrib,dR_contrib_min,cfg.choice_n_boot);

    CI_F = prctile(bootF,[2.5 97.5],1);
    CI_R = prctile(bootR,[2.5 97.5],1);

    % Assemble q table.
    q_table = table( ...
        cfg.q_list(:),Mb_q(:), ...
        [threshold_summary.P_CIR]', ...
        [threshold_summary.P_GOU]', ...
        [threshold_summary.DeltaF]', ...
        CI_F(1,1:nQ)',CI_F(2,1:nQ)', ...
        [threshold_summary.RMST_CIR_s]'/60, ...
        [threshold_summary.RMST_GOU_s]'/60, ...
        [threshold_summary.DeltaRMST_s]'/60, ...
        CI_R(1,1:nQ)',CI_R(2,1:nQ)', ...
        [threshold_summary.meanY_CIR]', ...
        [threshold_summary.meanY_GOU]', ...
        100*[threshold_summary.signedMeanDifference]', ...
        [threshold_summary.fano_CIR]', ...
        [threshold_summary.fano_GOU]', ...
        100*[threshold_summary.signedFanoDifference]', ...
        [threshold_summary.identityError_CIR]', ...
        [threshold_summary.identityError_GOU]', ...
        'VariableNames',{ ...
        'q','Mb', ...
        'F_CIR','F_GammaOU','DeltaF', ...
        'DeltaF_CI_low','DeltaF_CI_high', ...
        'RMST_CIR_min','RMST_GammaOU_min','DeltaRMST_min', ...
        'DeltaRMST_CI_low_min','DeltaRMST_CI_high_min', ...
        'meanY_CIR','meanY_GammaOU','signedMeanDifference_pct', ...
        'Fano_CIR','Fano_GammaOU','signedFanoDifference_pct', ...
        'endpoint_FPT_identity_error_CIR', ...
        'endpoint_FPT_identity_error_GammaOU'});

    % Assemble horizon table.
    offset = nQ;
    horizon_table = table( ...
        cfg.horizon_min(:),Mb_horizon(:), ...
        [horizon_summary.P_CIR]', ...
        [horizon_summary.P_GOU]', ...
        [horizon_summary.DeltaF]', ...
        CI_F(1,offset+(1:nH))',CI_F(2,offset+(1:nH))', ...
        [horizon_summary.RMST_CIR_s]'/60, ...
        [horizon_summary.RMST_GOU_s]'/60, ...
        [horizon_summary.DeltaRMST_s]'/60, ...
        CI_R(1,offset+(1:nH))',CI_R(2,offset+(1:nH))', ...
        100*([horizon_summary.DeltaRMST_s]'./horizon_s(:)), ...
        [horizon_summary.meanY_CIR]', ...
        [horizon_summary.meanY_GOU]', ...
        100*[horizon_summary.signedMeanDifference]', ...
        [horizon_summary.fano_CIR]', ...
        [horizon_summary.fano_GOU]', ...
        100*[horizon_summary.signedFanoDifference]', ...
        [horizon_summary.identityError_CIR]', ...
        [horizon_summary.identityError_GOU]', ...
        'VariableNames',{ ...
        'horizon_min','Mb', ...
        'F_CIR','F_GammaOU','DeltaF', ...
        'DeltaF_CI_low','DeltaF_CI_high', ...
        'RMST_CIR_min','RMST_GammaOU_min','DeltaRMST_min', ...
        'DeltaRMST_CI_low_min','DeltaRMST_CI_high_min', ...
        'DeltaRMST_over_horizon_pct', ...
        'meanY_CIR','meanY_GammaOU','signedMeanDifference_pct', ...
        'Fano_CIR','Fano_GammaOU','signedFanoDifference_pct', ...
        'endpoint_FPT_identity_error_CIR', ...
        'endpoint_FPT_identity_error_GammaOU'});

    RESULT.threshold_quantile.table = q_table;
    RESULT.horizon.table = horizon_table;

    if save_raw
        RESULT.threshold_horizon_raw.CAL = CAL;
        RESULT.threshold_horizon_raw.EVAL = EVAL;
        RESULT.threshold_horizon_raw.bootF = bootF;
        RESULT.threshold_horizon_raw.bootR = bootR;
    end

    fprintf('\nThreshold-quantile robustness summary:\n');
    disp(q_table(:,{ ...
        'q','Mb','DeltaF','DeltaF_CI_low','DeltaF_CI_high', ...
        'DeltaRMST_min','DeltaRMST_CI_low_min', ...
        'DeltaRMST_CI_high_min'}));

    fprintf('\nObservation-horizon robustness summary:\n');
    disp(horizon_table(:,{ ...
        'horizon_min','Mb','DeltaF','DeltaF_CI_low', ...
        'DeltaF_CI_high','DeltaRMST_min', ...
        'DeltaRMST_CI_low_min','DeltaRMST_CI_high_min'}));

    % ---------------- Figure 1 ------------------------------------------
    fig1 = figure( ...
        'Color','w', ...
        'Units','centimeters', ...
        'Position',[1 1 18.3 11.8], ...
        'Renderer','painters');

    tl1 = tiledlayout(fig1,2,3, ...
        'TileSpacing','compact','Padding','compact');

    % a | Mb versus q.
    ax = nexttile(tl1,1);
    plot(ax,q_table.q,q_table.Mb,'o-', ...
        'Color',style.dark_gray, ...
        'MarkerFaceColor','w', ...
        'MarkerSize',style.marker_size, ...
        'LineWidth',style.data_line);
    xline(ax,cfg.q_main,':', ...
        'Color',style.light_gray, ...
        'HandleVisibility','off');
    xlabel(ax,'Threshold quantile, q');
    ylabel(ax,'Calibrated M_b');
    title(ax,'Threshold calibration');
    panel_label(ax,'a',style);
    format_axes(ax,style);

    % b | DeltaF versus q.
    ax = nexttile(tl1,2);
    errorbar(ax,q_table.q,q_table.DeltaF, ...
        q_table.DeltaF-q_table.DeltaF_CI_low, ...
        q_table.DeltaF_CI_high-q_table.DeltaF, ...
        'o-','Color',style.gou, ...
        'MarkerFaceColor','w', ...
        'MarkerSize',style.marker_size, ...
        'LineWidth',style.data_line, ...
        'CapSize',3);
    yline(ax,0,':','Color',style.light_gray, ...
        'HandleVisibility','off');
    xline(ax,cfg.q_main,':','Color',style.light_gray, ...
        'HandleVisibility','off');
    xlabel(ax,'Threshold quantile, q');
    ylabel(ax,'\DeltaF_{60}');
    title(ax,'Crossing-probability robustness');
    panel_label(ax,'b',style);
    format_axes(ax,style);

    % c | DeltaRMST versus q.
    ax = nexttile(tl1,3);
    errorbar(ax,q_table.q,q_table.DeltaRMST_min, ...
        q_table.DeltaRMST_min-q_table.DeltaRMST_CI_low_min, ...
        q_table.DeltaRMST_CI_high_min-q_table.DeltaRMST_min, ...
        'o-','Color',style.cir, ...
        'MarkerFaceColor','w', ...
        'MarkerSize',style.marker_size, ...
        'LineWidth',style.data_line, ...
        'CapSize',3);
    yline(ax,0,':','Color',style.light_gray, ...
        'HandleVisibility','off');
    xline(ax,cfg.q_main,':','Color',style.light_gray, ...
        'HandleVisibility','off');
    xlabel(ax,'Threshold quantile, q');
    ylabel(ax,'\DeltaRMST_{60} (min)');
    title(ax,'First-passage timing robustness');
    panel_label(ax,'c',style);
    format_axes(ax,style);

    % d | Mb versus horizon.
    ax = nexttile(tl1,4);
    plot(ax,horizon_table.horizon_min,horizon_table.Mb,'s-', ...
        'Color',style.dark_gray, ...
        'MarkerFaceColor','w', ...
        'MarkerSize',style.marker_size, ...
        'LineWidth',style.data_line);
    xline(ax,cfg.horizon_main_min,':', ...
        'Color',style.light_gray, ...
        'HandleVisibility','off');
    xlabel(ax,'Observation horizon, T_{obs} (min)');
    ylabel(ax,'Calibrated M_b(T_{obs})');
    title(ax,'Time-matched threshold calibration');
    panel_label(ax,'d',style);
    format_axes(ax,style);

    % e | DeltaF versus horizon.
    ax = nexttile(tl1,5);
    errorbar(ax,horizon_table.horizon_min,horizon_table.DeltaF, ...
        horizon_table.DeltaF-horizon_table.DeltaF_CI_low, ...
        horizon_table.DeltaF_CI_high-horizon_table.DeltaF, ...
        's-','Color',style.gou, ...
        'MarkerFaceColor','w', ...
        'MarkerSize',style.marker_size, ...
        'LineWidth',style.data_line, ...
        'CapSize',3);
    yline(ax,0,':','Color',style.light_gray, ...
        'HandleVisibility','off');
    xline(ax,cfg.horizon_main_min,':', ...
        'Color',style.light_gray, ...
        'HandleVisibility','off');
    xlabel(ax,'Observation horizon, T_{obs} (min)');
    ylabel(ax,'\DeltaF_{T_{obs}}');
    title(ax,'Horizon robustness of crossing risk');
    panel_label(ax,'e',style);
    format_axes(ax,style);

    % f | DeltaRMST versus horizon.
    ax = nexttile(tl1,6);
    errorbar(ax,horizon_table.horizon_min, ...
        horizon_table.DeltaRMST_min, ...
        horizon_table.DeltaRMST_min- ...
            horizon_table.DeltaRMST_CI_low_min, ...
        horizon_table.DeltaRMST_CI_high_min- ...
            horizon_table.DeltaRMST_min, ...
        's-','Color',style.cir, ...
        'MarkerFaceColor','w', ...
        'MarkerSize',style.marker_size, ...
        'LineWidth',style.data_line, ...
        'CapSize',3);
    yline(ax,0,':','Color',style.light_gray, ...
        'HandleVisibility','off');
    xline(ax,cfg.horizon_main_min,':', ...
        'Color',style.light_gray, ...
        'HandleVisibility','off');
    xlabel(ax,'Observation horizon, T_{obs} (min)');
    ylabel(ax,'\DeltaRMST_{T_{obs}} (min)');
    title(ax,'Horizon robustness of timing shift');
    panel_label(ax,'f',style);
    format_axes(ax,style);

    if export_figures
        export_figure(fig1,fullfile(out_dir, ...
            ['TDCME_SI_threshold_horizon_robustness_' ...
            lower(run_level)]),cfg.export_dpi);
    end

    if export_tables
        writetable(q_table,fullfile(out_dir, ...
            ['TDCME_threshold_quantile_robustness_' ...
            lower(run_level) '.csv']));

        writetable(horizon_table,fullfile(out_dir, ...
            ['TDCME_horizon_robustness_' ...
            lower(run_level) '.csv']));
    end
end

%% =========================================================================
% 3. Numerical time-step convergence
% =========================================================================
if run_dt_convergence
    fprintf('\n------------------------------------------------------------\n');
    fprintf('Module 2: numerical time-step convergence\n');
    fprintf('------------------------------------------------------------\n');

    nDt = numel(cfg.dt_list);
    dt_summary = repmat(empty_effect_summary(),nDt,1);
    dt_Mb = zeros(nDt,1);
    dt_CI_F = zeros(nDt,2);
    dt_CI_R = zeros(nDt,2);

    T_s = cfg.horizon_main_min*60;

    for idt=1:nDt
        dt = cfg.dt_list(idt);
        validate_time_points(dt,T_s);

        seed0 = cfg.seed_dt_base+100*idt;

        fprintf('dt convergence %d/%d: dt=%.3f s\n', ...
            idt,nDt,dt);

        rng(seed0+1,'twister');

        CALdt = simulate_multi_horizon_tdcme( ...
            p,dt,T_s,T_s,cfg.dt_n_cal,[]);

        pooled_cal = [CALdt.CIR.Y(:,1);CALdt.GOU.Y(:,1)];

        dt_Mb(idt) = max(1,round(quantile( ...
            pooled_cal,cfg.q_main)));

        rng(seed0+2,'twister');

        EVALdt = simulate_multi_horizon_tdcme( ...
            p,dt,T_s,T_s,cfg.dt_n_eval,dt_Mb(idt));

        [S,dF_i,dR_i] = summarize_effect( ...
            EVALdt.CIR.Y(:,1), ...
            EVALdt.GOU.Y(:,1), ...
            EVALdt.CIR.FPT(:,1), ...
            EVALdt.GOU.FPT(:,1), ...
            dt_Mb(idt),T_s);

        dt_summary(idt) = S;

        rng(seed0+3,'twister');

        [bootFdt,bootRdt] = paired_bootstrap_matrix( ...
            dF_i,dR_i/60,cfg.dt_n_boot);

        dt_CI_F(idt,:) = prctile(bootFdt,[2.5 97.5]);
        dt_CI_R(idt,:) = prctile(bootRdt,[2.5 97.5]);
    end

    dt_table = table( ...
        cfg.dt_list(:),dt_Mb(:), ...
        [dt_summary.P_CIR]', ...
        [dt_summary.P_GOU]', ...
        [dt_summary.DeltaF]', ...
        dt_CI_F(:,1),dt_CI_F(:,2), ...
        [dt_summary.RMST_CIR_s]'/60, ...
        [dt_summary.RMST_GOU_s]'/60, ...
        [dt_summary.DeltaRMST_s]'/60, ...
        dt_CI_R(:,1),dt_CI_R(:,2), ...
        [dt_summary.meanY_CIR]', ...
        [dt_summary.meanY_GOU]', ...
        100*[dt_summary.signedMeanDifference]', ...
        [dt_summary.fano_CIR]', ...
        [dt_summary.fano_GOU]', ...
        100*[dt_summary.signedFanoDifference]', ...
        'VariableNames',{ ...
        'dt_s','Mb', ...
        'F_CIR','F_GammaOU','DeltaF', ...
        'DeltaF_CI_low','DeltaF_CI_high', ...
        'RMST_CIR_min','RMST_GammaOU_min','DeltaRMST_min', ...
        'DeltaRMST_CI_low_min','DeltaRMST_CI_high_min', ...
        'meanY_CIR','meanY_GammaOU','signedMeanDifference_pct', ...
        'Fano_CIR','Fano_GammaOU','signedFanoDifference_pct'});

    % Relative difference from the finest time step.
    [~,idx_finest] = min(dt_table.dt_s);

    dt_table.DeltaF_difference_from_finest = ...
        dt_table.DeltaF-dt_table.DeltaF(idx_finest);

    dt_table.DeltaRMST_difference_from_finest_min = ...
        dt_table.DeltaRMST_min- ...
        dt_table.DeltaRMST_min(idx_finest);

    RESULT.dt_convergence.table = dt_table;

    fprintf('\nNumerical convergence summary:\n');
    disp(dt_table(:,{ ...
        'dt_s','Mb','DeltaF','DeltaF_CI_low','DeltaF_CI_high', ...
        'DeltaRMST_min','DeltaRMST_CI_low_min', ...
        'DeltaRMST_CI_high_min', ...
        'DeltaRMST_difference_from_finest_min'}));

    if export_tables
        writetable(dt_table,fullfile(out_dir, ...
            ['TDCME_dt_convergence_' lower(run_level) '.csv']));
    end
end

%% =========================================================================
% 4. Independent-seed replication
% =========================================================================
if run_seed_replication
    fprintf('\n------------------------------------------------------------\n');
    fprintf('Module 3: independent-seed replication\n');
    fprintf('------------------------------------------------------------\n');

    nSeed = cfg.n_seed_batches;
    seed_summary = repmat(empty_effect_summary(),nSeed,1);
    seed_Mb = zeros(nSeed,1);
    seed_CI_F = zeros(nSeed,2);
    seed_CI_R = zeros(nSeed,2);

    T_s = cfg.horizon_main_min*60;
    dt = cfg.dt_main;

    for iseed=1:nSeed
        seed0 = cfg.seed_seed_base+1000*iseed;

        fprintf('Independent batch %d/%d\n',iseed,nSeed);

        rng(seed0+1,'twister');

        CALseed = simulate_multi_horizon_tdcme( ...
            p,dt,T_s,T_s,cfg.seed_n_cal,[]);

        pooled_cal = [ ...
            CALseed.CIR.Y(:,1); ...
            CALseed.GOU.Y(:,1)];

        seed_Mb(iseed) = max(1,round(quantile( ...
            pooled_cal,cfg.q_main)));

        rng(seed0+2,'twister');

        EVALseed = simulate_multi_horizon_tdcme( ...
            p,dt,T_s,T_s,cfg.seed_n_eval,seed_Mb(iseed));

        [S,dF_i,dR_i] = summarize_effect( ...
            EVALseed.CIR.Y(:,1), ...
            EVALseed.GOU.Y(:,1), ...
            EVALseed.CIR.FPT(:,1), ...
            EVALseed.GOU.FPT(:,1), ...
            seed_Mb(iseed),T_s);

        seed_summary(iseed) = S;

        rng(seed0+3,'twister');

        [bootFs,bootRs] = paired_bootstrap_matrix( ...
            dF_i,dR_i/60,cfg.seed_n_boot);

        seed_CI_F(iseed,:) = prctile(bootFs,[2.5 97.5]);
        seed_CI_R(iseed,:) = prctile(bootRs,[2.5 97.5]);
    end

    seed_table = table( ...
        (1:nSeed)',seed_Mb(:), ...
        [seed_summary.DeltaF]', ...
        seed_CI_F(:,1),seed_CI_F(:,2), ...
        [seed_summary.DeltaRMST_s]'/60, ...
        seed_CI_R(:,1),seed_CI_R(:,2), ...
        [seed_summary.meanY_CIR]', ...
        [seed_summary.meanY_GOU]', ...
        100*[seed_summary.signedMeanDifference]', ...
        [seed_summary.fano_CIR]', ...
        [seed_summary.fano_GOU]', ...
        100*[seed_summary.signedFanoDifference]', ...
        'VariableNames',{ ...
        'seed_batch','Mb','DeltaF', ...
        'DeltaF_CI_low','DeltaF_CI_high', ...
        'DeltaRMST_min','DeltaRMST_CI_low_min', ...
        'DeltaRMST_CI_high_min', ...
        'meanY_CIR','meanY_GammaOU','signedMeanDifference_pct', ...
        'Fano_CIR','Fano_GammaOU','signedFanoDifference_pct'});

    % Inverse-variance pooled effects across independent batches.
    [poolR,poolRlo,poolRhi] = inverse_variance_pool( ...
        seed_table.DeltaRMST_min, ...
        seed_table.DeltaRMST_CI_low_min, ...
        seed_table.DeltaRMST_CI_high_min);

    [poolF,poolFlo,poolFhi] = inverse_variance_pool( ...
        seed_table.DeltaF, ...
        seed_table.DeltaF_CI_low, ...
        seed_table.DeltaF_CI_high);

    pooled_seed_table = table( ...
        poolF,poolFlo,poolFhi, ...
        poolR,poolRlo,poolRhi, ...
        'VariableNames',{ ...
        'pooled_DeltaF','pooled_DeltaF_CI_low', ...
        'pooled_DeltaF_CI_high', ...
        'pooled_DeltaRMST_min', ...
        'pooled_DeltaRMST_CI_low_min', ...
        'pooled_DeltaRMST_CI_high_min'});

    RESULT.seed_replication.table = seed_table;
    RESULT.seed_replication.pooled = pooled_seed_table;

    fprintf('\nIndependent-seed summary:\n');
    disp(seed_table(:,{ ...
        'seed_batch','Mb','DeltaF','DeltaRMST_min', ...
        'DeltaRMST_CI_low_min','DeltaRMST_CI_high_min'}));

    fprintf('Pooled independent-seed effects:\n');
    disp(pooled_seed_table);

    if export_tables
        writetable(seed_table,fullfile(out_dir, ...
            ['TDCME_seed_replication_' lower(run_level) '.csv']));

        writetable(pooled_seed_table,fullfile(out_dir, ...
            ['TDCME_seed_replication_pooled_' ...
            lower(run_level) '.csv']));
    end
end

%% =========================================================================
% 5. Numerical-convergence and seed-replication figure
% =========================================================================
if run_dt_convergence || run_seed_replication
    fig2 = figure( ...
        'Color','w', ...
        'Units','centimeters', ...
        'Position',[1 1 18.3 11.5], ...
        'Renderer','painters');

    tl2 = tiledlayout(fig2,2,2, ...
        'TileSpacing','compact','Padding','compact');

    % a | DeltaF convergence.
    ax = nexttile(tl2,1);

    if run_dt_convergence
        errorbar(ax,dt_table.dt_s,dt_table.DeltaF, ...
            dt_table.DeltaF-dt_table.DeltaF_CI_low, ...
            dt_table.DeltaF_CI_high-dt_table.DeltaF, ...
            'o-','Color',style.gou, ...
            'MarkerFaceColor','w', ...
            'MarkerSize',style.marker_size, ...
            'LineWidth',style.data_line, ...
            'CapSize',3);

        set(ax,'XDir','reverse');
        xline(ax,cfg.dt_main,':', ...
            'Color',style.light_gray, ...
            'HandleVisibility','off');
        yline(ax,0,':','Color',style.light_gray, ...
            'HandleVisibility','off');
        xlabel(ax,'Time step, dt (s)');
        ylabel(ax,'\DeltaF_{60}');
        title(ax,'Numerical convergence of crossing risk');
    else
        axis(ax,'off');
    end

    panel_label(ax,'a',style);
    format_axes_if_visible(ax,style);

    % b | DeltaRMST convergence.
    ax = nexttile(tl2,2);

    if run_dt_convergence
        errorbar(ax,dt_table.dt_s,dt_table.DeltaRMST_min, ...
            dt_table.DeltaRMST_min- ...
                dt_table.DeltaRMST_CI_low_min, ...
            dt_table.DeltaRMST_CI_high_min- ...
                dt_table.DeltaRMST_min, ...
            'o-','Color',style.cir, ...
            'MarkerFaceColor','w', ...
            'MarkerSize',style.marker_size, ...
            'LineWidth',style.data_line, ...
            'CapSize',3);

        set(ax,'XDir','reverse');
        xline(ax,cfg.dt_main,':', ...
            'Color',style.light_gray, ...
            'HandleVisibility','off');
        yline(ax,0,':','Color',style.light_gray, ...
            'HandleVisibility','off');
        xlabel(ax,'Time step, dt (s)');
        ylabel(ax,'\DeltaRMST_{60} (min)');
        title(ax,'Numerical convergence of FPT timing');
    else
        axis(ax,'off');
    end

    panel_label(ax,'b',style);
    format_axes_if_visible(ax,style);

    % c | low-order output audit across dt.
    ax = nexttile(tl2,3);

    if run_dt_convergence
        hold(ax,'on');

        plot(ax,dt_table.dt_s, ...
            abs(dt_table.signedMeanDifference_pct), ...
            'o-','Color',style.cir, ...
            'MarkerFaceColor','w', ...
            'MarkerSize',style.marker_size, ...
            'LineWidth',style.data_line);

        plot(ax,dt_table.dt_s, ...
            abs(dt_table.signedFanoDifference_pct), ...
            's--','Color',style.gou, ...
            'MarkerFaceColor','w', ...
            'MarkerSize',style.marker_size, ...
            'LineWidth',style.data_line);

        set(ax,'XDir','reverse');
        xline(ax,cfg.dt_main,':', ...
            'Color',style.light_gray, ...
            'HandleVisibility','off');

        xlabel(ax,'Time step, dt (s)');
        ylabel(ax,'Absolute CIR-GammaOU difference (%)');
        title(ax,'Low-order output audit across dt');
        legend(ax,{'Mean cumulative ROS','Fano factor'}, ...
            'Location','best','Box','off');
    else
        axis(ax,'off');
    end

    panel_label(ax,'c',style);
    format_axes_if_visible(ax,style);

    % d | independent-seed forest plot.
    ax = nexttile(tl2,4);

    if run_seed_replication
        hold(ax,'on');

        ySeed = (1:nSeed)';

        for iseed=1:nSeed
            draw_horizontal_interval(ax, ...
                seed_table.DeltaRMST_CI_low_min(iseed), ...
                seed_table.DeltaRMST_CI_high_min(iseed), ...
                seed_table.DeltaRMST_min(iseed), ...
                ySeed(iseed),style.cir,style);
        end

        yPool = nSeed+1.3;

        draw_horizontal_interval(ax, ...
            poolRlo,poolRhi,poolR,yPool, ...
            style.dark_gray,style);

        xline(ax,0,':','Color',style.light_gray, ...
            'HandleVisibility','off');

        set(ax, ...
            'YTick',[ySeed;yPool], ...
            'YTickLabel',[ ...
                compose('Seed %d',1:nSeed)'; ...
                "Pooled"], ...
            'YDir','reverse');

        ylim(ax,[0.4 yPool+0.7]);
        xlabel(ax,'\DeltaRMST_{60} (min)');
        ylabel(ax,'');
        title(ax,'Independent-seed replication');
    else
        axis(ax,'off');
    end

    panel_label(ax,'d',style);
    format_axes_if_visible(ax,style);

    if export_figures
        export_figure(fig2,fullfile(out_dir, ...
            ['TDCME_SI_numerical_seed_robustness_' ...
            lower(run_level)]),cfg.export_dpi);
    end
end

%% =========================================================================
% 6. Save compact result structure
% =========================================================================
if save_results
    save(fullfile(out_dir, ...
        ['TDCME_core_robustness_results_' ...
        lower(run_level) '.mat']), ...
        'RESULT','-v7.3');
end

fprintf('\n============================================================\n');
fprintf('TD-CME robustness analysis completed successfully.\n');
fprintf('Results directory:\n%s\n',out_dir);
fprintf('============================================================\n');

diary off;

%% =========================================================================
% Local functions
% =========================================================================

function OUT = simulate_multi_horizon_tdcme( ...
        p,dt,T_max,horizon_s,n,thresholds)

    horizon_s = horizon_s(:)';
    horizon_steps = round(horizon_s/dt);
    n_steps = round(T_max/dt);

    if any(abs(horizon_steps*dt-horizon_s)>1e-10)
        error('Every horizon must lie exactly on the dt grid.');
    end

    thresholds = thresholds(:)';
    n_thresholds = numel(thresholds);
    n_horizon = numel(horizon_s);

    % Shared matched initialization.
    x1_initial = binornd(p.N,p.p_x1_ss,n,1);
    c1_initial = gamrnd(p.shape1,p.theta1,n,1);
    c2_initial = gamrnd(p.shape2,p.theta2,n,1);

    x1C = x1_initial;
    x2C = p.N-x1C;
    yC = zeros(n,1);
    c1C = c1_initial;
    c2C = c2_initial;

    x1G = x1_initial;
    x2G = p.N-x1G;
    yG = zeros(n,1);
    c1G = c1_initial;
    c2G = c2_initial;

    YC = zeros(n,n_horizon);
    YG = zeros(n,n_horizon);

    if n_thresholds>0
        FPTC = nan(n,n_thresholds);
        FPTG = nan(n,n_thresholds);
    else
        FPTC = zeros(n,0);
        FPTG = zeros(n,0);
    end

    horizon_pointer = 1;

    for it=1:n_steps
        t_now = it*dt;

        % ---------------- CIR -------------------------------------------
        c1_old = c1C;
        c2_old = c2C;

        c1C = cir_exact_transition( ...
            c1C,p.kappa1,p.shape1,p.theta1,dt);

        c2C = cir_exact_transition( ...
            c2C,p.kappa2,p.shape2,p.theta2,dt);

        I1C = 0.5*(c1_old+c1C)*dt;
        I2C = 0.5*(c2_old+c2C)*dt;

        pr1C = min(max(1-exp(-I1C),0),1);
        pr2C = min(max(1-exp(-I2C),0),1);

        r1C = binornd(x1C,pr1C);
        r2C = binornd(x2C,pr2C);

        x1C = x1C-r1C+r2C;
        x2C = p.N-x1C;
        yC = yC+r1C;

        % ---------------- GammaOU ---------------------------------------
        [c1G,I1G] = gamma_ou_step_with_integral( ...
            c1G,p.a1,p.kappa1,p.theta1,dt);

        [c2G,I2G] = gamma_ou_step_with_integral( ...
            c2G,p.a2,p.kappa2,p.theta2,dt);

        pr1G = min(max(1-exp(-I1G),0),1);
        pr2G = min(max(1-exp(-I2G),0),1);

        r1G = binornd(x1G,pr1G);
        r2G = binornd(x2G,pr2G);

        x1G = x1G-r1G+r2G;
        x2G = p.N-x1G;
        yG = yG+r1G;

        % First passage for all requested thresholds.
        for ith=1:n_thresholds
            hitC = isnan(FPTC(:,ith)) & yC>=thresholds(ith);
            hitG = isnan(FPTG(:,ith)) & yG>=thresholds(ith);

            FPTC(hitC,ith) = t_now;
            FPTG(hitG,ith) = t_now;
        end

        % Endpoint recording.
        if horizon_pointer<=n_horizon && ...
                it==horizon_steps(horizon_pointer)

            YC(:,horizon_pointer) = yC;
            YG(:,horizon_pointer) = yG;

            horizon_pointer = horizon_pointer+1;
        end
    end

    OUT.horizon_s = horizon_s;
    OUT.thresholds = thresholds;

    OUT.CIR.Y = YC;
    OUT.GOU.Y = YG;

    OUT.CIR.FPT = FPTC;
    OUT.GOU.FPT = FPTG;
end

function [S,dF_i,dR_i] = summarize_effect( ...
        YC,YG,fptC,fptG,threshold,tau)

    eventC = ~isnan(fptC) & fptC<=tau;
    eventG = ~isnan(fptG) & fptG<=tau;

    timeC = restricted_time(fptC,tau);
    timeG = restricted_time(fptG,tau);

    dF_i = double(eventG)-double(eventC);
    dR_i = timeC-timeG;

    meanC = mean(YC);
    meanG = mean(YG);

    varC = var(YC,0);
    varG = var(YG,0);

    fanoC = varC/max(meanC,eps);
    fanoG = varG/max(meanG,eps);

    endpointC = mean(YC>=threshold);
    endpointG = mean(YG>=threshold);

    S = empty_effect_summary();

    S.threshold = threshold;

    S.P_CIR = mean(eventC);
    S.P_GOU = mean(eventG);
    S.DeltaF = mean(dF_i);

    S.RMST_CIR_s = mean(timeC);
    S.RMST_GOU_s = mean(timeG);
    S.DeltaRMST_s = mean(dR_i);

    S.meanY_CIR = meanC;
    S.meanY_GOU = meanG;

    S.fano_CIR = fanoC;
    S.fano_GOU = fanoG;

    S.signedMeanDifference = ...
        (meanG-meanC)/max(0.5*(meanC+meanG),eps);

    S.signedFanoDifference = ...
        (fanoG-fanoC)/max(0.5*(fanoC+fanoG),eps);

    S.identityError_CIR = abs(endpointC-S.P_CIR);
    S.identityError_GOU = abs(endpointG-S.P_GOU);
end

function S = empty_effect_summary()
    S = struct( ...
        'threshold',nan, ...
        'P_CIR',nan, ...
        'P_GOU',nan, ...
        'DeltaF',nan, ...
        'RMST_CIR_s',nan, ...
        'RMST_GOU_s',nan, ...
        'DeltaRMST_s',nan, ...
        'meanY_CIR',nan, ...
        'meanY_GOU',nan, ...
        'fano_CIR',nan, ...
        'fano_GOU',nan, ...
        'signedMeanDifference',nan, ...
        'signedFanoDifference',nan, ...
        'identityError_CIR',nan, ...
        'identityError_GOU',nan);
end

function tRestricted = restricted_time(fpt,tau)
    tRestricted = fpt;
    tRestricted(isnan(tRestricted)) = tau;
    tRestricted = min(tRestricted,tau);
end

function [bootF,bootR] = paired_bootstrap_matrix( ...
        dF_contrib,dR_contrib,n_boot)

    if isvector(dF_contrib)
        dF_contrib = dF_contrib(:);
    end

    if isvector(dR_contrib)
        dR_contrib = dR_contrib(:);
    end

    n = size(dF_contrib,1);
    nCond = size(dF_contrib,2);

    if size(dR_contrib,1)~=n || ...
            size(dR_contrib,2)~=nCond
        error('Bootstrap contribution matrices must have equal dimensions.');
    end

    bootF = zeros(n_boot,nCond);
    bootR = zeros(n_boot,nCond);

    for ib=1:n_boot
        idx = randi(n,n,1);

        bootF(ib,:) = mean(dF_contrib(idx,:),1);
        bootR(ib,:) = mean(dR_contrib(idx,:),1);
    end
end

function validate_time_points(dt,time_points)
    steps = time_points/dt;

    if any(abs(steps-round(steps))>1e-10)
        error('All time points must be exact integer multiples of dt.');
    end
end

function cNext = cir_exact_transition( ...
        c,kappa,shape,theta,dt)

    if dt==0
        cNext = c;
        return;
    end

    decay = exp(-kappa*dt);
    oneMinusDecay = 1-decay;

    scale = 0.5*theta*oneMinusDecay;
    degreesFreedom = 2*shape;

    noncentrality = ...
        2*decay.*max(c,0)./ ...
        (theta*oneMinusDecay);

    if exist('ncx2rnd','file')==2
        cNext = scale.*ncx2rnd( ...
            degreesFreedom,noncentrality);
    else
        poissonMix = poissrnd(noncentrality/2);

        cNext = scale.*gamrnd( ...
            degreesFreedom/2+poissonMix,2);
    end

    cNext = max(cNext,0);
end

function [cNew,I] = gamma_ou_step_with_integral( ...
        cOld,a,kappa,theta,dt)

    if dt==0
        cNew = cOld;
        I = zeros(size(cOld));
        return;
    end

    n = numel(cOld);
    decay = exp(-kappa*dt);

    cNew = cOld.*decay;
    I = cOld.*(1-decay)/kappa;

    jumpNumber = poissrnd(a*dt,n,1);
    totalJumps = sum(jumpNumber);

    if totalJumps>0
        owner = repelem((1:n)',jumpNumber);

        jumpTime = dt*rand(totalJumps,1);
        jumpSize = gamrnd(1,theta,totalJumps,1);

        endpointContribution = ...
            jumpSize.*exp(-kappa*(dt-jumpTime));

        integralContribution = ...
            jumpSize.* ...
            (1-exp(-kappa*(dt-jumpTime)))/kappa;

        cNew = cNew+accumarray( ...
            owner,endpointContribution,[n 1],@sum,0);

        I = I+accumarray( ...
            owner,integralContribution,[n 1],@sum,0);
    end

    cNew = max(cNew,0);
    I = max(I,0);
end

function [pooled,ciLow,ciHigh] = inverse_variance_pool( ...
        estimate,ciLowInput,ciHighInput)

    estimate = estimate(:);
    ciLowInput = ciLowInput(:);
    ciHighInput = ciHighInput(:);

    se = (ciHighInput-ciLowInput)/(2*1.96);
    valid = isfinite(estimate) & isfinite(se) & se>0;

    if ~any(valid)
        pooled = nan;
        ciLow = nan;
        ciHigh = nan;
        return;
    end

    weight = 1./(se(valid).^2);
    pooled = sum(weight.*estimate(valid))/sum(weight);

    pooledSE = sqrt(1/sum(weight));

    ciLow = pooled-1.96*pooledSE;
    ciHigh = pooled+1.96*pooledSE;
end

function draw_horizontal_interval( ...
        ax,lo,hi,estimate,y,color,style)

    plot(ax,[lo hi],[y y],'-', ...
        'Color',color, ...
        'LineWidth',style.data_line);

    plot(ax,lo,y,'|', ...
        'Color',color, ...
        'MarkerSize',6, ...
        'LineWidth',style.data_line);

    plot(ax,hi,y,'|', ...
        'Color',color, ...
        'MarkerSize',6, ...
        'LineWidth',style.data_line);

    plot(ax,estimate,y,'o', ...
        'Color',color, ...
        'MarkerFaceColor','w', ...
        'MarkerSize',style.marker_size+1, ...
        'LineWidth',style.data_line);
end

function style = journal_style()
    style.font_name = 'Arial';

    style.axis_font = 7.0;
    style.label_font = 7.6;
    style.title_font = 8.2;
    style.panel_font = 10.5;

    style.axis_line = 0.75;
    style.data_line = 1.25;
    style.marker_size = 4.8;

    style.cir = hex2rgb('#3C5488');
    style.gou = hex2rgb('#E64B35');

    style.dark_gray = [0.20 0.20 0.20];
    style.light_gray = [0.72 0.72 0.72];
end

function format_axes(ax,style)
    set(ax, ...
        'Box','off', ...
        'TickDir','out', ...
        'TickLength',[0.018 0.018], ...
        'LineWidth',style.axis_line, ...
        'FontName',style.font_name, ...
        'FontSize',style.axis_font, ...
        'Layer','top', ...
        'XColor',[0.12 0.12 0.12], ...
        'YColor',[0.12 0.12 0.12]);

    ax.XLabel.FontName = style.font_name;
    ax.XLabel.FontSize = style.label_font;

    ax.YLabel.FontName = style.font_name;
    ax.YLabel.FontSize = style.label_font;

    ax.Title.FontName = style.font_name;
    ax.Title.FontSize = style.title_font;

    grid(ax,'off');

    try
        if ~isempty(ax.Toolbar)
            ax.Toolbar.Visible = 'off';
        end
    catch
    end

    try
        disableDefaultInteractivity(ax);
    catch
    end
end

function format_axes_if_visible(ax,style)
    if strcmp(ax.Visible,'on')
        format_axes(ax,style);
    end
end

function panel_label(ax,labelText,style)
    text(ax,-0.14,1.08,labelText, ...
        'Units','normalized', ...
        'FontName',style.font_name, ...
        'FontSize',style.panel_font, ...
        'FontWeight','bold', ...
        'Color',[0 0 0], ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','top', ...
        'Clipping','off');
end

function export_figure(fig,stem,dpi)
    drawnow;

    if exist('exportgraphics','file')==2
        exportgraphics(fig,[stem '.pdf'], ...
            'ContentType','vector', ...
            'BackgroundColor','white');

        exportgraphics(fig,[stem '.png'], ...
            'Resolution',dpi, ...
            'BackgroundColor','white');

        exportgraphics(fig,[stem '.tif'], ...
            'Resolution',dpi, ...
            'BackgroundColor','white');
    else
        set(fig,'PaperPositionMode','auto');

        print(fig,[stem '.pdf'], ...
            '-dpdf','-painters','-bestfit');

        print(fig,[stem '.png'], ...
            '-dpng',sprintf('-r%d',dpi));

        print(fig,[stem '.tif'], ...
            '-dtiff',sprintf('-r%d',dpi));
    end
end

function rgb = hex2rgb(hexString)
    hexString = strrep(hexString,'#','');
    rgb = sscanf(hexString,'%2x%2x%2x',[1 3])/255;
end
