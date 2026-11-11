% =========================================================================
% TDCME_RegimeMap.m
%
% Publication-oriented TD-CME condition analysis
% ------------------------------------------------
% Scientific question:
%   Under pointwise matched CIR and GammaOU stochastic inputs, in which
%   timescale-dispersion regimes does higher-order path architecture alter
%   full-state ROS-burden first-passage risk?
%
% Main coordinates:
%
%   Theta_1 = tau_1/T_obs = 1/(kappa_1*T_obs)
%
%   CV_1^2  = Var(c_1)/E(c_1)^2 = 1/shape_1
%
% At every sampled condition:
%   - mu_c1 is fixed;
%   - CIR and GammaOU have the same stationary Gamma marginal;
%   - CIR and GammaOU have the same exponential autocorrelation;
%   - only higher-order temporal path architecture differs.
%
% Publication-oriented improvements
% ---------------------------------
% 1. Regular main grid expanded from 6x6 to 6x8 = 48 conditions:
%
%      Theta_1 = [0.20, 0.35, 0.60, 0.9259, 1.50, 2.50]
%      CV_1^2  = [0.50, 1, 2, 3, 4, 5, 8, 10]
%
%    The added CV_1^2=3 and 4 conditions resolve the transition between
%    low-dispersion and high-dispersion regimes.
%
% 2. Optional local refinement at:
%
%      Theta_1 = [0.45, 0.75, 1.20]
%      CV_1^2  = [5, 8]
%
%    These six extra conditions refine the memory dependence in the
%    high-dispersion region without pretending that they form a regular
%    two-dimensional grid.
%
% 3. Discrete cell heatmaps are used. No contour interpolation is applied.
%    Every colored cell corresponds to an actually simulated condition.
%
% 4. Paired bootstrap confidence intervals are retained.
%
% 5. The primary RMST test is a paired test on restricted-time
%    contributions:
%
%      D_i = min(T_CIR,i,T_obs) - min(T_GOU,i,T_obs).
%
%    Benjamini-Hochberg FDR is applied across the regular-grid cells.
%
% 6. A prespecified mean-output audit tolerance is used:
%
%      |signed relative mean difference at 60 min| <= 3%.
%
%    Cells failing this audit are marked rather than silently interpreted
%    as pure path-architecture effects.
%
% 7. A baseline-fixed absolute threshold is evaluated in SI, together with
%    absolute-threshold event-rate maps.
%
% 8. Optional multi-seed, high-precision replication is implemented for
%    selected baseline, transition, and high-effect conditions.
%
% Main manuscript figure
% ----------------------
%   a. Discrete standardized DeltaRMST_60 regime map.
%      Circle: FDR q<0.05 and mean audit passed.
%      x: mean audit failed.
%      Star: manuscript baseline.
%
%   b. Memory slice at CV_1^2=5, including local refinement.
%
%   c. Dispersion slice at Theta_1=baseline.
%
% SI figure 1
% -----------
%   a. Standardized DeltaF_60 map
%   b. Mean-output audit map
%   c. Fano-output difference map
%   d. Cell-specific Q0.90 threshold map
%   e. Baseline-fixed-threshold DeltaRMST_60 map
%   f. Baseline-fixed-threshold pooled event-rate map
%
% SI figure 2
% -----------
%   a. DeltaRMST bootstrap CI half-width
%   b. FDR q-value map
%   c. Refined memory slices at CV_1^2=5 and 8
%   d. Multi-seed key-condition replication
%
% Strict TD-CME kernel retained
% -----------------------------
% N = 1200
%
% Baseline c1:
%   kappa1 = 3.0e-4 s^-1
%   shape1 = 0.20
%   mu_c1  = 1.5e-5 s^-1
%
% Fixed c2:
%   kappa2 = 4.0e-4 s^-1
%   shape2 = 3.50
%   mu_c2  = 1.7e-5 s^-1
%
% dt = 2 s
% T_obs = 3600 s = 60 min
%
% Numerical implementation:
%   - paired stationary-binomial Fe2+/Fe3+ initialization
%   - paired stationary-Gamma c1/c2 initialization
%   - exact CIR noncentral-chi-square transition
%   - exact GammaOU exponential decay and jump-time integral
%   - integrated-hazard binomial reaction updates
%   - independent calibration and evaluation samples
%   - paired bootstrap
%
% Required toolbox
% ----------------
% Statistics and Machine Learning Toolbox
%
% Optional toolbox
% ----------------
% Parallel Computing Toolbox
%
% Windows parallel-pool note
% --------------------------
% The safest default is serial execution. A local parallel pool may fail
% when MATLAB paths contain Chinese or other non-ASCII characters. Move the
% project to an ASCII-only folder and use an ASCII-only JobStorageLocation
% before enabling parallel execution.
%
% Recommended workflow
% --------------------
%   run_level = 'smoke' -> 'pilot' -> 'final'
%
% Final runs are computationally expensive.
% =========================================================================

clear; clc; close all;

%% =========================================================================
% 0. Run configuration
% =========================================================================
run_level = 'final';  % 'smoke' | 'pilot' | 'final'

make_main_figure = true;
make_SI_figure_1 = true;
make_SI_figure_2 = true;

export_figures   = true;
export_tables    = true;
save_results_mat = true;
save_matlab_fig  = false;

% Serial is deliberately the default after the previously observed Windows
% character-encoding error during parpool startup.
cfg.use_parallel = false;
cfg.max_workers  = 4;
cfg.parallel_job_storage = '';

cfg.export_dpi = 600;

% Mean-output equivalence audit. This is descriptive and prespecified.
cfg.mean_audit_tolerance = 0.03;

% Main statistical threshold.
cfg.fdr_alpha = 0.05;

% Fixed seeds.
cfg.seed_baseline_calibration = 20262001;
cfg.seed_cell_base            = 20262100;
cfg.seed_bootstrap_base       = 20262200;
cfg.seed_refinement_base      = 20262300;
cfg.seed_keypoint_base        = 20262400;

% Manuscript baseline.
baseline = struct();
baseline.kappa1 = 3.0e-4;
baseline.shape1 = 0.20;

T_obs = 3600;
baseline.Theta1 = 1/(baseline.kappa1*T_obs);
baseline.CV2_1  = 1/baseline.shape1;

switch lower(run_level)
    case 'smoke'
        % Fast structural test. Includes the transition region.
        cfg.Theta1_vec = [0.35, 0.60, baseline.Theta1, 1.50];
        cfg.CV2_vec    = [1, 2, 3, 4, baseline.CV2_1, 8];

        cfg.n_cal  = 250;
        cfg.n_eval = 500;
        cfg.n_boot = 150;

        cfg.run_local_refinement = false;
        cfg.run_keypoint_replication = false;

        cfg.n_cal_refine  = 250;
        cfg.n_eval_refine = 500;
        cfg.n_boot_refine = 150;

        cfg.key_n_seeds = 2;
        cfg.key_n_cal   = 400;
        cfg.key_n_eval  = 800;
        cfg.key_n_boot  = 200;

    case 'pilot'
        % Full 48-cell grid at moderate precision.
        cfg.Theta1_vec = [ ...
            0.20, 0.35, 0.60, baseline.Theta1, 1.50, 2.50];

        cfg.CV2_vec = [ ...
            0.50, 1.00, 2.00, 3.00, 4.00, ...
            baseline.CV2_1, 8.00, 10.0];

        cfg.n_cal  = 800;
        cfg.n_eval = 1600;
        cfg.n_boot = 400;

        cfg.run_local_refinement = true;
        cfg.run_keypoint_replication = false;

        cfg.n_cal_refine  = 800;
        cfg.n_eval_refine = 1600;
        cfg.n_boot_refine = 400;

        cfg.key_n_seeds = 3;
        cfg.key_n_cal   = 1200;
        cfg.key_n_eval  = 2500;
        cfg.key_n_boot  = 500;

    case 'final'
        % Publication-oriented 48-cell regular grid.
        cfg.Theta1_vec = [ ...
            0.20, 0.35, 0.60, baseline.Theta1, 1.50, 2.50];

        cfg.CV2_vec = [ ...
            0.50, 1.00, 2.00, 3.00, 4.00, ...
            baseline.CV2_1, 8.00, 10.0];

        cfg.n_cal  = 2000;
        cfg.n_eval = 4000;
        cfg.n_boot = 1000;

        cfg.run_local_refinement = true;
        cfg.run_keypoint_replication = true;

        cfg.n_cal_refine  = 2000;
        cfg.n_eval_refine = 4000;
        cfg.n_boot_refine = 1000;

        % Three independent high-precision batches for each key condition.
        cfg.key_n_seeds = 3;
        cfg.key_n_cal   = 3000;
        cfg.key_n_eval  = 6000;
        cfg.key_n_boot  = 800;

    otherwise
        error('run_level must be ''smoke'', ''pilot'', or ''final''.');
end

% Local off-grid refinement.
cfg.refine_Theta1_vec = [0.45, 0.75, 1.20];
cfg.refine_CV2_vec    = [baseline.CV2_1, 8.00];

style = journal_style();

%% =========================================================================
% 1. Fixed TD-CME parameters
% =========================================================================
p = struct();

p.N = 1200;

% c1 mean fixed across all conditions.
p.mu_c1 = 1.5e-5;

% c2 fixed at the manuscript baseline.
p.kappa2 = 4.0e-4;
p.shape2 = 3.50;
p.mu_c2  = 1.7e-5;
p.theta2 = p.mu_c2/p.shape2;
p.a2     = p.shape2*p.kappa2;

% Shared stationary iron initialization.
p.p_x1_ss = p.mu_c2/(p.mu_c1+p.mu_c2);
p.x1_ss   = p.N*p.p_x1_ss;
p.x2_ss   = p.N-p.x1_ss;

sim = struct();
sim.dt          = 2.0;
sim.T_total     = T_obs;
sim.q_threshold = 0.90;

validate_inputs(sim,cfg,baseline);

fprintf('============================================================\n');
fprintf('TD-CME publication-ready timescale-dispersion regime map\n');
fprintf('============================================================\n');
fprintf('Run level: %s\n',run_level);
fprintf('Regular grid: %d Theta values x %d CV^2 values = %d cells\n', ...
    numel(cfg.Theta1_vec),numel(cfg.CV2_vec), ...
    numel(cfg.Theta1_vec)*numel(cfg.CV2_vec));
fprintf('n_cal=%d, n_eval=%d, n_boot=%d per regular cell\n', ...
    cfg.n_cal,cfg.n_eval,cfg.n_boot);
fprintf('Mean audit tolerance = %.1f %%\n', ...
    100*cfg.mean_audit_tolerance);
fprintf('FDR alpha = %.3f\n',cfg.fdr_alpha);
fprintf('Baseline Theta_1 = %.6f\n',baseline.Theta1);
fprintf('Baseline CV_1^2  = %.6f\n',baseline.CV2_1);
fprintf('Observation horizon = %.1f min; dt = %.1f s\n', ...
    sim.T_total/60,sim.dt);
fprintf('Local refinement enabled: %d\n',cfg.run_local_refinement);
fprintf('Keypoint replication enabled: %d\n', ...
    cfg.run_keypoint_replication);
fprintf('============================================================\n\n');

%% =========================================================================
% 2. Robust parallel setup
% =========================================================================
useParallel = configure_parallel_execution(cfg);

%% =========================================================================
% 3. Baseline-fixed absolute threshold
% =========================================================================
baselineParams = make_cell_parameters( ...
    p,sim,baseline.Theta1,baseline.CV2_1);

nBaseCal = max(cfg.n_cal,2000);

rng(cfg.seed_baseline_calibration,'twister');

fprintf('Calibrating baseline-fixed absolute threshold...\n');

CAL_BASE = simulate_calibration( ...
    baselineParams,sim,nBaseCal);

M_abs = max(1,round(quantile( ...
    [CAL_BASE.CIR.Y60;CAL_BASE.GOU.Y60], ...
    sim.q_threshold)));

fprintf('Baseline-fixed absolute threshold M_abs = %d\n\n',M_abs);

%% =========================================================================
% 4. Regular-grid specifications
% =========================================================================
nTheta = numel(cfg.Theta1_vec);
nCV2   = numel(cfg.CV2_vec);
nCells = nTheta*nCV2;

regularSpec = build_cell_specifications( ...
    cfg.Theta1_vec,cfg.CV2_vec, ...
    cfg.seed_cell_base,cfg.seed_bootstrap_base);

%% =========================================================================
% 5. Regular-grid simulation
% =========================================================================
regularResults = cell(nCells,1);

if useParallel
    parfor icell=1:nCells
        regularResults{icell} = run_one_regime_cell( ...
            p,sim,cfg.n_cal,cfg.n_eval,cfg.n_boot, ...
            regularSpec(icell),M_abs,cfg.mean_audit_tolerance);
    end
else
    for icell=1:nCells
        fprintf('Regular cell %d/%d: Theta=%.4g, CV^2=%.4g\n', ...
            icell,nCells, ...
            regularSpec(icell).Theta1, ...
            regularSpec(icell).CV2_1);

        regularResults{icell} = run_one_regime_cell( ...
            p,sim,cfg.n_cal,cfg.n_eval,cfg.n_boot, ...
            regularSpec(icell),M_abs,cfg.mean_audit_tolerance);
    end
end

fprintf('\nRegular-grid simulation completed.\n\n');

%% =========================================================================
% 6. Assemble regular-grid matrices
% =========================================================================
M = assemble_regular_matrices( ...
    regularResults,nCV2,nTheta);

% Benjamini-Hochberg correction across the regular-grid RMST tests.
[M.q_RMST,M.reject_RMST] = ...
    benjamini_hochberg(M.p_RMST,cfg.fdr_alpha);

M.FDR_and_mean_pass = ...
    M.reject_RMST & M.mean_audit_pass;

M.nominal_CI_excludes_zero = ...
    (M.CIlo_RMST>0 & M.CIhi_RMST>0) | ...
    (M.CIlo_RMST<0 & M.CIhi_RMST<0);

%% =========================================================================
% 7. Optional local memory refinement
% =========================================================================
refinementResults = {};
refinementTable = table();

if cfg.run_local_refinement
    refineSpec = build_cell_specifications( ...
        cfg.refine_Theta1_vec,cfg.refine_CV2_vec, ...
        cfg.seed_refinement_base, ...
        cfg.seed_refinement_base+5000);

    nRefine = numel(refineSpec);
    refinementResults = cell(nRefine,1);

    fprintf('Running %d local-refinement conditions...\n',nRefine);

    if useParallel
        parfor icell=1:nRefine
            refinementResults{icell} = run_one_regime_cell( ...
                p,sim, ...
                cfg.n_cal_refine,cfg.n_eval_refine,cfg.n_boot_refine, ...
                refineSpec(icell),M_abs,cfg.mean_audit_tolerance);
        end
    else
        for icell=1:nRefine
            fprintf('Refinement cell %d/%d: Theta=%.4g, CV^2=%.4g\n', ...
                icell,nRefine, ...
                refineSpec(icell).Theta1, ...
                refineSpec(icell).CV2_1);

            refinementResults{icell} = run_one_regime_cell( ...
                p,sim, ...
                cfg.n_cal_refine,cfg.n_eval_refine,cfg.n_boot_refine, ...
                refineSpec(icell),M_abs,cfg.mean_audit_tolerance);
        end
    end

    refinementTable = results_to_table(refinementResults);

    fprintf('Local-refinement simulation completed.\n\n');
end

%% =========================================================================
% 8. Optional multi-seed key-condition replication
% =========================================================================
keyReplication = struct();
keyReplicationTable = table();

if cfg.run_keypoint_replication
    % Selected conditions:
    %   1. manuscript baseline;
    %   2. dispersion transition CV^2=3;
    %   3. dispersion transition CV^2=4;
    %   4. high-dispersion, short-memory candidate;
    %   5. high-dispersion extreme candidate.
    keyPoints = struct( ...
        'label',{ ...
            'Baseline', ...
            'Transition CV2=3', ...
            'Transition CV2=4', ...
            'High dispersion', ...
            'High dispersion extreme'}, ...
        'Theta1',{ ...
            baseline.Theta1, ...
            baseline.Theta1, ...
            baseline.Theta1, ...
            0.35, ...
            0.35}, ...
        'CV2_1',{ ...
            baseline.CV2_1, ...
            3.00, ...
            4.00, ...
            8.00, ...
            10.0});

    fprintf('Running multi-seed key-condition replication...\n');

    keyReplication = run_keypoint_replication( ...
        p,sim,cfg,keyPoints,M_abs);

    keyReplicationTable = keyReplication.seed_table;

    fprintf('Key-condition replication completed.\n\n');
end

%% =========================================================================
% 9. Baseline and transition summaries
% =========================================================================
[~,idxThetaBase] = min(abs(cfg.Theta1_vec-baseline.Theta1));
[~,idxCVBase]    = min(abs(cfg.CV2_vec-baseline.CV2_1));

fprintf('Regular-grid baseline cell:\n');
fprintf('  Theta_1 = %.6f\n',cfg.Theta1_vec(idxThetaBase));
fprintf('  CV_1^2  = %.6f\n',cfg.CV2_vec(idxCVBase));
fprintf('  M_b = %.0f\n',M.Mb(idxCVBase,idxThetaBase));
fprintf('  DeltaRMST_60 = %.4f min [%.4f, %.4f]\n', ...
    M.DeltaRMST(idxCVBase,idxThetaBase), ...
    M.CIlo_RMST(idxCVBase,idxThetaBase), ...
    M.CIhi_RMST(idxCVBase,idxThetaBase));
fprintf('  paired p = %.4g; FDR q = %.4g\n', ...
    M.p_RMST(idxCVBase,idxThetaBase), ...
    M.q_RMST(idxCVBase,idxThetaBase));
fprintf('  signed mean difference = %.3f %%\n\n', ...
    M.signedMeanPct(idxCVBase,idxThetaBase));

%% =========================================================================
% 10. Main manuscript figure
% =========================================================================
if make_main_figure
    figMain = figure( ...
        'Color','w', ...
        'Units','centimeters', ...
        'Position',[1.0 1.0 style.main_width_cm style.main_height_cm], ...
        'Renderer','painters');

    tl = tiledlayout(figMain,1,3, ...
        'TileSpacing','compact','Padding','compact');

    % ------------------------------------------------------------------ a
    ax1 = nexttile(tl,1);

    limRMST = symmetric_color_limit(M.DeltaRMST,0.05);

    plot_discrete_regime_map( ...
        ax1,cfg.Theta1_vec,cfg.CV2_vec,M.DeltaRMST, ...
        [-limRMST limRMST],style.diverging_map, ...
        '\DeltaRMST_{60} (min)', ...
        'Timescale-dispersion regime', ...
        baseline,style);

    overlay_regime_markers( ...
        ax1,cfg.Theta1_vec,cfg.CV2_vec, ...
        M.FDR_and_mean_pass, ...
        ~M.mean_audit_pass, ...
        baseline,style);

    add_panel_label(ax1,'a',style);

    text(ax1,0.02,0.025, ...
        sprintf(['Positive: earlier GammaOU crossing\n' ...
                 'circle: FDR q<%.2f + mean audit; ' ...
                 'x: |mean difference|>%.0f%%'], ...
        cfg.fdr_alpha,100*cfg.mean_audit_tolerance), ...
        'Units','normalized', ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','bottom', ...
        'FontName',style.font_name, ...
        'FontSize',style.annotation_font, ...
        'Color',style.dark_gray);

    % ------------------------------------------------------------------ b
    ax2 = nexttile(tl,2);
    hold(ax2,'on');

    [thetaSlice,rmstTheta,loTheta,hiTheta] = ...
        extract_memory_slice( ...
        cfg.Theta1_vec,cfg.CV2_vec,M, ...
        baseline.CV2_1,refinementTable);

    errorbar(ax2,thetaSlice,rmstTheta, ...
        rmstTheta-loTheta,hiTheta-rmstTheta, ...
        'o-','Color',style.cir, ...
        'MarkerFaceColor','w', ...
        'MarkerSize',style.marker_size, ...
        'LineWidth',style.data_line, ...
        'CapSize',3);

    yline(ax2,0,':', ...
        'Color',style.light_gray, ...
        'LineWidth',style.reference_line);

    xline(ax2,baseline.Theta1,':', ...
        'Color',style.light_gray, ...
        'LineWidth',style.reference_line);

    xlabel(ax2,'Memory ratio, \Theta_1');
    ylabel(ax2,'\DeltaRMST_{60} (min)');

    title(ax2,sprintf('Memory slice at CV_1^2=%.0f',baseline.CV2_1), ...
        'FontWeight','normal', ...
        'FontSize',style.title_font);

    apply_journal_axes(ax2,style);
    add_panel_label(ax2,'b',style);

    text(ax2,0.03,0.04, ...
        'Refined points: \Theta_1=0.45, 0.75 and 1.20', ...
        'Units','normalized', ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','bottom', ...
        'FontName',style.font_name, ...
        'FontSize',style.annotation_font, ...
        'Color',style.dark_gray);

    % ------------------------------------------------------------------ c
    ax3 = nexttile(tl,3);
    hold(ax3,'on');

    [cvSlice,rmstCV,loCV,hiCV] = ...
        extract_dispersion_slice( ...
        cfg.Theta1_vec,cfg.CV2_vec,M,baseline.Theta1);

    errorbar(ax3,cvSlice,rmstCV, ...
        rmstCV-loCV,hiCV-rmstCV, ...
        's--','Color',style.gou, ...
        'MarkerFaceColor','w', ...
        'MarkerSize',style.marker_size, ...
        'LineWidth',style.data_line, ...
        'CapSize',3);

    yline(ax3,0,':', ...
        'Color',style.light_gray, ...
        'LineWidth',style.reference_line);

    xline(ax3,baseline.CV2_1,':', ...
        'Color',style.light_gray, ...
        'LineWidth',style.reference_line);

    xlabel(ax3,'Stationary dispersion, CV_1^2');
    ylabel(ax3,'\DeltaRMST_{60} (min)');

    title(ax3,sprintf('Dispersion slice at \\Theta_1=%.2f', ...
        baseline.Theta1), ...
        'FontWeight','normal', ...
        'FontSize',style.title_font);

    apply_journal_axes(ax3,style);
    add_panel_label(ax3,'c',style);

    text(ax3,0.03,0.04, ...
        'Transition region resolved by CV_1^2=3 and 4', ...
        'Units','normalized', ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','bottom', ...
        'FontName',style.font_name, ...
        'FontSize',style.annotation_font, ...
        'Color',style.dark_gray);

    drawnow;

    if export_figures
        export_publication_figure( ...
            figMain, ...
            sprintf('TDCME_TimescaleDispersion_PUBLICATION_MAIN_%s', ...
            lower(run_level)), ...
            cfg.export_dpi,save_matlab_fig);
    end
end

%% =========================================================================
% 11. SI figure 1: complete audits
% =========================================================================
if make_SI_figure_1
    figSI1 = figure( ...
        'Color','w', ...
        'Units','centimeters', ...
        'Position',[1.0 1.0 style.si_width_cm style.si_height_cm], ...
        'Renderer','painters');

    tlSI1 = tiledlayout(figSI1,2,3, ...
        'TileSpacing','compact','Padding','compact');

    % ------------------------------------------------------------------ a
    axS1 = nexttile(tlSI1,1);

    limF = symmetric_color_limit(M.DeltaF,0.005);

    plot_discrete_regime_map( ...
        axS1,cfg.Theta1_vec,cfg.CV2_vec,M.DeltaF, ...
        [-limF limF],style.diverging_map, ...
        '\DeltaF_{60}', ...
        'Crossing-probability shift', ...
        baseline,style);

    overlay_regime_markers( ...
        axS1,cfg.Theta1_vec,cfg.CV2_vec, ...
        M.nominal_CI_excludes_zero, ...
        false(size(M.mean_audit_pass)), ...
        baseline,style);

    add_panel_label(axS1,'a',style);

    % ------------------------------------------------------------------ b
    axS2 = nexttile(tlSI1,2);

    limMean = symmetric_color_limit(M.signedMeanPct,0.5);

    plot_discrete_regime_map( ...
        axS2,cfg.Theta1_vec,cfg.CV2_vec,M.signedMeanPct, ...
        [-limMean limMean],style.diverging_map, ...
        'Signed mean difference at 60 min (%)', ...
        'Mean-output audit', ...
        baseline,style);

    overlay_regime_markers( ...
        axS2,cfg.Theta1_vec,cfg.CV2_vec, ...
        false(size(M.mean_audit_pass)), ...
        ~M.mean_audit_pass, ...
        baseline,style);

    add_panel_label(axS2,'b',style);

    % ------------------------------------------------------------------ c
    axS3 = nexttile(tlSI1,3);

    limFano = symmetric_color_limit(M.signedFanoPct,2);

    plot_discrete_regime_map( ...
        axS3,cfg.Theta1_vec,cfg.CV2_vec,M.signedFanoPct, ...
        [-limFano limFano],style.diverging_map, ...
        'Signed Fano difference at 60 min (%)', ...
        'Second-order output divergence', ...
        baseline,style);

    add_panel_label(axS3,'c',style);

    % ------------------------------------------------------------------ d
    axS4 = nexttile(tlSI1,4);

    plot_discrete_regime_map( ...
        axS4,cfg.Theta1_vec,cfg.CV2_vec,M.Mb, ...
        [min(M.Mb(:)) max(M.Mb(:))],style.sequential_map, ...
        'Cell-specific M_b', ...
        'Cell-specific pooled Q_{0.90} threshold', ...
        baseline,style);

    add_panel_label(axS4,'d',style);

    % ------------------------------------------------------------------ e
    axS5 = nexttile(tlSI1,5);

    limAbs = symmetric_color_limit(M.DeltaRMST_abs,0.05);

    plot_discrete_regime_map( ...
        axS5,cfg.Theta1_vec,cfg.CV2_vec,M.DeltaRMST_abs, ...
        [-limAbs limAbs],style.diverging_map, ...
        '\DeltaRMST_{60} at fixed M_{abs} (min)', ...
        'Baseline-fixed absolute-threshold timing', ...
        baseline,style);

    add_panel_label(axS5,'e',style);

    text(axS5,0.02,0.035, ...
        sprintf('Baseline-fixed M_{abs}=%d',M_abs), ...
        'Units','normalized', ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','bottom', ...
        'FontName',style.font_name, ...
        'FontSize',style.annotation_font, ...
        'Color',style.dark_gray);

    % ------------------------------------------------------------------ f
    axS6 = nexttile(tlSI1,6);

    plot_discrete_regime_map( ...
        axS6,cfg.Theta1_vec,cfg.CV2_vec,M.eventRate_abs_pooled, ...
        [0 max(M.eventRate_abs_pooled(:))],style.sequential_map, ...
        'Pooled event probability at fixed M_{abs}', ...
        'Absolute-threshold event-rate audit', ...
        baseline,style);

    add_panel_label(axS6,'f',style);

    drawnow;

    if export_figures
        export_publication_figure( ...
            figSI1, ...
            sprintf('TDCME_TimescaleDispersion_PUBLICATION_SI1_%s', ...
            lower(run_level)), ...
            cfg.export_dpi,save_matlab_fig);
    end
end

%% =========================================================================
% 12. SI figure 2: uncertainty, FDR, refinement, replication
% =========================================================================
if make_SI_figure_2
    figSI2 = figure( ...
        'Color','w', ...
        'Units','centimeters', ...
        'Position',[1.0 1.0 style.si_width_cm style.si_height_cm], ...
        'Renderer','painters');

    tlSI2 = tiledlayout(figSI2,2,2, ...
        'TileSpacing','compact','Padding','compact');

    % ------------------------------------------------------------------ a
    axU1 = nexttile(tlSI2,1);

    plot_discrete_regime_map( ...
        axU1,cfg.Theta1_vec,cfg.CV2_vec,M.CIhalf_RMST, ...
        [0 max(M.CIhalf_RMST(:))],style.sequential_map, ...
        '95% CI half-width (min)', ...
        'Uncertainty of \DeltaRMST_{60}', ...
        baseline,style);

    add_panel_label(axU1,'a',style);

    % ------------------------------------------------------------------ b
    axU2 = nexttile(tlSI2,2);

    qDisplay = -log10(max(M.q_RMST,1e-6));

    plot_discrete_regime_map( ...
        axU2,cfg.Theta1_vec,cfg.CV2_vec,qDisplay, ...
        [0 max(qDisplay(:))],style.sequential_map, ...
        '-log_{10}(FDR q)', ...
        'Multiple-comparison-adjusted evidence', ...
        baseline,style);

    overlay_regime_markers( ...
        axU2,cfg.Theta1_vec,cfg.CV2_vec, ...
        M.reject_RMST, ...
        false(size(M.reject_RMST)), ...
        baseline,style);

    add_panel_label(axU2,'b',style);

    % ------------------------------------------------------------------ c
    axU3 = nexttile(tlSI2,3);
    hold(axU3,'on');

    plot_refined_memory_slices( ...
        axU3,cfg,M,refinementTable,baseline,style);

    yline(axU3,0,':', ...
        'Color',style.light_gray, ...
        'LineWidth',style.reference_line);

    xlabel(axU3,'Memory ratio, \Theta_1');
    ylabel(axU3,'\DeltaRMST_{60} (min)');
    title(axU3,'Locally refined high-dispersion memory slices', ...
        'FontWeight','normal', ...
        'FontSize',style.title_font);

    legend(axU3,'Location','best','Box','off', ...
        'FontSize',style.legend_font);

    apply_journal_axes(axU3,style);
    add_panel_label(axU3,'c',style);

    % ------------------------------------------------------------------ d
    axU4 = nexttile(tlSI2,4);
    hold(axU4,'on');

    if cfg.run_keypoint_replication
        plot_keypoint_replication( ...
            axU4,keyReplication,style);
    else
        axis(axU4,'off');

        text(axU4,0.5,0.5, ...
            ['Multi-seed key-condition replication is disabled in this ' ...
             'run level. Enable cfg.run_keypoint_replication for final.'], ...
            'Units','normalized', ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','middle', ...
            'FontName',style.font_name, ...
            'FontSize',style.label_font, ...
            'Color',style.dark_gray);
    end

    title(axU4,'Independent-seed key-condition replication', ...
        'FontWeight','normal', ...
        'FontSize',style.title_font);

    add_panel_label(axU4,'d',style);

    drawnow;

    if export_figures
        export_publication_figure( ...
            figSI2, ...
            sprintf('TDCME_TimescaleDispersion_PUBLICATION_SI2_%s', ...
            lower(run_level)), ...
            cfg.export_dpi,save_matlab_fig);
    end
end

%% =========================================================================
% 13. Export tables
% =========================================================================
regularTable = results_to_table(regularResults);

% Add FDR and matrix-derived flags using regular cell ordering.
qVector = nan(height(regularTable),1);
rejectVector = false(height(regularTable),1);
meanPassVector = false(height(regularTable),1);

for ir=1:height(regularTable)
    iCV = regularTable.iCV2(ir);
    iTh = regularTable.iTheta(ir);

    qVector(ir) = M.q_RMST(iCV,iTh);
    rejectVector(ir) = M.reject_RMST(iCV,iTh);
    meanPassVector(ir) = M.mean_audit_pass(iCV,iTh);
end

regularTable.FDR_q_RMST = qVector;
regularTable.FDR_reject_RMST = rejectVector;
regularTable.mean_audit_pass = meanPassVector;

parameterTable = table( ...
    string(run_level), ...
    p.N,p.mu_c1, ...
    p.kappa2,p.shape2,p.mu_c2,p.theta2,p.a2, ...
    sim.dt,sim.T_total,sim.q_threshold, ...
    baseline.Theta1,baseline.CV2_1,M_abs, ...
    cfg.mean_audit_tolerance,cfg.fdr_alpha, ...
    cfg.n_cal,cfg.n_eval,cfg.n_boot, ...
    cfg.run_local_refinement,cfg.run_keypoint_replication, ...
    'VariableNames',{ ...
    'run_level','N','mu_c1', ...
    'kappa2','shape2','mu_c2','theta2','a2', ...
    'dt_s','T_total_s','quantile_threshold', ...
    'baseline_Theta1','baseline_CV1_squared','baseline_fixed_Mabs', ...
    'mean_audit_tolerance','FDR_alpha', ...
    'n_cal_per_regular_cell','n_eval_per_regular_cell', ...
    'n_boot_per_regular_cell', ...
    'local_refinement_enabled','keypoint_replication_enabled'});

if export_tables
    writetable(regularTable, ...
        sprintf('TDCME_TimescaleDispersion_PUBLICATION_regular_cells_%s.csv', ...
        lower(run_level)));

    writetable(parameterTable, ...
        sprintf('TDCME_TimescaleDispersion_PUBLICATION_parameters_%s.csv', ...
        lower(run_level)));

    if ~isempty(refinementTable)
        writetable(refinementTable, ...
            sprintf('TDCME_TimescaleDispersion_PUBLICATION_refinement_%s.csv', ...
            lower(run_level)));
    end

    if ~isempty(keyReplicationTable)
        writetable(keyReplicationTable, ...
            sprintf('TDCME_TimescaleDispersion_PUBLICATION_key_seeds_%s.csv', ...
            lower(run_level)));

        writetable(keyReplication.summary_table, ...
            sprintf('TDCME_TimescaleDispersion_PUBLICATION_key_summary_%s.csv', ...
            lower(run_level)));
    end
end

%% =========================================================================
% 14. Save MAT results
% =========================================================================
if save_results_mat
    SAVE = struct();

    SAVE.parameters = p;
    SAVE.simulation = sim;
    SAVE.configuration = cfg;
    SAVE.baseline = baseline;
    SAVE.baseline_fixed_threshold = M_abs;

    SAVE.regular_spec = regularSpec;
    SAVE.regular_results = regularResults;
    SAVE.regular_matrices = M;
    SAVE.regular_table = regularTable;

    SAVE.refinement_results = refinementResults;
    SAVE.refinement_table = refinementTable;

    SAVE.key_replication = keyReplication;
    SAVE.parameter_table = parameterTable;

    save(sprintf( ...
        'TDCME_TimescaleDispersion_PUBLICATION_results_%s.mat', ...
        lower(run_level)), ...
        'SAVE','-v7.3');
end

fprintf('\nPublication-oriented regime-map analysis completed.\n');
fprintf(['Interpret precise boundaries only after final multi-seed checks. ' ...
    'The discrete map supports regional, not exact phase-boundary, claims.\n']);

%% =========================================================================
% Local functions
% =========================================================================

function validate_inputs(sim,cfg,baseline)
    if sim.dt<=0 || sim.T_total<=0
        error('dt and T_total must be positive.');
    end

    if abs(sim.T_total/sim.dt-round(sim.T_total/sim.dt))>1e-12
        error('T_total must be divisible by dt.');
    end

    if any(diff(cfg.Theta1_vec)<=0)
        error('Theta1_vec must be strictly increasing.');
    end

    if any(diff(cfg.CV2_vec)<=0)
        error('CV2_vec must be strictly increasing.');
    end

    if ~any(abs(cfg.Theta1_vec-baseline.Theta1)<1e-10)
        error('The manuscript baseline Theta1 must be included in the grid.');
    end

    if ~any(abs(cfg.CV2_vec-baseline.CV2_1)<1e-10)
        error('The manuscript baseline CV2 must be included in the grid.');
    end

    if cfg.n_cal<100 || cfg.n_eval<100
        error('Calibration/evaluation sample sizes are too small.');
    end
end

function useParallel = configure_parallel_execution(cfg)
    useParallel = false;

    if ~cfg.use_parallel
        fprintf('Serial execution enabled.\n\n');
        return;
    end

    try
        if ~license('test','Distrib_Computing_Toolbox')
            fprintf(['Parallel Computing Toolbox unavailable; ' ...
                'serial execution enabled.\n\n']);
            return;
        end

        if any(double(pwd)>127)
            fprintf(2,[ ...
                'Parallel execution disabled because the current folder ' ...
                'contains non-ASCII characters:\n  %s\n\n'],pwd);
            return;
        end

        pool = gcp('nocreate');

        if isempty(pool)
            cluster = parcluster('local');

            if ~isempty(cfg.parallel_job_storage)
                storageFolder = cfg.parallel_job_storage;
            else
                storageFolder = fullfile(tempdir,'MATLAB_TDCME_parallel_jobs');
            end

            if any(double(storageFolder)>127)
                error('Parallel JobStorageLocation must contain only ASCII characters.');
            end

            if ~exist(storageFolder,'dir')
                [ok,msg] = mkdir(storageFolder);

                if ~ok
                    error('Cannot create JobStorageLocation: %s',msg);
                end
            end

            cluster.JobStorageLocation = storageFolder;
            pool = parpool(cluster,cfg.max_workers);
        end

        useParallel = true;

        fprintf('Parallel execution enabled with %d workers.\n\n', ...
            pool.NumWorkers);

    catch ME
        fprintf(2,[ ...
            'Parallel startup failed; serial execution will be used.\n' ...
            'MATLAB message: %s\n\n'],ME.message);

        useParallel = false;
    end
end

function spec = build_cell_specifications( ...
        ThetaVec,CV2Vec,seedCellBase,seedBootBase)

    nTheta = numel(ThetaVec);
    nCV = numel(CV2Vec);
    n = nTheta*nCV;

    spec = repmat(struct( ...
        'linear_index',[], ...
        'iTheta',[], ...
        'iCV2',[], ...
        'Theta1',[], ...
        'CV2_1',[], ...
        'seed_cal',[], ...
        'seed_eval',[], ...
        'seed_boot',[]),n,1);

    kk = 0;

    for iCV=1:nCV
        for iTh=1:nTheta
            kk = kk+1;

            spec(kk).linear_index = kk;
            spec(kk).iTheta = iTh;
            spec(kk).iCV2 = iCV;

            spec(kk).Theta1 = ThetaVec(iTh);
            spec(kk).CV2_1 = CV2Vec(iCV);

            spec(kk).seed_cal  = seedCellBase+1000+kk;
            spec(kk).seed_eval = seedCellBase+2000+kk;
            spec(kk).seed_boot = seedBootBase+kk;
        end
    end
end

function pc = make_cell_parameters(p,sim,Theta1,CV2_1)
    pc = p;

    pc.Theta1 = Theta1;
    pc.CV2_1  = CV2_1;

    pc.kappa1 = 1/(Theta1*sim.T_total);
    pc.shape1 = 1/CV2_1;
    pc.theta1 = pc.mu_c1/pc.shape1;
    pc.a1     = pc.shape1*pc.kappa1;

    pc.var_c1 = pc.mu_c1^2*CV2_1;
    pc.tau1   = 1/pc.kappa1;

    pc.expected_GammaOU_jumps_60min = ...
        pc.a1*sim.T_total;
end

function R = run_one_regime_cell( ...
        p,sim,nCal,nEval,nBoot,spec,M_abs,meanTolerance)

    pc = make_cell_parameters( ...
        p,sim,spec.Theta1,spec.CV2_1);

    % Independent calibration.
    rng(spec.seed_cal,'twister');

    CAL = simulate_calibration(pc,sim,nCal);

    Mb_std = max(1,round(quantile( ...
        [CAL.CIR.Y60;CAL.GOU.Y60], ...
        sim.q_threshold)));

    % Independent paired evaluation.
    rng(spec.seed_eval,'twister');

    EVAL = simulate_evaluation( ...
        pc,sim,nEval,Mb_std,M_abs);

    % Full-state output summaries.
    meanC = mean(EVAL.CIR.Y60);
    meanG = mean(EVAL.GOU.Y60);

    varC = var(EVAL.CIR.Y60,0);
    varG = var(EVAL.GOU.Y60,0);

    fanoC = varC/max(meanC,eps);
    fanoG = varG/max(meanG,eps);

    signedMean = ...
        (meanG-meanC)/max(0.5*(meanC+meanG),eps);

    signedFano = ...
        (fanoG-fanoC)/max(0.5*(fanoC+fanoG),eps);

    stdStats = fpt_effect_summary( ...
        EVAL.CIR.FPT_std,EVAL.GOU.FPT_std,sim.T_total);

    absStats = fpt_effect_summary( ...
        EVAL.CIR.FPT_abs,EVAL.GOU.FPT_abs,sim.T_total);

    % Paired restricted-time contributions for the primary RMST test.
    timeC = restricted_time(EVAL.CIR.FPT_std,sim.T_total);
    timeG = restricted_time(EVAL.GOU.FPT_std,sim.T_total);
    pairedDifference = timeC-timeG;

    p_RMST = paired_mean_pvalue(pairedDifference);

    % Paired bootstrap.
    rng(spec.seed_boot,'twister');

    BOOT = bootstrap_cell_effects( ...
        EVAL,sim,nBoot);

    CI_DeltaF_std = ...
        prctile(BOOT.DeltaF_std,[2.5 97.5]);

    CI_DeltaRMST_std_min = ...
        prctile(BOOT.DeltaRMST_std/60,[2.5 97.5]);

    CI_DeltaF_abs = ...
        prctile(BOOT.DeltaF_abs,[2.5 97.5]);

    CI_DeltaRMST_abs_min = ...
        prctile(BOOT.DeltaRMST_abs/60,[2.5 97.5]);

    SE_DeltaRMST_std_min = ...
        std(BOOT.DeltaRMST_std/60,0);

    R = struct();

    R.linear_index = spec.linear_index;
    R.iTheta = spec.iTheta;
    R.iCV2 = spec.iCV2;

    R.Theta1 = pc.Theta1;
    R.CV2_1 = pc.CV2_1;

    R.kappa1 = pc.kappa1;
    R.shape1 = pc.shape1;
    R.theta1 = pc.theta1;
    R.a1 = pc.a1;

    R.var_c1 = pc.var_c1;
    R.tau1_s = pc.tau1;
    R.expected_GammaOU_jumps_60min = ...
        pc.expected_GammaOU_jumps_60min;

    R.Mb_std = Mb_std;
    R.M_abs = M_abs;

    R.meanY60_CIR = meanC;
    R.meanY60_GOU = meanG;

    R.varY60_CIR = varC;
    R.varY60_GOU = varG;

    R.fano60_CIR = fanoC;
    R.fano60_GOU = fanoG;

    R.RelativeMeanSigned = signedMean;
    R.RelativeFanoSigned = signedFano;
    R.mean_audit_pass = abs(signedMean)<=meanTolerance;

    R.F60_std_CIR = stdStats.Ftau_CIR;
    R.F60_std_GOU = stdStats.Ftau_GOU;
    R.DeltaF_std = stdStats.DeltaF;
    R.DeltaRMST_std_min = stdStats.DeltaRMST/60;
    R.CI_DeltaF_std = CI_DeltaF_std;
    R.CI_DeltaRMST_std_min = CI_DeltaRMST_std_min;
    R.SE_DeltaRMST_std_min = SE_DeltaRMST_std_min;
    R.p_RMST = p_RMST;

    R.F60_abs_CIR = absStats.Ftau_CIR;
    R.F60_abs_GOU = absStats.Ftau_GOU;
    R.eventRate_abs_pooled = ...
        0.5*(absStats.Ftau_CIR+absStats.Ftau_GOU);

    R.DeltaF_abs = absStats.DeltaF;
    R.DeltaRMST_abs_min = absStats.DeltaRMST/60;
    R.CI_DeltaF_abs = CI_DeltaF_abs;
    R.CI_DeltaRMST_abs_min = CI_DeltaRMST_abs_min;
end

function OUT = simulate_calibration(p,sim,n)
    nSteps = round(sim.T_total/sim.dt);

    x1Initial = binornd(p.N,p.p_x1_ss,n,1);
    c1Initial = gamrnd(p.shape1,p.theta1,n,1);
    c2Initial = gamrnd(p.shape2,p.theta2,n,1);

    x1C = x1Initial;
    x2C = p.N-x1C;
    yC  = zeros(n,1);
    c1C = c1Initial;
    c2C = c2Initial;

    x1G = x1Initial;
    x2G = p.N-x1G;
    yG  = zeros(n,1);
    c1G = c1Initial;
    c2G = c2Initial;

    for it=1:nSteps
        % CIR.
        c1Old = c1C;
        c2Old = c2C;

        c1C = cir_exact_transition( ...
            c1C,p.kappa1,p.shape1,p.theta1,sim.dt);

        c2C = cir_exact_transition( ...
            c2C,p.kappa2,p.shape2,p.theta2,sim.dt);

        I1C = 0.5*(c1Old+c1C)*sim.dt;
        I2C = 0.5*(c2Old+c2C)*sim.dt;

        pr1C = min(max(1-exp(-I1C),0),1);
        pr2C = min(max(1-exp(-I2C),0),1);

        r1C = binornd(x1C,pr1C);
        r2C = binornd(x2C,pr2C);

        x1C = x1C-r1C+r2C;
        x2C = p.N-x1C;
        yC  = yC+r1C;

        % GammaOU.
        [c1G,I1G] = gamma_ou_step_with_integral( ...
            c1G,p.a1,p.kappa1,p.theta1,sim.dt);

        [c2G,I2G] = gamma_ou_step_with_integral( ...
            c2G,p.a2,p.kappa2,p.theta2,sim.dt);

        pr1G = min(max(1-exp(-I1G),0),1);
        pr2G = min(max(1-exp(-I2G),0),1);

        r1G = binornd(x1G,pr1G);
        r2G = binornd(x2G,pr2G);

        x1G = x1G-r1G+r2G;
        x2G = p.N-x1G;
        yG  = yG+r1G;
    end

    OUT.CIR.Y60 = yC;
    OUT.GOU.Y60 = yG;
end

function OUT = simulate_evaluation( ...
        p,sim,n,Mb_std,M_abs)

    nSteps = round(sim.T_total/sim.dt);

    x1Initial = binornd(p.N,p.p_x1_ss,n,1);
    c1Initial = gamrnd(p.shape1,p.theta1,n,1);
    c2Initial = gamrnd(p.shape2,p.theta2,n,1);

    x1C = x1Initial;
    x2C = p.N-x1C;
    yC  = zeros(n,1);
    c1C = c1Initial;
    c2C = c2Initial;

    x1G = x1Initial;
    x2G = p.N-x1G;
    yG  = zeros(n,1);
    c1G = c1Initial;
    c2G = c2Initial;

    FPTstdC = nan(n,1);
    FPTstdG = nan(n,1);

    FPTabsC = nan(n,1);
    FPTabsG = nan(n,1);

    for it=1:nSteps
        tNow = it*sim.dt;

        % CIR.
        c1Old = c1C;
        c2Old = c2C;

        c1C = cir_exact_transition( ...
            c1C,p.kappa1,p.shape1,p.theta1,sim.dt);

        c2C = cir_exact_transition( ...
            c2C,p.kappa2,p.shape2,p.theta2,sim.dt);

        I1C = 0.5*(c1Old+c1C)*sim.dt;
        I2C = 0.5*(c2Old+c2C)*sim.dt;

        pr1C = min(max(1-exp(-I1C),0),1);
        pr2C = min(max(1-exp(-I2C),0),1);

        r1C = binornd(x1C,pr1C);
        r2C = binornd(x2C,pr2C);

        x1C = x1C-r1C+r2C;
        x2C = p.N-x1C;
        yC  = yC+r1C;

        FPTstdC = update_fpt(FPTstdC,yC,Mb_std,tNow);
        FPTabsC = update_fpt(FPTabsC,yC,M_abs,tNow);

        % GammaOU.
        [c1G,I1G] = gamma_ou_step_with_integral( ...
            c1G,p.a1,p.kappa1,p.theta1,sim.dt);

        [c2G,I2G] = gamma_ou_step_with_integral( ...
            c2G,p.a2,p.kappa2,p.theta2,sim.dt);

        pr1G = min(max(1-exp(-I1G),0),1);
        pr2G = min(max(1-exp(-I2G),0),1);

        r1G = binornd(x1G,pr1G);
        r2G = binornd(x2G,pr2G);

        x1G = x1G-r1G+r2G;
        x2G = p.N-x1G;
        yG  = yG+r1G;

        FPTstdG = update_fpt(FPTstdG,yG,Mb_std,tNow);
        FPTabsG = update_fpt(FPTabsG,yG,M_abs,tNow);
    end

    OUT.CIR.Y60 = yC;
    OUT.CIR.FPT_std = FPTstdC;
    OUT.CIR.FPT_abs = FPTabsC;

    OUT.GOU.Y60 = yG;
    OUT.GOU.FPT_std = FPTstdG;
    OUT.GOU.FPT_abs = FPTabsG;
end

function fpt = update_fpt( ...
        fpt,currentValue,threshold,currentTime)

    crossing = isnan(fpt) & currentValue>=threshold;
    fpt(crossing) = currentTime;
end

function tRestricted = restricted_time(fpt,tau)
    tRestricted = fpt;
    tRestricted(isnan(tRestricted)) = tau;
    tRestricted = min(tRestricted,tau);
end

function S = fpt_effect_summary(fptC,fptG,tau)
    eventC = ~isnan(fptC) & fptC<=tau;
    eventG = ~isnan(fptG) & fptG<=tau;

    timeC = restricted_time(fptC,tau);
    timeG = restricted_time(fptG,tau);

    S.Ftau_CIR = mean(eventC);
    S.Ftau_GOU = mean(eventG);

    S.DeltaF = S.Ftau_GOU-S.Ftau_CIR;

    S.RMST_CIR = mean(timeC);
    S.RMST_GOU = mean(timeG);

    S.DeltaRMST = S.RMST_CIR-S.RMST_GOU;
end

function pValue = paired_mean_pvalue(differenceVector)
    differenceVector = differenceVector(isfinite(differenceVector));
    n = numel(differenceVector);

    if n<2
        pValue = nan;
        return;
    end

    sd = std(differenceVector,0);

    if sd==0
        if mean(differenceVector)==0
            pValue = 1;
        else
            pValue = 0;
        end
        return;
    end

    tStatistic = mean(differenceVector)/(sd/sqrt(n));
    pValue = 2*tcdf(-abs(tStatistic),n-1);
end

function B = bootstrap_cell_effects(EVAL,sim,nBoot)
    n = numel(EVAL.CIR.Y60);

    B.DeltaF_std = nan(nBoot,1);
    B.DeltaRMST_std = nan(nBoot,1);

    B.DeltaF_abs = nan(nBoot,1);
    B.DeltaRMST_abs = nan(nBoot,1);

    for ib=1:nBoot
        idx = randi(n,n,1);

        stdStats = fpt_effect_summary( ...
            EVAL.CIR.FPT_std(idx), ...
            EVAL.GOU.FPT_std(idx), ...
            sim.T_total);

        absStats = fpt_effect_summary( ...
            EVAL.CIR.FPT_abs(idx), ...
            EVAL.GOU.FPT_abs(idx), ...
            sim.T_total);

        B.DeltaF_std(ib) = stdStats.DeltaF;
        B.DeltaRMST_std(ib) = stdStats.DeltaRMST;

        B.DeltaF_abs(ib) = absStats.DeltaF;
        B.DeltaRMST_abs(ib) = absStats.DeltaRMST;
    end
end

function cNext = cir_exact_transition( ...
        c,kappa,shape,theta,dt)

    if dt==0
        cNext=c;
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
        cNew=cOld;
        I=zeros(size(cOld));
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

function M = assemble_regular_matrices( ...
        results,nCV,nTheta)

    fieldNames = { ...
        'DeltaRMST','DeltaF','DeltaRMST_abs','DeltaF_abs', ...
        'signedMeanPct','signedFanoPct','Mb', ...
        'CIlo_RMST','CIhi_RMST','CIhalf_RMST', ...
        'CIlo_F','CIhi_F','p_RMST', ...
        'mean_audit_pass','eventRate_abs_pooled', ...
        'F60_std_CIR','F60_std_GOU'};

    for k=1:numel(fieldNames)
        if strcmp(fieldNames{k},'mean_audit_pass')
            M.(fieldNames{k}) = false(nCV,nTheta);
        else
            M.(fieldNames{k}) = nan(nCV,nTheta);
        end
    end

    for i=1:numel(results)
        r = results{i};

        iCV = r.iCV2;
        iTh = r.iTheta;

        M.DeltaRMST(iCV,iTh) = r.DeltaRMST_std_min;
        M.DeltaF(iCV,iTh) = r.DeltaF_std;

        M.DeltaRMST_abs(iCV,iTh) = r.DeltaRMST_abs_min;
        M.DeltaF_abs(iCV,iTh) = r.DeltaF_abs;

        M.signedMeanPct(iCV,iTh) = ...
            100*r.RelativeMeanSigned;

        M.signedFanoPct(iCV,iTh) = ...
            100*r.RelativeFanoSigned;

        M.Mb(iCV,iTh) = r.Mb_std;

        M.CIlo_RMST(iCV,iTh) = ...
            r.CI_DeltaRMST_std_min(1);

        M.CIhi_RMST(iCV,iTh) = ...
            r.CI_DeltaRMST_std_min(2);

        M.CIhalf_RMST(iCV,iTh) = ...
            0.5*(r.CI_DeltaRMST_std_min(2)- ...
                 r.CI_DeltaRMST_std_min(1));

        M.CIlo_F(iCV,iTh) = r.CI_DeltaF_std(1);
        M.CIhi_F(iCV,iTh) = r.CI_DeltaF_std(2);

        M.p_RMST(iCV,iTh) = r.p_RMST;
        M.mean_audit_pass(iCV,iTh) = r.mean_audit_pass;

        M.eventRate_abs_pooled(iCV,iTh) = ...
            r.eventRate_abs_pooled;

        M.F60_std_CIR(iCV,iTh) = r.F60_std_CIR;
        M.F60_std_GOU(iCV,iTh) = r.F60_std_GOU;
    end
end

function [qMatrix,rejectMatrix] = ...
        benjamini_hochberg(pMatrix,alpha)

    p = pMatrix(:);
    valid = isfinite(p);

    pValid = p(valid);
    m = numel(pValid);

    qValid = nan(size(pValid));
    rejectValid = false(size(pValid));

    if m>0
        [pSort,order] = sort(pValid);
        rankVector = (1:m)';

        qSort = pSort.*m./rankVector;

        for k=m-1:-1:1
            qSort(k) = min(qSort(k),qSort(k+1));
        end

        qSort = min(qSort,1);

        inverseOrder = zeros(m,1);
        inverseOrder(order) = 1:m;

        qValid = qSort(inverseOrder);
        rejectValid = qValid<=alpha;
    end

    q = nan(size(p));
    reject = false(size(p));

    q(valid) = qValid;
    reject(valid) = rejectValid;

    qMatrix = reshape(q,size(pMatrix));
    rejectMatrix = reshape(reject,size(pMatrix));
end

function T = results_to_table(results)
    n = numel(results);

    linear_index = zeros(n,1);
    iTheta = zeros(n,1);
    iCV2 = zeros(n,1);

    Theta1 = zeros(n,1);
    CV1_squared = zeros(n,1);

    kappa1 = zeros(n,1);
    shape1 = zeros(n,1);
    theta1 = zeros(n,1);
    a1 = zeros(n,1);
    expectedJumps = zeros(n,1);

    Mb = zeros(n,1);

    meanC = zeros(n,1);
    meanG = zeros(n,1);
    fanoC = zeros(n,1);
    fanoG = zeros(n,1);

    signedMean = zeros(n,1);
    signedFano = zeros(n,1);

    dF = zeros(n,1);
    dFLo = zeros(n,1);
    dFHi = zeros(n,1);

    dRMST = zeros(n,1);
    dRMSTLo = zeros(n,1);
    dRMSTHi = zeros(n,1);
    seRMST = zeros(n,1);

    pRMST = zeros(n,1);
    meanPass = false(n,1);

    dFabs = zeros(n,1);
    dRMSTabs = zeros(n,1);
    eventRateAbs = zeros(n,1);

    for i=1:n
        r = results{i};

        linear_index(i) = r.linear_index;
        iTheta(i) = r.iTheta;
        iCV2(i) = r.iCV2;

        Theta1(i) = r.Theta1;
        CV1_squared(i) = r.CV2_1;

        kappa1(i) = r.kappa1;
        shape1(i) = r.shape1;
        theta1(i) = r.theta1;
        a1(i) = r.a1;
        expectedJumps(i) = ...
            r.expected_GammaOU_jumps_60min;

        Mb(i) = r.Mb_std;

        meanC(i) = r.meanY60_CIR;
        meanG(i) = r.meanY60_GOU;
        fanoC(i) = r.fano60_CIR;
        fanoG(i) = r.fano60_GOU;

        signedMean(i) = 100*r.RelativeMeanSigned;
        signedFano(i) = 100*r.RelativeFanoSigned;

        dF(i) = r.DeltaF_std;
        dFLo(i) = r.CI_DeltaF_std(1);
        dFHi(i) = r.CI_DeltaF_std(2);

        dRMST(i) = r.DeltaRMST_std_min;
        dRMSTLo(i) = r.CI_DeltaRMST_std_min(1);
        dRMSTHi(i) = r.CI_DeltaRMST_std_min(2);
        seRMST(i) = r.SE_DeltaRMST_std_min;

        pRMST(i) = r.p_RMST;
        meanPass(i) = r.mean_audit_pass;

        dFabs(i) = r.DeltaF_abs;
        dRMSTabs(i) = r.DeltaRMST_abs_min;
        eventRateAbs(i) = r.eventRate_abs_pooled;
    end

    T = table( ...
        linear_index,iTheta,iCV2, ...
        Theta1,CV1_squared, ...
        kappa1,shape1,theta1,a1,expectedJumps, ...
        Mb, ...
        meanC,meanG,fanoC,fanoG, ...
        signedMean,signedFano, ...
        dF,dFLo,dFHi, ...
        dRMST,dRMSTLo,dRMSTHi,seRMST,pRMST,meanPass, ...
        dFabs,dRMSTabs,eventRateAbs, ...
        'VariableNames',{ ...
        'linear_index','iTheta','iCV2', ...
        'Theta1','CV1_squared', ...
        'kappa1_per_s','shape1','theta1_scale', ...
        'GammaOU_jump_rate_per_s', ...
        'expected_GammaOU_jumps_60min', ...
        'Mb_cell_specific', ...
        'meanY60_CIR','meanY60_GammaOU', ...
        'Fano60_CIR','Fano60_GammaOU', ...
        'signed_relative_mean_difference_percent', ...
        'signed_relative_Fano_difference_percent', ...
        'DeltaF60_standardized', ...
        'DeltaF60_standardized_CI_low', ...
        'DeltaF60_standardized_CI_high', ...
        'DeltaRMST60_standardized_min', ...
        'DeltaRMST60_standardized_CI_low_min', ...
        'DeltaRMST60_standardized_CI_high_min', ...
        'SE_DeltaRMST60_standardized_min', ...
        'paired_p_RMST', ...
        'mean_audit_pass', ...
        'DeltaF60_baseline_fixed_threshold', ...
        'DeltaRMST60_baseline_fixed_threshold_min', ...
        'pooled_event_rate_baseline_fixed_threshold'});
end

function [theta,estimate,lo,hi] = extract_memory_slice( ...
        ThetaVec,CV2Vec,M,targetCV,refinementTable)

    [~,iCV] = min(abs(CV2Vec-targetCV));

    theta = ThetaVec(:);
    estimate = M.DeltaRMST(iCV,:)';
    lo = M.CIlo_RMST(iCV,:)';
    hi = M.CIhi_RMST(iCV,:)';

    if ~isempty(refinementTable)
        select = abs(refinementTable.CV1_squared-targetCV)<1e-10;

        theta = [theta;refinementTable.Theta1(select)];
        estimate = [estimate; ...
            refinementTable.DeltaRMST60_standardized_min(select)];

        lo = [lo; ...
            refinementTable.DeltaRMST60_standardized_CI_low_min(select)];

        hi = [hi; ...
            refinementTable.DeltaRMST60_standardized_CI_high_min(select)];

        [theta,order] = sort(theta);
        estimate = estimate(order);
        lo = lo(order);
        hi = hi(order);
    end
end

function [cv,estimate,lo,hi] = extract_dispersion_slice( ...
        ThetaVec,CV2Vec,M,targetTheta)

    [~,iTheta] = min(abs(ThetaVec-targetTheta));

    cv = CV2Vec(:);
    estimate = M.DeltaRMST(:,iTheta);
    lo = M.CIlo_RMST(:,iTheta);
    hi = M.CIhi_RMST(:,iTheta);
end

function xScaled = rescale_to_interval(x,newMin,newMax)
    if max(x)==min(x)
        xScaled = repmat(0.5*(newMin+newMax),size(x));
        return;
    end

    xScaled = newMin + ...
        (x-min(x))./(max(x)-min(x)).*(newMax-newMin);
end

function plot_refined_memory_slices( ...
        ax,cfg,M,refinementTable,baseline,style)

    hold(ax,'on');

    targetCV = [baseline.CV2_1,8];

    lineColors = [style.cir;style.gou];
    markerList = {'o','s'};

    for k=1:numel(targetCV)
        [theta,estimate,lo,hi] = ...
            extract_memory_slice( ...
            cfg.Theta1_vec,cfg.CV2_vec,M, ...
            targetCV(k),refinementTable);

        errorbar(ax,theta,estimate, ...
            estimate-lo,hi-estimate, ...
            [markerList{k} '-'], ...
            'Color',lineColors(k,:), ...
            'MarkerFaceColor','w', ...
            'MarkerSize',style.marker_size, ...
            'LineWidth',style.data_line, ...
            'CapSize',3, ...
            'DisplayName',sprintf('CV_1^2=%.0f',targetCV(k)));
    end
end

function K = run_keypoint_replication( ...
        p,sim,cfg,keyPoints,M_abs)

    nPoints = numel(keyPoints);
    nSeeds = cfg.key_n_seeds;

    seedResults = cell(nPoints,nSeeds);

    for ip=1:nPoints
        for iseed=1:nSeeds
            spec = struct();

            spec.linear_index = iseed;
            spec.iTheta = 1;
            spec.iCV2 = ip;

            spec.Theta1 = keyPoints(ip).Theta1;
            spec.CV2_1 = keyPoints(ip).CV2_1;

            seedOffset = 1000*ip+10*iseed;

            spec.seed_cal = ...
                cfg.seed_keypoint_base+seedOffset;

            spec.seed_eval = ...
                cfg.seed_keypoint_base+10000+seedOffset;

            spec.seed_boot = ...
                cfg.seed_keypoint_base+20000+seedOffset;

            fprintf(['  key point %d/%d, seed %d/%d: ' ...
                '%s, Theta=%.4g, CV^2=%.4g\n'], ...
                ip,nPoints,iseed,nSeeds, ...
                keyPoints(ip).label, ...
                keyPoints(ip).Theta1, ...
                keyPoints(ip).CV2_1);

            seedResults{ip,iseed} = run_one_regime_cell( ...
                p,sim, ...
                cfg.key_n_cal,cfg.key_n_eval,cfg.key_n_boot, ...
                spec,M_abs,cfg.mean_audit_tolerance);
        end
    end

    label = strings(nPoints*nSeeds,1);
    pointIndex = zeros(nPoints*nSeeds,1);
    seedIndex = zeros(nPoints*nSeeds,1);

    Theta1 = zeros(nPoints*nSeeds,1);
    CV2 = zeros(nPoints*nSeeds,1);

    estimate = zeros(nPoints*nSeeds,1);
    ciLow = zeros(nPoints*nSeeds,1);
    ciHigh = zeros(nPoints*nSeeds,1);
    se = zeros(nPoints*nSeeds,1);

    row = 0;

    for ip=1:nPoints
        for iseed=1:nSeeds
            row = row+1;
            r = seedResults{ip,iseed};

            label(row) = string(keyPoints(ip).label);
            pointIndex(row) = ip;
            seedIndex(row) = iseed;

            Theta1(row) = r.Theta1;
            CV2(row) = r.CV2_1;

            estimate(row) = r.DeltaRMST_std_min;
            ciLow(row) = r.CI_DeltaRMST_std_min(1);
            ciHigh(row) = r.CI_DeltaRMST_std_min(2);
            se(row) = max(r.SE_DeltaRMST_std_min,1e-8);
        end
    end

    seedTable = table( ...
        label,pointIndex,seedIndex,Theta1,CV2, ...
        estimate,ciLow,ciHigh,se, ...
        'VariableNames',{ ...
        'label','point_index','seed_index','Theta1','CV1_squared', ...
        'DeltaRMST60_min','CI_low_min','CI_high_min','SE_min'});

    summaryLabel = strings(nPoints,1);
    summaryTheta = zeros(nPoints,1);
    summaryCV2 = zeros(nPoints,1);

    pooledEstimate = zeros(nPoints,1);
    pooledSE = zeros(nPoints,1);
    pooledLow = zeros(nPoints,1);
    pooledHigh = zeros(nPoints,1);

    heterogeneityQ = zeros(nPoints,1);

    for ip=1:nPoints
        select = pointIndex==ip;

        est = estimate(select);
        seLocal = se(select);

        weights = 1./(seLocal.^2);

        pooledEstimate(ip) = ...
            sum(weights.*est)/sum(weights);

        pooledSE(ip) = sqrt(1/sum(weights));

        pooledLow(ip) = pooledEstimate(ip)-1.96*pooledSE(ip);
        pooledHigh(ip) = pooledEstimate(ip)+1.96*pooledSE(ip);

        heterogeneityQ(ip) = ...
            sum(weights.*(est-pooledEstimate(ip)).^2);

        summaryLabel(ip) = string(keyPoints(ip).label);
        summaryTheta(ip) = keyPoints(ip).Theta1;
        summaryCV2(ip) = keyPoints(ip).CV2_1;
    end

    summaryTable = table( ...
        summaryLabel,summaryTheta,summaryCV2, ...
        pooledEstimate,pooledSE,pooledLow,pooledHigh,heterogeneityQ, ...
        'VariableNames',{ ...
        'label','Theta1','CV1_squared', ...
        'pooled_DeltaRMST60_min','pooled_SE_min', ...
        'pooled_CI_low_min','pooled_CI_high_min', ...
        'fixed_effect_heterogeneity_Q'});

    K = struct();
    K.key_points = keyPoints;
    K.seed_results = seedResults;
    K.seed_table = seedTable;
    K.summary_table = summaryTable;
end

function plot_keypoint_replication(ax,K,style)
    summary = K.summary_table;
    seedTable = K.seed_table;

    nPoints = height(summary);
    y = (1:nPoints)';

    hold(ax,'on');

    for ip=1:nPoints
        % Fixed-effect pooled 95% interval.
        plot(ax, ...
            [summary.pooled_CI_low_min(ip), ...
             summary.pooled_CI_high_min(ip)], ...
            [ip ip], ...
            '-','Color',style.theory, ...
            'LineWidth',style.data_line);

        plot(ax, ...
            summary.pooled_CI_low_min(ip),ip, ...
            '|','Color',style.theory, ...
            'MarkerSize',6, ...
            'LineWidth',style.data_line);

        plot(ax, ...
            summary.pooled_CI_high_min(ip),ip, ...
            '|','Color',style.theory, ...
            'MarkerSize',6, ...
            'LineWidth',style.data_line);

        plot(ax, ...
            summary.pooled_DeltaRMST60_min(ip),ip, ...
            'o','Color',style.theory, ...
            'MarkerFaceColor','w', ...
            'MarkerSize',style.marker_size+1, ...
            'LineWidth',style.data_line);

        % Individual independent-seed estimates.
        select = seedTable.point_index==ip;
        nSeed = sum(select);

        scatter(ax, ...
            seedTable.DeltaRMST60_min(select), ...
            ip+0.12*linspace(-1,1,nSeed)', ...
            18,style.gou,'filled', ...
            'MarkerFaceAlpha',0.65);
    end

    xline(ax,0,':', ...
        'Color',style.light_gray, ...
        'LineWidth',style.reference_line);

    set(ax,'YTick',y, ...
        'YTickLabel',cellstr(summary.label));

    ylim(ax,[0.4 nPoints+0.6]);
    set(ax,'YDir','reverse');

    xlabel(ax,'\DeltaRMST_{60} (min)');
    ylabel(ax,'');

    apply_journal_axes(ax,style);
end

function plot_discrete_regime_map( ...
        ax,ThetaVec,CV2Vec,Z,colorLimits,cmap, ...
        colorbarLabel,panelTitle,baseline,style)

    hold(ax,'on');

    colorLimits = sanitize_color_limits(colorLimits);

    xCenter = log10(ThetaVec);
    yCenter = log10(CV2Vec);

    xEdge = centers_to_edges(xCenter);
    yEdge = centers_to_edges(yCenter);

    for i=1:numel(CV2Vec)
        for j=1:numel(ThetaVec)
            xPatch = [ ...
                xEdge(j),xEdge(j+1),xEdge(j+1),xEdge(j)];

            yPatch = [ ...
                yEdge(i),yEdge(i),yEdge(i+1),yEdge(i+1)];

            patch(ax, ...
                'XData',xPatch, ...
                'YData',yPatch, ...
                'FaceVertexCData',Z(i,j), ...
                'FaceColor','flat', ...
                'EdgeColor',[1 1 1], ...
                'LineWidth',0.45);
        end
    end

    colormap(ax,cmap);
    caxis(ax,colorLimits);

    cb = colorbar(ax);
    cb.Label.String = colorbarLabel;
    cb.Label.FontName = style.font_name;
    cb.Label.FontSize = style.label_font;

    set(ax, ...
        'XTick',xCenter, ...
        'XTickLabel',compose('%.2g',ThetaVec), ...
        'YTick',yCenter, ...
        'YTickLabel',compose('%.2g',CV2Vec));

    xlabel(ax,'Memory ratio, \Theta_1=1/(\kappa_1T_{obs})');
    ylabel(ax,'Stationary dispersion, CV_1^2=1/shape_1');

    title(ax,panelTitle, ...
        'FontWeight','normal', ...
        'FontSize',style.title_font);

    xlim(ax,[xEdge(1) xEdge(end)]);
    ylim(ax,[yEdge(1) yEdge(end)]);

    apply_journal_axes(ax,style);
end

function overlay_regime_markers( ...
        ax,ThetaVec,CV2Vec,positiveMask,failMask,baseline,style)

    x = log10(ThetaVec);
    y = log10(CV2Vec);
    [X,Y] = meshgrid(x,y);

    if any(positiveMask(:))
        plot(ax,X(positiveMask),Y(positiveMask), ...
            'o','Color',[0.05 0.05 0.05], ...
            'MarkerFaceColor','w', ...
            'MarkerSize',4.6, ...
            'LineWidth',0.85);
    end

    if any(failMask(:))
        plot(ax,X(failMask),Y(failMask), ...
            'x','Color',[0.05 0.05 0.05], ...
            'MarkerSize',5.5, ...
            'LineWidth',1.0);
    end

    plot(ax,log10(baseline.Theta1),log10(baseline.CV2_1), ...
        'p','Color',[0.05 0.05 0.05], ...
        'MarkerFaceColor',[1.00 0.82 0.18], ...
        'MarkerSize',10, ...
        'LineWidth',1.1);
end

function edges = centers_to_edges(centers)
    centers = centers(:)';

    if numel(centers)==1
        edges = [centers-0.5 centers+0.5];
        return;
    end

    midpoint = 0.5*(centers(1:end-1)+centers(2:end));

    firstEdge = centers(1)-0.5*(centers(2)-centers(1));
    lastEdge = centers(end)+0.5*(centers(end)-centers(end-1));

    edges = [firstEdge midpoint lastEdge];
end

function limits = sanitize_color_limits(limits)
    limits = double(limits(:)');

    if numel(limits)~=2 || any(~isfinite(limits))
        limits = [0 1];
        return;
    end

    if limits(2)<=limits(1)
        center = limits(1);

        if center==0
            limits = [-0.5 0.5];
        else
            pad = max(abs(center)*0.05,1e-6);
            limits = [center-pad center+pad];
        end
    end
end

function lim = symmetric_color_limit(Z,minLimit)
    values = abs(Z(isfinite(Z)));

    if isempty(values)
        lim = minLimit;
    else
        lim = max(max(values),minLimit);
    end
end

function style = journal_style()
    style.font_name = 'Arial';

    style.main_width_cm  = 18.3;
    style.main_height_cm = 7.8;

    style.si_width_cm  = 18.3;
    style.si_height_cm = 12.5;

    style.axis_font       = 7.0;
    style.label_font      = 7.5;
    style.title_font      = 8.0;
    style.legend_font     = 6.4;
    style.annotation_font = 5.9;
    style.small_font      = 5.6;
    style.panel_font      = 10.5;

    style.axis_line      = 0.75;
    style.data_line      = 1.20;
    style.reference_line = 0.75;
    style.marker_size    = 4.5;

    style.cir = hex2rgb('#3C5488');
    style.gou = hex2rgb('#E64B35');

    style.theory     = [0.12 0.12 0.12];
    style.dark_gray  = [0.27 0.27 0.27];
    style.light_gray = [0.72 0.72 0.72];

    style.diverging_map = blue_white_red(256);
    style.sequential_map = parula(256);
end

function apply_journal_axes(ax,style)
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
            ax.Toolbar.Visible='off';
        end
    catch
    end

    try
        disableDefaultInteractivity(ax);
    catch
    end
end

function add_panel_label(ax,labelText,style)
    text(ax,-0.16,1.10,labelText, ...
        'Units','normalized', ...
        'FontName',style.font_name, ...
        'FontSize',style.panel_font, ...
        'FontWeight','bold', ...
        'Color',[0 0 0], ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','top', ...
        'Clipping','off');
end

function export_publication_figure( ...
        fig,stem,dpi,saveMatlabFigure)

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

    if saveMatlabFigure
        savefig(fig,[stem '.fig']);
    end
end

function map = blue_white_red(n)
    if nargin<1
        n=256;
    end

    blue  = [0.20 0.35 0.65];
    white = [1.00 1.00 1.00];
    red   = [0.78 0.20 0.16];

    n1 = floor(n/2);
    n2 = n-n1;

    left = [ ...
        linspace(blue(1),white(1),n1)', ...
        linspace(blue(2),white(2),n1)', ...
        linspace(blue(3),white(3),n1)'];

    right = [ ...
        linspace(white(1),red(1),n2)', ...
        linspace(white(2),red(2),n2)', ...
        linspace(white(3),red(3),n2)'];

    map = [left;right];
end

function rgb = hex2rgb(hexString)
    hexString = strrep(hexString,'#','');
    rgb = sscanf(hexString,'%2x%2x%2x',[1 3])/255;
end
