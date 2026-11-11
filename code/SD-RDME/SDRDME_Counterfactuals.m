%% SDRDME_Counterfactuals.m
% -------------------------------------------------------------------------
% Standalone U/H/B/F counterfactual decomposition for mitochondrial SD-RDME.
%
% U: uniform reaction + uniform diffusion
% H: heterogeneous reaction c1(r) + uniform diffusion
% B: uniform reaction + heterogeneous diffusion D(r)
% F: heterogeneous reaction c1(r) + heterogeneous diffusion D(r)
%
% Design:
%   1) Geometry, reaction parameters, diffusion parameters, initial iron
%      state, observation window, and exact next-subvolume SSA are identical
%      to the finalized predictor-independent Fig. 3B-C analysis.
%   2) U/H/B are matched to F by the exact closed first-moment system at
%      Tobs = 60 min. Matching scales c1 only and preserves each model's
%      spatial pattern.
%   3) The local region A_0.85 is fixed a priori as the top 15% of the
%      window-mean c1-hat score (radius 0.10), exactly matching the
%      candidate-screening and held-out-validation definition.
%   4) Main Fig. 3E shows the counterfactual phenotype decomposition.
%   5) SI shows event-resolved Fe2+ supply, loss, reduction, Fenton
%      consumption, inventory, and exact regional mass-balance closure.
%
% Outputs:
%   Fig3E_UHBF_counterfactual_main.{pdf,png,tif,fig}
%   FigS_UHBF_detailedMechanisms.{pdf,png,tif,fig}
%   Fig3E_UHBF_model_summary.csv
%   Fig3E_UHBF_contrast_summary.csv
%   Fig3E_UHBF_counterfactual_data.mat
%
% v2 corrections:
%   - A_0.85 now uses the locked radius-0.10 window-mean c1-hat score.
%   - The burden interaction is dimensionless and matches the ratio axis.
% -------------------------------------------------------------------------

clear; close all; clc;
clear functions;

%% --------------------------- Run configuration --------------------------
RUN_MODE = 'final';          % 'smoke' | 'pilot' | 'final'
USE_PARALLEL = true;
AUTO_START_POOL = true;
N_WORKERS = 8;
SAVE_FIG = true;
SAVE_DATA = true;

switch lower(RUN_MODE)
    case 'smoke'
        nRepBurden = 20;
        nRepFPT    = 30;
        nBoot      = 200;
    case 'pilot'
        nRepBurden = 100;
        nRepFPT    = 150;
        nBoot      = 500;
    case 'final'
        nRepBurden = 300;
        nRepFPT    = 300;
        nBoot      = 1000;
    otherwise
        error('Unknown RUN_MODE: %s', RUN_MODE);
end

seedBase = 20260729;

outMain     = 'Fig3E_UHBF_counterfactual_main';
outSI       = 'FigS_UHBF_detailedMechanisms';
outModelCSV = 'Fig3E_UHBF_model_summary.csv';
outContrCSV = 'Fig3E_UHBF_contrast_summary.csv';
outData     = 'Fig3E_UHBF_counterfactual_data.mat';

%% ------------------------------ Style -----------------------------------
set(groot,'defaultFigureColor','w');
set(groot,'defaultAxesFontName','Arial');
set(groot,'defaultTextFontName','Arial');
set(groot,'defaultAxesFontSize',7);
set(groot,'defaultTextFontSize',7);
set(groot,'defaultAxesLineWidth',0.70);
set(groot,'defaultLineLineWidth',1.05);
set(groot,'defaultAxesTickDir','out');
set(groot,'defaultAxesBox','on');
set(groot,'defaultLegendInterpreter','tex');
set(groot,'defaultAxesTickLabelInterpreter','tex');

colU = [0.40 0.43 0.48];
colH = [0.85 0.33 0.10];
colB = [0.15 0.43 0.76];
colF = [0.49 0.23 0.67];

modelColors = [colU; colH; colB; colF];
modelNames  = {'U','H','B','F'};
modelLabels = {'Uniform','Hotspot only','Bottleneck only','Full'};

%% ---------------------- Locked Fig. 3 parameters ------------------------
% Grid and crista-inspired geometry
nx = 60;
ny = 36;
aOut = 1.12;  bOut = 0.62;
aIn  = 0.94;  bIn  = 0.49;
ellG = 0.085;
ellS = 0.042;
connMatrix = 1.00;
connLumen  = 0.82;
connNeck   = 0.18;

% SD-RDME reaction and diffusion fields, min-based units
D0_Fe2   = 3.0e-4;
D2_scale = 0.70;
D0_Fe3   = D2_scale * D0_Fe2;
c1_0     = 6.0e-3;
etaS     = 3.40;
c2_0     = 1.8e-2;

% Observation window and initial iron state
Tobs     = 60;
dtRecord = 0.50;
tvec     = 0:dtRecord:Tobs;
FeTot_per_voxel = 6;
fracFe2_0       = 0.65;  % actual rounded state is Fe2+/Fe3+ = 4/2

% Closed first-moment global matching
dtMoment = 0.05;
nBisect  = 14;

% Fixed reaction-field region and control-anchored threshold
qA = 0.85;
windowRadius = 0.10;
thresholdFactor = 1.00;
alphaCI = 0.05;

%% --------------------------- Parallel pool ------------------------------
usePar = configureParallel(USE_PARALLEL,AUTO_START_POOL,N_WORKERS);

%% -------------------------- Geometry and fields -------------------------
geom = buildCurrentCristaGeometry(nx,ny,aOut,bOut,aIn,bIn, ...
    ellG,ellS,connMatrix,connLumen,connNeck);

mask   = geom.OmegaR;
x      = geom.x;
y      = geom.y;
Xg     = geom.X;
Yg     = geom.Y;
S      = geom.S;
G      = geom.G;
K      = nnz(mask);
idxMap = buildMaskIndexMap(mask);

D1Full = D0_Fe2 * G;
D2Full = D0_Fe3 * G;
D1Full(~mask) = NaN;
D2Full(~mask) = NaN;

c1Full = c1_0 * (1 + etaS*S);
c2Full = c2_0 * ones(size(c1Full));
c1Full(~mask) = NaN;
c2Full(~mask) = NaN;

c1F = c1Full(mask);
c2F = c2Full(mask);
D1F = D1Full(mask);
D2F = D2Full(mask);

c1Mean = mean(c1F);
D1Mean = mean(D1F);
D2Mean = mean(D2F);

N0   = FeTot_per_voxel * ones(K,1);
x1_0 = round(fracFe2_0*N0);
x2_0 = N0-x1_0;

% A_0.85 is selected from the same window-mean c1-hat score used in the
% candidate-field screening and held-out validation, then frozen before
% any U/H/B/F outcome simulations are generated.
coordX = Xg(mask);
coordY = Yg(mask);
win = buildCircularWindowsSparse(coordX,coordY,windowRadius);

c1Center = normalizeVector(c1F);
c1Window = (win.W*c1Center)./win.size;
c1Window = normalizeVector(c1Window);

A = topFractionMask(c1Window,1-qA);
nA = nnz(A);

Amap = false(size(mask));
Amap(mask) = A;

fprintf('============================================================\n');
fprintf('U/H/B/F counterfactual decomposition (%s)\n',upper(RUN_MODE));
fprintf('Grid                          = %d x %d\n',nx,ny);
fprintf('Omega_R voxels K              = %d\n',K);
fprintf('|A_%.2f|/K                    = %.4f (%d voxels)\n',qA,nA/K,nA);
fprintf('Window radius                 = %.3f\n',windowRadius);
fprintf('Window-size range             = [%d, %d]\n',min(win.size),max(win.size));
fprintf('Median window size            = %.1f\n',median(win.size));
fprintf('Tobs / dtRecord               = %.1f / %.2f min\n',Tobs,dtRecord);
fprintf('Initial Fe2+/Fe3+ per voxel   = %d / %d\n',x1_0(1),x2_0(1));
fprintf('c1 range                      = [%.5g, %.5g] min^-1\n',min(c1F),max(c1F));
fprintf('D_Fe2 range                   = [%.5g, %.5g]\n',min(D1F),max(D1F));
fprintf('Burden trajectories/model     = %d\n',nRepBurden);
fprintf('FPT trajectories/model        = %d\n',nRepFPT);
fprintf('Bootstrap replicates          = %d\n',nBoot);
fprintf('============================================================\n');

%% -------------------------- Diffusion edges -----------------------------
dx = mean(diff(x));
dy = mean(diff(y));

% Heterogeneous diffusion, used in B and F.
[fromHet1,toHet1,qHet1] = buildDirectedEdgesHarmonic( ...
    mask,idxMap,D1Full,dx,dy);
[fromHet2,toHet2,qHet2] = buildDirectedEdgesHarmonic( ...
    mask,idxMap,D2Full,dx,dy);

% Uniform diffusion, used in U and H.
D1Uniform = NaN(size(D1Full));
D2Uniform = NaN(size(D2Full));
D1Uniform(mask) = D1Mean;
D2Uniform(mask) = D2Mean;

[fromUni1,toUni1,qUni1] = buildDirectedEdgesHarmonic( ...
    mask,idxMap,D1Uniform,dx,dy);
[fromUni2,toUni2,qUni2] = buildDirectedEdgesHarmonic( ...
    mask,idxMap,D2Uniform,dx,dy);

%% ------------------------ Define U/H/B/F models -------------------------
M = repmat(struct(),1,4);

% U: uniform reaction, uniform diffusion
M(1).name   = 'U';
M(1).label  = modelLabels{1};
M(1).color  = colU;
M(1).c1Base = c1Mean*ones(K,1);
M(1).c2     = c2_0*ones(K,1);
M(1).from1  = fromUni1; M(1).to1 = toUni1; M(1).q1 = qUni1;
M(1).from2  = fromUni2; M(1).to2 = toUni2; M(1).q2 = qUni2;

% H: heterogeneous reaction, uniform diffusion
M(2).name   = 'H';
M(2).label  = modelLabels{2};
M(2).color  = colH;
M(2).c1Base = c1F;
M(2).c2     = c2_0*ones(K,1);
M(2).from1  = fromUni1; M(2).to1 = toUni1; M(2).q1 = qUni1;
M(2).from2  = fromUni2; M(2).to2 = toUni2; M(2).q2 = qUni2;

% B: uniform reaction, heterogeneous diffusion
M(3).name   = 'B';
M(3).label  = modelLabels{3};
M(3).color  = colB;
M(3).c1Base = c1Mean*ones(K,1);
M(3).c2     = c2_0*ones(K,1);
M(3).from1  = fromHet1; M(3).to1 = toHet1; M(3).q1 = qHet1;
M(3).from2  = fromHet2; M(3).to2 = toHet2; M(3).q2 = qHet2;

% F: heterogeneous reaction, heterogeneous diffusion
M(4).name   = 'F';
M(4).label  = modelLabels{4};
M(4).color  = colF;
M(4).c1Base = c1F;
M(4).c2     = c2_0*ones(K,1);
M(4).from1  = fromHet1; M(4).to1 = toHet1; M(4).q1 = qHet1;
M(4).from2  = fromHet2; M(4).to2 = toHet2; M(4).q2 = qHet2;

%% ---------------- Exact closed-moment global matching -------------------
fprintf('\nSolving the closed first-moment global matching problem...\n');

M(4).alpha = 1.0;
M(4).moment = momentEulerAll( ...
    M(4).alpha*M(4).c1Base,M(4).c2, ...
    M(4).from1,M(4).to1,M(4).q1, ...
    M(4).from2,M(4).to2,M(4).q2, ...
    x1_0,x2_0,tvec,dtMoment);

targetGlobal = M(4).moment.globalM3(end);

for m = 1:3
    lo = 0.15;
    hi = 5.50;

    % Bracketing audit.
    momLo = momentEulerAll(lo*M(m).c1Base,M(m).c2, ...
        M(m).from1,M(m).to1,M(m).q1, ...
        M(m).from2,M(m).to2,M(m).q2, ...
        x1_0,x2_0,tvec,dtMoment);
    momHi = momentEulerAll(hi*M(m).c1Base,M(m).c2, ...
        M(m).from1,M(m).to1,M(m).q1, ...
        M(m).from2,M(m).to2,M(m).q2, ...
        x1_0,x2_0,tvec,dtMoment);

    assert(momLo.globalM3(end) <= targetGlobal && ...
           targetGlobal <= momHi.globalM3(end), ...
        'Bisection interval does not bracket the Full target for model %s.', ...
        M(m).name);

    for it = 1:nBisect
        alpha = 0.5*(lo+hi);
        trial = momentEulerAll(alpha*M(m).c1Base,M(m).c2, ...
            M(m).from1,M(m).to1,M(m).q1, ...
            M(m).from2,M(m).to2,M(m).q2, ...
            x1_0,x2_0,tvec,dtMoment);

        if trial.globalM3(end) < targetGlobal
            lo = alpha;
        else
            hi = alpha;
        end
    end

    M(m).alpha = 0.5*(lo+hi);
    M(m).moment = momentEulerAll( ...
        M(m).alpha*M(m).c1Base,M(m).c2, ...
        M(m).from1,M(m).to1,M(m).q1, ...
        M(m).from2,M(m).to2,M(m).q2, ...
        x1_0,x2_0,tvec,dtMoment);
end

for m = 1:4
    M(m).c1 = M(m).alpha*M(m).c1Base;
    M(m).out = buildOutgoingLists(K, ...
        M(m).from1,M(m).to1,M(m).q1, ...
        M(m).from2,M(m).to2,M(m).q2);
    M(m).momentRatio = M(m).moment.globalM3(end)/targetGlobal;

    fprintf('%s: alpha = %.6f, moment ratio to Full = %.6f\n', ...
        M(m).name,M(m).alpha,M(m).momentRatio);
end

% Control-anchored regional FPT threshold.
mStar = M(1).moment.globalM3(end)/K;
M_A = max(1,round(thresholdFactor*mStar*nA));

fprintf('Regional threshold M_A       = %d\n',M_A);
fprintf('M_A / |A|                    = %.4f events/voxel\n',M_A/nA);

%% ------------------- Burden and mechanism ensembles --------------------
for m = 1:4
    fprintf('\nRunning burden/mechanism ensemble for %s...\n',M(m).name);
    M(m).burden = runMechanismEnsembleNSM( ...
        M(m).c1,M(m).c2,M(m).out,A,x1_0,x2_0,tvec, ...
        nRepBurden,seedBase+1000*m,usePar);
end

%% -------------------------- Regional FPT ensembles ----------------------
for m = 1:4
    fprintf('Running regional FPT ensemble for %s...\n',M(m).name);
    M(m).FPT = runRegionFPTensembleNSM( ...
        M(m).c1,M(m).c2,M(m).out,A,x1_0,x2_0, ...
        M_A,Tobs,nRepFPT,seedBase+10000+1000*m,usePar);
end

%% ------------------------ Statistical summaries -------------------------
idxU = 1;
idxH = 2;
idxB = 3;
idxF = 4;

for m = 1:4
    gFinal = M(m).burden.globalDensityMat(:,end);
    lFinal = M(m).burden.localBurdenMat(:,end);
    fptTrunc = min(M(m).FPT,Tobs);

    [M(m).summary.globalMean,M(m).summary.globalCI] = ...
        bootstrapMeanCI(gFinal,nBoot,alphaCI,seedBase+20000+m);

    [M(m).summary.localMean,M(m).summary.localCI] = ...
        bootstrapMeanCI(lFinal,nBoot,alphaCI,seedBase+21000+m);

    [M(m).summary.RMST,M(m).summary.RMST_CI] = ...
        bootstrapMeanCI(fptTrunc,nBoot,alphaCI,seedBase+22000+m);

    M(m).summary.P60 = mean(M(m).FPT<=Tobs);
end

% Local burden ratios relative to U.
localRatio = zeros(4,1);
localRatioCI = zeros(4,2);
localRatio(1) = 1;
localRatioCI(1,:) = [1 1];

for m = 2:4
    [localRatio(m),localRatioCI(m,:)] = bootstrapRatioIndependent( ...
        M(m).burden.localBurdenMat(:,end), ...
        M(idxU).burden.localBurdenMat(:,end), ...
        nBoot,alphaCI,seedBase+23000+m);
end

% RMST effects used in the causal decomposition.
[EH,EH_CI] = bootstrapDifferenceIndependent( ...
    min(M(idxU).FPT,Tobs),min(M(idxH).FPT,Tobs), ...
    nBoot,alphaCI,seedBase+24001);

[EB,EB_CI] = bootstrapDifferenceIndependent( ...
    min(M(idxU).FPT,Tobs),min(M(idxB).FPT,Tobs), ...
    nBoot,alphaCI,seedBase+24002);

[EBgivenH,EBgivenH_CI] = bootstrapDifferenceIndependent( ...
    min(M(idxH).FPT,Tobs),min(M(idxF).FPT,Tobs), ...
    nBoot,alphaCI,seedBase+24003);

[EF,EF_CI] = bootstrapDifferenceIndependent( ...
    min(M(idxU).FPT,Tobs),min(M(idxF).FPT,Tobs), ...
    nBoot,alphaCI,seedBase+24004);

% Additive interaction terms, retained for SI / exported tables.
[IBurden,IBurden_CI] = bootstrapNormalizedBurdenInteraction( ...
    M(idxU).burden.localBurdenMat(:,end), ...
    M(idxH).burden.localBurdenMat(:,end), ...
    M(idxB).burden.localBurdenMat(:,end), ...
    M(idxF).burden.localBurdenMat(:,end), ...
    nBoot,alphaCI,seedBase+25001);

[IFPT,IFPT_CI] = bootstrapAdditiveFPTInteraction( ...
    min(M(idxU).FPT,Tobs),min(M(idxH).FPT,Tobs), ...
    min(M(idxB).FPT,Tobs),min(M(idxF).FPT,Tobs), ...
    nBoot,alphaCI,seedBase+25002);

% Mechanism quantities and exact mass-balance audit.
for m = 1:4
    mech = M(m).burden;

    M(m).mechanism.netDiff60 = ...
        mech.NinMat(:,end)-mech.NoutMat(:,end);

    M(m).mechanism.netReaction60 = ...
        mech.NredMat(:,end)-mech.NfenMat(:,end);

    M(m).mechanism.deltaFe260 = ...
        mech.Fe2A_mat(:,end)-mech.Fe2A_mat(:,1);

    M(m).mechanism.balanceResidual60 = ...
        M(m).mechanism.deltaFe260 - ...
        (M(m).mechanism.netDiff60+M(m).mechanism.netReaction60);
end

maxBalanceResidual = max(arrayfun(@(s) ...
    max(abs(s.mechanism.balanceResidual60)),M));

%% ------------------------------ Console ---------------------------------
fprintf('\n================ Counterfactual summary ====================\n');
for m = 1:4
    fprintf('%s: SSA global ratio to Full = %.5f; ', ...
        M(m).name,M(m).summary.globalMean/M(idxF).summary.globalMean);
    fprintf('local burden/U = %.4f; RMST60 = %.3f; P60 = %.3f\n', ...
        localRatio(m),M(m).summary.RMST,M(m).summary.P60);
end
fprintf('E_H       = %.3f [%.3f, %.3f] min\n',EH,EH_CI(1),EH_CI(2));
fprintf('E_B       = %.3f [%.3f, %.3f] min\n',EB,EB_CI(1),EB_CI(2));
fprintf('E_B|H     = %.3f [%.3f, %.3f] min\n', ...
    EBgivenH,EBgivenH_CI(1),EBgivenH_CI(2));
fprintf('E_F       = %.3f [%.3f, %.3f] min\n',EF,EF_CI(1),EF_CI(2));
fprintf('I_HB burden ratio = %.4f [%.4f, %.4f]\n', ...
    IBurden,IBurden_CI(1),IBurden_CI(2));
fprintf('I_HB FPT    = %.4f [%.4f, %.4f] min\n', ...
    IFPT,IFPT_CI(1),IFPT_CI(2));
fprintf('Max Fe2 mass-balance residual = %.3e events/voxel\n', ...
    maxBalanceResidual);
fprintf('============================================================\n');

%% ============================== Main Fig. 3E ============================
figMain = figure('Units','centimeters','Position',[1.0 1.5 18.4 9.0], ...
    'Color','w');
tl = tiledlayout(figMain,2,2,'TileSpacing','compact','Padding','compact');

% A. Global matching trajectories.
ax = nexttile(tl,1); hold(ax,'on');
hMainModels = gobjects(4,1);
for m = 1:4
    hMainModels(m) = plotMeanSEMWithHandle( ...
        ax,tvec,M(m).burden.globalDensityMat,M(m).color,0.10);
end
xlabel(ax,'Time (min)');
ylabel(ax,'Global burden per voxel');
title(ax,'Global burden matching','FontSize',8.5,'FontWeight','normal');
legend(ax,hMainModels,modelLabels,'Location','northwest', ...
    'Box','off','FontSize',5.7);
text(ax,0.02,0.97,'A','Units','normalized','FontWeight','bold', ...
    'FontSize',12,'VerticalAlignment','top');
text(ax,0.98,0.07,sprintf('max final mismatch = %.2f%%', ...
    100*max(abs(arrayfun(@(s) ...
    s.summary.globalMean/M(idxF).summary.globalMean-1,M)))), ...
    'Units','normalized','HorizontalAlignment','right','FontSize',5.8);
xlim(ax,[0 Tobs]);
formatAxes(ax);

% B. Local burden decomposition in A_0.85.
ax = nexttile(tl,2); hold(ax,'on');
for m = 1:4
    errorbar(ax,m,localRatio(m), ...
        max(localRatio(m)-localRatioCI(m,1),0), ...
        max(localRatioCI(m,2)-localRatio(m),0), ...
        'o','Color',M(m).color,'MarkerFaceColor',M(m).color, ...
        'MarkerSize',4.8,'LineWidth',0.95,'CapSize',5);
end
plot(ax,[0.5 4.5],[1 1],':','Color',[0.60 0.60 0.60], ...
    'LineWidth',0.75,'HandleVisibility','off');
set(ax,'XLim',[0.5 4.5],'XTick',1:4,'XTickLabel',modelNames);
ylabel(ax,'Local burden at 60 min / U');
title(ax,'Local burden decomposition','FontSize',8.5,'FontWeight','normal');
text(ax,0.02,0.97,'B','Units','normalized','FontWeight','bold', ...
    'FontSize',12,'VerticalAlignment','top');
text(ax,0.98,0.07,sprintf('$I_{HB}^{ratio}=%.2f$ [%.2f, %.2f]', ...
    IBurden,IBurden_CI(1),IBurden_CI(2)), ...
    'Units','normalized','HorizontalAlignment','right', ...
    'Interpreter','latex','FontSize',5.8);
formatAxes(ax);

% C. Regional FPT cumulative crossing probabilities with censoring.
ax = nexttile(tl,3); hold(ax,'on');
for m = 1:4
    [tx,Fx] = censoredCrossingCDF(M(m).FPT,Tobs);
    stairs(ax,tx,Fx,'Color',M(m).color,'LineWidth',1.20);
end
xlabel(ax,'Regional FPT (min)');
ylabel(ax,'Cumulative crossing probability');
title(ax,'Local threshold crossing','FontSize',8.5,'FontWeight','normal');
legend(ax,modelLabels,'Location','northwest','Box','off','FontSize',5.7);
text(ax,0.02,0.97,'C','Units','normalized','FontWeight','bold', ...
    'FontSize',12,'VerticalAlignment','top');
text(ax,0.98,0.08,sprintf('$P_{60}$: U %.2f, H %.2f, B %.2f, F %.2f', ...
    M(1).summary.P60,M(2).summary.P60,M(3).summary.P60,M(4).summary.P60), ...
    'Units','normalized','HorizontalAlignment','right', ...
    'Interpreter','latex','FontSize',5.5);
xlim(ax,[0 Tobs]); ylim(ax,[0 1]);
formatAxes(ax);

% D. Factorial RMST effects.
ax = nexttile(tl,4); hold(ax,'on');
effectVal = [EH,EB,EBgivenH];
effectLo  = [EH_CI(1),EB_CI(1),EBgivenH_CI(1)];
effectHi  = [EH_CI(2),EB_CI(2),EBgivenH_CI(2)];
effectCol = [colH;colB;colF];

for k = 1:3
    errorbar(ax,k,effectVal(k), ...
        max(effectVal(k)-effectLo(k),0), ...
        max(effectHi(k)-effectVal(k),0), ...
        'o','Color',effectCol(k,:), ...
        'MarkerFaceColor',effectCol(k,:), ...
        'MarkerSize',4.8,'LineWidth',0.95,'CapSize',5);
end
plot(ax,[0.5 3.5],[0 0],':','Color',[0.60 0.60 0.60], ...
    'LineWidth',0.75,'HandleVisibility','off');
set(ax,'XLim',[0.5 3.5],'XTick',1:3, ...
    'XTickLabel',{'hotspot alone','bottleneck alone','bottleneck | hotspot'});
ylabel(ax,'RMST advancement (min)');
title(ax,'Counterfactual FPT effects','FontSize',8.5,'FontWeight','normal');
text(ax,0.02,0.97,'D','Units','normalized','FontWeight','bold', ...
    'FontSize',12,'VerticalAlignment','top');
text(ax,0.98,0.08,sprintf('$E_F=%.2f$ min',EF), ...
    'Units','normalized','HorizontalAlignment','right', ...
    'Interpreter','latex','FontSize',5.8);
formatAxes(ax);

if SAVE_FIG
    exportgraphics(figMain,[outMain '.pdf'],'ContentType','vector');
    exportgraphics(figMain,[outMain '.png'],'Resolution',600);
    exportgraphics(figMain,[outMain '.tif'],'Resolution',600);
    savefig(figMain,[outMain '.fig']);
end

%% ========================= SI detailed mechanisms =======================
figSI = figure('Units','centimeters','Position',[1.0 1.0 19.0 11.8], ...
    'Color','w');
tlSI = tiledlayout(figSI,2,3,'TileSpacing','compact','Padding','compact');

% A. Fe2+ inventory.
ax = nexttile(tlSI,1); hold(ax,'on');
hSIModels = gobjects(4,1);
for m = 1:4
    hSIModels(m) = plotMeanSEMWithHandle( ...
        ax,tvec,M(m).burden.Fe2A_mat,M(m).color,0.10);
end
xlabel(ax,'Time (min)');
ylabel(ax,'Fe^{2+} per voxel in A');
title(ax,'Local Fe^{2+} inventory','FontSize',8.3,'FontWeight','normal');
legend(ax,hSIModels,modelNames,'Location','best', ...
    'Box','off','FontSize',5.6);
text(ax,0.02,0.97,'A','Units','normalized','FontWeight','bold', ...
    'FontSize',11.5,'VerticalAlignment','top');
formatAxes(ax);

% B. Cumulative Fe2+ influx.
ax = nexttile(tlSI,2); hold(ax,'on');
for m = 1:4
    plotMeanSEM(ax,tvec,M(m).burden.NinMat,M(m).color,0.10);
end
xlabel(ax,'Time (min)');
ylabel(ax,'Cumulative influx per voxel in A');
title(ax,'Boundary Fe^{2+} influx','FontSize',8.3,'FontWeight','normal');
text(ax,0.02,0.97,'B','Units','normalized','FontWeight','bold', ...
    'FontSize',11.5,'VerticalAlignment','top');
formatAxes(ax);

% C. Cumulative Fe2+ efflux.
ax = nexttile(tlSI,3); hold(ax,'on');
for m = 1:4
    plotMeanSEM(ax,tvec,M(m).burden.NoutMat,M(m).color,0.10);
end
xlabel(ax,'Time (min)');
ylabel(ax,'Cumulative efflux per voxel in A');
title(ax,'Boundary Fe^{2+} efflux','FontSize',8.3,'FontWeight','normal');
text(ax,0.02,0.97,'C','Units','normalized','FontWeight','bold', ...
    'FontSize',11.5,'VerticalAlignment','top');
formatAxes(ax);

% D. Local R2 reduction.
ax = nexttile(tlSI,4); hold(ax,'on');
for m = 1:4
    plotMeanSEM(ax,tvec,M(m).burden.NredMat,M(m).color,0.10);
end
xlabel(ax,'Time (min)');
ylabel(ax,'Cumulative R2 events per voxel in A');
title(ax,'Local Fe^{3+} reduction','FontSize',8.3,'FontWeight','normal');
text(ax,0.02,0.97,'D','Units','normalized','FontWeight','bold', ...
    'FontSize',11.5,'VerticalAlignment','top');
formatAxes(ax);

% E. Local R1 / Fenton consumption.
ax = nexttile(tlSI,5); hold(ax,'on');
for m = 1:4
    plotMeanSEM(ax,tvec,M(m).burden.NfenMat,M(m).color,0.10);
end
xlabel(ax,'Time (min)');
ylabel(ax,'Cumulative R1 events per voxel in A');
title(ax,'Local Fenton consumption','FontSize',8.3,'FontWeight','normal');
text(ax,0.02,0.97,'E','Units','normalized','FontWeight','bold', ...
    'FontSize',11.5,'VerticalAlignment','top');
formatAxes(ax);

% F. Exact Fe2+ regional mass-balance closure.
ax = nexttile(tlSI,6); hold(ax,'on');
balanceNames = {'net diffusion','net reaction','observed \DeltaFe^{2+}'};
offsets = [-0.18 -0.06 0.06 0.18];

for m = 1:4
    vals = [mean(M(m).mechanism.netDiff60), ...
            mean(M(m).mechanism.netReaction60), ...
            mean(M(m).mechanism.deltaFe260)];
    plot(ax,(1:3)+offsets(m),vals,'-o','Color',M(m).color, ...
        'MarkerFaceColor',M(m).color,'MarkerSize',3.8, ...
        'LineWidth',0.90);
end
plot(ax,[0.5 3.5],[0 0],':','Color',[0.60 0.60 0.60], ...
    'LineWidth',0.75,'HandleVisibility','off');
set(ax,'XLim',[0.5 3.5],'XTick',1:3,'XTickLabel',balanceNames);
ylabel(ax,'Events per voxel in A at 60 min');
title(ax,'Regional Fe^{2+} mass-balance closure', ...
    'FontSize',8.3,'FontWeight','normal');
legend(ax,modelNames,'Location','best','Box','off','FontSize',5.5);
text(ax,0.02,0.97,'F','Units','normalized','FontWeight','bold', ...
    'FontSize',11.5,'VerticalAlignment','top');
text(ax,0.98,0.06,sprintf('max residual = %.1e',maxBalanceResidual), ...
    'Units','normalized','HorizontalAlignment','right','FontSize',5.5);
formatAxes(ax);

if SAVE_FIG
    exportgraphics(figSI,[outSI '.pdf'],'ContentType','vector');
    exportgraphics(figSI,[outSI '.png'],'Resolution',600);
    exportgraphics(figSI,[outSI '.tif'],'Resolution',600);
    savefig(figSI,[outSI '.fig']);
end

%% --------------------------- Export summaries ---------------------------
% Numeric summary table for CSV export.
modelTable = table( ...
    string(modelNames(:)), ...
    [M.alpha].', ...
    [M.momentRatio].', ...
    arrayfun(@(s)s.summary.globalMean/M(idxF).summary.globalMean,M).', ...
    arrayfun(@(s)s.summary.localMean,M).', ...
    localRatio, ...
    arrayfun(@(s)s.summary.RMST,M).', ...
    arrayfun(@(s)s.summary.P60,M).', ...
    'VariableNames',{'Model','AlphaMatch','MomentGlobalRatio', ...
    'SSAGlobalRatio','LocalBurdenMean','LocalBurdenRatioToU', ...
    'RMST60','Pcross60'});

contrastTable = table( ...
    ["Hotspot_alone";"Bottleneck_alone";"Bottleneck_given_hotspot"; ...
     "Full_total";"Burden_ratio_interaction";"FPT_interaction"], ...
    [EH;EB;EBgivenH;EF;IBurden;IFPT], ...
    [EH_CI(1);EB_CI(1);EBgivenH_CI(1);EF_CI(1); ...
     IBurden_CI(1);IFPT_CI(1)], ...
    [EH_CI(2);EB_CI(2);EBgivenH_CI(2);EF_CI(2); ...
     IBurden_CI(2);IFPT_CI(2)], ...
    'VariableNames',{'Contrast','Estimate','CI_L','CI_U'});

if SAVE_DATA
    writetable(modelTable,outModelCSV);
    writetable(contrastTable,outContrCSV);

    metadata = struct();
    metadata.RUN_MODE = RUN_MODE;
    metadata.seedBase = seedBase;
    metadata.nRepBurden = nRepBurden;
    metadata.nRepFPT = nRepFPT;
    metadata.nBoot = nBoot;
    metadata.qA = qA;
    metadata.nA = nA;
    metadata.K = K;
    metadata.M_A = M_A;
    metadata.Tobs = Tobs;
    metadata.dtRecord = dtRecord;
    metadata.parameters = struct( ...
        'D0_Fe2',D0_Fe2,'D0_Fe3',D0_Fe3, ...
        'c1_0',c1_0,'etaS',etaS,'c2_0',c2_0, ...
        'FeTot_per_voxel',FeTot_per_voxel, ...
        'fracFe2_0',fracFe2_0);

    save(outData,'metadata','M','A','Amap','c1Center','c1Window', ...
        'win','mask','x','y', ...
        'Xg','Yg','tvec','targetGlobal','mStar','M_A', ...
        'localRatio','localRatioCI', ...
        'EH','EH_CI','EB','EB_CI','EBgivenH','EBgivenH_CI', ...
        'EF','EF_CI','IBurden','IBurden_CI','IFPT','IFPT_CI', ...
        'modelTable','contrastTable','-v7.3');
end

fprintf('\nCompleted. Outputs:\n');
fprintf('  %s.{pdf,png,tif,fig}\n',outMain);
fprintf('  %s.{pdf,png,tif,fig}\n',outSI);
fprintf('  %s\n',outModelCSV);
fprintf('  %s\n',outContrCSV);
fprintf('  %s\n',outData);

%% ======================== Analysis-specific functions ===================

function out = runMechanismEnsembleNSM(c1,c2,outgoing,A,x1_0,x2_0, ...
    tvec,nRep,seedBase,usePar)

    K = numel(c1);
    nT = numel(tvec);

    globalDensityMat = zeros(nRep,nT);
    localBurdenMat   = zeros(nRep,nT);
    Fe2A_mat         = zeros(nRep,nT);
    NinMat           = zeros(nRep,nT);
    NoutMat          = zeros(nRep,nT);
    NredMat          = zeros(nRep,nT);
    NfenMat          = zeros(nRep,nT);

    if usePar
        parfor r = 1:nRep
            rng(seedBase+r,'twister');
            one = simulateMechanismOneNSM( ...
                c1,c2,outgoing,A,x1_0,x2_0,tvec);

            globalDensityMat(r,:) = one.globalDensityRec;
            localBurdenMat(r,:)   = one.localBurdenRec;
            Fe2A_mat(r,:)         = one.Fe2Arec;
            NinMat(r,:)           = one.NinRec;
            NoutMat(r,:)          = one.NoutRec;
            NredMat(r,:)          = one.NredRec;
            NfenMat(r,:)          = one.NfenRec;
        end
    else
        for r = 1:nRep
            rng(seedBase+r,'twister');
            one = simulateMechanismOneNSM( ...
                c1,c2,outgoing,A,x1_0,x2_0,tvec);

            globalDensityMat(r,:) = one.globalDensityRec;
            localBurdenMat(r,:)   = one.localBurdenRec;
            Fe2A_mat(r,:)         = one.Fe2Arec;
            NinMat(r,:)           = one.NinRec;
            NoutMat(r,:)          = one.NoutRec;
            NredMat(r,:)          = one.NredRec;
            NfenMat(r,:)          = one.NfenRec;

            if mod(r,max(1,round(nRep/10)))==0
                fprintf('  completed %d / %d burden trajectories\n',r,nRep);
            end
        end
    end

    out.globalDensityMat = globalDensityMat;
    out.localBurdenMat   = localBurdenMat;
    out.Fe2A_mat         = Fe2A_mat;
    out.NinMat           = NinMat;
    out.NoutMat          = NoutMat;
    out.NredMat          = NredMat;
    out.NfenMat          = NfenMat;
end

function one = simulateMechanismOneNSM(c1,c2,out,A,x1_0,x2_0,tvec)

    K = numel(c1);
    nA = nnz(A);
    nT = numel(tvec);

    x1 = x1_0(:);
    x2 = x2_0(:);

    voxelRate = x1.*(c1+out.q1sum)+x2.*(c2+out.q2sum);
    tree = fenwickBuild(voxelRate);

    globalROS = 0;
    localROS  = 0;
    Nin  = 0;
    Nout = 0;
    Nred = 0;
    Nfen = 0;

    globalDensityRec = zeros(1,nT);
    localBurdenRec   = zeros(1,nT);
    Fe2Arec          = zeros(1,nT);
    NinRec           = zeros(1,nT);
    NoutRec          = zeros(1,nT);
    NredRec          = zeros(1,nT);
    NfenRec          = zeros(1,nT);

    Fe2Arec(1) = sum(x1(A))/nA;

    t = tvec(1);
    outIdx = 1;
    Tend = tvec(end);

    while t<Tend
        a0 = fenwickTotal(tree);
        if ~(a0>0) || ~isfinite(a0)
            break;
        end

        tau = -log(rand)/a0;
        tNext = t+tau;

        while outIdx<nT && tvec(outIdx+1)<=min(tNext,Tend)
            outIdx = outIdx+1;
            globalDensityRec(outIdx) = globalROS/K;
            localBurdenRec(outIdx)   = localROS/nA;
            Fe2Arec(outIdx)          = sum(x1(A))/nA;
            NinRec(outIdx)           = Nin/nA;
            NoutRec(outIdx)          = Nout/nA;
            NredRec(outIdx)          = Nred/nA;
            NfenRec(outIdx)          = Nfen/nA;
        end

        if tNext>Tend
            break;
        end

        k = fenwickFind(tree,rand*a0);

        aR1 = c1(k)*x1(k);
        aR2 = c2(k)*x2(k);
        aD1 = out.q1sum(k)*x1(k);
        aD2 = out.q2sum(k)*x2(k);
        aLocal = aR1+aR2+aD1+aD2;

        if aLocal<=0
            old = voxelRate(k);
            voxelRate(k) = 0;
            tree = fenwickAdd(tree,k,-old);
            t = tNext;
            continue;
        end

        u = rand*aLocal;
        touched = k;

        if u<=aR1
            x1(k) = x1(k)-1;
            x2(k) = x2(k)+1;
            globalROS = globalROS+1;

            if A(k)
                localROS = localROS+1;
                Nfen = Nfen+1;
            end

        elseif u<=aR1+aR2
            x2(k) = x2(k)-1;
            x1(k) = x1(k)+1;

            if A(k)
                Nred = Nred+1;
            end

        elseif u<=aR1+aR2+aD1
            j = chooseDestination(out.to1{k},out.q1{k});
            x1(k) = x1(k)-1;
            x1(j) = x1(j)+1;
            touched = [k;j]; %#ok<AGROW>

            if ~A(k) && A(j)
                Nin = Nin+1;
            elseif A(k) && ~A(j)
                Nout = Nout+1;
            end

        else
            j = chooseDestination(out.to2{k},out.q2{k});
            x2(k) = x2(k)-1;
            x2(j) = x2(j)+1;
            touched = [k;j]; %#ok<AGROW>
        end

        [voxelRate,tree] = refreshVoxelRates( ...
            unique(touched),voxelRate,tree, ...
            x1,x2,c1,c2,out.q1sum,out.q2sum);

        t = tNext;
    end

    while outIdx<nT
        outIdx = outIdx+1;
        globalDensityRec(outIdx) = globalROS/K;
        localBurdenRec(outIdx)   = localROS/nA;
        Fe2Arec(outIdx)          = sum(x1(A))/nA;
        NinRec(outIdx)           = Nin/nA;
        NoutRec(outIdx)          = Nout/nA;
        NredRec(outIdx)          = Nred/nA;
        NfenRec(outIdx)          = Nfen/nA;
    end

    one.globalDensityRec = globalDensityRec;
    one.localBurdenRec   = localBurdenRec;
    one.Fe2Arec          = Fe2Arec;
    one.NinRec           = NinRec;
    one.NoutRec          = NoutRec;
    one.NredRec          = NredRec;
    one.NfenRec          = NfenRec;
end

function FPT = runRegionFPTensembleNSM(c1,c2,outgoing,A,x1_0,x2_0, ...
    M_A,Tobs,nRep,seedBase,usePar)

    FPT = inf(nRep,1);

    if usePar
        parfor r = 1:nRep
            rng(seedBase+r,'twister');
            FPT(r) = simulateRegionFPTOneNSM( ...
                c1,c2,outgoing,A,x1_0,x2_0,M_A,Tobs);
        end
    else
        for r = 1:nRep
            rng(seedBase+r,'twister');
            FPT(r) = simulateRegionFPTOneNSM( ...
                c1,c2,outgoing,A,x1_0,x2_0,M_A,Tobs);

            if mod(r,max(1,round(nRep/10)))==0
                fprintf('  completed %d / %d FPT trajectories\n',r,nRep);
            end
        end
    end
end

function fpt = simulateRegionFPTOneNSM( ...
    c1,c2,out,A,x1_0,x2_0,M_A,Tobs)

    x1 = x1_0(:);
    x2 = x2_0(:);

    localROS = 0;
    fpt = inf;

    voxelRate = x1.*(c1+out.q1sum)+x2.*(c2+out.q2sum);
    tree = fenwickBuild(voxelRate);

    t = 0;

    while t<Tobs
        a0 = fenwickTotal(tree);
        if ~(a0>0) || ~isfinite(a0)
            break;
        end

        tau = -log(rand)/a0;
        tNext = t+tau;

        if tNext>Tobs
            break;
        end

        k = fenwickFind(tree,rand*a0);

        aR1 = c1(k)*x1(k);
        aR2 = c2(k)*x2(k);
        aD1 = out.q1sum(k)*x1(k);
        aD2 = out.q2sum(k)*x2(k);
        aLocal = aR1+aR2+aD1+aD2;

        if aLocal<=0
            old = voxelRate(k);
            voxelRate(k) = 0;
            tree = fenwickAdd(tree,k,-old);
            t = tNext;
            continue;
        end

        u = rand*aLocal;
        touched = k;

        if u<=aR1
            x1(k) = x1(k)-1;
            x2(k) = x2(k)+1;

            if A(k)
                localROS = localROS+1;
                if localROS>=M_A
                    fpt = tNext;
                    return;
                end
            end

        elseif u<=aR1+aR2
            x2(k) = x2(k)-1;
            x1(k) = x1(k)+1;

        elseif u<=aR1+aR2+aD1
            j = chooseDestination(out.to1{k},out.q1{k});
            x1(k) = x1(k)-1;
            x1(j) = x1(j)+1;
            touched = [k;j]; %#ok<AGROW>

        else
            j = chooseDestination(out.to2{k},out.q2{k});
            x2(k) = x2(k)-1;
            x2(j) = x2(j)+1;
            touched = [k;j]; %#ok<AGROW>
        end

        [voxelRate,tree] = refreshVoxelRates( ...
            unique(touched),voxelRate,tree, ...
            x1,x2,c1,c2,out.q1sum,out.q2sum);

        t = tNext;
    end
end

function hLine = plotMeanSEMWithHandle(ax,t,mat,lineColor,alphaVal)
    mu = mean(mat,1);
    se = std(mat,0,1)/sqrt(size(mat,1));

    fill(ax,[t,fliplr(t)],[mu-se,fliplr(mu+se)],lineColor, ...
        'FaceAlpha',alphaVal,'EdgeColor','none', ...
        'HandleVisibility','off');

    hLine = plot(ax,t,mu,'Color',lineColor,'LineWidth',1.10);
end

function [mu,ci] = bootstrapMeanCI(x,nBoot,alpha,seed)
    x = x(:);
    n = numel(x);
    mu = mean(x);
    boot = zeros(nBoot,1);

    rng(seed,'twister');
    for b = 1:nBoot
        boot(b) = mean(x(randi(n,n,1)));
    end

    ci = [quantileLocal(boot,alpha/2), ...
          quantileLocal(boot,1-alpha/2)];
end

function [ratio,ci] = bootstrapRatioIndependent(x,y,nBoot,alpha,seed)
    x = x(:);
    y = y(:);
    nx = numel(x);
    ny = numel(y);

    ratio = mean(x)/max(mean(y),eps);
    boot = zeros(nBoot,1);

    rng(seed,'twister');
    for b = 1:nBoot
        bx = x(randi(nx,nx,1));
        by = y(randi(ny,ny,1));
        boot(b) = mean(bx)/max(mean(by),eps);
    end

    ci = [quantileLocal(boot,alpha/2), ...
          quantileLocal(boot,1-alpha/2)];
end

function [delta,ci] = bootstrapDifferenceIndependent(a,b,nBoot,alpha,seed)
    a = a(:);
    b = b(:);
    na = numel(a);
    nb = numel(b);

    delta = mean(a)-mean(b);
    boot = zeros(nBoot,1);

    rng(seed,'twister');
    for k = 1:nBoot
        aa = a(randi(na,na,1));
        bb = b(randi(nb,nb,1));
        boot(k) = mean(aa)-mean(bb);
    end

    ci = [quantileLocal(boot,alpha/2), ...
          quantileLocal(boot,1-alpha/2)];
end

function [interaction,ci] = bootstrapNormalizedBurdenInteraction( ...
    U,H,B,F,nBoot,alpha,seed)

    U = U(:); H = H(:); B = B(:); F = F(:);
    nU = numel(U); nH = numel(H);
    nB = numel(B); nF = numel(F);

    % Dimensionless interaction on the same ratio scale as panel B:
    % I_HB = F/U - H/U - B/U + 1.
    meanU = max(mean(U),eps);
    interaction = mean(F)/meanU ...
                - mean(H)/meanU ...
                - mean(B)/meanU + 1;

    boot = zeros(nBoot,1);

    rng(seed,'twister');
    for k = 1:nBoot
        UU = U(randi(nU,nU,1));
        HH = H(randi(nH,nH,1));
        BB = B(randi(nB,nB,1));
        FF = F(randi(nF,nF,1));

        meanUU = max(mean(UU),eps);
        boot(k) = mean(FF)/meanUU ...
                - mean(HH)/meanUU ...
                - mean(BB)/meanUU + 1;
    end

    ci = [quantileLocal(boot,alpha/2), ...
          quantileLocal(boot,1-alpha/2)];
end

function [interaction,ci] = bootstrapAdditiveFPTInteraction( ...
    U,H,B,F,nBoot,alpha,seed)

    U = U(:); H = H(:); B = B(:); F = F(:);
    nU = numel(U); nH = numel(H);
    nB = numel(B); nF = numel(F);

    % Effects are RMST advancement: E_M = RMST_U - RMST_M.
    interaction = (mean(U)-mean(F)) ...
                - (mean(U)-mean(H)) ...
                - (mean(U)-mean(B));

    boot = zeros(nBoot,1);

    rng(seed,'twister');
    for k = 1:nBoot
        UU = U(randi(nU,nU,1));
        HH = H(randi(nH,nH,1));
        BB = B(randi(nB,nB,1));
        FF = F(randi(nF,nF,1));

        boot(k) = (mean(UU)-mean(FF)) ...
                - (mean(UU)-mean(HH)) ...
                - (mean(UU)-mean(BB));
    end

    ci = [quantileLocal(boot,alpha/2), ...
          quantileLocal(boot,1-alpha/2)];
end

function win = buildCircularWindowsSparse(coordX,coordY,radius)
    coordX = coordX(:);
    coordY = coordY(:);
    K = numel(coordX);

    I = cell(K,1);
    J = cell(K,1);

    radius2 = radius^2 + 1e-14;

    for k = 1:K
        idx = find((coordX-coordX(k)).^2 + ...
                   (coordY-coordY(k)).^2 <= radius2);

        I{k} = repmat(k,numel(idx),1);
        J{k} = idx(:);
    end

    I = vertcat(I{:});
    J = vertcat(J{:});

    win.W = sparse(I,J,1,K,K);
    win.size = full(sum(win.W,2));

    if any(win.size<1)
        error('At least one local window is empty.');
    end
end

function v = normalizeVector(v)
    v = v(:);
    lo = min(v);
    hi = max(v);
    v = (v-lo)/(hi-lo+eps);
end

function maskTop = topFractionMask(score,topFraction)
    score = score(:);
    K = numel(score);
    nTop = max(1,ceil(topFraction*K));

    [~,order] = sortrows([-score,(1:K).'],[1 2]);

    maskTop = false(K,1);
    maskTop(order(1:nTop)) = true;
end

function [xCDF,FCDF] = censoredCrossingCDF(t,Tobs)
    t = t(:);
    nTotal = numel(t);

    events = sort(t(isfinite(t) & t<=Tobs));
    nEvents = numel(events);

    if nEvents==0
        xCDF = [0;Tobs];
        FCDF = [0;0];
    else
        xCDF = [0;events;Tobs];
        FCDF = [0;(1:nEvents)'/nTotal;nEvents/nTotal];
    end
end

%% ================= Reused validated Fig. 3B-C helpers ==================

function usePar = configureParallel(requested,autoStart,nWorkers)
    usePar = false;
    if ~requested
        return;
    end
    try
        hasPCT = license('test','Distrib_Computing_Toolbox');
        if ~hasPCT
            warning('Parallel Computing Toolbox unavailable; using serial loops.');
            return;
        end
        pool = gcp('nocreate');
        if isempty(pool) && autoStart
            pool = parpool('local',nWorkers);
        end
        usePar = ~isempty(pool);
    catch ME
        warning('Parallel setup failed (%s); using serial loops.',ME.message);
        usePar = false;
    end
end

function geom = buildCurrentCristaGeometry(nx,ny,aOut,bOut,aIn,bIn, ...
    ellG,ellS,connMatrix,connLumen,connNeck)

    x = linspace(-1.23,1.23,nx);
    y = linspace(-0.70,0.70,ny);
    [X,Y] = meshgrid(x,y);

    rhoOut = sqrt((X./aOut).^2 + (Y./bOut).^2);
    rhoIn  = sqrt((X./aIn ).^2 + (Y./bIn ).^2);
    maskOut = rhoOut <= 1;
    OmegaR  = rhoIn <= 1;

    cristaSpec = defaultCristaSpec();
    nCrista = numel(cristaSpec);
    nPtsCrista = 200;

    crista = cell(1,nCrista);
    cristaLumenMask = false(size(X));
    for j = 1:nCrista
        crTemp = generateCristaAnchored(cristaSpec(j),aIn,bIn,nPtsCrista);
        thisLumen = inpolygon(X,Y,crTemp.polyX,crTemp.polyY) & OmegaR;
        crTemp.lumenMask = thisLumen;
        crista{j} = crTemp;
        cristaLumenMask = cristaLumenMask | thisLumen;
    end
    matrixMask = OmegaR & ~cristaLumenMask;

    thetaAll = linspace(-pi,pi,1200);
    xOuter = aOut*cos(thetaAll);
    yOuter = bOut*sin(thetaAll);

    gapMask = true(size(thetaAll));
    for j = 1:nCrista
        cr = crista{j};
        dth = angdiff_local(thetaAll,cr.theta0);
        gapMask = gapMask & (abs(dth) > cr.gapHalfTheta);
    end

    xInnerPlot = aIn*cos(thetaAll);
    yInnerPlot = bIn*sin(thetaAll);
    xInnerPlot(~gapMask) = NaN;
    yInnerPlot(~gapMask) = NaN;

    memX = xInnerPlot(gapMask);
    memY = yInnerPlot(gapMask);
    for j = 1:nCrista
        cr = crista{j};
        memX = [memX(:); cr.xMem(:)]; %#ok<AGROW>
        memY = [memY(:); cr.yMem(:)]; %#ok<AGROW>
    end
    sampleStep = 3;
    memX = memX(1:sampleStep:end);
    memY = memY(1:sampleStep:end);

    dR = inf(size(X));
    for p = 1:numel(memX)
        dR = min(dR,hypot(X-memX(p),Y-memY(p)));
    end
    dR(~OmegaR) = NaN;

    chiBase = connMatrix*ones(size(X));
    chiBase(cristaLumenMask) = connLumen;
    chiBase(~OmegaR) = NaN;

    Kneck = zeros(size(X));
    for j = 1:nCrista
        cr = crista{j};
        dNeck = inf(size(X));
        for p = 1:numel(cr.neckX)
            dNeck = min(dNeck,hypot(X-cr.neckX(p),Y-cr.neckY(p)));
        end
        sigmaNeck = 0.75*cr.neckRadius;
        thisK = exp(-(dNeck.^2)/(2*sigmaNeck^2));
        Kneck = max(Kneck,thisK);
    end
    Kneck(~OmegaR) = NaN;

    chiG = chiBase - (chiBase-connNeck).*Kneck;
    chiG(~OmegaR) = NaN;
    chiG = min(max(chiG,connNeck),connMatrix);
    chiG(~OmegaR) = NaN;

    G = (1-exp(-dR./ellG)).*chiG;
    G(~OmegaR) = NaN;
    G = min(max(G,0),1);

    S0 = zeros(size(X));
    for p = 1:numel(memX)
        S0 = S0 + exp(-((X-memX(p)).^2 + (Y-memY(p)).^2)/(2*ellS^2));
    end
    S0(~OmegaR) = NaN;
    S = normalizeInside(S0,OmegaR);

    geom.x=x; geom.y=y; geom.X=X; geom.Y=Y;
    geom.maskOut=maskOut; geom.OmegaR=OmegaR;
    geom.matrixMask=matrixMask; geom.cristaLumenMask=cristaLumenMask;
    geom.crista=crista; geom.xOuter=xOuter; geom.yOuter=yOuter;
    geom.xInnerPlot=xInnerPlot; geom.yInnerPlot=yInnerPlot;
    geom.G=G; geom.S=S;
end

function sp = defaultCristaSpec()
    sp = repmat(struct('x0',0,'fromTop',true,'L',0,'w',0, ...
        'amp',0,'phase',0,'bias',0),1,7);
    sp(1).x0=-0.68; sp(1).fromTop=true;  sp(1).L=0.34; sp(1).w=0.043; sp(1).amp=0.026; sp(1).phase=0.30; sp(1).bias=-0.030;
    sp(2).x0=-0.48; sp(2).fromTop=false; sp(2).L=0.42; sp(2).w=0.050; sp(2).amp=0.032; sp(2).phase=1.40; sp(2).bias= 0.020;
    sp(3).x0=-0.04; sp(3).fromTop=true;  sp(3).L=0.68; sp(3).w=0.060; sp(3).amp=0.018; sp(3).phase=2.20; sp(3).bias= 0.000;
    sp(4).x0= 0.39; sp(4).fromTop=false; sp(4).L=0.42; sp(4).w=0.050; sp(4).amp=0.026; sp(4).phase=0.70; sp(4).bias=-0.020;
    sp(5).x0= 0.55; sp(5).fromTop=true;  sp(5).L=0.54; sp(5).w=0.056; sp(5).amp=0.020; sp(5).phase=2.10; sp(5).bias= 0.010;
    sp(6).x0= 0.67; sp(6).fromTop=false; sp(6).L=0.37; sp(6).w=0.046; sp(6).amp=0.024; sp(6).phase=3.00; sp(6).bias= 0.030;
    sp(7).x0= 0.79; sp(7).fromTop=true;  sp(7).L=0.30; sp(7).w=0.040; sp(7).amp=0.020; sp(7).phase=1.10; sp(7).bias=-0.012;
end

function cr = generateCristaAnchored(sp,aIn,bIn,nPts)
    x0 = sp.x0;
    yb = bIn*sqrt(max(0,1-(x0/aIn)^2));
    if sp.fromTop, y0=+yb; else, y0=-yb; end

    theta0 = atan2(y0/bIn,x0/aIn);
    nIn = -[x0/(aIn^2),y0/(bIn^2)];
    nIn = nIn/norm(nIn);
    tVec = [-aIn*sin(theta0),bIn*cos(theta0)];
    tVec = tVec/norm(tVec);

    u = linspace(0,1,nPts);
    shift = sp.amp*sin(2*pi*1.03*u+sp.phase) + sp.bias*sin(pi*u);
    xc = x0 + nIn(1)*sp.L*u + tVec(1)*shift;
    yc = y0 + nIn(2)*sp.L*u + tVec(2)*shift;

    inside = (xc./aIn).^2 + (yc./bIn).^2 <= 0.985;
    lastValid = find(inside,1,'last');
    if isempty(lastValid) || lastValid < 10
        lastValid = max(10,round(0.65*nPts));
    end
    xc=xc(1:lastValid); yc=yc(1:lastValid);
    u=linspace(0,1,numel(xc));

    wBody=sp.w; wOpen=0.32*sp.w; wTip=0.68*sp.w;
    w=widthProfileOpenToClosed(u,wBody,wOpen,wTip,0.22,0.18);

    dxc=gradient(xc); dyc=gradient(yc);
    L=hypot(dxc,dyc)+eps;
    qx=-dyc./L; qy=dxc./L;

    xL=xc+w.*qx; yL=yc+w.*qy;
    xR=xc-w.*qx; yR=yc-w.*qy;
    xL(1)=x0+wOpen*tVec(1); yL(1)=y0+wOpen*tVec(2);
    xR(1)=x0-wOpen*tVec(1); yR(1)=y0-wOpen*tVec(2);

    txEnd=dxc(end)/L(end); tyEnd=dyc(end)/L(end);
    qxEnd=-tyEnd; qyEnd=txEnd;
    phi=linspace(0,pi,22);
    xTip=xc(end)+w(end)*(cos(phi).*qxEnd+sin(phi).*txEnd);
    yTip=yc(end)+w(end)*(cos(phi).*qyEnd+sin(phi).*tyEnd);

    polyX=[xL,xTip(2:end-1),fliplr(xR)];
    polyY=[yL,yTip(2:end-1),fliplr(yR)];

    metric=sqrt((aIn*sin(theta0))^2+(bIn*cos(theta0))^2);
    gapHalfTheta=(1.35*wOpen)/metric;

    cr.xc=xc; cr.yc=yc; cr.w=w; cr.wOpen=wOpen;
    cr.polyX=polyX; cr.polyY=polyY; cr.xMem=polyX; cr.yMem=polyY;
    cr.theta0=theta0; cr.gapHalfTheta=gapHalfTheta;
    nNeckField=max(8,round(0.24*numel(xc)));
    cr.neckX=xc(1:nNeckField); cr.neckY=yc(1:nNeckField);
    cr.neckRadius=max(1.15*wBody,2.6*wOpen);
    cr.lumenMask=[];
end

function w = widthProfileOpenToClosed(u,wBody,wOpen,wTip,openFrac,tipFrac)
    w=wBody*ones(size(u));
    idxO=u<=openFrac;
    if any(idxO)
        q=u(idxO)/max(openFrac,eps); q=q.^2.*(3-2*q);
        w(idxO)=wOpen+(wBody-wOpen).*q;
    end
    idxT=u>=(1-tipFrac);
    if any(idxT)
        q=(u(idxT)-(1-tipFrac))/max(tipFrac,eps); q=q.^2.*(3-2*q);
        w(idxT)=wBody+(wTip-wBody).*q;
    end
end

function d = angdiff_local(a,b)
    d=angle(exp(1i*(a-b)));
end

function idxMap = buildMaskIndexMap(mask)
    idxMap=zeros(size(mask));
    idxMap(mask)=1:nnz(mask);
end

function [from,to,q] = buildDirectedEdgesHarmonic(mask,idxMap,Dfield,dx,dy)
    [nRows,nCols]=size(mask);
    from=[]; to=[]; q=[];
    dirs=[0 1;0 -1;1 0;-1 0];
    for r=1:nRows
        for c=1:nCols
            if ~mask(r,c), continue; end
            k=idxMap(r,c);
            for d=1:4
                rr=r+dirs(d,1); cc=c+dirs(d,2);
                if rr<1 || rr>nRows || cc<1 || cc>nCols || ~mask(rr,cc)
                    continue;
                end
                l=idxMap(rr,cc);
                Dk=Dfield(r,c); Dl=Dfield(rr,cc);
                Dij=2*Dk*Dl/max(Dk+Dl,eps);
                if dirs(d,1)==0, h2=dx^2; else, h2=dy^2; end
                from(end+1,1)=k; %#ok<AGROW>
                to(end+1,1)=l; %#ok<AGROW>
                q(end+1,1)=Dij/h2; %#ok<AGROW>
            end
        end
    end
end

function out = momentEulerAll(c1,c2,from1,to1,q1,from2,to2,q2, ...
    x1_0,x2_0,tvec,dt)

    K=numel(c1); nT=numel(tvec);
    m1=x1_0(:); m2=x2_0(:); m3=zeros(K,1);
    out.globalM3=zeros(nT,1);
    out.globalM3(1)=0;
    t=0;

    for it=2:nT
        tTarget=tvec(it);
        while t<tTarget-1e-12
            h=min(dt,tTarget-t);
            dm1=-c1.*m1+c2.*m2;
            dm2= c1.*m1-c2.*m2;
            dm3= c1.*m1;

            flux1=q1.*m1(from1);
            dm1=dm1+accumarray(to1,flux1,[K,1])-accumarray(from1,flux1,[K,1]);
            flux2=q2.*m2(from2);
            dm2=dm2+accumarray(to2,flux2,[K,1])-accumarray(from2,flux2,[K,1]);

            m1=max(m1+h*dm1,0);
            m2=max(m2+h*dm2,0);
            m3=max(m3+h*dm3,0);
            t=t+h;
        end
        out.globalM3(it)=sum(m3);
    end
    out.finalM1=m1; out.finalM2=m2; out.finalM3=m3;
end

function out = buildOutgoingLists(K,from1,to1,q1,from2,to2,q2)
    out.to1=cell(K,1); out.q1=cell(K,1);
    out.to2=cell(K,1); out.q2=cell(K,1);
    for e=1:numel(q1)
        k=from1(e);
        out.to1{k}(end+1,1)=to1(e); %#ok<AGROW>
        out.q1{k}(end+1,1)=q1(e); %#ok<AGROW>
    end
    for e=1:numel(q2)
        k=from2(e);
        out.to2{k}(end+1,1)=to2(e); %#ok<AGROW>
        out.q2{k}(end+1,1)=q2(e); %#ok<AGROW>
    end
    out.q1sum=zeros(K,1); out.q2sum=zeros(K,1);
    for k=1:K
        out.q1sum(k)=sum(out.q1{k});
        out.q2sum(k)=sum(out.q2{k});
    end
end

function [voxelRate,tree] = refreshVoxelRates(touched,voxelRate,tree, ...
    x1,x2,c1,c2,q1sum,q2sum)
    for ii=1:numel(touched)
        k=touched(ii);
        newRate=max(x1(k)*(c1(k)+q1sum(k))+x2(k)*(c2(k)+q2sum(k)),0);
        delta=newRate-voxelRate(k);
        voxelRate(k)=newRate;
        tree=fenwickAdd(tree,k,delta);
    end
end

function j = chooseDestination(to,q)
    if isempty(to), error('No outgoing destination available.'); end
    if numel(to)==1, j=to(1); return; end
    cs=cumsum(q);
    idx=find(cs>=rand*cs(end),1,'first');
    j=to(idx);
end

function tree = fenwickBuild(values)
    n=numel(values); tree=zeros(n,1);
    for i=1:n
        tree=fenwickAdd(tree,i,values(i));
    end
end

function tree = fenwickAdd(tree,idx,delta)
    n=numel(tree);
    while idx<=n
        tree(idx)=tree(idx)+delta;
        lowbit=idx-bitand(idx,idx-1);
        idx=idx+lowbit;
    end
end

function s = fenwickTotal(tree)
    idx=numel(tree); s=0;
    while idx>0
        s=s+tree(idx);
        lowbit=idx-bitand(idx,idx-1);
        idx=idx-lowbit;
    end
end

function idx = fenwickFind(tree,target)
    n=numel(tree);
    if target<=0, target=realmin; end
    idx=0;
    bit=2^floor(log2(n));
    while bit>0
        nxt=idx+bit;
        if nxt<=n && tree(nxt)<target
            idx=nxt;
            target=target-tree(nxt);
        end
        bit=floor(bit/2);
    end
    idx=min(idx+1,n);
end

function plotMeanSEM(ax,t,mat,lineColor,alphaVal)
    mu=mean(mat,1);
    se=std(mat,0,1)/sqrt(size(mat,1));
    patch(ax,[t fliplr(t)],[mu-se fliplr(mu+se)],lineColor, ...
        'FaceAlpha',alphaVal,'EdgeColor','none','HandleVisibility','off');
    plot(ax,t,mu,'-','Color',lineColor,'LineWidth',1.35);
end

function formatAxes(ax)
    set(ax,'Box','on','LineWidth',0.70,'TickDir','out','Layer','top', ...
        'FontSize',7,'XColor',[0.08 0.08 0.08], ...
        'YColor',[0.08 0.08 0.08]);
end

function q = quantileLocal(x,p)
    x=sort(x(:)); n=numel(x);
    if n==0, q=NaN; return; end
    p=min(max(p,0),1);
    pos=1+p*(n-1); lo=floor(pos); hi=ceil(pos);
    if lo==hi, q=x(lo); else, q=x(lo)+(pos-lo)*(x(hi)-x(lo)); end
end

function Y = normalizeInside(Y,mask)
    vals=Y(mask & ~isnan(Y));
    ymin=min(vals); ymax=max(vals);
    Y=(Y-ymin)/(ymax-ymin+eps);
    Y(~mask)=NaN;
end

function cmap = cmapBlueWhiteRed(n)
    anchors=[
        0.10 0.22 0.55
        0.38 0.65 0.84
        0.98 0.98 0.97
        0.94 0.56 0.34
        0.66 0.08 0.06];
    cmap=interpCmap(anchors,n);
end

function cmap = interpCmap(anchors,n)
    x0=linspace(0,1,size(anchors,1));
    xq=linspace(0,1,n);
    cmap=zeros(n,3);
    for c=1:3
        cmap(:,c)=interp1(x0,anchors(:,c),xq,'pchip');
    end
    cmap=max(min(cmap,1),0);
end
