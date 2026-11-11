%% Fig3BC_predictorIndependent_localBurden_FPT.m
% -------------------------------------------------------------------------
% Predictor-independent Fig. 3B--C for mitochondrial SD-RDME.
%
% Scientific purpose
%   Step 2 / Fig. 3B:
%     Demonstrate that a globally burden-matched SD-RDME produces spatially
%     heterogeneous local ROS burden, without defining any region from
%     c1-hat, B_D, X_HB, or any other candidate predictor.
%
%   Step 3 / Fig. 3C:
%     Demonstrate that the predictor-independent local burden redistribution
%     translates into spatially heterogeneous local first-passage times
%     (FPTs), again without defining a predictor-based region.
%
% Core design
%   1) Full SD-RDME and uniform RDME are matched by the ensemble-mean final
%      global ROS burden using exact first-moment equations.
%   2) Locality is defined only by fixed-radius windows W_k centered at every
%      voxel k. The same windows are used for burden and FPT analyses.
%   3) The local FPT threshold is anchored to the globally matched Uniform
%      control:
%          m* = E_U[M3_global(Tobs)] / K,
%          M_Wk = round(m* |W_k|).
%      Thus no c1-hat-, B_D-, X_HB-, or output-selected region enters the
%      threshold definition.
%   4) Exact continuous-time next-subvolume SSA with Fenwick-tree event
%      selection is used. No tau-leaping is used.
%
% Model events
%   R1: Fe2+ -> Fe3+ + cumulative ROS burden, propensity c1_k X1_k
%   R2: Fe3+ -> Fe2+, propensity c2_k X2_k
%   D1: Fe2+ diffusion k -> l, propensity lambda1_kl X1_k
%   D2: Fe3+ diffusion k -> l, propensity lambda2_kl X2_k
%
% Outputs
%   Fig3B_predictorIndependent_localBurden.{pdf,png,tif,fig}
%   Fig3C_predictorIndependent_localFPT.{pdf,png,tif,fig}
%   Fig3BC_predictorIndependent_data.mat
%
% The geometry and kinetic parameters are inherited from
% SDRDME_localburden4(3).m / the current Fig. 3A--C configuration.
% -------------------------------------------------------------------------

clear; close all; clc;
clear functions;

%% --------------------------- Run configuration --------------------------
RUN_MODE = 'final';           % 'smoke' | 'pilot' | 'final'
USE_PARALLEL = true;
AUTO_START_POOL = true;
N_WORKERS = 8;
SAVE_FIG = true;
SAVE_DATA = true;

switch lower(RUN_MODE)
    case 'smoke'
        nRepBurden = 20;
        nRepFPT    = 20;
        nBoot      = 200;
    case 'pilot'
        nRepBurden = 100;
        nRepFPT    = 120;
        nBoot      = 500;
    case 'final'
        nRepBurden = 300;
        nRepFPT    = 400;
        nBoot      = 1000;
    otherwise
        error('Unknown RUN_MODE: %s', RUN_MODE);
end

seedBase = 20260724;
outPrefixB = 'Fig3B_predictorIndependent_localBurden';
outPrefixC = 'Fig3C_predictorIndependent_localFPT';
dataFile   = 'Fig3BC_predictorIndependent_data.mat';

%% --------------------------- Figure style -------------------------------
set(groot,'defaultFigureColor','w');
set(groot,'defaultAxesFontName','Arial');
set(groot,'defaultTextFontName','Arial');
set(groot,'defaultAxesFontSize',7);
set(groot,'defaultTextFontSize',7);
set(groot,'defaultAxesLineWidth',0.70);
set(groot,'defaultLineLineWidth',1.05);
set(groot,'defaultAxesTickDir','out');
set(groot,'defaultAxesBox','on');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','tex');

%% ---------------------- Parameters inherited from prior code ------------
% Grid and crista-inspired geometry.
nx = 60;
ny = 36;
aOut = 1.12;  bOut = 0.62;
aIn  = 0.94;  bIn  = 0.49;
ellG = 0.085;
ellS = 0.042;
connMatrix = 1.00;
connLumen  = 0.82;
connNeck   = 0.18;

% SD-RDME reaction and diffusion fields, min-based units.
D0_Fe2   = 3.0e-4;
D2_scale = 0.70;
D0_Fe3   = D2_scale * D0_Fe2;
c1_0     = 6.0e-3;
etaS     = 3.40;
c2_0     = 1.8e-2;

% Observation window and initial iron state.
Tobs     = 60;          % min
dtRecord = 0.50;        % min
tvec     = 0:dtRecord:Tobs;
FeTot_per_voxel = 6;
fracFe2_0       = 0.65;

% Exact first-moment global matching.
dtMoment = 0.05;
nBisect  = 14;

% Predictor-independent local windows.
windowRadius = 0.10;

% Predictor-independent control-anchored FPT threshold density.
thresholdFactor = 1.00;

% Statistical summaries.
alphaCI = 0.05;
nBinsBurdenFPT = 8;

%% --------------------------- Parallel pool ------------------------------
usePar = configureParallel(USE_PARALLEL, AUTO_START_POOL, N_WORKERS);

%% ---------------------- Build current SD-RDME fields --------------------
geom = buildCurrentCristaGeometry(nx,ny,aOut,bOut,aIn,bIn, ...
    ellG,ellS,connMatrix,connLumen,connNeck);

mask = geom.OmegaR;
x    = geom.x;
y    = geom.y;
Xg   = geom.X;
Yg   = geom.Y;
K    = nnz(mask);

S = geom.S;
G = geom.G;

D_Fe2_sd = D0_Fe2 * G;
D_Fe3_sd = D0_Fe3 * G;
D_Fe2_sd(~mask) = NaN;
D_Fe3_sd(~mask) = NaN;

c1_sd = c1_0 * (1 + etaS * S);
c2_sd = c2_0 * ones(size(c1_sd));
c1_sd(~mask) = NaN;
c2_sd(~mask) = NaN;

idxMap = buildMaskIndexMap(mask);
c1_sd_vec = c1_sd(mask);
c2_sd_vec = c2_sd(mask);
D1_sd_vec = D_Fe2_sd(mask);
D2_sd_vec = D_Fe3_sd(mask);

N0 = FeTot_per_voxel * ones(K,1);
x1_0 = round(fracFe2_0 * N0);
x2_0 = N0 - x1_0;

fprintf('============================================================\n');
fprintf('Predictor-independent SD-RDME burden/FPT analysis (%s)\n', upper(RUN_MODE));
fprintf('Grid                         = %d x %d\n', nx, ny);
fprintf('Omega_R voxels K             = %d\n', K);
fprintf('Tobs                         = %.1f min\n', Tobs);
fprintf('Window radius                = %.3f\n', windowRadius);
fprintf('Burden trajectories/model    = %d\n', nRepBurden);
fprintf('FPT trajectories/model       = %d\n', nRepFPT);
fprintf('c1 range                     = [%.5g, %.5g] min^-1\n', ...
    min(c1_sd_vec), max(c1_sd_vec));
fprintf('D_Fe2 range                  = [%.5g, %.5g]\n', ...
    min(D1_sd_vec), max(D1_sd_vec));
fprintf('============================================================\n');

%% ---------------------- Directed diffusion edges ------------------------
dx = mean(diff(x));
dy = mean(diff(y));

[fromSD1,toSD1,qSD1] = buildDirectedEdgesHarmonic(mask,idxMap,D_Fe2_sd,dx,dy);
[fromSD2,toSD2,qSD2] = buildDirectedEdgesHarmonic(mask,idxMap,D_Fe3_sd,dx,dy);

% Uniform RDME uses spatially constant diffusion equal to the spatial mean.
D1_uni = NaN(size(D_Fe2_sd));
D2_uni = NaN(size(D_Fe3_sd));
D1_uni(mask) = mean(D1_sd_vec);
D2_uni(mask) = mean(D2_sd_vec);

[fromUNI1,toUNI1,qUNI1] = buildDirectedEdgesHarmonic(mask,idxMap,D1_uni,dx,dy);
[fromUNI2,toUNI2,qUNI2] = buildDirectedEdgesHarmonic(mask,idxMap,D2_uni,dx,dy);

%% ---------------------- Exact first-moment global matching --------------
fprintf('\nSolving exact first-moment equations for the SD-RDME target...\n');
SDmom = momentEulerAll(c1_sd_vec,c2_sd_vec,fromSD1,toSD1,qSD1, ...
    fromSD2,toSD2,qSD2,x1_0,x2_0,tvec,dtMoment);
targetGlobal = SDmom.globalM3(end);

c1_sd_mean = mean(c1_sd_vec);
lo = 0.15;
hi = 5.50;

fprintf('Matching Uniform RDME c1 by bisection...\n');
for it = 1:nBisect
    alpha = 0.5 * (lo + hi);
    c1_trial = alpha * c1_sd_mean * ones(K,1);
    trialMom = momentEulerAll(c1_trial,c2_0*ones(K,1), ...
        fromUNI1,toUNI1,qUNI1,fromUNI2,toUNI2,qUNI2, ...
        x1_0,x2_0,tvec,dtMoment);

    if trialMom.globalM3(end) < targetGlobal
        lo = alpha;
    else
        hi = alpha;
    end
end

alphaMatch = 0.5 * (lo + hi);
c1_uni_vec = alphaMatch * c1_sd_mean * ones(K,1);
c2_uni_vec = c2_0 * ones(K,1);

UNImom = momentEulerAll(c1_uni_vec,c2_uni_vec, ...
    fromUNI1,toUNI1,qUNI1,fromUNI2,toUNI2,qUNI2, ...
    x1_0,x2_0,tvec,dtMoment);

momentMismatch = abs(UNImom.globalM3(end)-SDmom.globalM3(end)) / ...
    max(SDmom.globalM3(end),eps);

fprintf('alphaMatch                   = %.6f\n', alphaMatch);
fprintf('Moment global mismatch       = %.5f %%\n', 100*momentMismatch);

%% ---------------------- Predictor-independent local windows -------------
coordX = Xg(mask);
coordY = Yg(mask);
win = buildCircularWindowsSparse(coordX,coordY,windowRadius);

fprintf('Window size range            = [%d, %d] voxels\n', ...
    min(win.size), max(win.size));
fprintf('Median window size           = %.1f voxels\n', median(win.size));

%% ========================================================================
% STEP 2 / FIG. 3B: predictor-independent local ROS burden heterogeneity
% ========================================================================

outSD = buildOutgoingLists(K,fromSD1,toSD1,qSD1,fromSD2,toSD2,qSD2);
outUNI = buildOutgoingLists(K,fromUNI1,toUNI1,qUNI1,fromUNI2,toUNI2,qUNI2);

fprintf('\nRunning exact burden SSA ensemble for SD-RDME...\n');
burdenSD = runBurdenEnsembleNSM(c1_sd_vec,c2_sd_vec,outSD, ...
    x1_0,x2_0,tvec,nRepBurden,seedBase+1000,usePar);

fprintf('Running exact burden SSA ensemble for global-matched Uniform RDME...\n');
burdenUNI = runBurdenEnsembleNSM(c1_uni_vec,c2_uni_vec,outUNI, ...
    x1_0,x2_0,tvec,nRepBurden,seedBase+2000,usePar);

% Global burden per voxel.
gSD  = burdenSD.globalMat  ./ K;
gUNI = burdenUNI.globalMat ./ K;

% Fixed-radius window burden density for each trajectory and window.
windowSumSD  = burdenSD.finalVoxelMat  * win.W.';
windowSumUNI = burdenUNI.finalVoxelMat * win.W.';
windowDensitySD  = windowSumSD  ./ win.size.';
windowDensityUNI = windowSumUNI ./ win.size.';

meanWinSD  = mean(windowDensitySD,1).';
meanWinUNI = mean(windowDensityUNI,1).';
deltaWinBurden = meanWinSD - meanWinUNI;

meanWinSD_map = vectorToMaskMap(meanWinSD,mask);
meanWinUNI_map = vectorToMaskMap(meanWinUNI,mask);
deltaWinBurden_map = vectorToMaskMap(deltaWinBurden,mask);

% Exact first-moment window fields; retained as an audit and for Fig. 3C.
momWinSD  = (win.W * SDmom.finalM3)  ./ win.size;
momWinUNI = (win.W * UNImom.finalM3) ./ win.size;
deltaWinBurdenMoment = momWinSD - momWinUNI;

% Audits and predictor-independent spatial heterogeneity summaries.
globalMismatchSSA = abs(mean(gSD(:,end))-mean(gUNI(:,end))) / ...
    max(mean(gSD(:,end)),eps);
spatialSD_UNI = std(meanWinUNI);
spatialSD_SD  = std(meanWinSD);
spatialSDRatio = spatialSD_SD / max(spatialSD_UNI,eps);
spread_UNI = quantileLocal(meanWinUNI,0.90)-quantileLocal(meanWinUNI,0.10);
spread_SD  = quantileLocal(meanWinSD,0.90)-quantileLocal(meanWinSD,0.10);
spreadRatio = spread_SD / max(spread_UNI,eps);
fractionPositiveBurden = mean(deltaWinBurden > 0);

momentVsSSA_R = pearsonLocal(deltaWinBurdenMoment,deltaWinBurden);

fprintf('\n================ Step 2 / local burden summary ================\n');
fprintf('Final global burden per voxel, Uniform = %.4f\n', mean(gUNI(:,end)));
fprintf('Final global burden per voxel, SD      = %.4f\n', mean(gSD(:,end)));
fprintf('SSA final global mismatch              = %.3f %%\n', 100*globalMismatchSSA);
fprintf('Spatial SD of local burden, Uniform    = %.4f\n', spatialSD_UNI);
fprintf('Spatial SD of local burden, SD         = %.4f\n', spatialSD_SD);
fprintf('Spatial-SD ratio SD/Uniform            = %.3f\n', spatialSDRatio);
fprintf('P90-P10 spread ratio SD/Uniform        = %.3f\n', spreadRatio);
fprintf('Fraction of windows with Delta burden>0= %.3f\n', fractionPositiveBurden);
fprintf('Moment-vs-SSA Delta-burden correlation = %.4f\n', momentVsSSA_R);
fprintf('===============================================================\n');

%% ---------------------- Plot predictor-independent Fig. 3B -------------
colSD   = [0.76 0.19 0.13];
colUNI  = [0.38 0.40 0.43];
colText = [0.08 0.08 0.08];

figB = figure('Units','centimeters','Position',[1.2 2.0 18.3 5.3],'Color','w');
tlB = tiledlayout(figB,1,4,'TileSpacing','compact','Padding','compact');

% B1. Global matching.
ax = nexttile(tlB,1); hold(ax,'on');
plotMeanSEM(ax,tvec,gUNI,colUNI,0.12);
plotMeanSEM(ax,tvec,gSD,colSD,0.12);
text(ax,0.02,0.98,'B','Units','normalized','FontSize',12, ...
    'FontWeight','bold','VerticalAlignment','top','Color',colText);
title(ax,'Global burden matching','FontSize',8.2,'FontWeight','normal');
xlabel(ax,'Time (min)');
ylabel(ax,'ROS burden per voxel');
legend(ax,{'Uniform RDME','SD-RDME'},'Box','off','Location','northwest', ...
    'FontSize',6.2,'Interpreter','tex');
text(ax,0.97,0.08,sprintf('mismatch = %.2f%%',100*globalMismatchSSA), ...
    'Units','normalized','HorizontalAlignment','right','FontSize',6.2, ...
    'Color',colText);
xlim(ax,[0 Tobs]);
formatAxes(ax);

% Shared sequential color scale for B2--B3.
allBurdenVals = [meanWinUNI; meanWinSD];
seqMax = max(quantileLocal(allBurdenVals,0.99),eps);
seqMin = max(0,quantileLocal(allBurdenVals,0.01));
if seqMax <= seqMin
    seqMin = min(allBurdenVals);
    seqMax = max(allBurdenVals) + eps;
end

% B2. Uniform local burden field.
ax = nexttile(tlB,2); hold(ax,'on');
plotMaskedMap(ax,x,y,meanWinUNI_map,mask,cmapSequentialBlue(256),[seqMin seqMax]);
contour(ax,x,y,double(mask),[0.5 0.5],'Color',[0.18 0.18 0.18], ...
    'LineWidth',0.55);
cb = colorbar(ax); styleColorbar(cb,'Local burden density');
title(ax,'Uniform RDME','FontSize',8.2,'FontWeight','normal');
formatMapAxes(ax);

% B3. SD-RDME local burden field.
ax = nexttile(tlB,3); hold(ax,'on');
plotMaskedMap(ax,x,y,meanWinSD_map,mask,cmapSequentialBlue(256),[seqMin seqMax]);
contour(ax,x,y,double(mask),[0.5 0.5],'Color',[0.18 0.18 0.18], ...
    'LineWidth',0.55);
cb = colorbar(ax); styleColorbar(cb,'Local burden density');
title(ax,'SD-RDME','FontSize',8.2,'FontWeight','normal');
formatMapAxes(ax);

% B4. Predictor-independent difference field.
ax = nexttile(tlB,4); hold(ax,'on');
divMax = max(quantileLocal(abs(deltaWinBurden),0.98),eps);
plotMaskedMap(ax,x,y,deltaWinBurden_map,mask,cmapBlueWhiteRed(256),[-divMax divMax]);
contour(ax,x,y,double(mask),[0.5 0.5],'Color',[0.18 0.18 0.18], ...
    'LineWidth',0.55);
cb = colorbar(ax); styleColorbar(cb,'$\Delta$ local burden','latex');
title(ax,'Local burden redistribution','FontSize',8.2,'FontWeight','normal');
text(ax,0.04,0.96,sprintf('spatial SD ratio = %.2f',spatialSDRatio), ...
    'Units','normalized','VerticalAlignment','top','FontSize',6.1, ...
    'Color',colText,'BackgroundColor','w','Margin',1.0);
text(ax,0.04,0.84,sprintf('P90-P10 ratio = %.2f',spreadRatio), ...
    'Units','normalized','VerticalAlignment','top','FontSize',6.1, ...
    'Color',colText,'BackgroundColor','w','Margin',1.0);
formatMapAxes(ax);

if SAVE_FIG
    exportgraphics(figB,[outPrefixB '.pdf'],'ContentType','vector');
    exportgraphics(figB,[outPrefixB '.png'],'Resolution',600);
    exportgraphics(figB,[outPrefixB '.tif'],'Resolution',600);
    savefig(figB,[outPrefixB '.fig']);
end

%% ========================================================================
% STEP 3 / FIG. 3C: predictor-independent local FPT heterogeneity
% ========================================================================

% Control-anchored threshold density based only on the globally matched
% Uniform RDME mean burden density at Tobs.
mStar = thresholdFactor * UNImom.globalM3(end) / K;
Mwindow = max(1,round(mStar * win.size));

fprintf('\nPredictor-independent FPT threshold definition:\n');
fprintf('m* = E_U[M3_global(Tobs)]/K = %.6f counts/voxel\n',mStar);
fprintf('Window threshold range            = [%d, %d] counts\n', ...
    min(Mwindow),max(Mwindow));

fprintf('\nRunning exact all-window FPT ensemble for SD-RDME...\n');
FPT_SD = runWindowFPTensembleNSM(c1_sd_vec,c2_sd_vec,outSD,win, ...
    x1_0,x2_0,Mwindow,Tobs,nRepFPT,seedBase+3000,usePar);

fprintf('Running exact all-window FPT ensemble for global-matched Uniform RDME...\n');
FPT_UNI = runWindowFPTensembleNSM(c1_uni_vec,c2_uni_vec,outUNI,win, ...
    x1_0,x2_0,Mwindow,Tobs,nRepFPT,seedBase+4000,usePar);

RMST_SD  = mean(min(FPT_SD,Tobs),1).';
RMST_UNI = mean(min(FPT_UNI,Tobs),1).';
DeltaRMST = RMST_UNI - RMST_SD;   % positive means earlier under SD-RDME

Pcross_SD  = mean(FPT_SD <= Tobs,1).';
Pcross_UNI = mean(FPT_UNI <= Tobs,1).';
DeltaPcross = Pcross_SD - Pcross_UNI;

DeltaRMST_map = vectorToMaskMap(DeltaRMST,mask);
DeltaPcross_map = vectorToMaskMap(DeltaPcross,mask);

% Predictor-independent burden-to-FPT relation.
[rhoBurdenFPT,rankBurden,rankFPT] = spearmanLocal(deltaWinBurdenMoment,DeltaRMST); %#ok<ASGLU>
[binX,binMed,binLo,binHi] = binnedQuantileSummary( ...
    deltaWinBurdenMoment,DeltaRMST,nBinsBurdenFPT);

fractionEarlier = mean(DeltaRMST > 0);
medianDeltaRMST = median(DeltaRMST);
meanDeltaRMST = mean(DeltaRMST);
q90DeltaRMST = quantileLocal(DeltaRMST,0.90);

% Trajectory-level mean local RMST provides a Monte-Carlo CI without
% treating overlapping windows as independent replicates.
trajMeanRMST_SD  = mean(min(FPT_SD,Tobs),2);
trajMeanRMST_UNI = mean(min(FPT_UNI,Tobs),2);
[meanDeltaRMSTTraj,ciDeltaRMST] = bootstrapDifferenceIndependent( ...
    trajMeanRMST_UNI,trajMeanRMST_SD,nBoot,alphaCI,seedBase+5000);

fprintf('\n================ Step 3 / local FPT summary ==================\n');
fprintf('Mean window DeltaRMST                = %.4f min\n',meanDeltaRMST);
fprintf('Median window DeltaRMST              = %.4f min\n',medianDeltaRMST);
fprintf('90th percentile DeltaRMST            = %.4f min\n',q90DeltaRMST);
fprintf('Fraction of windows with earlier FPT = %.3f\n',fractionEarlier);
fprintf('Trajectory-averaged DeltaRMST        = %.4f min\n',meanDeltaRMSTTraj);
fprintf('Bootstrap %.0f%% CI                   = [%.4f, %.4f] min\n', ...
    100*(1-alphaCI),ciDeltaRMST(1),ciDeltaRMST(2));
fprintf('Spearman rho(Delta burden,DeltaRMST) = %.4f\n',rhoBurdenFPT);
fprintf('===============================================================\n');

%% ---------------------- Plot predictor-independent Fig. 3C -------------
figC = figure('Units','centimeters','Position',[1.2 2.0 18.3 5.5],'Color','w');
tlC = tiledlayout(figC,1,3,'TileSpacing','compact','Padding','compact');

% C1. Spatial map of FPT advancement.
ax = nexttile(tlC,1); hold(ax,'on');
rmstMax = max(quantileLocal(abs(DeltaRMST),0.98),eps);
plotMaskedMap(ax,x,y,DeltaRMST_map,mask,cmapBlueWhiteRed(256),[-rmstMax rmstMax]);
contour(ax,x,y,double(mask),[0.5 0.5],'Color',[0.18 0.18 0.18], ...
    'LineWidth',0.55);
cb = colorbar(ax); styleColorbar(cb,'$\Delta\mathrm{RMST}_{60}$ (min)','latex');
text(ax,0.02,0.98,'C','Units','normalized','FontSize',12, ...
    'FontWeight','bold','VerticalAlignment','top','Color',colText);
title(ax,'Observed local FPT advancement','FontSize',8.2,'FontWeight','normal');
text(ax,0.04,0.95,sprintf('mean = %.2f min',meanDeltaRMST), ...
    'Units','normalized','VerticalAlignment','top','FontSize',6.1, ...
    'BackgroundColor','w','Margin',1.0,'Color',colText);
formatMapAxes(ax);

% C2. Distribution across all local windows.
ax = nexttile(tlC,2); hold(ax,'on');
[xECDF,FECDF] = simpleECDF(DeltaRMST);
plot(ax,xECDF,FECDF,'-','Color',colSD,'LineWidth',1.35);
plot(ax,[0 0],[0 1],'--','Color',[0.50 0.50 0.50], ...
    'LineWidth',0.75,'HandleVisibility','off');
title(ax,'Distribution across local windows','FontSize',8.2,'FontWeight','normal');
xlabel(ax,'$\Delta\mathrm{RMST}_{60}$ (min)','Interpreter','latex');
ylabel(ax,'Cumulative fraction');
text(ax,0.05,0.90,sprintf('fraction > 0 = %.2f',fractionEarlier), ...
    'Units','normalized','FontSize',6.3,'Color',colText);
text(ax,0.05,0.81,sprintf('median = %.2f min',medianDeltaRMST), ...
    'Units','normalized','FontSize',6.3,'Color',colText);
formatAxes(ax);
ylim(ax,[0 1]);

% C3. Local burden redistribution versus FPT advancement.
ax = nexttile(tlC,3); hold(ax,'on');
scatter(ax,deltaWinBurdenMoment,DeltaRMST,8,[0.62 0.64 0.67], ...
    'filled','MarkerFaceAlpha',0.42,'MarkerEdgeAlpha',0.20);
errorbar(ax,binX,binMed,binMed-binLo,binHi-binMed,'o-', ...
    'Color',colSD,'MarkerFaceColor','w','MarkerSize',3.8, ...
    'LineWidth',1.10,'CapSize',2.5);
plot(ax,xlim(ax),[0 0],'--','Color',[0.72 0.72 0.72], ...
    'LineWidth',0.65,'HandleVisibility','off');
plot(ax,[0 0],ylim(ax),'--','Color',[0.72 0.72 0.72], ...
    'LineWidth',0.65,'HandleVisibility','off');
title(ax,'Local burden--FPT correspondence','FontSize',8.2,'FontWeight','normal');
xlabel(ax,'$\Delta$ local ROS burden density','Interpreter','latex');
ylabel(ax,'$\Delta\mathrm{RMST}_{60}$ (min)','Interpreter','latex');
text(ax,0.05,0.92,sprintf('$\rho_s = %.2f$',rhoBurdenFPT), ...
    'Units','normalized','Interpreter','latex','FontSize',6.7,'Color',colText);
formatAxes(ax);

if SAVE_FIG
    exportgraphics(figC,[outPrefixC '.pdf'],'ContentType','vector');
    exportgraphics(figC,[outPrefixC '.png'],'Resolution',600);
    exportgraphics(figC,[outPrefixC '.tif'],'Resolution',600);
    savefig(figC,[outPrefixC '.fig']);
end

%% ---------------------- Save data ---------------------------------------
if SAVE_DATA
    meta = struct();
    meta.RUN_MODE = RUN_MODE;
    meta.nx = nx; meta.ny = ny; meta.K = K;
    meta.Tobs = Tobs; meta.dtRecord = dtRecord;
    meta.windowRadius = windowRadius;
    meta.nRepBurden = nRepBurden; meta.nRepFPT = nRepFPT;
    meta.alphaMatch = alphaMatch;
    meta.momentMismatch = momentMismatch;
    meta.globalMismatchSSA = globalMismatchSSA;
    meta.thresholdDensity = mStar;
    meta.parameters = struct('D0_Fe2',D0_Fe2,'D0_Fe3',D0_Fe3, ...
        'c1_0',c1_0,'etaS',etaS,'c2_0',c2_0, ...
        'FeTot_per_voxel',FeTot_per_voxel,'fracFe2_0',fracFe2_0);

    save(dataFile,'meta','x','y','mask','S','G','c1_sd','D_Fe2_sd', ...
        'tvec','gSD','gUNI','meanWinSD','meanWinUNI','deltaWinBurden', ...
        'deltaWinBurdenMoment','meanWinSD_map','meanWinUNI_map', ...
        'deltaWinBurden_map','windowDensitySD','windowDensityUNI', ...
        'Mwindow','mStar','FPT_SD','FPT_UNI','RMST_SD','RMST_UNI', ...
        'DeltaRMST','DeltaPcross','DeltaRMST_map','DeltaPcross_map', ...
        'rhoBurdenFPT','meanDeltaRMSTTraj','ciDeltaRMST','-v7.3');
end

fprintf('\nCompleted. Outputs:\n');
fprintf('  %s.*\n',outPrefixB);
fprintf('  %s.*\n',outPrefixC);
if SAVE_DATA
    fprintf('  %s\n',dataFile);
end

%% ============================ Local functions ===========================

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

function win = buildCircularWindowsSparse(coordX,coordY,radius)
    K=numel(coordX);
    I=cell(K,1); J=cell(K,1);
    members=cell(K,1);
    for k=1:K
        idx=find((coordX-coordX(k)).^2+(coordY-coordY(k)).^2 <= radius^2+1e-14);
        members{k}=idx(:);
        I{k}=repmat(k,numel(idx),1);
        J{k}=idx(:);
    end
    I=vertcat(I{:}); J=vertcat(J{:});
    W=sparse(I,J,1,K,K);
    sizeVec=full(sum(W,2));

    centersByVoxel=cell(K,1);
    for j=1:K
        centersByVoxel{j}=find(W(:,j));
    end

    win.W=W;
    win.size=sizeVec;
    win.members=members;
    win.centersByVoxel=centersByVoxel;
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

function out = runBurdenEnsembleNSM(c1,c2,outgoing,x1_0,x2_0,tvec,nRep,seedBase,usePar)
    K=numel(c1); nT=numel(tvec);
    globalMat=zeros(nRep,nT);
    finalVoxelMat=zeros(nRep,K);

    if usePar
        parfor r=1:nRep
            rng(seedBase+r,'twister');
            [g,v]=simulateBurdenOneNSM(c1,c2,outgoing,x1_0,x2_0,tvec);
            globalMat(r,:)=g;
            finalVoxelMat(r,:)=v.';
        end
    else
        for r=1:nRep
            rng(seedBase+r,'twister');
            [g,v]=simulateBurdenOneNSM(c1,c2,outgoing,x1_0,x2_0,tvec);
            globalMat(r,:)=g;
            finalVoxelMat(r,:)=v.';
            if mod(r,max(1,round(nRep/10)))==0
                fprintf('  finished %d / %d trajectories\n',r,nRep);
            end
        end
    end
    out.globalMat=globalMat;
    out.finalVoxelMat=finalVoxelMat;
end

function [globalRec,x3] = simulateBurdenOneNSM(c1,c2,out,x1_0,x2_0,tvec)
    K=numel(c1); nT=numel(tvec);
    x1=x1_0(:); x2=x2_0(:); x3=zeros(K,1);
    voxelRate=x1.*(c1+out.q1sum)+x2.*(c2+out.q2sum);
    tree=fenwickBuild(voxelRate);

    globalRec=zeros(1,nT);
    globalROS=0;
    t=tvec(1); outIdx=1; Tend=tvec(end);

    while t<Tend
        a0=fenwickTotal(tree);
        if ~(a0>0) || ~isfinite(a0), break; end
        tau=-log(rand)/a0;
        tNext=t+tau;

        while outIdx<nT && tvec(outIdx+1)<=min(tNext,Tend)
            outIdx=outIdx+1;
            globalRec(outIdx)=globalROS;
        end
        if tNext>Tend, break; end

        k=fenwickFind(tree,rand*a0);
        aR1=c1(k)*x1(k); aR2=c2(k)*x2(k);
        aD1=out.q1sum(k)*x1(k); aD2=out.q2sum(k)*x2(k);
        aLocal=aR1+aR2+aD1+aD2;
        if aLocal<=0
            old=voxelRate(k); voxelRate(k)=0;
            tree=fenwickAdd(tree,k,-old); t=tNext; continue;
        end

        u=rand*aLocal; touched=k;
        if u<=aR1
            x1(k)=x1(k)-1; x2(k)=x2(k)+1;
            x3(k)=x3(k)+1; globalROS=globalROS+1;
        elseif u<=aR1+aR2
            x2(k)=x2(k)-1; x1(k)=x1(k)+1;
        elseif u<=aR1+aR2+aD1
            j=chooseDestination(out.to1{k},out.q1{k});
            x1(k)=x1(k)-1; x1(j)=x1(j)+1; touched=[k;j]; %#ok<AGROW>
        else
            j=chooseDestination(out.to2{k},out.q2{k});
            x2(k)=x2(k)-1; x2(j)=x2(j)+1; touched=[k;j]; %#ok<AGROW>
        end

        [voxelRate,tree]=refreshVoxelRates(unique(touched),voxelRate,tree, ...
            x1,x2,c1,c2,out.q1sum,out.q2sum);
        t=tNext;
    end

    while outIdx<nT
        outIdx=outIdx+1; globalRec(outIdx)=globalROS;
    end
end

function FPT = runWindowFPTensembleNSM(c1,c2,outgoing,win,x1_0,x2_0, ...
    Mwindow,Tobs,nRep,seedBase,usePar)

    K=numel(c1);
    FPT=inf(nRep,K);
    if usePar
        parfor r=1:nRep
            rng(seedBase+r,'twister');
            FPT(r,:)=simulateWindowFPTOneNSM(c1,c2,outgoing,win, ...
                x1_0,x2_0,Mwindow,Tobs).';
        end
    else
        for r=1:nRep
            rng(seedBase+r,'twister');
            FPT(r,:)=simulateWindowFPTOneNSM(c1,c2,outgoing,win, ...
                x1_0,x2_0,Mwindow,Tobs).';
            if mod(r,max(1,round(nRep/10)))==0
                fprintf('  finished %d / %d FPT trajectories\n',r,nRep);
            end
        end
    end
end

function fpt = simulateWindowFPTOneNSM(c1,c2,out,win,x1_0,x2_0,Mwindow,Tobs)
    K=numel(c1);
    x1=x1_0(:); x2=x2_0(:);
    localROS=zeros(K,1);
    fpt=inf(K,1);
    active=true(K,1);

    voxelRate=x1.*(c1+out.q1sum)+x2.*(c2+out.q2sum);
    tree=fenwickBuild(voxelRate);
    t=0;

    while t<Tobs && any(active)
        a0=fenwickTotal(tree);
        if ~(a0>0) || ~isfinite(a0), break; end
        tau=-log(rand)/a0;
        tNext=t+tau;
        if tNext>Tobs, break; end

        k=fenwickFind(tree,rand*a0);
        aR1=c1(k)*x1(k); aR2=c2(k)*x2(k);
        aD1=out.q1sum(k)*x1(k); aD2=out.q2sum(k)*x2(k);
        aLocal=aR1+aR2+aD1+aD2;
        if aLocal<=0
            old=voxelRate(k); voxelRate(k)=0;
            tree=fenwickAdd(tree,k,-old); t=tNext; continue;
        end

        u=rand*aLocal; touched=k;
        if u<=aR1
            x1(k)=x1(k)-1; x2(k)=x2(k)+1;
            centers=win.centersByVoxel{k};
            if ~isempty(centers)
                centers=centers(active(centers));
                localROS(centers)=localROS(centers)+1;
                crossed=centers(localROS(centers)>=Mwindow(centers));
                if ~isempty(crossed)
                    fpt(crossed)=tNext;
                    active(crossed)=false;
                end
            end
        elseif u<=aR1+aR2
            x2(k)=x2(k)-1; x1(k)=x1(k)+1;
        elseif u<=aR1+aR2+aD1
            j=chooseDestination(out.to1{k},out.q1{k});
            x1(k)=x1(k)-1; x1(j)=x1(j)+1; touched=[k;j]; %#ok<AGROW>
        else
            j=chooseDestination(out.to2{k},out.q2{k});
            x2(k)=x2(k)-1; x2(j)=x2(j)+1; touched=[k;j]; %#ok<AGROW>
        end

        [voxelRate,tree]=refreshVoxelRates(unique(touched),voxelRate,tree, ...
            x1,x2,c1,c2,out.q1sum,out.q2sum);
        t=tNext;
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

function map = vectorToMaskMap(v,mask)
    map=NaN(size(mask)); map(mask)=v;
end

function plotMeanSEM(ax,t,mat,lineColor,alphaVal)
    mu=mean(mat,1);
    se=std(mat,0,1)/sqrt(size(mat,1));
    patch(ax,[t fliplr(t)],[mu-se fliplr(mu+se)],lineColor, ...
        'FaceAlpha',alphaVal,'EdgeColor','none','HandleVisibility','off');
    plot(ax,t,mu,'-','Color',lineColor,'LineWidth',1.35);
end

function plotMaskedMap(ax,x,y,map,mask,cmap,clim)
    h=imagesc(ax,x,y,map);
    set(ax,'YDir','normal');
    set(h,'AlphaData',mask,'Interpolation','nearest');
    axis(ax,'image');
    colormap(ax,cmap); caxis(ax,clim);
end

function styleColorbar(cb,labelText,varargin)
    cb.Box='off'; cb.TickDirection='out'; cb.FontSize=6.0;
    cb.Label.String=labelText; cb.Label.FontSize=6.8;
    if nargin>=3 && strcmpi(varargin{1},'latex')
        cb.Label.Interpreter='latex';
    end
end

function formatMapAxes(ax)
    formatAxes(ax);
    set(ax,'XTick',[-1 0 1],'YTick',[-0.5 0 0.5]);
    xlabel(ax,'x'); ylabel(ax,'y');
end

function formatAxes(ax)
    set(ax,'Box','on','LineWidth',0.70,'TickDir','out','Layer','top', ...
        'FontSize',7,'XColor',[0.08 0.08 0.08], ...
        'YColor',[0.08 0.08 0.08]);
end

function [x,F] = simpleECDF(v)
    x=sort(v(:)); F=(1:numel(x))'/numel(x);
end

function [rho,rx,ry] = spearmanLocal(x,y)
    rx=tiedRankLocal(x(:)); ry=tiedRankLocal(y(:));
    rho=pearsonLocal(rx,ry);
end

function r = tiedRankLocal(x)
    [xs,ord]=sort(x);
    r=zeros(size(x));
    i=1; n=numel(x);
    while i<=n
        j=i;
        while j<n && xs(j+1)==xs(i), j=j+1; end
        r(ord(i:j))=(i+j)/2;
        i=j+1;
    end
end

function r = pearsonLocal(x,y)
    x=x(:)-mean(x(:)); y=y(:)-mean(y(:));
    r=sum(x.*y)/sqrt(max(sum(x.^2)*sum(y.^2),eps));
end

function [binX,binMed,binLo,binHi] = binnedQuantileSummary(x,y,nBins)
    x=x(:); y=y(:);
    [~,ord]=sort(x); n=numel(x);
    edges=round(linspace(0,n,nBins+1));
    binX=NaN(nBins,1); binMed=binX; binLo=binX; binHi=binX;
    for b=1:nBins
        idx=ord(edges(b)+1:edges(b+1));
        binX(b)=median(x(idx));
        binMed(b)=median(y(idx));
        binLo(b)=quantileLocal(y(idx),0.25);
        binHi(b)=quantileLocal(y(idx),0.75);
    end
end

function [delta,ci] = bootstrapDifferenceIndependent(a,b,nBoot,alpha,seed)
    rng(seed,'twister');
    a=a(:); b=b(:);
    delta=mean(a)-mean(b);
    boot=zeros(nBoot,1);
    na=numel(a); nb=numel(b);
    for k=1:nBoot
        boot(k)=mean(a(randi(na,na,1)))-mean(b(randi(nb,nb,1)));
    end
    ci=[quantileLocal(boot,alpha/2),quantileLocal(boot,1-alpha/2)];
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

function cmap = cmapSequentialBlue(n)
    anchors=[
        0.97 0.98 0.98
        0.78 0.88 0.92
        0.45 0.70 0.78
        0.15 0.43 0.58
        0.03 0.18 0.31];
    cmap=interpCmap(anchors,n);
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
