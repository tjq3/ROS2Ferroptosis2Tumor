%% Fig3G_c1hat_heldoutSpatialValidation.m
% -------------------------------------------------------------------------
% Held-out spatial validation of c1-hat after three-candidate screening.
%
% Scientific role
%   Discovery stage (completed elsewhere):
%       Compare c1-hat, B_D and X_HB using the Fig. 3B-C final ensemble.
%
%   Validation stage (this script):
%       Freeze the selected descriptor c1-hat and its top-15% region before
%       generating any new outcome trajectories. Use independent SSA seeds
%       to test whether c1-hat localizes:
%         1) held-out local ROS-burden redistribution;
%         2) held-out local FPT advancement;
%         3) a predicted-zone burden increase and earlier threshold crossing.
%
% Strict anti-circularity design
%   - Only locked model fields/metadata are loaded from the Fig. 3B-C file.
%   - The discovery burden and FPT trajectories are NOT loaded or reused.
%   - The predictor, q=0.85, window radius, global matching coefficient and
%     threshold density are fixed before the held-out simulations.
%   - New non-overlapping random seeds are used for all validation ensembles.
%
% Locked parameters inherited from Fig. 3B-C
%   Grid: 60 x 36; K = 868 for the current geometry
%   Tobs = 60 min; dtRecord = 0.5 min
%   D0_Fe2 = 3.0e-4; D0_Fe3 = 0.70 D0_Fe2
%   c1_0 = 6.0e-3 min^-1; etaS = 3.40; c2_0 = 1.8e-2 min^-1
%   Initial iron per voxel: Fe2+/Fe3+ = 4/2
%   Fixed-radius windows: radius = 0.10
%   Window FPT thresholds and uniform-RDME matching are loaded unchanged.
%
% Main outputs
%   Fig3G_c1hat_heldoutValidation_main.{pdf,png,tif,fig}
%   FigS_c1hat_heldoutValidation_details.{pdf,png,tif,fig}
%   Fig3G_c1hat_heldoutValidation_summary.csv
%   Fig3G_c1hat_heldoutValidation_data.mat
%
% Notes
%   * Window-level burden and FPT maps use the same predictor-independent
%     windows as Fig. 3B-C.
%   * The predicted-zone burden is the mean local burden density across all
%     c1-hat top-15% window centers. It is a linear weighted sum of voxel
%     burden and is updated exactly at every R1 event.
%   * The predicted-zone FPT is the first time this mean burden density
%     reaches the locked control-anchored density threshold mStar.
% -------------------------------------------------------------------------

clear; close all; clc;
clear functions;

%% --------------------------- Configuration ------------------------------
RUN_MODE = 'final';                  % 'smoke' | 'pilot' | 'final'
LOCKED_FILE = 'Fig3BC_predictorIndependent_data.mat';
USE_PARALLEL = true;
AUTO_START_POOL = true;
N_WORKERS = 8;
SAVE_FIG = true;
SAVE_DATA = true;

qMain = 0.85;
qList = 0.75:0.025:0.95;
alphaCI = 0.05;

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
        error('Unknown RUN_MODE: %s',RUN_MODE);
end

% Completely separate from the Fig. 3B-C and candidate-screening seeds.
seedBase = 20260728;

outMain = 'Fig3G_c1hat_heldoutValidation_main';
outSI   = 'FigS_c1hat_heldoutValidation_details';
outCSV  = 'Fig3G_c1hat_heldoutValidation_summary.csv';
outMAT  = 'Fig3G_c1hat_heldoutValidation_data.mat';

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

colSD   = [0.76 0.19 0.13];
colUNI  = [0.38 0.40 0.43];
colPred = [0.88 0.36 0.10];
colBurden = [0.12 0.12 0.12];
colFPT    = [0.00 0.53 0.53];
colText   = [0.08 0.08 0.08];

%% ---------------------- Load locked model state -------------------------
if ~isfile(LOCKED_FILE)
    error(['Locked Fig. 3B-C data file not found: %s\n', ...
        'Run the finalized predictor-independent Fig. 3B-C script first.'], ...
        LOCKED_FILE);
end

% Deliberately do not load discovery outcome trajectories.
L = load(LOCKED_FILE,'meta','x','y','mask','S','G','c1_sd', ...
    'D_Fe2_sd','mStar','Mwindow');

required = {'meta','x','y','mask','S','G','c1_sd','D_Fe2_sd', ...
    'mStar','Mwindow'};
for i = 1:numel(required)
    if ~isfield(L,required{i})
        error('Locked data file is missing required variable: %s',required{i});
    end
end

meta = L.meta;
x = L.x(:).';
y = L.y(:).';
mask = logical(L.mask);
S = L.S;
G = L.G;
c1_sd = L.c1_sd;
D_Fe2_sd = L.D_Fe2_sd;
mStar = L.mStar;
MwindowLocked = L.Mwindow(:);

K = nnz(mask);
[nRows,nCols] = size(mask);
assert(numel(x)==nCols && numel(y)==nRows, ...
    'Coordinate vectors do not match the saved mask.');
[Xg,Yg] = meshgrid(x,y);

% Lock all kinetic and observation parameters to Fig. 3B-C metadata.
Tobs = meta.Tobs;
dtRecord = meta.dtRecord;
tvec = 0:dtRecord:Tobs;
windowRadius = meta.windowRadius;
alphaMatch = meta.alphaMatch;

p = meta.parameters;
D0_Fe2 = p.D0_Fe2;
D0_Fe3 = p.D0_Fe3;
c1_0 = p.c1_0;
etaS = p.etaS;
c2_0 = p.c2_0;
FeTot_per_voxel = p.FeTot_per_voxel;
fracFe2_0 = p.fracFe2_0;

% The actual rounded state inherited from Fig. 3B-C is 4/2 per voxel.
N0 = FeTot_per_voxel*ones(K,1);
x1_0 = round(fracFe2_0*N0);
x2_0 = N0-x1_0;

% Reconstruct the fields that were not saved separately.
D_Fe3_sd = (D0_Fe3/D0_Fe2)*D_Fe2_sd;
D_Fe3_sd(~mask) = NaN;
c2_sd = c2_0*ones(size(c1_sd));
c2_sd(~mask) = NaN;

idxMap = buildMaskIndexMap(mask);
c1_sd_vec = c1_sd(mask);
c2_sd_vec = c2_sd(mask);
D1_sd_vec = D_Fe2_sd(mask);
D2_sd_vec = D_Fe3_sd(mask);

fprintf('============================================================\n');
fprintf('Held-out c1-hat spatial validation (%s)\n',upper(RUN_MODE));
fprintf('Locked source file             = %s\n',LOCKED_FILE);
fprintf('Grid                           = %d x %d\n',numel(x),numel(y));
fprintf('Omega_R voxels K               = %d\n',K);
fprintf('Tobs / dtRecord                = %.1f / %.2f min\n',Tobs,dtRecord);
fprintf('Window radius                  = %.3f\n',windowRadius);
fprintf('qMain / predicted fraction     = %.2f / %.2f\n',qMain,1-qMain);
fprintf('Held-out burden trajectories   = %d per model\n',nRepBurden);
fprintf('Held-out FPT trajectories      = %d per model\n',nRepFPT);
fprintf('Bootstrap replicates           = %d\n',nBoot);
fprintf('Initial Fe2+/Fe3+ per voxel    = %d / %d\n',x1_0(1),x2_0(1));
fprintf('============================================================\n');

%% ---------------------- Build locked local windows ----------------------
coordX = Xg(mask);
coordY = Yg(mask);
win = buildCircularWindowsSparse(coordX,coordY,windowRadius);

assert(numel(MwindowLocked)==K, ...
    'Locked Mwindow length does not match the reconstructed windows.');
MwindowExpected = max(1,round(mStar*win.size));
if any(MwindowExpected~=MwindowLocked)
    error(['Reconstructed window thresholds differ from the locked Fig. 3B-C ', ...
        'thresholds. Check mask ordering or windowRadius.']);
end
Mwindow = MwindowLocked;

fprintf('Window-size range              = [%d, %d]\n', ...
    min(win.size),max(win.size));
fprintf('Median window size             = %.1f\n',median(win.size));
fprintf('Locked threshold density mStar = %.6f counts/voxel\n',mStar);
fprintf('Window-threshold range         = [%d, %d] counts\n', ...
    min(Mwindow),max(Mwindow));

%% ---------------------- Freeze the c1-hat prediction --------------------
% Use exactly the score definition selected in the screening analysis:
% center c1-hat normalized first, then averaged over the same local window.
c1Center = normalizeVector(c1_sd_vec);
c1Window = (win.W*c1Center)./win.size;
c1Window = normalizeVector(c1Window);

nPred = max(1,ceil((1-qMain)*K));
predMaskVec = topNMask(c1Window,nPred);
restMaskVec = ~predMaskVec;

c1WindowMap = vectorToMaskMap(c1Window,mask);
predMaskMap = false(size(mask));
predMaskMap(mask) = predMaskVec;

% Linear weights for the average local burden density across predicted or
% remaining window centers.  Each weight vector sums to one.
wPred = buildAverageWindowWeight(win,predMaskVec);
wRest = buildAverageWindowWeight(win,restMaskVec);

fprintf('Predicted top-q windows        = %d / %d (%.3f)\n', ...
    nnz(predMaskVec),K,nnz(predMaskVec)/K);
fprintf('sum(wPred) / sum(wRest)        = %.6f / %.6f\n', ...
    sum(wPred),sum(wRest));

%% ---------------------- Build SD and matched uniform RDME ---------------
dx = mean(diff(x));
dy = mean(diff(y));

[fromSD1,toSD1,qSD1] = buildDirectedEdgesHarmonic( ...
    mask,idxMap,D_Fe2_sd,dx,dy);
[fromSD2,toSD2,qSD2] = buildDirectedEdgesHarmonic( ...
    mask,idxMap,D_Fe3_sd,dx,dy);

D1_uni = NaN(size(D_Fe2_sd));
D2_uni = NaN(size(D_Fe3_sd));
D1_uni(mask) = mean(D1_sd_vec);
D2_uni(mask) = mean(D2_sd_vec);

[fromUNI1,toUNI1,qUNI1] = buildDirectedEdgesHarmonic( ...
    mask,idxMap,D1_uni,dx,dy);
[fromUNI2,toUNI2,qUNI2] = buildDirectedEdgesHarmonic( ...
    mask,idxMap,D2_uni,dx,dy);

c1_uni_vec = alphaMatch*mean(c1_sd_vec)*ones(K,1);
c2_uni_vec = c2_0*ones(K,1);

outSD = buildOutgoingLists(K,fromSD1,toSD1,qSD1,fromSD2,toSD2,qSD2);
outUNI = buildOutgoingLists(K,fromUNI1,toUNI1,qUNI1, ...
    fromUNI2,toUNI2,qUNI2);

usePar = configureParallel(USE_PARALLEL,AUTO_START_POOL,N_WORKERS);

%% ========================================================================
% Independent burden validation ensemble
% ========================================================================
fprintf('\nRunning held-out burden ensemble for SD-RDME...\n');
burdenSD = runHeldoutBurdenEnsembleNSM(c1_sd_vec,c2_sd_vec,outSD, ...
    x1_0,x2_0,tvec,wPred,wRest,nRepBurden,seedBase+1000,usePar);

fprintf('Running held-out burden ensemble for matched uniform RDME...\n');
burdenUNI = runHeldoutBurdenEnsembleNSM(c1_uni_vec,c2_uni_vec,outUNI, ...
    x1_0,x2_0,tvec,wPred,wRest,nRepBurden,seedBase+2000,usePar);

globalSD  = burdenSD.globalMat/K;
globalUNI = burdenUNI.globalMat/K;

windowSumSD = burdenSD.finalVoxelMat*win.W.';
windowSumUNI = burdenUNI.finalVoxelMat*win.W.';
windowDensitySD = windowSumSD./win.size.';
windowDensityUNI = windowSumUNI./win.size.';

meanWinSD = mean(windowDensitySD,1).';
meanWinUNI = mean(windowDensityUNI,1).';
deltaBurden = meanWinSD-meanWinUNI;
deltaBurdenMap = vectorToMaskMap(deltaBurden,mask);

globalMismatch = abs(mean(globalSD(:,end))-mean(globalUNI(:,end)))/ ...
    max(mean(globalSD(:,end)),eps);

%% ========================================================================
% Independent FPT validation ensemble
% ========================================================================
fprintf('\nRunning held-out all-window FPT ensemble for SD-RDME...\n');
[FPT_SD,predFPT_SD,restFPT_SD] = runHeldoutFPTEnsembleNSM( ...
    c1_sd_vec,c2_sd_vec,outSD,win,x1_0,x2_0,Mwindow,Tobs, ...
    wPred,wRest,mStar,nRepFPT,seedBase+3000,usePar);

fprintf('Running held-out all-window FPT ensemble for matched uniform RDME...\n');
[FPT_UNI,predFPT_UNI,restFPT_UNI] = runHeldoutFPTEnsembleNSM( ...
    c1_uni_vec,c2_uni_vec,outUNI,win,x1_0,x2_0,Mwindow,Tobs, ...
    wPred,wRest,mStar,nRepFPT,seedBase+4000,usePar);

FPT_SD_c = min(FPT_SD,Tobs);
FPT_UNI_c = min(FPT_UNI,Tobs);
RMST_SD = mean(FPT_SD_c,1).';
RMST_UNI = mean(FPT_UNI_c,1).';
deltaRMST = RMST_UNI-RMST_SD;
deltaRMSTMap = vectorToMaskMap(deltaRMST,mask);

%% ---------------------- Held-out spatial validation ---------------------
nTop = nnz(predMaskVec);
obsBurdenTop = topNMask(deltaBurden,nTop);
obsFPTTop = topNMask(deltaRMST,nTop);

pointBurden = computeSpatialMetrics(c1Window,deltaBurden, ...
    predMaskVec,obsBurdenTop);
pointFPT = computeSpatialMetrics(c1Window,deltaRMST, ...
    predMaskVec,obsFPTTop);

% Trajectory bootstrap retains the overlapping-window dependence structure.
fprintf('\nBootstrapping held-out burden localization metrics...\n');
bootBurdenFields = bootstrapMeanDifferenceFields( ...
    windowDensitySD,windowDensityUNI,nBoot,seedBase+5000,usePar);

fprintf('Bootstrapping held-out FPT localization metrics...\n');
bootFPTFields = bootstrapMeanDifferenceFields( ...
    FPT_UNI_c,FPT_SD_c,nBoot,seedBase+6000,usePar);

bootMetricsBurden = bootstrapSpatialMetrics( ...
    c1Window,predMaskVec,bootBurdenFields);
bootMetricsFPT = bootstrapSpatialMetrics( ...
    c1Window,predMaskVec,bootFPTFields);

ciBurden = metricCI(bootMetricsBurden,alphaCI);
ciFPT = metricCI(bootMetricsFPT,alphaCI);

%% ---------------------- Predicted-zone validation -----------------------
predEndSD = burdenSD.predMat(:,end);
predEndUNI = burdenUNI.predMat(:,end);
restEndSD = burdenSD.restMat(:,end);
restEndUNI = burdenUNI.restMat(:,end);

[deltaPredBurden,ciPredBurden] = bootstrapDifferenceIndependent( ...
    predEndSD,predEndUNI,nBoot,alphaCI,seedBase+7000);
[deltaRestBurden,ciRestBurden] = bootstrapDifferenceIndependent( ...
    restEndSD,restEndUNI,nBoot,alphaCI,seedBase+7100);

predRatioBurden = mean(predEndSD)/max(mean(predEndUNI),eps);
restRatioBurden = mean(restEndSD)/max(mean(restEndUNI),eps);

predRMST_SD = min(predFPT_SD,Tobs);
predRMST_UNI = min(predFPT_UNI,Tobs);
restRMST_SD = min(restFPT_SD,Tobs);
restRMST_UNI = min(restFPT_UNI,Tobs);

[deltaPredRMST,ciPredRMST] = bootstrapDifferenceIndependent( ...
    predRMST_UNI,predRMST_SD,nBoot,alphaCI,seedBase+7200);
[deltaRestRMST,ciRestRMST] = bootstrapDifferenceIndependent( ...
    restRMST_UNI,restRMST_SD,nBoot,alphaCI,seedBase+7300);

PcrossPredSD = mean(predFPT_SD<=Tobs);
PcrossPredUNI = mean(predFPT_UNI<=Tobs);
PcrossRestSD = mean(restFPT_SD<=Tobs);
PcrossRestUNI = mean(restFPT_UNI<=Tobs);

%% ---------------------- Quantile robustness -----------------------------
Jq = zeros(numel(qList),2);
Eq = zeros(numel(qList),2);
for iq = 1:numel(qList)
    nQ = max(1,ceil((1-qList(iq))*K));
    predQ = topNMask(c1Window,nQ);
    obsBQ = topNMask(deltaBurden,nQ);
    obsFQ = topNMask(deltaRMST,nQ);
    [Jq(iq,1),Eq(iq,1)] = overlapMetricsEqualSize(predQ,obsBQ);
    [Jq(iq,2),Eq(iq,2)] = overlapMetricsEqualSize(predQ,obsFQ);
end

%% ---------------------- Audits and console output -----------------------
sizeRhoBurden = spearmanLocal(win.size,deltaBurden);
sizeRhoFPT = spearmanLocal(win.size,deltaRMST);

fprintf('\n================ Held-out validation summary ================\n');
fprintf('Independent SSA global mismatch          = %.3f %%\n', ...
    100*globalMismatch);
fprintf('\nSpatial localization of held-out burden:\n');
printMetricLine(pointBurden,ciBurden);
fprintf('\nSpatial localization of held-out FPT:\n');
printMetricLine(pointFPT,ciFPT);

fprintf('\nPredicted top-15%% zone:\n');
fprintf('Final burden SD / U                      = %.4f / %.4f\n', ...
    mean(predEndSD),mean(predEndUNI));
fprintf('Delta burden [95%% CI]                    = %.4f [%.4f, %.4f]\n', ...
    deltaPredBurden,ciPredBurden(1),ciPredBurden(2));
fprintf('Burden ratio SD/U                        = %.4f\n',predRatioBurden);
fprintf('Delta RMST60 [95%% CI]                   = %.4f [%.4f, %.4f] min\n', ...
    deltaPredRMST,ciPredRMST(1),ciPredRMST(2));
fprintf('Pcross60 SD / U                          = %.3f / %.3f\n', ...
    PcrossPredSD,PcrossPredUNI);

fprintf('\nRemaining windows:\n');
fprintf('Delta burden [95%% CI]                    = %.4f [%.4f, %.4f]\n', ...
    deltaRestBurden,ciRestBurden(1),ciRestBurden(2));
fprintf('Burden ratio SD/U                        = %.4f\n',restRatioBurden);
fprintf('Delta RMST60 [95%% CI]                   = %.4f [%.4f, %.4f] min\n', ...
    deltaRestRMST,ciRestRMST(1),ciRestRMST(2));
fprintf('Pcross60 SD / U                          = %.3f / %.3f\n', ...
    PcrossRestSD,PcrossRestUNI);

fprintf('\nWindow-size audit:\n');
fprintf('rho(size, Delta burden)                  = %.4f\n',sizeRhoBurden);
fprintf('rho(size, Delta RMST60)                  = %.4f\n',sizeRhoFPT);
fprintf('=============================================================\n');

%% ========================================================================
% Main figure: held-out spatial and regional validation
% ========================================================================
figMain = figure('Units','centimeters','Position',[1.0 1.0 18.8 10.2], ...
    'Color','w');
tl = tiledlayout(figMain,2,3,'TileSpacing','compact','Padding','compact');

% A. Frozen c1-hat prediction.
ax = nexttile(tl,1); hold(ax,'on');
plotMaskedMap(ax,x,y,c1WindowMap,mask,natureBlueYellow(256),[0 1]);
contour(ax,Xg,Yg,double(predMaskMap),[0.5 0.5], ...
    'Color',[1 1 1],'LineStyle','--','LineWidth',1.25);
contour(ax,Xg,Yg,double(mask),[0.5 0.5], ...
    'Color',[0.18 0.18 0.18],'LineWidth',0.55);
cb = colorbar(ax); styleColorbar(cb,'Window-mean $\widehat c_1$','latex');
title(ax,'Pre-specified reaction-field prediction', ...
    'FontSize',8.2,'FontWeight','normal');
text(ax,0.02,0.98,'A','Units','normalized','FontSize',12, ...
    'FontWeight','bold','VerticalAlignment','top','Color',colText);
text(ax,0.04,0.08,'white dashed: predicted top 15%', ...
    'Units','normalized','FontSize',5.7,'BackgroundColor','w', ...
    'Margin',1,'Color',colText);
formatMapAxes(ax);

% B. Held-out burden redistribution.
ax = nexttile(tl,2); hold(ax,'on');
bMax = max(quantileLocal(abs(deltaBurden),0.98),eps);
plotMaskedMap(ax,x,y,deltaBurdenMap,mask,cmapBlueWhiteRed(256),[-bMax bMax]);
contour(ax,Xg,Yg,double(predMaskMap),[0.5 0.5], ...
    'Color',[0.08 0.08 0.08],'LineStyle','--','LineWidth',1.15);
contour(ax,Xg,Yg,double(mask),[0.5 0.5], ...
    'Color',[0.18 0.18 0.18],'LineWidth',0.55);
cb = colorbar(ax); styleColorbar(cb,'$\Delta$ local burden','latex');
title(ax,'Held-out local-burden redistribution', ...
    'FontSize',8.2,'FontWeight','normal');
text(ax,0.02,0.98,'B','Units','normalized','FontSize',12, ...
    'FontWeight','bold','VerticalAlignment','top','Color',colText);
text(ax,0.04,0.08,sprintf('$J=%.2f,\\ \\rho_s=%.2f$', ...
    pointBurden.J,pointBurden.rho),'Units','normalized', ...
    'Interpreter','latex','FontSize',6.1,'BackgroundColor','w', ...
    'Margin',1,'Color',colText);
formatMapAxes(ax);

% C. Held-out FPT redistribution.
ax = nexttile(tl,3); hold(ax,'on');
fMax = max(quantileLocal(abs(deltaRMST),0.98),eps);
plotMaskedMap(ax,x,y,deltaRMSTMap,mask,cmapBlueWhiteRed(256),[-fMax fMax]);
contour(ax,Xg,Yg,double(predMaskMap),[0.5 0.5], ...
    'Color',[0.08 0.08 0.08],'LineStyle','--','LineWidth',1.15);
contour(ax,Xg,Yg,double(mask),[0.5 0.5], ...
    'Color',[0.18 0.18 0.18],'LineWidth',0.55);
cb = colorbar(ax); styleColorbar(cb,'$\Delta\mathrm{RMST}_{60}$ (min)','latex');
title(ax,'Held-out local-FPT redistribution', ...
    'FontSize',8.2,'FontWeight','normal');
text(ax,0.02,0.98,'C','Units','normalized','FontSize',12, ...
    'FontWeight','bold','VerticalAlignment','top','Color',colText);
text(ax,0.04,0.08,sprintf('$J=%.2f,\\ \\rho_s=%.2f$', ...
    pointFPT.J,pointFPT.rho),'Units','normalized', ...
    'Interpreter','latex','FontSize',6.1,'BackgroundColor','w', ...
    'Margin',1,'Color',colText);
formatMapAxes(ax);

% D. Predicted-zone burden trajectory.
ax = nexttile(tl,4); hold(ax,'on');
plotMeanSEM(ax,tvec,burdenUNI.predMat,colUNI,0.13);
plotMeanSEM(ax,tvec,burdenSD.predMat,colSD,0.13);
title(ax,'Predicted-zone ROS burden','FontSize',8.2,'FontWeight','normal');
text(ax,0.02,0.98,'D','Units','normalized','FontSize',12, ...
    'FontWeight','bold','VerticalAlignment','top','Color',colText);
xlabel(ax,'Time (min)');
ylabel(ax,'Mean local burden density');
legend(ax,{'Uniform RDME','SD-RDME'},'Box','off', ...
    'Location','northwest','FontSize',5.8,'Interpreter','tex');
text(ax,0.97,0.08,sprintf('$\\Delta B_{60}=%.2f$\\newline$SD/U=%.2f$', ...
    deltaPredBurden,predRatioBurden),'Units','normalized', ...
    'HorizontalAlignment','right','Interpreter','latex', ...
    'FontSize',6.0,'Color',colText);
xlim(ax,[0 Tobs]);
formatAxes(ax);

% E. Predicted-zone FPT CDF.
ax = nexttile(tl,5); hold(ax,'on');
[tU,FU] = finiteECDF(predFPT_UNI,Tobs);
[tS,FS] = finiteECDF(predFPT_SD,Tobs);
stairs(ax,tU,FU,'-','Color',colUNI,'LineWidth',1.25);
stairs(ax,tS,FS,'-','Color',colSD,'LineWidth',1.25);
title(ax,'Predicted-zone threshold crossing', ...
    'FontSize',8.2,'FontWeight','normal');
text(ax,0.02,0.98,'E','Units','normalized','FontSize',12, ...
    'FontWeight','bold','VerticalAlignment','top','Color',colText);
xlabel(ax,'FPT (min)');
ylabel(ax,'Cumulative crossing probability');
legend(ax,{'Uniform RDME','SD-RDME'},'Box','off', ...
    'Location','northwest','FontSize',5.8,'Interpreter','tex');
text(ax,0.97,0.09,sprintf('$\\Delta\\mathrm{RMST}_{60}=%.2f$ min\\newline$P_{60}: %.2f\\rightarrow%.2f$', ...
    deltaPredRMST,PcrossPredUNI,PcrossPredSD), ...
    'Units','normalized','HorizontalAlignment','right', ...
    'Interpreter','latex','FontSize',5.9,'Color',colText);
xlim(ax,[0 Tobs]); ylim(ax,[0 1]);
formatAxes(ax);

% F. Held-out validation metric summary.
ax = nexttile(tl,6); hold(ax,'on');
metricNames = {'$\rho_s$','Jaccard','AUROC'};
valB = [pointBurden.rho,pointBurden.J,pointBurden.AUC];
loB  = [ciBurden.rho(1),ciBurden.J(1),ciBurden.AUC(1)];
hiB  = [ciBurden.rho(2),ciBurden.J(2),ciBurden.AUC(2)];
valF = [pointFPT.rho,pointFPT.J,pointFPT.AUC];
loF  = [ciFPT.rho(1),ciFPT.J(1),ciFPT.AUC(1)];
hiF  = [ciFPT.rho(2),ciFPT.J(2),ciFPT.AUC(2)];

xx = 1:3;
errorbar(ax,xx-0.08,valB,valB-loB,hiB-valB,'o', ...
    'Color',colBurden,'MarkerFaceColor',colBurden, ...
    'MarkerSize',4.2,'LineWidth',0.9,'CapSize',4);
errorbar(ax,xx+0.08,valF,valF-loF,hiF-valF,'s', ...
    'Color',colFPT,'MarkerFaceColor','w', ...
    'MarkerSize',4.2,'LineWidth',0.9,'CapSize',4);
chanceJ = (1-qMain)/(2-(1-qMain));
plot(ax,[0.5 3.5],[chanceJ chanceJ],':', ...
    'Color',[0.70 0.70 0.70],'LineWidth',0.7, ...
    'HandleVisibility','off');
set(ax,'XTick',xx,'XTickLabel',metricNames, ...
    'TickLabelInterpreter','latex');
xlim(ax,[0.5 3.5]); ylim(ax,[0 1.02]);
ylabel(ax,'Held-out performance');
title(ax,'Independent localization performance', ...
    'FontSize',8.2,'FontWeight','normal');
text(ax,0.02,0.98,'F','Units','normalized','FontSize',12, ...
    'FontWeight','bold','VerticalAlignment','top','Color',colText);
legend(ax,{'local burden','local FPT'},'Box','off', ...
    'Location','southeast','FontSize',5.7);
formatAxes(ax);

if SAVE_FIG
    exportgraphics(figMain,[outMain '.pdf'],'ContentType','vector');
    exportgraphics(figMain,[outMain '.png'],'Resolution',600);
    exportgraphics(figMain,[outMain '.tif'],'Resolution',600);
    savefig(figMain,[outMain '.fig']);
end

%% ========================================================================
% Supplementary validation audits
% ========================================================================
figSI = figure('Units','centimeters','Position',[1.0 2.0 18.5 5.4], ...
    'Color','w');
tlSI = tiledlayout(figSI,1,4,'TileSpacing','compact','Padding','compact');

% S1. Independent global matching audit.
ax = nexttile(tlSI,1); hold(ax,'on');
plotMeanSEM(ax,tvec,globalUNI,colUNI,0.12);
plotMeanSEM(ax,tvec,globalSD,colSD,0.12);
title(ax,'Held-out global matching','FontSize',8.1,'FontWeight','normal');
xlabel(ax,'Time (min)'); ylabel(ax,'Global burden per voxel');
text(ax,0.02,0.98,'A','Units','normalized','FontSize',11.5, ...
    'FontWeight','bold','VerticalAlignment','top');
text(ax,0.97,0.08,sprintf('mismatch = %.2f%%',100*globalMismatch), ...
    'Units','normalized','HorizontalAlignment','right','FontSize',5.8);
formatAxes(ax);

% S2. Predicted versus remaining burden effects.
ax = nexttile(tlSI,2); hold(ax,'on');
vals = [deltaPredBurden,deltaRestBurden];
lo = [ciPredBurden(1),ciRestBurden(1)];
hi = [ciPredBurden(2),ciRestBurden(2)];
errorbar(ax,1:2,vals,vals-lo,hi-vals,'o', ...
    'Color',colSD,'MarkerFaceColor','w','MarkerSize',4.5, ...
    'LineWidth',1.0,'CapSize',5);
plot(ax,[0.5 2.5],[0 0],':','Color',[0.65 0.65 0.65], ...
    'LineWidth',0.7);
set(ax,'XTick',1:2,'XTickLabel',{'predicted top 15%','remaining'});
ylabel(ax,'$\Delta$ burden at 60 min','Interpreter','latex');
title(ax,'Spatial selectivity of burden effect', ...
    'FontSize',8.1,'FontWeight','normal');
text(ax,0.02,0.98,'B','Units','normalized','FontSize',11.5, ...
    'FontWeight','bold','VerticalAlignment','top');
formatAxes(ax);

% S3. Predicted versus remaining FPT effects.
ax = nexttile(tlSI,3); hold(ax,'on');
vals = [deltaPredRMST,deltaRestRMST];
lo = [ciPredRMST(1),ciRestRMST(1)];
hi = [ciPredRMST(2),ciRestRMST(2)];
errorbar(ax,1:2,vals,vals-lo,hi-vals,'o', ...
    'Color',colFPT,'MarkerFaceColor','w','MarkerSize',4.5, ...
    'LineWidth',1.0,'CapSize',5);
plot(ax,[0.5 2.5],[0 0],':','Color',[0.65 0.65 0.65], ...
    'LineWidth',0.7);
set(ax,'XTick',1:2,'XTickLabel',{'predicted top 15%','remaining'});
ylabel(ax,'$\Delta\mathrm{RMST}_{60}$ (min)','Interpreter','latex');
title(ax,'Spatial selectivity of FPT effect', ...
    'FontSize',8.1,'FontWeight','normal');
text(ax,0.02,0.98,'C','Units','normalized','FontSize',11.5, ...
    'FontWeight','bold','VerticalAlignment','top');
formatAxes(ax);

% S4. q robustness on held-out trajectories.
ax = nexttile(tlSI,4); hold(ax,'on');
plot(ax,qList,Jq(:,1),'-o','Color',colBurden, ...
    'MarkerFaceColor','w','MarkerSize',3.3,'LineWidth',1.0);
plot(ax,qList,Jq(:,2),'--s','Color',colFPT, ...
    'MarkerFaceColor','w','MarkerSize',3.2,'LineWidth',1.0);
plot(ax,[qMain qMain],[0 1],':','Color',[0.60 0.60 0.60], ...
    'LineWidth',0.7,'HandleVisibility','off');
xlabel(ax,'Quantile $q$','Interpreter','latex');
ylabel(ax,'Jaccard overlap');
title(ax,'Held-out quantile robustness', ...
    'FontSize',8.1,'FontWeight','normal');
text(ax,0.02,0.98,'D','Units','normalized','FontSize',11.5, ...
    'FontWeight','bold','VerticalAlignment','top');
legend(ax,{'local burden','local FPT'},'Box','off', ...
    'Location','southwest','FontSize',5.5);
xlim(ax,[min(qList) max(qList)]); ylim(ax,[0 1]);
formatAxes(ax);

if SAVE_FIG
    exportgraphics(figSI,[outSI '.pdf'],'ContentType','vector');
    exportgraphics(figSI,[outSI '.png'],'Resolution',600);
    exportgraphics(figSI,[outSI '.tif'],'Resolution',600);
    savefig(figSI,[outSI '.fig']);
end

%% ---------------------- Save data and summary ---------------------------
summaryTable = table( ...
    pointBurden.rho,ciBurden.rho(1),ciBurden.rho(2), ...
    pointBurden.J,ciBurden.J(1),ciBurden.J(2), ...
    pointBurden.E,pointBurden.AUC,ciBurden.AUC(1),ciBurden.AUC(2), ...
    pointBurden.AP, ...
    pointFPT.rho,ciFPT.rho(1),ciFPT.rho(2), ...
    pointFPT.J,ciFPT.J(1),ciFPT.J(2), ...
    pointFPT.E,pointFPT.AUC,ciFPT.AUC(1),ciFPT.AUC(2), ...
    pointFPT.AP, ...
    deltaPredBurden,ciPredBurden(1),ciPredBurden(2),predRatioBurden, ...
    deltaPredRMST,ciPredRMST(1),ciPredRMST(2), ...
    PcrossPredUNI,PcrossPredSD,globalMismatch, ...
    'VariableNames',{ ...
    'Burden_Spearman','Burden_Spearman_CI_L','Burden_Spearman_CI_U', ...
    'Burden_Jaccard','Burden_Jaccard_CI_L','Burden_Jaccard_CI_U', ...
    'Burden_Enrichment','Burden_AUROC','Burden_AUROC_CI_L', ...
    'Burden_AUROC_CI_U','Burden_AP', ...
    'FPT_Spearman','FPT_Spearman_CI_L','FPT_Spearman_CI_U', ...
    'FPT_Jaccard','FPT_Jaccard_CI_L','FPT_Jaccard_CI_U', ...
    'FPT_Enrichment','FPT_AUROC','FPT_AUROC_CI_L', ...
    'FPT_AUROC_CI_U','FPT_AP', ...
    'Predicted_DeltaBurden','Predicted_DeltaBurden_CI_L', ...
    'Predicted_DeltaBurden_CI_U','Predicted_BurdenRatio', ...
    'Predicted_DeltaRMST60','Predicted_DeltaRMST60_CI_L', ...
    'Predicted_DeltaRMST60_CI_U','Predicted_Pcross60_UNI', ...
    'Predicted_Pcross60_SD','GlobalMismatch'});

if SAVE_DATA
    writetable(summaryTable,outCSV);
    validationMeta = struct();
    validationMeta.RUN_MODE = RUN_MODE;
    validationMeta.lockedFile = LOCKED_FILE;
    validationMeta.seedBase = seedBase;
    validationMeta.qMain = qMain;
    validationMeta.qList = qList;
    validationMeta.nRepBurden = nRepBurden;
    validationMeta.nRepFPT = nRepFPT;
    validationMeta.nBoot = nBoot;
    validationMeta.windowRadius = windowRadius;
    validationMeta.mStar = mStar;
    validationMeta.alphaMatch = alphaMatch;

    save(outMAT,'validationMeta','x','y','mask','c1Window', ...
        'c1WindowMap','predMaskVec','predMaskMap','wPred','wRest', ...
        'globalSD','globalUNI','windowDensitySD','windowDensityUNI', ...
        'deltaBurden','deltaBurdenMap','FPT_SD','FPT_UNI', ...
        'deltaRMST','deltaRMSTMap','predFPT_SD','predFPT_UNI', ...
        'restFPT_SD','restFPT_UNI','pointBurden','pointFPT', ...
        'ciBurden','ciFPT','Jq','Eq','qList','summaryTable','-v7.3');
end

fprintf('\nCompleted. Outputs:\n');
fprintf('  %s.{pdf,png,tif,fig}\n',outMain);
fprintf('  %s.{pdf,png,tif,fig}\n',outSI);
fprintf('  %s\n',outCSV);
fprintf('  %s\n',outMAT);

%% ==================== Validation-specific functions =====================

function w = buildAverageWindowWeight(win,centerMask)
    centerMask = logical(centerMask(:));
    n = nnz(centerMask);
    if n<1
        error('The requested center set is empty.');
    end
    rowWeight = zeros(numel(centerMask),1);
    rowWeight(centerMask) = 1./(n*win.size(centerMask));
    w = full(win.W.'*rowWeight);
    w = w(:);
    if abs(sum(w)-1)>1e-10
        error('Average-window weight does not sum to one.');
    end
end

function out = runHeldoutBurdenEnsembleNSM(c1,c2,outgoing,x1_0,x2_0, ...
    tvec,wPred,wRest,nRep,seedBase,usePar)

    K = numel(c1);
    nT = numel(tvec);
    globalMat = zeros(nRep,nT);
    predMat = zeros(nRep,nT);
    restMat = zeros(nRep,nT);
    finalVoxelMat = zeros(nRep,K);

    if usePar
        parfor r = 1:nRep
            rng(seedBase+r,'twister');
            [g,pd,rd,v] = simulateHeldoutBurdenOneNSM( ...
                c1,c2,outgoing,x1_0,x2_0,tvec,wPred,wRest);
            globalMat(r,:) = g;
            predMat(r,:) = pd;
            restMat(r,:) = rd;
            finalVoxelMat(r,:) = v.';
        end
    else
        for r = 1:nRep
            rng(seedBase+r,'twister');
            [g,pd,rd,v] = simulateHeldoutBurdenOneNSM( ...
                c1,c2,outgoing,x1_0,x2_0,tvec,wPred,wRest);
            globalMat(r,:) = g;
            predMat(r,:) = pd;
            restMat(r,:) = rd;
            finalVoxelMat(r,:) = v.';
            if mod(r,max(1,round(nRep/10)))==0
                fprintf('  burden trajectories: %d / %d\n',r,nRep);
            end
        end
    end

    out.globalMat = globalMat;
    out.predMat = predMat;
    out.restMat = restMat;
    out.finalVoxelMat = finalVoxelMat;
end

function [globalRec,predRec,restRec,x3] = ...
    simulateHeldoutBurdenOneNSM(c1,c2,out,x1_0,x2_0,tvec,wPred,wRest)

    K = numel(c1);
    nT = numel(tvec);
    x1 = x1_0(:);
    x2 = x2_0(:);
    x3 = zeros(K,1);

    voxelRate = x1.*(c1+out.q1sum)+x2.*(c2+out.q2sum);
    tree = fenwickBuild(voxelRate);

    globalRec = zeros(1,nT);
    predRec = zeros(1,nT);
    restRec = zeros(1,nT);

    globalROS = 0;
    predBurden = 0;
    restBurden = 0;

    t = tvec(1);
    outIdx = 1;
    Tend = tvec(end);

    while t<Tend
        a0 = fenwickTotal(tree);
        if ~(a0>0) || ~isfinite(a0), break; end
        tau = -log(rand)/a0;
        tNext = t+tau;

        while outIdx<nT && tvec(outIdx+1)<=min(tNext,Tend)
            outIdx = outIdx+1;
            globalRec(outIdx) = globalROS;
            predRec(outIdx) = predBurden;
            restRec(outIdx) = restBurden;
        end
        if tNext>Tend, break; end

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
            x3(k) = x3(k)+1;
            globalROS = globalROS+1;
            predBurden = predBurden+wPred(k);
            restBurden = restBurden+wRest(k);
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

        [voxelRate,tree] = refreshVoxelRates(unique(touched), ...
            voxelRate,tree,x1,x2,c1,c2,out.q1sum,out.q2sum);
        t = tNext;
    end

    while outIdx<nT
        outIdx = outIdx+1;
        globalRec(outIdx) = globalROS;
        predRec(outIdx) = predBurden;
        restRec(outIdx) = restBurden;
    end
end

function [FPT,predFPT,restFPT] = runHeldoutFPTEnsembleNSM( ...
    c1,c2,outgoing,win,x1_0,x2_0,Mwindow,Tobs,wPred,wRest,mStar, ...
    nRep,seedBase,usePar)

    K = numel(c1);
    FPT = inf(nRep,K);
    predFPT = inf(nRep,1);
    restFPT = inf(nRep,1);

    if usePar
        parfor r = 1:nRep
            rng(seedBase+r,'twister');
            [f,pf,rf] = simulateHeldoutFPTOneNSM( ...
                c1,c2,outgoing,win,x1_0,x2_0,Mwindow,Tobs, ...
                wPred,wRest,mStar);
            FPT(r,:) = f.';
            predFPT(r) = pf;
            restFPT(r) = rf;
        end
    else
        for r = 1:nRep
            rng(seedBase+r,'twister');
            [f,pf,rf] = simulateHeldoutFPTOneNSM( ...
                c1,c2,outgoing,win,x1_0,x2_0,Mwindow,Tobs, ...
                wPred,wRest,mStar);
            FPT(r,:) = f.';
            predFPT(r) = pf;
            restFPT(r) = rf;
            if mod(r,max(1,round(nRep/10)))==0
                fprintf('  FPT trajectories: %d / %d\n',r,nRep);
            end
        end
    end
end

function [fpt,predFPT,restFPT] = simulateHeldoutFPTOneNSM( ...
    c1,c2,out,win,x1_0,x2_0,Mwindow,Tobs,wPred,wRest,mStar)

    K = numel(c1);
    x1 = x1_0(:);
    x2 = x2_0(:);

    localROS = zeros(K,1);
    fpt = inf(K,1);
    active = true(K,1);

    predBurden = 0;
    restBurden = 0;
    predFPT = inf;
    restFPT = inf;

    voxelRate = x1.*(c1+out.q1sum)+x2.*(c2+out.q2sum);
    tree = fenwickBuild(voxelRate);
    t = 0;

    while t<Tobs && (any(active) || isinf(predFPT) || isinf(restFPT))
        a0 = fenwickTotal(tree);
        if ~(a0>0) || ~isfinite(a0), break; end
        tau = -log(rand)/a0;
        tNext = t+tau;
        if tNext>Tobs, break; end

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

            centers = win.centersByVoxel{k};
            if ~isempty(centers)
                centers = centers(active(centers));
                localROS(centers) = localROS(centers)+1;
                crossed = centers(localROS(centers)>=Mwindow(centers));
                if ~isempty(crossed)
                    fpt(crossed) = tNext;
                    active(crossed) = false;
                end
            end

            predBurden = predBurden+wPred(k);
            restBurden = restBurden+wRest(k);

            if isinf(predFPT) && predBurden>=mStar
                predFPT = tNext;
            end
            if isinf(restFPT) && restBurden>=mStar
                restFPT = tNext;
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

        [voxelRate,tree] = refreshVoxelRates(unique(touched), ...
            voxelRate,tree,x1,x2,c1,c2,out.q1sum,out.q2sum);
        t = tNext;
    end
end

function M = computeSpatialMetrics(score,target,predTop,obsTop)
    M = struct();
    M.rho = spearmanLocal(score,target);
    [M.J,M.E] = overlapMetricsEqualSize(predTop,obsTop);
    [M.AUC,M.AP] = binaryDiscrimination(score,obsTop);
end

function B = bootstrapSpatialMetrics(score,predTop,bootFields)
    nBoot = size(bootFields,1);
    nTop = nnz(predTop);
    B.rho = zeros(nBoot,1);
    B.J = zeros(nBoot,1);
    B.E = zeros(nBoot,1);
    B.AUC = zeros(nBoot,1);
    B.AP = zeros(nBoot,1);

    for b = 1:nBoot
        target = bootFields(b,:).';
        obsTop = topNMask(target,nTop);
        B.rho(b) = spearmanLocal(score,target);
        [B.J(b),B.E(b)] = overlapMetricsEqualSize(predTop,obsTop);
        [B.AUC(b),B.AP(b)] = binaryDiscrimination(score,obsTop);
    end
end

function C = metricCI(B,alpha)
    C.rho = percentileCI(B.rho,alpha);
    C.J = percentileCI(B.J,alpha);
    C.E = percentileCI(B.E,alpha);
    C.AUC = percentileCI(B.AUC,alpha);
    C.AP = percentileCI(B.AP,alpha);
end

function printMetricLine(M,C)
    fprintf('rho_s [95%% CI]                         = %.4f [%.4f, %.4f]\n', ...
        M.rho,C.rho(1),C.rho(2));
    fprintf('Jaccard [95%% CI]                      = %.4f [%.4f, %.4f]\n', ...
        M.J,C.J(1),C.J(2));
    fprintf('Enrichment                             = %.4f\n',M.E);
    fprintf('AUROC [95%% CI]                        = %.4f [%.4f, %.4f]\n', ...
        M.AUC,C.AUC(1),C.AUC(2));
    fprintf('Average precision                      = %.4f\n',M.AP);
end

function bootFields = bootstrapMeanDifferenceFields(A,B,nBoot,seed,usePar)
    % Returns mean(A*) - mean(B*) for each complete-trajectory bootstrap.
    nA = size(A,1);
    nB = size(B,1);
    K = size(A,2);
    bootFields = zeros(nBoot,K);

    if usePar
        parfor b = 1:nBoot
            rng(seed+b,'twister');
            iA = randi(nA,nA,1);
            iB = randi(nB,nB,1);
            bootFields(b,:) = mean(A(iA,:),1)-mean(B(iB,:),1);
        end
    else
        for b = 1:nBoot
            rng(seed+b,'twister');
            iA = randi(nA,nA,1);
            iB = randi(nB,nB,1);
            bootFields(b,:) = mean(A(iA,:),1)-mean(B(iB,:),1);
            if mod(b,max(1,round(nBoot/10)))==0
                fprintf('  bootstrap fields: %d / %d\n',b,nBoot);
            end
        end
    end
end

function m = topNMask(v,nTop)
    v = v(:);
    K = numel(v);
    nTop = min(max(round(nTop),1),K);
    [~,ord] = sortrows([-v,(1:K).'],[1 2]);
    m = false(K,1);
    m(ord(1:nTop)) = true;
end

function [J,E] = overlapMetricsEqualSize(a,b)
    a = logical(a(:));
    b = logical(b(:));
    inter = nnz(a&b);
    uni = nnz(a|b);
    J = inter/max(uni,1);
    expected = nnz(a)*nnz(b)/numel(a);
    E = inter/max(expected,eps);
end

function [AUC,AP] = binaryDiscrimination(score,label)
    score = score(:);
    label = logical(label(:));
    nPos = nnz(label);
    nNeg = numel(label)-nPos;
    if nPos==0 || nNeg==0
        AUC = NaN;
        AP = NaN;
        return;
    end

    r = tiedRankLocal(score);
    AUC = (sum(r(label))-nPos*(nPos+1)/2)/(nPos*nNeg);

    [~,ord] = sortrows([-score,(1:numel(score)).'],[1 2]);
    yy = label(ord);
    tp = cumsum(yy);
    precision = tp./(1:numel(yy)).';
    AP = sum(precision(yy))/nPos;
end

function ci = percentileCI(v,alpha)
    ci = [quantileLocal(v,alpha/2),quantileLocal(v,1-alpha/2)];
end

function v = normalizeVector(v)
    v = v(:);
    lo = min(v);
    hi = max(v);
    v = (v-lo)/(hi-lo+eps);
end

function [xF,F] = finiteECDF(t,Tobs)

    t = t(:);
    nTotal = numel(t);

    % Only actual threshold-crossing events contribute CDF jumps.
    events = sort(t(isfinite(t) & t <= Tobs));

    if isempty(events)
        xF = [0; Tobs];
        F  = [0; 0];
    else
        xF = [0; events; Tobs];
        F  = [0; (1:numel(events))'/nTotal; numel(events)/nTotal];
    end
end

function cmap = natureBlueYellow(n)
    anchors = [0.06 0.13 0.34
               0.10 0.36 0.57
               0.23 0.61 0.61
               0.65 0.74 0.38
               0.96 0.82 0.28];
    cmap = interpCmap(anchors,n);
end


%% ================= Reused validated Fig. 3B-C functions ================


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
