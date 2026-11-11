%% SDRDME_CoreRobustness.m
% -------------------------------------------------------------------------
% FINAL targeted robustness analysis for the mitochondrial SD-RDME.
%
% PURPOSE
%   Test whether the two load-bearing main-text phenotypes—local ROS
%   burden redistribution and Full-vs-Uniform local FPT advancement—
%   depend on four analyst/numerical choices:
%       A. high-risk region quantile q;
%       B. fixed-radius local window r_W;
%       C. regional FPT threshold multiplier;
%       D. spatial-grid resolution.
%
% THE ANALYSIS DELIBERATELY DOES NOT SCAN ALL BIOCHEMICAL PARAMETERS.
% It preserves the finalized Fig. 3 geometry, fields, reaction/diffusion
% parameters, initial iron state, exact first-moment global matching,
% harmonic interface rates, and exact next-subvolume SSA.
%
% LOCKED MODEL
%   U: uniform reaction + uniform diffusion
%   H: heterogeneous reaction c1(r) + uniform diffusion
%   B: uniform reaction + heterogeneous diffusion D(r)
%   F: heterogeneous reaction c1(r) + heterogeneous diffusion D(r)
%
% PRIMARY DISPLAYED ENDPOINTS
%   Region quantile:       local burden ratio B_A^F/B_A^U and
%                          Delta_UF = RMST_U - RMST_F.
%   Local-window radius:   Delta H_B = H_B^F-H_B^U and Delta_UF.
%   FPT threshold:         Delta_UF and P(T_A<=60 min) for U and F.
%   Grid resolution:       Delta H_B and Delta_UF.
%
% U/H/B/F mechanism effects (E_H, E_B and E_B|H) are still computed for
% every condition and exported as an audit table and a compact secondary
% figure, but they no longer dominate all four primary robustness panels.
%
% FOUR ROBUSTNESS AXES
%   qList                 = [0.80 0.85 0.90]
%   windowRadiusList      = [0.08 0.10 0.12]
%   thresholdMultiplier   = [0.80 1.00 1.20 1.40]
%   gridList              = [48x29, 60x36, 72x43]
%
% GRID-RESOLUTION RULE
%   Physical coordinates and all physical length scales are fixed.
%   Diffusion rates are recomputed as D_kl / h^2.
%   The 60x36 reference retains the exact 4/2 Fe2+/Fe3+ state per voxel.
%   Coarse/fine grids preserve the corresponding initial iron density by
%   scaling total molecule counts with voxel area and distributing the
%   resulting integer counts uniformly and deterministically.
%
% COMPUTATIONAL EFFICIENCY
%   All q-, radius-, and threshold-robustness endpoints on the reference
%   grid are recorded in one multi-region exact SSA ensemble per model.
%
% OUTPUTS
%   SDRDME_CoreRobustness.{pdf,png,tif,fig}
%   FigS_SDRDME_fourRobustness_mechanismAudit_v2.{pdf,png,tif,fig}
%   FigS_SDRDME_fourRobustness_audit_v2.csv
%   FigS_SDRDME_fourRobustness_gridAudit_v2.csv
%   FigS_SDRDME_fourRobustness_data_v2.mat
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
        nRepFPT = 20;
        nBoot = 200;
    case 'pilot'
        nRepFPT = 100;
        nBoot = 500;
    case 'final'
        nRepFPT = 300;
        nBoot = 1000;
    otherwise
        error('Unknown RUN_MODE: %s',RUN_MODE);
end

seedBase = 20260806;
alphaCI = 0.05;

outFigure = 'FigS_SDRDME_fourRobustness_FINAL_v2';
outMechanismFigure = 'FigS_SDRDME_fourRobustness_mechanismAudit_v2';
outAuditCSV = 'FigS_SDRDME_fourRobustness_audit_v2.csv';
outGridCSV = 'FigS_SDRDME_fourRobustness_gridAudit_v2.csv';
outMAT = 'FigS_SDRDME_fourRobustness_data_v2.mat';

%% ----------------------------- Figure style -----------------------------
set(groot,'defaultFigureColor','w');
set(groot,'defaultAxesFontName','Arial');
set(groot,'defaultTextFontName','Arial');
set(groot,'defaultAxesFontSize',7);
set(groot,'defaultTextFontSize',7);
set(groot,'defaultAxesLineWidth',0.72);
set(groot,'defaultLineLineWidth',1.05);
set(groot,'defaultAxesTickDir','out');
set(groot,'defaultAxesBox','on');
set(groot,'defaultLegendInterpreter','tex');
set(groot,'defaultAxesTickLabelInterpreter','tex');

% Main-figure endpoint colors.
colBurden = [0.08 0.45 0.62];
colFPT = [0.86 0.34 0.10];
colP60U = [0.36 0.38 0.42];
colP60F = [0.49 0.23 0.67];
colReference = [0.48 0.48 0.48];

% Secondary U/H/B/F mechanism-audit colors.
colDelta = [0.20 0.21 0.23];
colH = [0.85 0.33 0.10];
colB = [0.15 0.43 0.76];
colBH = [0.49 0.23 0.67];
effectColors = [colDelta;colH;colB;colBH];
effectLabels = {'Full vs Uniform','hotspot alone', ...
    'bottleneck alone','bottleneck | hotspot'};

%% ----------------------- Locked finalized parameters --------------------
P = struct();

% Reference grid and crista-inspired geometry.
P.nx = 60;
P.ny = 36;
P.aOut = 1.12;
P.bOut = 0.62;
P.aIn = 0.94;
P.bIn = 0.49;
P.ellG = 0.085;
P.ellS = 0.042;
P.connMatrix = 1.00;
P.connLumen = 0.82;
P.connNeck = 0.18;

% Reaction and diffusion fields, min-based units.
P.D0_Fe2 = 3.0e-4;
P.D0_Fe3 = 0.70*P.D0_Fe2;
P.c1_0 = 6.0e-3;
P.etaS = 3.40;
P.c2_0 = 1.8e-2;

% Observation and exact first-moment matching.
P.Tobs = 60;
P.dtRecord = 0.50;
P.tvec = 0:P.dtRecord:P.Tobs;
P.dtMoment = 0.05;
P.nBisect = 14;

% Exact locked reference initial state.
P.Fe2_ref_per_voxel = 4;
P.Fe3_ref_per_voxel = 2;

% Baseline local-analysis choices.
P.qRef = 0.85;
P.windowRadiusRef = 0.10;
P.thresholdMultiplierRef = 1.00;

% Required robustness values.
qList = [0.80 0.85 0.90];
windowRadiusList = [0.08 0.10 0.12];
thresholdMultiplierList = [0.80 1.00 1.20 1.40];
gridList = [48 29; 60 36; 72 43];
gridLabels = {'coarse 48x29','reference 60x36','fine 72x43'};

assert(any(abs(qList-P.qRef)<1e-12),'qList must contain q=0.85.');
assert(any(abs(windowRadiusList-P.windowRadiusRef)<1e-12), ...
    'windowRadiusList must contain r_W=0.10.');
assert(any(abs(thresholdMultiplierList-P.thresholdMultiplierRef)<1e-12), ...
    'thresholdMultiplierList must contain multiplier 1.0.');
assert(any(gridList(:,1)==P.nx & gridList(:,2)==P.ny), ...
    'gridList must contain the 60x36 reference grid.');

%% --------------------------- Parallel pool ------------------------------
usePar = configureParallel(USE_PARALLEL,AUTO_START_POOL,N_WORKERS);

%% =======================================================================
% PART I. Build the exact finalized 60x36 reference case
% ========================================================================
fprintf('============================================================\n');
fprintf('Four targeted SD-RDME robustness analyses (%s)\n',upper(RUN_MODE));
fprintf('FPT trajectories/model ensemble = %d\n',nRepFPT);
fprintf('Bootstrap replicates            = %d\n',nBoot);
fprintf('============================================================\n');

Base = buildMatchedCase(P.nx,P.ny,P, ...
    'locked-reference',seedBase+100);

assert(Base.K==868, ...
    ['The reference geometry no longer has K=868 voxels. ' ...
     'The robustness analysis must use the finalized Fig. 3 geometry.']);
assert(all(Base.x1_0==P.Fe2_ref_per_voxel) && ...
       all(Base.x2_0==P.Fe3_ref_per_voxel), ...
    'The reference grid must retain the exact 4/2 initial iron state.');

fprintf('\nReference case audit\n');
fprintf('K                               = %d\n',Base.K);
fprintf('Initial Fe2+/Fe3+ per voxel     = %d / %d\n', ...
    Base.x1_0(1),Base.x2_0(1));
fprintf('Full target at 60 min           = %.6f events\n',Base.targetGlobal);
fprintf('Maximum exact-moment mismatch   = %.5f %%\n', ...
    100*Base.globalMismatchMax);
fprintf('Matching alpha U/H/B            = %.6f / %.6f / %.6f\n', ...
    Base.M(1).alpha,Base.M(2).alpha,Base.M(3).alpha);

%% =======================================================================
% PART II. Construct all reference-grid robustness definitions
% ========================================================================
% The reference c1 field is fixed; only region construction and threshold
% definition change.
nQ = numel(qList);
nR = numel(windowRadiusList);
nT = numel(thresholdMultiplierList);

Aq = false(Base.K,nQ);
Mq = zeros(1,nQ);
qRegionFraction = zeros(1,nQ);
qLocalBurdenRatio = zeros(1,nQ);

for iq=1:nQ
    Def = defineRegionAndBurden(Base,P.windowRadiusRef,qList(iq));
    Aq(:,iq)=Def.A;
    Mq(iq)=Def.M_A;
    qRegionFraction(iq)=Def.nA/Base.K;
    qLocalBurdenRatio(iq)=Def.localBurdenRatioFU;
end

Ar = false(Base.K,nR);
Mr = zeros(1,nR);
radiusRegionFraction = zeros(1,nR);
radiusLocalBurdenRatio = zeros(1,nR);
radiusDeltaHB = zeros(1,nR);

for ir=1:nR
    Def = defineRegionAndBurden(Base,windowRadiusList(ir),P.qRef);
    Ar(:,ir)=Def.A;
    Mr(ir)=Def.M_A;
    radiusRegionFraction(ir)=Def.nA/Base.K;
    radiusLocalBurdenRatio(ir)=Def.localBurdenRatioFU;
    radiusDeltaHB(ir)=Def.deltaHB;
end

DefRef = defineRegionAndBurden(Base,P.windowRadiusRef,P.qRef);
At = repmat(DefRef.A,1,nT);
Mt = max(1,round(DefRef.M_A*thresholdMultiplierList));

% One exact multi-region ensemble per U/H/B/F model records all reference
% q-, radius-, and threshold-robustness endpoints.
Aall = [Aq Ar At];
Mall = [Mq Mr Mt];

idxQ = 1:nQ;
idxR = nQ+(1:nR);
idxT = nQ+nR+(1:nT);

FPTBase = cell(4,1);
baseEnsembleSeed = seedBase+100000;

fprintf('\nRunning reference-grid multi-region exact SSA ensembles...\n');
for m=1:4
    FPTBase{m}=runMultiRegionFPTensembleNSM( ...
        Base.M(m).c1,Base.c2vec,Base.M(m).out,Aall,Mall, ...
        Base.x1_0,Base.x2_0,P.Tobs,nRepFPT, ...
        baseEnsembleSeed,usePar);
end

% Identical definitions appear in q, radius, and threshold groups. They must
% produce exactly identical FPT columns because they were recorded from the
% same trajectories.
idxQRef = find(abs(qList-P.qRef)<1e-12,1);
idxRRef = find(abs(windowRadiusList-P.windowRadiusRef)<1e-12,1);
idxTRef = find(abs(thresholdMultiplierList-1)<1e-12,1);

for m=1:4
    assert(isequal(FPTBase{m}(:,idxQ(idxQRef)), ...
                   FPTBase{m}(:,idxR(idxRRef))) && ...
           isequal(FPTBase{m}(:,idxQ(idxQRef)), ...
                   FPTBase{m}(:,idxT(idxTRef))), ...
        'Duplicated reference FPT definitions are not identical.');
end

SummaryQ = summarizeFPTContrasts( ...
    cellfun(@(x)x(:,idxQ),FPTBase,'UniformOutput',false), ...
    P.Tobs,nBoot,alphaCI,seedBase+200000);

SummaryR = summarizeFPTContrasts( ...
    cellfun(@(x)x(:,idxR),FPTBase,'UniformOutput',false), ...
    P.Tobs,nBoot,alphaCI,seedBase+210000);

SummaryT = summarizeFPTContrasts( ...
    cellfun(@(x)x(:,idxT),FPTBase,'UniformOutput',false), ...
    P.Tobs,nBoot,alphaCI,seedBase+220000);

%% =======================================================================
% PART III. Grid-resolution robustness
% ========================================================================
nG = size(gridList,1);
GridCase = cell(nG,1);
GridDef = cell(nG,1);
GridSummary = cell(nG,1);

gridK = zeros(nG,1);
gridDx = zeros(nG,1);
gridDy = zeros(nG,1);
gridVoxelArea = zeros(nG,1);
gridFe2Total = zeros(nG,1);
gridFe3Total = zeros(nG,1);
gridMeanFe2 = zeros(nG,1);
gridMeanFe3 = zeros(nG,1);
gridMismatch = zeros(nG,1);
gridDeltaHB = zeros(nG,1);
gridLocalBurdenRatio = zeros(nG,1);
gridRegionFraction = zeros(nG,1);
gridThresholdDensity = zeros(nG,1);

gridEffect = zeros(nG,4);
gridEffectL = zeros(nG,4);
gridEffectU = zeros(nG,4);
gridP60 = zeros(nG,4);

idxGridRef = find(gridList(:,1)==P.nx & gridList(:,2)==P.ny,1);

for ig=1:nG
    nxg=gridList(ig,1);
    nyg=gridList(ig,2);

    if ig==idxGridRef
        GCase=Base;
        GDef=DefRef;
        GSummary=extractSingleConditionSummary(SummaryQ,idxQRef);
    else
        fprintf('\nBuilding and simulating grid %dx%d...\n',nxg,nyg);
        GCase=buildMatchedCase(nxg,nyg,P, ...
            'density-preserving',seedBase+300+ig);
        GDef=defineRegionAndBurden( ...
            GCase,P.windowRadiusRef,P.qRef);

        FPTGrid=cell(4,1);
        gridSeed=seedBase+300000+1000*ig;
        for m=1:4
            FPTGrid{m}=runMultiRegionFPTensembleNSM( ...
                GCase.M(m).c1,GCase.c2vec,GCase.M(m).out, ...
                GDef.A,GDef.M_A,GCase.x1_0,GCase.x2_0, ...
                P.Tobs,nRepFPT,gridSeed,usePar);
        end
        GSummary=summarizeFPTContrasts( ...
            FPTGrid,P.Tobs,nBoot,alphaCI,seedBase+310000+1000*ig);
    end

    GridCase{ig}=GCase;
    GridDef{ig}=GDef;
    GridSummary{ig}=GSummary;

    gridK(ig)=GCase.K;
    gridDx(ig)=GCase.dx;
    gridDy(ig)=GCase.dy;
    gridVoxelArea(ig)=GCase.dx*GCase.dy;
    gridFe2Total(ig)=sum(GCase.x1_0);
    gridFe3Total(ig)=sum(GCase.x2_0);
    gridMeanFe2(ig)=mean(GCase.x1_0);
    gridMeanFe3(ig)=mean(GCase.x2_0);
    gridMismatch(ig)=GCase.globalMismatchMax;
    gridDeltaHB(ig)=GDef.deltaHB;
    gridLocalBurdenRatio(ig)=GDef.localBurdenRatioFU;
    gridRegionFraction(ig)=GDef.nA/GCase.K;
    gridThresholdDensity(ig)=GDef.M_A/GDef.nA;

    gridEffect(ig,:)=GSummary.point(1,:);
    gridEffectL(ig,:)=GSummary.ciL(1,:);
    gridEffectU(ig,:)=GSummary.ciU(1,:);
    gridP60(ig,:)=GSummary.P60(1,:);
end

%% =======================================================================
% PART IV. Assemble long-form audit table
% ========================================================================
rows = {};

for iq=1:nQ
    rows(end+1,:)=makeAuditRow( ...
        'Region quantile',sprintf('q=%.2f',qList(iq)),qList(iq), ...
        P.nx,P.ny,Base.K,P.windowRadiusRef,qList(iq),1.0, ...
        qRegionFraction(iq),Mq(iq),Mq(iq)/(qRegionFraction(iq)*Base.K), ...
        Base.globalMismatchMax,DefRef.deltaHB,qLocalBurdenRatio(iq), ...
        SummaryQ,iq); %#ok<AGROW>
end

for ir=1:nR
    rows(end+1,:)=makeAuditRow( ...
        'Window radius',sprintf('r=%.2f',windowRadiusList(ir)), ...
        windowRadiusList(ir),P.nx,P.ny,Base.K, ...
        windowRadiusList(ir),P.qRef,1.0, ...
        radiusRegionFraction(ir),Mr(ir), ...
        Mr(ir)/(radiusRegionFraction(ir)*Base.K), ...
        Base.globalMismatchMax,radiusDeltaHB(ir), ...
        radiusLocalBurdenRatio(ir),SummaryR,ir); %#ok<AGROW>
end

for it=1:nT
    rows(end+1,:)=makeAuditRow( ...
        'Threshold multiplier', ...
        sprintf('M=%.1fx',thresholdMultiplierList(it)), ...
        thresholdMultiplierList(it),P.nx,P.ny,Base.K, ...
        P.windowRadiusRef,P.qRef,thresholdMultiplierList(it), ...
        DefRef.nA/Base.K,Mt(it),Mt(it)/DefRef.nA, ...
        Base.globalMismatchMax,DefRef.deltaHB, ...
        DefRef.localBurdenRatioFU,SummaryT,it); %#ok<AGROW>
end

for ig=1:nG
    GS=GridSummary{ig};
    GD=GridDef{ig};
    GC=GridCase{ig};
    rows(end+1,:)=makeAuditRow( ...
        'Grid resolution',gridLabels{ig},ig, ...
        gridList(ig,1),gridList(ig,2),GC.K, ...
        P.windowRadiusRef,P.qRef,1.0,GD.nA/GC.K, ...
        GD.M_A,GD.M_A/GD.nA,GC.globalMismatchMax, ...
        GD.deltaHB,GD.localBurdenRatioFU,GS,1); %#ok<AGROW>
end

auditTable=cell2table(rows,'VariableNames',{ ...
    'RobustnessType','ConditionLabel','ConditionValue', ...
    'Nx','Ny','K','WindowRadius','RegionQuantile', ...
    'ThresholdMultiplier','RegionFraction','RegionalThreshold', ...
    'ThresholdDensity','MaximumGlobalMismatch', ...
    'DeltaBurdenHeterogeneity','LocalBurdenRatio_FU', ...
    'DeltaRMST_UF','DeltaRMST_UF_CI_L','DeltaRMST_UF_CI_U', ...
    'EffectH','EffectH_CI_L','EffectH_CI_U', ...
    'EffectB','EffectB_CI_L','EffectB_CI_U', ...
    'EffectBgivenH','EffectBgivenH_CI_L','EffectBgivenH_CI_U', ...
    'P60_U','P60_H','P60_B','P60_F'});

gridAuditTable=table( ...
    string(gridLabels(:)),gridList(:,1),gridList(:,2),gridK, ...
    gridDx,gridDy,gridVoxelArea,gridFe2Total,gridFe3Total, ...
    gridMeanFe2,gridMeanFe3,gridMismatch,gridRegionFraction, ...
    gridThresholdDensity,gridDeltaHB,gridLocalBurdenRatio, ...
    'VariableNames',{ ...
    'GridLabel','Nx','Ny','K','Dx','Dy','VoxelArea', ...
    'InitialFe2Total','InitialFe3Total','MeanFe2PerVoxel', ...
    'MeanFe3PerVoxel','MaximumGlobalMismatch','RegionFraction', ...
    'ThresholdDensity','DeltaBurdenHeterogeneity', ...
    'LocalBurdenRatio_FU'});

%% =======================================================================
% PART V. Phenotype-centered four-panel SI robustness figure
% ========================================================================
% The primary figure displays the main-text phenomena:
%   local burden redistribution and Full-vs-Uniform local FPT advancement.
% U/H/B/F causal effects are retained in the audit table and in one compact
% secondary mechanism-audit figure below.

fig=figure('Units','centimeters','Position',[1 1 18.8 11.6], ...
    'Color','w','Renderer','painters');
tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

% A. Region quantile: local burden enrichment + Full-vs-Uniform FPT.
ax=nexttile(tl,1);
[hBurdenA,hFPTA]=plotBurdenAndFPTDualAxis( ...
    ax,qList,qLocalBurdenRatio,1, ...
    SummaryQ.point(:,1),SummaryQ.ciL(:,1),SummaryQ.ciU(:,1), ...
    colBurden,colFPT, ...
    'Local burden ratio, B_A^F/B_A^U', ...
    '\DeltaRMST_{60}, U-F (min)');
xline(ax,P.qRef,':','Color',colReference, ...
    'LineWidth',0.85,'HandleVisibility','off');
xlabel(ax,'Region quantile, q');
title(ax,'Region-definition robustness', ...
    'FontWeight','normal','FontSize',8.4);
panelLabelLocal(ax,'A');
legend(ax,[hBurdenA hFPTA], ...
    {'local ROS burden','local FPT advancement'}, ...
    'Location','northwest','Box','off','FontSize',5.8);
text(ax,0.98,0.06,sprintf( ...
    'burden > 1: %d/%d; FPT > 0: %d/%d', ...
    nnz(qLocalBurdenRatio>1),nQ,nnz(SummaryQ.point(:,1)>0),nQ), ...
    'Units','normalized','HorizontalAlignment','right','FontSize',5.5);

% B. Local-window radius: burden heterogeneity + Full-vs-Uniform FPT.
ax=nexttile(tl,2);
[hBurdenB,hFPTB]=plotBurdenAndFPTDualAxis( ...
    ax,windowRadiusList,radiusDeltaHB,0, ...
    SummaryR.point(:,1),SummaryR.ciL(:,1),SummaryR.ciU(:,1), ...
    colBurden,colFPT, ...
    '\DeltaH_B = H_B^F-H_B^U', ...
    '\DeltaRMST_{60}, U-F (min)');
xline(ax,P.windowRadiusRef,':','Color',colReference, ...
    'LineWidth',0.85,'HandleVisibility','off');
xlabel(ax,'Local-window radius, r_W');
title(ax,'Local-scale robustness', ...
    'FontWeight','normal','FontSize',8.4);
panelLabelLocal(ax,'B');
legend(ax,[hBurdenB hFPTB], ...
    {'burden heterogeneity','local FPT advancement'}, ...
    'Location','northwest','Box','off','FontSize',5.8);
text(ax,0.98,0.06,sprintf( ...
    '\DeltaH_B > 0: %d/%d; FPT > 0: %d/%d', ...
    nnz(radiusDeltaHB>0),nR,nnz(SummaryR.point(:,1)>0),nR), ...
    'Units','normalized','HorizontalAlignment','right','FontSize',5.5);

% C. Threshold: Full-vs-Uniform RMST effect + crossing probabilities.
ax=nexttile(tl,3);
[hFPTC,hP60U,hP60F]=plotThresholdPhenotypes( ...
    ax,thresholdMultiplierList, ...
    SummaryT.point(:,1),SummaryT.ciL(:,1),SummaryT.ciU(:,1), ...
    SummaryT.P60(:,1),SummaryT.P60(:,4), ...
    colFPT,colP60U,colP60F);
xline(ax,1,':','Color',colReference, ...
    'LineWidth',0.85,'HandleVisibility','off');
xlabel(ax,'Regional threshold multiplier');
title(ax,'FPT-threshold robustness', ...
    'FontWeight','normal','FontSize',8.4);
panelLabelLocal(ax,'C');
legend(ax,[hFPTC hP60U hP60F], ...
    {'\DeltaRMST_{60}, U-F','P_{60}, Uniform','P_{60}, Full'}, ...
    'Location','southwest','Box','off','FontSize',5.7);
text(ax,0.98,0.93,sprintf('FPT effect > 0: %d/%d thresholds', ...
    nnz(SummaryT.point(:,1)>0),nT), ...
    'Units','normalized','HorizontalAlignment','right', ...
    'VerticalAlignment','top','FontSize',5.5);

% D. Grid: burden heterogeneity + Full-vs-Uniform FPT.
ax=nexttile(tl,4);
[hBurdenD,hFPTD]=plotBurdenAndFPTDualAxis( ...
    ax,1:nG,gridDeltaHB,0, ...
    gridEffect(:,1),gridEffectL(:,1),gridEffectU(:,1), ...
    colBurden,colFPT, ...
    '\DeltaH_B = H_B^F-H_B^U', ...
    '\DeltaRMST_{60}, U-F (min)');
xline(ax,idxGridRef,':','Color',colReference, ...
    'LineWidth',0.85,'HandleVisibility','off');
set(ax,'XTick',1:nG,'XTickLabel',gridLabels);
xtickangle(ax,18);
xlabel(ax,'Spatial grid');
title(ax,'Grid-resolution robustness', ...
    'FontWeight','normal','FontSize',8.4);
panelLabelLocal(ax,'D');
legend(ax,[hBurdenD hFPTD], ...
    {'burden heterogeneity','local FPT advancement'}, ...
    'Location','northwest','Box','off','FontSize',5.8);
text(ax,0.98,0.06,sprintf( ...
    '\DeltaH_B > 0: %d/%d; FPT > 0: %d/%d; mismatch %.2f%%', ...
    nnz(gridDeltaHB>0),nG,nnz(gridEffect(:,1)>0),nG, ...
    100*max([Base.globalMismatchMax;gridMismatch])), ...
    'Units','normalized','HorizontalAlignment','right','FontSize',5.3);

if SAVE_FIG
    exportgraphics(fig,[outFigure '.pdf'],'ContentType','vector');
    exportgraphics(fig,[outFigure '.png'],'Resolution',600);
    exportgraphics(fig,[outFigure '.tif'],'Resolution',600);
    savefig(fig,[outFigure '.fig']);
end

%% ------------------ Compact U/H/B/F mechanism audit --------------------
% Deduplicate the identical reference condition: retain all q conditions,
% then add only non-reference radius, threshold and grid conditions.
idxRNonRef=find(abs(windowRadiusList-P.windowRadiusRef)>1e-12);
idxTNonRef=find(abs(thresholdMultiplierList-1)>1e-12);
idxGNonRef=setdiff(1:nG,idxGridRef);

mechanismEffects=[ ...
    SummaryQ.point(:,2:4); ...
    SummaryR.point(idxRNonRef,2:4); ...
    SummaryT.point(idxTNonRef,2:4); ...
    gridEffect(idxGNonRef,2:4)];

mechanismGroup=[ ...
    repmat("region q",nQ,1); ...
    repmat("window radius",numel(idxRNonRef),1); ...
    repmat("threshold",numel(idxTNonRef),1); ...
    repmat("grid",numel(idxGNonRef),1)];

mechanismAuditTable=table( ...
    mechanismGroup,mechanismEffects(:,1),mechanismEffects(:,2), ...
    mechanismEffects(:,3), ...
    'VariableNames',{'RobustnessAxis','EffectH','EffectB', ...
    'EffectBgivenH'});

figM=figure('Units','centimeters','Position',[1 1 9.4 7.2], ...
    'Color','w','Renderer','painters');
axM=axes(figM); hold(axM,'on');
yline(axM,0,':','Color',[0.52 0.52 0.52], ...
    'LineWidth',0.8,'HandleVisibility','off');

mechanismColors=[colH;colB;colBH];
mechanismNames={'hotspot alone','bottleneck alone','bottleneck | hotspot'};

for e=1:3
    nPoint=size(mechanismEffects,1);
    jitter=linspace(-0.18,0.18,nPoint).';
    lightColor=lightenColorLocal(mechanismColors(e,:),0.55);

    scatter(axM,e+jitter,mechanismEffects(:,e),20,lightColor, ...
        'filled','MarkerEdgeColor','none');

    med=median(mechanismEffects(:,e));
    plot(axM,[e-0.18 e+0.18],[med med],'-', ...
        'Color',mechanismColors(e,:),'LineWidth',1.7);
    plot(axM,e,med,'o','Color',mechanismColors(e,:), ...
        'MarkerFaceColor','w','MarkerSize',5.2,'LineWidth',1.1);
end

set(axM,'XTick',1:3,'XTickLabel',mechanismNames);
xtickangle(axM,16);
xlim(axM,[0.5 3.5]);
ylabel(axM,'RMST advancement (min)');
title(axM,'U/H/B/F ordering across all robustness choices', ...
    'FontWeight','normal','FontSize',8.3);
text(axM,0.02,0.97,'Audit','Units','normalized', ...
    'FontWeight','bold','FontSize',9.5,'VerticalAlignment','top');
formatAxes(axM);

if SAVE_FIG
    exportgraphics(figM,[outMechanismFigure '.pdf'],'ContentType','vector');
    exportgraphics(figM,[outMechanismFigure '.png'],'Resolution',600);
    exportgraphics(figM,[outMechanismFigure '.tif'],'Resolution',600);
    savefig(figM,[outMechanismFigure '.fig']);
end

%% ------------------------------ Save data -------------------------------
if SAVE_DATA
    writetable(auditTable,outAuditCSV);
    writetable(gridAuditTable,outGridCSV);

    metadata=struct();
    metadata.RUN_MODE=RUN_MODE;
    metadata.nRepFPT=nRepFPT;
    metadata.nBoot=nBoot;
    metadata.parameters=P;
    metadata.qList=qList;
    metadata.windowRadiusList=windowRadiusList;
    metadata.thresholdMultiplierList=thresholdMultiplierList;
    metadata.gridList=gridList;
    metadata.effectLabels=effectLabels;
    metadata.referenceAlpha=[Base.M.alpha];

    save(outMAT,'metadata','auditTable','gridAuditTable', ...
        'mechanismAuditTable','SummaryQ','SummaryR','SummaryT', ...
        'gridEffect','gridEffectL','gridEffectU','gridP60', ...
        'qLocalBurdenRatio','radiusDeltaHB','gridDeltaHB', ...
        '-v7.3');
end

fprintf('\n================ Final robustness summary ================\n');
fprintf('q robustness: Delta U-F positive              = %d/%d\n', ...
    nnz(SummaryQ.point(:,1)>0),nQ);
fprintf('radius robustness: Delta U-F positive         = %d/%d\n', ...
    nnz(SummaryR.point(:,1)>0),nR);
fprintf('threshold robustness: Delta U-F positive      = %d/%d\n', ...
    nnz(SummaryT.point(:,1)>0),nT);
fprintf('grid robustness: Delta U-F positive           = %d/%d\n', ...
    nnz(gridEffect(:,1)>0),nG);
fprintf('radius robustness: Delta H_B positive         = %d/%d\n', ...
    nnz(radiusDeltaHB>0),nR);
fprintf('grid robustness: Delta H_B positive           = %d/%d\n', ...
    nnz(gridDeltaHB>0),nG);
fprintf('maximum exact-moment mismatch                 = %.5f %%\n', ...
    100*max([Base.globalMismatchMax;gridMismatch]));
fprintf('==========================================================\n');

fprintf('\nCompleted. Outputs:\n');
fprintf('  %s.{pdf,png,tif,fig}\n',outFigure);
fprintf('  %s.{pdf,png,tif,fig}\n',outMechanismFigure);
fprintf('  %s\n',outAuditCSV);
fprintf('  %s\n',outGridCSV);
fprintf('  %s\n',outMAT);

%% ======================== Analysis-specific functions ===================

function Case = buildMatchedCase(nx,ny,P,initialMode,allocationSeed)
    geom=buildCurrentCristaGeometry(nx,ny, ...
        P.aOut,P.bOut,P.aIn,P.bIn,P.ellG,P.ellS, ...
        P.connMatrix,P.connLumen,P.connNeck);

    mask=geom.OmegaR;
    K=nnz(mask);
    idxMap=buildMaskIndexMap(mask);

    x=geom.x;
    y=geom.y;
    Xg=geom.X;
    Yg=geom.Y;
    S=geom.S;
    G=geom.G;
    dx=mean(diff(x));
    dy=mean(diff(y));

    D1Full=P.D0_Fe2*G;
    D2Full=P.D0_Fe3*G;
    D1Full(~mask)=NaN;
    D2Full(~mask)=NaN;

    c1Full=P.c1_0*(1+P.etaS*S);
    c1Full(~mask)=NaN;

    c1F=c1Full(mask);
    D1F=D1Full(mask);
    D2F=D2Full(mask);

    c1Mean=mean(c1F);
    D1Mean=mean(D1F);
    D2Mean=mean(D2F);
    c2vec=P.c2_0*ones(K,1);

    [x1_0,x2_0,initAudit]=makeInitialIronState( ...
        K,dx,dy,P,initialMode,allocationSeed);

    [fromHet1,toHet1,qHet1]=buildDirectedEdgesHarmonic( ...
        mask,idxMap,D1Full,dx,dy);
    [fromHet2,toHet2,qHet2]=buildDirectedEdgesHarmonic( ...
        mask,idxMap,D2Full,dx,dy);

    D1Uniform=NaN(size(D1Full));
    D2Uniform=NaN(size(D2Full));
    D1Uniform(mask)=D1Mean;
    D2Uniform(mask)=D2Mean;

    [fromUni1,toUni1,qUni1]=buildDirectedEdgesHarmonic( ...
        mask,idxMap,D1Uniform,dx,dy);
    [fromUni2,toUni2,qUni2]=buildDirectedEdgesHarmonic( ...
        mask,idxMap,D2Uniform,dx,dy);

    modelTemplate=struct( ...
        'name','','c1Base',[],'c1',[],'alpha',NaN, ...
        'from1',[],'to1',[],'q1',[], ...
        'from2',[],'to2',[],'q2',[], ...
        'moment',[],'out',[],'momentRatio',NaN);
    M=repmat(modelTemplate,1,4);

    M(1).name='U';
    M(1).c1Base=c1Mean*ones(K,1);
    M(1).from1=fromUni1; M(1).to1=toUni1; M(1).q1=qUni1;
    M(1).from2=fromUni2; M(1).to2=toUni2; M(1).q2=qUni2;

    M(2).name='H';
    M(2).c1Base=c1F;
    M(2).from1=fromUni1; M(2).to1=toUni1; M(2).q1=qUni1;
    M(2).from2=fromUni2; M(2).to2=toUni2; M(2).q2=qUni2;

    M(3).name='B';
    M(3).c1Base=c1Mean*ones(K,1);
    M(3).from1=fromHet1; M(3).to1=toHet1; M(3).q1=qHet1;
    M(3).from2=fromHet2; M(3).to2=toHet2; M(3).q2=qHet2;

    M(4).name='F';
    M(4).c1Base=c1F;
    M(4).from1=fromHet1; M(4).to1=toHet1; M(4).q1=qHet1;
    M(4).from2=fromHet2; M(4).to2=toHet2; M(4).q2=qHet2;

    M(4).alpha=1;
    M(4).moment=momentEulerAll( ...
        M(4).c1Base,c2vec, ...
        M(4).from1,M(4).to1,M(4).q1, ...
        M(4).from2,M(4).to2,M(4).q2, ...
        x1_0,x2_0,P.tvec,P.dtMoment);

    targetGlobal=M(4).moment.globalM3(end);

    for m=1:3
        [M(m).alpha,M(m).moment]=matchModelToTarget( ...
            M(m).c1Base,c2vec, ...
            M(m).from1,M(m).to1,M(m).q1, ...
            M(m).from2,M(m).to2,M(m).q2, ...
            x1_0,x2_0,P.tvec,P.dtMoment, ...
            targetGlobal,P.nBisect);
    end

    for m=1:4
        M(m).c1=M(m).alpha*M(m).c1Base;
        M(m).out=buildOutgoingLists(K, ...
            M(m).from1,M(m).to1,M(m).q1, ...
            M(m).from2,M(m).to2,M(m).q2);
        M(m).momentRatio=M(m).moment.globalM3(end)/targetGlobal;
    end

    Case=struct();
    Case.nx=nx; Case.ny=ny;
    Case.geom=geom; Case.mask=mask; Case.K=K;
    Case.x=x; Case.y=y; Case.X=Xg; Case.Y=Yg;
    Case.S=S; Case.G=G; Case.c1Full=c1Full;
    Case.c1F=c1F; Case.c2vec=c2vec;
    Case.dx=dx; Case.dy=dy;
    Case.x1_0=x1_0; Case.x2_0=x2_0;
    Case.initialAudit=initAudit;
    Case.M=M;
    Case.targetGlobal=targetGlobal;
    Case.globalMismatchMax=max(abs([M.momentRatio]-1));
    Case.mStar=M(1).moment.globalM3(end)/K;
end

function [x1,x2,audit] = makeInitialIronState( ...
    K,dx,dy,P,mode,seed)

    baseDx=2.46/(P.nx-1);
    baseDy=1.40/(P.ny-1);

    switch lower(mode)
        case 'locked-reference'
            x1=P.Fe2_ref_per_voxel*ones(K,1);
            x2=P.Fe3_ref_per_voxel*ones(K,1);

        case 'density-preserving'
            densityFe2=P.Fe2_ref_per_voxel/(baseDx*baseDy);
            densityFe3=P.Fe3_ref_per_voxel/(baseDx*baseDy);

            totalFe2=max(1,round(densityFe2*K*dx*dy));
            totalFe3=max(1,round(densityFe3*K*dx*dy));

            rng(seed,'twister');
            perm=randperm(K);
            x1=balancedIntegerAllocation(totalFe2,K,perm);
            x2=balancedIntegerAllocation(totalFe3,K,perm);

        otherwise
            error('Unknown initial-state mode: %s',mode);
    end

    audit=struct();
    audit.totalFe2=sum(x1);
    audit.totalFe3=sum(x2);
    audit.meanFe2=mean(x1);
    audit.meanFe3=mean(x2);
    audit.dx=dx;
    audit.dy=dy;
    audit.voxelArea=dx*dy;
end

function x = balancedIntegerAllocation(totalCount,K,perm)
    base=floor(totalCount/K);
    remainder=totalCount-base*K;
    x=base*ones(K,1);
    if remainder>0
        x(perm(1:remainder))=x(perm(1:remainder))+1;
    end
end

function [alpha,momentOut] = matchModelToTarget( ...
    c1Base,c2,from1,to1,q1,from2,to2,q2, ...
    x1_0,x2_0,tvec,dtMoment,targetGlobal,nBisect)

    lo=0.10;
    hi=7.00;

    momLo=momentEulerAll(lo*c1Base,c2, ...
        from1,to1,q1,from2,to2,q2, ...
        x1_0,x2_0,tvec,dtMoment);
    momHi=momentEulerAll(hi*c1Base,c2, ...
        from1,to1,q1,from2,to2,q2, ...
        x1_0,x2_0,tvec,dtMoment);

    while momLo.globalM3(end)>targetGlobal
        hi=lo;
        lo=0.5*lo;
        momLo=momentEulerAll(lo*c1Base,c2, ...
            from1,to1,q1,from2,to2,q2, ...
            x1_0,x2_0,tvec,dtMoment);
    end

    while momHi.globalM3(end)<targetGlobal
        lo=hi;
        hi=2*hi;
        momHi=momentEulerAll(hi*c1Base,c2, ...
            from1,to1,q1,from2,to2,q2, ...
            x1_0,x2_0,tvec,dtMoment);
        if hi>100
            error('Could not bracket the Full global target.');
        end
    end

    for it=1:nBisect
        mid=0.5*(lo+hi);
        trial=momentEulerAll(mid*c1Base,c2, ...
            from1,to1,q1,from2,to2,q2, ...
            x1_0,x2_0,tvec,dtMoment);

        if trial.globalM3(end)<targetGlobal
            lo=mid;
        else
            hi=mid;
        end
    end

    alpha=0.5*(lo+hi);
    momentOut=momentEulerAll(alpha*c1Base,c2, ...
        from1,to1,q1,from2,to2,q2, ...
        x1_0,x2_0,tvec,dtMoment);
end

function Def = defineRegionAndBurden(Case,radius,q)
    coordX=Case.X(Case.mask);
    coordY=Case.Y(Case.mask);
    win=buildCircularWindowsSparse(coordX,coordY,radius);

    c1Center=normalizeVector(Case.c1F);
    c1Window=(win.W*c1Center)./win.size;
    c1Window=normalizeVector(c1Window);

    A=topFractionMask(c1Window,1-q);
    nA=nnz(A);
    M_A=max(1,round(Case.mStar*nA));

    burdenU=(win.W*Case.M(1).moment.finalM3)./win.size;
    burdenF=(win.W*Case.M(4).moment.finalM3)./win.size;

    HBU=normalizedP90P10(burdenU);
    HBF=normalizedP90P10(burdenF);

    Def=struct();
    Def.radius=radius;
    Def.q=q;
    Def.win=win;
    Def.A=A;
    Def.nA=nA;
    Def.M_A=M_A;
    Def.HBU=HBU;
    Def.HBF=HBF;
    Def.deltaHB=HBF-HBU;
    Def.localBurdenRatioFU= ...
        mean(Case.M(4).moment.finalM3(A))/ ...
        max(mean(Case.M(1).moment.finalM3(A)),eps);
end

function H = normalizedP90P10(x)
    x=x(:);
    H=(quantileLocal(x,0.90)-quantileLocal(x,0.10))/ ...
        max(mean(x),eps);
end

function win = buildCircularWindowsSparse(coordX,coordY,radius)
    coordX=coordX(:);
    coordY=coordY(:);
    K=numel(coordX);

    I=cell(K,1);
    J=cell(K,1);
    radius2=radius^2+1e-14;

    for k=1:K
        idx=find((coordX-coordX(k)).^2+ ...
            (coordY-coordY(k)).^2<=radius2);
        I{k}=repmat(k,numel(idx),1);
        J{k}=idx(:);
    end

    I=vertcat(I{:});
    J=vertcat(J{:});

    win.W=sparse(I,J,1,K,K);
    win.size=full(sum(win.W,2));

    if any(win.size<1)
        error('At least one local window is empty.');
    end
end

function v = normalizeVector(v)
    v=v(:);
    v=(v-min(v))/(max(v)-min(v)+eps);
end

function maskTop = topFractionMask(score,topFraction)
    score=score(:);
    K=numel(score);
    nTop=max(1,ceil(topFraction*K));

    [~,ord]=sortrows([-score,(1:K).'],[1 2]);
    maskTop=false(K,1);
    maskTop(ord(1:nTop))=true;
end

function FPT = runMultiRegionFPTensembleNSM( ...
    c1,c2,outgoing,Aset,thresholds,x1_0,x2_0, ...
    Tobs,nRep,seedBase,usePar)

    Aset=logical(Aset);
    if isvector(Aset)
        Aset=Aset(:);
    end

    thresholds=thresholds(:).';
    nCond=size(Aset,2);

    assert(size(Aset,1)==numel(c1), ...
        'Aset must have one row per voxel.');
    assert(numel(thresholds)==nCond, ...
        'One threshold is required per region condition.');

    FPT=inf(nRep,nCond);

    if usePar
        parfor r=1:nRep
            rng(seedBase+r,'twister');
            FPT(r,:)=simulateMultiRegionFPTOneNSM( ...
                c1,c2,outgoing,Aset,thresholds, ...
                x1_0,x2_0,Tobs);
        end
    else
        for r=1:nRep
            rng(seedBase+r,'twister');
            FPT(r,:)=simulateMultiRegionFPTOneNSM( ...
                c1,c2,outgoing,Aset,thresholds, ...
                x1_0,x2_0,Tobs);
        end
    end
end

function fpt = simulateMultiRegionFPTOneNSM( ...
    c1,c2,out,Aset,thresholds,x1_0,x2_0,Tobs)

    x1=x1_0(:);
    x2=x2_0(:);

    nCond=size(Aset,2);
    localROS=zeros(1,nCond);
    fpt=inf(1,nCond);
    active=true(1,nCond);

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

        aR1=c1(k)*x1(k);
        aR2=c2(k)*x2(k);
        aD1=out.q1sum(k)*x1(k);
        aD2=out.q2sum(k)*x2(k);
        aLocal=aR1+aR2+aD1+aD2;

        if aLocal<=0
            oldRate=voxelRate(k);
            voxelRate(k)=0;
            tree=fenwickAdd(tree,k,-oldRate);
            t=tNext;
            continue;
        end

        u=rand*aLocal;
        touched=k;

        if u<=aR1
            x1(k)=x1(k)-1;
            x2(k)=x2(k)+1;

            hit=active & Aset(k,:);
            localROS(hit)=localROS(hit)+1;

            crossed=hit & (localROS>=thresholds);
            fpt(crossed)=tNext;
            active(crossed)=false;

        elseif u<=aR1+aR2
            x2(k)=x2(k)-1;
            x1(k)=x1(k)+1;

        elseif u<=aR1+aR2+aD1
            j=chooseDestination(out.to1{k},out.q1{k});
            x1(k)=x1(k)-1;
            x1(j)=x1(j)+1;
            touched=[k;j]; %#ok<AGROW>

        else
            j=chooseDestination(out.to2{k},out.q2{k});
            x2(k)=x2(k)-1;
            x2(j)=x2(j)+1;
            touched=[k;j]; %#ok<AGROW>
        end

        [voxelRate,tree]=refreshVoxelRates( ...
            unique(touched),voxelRate,tree, ...
            x1,x2,c1,c2,out.q1sum,out.q2sum);

        t=tNext;
    end
end

function Summary = summarizeFPTContrasts( ...
    FPT,Tobs,nBoot,alpha,seedBase)

    nCond=size(FPT{1},2);
    point=zeros(nCond,4);
    ciL=zeros(nCond,4);
    ciU=zeros(nCond,4);
    P60=zeros(nCond,4);
    RMST=zeros(nCond,4);

    C=cell(4,1);
    for m=1:4
        C{m}=min(FPT{m},Tobs);
        RMST(:,m)=mean(C{m},1).';
        P60(:,m)=mean(FPT{m}<=Tobs,1).';
    end

    for j=1:nCond
        contrasts=[
            C{1}(:,j)-C{4}(:,j), ...
            C{1}(:,j)-C{2}(:,j), ...
            C{1}(:,j)-C{3}(:,j), ...
            C{2}(:,j)-C{4}(:,j)];

        point(j,:)=mean(contrasts,1);

        for e=1:4
            ci=bootstrapMeanVectorCI( ...
                contrasts(:,e),nBoot,alpha, ...
                seedBase+10000*j+100*e);
            ciL(j,e)=ci(1);
            ciU(j,e)=ci(2);
        end
    end

    Summary=struct();
    Summary.point=point;
    Summary.ciL=ciL;
    Summary.ciU=ciU;
    Summary.P60=P60;
    Summary.RMST=RMST;
end

function ci = bootstrapMeanVectorCI(x,nBoot,alpha,seed)
    x=x(:);
    n=numel(x);
    boot=zeros(nBoot,1);

    rng(seed,'twister');
    for b=1:nBoot
        idx=randi(n,n,1);
        boot(b)=mean(x(idx));
    end

    ci=[quantileLocal(boot,alpha/2), ...
        quantileLocal(boot,1-alpha/2)];
end

function S = extractSingleConditionSummary(Summary,index)
    S=struct();
    S.point=Summary.point(index,:);
    S.ciL=Summary.ciL(index,:);
    S.ciU=Summary.ciU(index,:);
    S.P60=Summary.P60(index,:);
    S.RMST=Summary.RMST(index,:);
end

function row = makeAuditRow( ...
    type,label,value,nx,ny,K,radius,q,thresholdMultiplier, ...
    regionFraction,M_A,thresholdDensity,mismatch,deltaHB, ...
    localRatio,Summary,index)

    row={ ...
        type,label,value,nx,ny,K,radius,q,thresholdMultiplier, ...
        regionFraction,M_A,thresholdDensity,mismatch,deltaHB, ...
        localRatio, ...
        Summary.point(index,1),Summary.ciL(index,1),Summary.ciU(index,1), ...
        Summary.point(index,2),Summary.ciL(index,2),Summary.ciU(index,2), ...
        Summary.point(index,3),Summary.ciL(index,3),Summary.ciU(index,3), ...
        Summary.point(index,4),Summary.ciL(index,4),Summary.ciU(index,4), ...
        Summary.P60(index,1),Summary.P60(index,2), ...
        Summary.P60(index,3),Summary.P60(index,4)};
end

function [hBurden,hFPT] = plotBurdenAndFPTDualAxis( ...
    ax,x,burdenMetric,burdenNull,fpt,fptL,fptU, ...
    burdenColor,fptColor,burdenLabel,fptLabel)

    x=x(:);
    burdenMetric=burdenMetric(:);
    fpt=fpt(:);
    fptL=fptL(:);
    fptU=fptU(:);

    yyaxis(ax,'left');
    hold(ax,'on');
    hBurden=plot(ax,x,burdenMetric,'-o', ...
        'Color',burdenColor,'MarkerFaceColor',burdenColor, ...
        'MarkerEdgeColor','w','MarkerSize',4.4,'LineWidth',1.2);
    yline(ax,burdenNull,':','Color',lightenColorLocal(burdenColor,0.35), ...
        'LineWidth',0.8,'HandleVisibility','off');
    ylabel(ax,burdenLabel);
    ax.YAxis(1).Color=burdenColor;
    ylim(ax,paddedLimitsWithReferenceLocal(burdenMetric,burdenNull,0.14));

    yyaxis(ax,'right');
    hFPT=errorbar(ax,x,fpt,max(fpt-fptL,0),max(fptU-fpt,0), ...
        '-s','Color',fptColor,'MarkerFaceColor','w', ...
        'MarkerEdgeColor',fptColor,'MarkerSize',4.2, ...
        'LineWidth',1.15,'CapSize',3);
    yline(ax,0,':','Color',lightenColorLocal(fptColor,0.40), ...
        'LineWidth',0.8,'HandleVisibility','off');
    ylabel(ax,fptLabel);
    ax.YAxis(2).Color=fptColor;
    ylim(ax,paddedLimitsWithReferenceLocal([fptL;fptU],0,0.12));

    formatAxes(ax);
end

function [hFPT,hP60U,hP60F] = plotThresholdPhenotypes( ...
    ax,x,fpt,fptL,fptU,p60U,p60F,fptColor,p60UColor,p60FColor)

    x=x(:);
    fpt=fpt(:);
    fptL=fptL(:);
    fptU=fptU(:);
    p60U=p60U(:);
    p60F=p60F(:);

    yyaxis(ax,'left');
    hold(ax,'on');
    hFPT=errorbar(ax,x,fpt,max(fpt-fptL,0),max(fptU-fpt,0), ...
        '-s','Color',fptColor,'MarkerFaceColor','w', ...
        'MarkerEdgeColor',fptColor,'MarkerSize',4.3, ...
        'LineWidth',1.18,'CapSize',3);
    yline(ax,0,':','Color',lightenColorLocal(fptColor,0.40), ...
        'LineWidth',0.8,'HandleVisibility','off');
    ylabel(ax,'\DeltaRMST_{60}, U-F (min)');
    ax.YAxis(1).Color=fptColor;
    ylim(ax,paddedLimitsWithReferenceLocal([fptL;fptU],0,0.12));

    yyaxis(ax,'right');
    hP60U=plot(ax,x,p60U,'--o','Color',p60UColor, ...
        'MarkerFaceColor','w','MarkerEdgeColor',p60UColor, ...
        'MarkerSize',4.0,'LineWidth',1.05);
    hP60F=plot(ax,x,p60F,'-d','Color',p60FColor, ...
        'MarkerFaceColor','w','MarkerEdgeColor',p60FColor, ...
        'MarkerSize',4.1,'LineWidth',1.08);
    ylabel(ax,'Crossing probability by 60 min');
    ax.YAxis(2).Color=[0.27 0.27 0.27];
    ylim(ax,[0 1.04]);

    formatAxes(ax);
end

function lim = paddedLimitsWithReferenceLocal(values,reference,padFraction)
    values=values(:);
    values=values(isfinite(values));
    values=[values;reference];

    lo=min(values);
    hi=max(values);
    span=hi-lo;

    if span<1e-10
        span=max(abs(hi),1);
    end

    lim=[lo-padFraction*span,hi+padFraction*span];

    if reference==0
        lim(1)=min(lim(1),-0.04*max(abs(lim(2)),1));
    end
end

function panelLabelLocal(ax,label)
    text(ax,0.02,0.97,label,'Units','normalized', ...
        'FontSize',11.5,'FontWeight','bold','VerticalAlignment','top');
end

function c = lightenColorLocal(c0,amount)
    amount=min(max(amount,0),1);
    c=c0+(1-c0)*amount;
end

%% ================= Reused validated SD-RDME helpers ====================

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
