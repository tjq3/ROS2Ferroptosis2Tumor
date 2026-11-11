%% SDRDME_GeometryRobustness.m
% -------------------------------------------------------------------------
% FINAL SI geometry-generalization analysis for the mitochondrial SD-RDME.
%
% PURPOSE
%   Test whether the load-bearing conclusions of the main-text SD-RDME
%   analysis persist after changing mitochondrial outline and crista
%   architecture, without repeating the entire main-text workflow.
%
% FOUR-PANEL SI FIGURE
%   A. Five representative geometries. Reference is exactly the finalized
%      main-text geometry; alternatives perturb either overall outline or
%      crista density.
%   B. Predictor-independent local-burden heterogeneity:
%         Delta H_B = H_B^F - H_B^U,
%      where H_B is the normalized P90-P10 spread across all fixed-radius
%      local windows. Positive values mean that globally matched SD-RDME
%      retains stronger local burden redistribution.
%   C. Regional FPT advancement in the pre-specified c1-hat top-15% region:
%         Delta RMST_60 = RMST_U - RMST_F.
%   D. Geometry robustness of the U/H/B/F causal ordering:
%         E_H   = RMST_U - RMST_H
%         E_B   = RMST_U - RMST_B
%         E_B|H = RMST_H - RMST_F.
%
% WORKLOAD CONTROL
%   The analysis deliberately does NOT repeat candidate-field screening,
%   held-out spatial prediction, detailed flux ledgers, rescue experiments,
%   hotspot topology, or parameter sensitivity in every geometry.
%
% GEOMETRY SET
%   1. Reference-like: exact finalized main-text structure.
%   2. Elongated: area-preserving increase in aspect ratio; crista count
%      retained at seven.
%   3. Rounded: area-preserving decrease in aspect ratio; crista count
%      retained at seven.
%   4. Sparse cristae: reference outline with four cristae.
%   5. Dense cristae: reference outline with ten narrower cristae.
%
% CONTROL LOGIC
%   For every geometry, Full defines its own final global ROS target.
%   U/H/B are independently matched to that target on the SAME geometry.
%   All kinetics, initial iron state, window radius, q=0.85, threshold rule,
%   harmonic interface rates, and exact next-subvolume SSA are unchanged.
%
% OUTPUTS
%   FigS_geometryGeneralization_FINAL.{pdf,png,tif,fig}
%   FigS_geometryGeneralization_audit.csv
%   FigS_geometryGeneralization_family_summary.csv
%   FigS_geometryGeneralization_data.mat
% -------------------------------------------------------------------------

clear; close all; clc;
clear functions;

%% --------------------------- Run configuration --------------------------
RUN_MODE = 'final';       % 'smoke' | 'pilot' | 'final'
USE_PARALLEL = true;
AUTO_START_POOL = true;
N_WORKERS = 8;
SAVE_FIG = true;
SAVE_DATA = true;

switch lower(RUN_MODE)
    case 'smoke'
        nAlternativeSeeds = 1;   % 1 exact reference + 4 alternatives
        nRepFPT = 20;
        nBootTrajectory = 100;
        nBootGeometry = 300;
    case 'pilot'
        nAlternativeSeeds = 2;   % 1 exact reference + 8 alternatives
        nRepFPT = 60;
        nBootTrajectory = 500;
        nBootGeometry = 1000;
    case 'final'
        nAlternativeSeeds = 2;   % final SI workload: 9 geometries total
        nRepFPT = 150;
        nBootTrajectory = 1000;
        nBootGeometry = 2000;
    otherwise
        error('Unknown RUN_MODE: %s',RUN_MODE);
end

seedBase = 20260805;
alphaCI = 0.05;

outFigure = 'FigS_geometryGeneralization_FINAL';
outAuditCSV = 'FigS_geometryGeneralization_audit.csv';
outFamilyCSV = 'FigS_geometryGeneralization_family_summary.csv';
outMAT = 'FigS_geometryGeneralization_data.mat';

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

familyColors = [
    0.34 0.37 0.42
    0.12 0.47 0.71
    0.20 0.60 0.38
    0.93 0.55 0.14
    0.78 0.25 0.20
];

%% ---------------- Locked finalized Fig. 3 parameters -------------------
% Exact reference grid and geometry parameters.
nx = 60;
ny = 36;
aOutRef = 1.12;  bOutRef = 0.62;
aInRef  = 0.94;  bInRef  = 0.49;

ellG = 0.085;
ellS = 0.042;

connMatrix = 1.00;
connLumen  = 0.82;
connNeck   = 0.18;

% Reaction and diffusion parameters.
D0_Fe2   = 3.0e-4;
D2_scale = 0.70;
D0_Fe3   = D2_scale*D0_Fe2;

c1_0 = 6.0e-3;
etaS = 3.40;
c2_0 = 1.8e-2;

% Observation and initial state.
Tobs = 60;
dtRecord = 0.50;
tvec = 0:dtRecord:Tobs;

FeTot_per_voxel = 6;
fracFe2_0 = 0.65;   % rounded initial state = Fe2+/Fe3+ = 4/2

% Exact first-moment matching.
dtMoment = 0.05;
nBisect = 14;

% Local definition inherited from the finalized c1-hat workflow.
windowRadius = 0.10;
qA = 0.85;
topFraction = 1-qA;
thresholdFactor = 1.00;

%% --------------------------- Parallel pool ------------------------------
usePar = configureParallel(USE_PARALLEL,AUTO_START_POOL,N_WORKERS);

%% ------------------------- Geometry definitions -------------------------
families = geometryFamilyDefinitionsFinal( ...
    aOutRef,bOutRef,aInRef,bInRef,nAlternativeSeeds);

nFamily = numel(families);
nGeometry = sum([families.nInstances]);

familyName = strings(nGeometry,1);
familyShort = strings(nGeometry,1);
familyIndex = zeros(nGeometry,1);
instanceWithinFamily = zeros(nGeometry,1);
geometrySeed = zeros(nGeometry,1);

Kvec = zeros(nGeometry,1);
aspectRatio = zeros(nGeometry,1);
nCristaVec = zeros(nGeometry,1);
lumenFraction = zeros(nGeometry,1);
cristaOverlapFraction = zeros(nGeometry,1);
meanS = zeros(nGeometry,1);
meanG = zeros(nGeometry,1);
cvG = zeros(nGeometry,1);

globalMismatchMax = zeros(nGeometry,1);
regionFraction = zeros(nGeometry,1);
thresholdDensity = zeros(nGeometry,1);

heterogeneityU = zeros(nGeometry,1);
heterogeneityF = zeros(nGeometry,1);
deltaBurdenHeterogeneity = zeros(nGeometry,1);
fractionDeltaBurdenPositive = zeros(nGeometry,1);
localBurdenRatioFU = zeros(nGeometry,1);
predictorJaccard = zeros(nGeometry,1);

RMST_U = zeros(nGeometry,1);
RMST_H = zeros(nGeometry,1);
RMST_B = zeros(nGeometry,1);
RMST_F = zeros(nGeometry,1);

DeltaRMST_UF = zeros(nGeometry,1);
DeltaRMST_UF_L = zeros(nGeometry,1);
DeltaRMST_UF_U = zeros(nGeometry,1);

EffectH = zeros(nGeometry,1);
EffectH_L = zeros(nGeometry,1);
EffectH_U = zeros(nGeometry,1);

EffectB = zeros(nGeometry,1);
EffectB_L = zeros(nGeometry,1);
EffectB_U = zeros(nGeometry,1);

EffectBgivenH = zeros(nGeometry,1);
EffectBgivenH_L = zeros(nGeometry,1);
EffectBgivenH_U = zeros(nGeometry,1);

P60_U = zeros(nGeometry,1);
P60_H = zeros(nGeometry,1);
P60_B = zeros(nGeometry,1);
P60_F = zeros(nGeometry,1);

RepresentativeGeometry = cell(nFamily,1);

fprintf('============================================================\n');
fprintf('FINAL SI geometry generalization (%s)\n',upper(RUN_MODE));
fprintf('Geometry families             = %d\n',nFamily);
fprintf('Total geometry instances      = %d\n',nGeometry);
fprintf('FPT trajectories/model        = %d\n',nRepFPT);
fprintf('Trajectory bootstrap          = %d\n',nBootTrajectory);
fprintf('Geometry bootstrap            = %d\n',nBootGeometry);
fprintf('============================================================\n');

%% ----------------------------- Main loop --------------------------------
gCounter = 0;

for f = 1:nFamily
    for r = 1:families(f).nInstances
        gCounter = gCounter+1;
        thisSeed = seedBase + 10000*f + 100*r;

        fprintf('\n------------------------------------------------------------\n');
        fprintf('Geometry %d/%d: %s, instance %d, seed %d\n', ...
            gCounter,nGeometry,families(f).label,r,thisSeed);
        fprintf('------------------------------------------------------------\n');

        [geom,geomSpec] = generateAcceptedGeometryFinal( ...
            families(f),r,thisSeed,nx,ny, ...
            ellG,ellS,connMatrix,connLumen,connNeck);

        mask = geom.OmegaR;
        x = geom.x;
        y = geom.y;
        Xg = geom.X;
        Yg = geom.Y;
        S = geom.S;
        G = geom.G;
        K = nnz(mask);
        idxMap = buildMaskIndexMap(mask);

        % The first geometry is required to be the exact finalized reference.
        if families(f).isExactReference
            assert(abs(geomSpec.aOut-aOutRef)<1e-12 && ...
                   abs(geomSpec.bOut-bOutRef)<1e-12 && ...
                   abs(geomSpec.aIn-aInRef)<1e-12 && ...
                   abs(geomSpec.bIn-bInRef)<1e-12, ...
                'Reference outline no longer matches the finalized Fig. 3.');
            assert(numel(geomSpec.cristaSpec)==7, ...
                'Reference crista number no longer matches the finalized Fig. 3.');
            assert(K==868, ...
                ['Reference geometry no longer has K=868 voxels. ' ...
                 'Check the reference geometry helper before proceeding.']);
        end

        % Full fields: identical definitions to finalized Fig. 3.
        D1Full = D0_Fe2*G;
        D2Full = D0_Fe3*G;
        D1Full(~mask) = NaN;
        D2Full(~mask) = NaN;

        c1Full = c1_0*(1+etaS*S);
        c1Full(~mask) = NaN;

        c1F = c1Full(mask);
        D1F = D1Full(mask);
        D2F = D2Full(mask);

        c1Mean = mean(c1F);
        D1Mean = mean(D1F);
        D2Mean = mean(D2F);
        c2vec = c2_0*ones(K,1);

        N0 = FeTot_per_voxel*ones(K,1);
        x1_0 = round(fracFe2_0*N0);
        x2_0 = N0-x1_0;

        % Geometry-specific, pre-specified c1-hat top-15% region.
        coordX = Xg(mask);
        coordY = Yg(mask);
        win = buildCircularWindowsSparse(coordX,coordY,windowRadius);

        c1Center = normalizeVector(c1F);
        c1Window = (win.W*c1Center)./win.size;
        c1Window = normalizeVector(c1Window);

        A = topFractionMask(c1Window,topFraction);
        nA = nnz(A);
        Amap = false(size(mask));
        Amap(mask) = A;

        % Diffusion graphs.
        dx = mean(diff(x));
        dy = mean(diff(y));

        [fromHet1,toHet1,qHet1] = buildDirectedEdgesHarmonic( ...
            mask,idxMap,D1Full,dx,dy);
        [fromHet2,toHet2,qHet2] = buildDirectedEdgesHarmonic( ...
            mask,idxMap,D2Full,dx,dy);

        D1Uniform = NaN(size(D1Full));
        D2Uniform = NaN(size(D2Full));
        D1Uniform(mask) = D1Mean;
        D2Uniform(mask) = D2Mean;

        [fromUni1,toUni1,qUni1] = buildDirectedEdgesHarmonic( ...
            mask,idxMap,D1Uniform,dx,dy);
        [fromUni2,toUni2,qUni2] = buildDirectedEdgesHarmonic( ...
            mask,idxMap,D2Uniform,dx,dy);

        % U/H/B/F on the same geometry.
        M = repmat(struct(),1,4);

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

        % Full defines the within-geometry target.
        M(4).alpha = 1;
        M(4).moment = momentEulerAll( ...
            M(4).c1Base,c2vec, ...
            M(4).from1,M(4).to1,M(4).q1, ...
            M(4).from2,M(4).to2,M(4).q2, ...
            x1_0,x2_0,tvec,dtMoment);

        targetGlobal = M(4).moment.globalM3(end);

        for m = 1:3
            [M(m).alpha,M(m).moment] = matchModelToTargetFinal( ...
                M(m).c1Base,c2vec, ...
                M(m).from1,M(m).to1,M(m).q1, ...
                M(m).from2,M(m).to2,M(m).q2, ...
                x1_0,x2_0,tvec,dtMoment,targetGlobal,nBisect);
        end

        for m = 1:4
            M(m).c1 = M(m).alpha*M(m).c1Base;
            M(m).out = buildOutgoingLists(K, ...
                M(m).from1,M(m).to1,M(m).q1, ...
                M(m).from2,M(m).to2,M(m).q2);
            M(m).momentRatio = M(m).moment.globalM3(end)/targetGlobal;
        end

        % Predictor-independent burden heterogeneity across ALL windows.
        windowBurdenU = (win.W*M(1).moment.finalM3)./win.size;
        windowBurdenF = (win.W*M(4).moment.finalM3)./win.size;

        HBU = normalizedP90P10(windowBurdenU);
        HBF = normalizedP90P10(windowBurdenF);
        deltaHB = HBF-HBU;

        deltaWindowBurden = windowBurdenF-windowBurdenU;
        fracPositive = mean(deltaWindowBurden>0);

        actualTop = topFractionMask(deltaWindowBurden,topFraction);
        jac = nnz(A&actualTop)/max(nnz(A|actualTop),1);

        localU = mean(M(1).moment.finalM3(A));
        localF = mean(M(4).moment.finalM3(A));
        localRatioFU = localF/max(localU,eps);

        % Control-anchored, geometry-specific regional threshold.
        mStar = M(1).moment.globalM3(end)/K;
        M_A = max(1,round(thresholdFactor*mStar*nA));

        % Exact regional FPT ensembles. The same random-number stream index
        % is used across U/H/B/F to reduce Monte Carlo contrast variance.
        FPT = cell(4,1);
        commonFPTSeed = seedBase + 1000000 + 10000*gCounter;

        for m = 1:4
            FPT{m} = runRegionFPTensembleNSM( ...
                M(m).c1,c2vec,M(m).out,A, ...
                x1_0,x2_0,M_A,Tobs,nRepFPT,commonFPTSeed,usePar);
        end

        rmst = zeros(4,1);
        p60 = zeros(4,1);
        for m = 1:4
            rmst(m) = mean(min(FPT{m},Tobs));
            p60(m) = mean(FPT{m}<=Tobs);
        end

        [dUF,ciUF] = bootstrapPairedDifference( ...
            min(FPT{1},Tobs),min(FPT{4},Tobs), ...
            nBootTrajectory,alphaCI,thisSeed+5001);

        [eH,ciH] = bootstrapPairedDifference( ...
            min(FPT{1},Tobs),min(FPT{2},Tobs), ...
            nBootTrajectory,alphaCI,thisSeed+5002);

        [eB,ciB] = bootstrapPairedDifference( ...
            min(FPT{1},Tobs),min(FPT{3},Tobs), ...
            nBootTrajectory,alphaCI,thisSeed+5003);

        [eBH,ciBH] = bootstrapPairedDifference( ...
            min(FPT{2},Tobs),min(FPT{4},Tobs), ...
            nBootTrajectory,alphaCI,thisSeed+5004);

        % Store audit and endpoints.
        familyName(gCounter)=string(families(f).label);
        familyShort(gCounter)=string(families(f).shortLabel);
        familyIndex(gCounter)=f;
        instanceWithinFamily(gCounter)=r;
        geometrySeed(gCounter)=thisSeed;

        Kvec(gCounter)=K;
        aspectRatio(gCounter)=geomSpec.aIn/geomSpec.bIn;
        nCristaVec(gCounter)=numel(geomSpec.cristaSpec);
        lumenFraction(gCounter)=nnz(geom.cristaLumenMask)/K;
        cristaOverlapFraction(gCounter)=geom.cristaOverlapFraction;
        meanS(gCounter)=mean(S(mask));
        meanG(gCounter)=mean(G(mask));
        cvG(gCounter)=std(G(mask))/max(mean(G(mask)),eps);

        globalMismatchMax(gCounter)=max(abs([M.momentRatio]-1));
        regionFraction(gCounter)=nA/K;
        thresholdDensity(gCounter)=M_A/nA;

        heterogeneityU(gCounter)=HBU;
        heterogeneityF(gCounter)=HBF;
        deltaBurdenHeterogeneity(gCounter)=deltaHB;
        fractionDeltaBurdenPositive(gCounter)=fracPositive;
        localBurdenRatioFU(gCounter)=localRatioFU;
        predictorJaccard(gCounter)=jac;

        RMST_U(gCounter)=rmst(1);
        RMST_H(gCounter)=rmst(2);
        RMST_B(gCounter)=rmst(3);
        RMST_F(gCounter)=rmst(4);

        DeltaRMST_UF(gCounter)=dUF;
        DeltaRMST_UF_L(gCounter)=ciUF(1);
        DeltaRMST_UF_U(gCounter)=ciUF(2);

        EffectH(gCounter)=eH;
        EffectH_L(gCounter)=ciH(1);
        EffectH_U(gCounter)=ciH(2);

        EffectB(gCounter)=eB;
        EffectB_L(gCounter)=ciB(1);
        EffectB_U(gCounter)=ciB(2);

        EffectBgivenH(gCounter)=eBH;
        EffectBgivenH_L(gCounter)=ciBH(1);
        EffectBgivenH_U(gCounter)=ciBH(2);

        P60_U(gCounter)=p60(1);
        P60_H(gCounter)=p60(2);
        P60_B(gCounter)=p60(3);
        P60_F(gCounter)=p60(4);

        if r==1
            rep=struct();
            rep.x=x; rep.y=y; rep.mask=mask;
            rep.c1Map=normalizeInside(c1Full,mask);
            rep.Amap=Amap;
            rep.family=families(f).label;
            rep.shortLabel=families(f).shortLabel;
            rep.aspect=aspectRatio(gCounter);
            rep.nCrista=nCristaVec(gCounter);
            RepresentativeGeometry{f}=rep;
        end

        fprintf('K=%d; AR=%.2f; nCrista=%d; lumen frac=%.3f\n', ...
            K,aspectRatio(gCounter),nCristaVec(gCounter), ...
            lumenFraction(gCounter));
        fprintf('max global mismatch=%.4f%%\n', ...
            100*globalMismatchMax(gCounter));
        fprintf('Delta H_B=%.4f; F/U local ratio=%.3f; Jaccard=%.3f\n', ...
            deltaHB,localRatioFU,jac);
        fprintf('Delta RMST U-F=%.3f [%.3f, %.3f] min\n', ...
            dUF,ciUF(1),ciUF(2));
        fprintf('E_H=%.3f; E_B=%.3f; E_B|H=%.3f min\n',eH,eB,eBH);
    end
end

%% -------------------- Across-geometry family-balanced CI ----------------
[medianDeltaHB,ciDeltaHB] = bootstrapFamilyBalancedMedianCI( ...
    deltaBurdenHeterogeneity,familyIndex,nFamily, ...
    nBootGeometry,alphaCI,seedBase+8001);

[medianDeltaUF,ciDeltaUF] = bootstrapFamilyBalancedMedianCI( ...
    DeltaRMST_UF,familyIndex,nFamily, ...
    nBootGeometry,alphaCI,seedBase+8002);

[medianEH,ciEH] = bootstrapFamilyBalancedMedianCI( ...
    EffectH,familyIndex,nFamily,nBootGeometry,alphaCI,seedBase+8003);
[medianEB,ciEB] = bootstrapFamilyBalancedMedianCI( ...
    EffectB,familyIndex,nFamily,nBootGeometry,alphaCI,seedBase+8004);
[medianEBH,ciEBH] = bootstrapFamilyBalancedMedianCI( ...
    EffectBgivenH,familyIndex,nFamily,nBootGeometry,alphaCI,seedBase+8005);

fprintf('\n================ Geometry-generalization summary ================\n');
fprintf('Delta H_B > 0                         = %d/%d geometries\n', ...
    nnz(deltaBurdenHeterogeneity>0),nGeometry);
fprintf('Median Delta H_B                      = %.4f [%.4f, %.4f]\n', ...
    medianDeltaHB,ciDeltaHB(1),ciDeltaHB(2));
fprintf('Delta RMST U-F > 0                    = %d/%d geometries\n', ...
    nnz(DeltaRMST_UF>0),nGeometry);
fprintf('Median Delta RMST U-F                 = %.3f [%.3f, %.3f]\n', ...
    medianDeltaUF,ciDeltaUF(1),ciDeltaUF(2));
fprintf('Median E_H                            = %.3f [%.3f, %.3f]\n', ...
    medianEH,ciEH(1),ciEH(2));
fprintf('Median E_B                            = %.3f [%.3f, %.3f]\n', ...
    medianEB,ciEB(1),ciEB(2));
fprintf('Median E_B|H                          = %.3f [%.3f, %.3f]\n', ...
    medianEBH,ciEBH(1),ciEBH(2));
fprintf('Maximum global mismatch               = %.4f%%\n', ...
    100*max(globalMismatchMax));
fprintf('=================================================================\n');

%% ------------------------------ SI figure -------------------------------
fig = figure('Units','centimeters','Position',[1.0 1.0 18.8 12.0], ...
    'Color','w');
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

% Reserve tile A and replace it with a nested five-geometry panel.
axPlaceholder = nexttile(tl,1);
drawnow;
panelPos = axPlaceholder.Position;
delete(axPlaceholder);

panelA = uipanel(fig,'Units','normalized','Position',panelPos, ...
    'BorderType','none','BackgroundColor','w');
tlA = tiledlayout(panelA,1,nFamily,'TileSpacing','compact','Padding','compact');
title(tlA,'Alternative mitochondrial geometries', ...
    'FontSize',8.3,'FontWeight','normal');

for f=1:nFamily
    ax=nexttile(tlA,f);
    plotCompactGeometry(ax,RepresentativeGeometry{f},cmapSequentialBlue(256));
    title(ax,RepresentativeGeometry{f}.shortLabel, ...
        'FontSize',6.4,'FontWeight','normal');
    if f==1
        text(ax,0.01,1.07,'A','Units','normalized', ...
            'FontSize',11.5,'FontWeight','bold', ...
            'VerticalAlignment','top');
    end
end

% B. Predictor-independent burden heterogeneity.
ax = nexttile(tl,2); hold(ax,'on');
plot(ax,[0.5 nFamily+0.5],[0 0],':','Color',[0.55 0.55 0.55], ...
    'LineWidth',0.8);

for g=1:nGeometry
    xg=familyIndex(g)+deterministicJitter( ...
        instanceWithinFamily(g),families(familyIndex(g)).nInstances,0.11);
    scatter(ax,xg,deltaBurdenHeterogeneity(g),27, ...
        familyColors(familyIndex(g),:),'filled', ...
        'MarkerEdgeColor','w','LineWidth',0.45);
end

for f=1:nFamily
    idx=familyIndex==f;
    medf=median(deltaBurdenHeterogeneity(idx));
    plot(ax,[f-0.16 f+0.16],[medf medf],'-', ...
        'Color',familyColors(f,:),'LineWidth',1.5);
end

set(ax,'XTick',1:nFamily,'XTickLabel',{families.shortLabel});
xtickangle(ax,24);
ylabel(ax,'\Delta H_B = H_B^F - H_B^U');
title(ax,'Local burden heterogeneity generalizes', ...
    'FontSize',8.3,'FontWeight','normal');
text(ax,0.02,0.97,'B','Units','normalized','FontSize',11.5, ...
    'FontWeight','bold','VerticalAlignment','top');
text(ax,0.98,0.07,sprintf('max global mismatch = %.2f%%', ...
    100*max(globalMismatchMax)), ...
    'Units','normalized','HorizontalAlignment','right','FontSize',5.6);
xlim(ax,[0.5 nFamily+0.5]);
formatAxes(ax);

% C. Full-versus-Uniform regional FPT effect.
ax = nexttile(tl,3); hold(ax,'on');
plot(ax,[0.5 nFamily+0.5],[0 0],':','Color',[0.55 0.55 0.55], ...
    'LineWidth',0.8);

for g=1:nGeometry
    xg=familyIndex(g)+deterministicJitter( ...
        instanceWithinFamily(g),families(familyIndex(g)).nInstances,0.11);
    errorbar(ax,xg,DeltaRMST_UF(g), ...
        max(DeltaRMST_UF(g)-DeltaRMST_UF_L(g),0), ...
        max(DeltaRMST_UF_U(g)-DeltaRMST_UF(g),0), ...
        'o','Color',familyColors(familyIndex(g),:), ...
        'MarkerFaceColor','w','MarkerSize',3.8, ...
        'LineWidth',0.72,'CapSize',3);
end

set(ax,'XTick',1:nFamily,'XTickLabel',{families.shortLabel});
xtickangle(ax,24);
ylabel(ax,'\DeltaRMST_{60}, U-F (min)');
title(ax,'Local FPT advancement across geometries', ...
    'FontSize',8.3,'FontWeight','normal');
text(ax,0.02,0.97,'C','Units','normalized','FontSize',11.5, ...
    'FontWeight','bold','VerticalAlignment','top');
text(ax,0.98,0.07,sprintf('%d/%d positive; median %.2f min', ...
    nnz(DeltaRMST_UF>0),nGeometry,medianDeltaUF), ...
    'Units','normalized','HorizontalAlignment','right','FontSize',5.6);
xlim(ax,[0.5 nFamily+0.5]);
formatAxes(ax);

% D. U/H/B/F causal ordering.
ax = nexttile(tl,4); hold(ax,'on');
plot(ax,[0.5 3.5],[0 0],':','Color',[0.55 0.55 0.55], ...
    'LineWidth',0.8);

effectData=[EffectH EffectB EffectBgivenH];
effectColors=[colH;colB;colF];
effectMedian=[medianEH medianEB medianEBH];
effectCI=[ciEH;ciEB;ciEBH];

for e=1:3
    lightCol=lightenColor(effectColors(e,:),0.55);
    for g=1:nGeometry
        xg=e+deterministicJitter(g,nGeometry,0.18);
        scatter(ax,xg,effectData(g,e),17,lightCol,'filled', ...
            'MarkerEdgeColor','none');
    end

    errorbar(ax,e,effectMedian(e), ...
        max(effectMedian(e)-effectCI(e,1),0), ...
        max(effectCI(e,2)-effectMedian(e),0), ...
        'o','Color',effectColors(e,:), ...
        'MarkerFaceColor','w','MarkerSize',5.1, ...
        'LineWidth',1.15,'CapSize',5);
end

set(ax,'XTick',1:3, ...
    'XTickLabel',{'hotspot alone','bottleneck alone','bottleneck | hotspot'});
xtickangle(ax,18);
ylabel(ax,'RMST advancement (min)');
title(ax,'Reaction-hotspot dominance is geometry-robust', ...
    'FontSize',8.3,'FontWeight','normal');
text(ax,0.02,0.97,'D','Units','normalized','FontSize',11.5, ...
    'FontWeight','bold','VerticalAlignment','top');
xlim(ax,[0.5 3.5]);
formatAxes(ax);

if SAVE_FIG
    exportgraphics(fig,[outFigure '.pdf'],'ContentType','vector');
    exportgraphics(fig,[outFigure '.png'],'Resolution',600);
    exportgraphics(fig,[outFigure '.tif'],'Resolution',600);
    savefig(fig,[outFigure '.fig']);
end

%% ------------------------------ Audit table -----------------------------
auditTable = table( ...
    familyName,familyShort,familyIndex,instanceWithinFamily,geometrySeed, ...
    Kvec,aspectRatio,nCristaVec,lumenFraction,cristaOverlapFraction, ...
    meanS,meanG,cvG,globalMismatchMax,regionFraction,thresholdDensity, ...
    heterogeneityU,heterogeneityF,deltaBurdenHeterogeneity, ...
    fractionDeltaBurdenPositive,localBurdenRatioFU,predictorJaccard, ...
    RMST_U,RMST_H,RMST_B,RMST_F, ...
    DeltaRMST_UF,DeltaRMST_UF_L,DeltaRMST_UF_U, ...
    EffectH,EffectH_L,EffectH_U, ...
    EffectB,EffectB_L,EffectB_U, ...
    EffectBgivenH,EffectBgivenH_L,EffectBgivenH_U, ...
    P60_U,P60_H,P60_B,P60_F);

familyRows = strings(nFamily,1);
familyN = zeros(nFamily,1);
familyMedianDeltaHB = zeros(nFamily,1);
familyMedianDeltaUF = zeros(nFamily,1);
familyMedianEH = zeros(nFamily,1);
familyMedianEB = zeros(nFamily,1);
familyMedianEBH = zeros(nFamily,1);
familyMaxMismatch = zeros(nFamily,1);

for f=1:nFamily
    idx=familyIndex==f;
    familyRows(f)=string(families(f).label);
    familyN(f)=nnz(idx);
    familyMedianDeltaHB(f)=median(deltaBurdenHeterogeneity(idx));
    familyMedianDeltaUF(f)=median(DeltaRMST_UF(idx));
    familyMedianEH(f)=median(EffectH(idx));
    familyMedianEB(f)=median(EffectB(idx));
    familyMedianEBH(f)=median(EffectBgivenH(idx));
    familyMaxMismatch(f)=max(globalMismatchMax(idx));
end

familySummaryTable = table( ...
    familyRows,familyN,familyMedianDeltaHB,familyMedianDeltaUF, ...
    familyMedianEH,familyMedianEB,familyMedianEBH,familyMaxMismatch, ...
    'VariableNames',{'GeometryFamily','Ninstances', ...
    'MedianDeltaBurdenHeterogeneity','MedianDeltaRMST_UF', ...
    'MedianEffectH','MedianEffectB','MedianEffectBgivenH', ...
    'MaximumGlobalMismatch'});

if SAVE_DATA
    writetable(auditTable,outAuditCSV);
    writetable(familySummaryTable,outFamilyCSV);

    metadata=struct();
    metadata.RUN_MODE=RUN_MODE;
    metadata.nAlternativeSeeds=nAlternativeSeeds;
    metadata.nRepFPT=nRepFPT;
    metadata.nBootTrajectory=nBootTrajectory;
    metadata.nBootGeometry=nBootGeometry;
    metadata.referenceK=868;
    metadata.windowRadius=windowRadius;
    metadata.qA=qA;
    metadata.families=families;
    metadata.kinetics=struct( ...
        'D0_Fe2',D0_Fe2,'D0_Fe3',D0_Fe3, ...
        'c1_0',c1_0,'etaS',etaS,'c2_0',c2_0, ...
        'Tobs',Tobs,'FeTot_per_voxel',FeTot_per_voxel, ...
        'fracFe2_0',fracFe2_0);

    save(outMAT,'metadata','auditTable','familySummaryTable', ...
        'RepresentativeGeometry', ...
        'medianDeltaHB','ciDeltaHB','medianDeltaUF','ciDeltaUF', ...
        'medianEH','ciEH','medianEB','ciEB','medianEBH','ciEBH', ...
        '-v7.3');
end

fprintf('\nCompleted. Outputs:\n');
fprintf('  %s.{pdf,png,tif,fig}\n',outFigure);
fprintf('  %s\n',outAuditCSV);
fprintf('  %s\n',outFamilyCSV);
fprintf('  %s\n',outMAT);

%% ======================== Analysis-specific functions ===================

function families = geometryFamilyDefinitionsFinal( ...
    aOutRef,bOutRef,aInRef,bInRef,nAlternativeSeeds)

    % MATLAB does not permit assigning a populated structure into one
    % element of a structure array that was created from struct() with no
    % fields.  Initialize a full-field template first so every element has
    % exactly the same field set.
    familyTemplate=struct( ...
        'label','', ...
        'shortLabel','', ...
        'aOut',NaN, ...
        'bOut',NaN, ...
        'aIn',NaN, ...
        'bIn',NaN, ...
        'nCrista',NaN, ...
        'lengthScale',NaN, ...
        'widthScale',NaN, ...
        'ampScale',NaN, ...
        'nInstances',NaN, ...
        'isExactReference',false);

    families=repmat(familyTemplate,1,5);

    families(1)=struct( ...
        'label','Reference-like','shortLabel','reference', ...
        'aOut',aOutRef,'bOut',bOutRef,'aIn',aInRef,'bIn',bInRef, ...
        'nCrista',7,'lengthScale',1.00,'widthScale',1.00, ...
        'ampScale',1.00,'nInstances',1,'isExactReference',true);

    % Area-preserving elongation of outer and inner ellipses.
    families(2)=struct( ...
        'label','Elongated','shortLabel','elongated', ...
        'aOut',1.20,'bOut',(aOutRef*bOutRef)/1.20, ...
        'aIn',1.06,'bIn',(aInRef*bInRef)/1.06, ...
        'nCrista',7,'lengthScale',0.98,'widthScale',0.92, ...
        'ampScale',1.00,'nInstances',nAlternativeSeeds, ...
        'isExactReference',false);

    % Area-preserving rounded outline; crista number retained.
    families(3)=struct( ...
        'label','Rounded','shortLabel','rounded', ...
        'aOut',1.00,'bOut',(aOutRef*bOutRef)/1.00, ...
        'aIn',0.84,'bIn',(aInRef*bInRef)/0.84, ...
        'nCrista',7,'lengthScale',1.02,'widthScale',1.00, ...
        'ampScale',0.95,'nInstances',nAlternativeSeeds, ...
        'isExactReference',false);

    % Internal-structure perturbations on the exact reference outline.
    families(4)=struct( ...
        'label','Sparse cristae','shortLabel','sparse', ...
        'aOut',aOutRef,'bOut',bOutRef,'aIn',aInRef,'bIn',bInRef, ...
        'nCrista',4,'lengthScale',1.00,'widthScale',1.08, ...
        'ampScale',1.00,'nInstances',nAlternativeSeeds, ...
        'isExactReference',false);

    families(5)=struct( ...
        'label','Dense cristae','shortLabel','dense', ...
        'aOut',aOutRef,'bOut',bOutRef,'aIn',aInRef,'bIn',bInRef, ...
        'nCrista',10,'lengthScale',0.82,'widthScale',0.74, ...
        'ampScale',0.92,'nInstances',nAlternativeSeeds, ...
        'isExactReference',false);
end

function [geom,spec] = generateAcceptedGeometryFinal( ...
    family,replicateIndex,seed,nx,ny, ...
    ellG,ellS,connMatrix,connLumen,connNeck)

    if family.isExactReference
        spec=struct();
        spec.aOut=family.aOut; spec.bOut=family.bOut;
        spec.aIn=family.aIn; spec.bIn=family.bIn;
        spec.cristaSpec=referenceCristaSpecFinal();

        geom=buildGeometryFromSpecFinal( ...
            nx,ny,spec.aOut,spec.bOut,spec.aIn,spec.bIn, ...
            ellG,ellS,connMatrix,connLumen,connNeck,spec.cristaSpec);
        return;
    end

    maxAttempts=60;

    for attempt=1:maxAttempts
        thisSeed=seed+7919*(attempt-1);
        spec=makeGeometryInstanceFinal(family,replicateIndex,thisSeed);

        geom=buildGeometryFromSpecFinal( ...
            nx,ny,spec.aOut,spec.bOut,spec.aIn,spec.bIn, ...
            ellG,ellS,connMatrix,connLumen,connNeck,spec.cristaSpec);

        K=nnz(geom.OmegaR);
        lumenFrac=nnz(geom.cristaLumenMask)/max(K,1);

        valid=K>=790 && K<=940 ...
            && lumenFrac>=0.012 && lumenFrac<=0.34 ...
            && geom.cristaOverlapFraction<=0.035 ...
            && range(geom.S(geom.OmegaR))>0.50 ...
            && range(geom.G(geom.OmegaR))>0.15;

        if valid
            return;
        end
    end

    error('Could not generate a valid %s geometry after %d attempts.', ...
        family.label,maxAttempts);
end

function spec = makeGeometryInstanceFinal(family,replicateIndex,seed)
    rng(seed,'twister');

    spec.aOut=family.aOut;
    spec.bOut=family.bOut;
    spec.aIn=family.aIn;
    spec.bIn=family.bIn;

    % Independent geometry seeds perturb crista placement and curvature,
    % not the defining family-level outline or crista count.
    spec.cristaSpec=proceduralCristaSpecFinal( ...
        family.nCrista,spec.aIn, ...
        family.lengthScale,family.widthScale,family.ampScale,seed);
end

function sp = referenceCristaSpecFinal()
    % Exact finalized main-text crista specification.
    sp=repmat(struct('x0',0,'fromTop',true,'L',0,'w',0, ...
        'amp',0,'phase',0,'bias',0),1,7);

    sp(1).x0=-0.68; sp(1).fromTop=true;  sp(1).L=0.34; sp(1).w=0.043; sp(1).amp=0.026; sp(1).phase=0.30; sp(1).bias=-0.030;
    sp(2).x0=-0.48; sp(2).fromTop=false; sp(2).L=0.42; sp(2).w=0.050; sp(2).amp=0.032; sp(2).phase=1.40; sp(2).bias= 0.020;
    sp(3).x0=-0.04; sp(3).fromTop=true;  sp(3).L=0.68; sp(3).w=0.060; sp(3).amp=0.018; sp(3).phase=2.20; sp(3).bias= 0.000;
    sp(4).x0= 0.39; sp(4).fromTop=false; sp(4).L=0.42; sp(4).w=0.050; sp(4).amp=0.026; sp(4).phase=0.70; sp(4).bias=-0.020;
    sp(5).x0= 0.55; sp(5).fromTop=true;  sp(5).L=0.54; sp(5).w=0.056; sp(5).amp=0.020; sp(5).phase=2.10; sp(5).bias= 0.010;
    sp(6).x0= 0.67; sp(6).fromTop=false; sp(6).L=0.37; sp(6).w=0.046; sp(6).amp=0.024; sp(6).phase=3.00; sp(6).bias= 0.030;
    sp(7).x0= 0.79; sp(7).fromTop=true;  sp(7).L=0.30; sp(7).w=0.040; sp(7).amp=0.020; sp(7).phase=1.10; sp(7).bias=-0.012;
end

function sp = proceduralCristaSpecFinal( ...
    nCrista,aIn,lengthScale,widthScale,ampScale,seed)

    rng(seed+1729,'twister');

    sp=repmat(struct('x0',0,'fromTop',true,'L',0,'w',0, ...
        'amp',0,'phase',0,'bias',0),1,nCrista);

    xBase=linspace(-0.79*aIn,0.79*aIn,nCrista);

    if nCrista>1
        spacing=median(diff(xBase));
    else
        spacing=0.20*aIn;
    end

    xJitter=0.11*spacing*randn(1,nCrista);
    xList=sort(xBase+xJitter);
    xList=min(max(xList,-0.84*aIn),0.84*aIn);

    orientationOffset=randi([0 1]);

    for j=1:nCrista
        radial=abs(xList(j))/max(aIn,eps);
        centrality=max(0,1-radial);

        sp(j).x0=xList(j);
        sp(j).fromTop=mod(j+orientationOffset,2)==1;
        sp(j).L=lengthScale*(0.30+0.37*centrality)*exp(0.05*randn);
        sp(j).w=widthScale*(0.041+0.015*centrality)*exp(0.04*randn);
        sp(j).amp=ampScale*(0.016+0.012*rand)*exp(0.09*randn);
        sp(j).phase=2*pi*rand;
        sp(j).bias=0.018*randn;
    end
end

function geom = buildGeometryFromSpecFinal( ...
    nx,ny,aOut,bOut,aIn,bIn, ...
    ellG,ellS,connMatrix,connLumen,connNeck,cristaSpec)

    x=linspace(-1.23,1.23,nx);
    y=linspace(-0.70,0.70,ny);
    [X,Y]=meshgrid(x,y);

    rhoOut=sqrt((X./aOut).^2+(Y./bOut).^2);
    rhoIn=sqrt((X./aIn).^2+(Y./bIn).^2);

    maskOut=rhoOut<=1;
    OmegaR=rhoIn<=1;

    nCrista=numel(cristaSpec);
    nPtsCrista=200;

    crista=cell(1,nCrista);
    lumenCount=zeros(size(X));
    cristaLumenMask=false(size(X));

    for j=1:nCrista
        crTemp=generateCristaAnchored(cristaSpec(j),aIn,bIn,nPtsCrista);
        thisLumen=inpolygon(X,Y,crTemp.polyX,crTemp.polyY)&OmegaR;
        crTemp.lumenMask=thisLumen;
        crista{j}=crTemp;
        lumenCount=lumenCount+double(thisLumen);
        cristaLumenMask=cristaLumenMask|thisLumen;
    end

    matrixMask=OmegaR&~cristaLumenMask;
    cristaOverlapFraction=nnz(lumenCount>1)/max(nnz(OmegaR),1);

    thetaAll=linspace(-pi,pi,1200);
    xOuter=aOut*cos(thetaAll);
    yOuter=bOut*sin(thetaAll);

    gapMask=true(size(thetaAll));
    for j=1:nCrista
        cr=crista{j};
        dth=angdiff_local(thetaAll,cr.theta0);
        gapMask=gapMask&(abs(dth)>cr.gapHalfTheta);
    end

    xInnerPlot=aIn*cos(thetaAll);
    yInnerPlot=bIn*sin(thetaAll);
    xInnerPlot(~gapMask)=NaN;
    yInnerPlot(~gapMask)=NaN;

    memX=xInnerPlot(gapMask);
    memY=yInnerPlot(gapMask);

    for j=1:nCrista
        cr=crista{j};
        memX=[memX(:);cr.xMem(:)]; %#ok<AGROW>
        memY=[memY(:);cr.yMem(:)]; %#ok<AGROW>
    end

    sampleStep=3;
    memX=memX(1:sampleStep:end);
    memY=memY(1:sampleStep:end);

    dR=inf(size(X));
    for p=1:numel(memX)
        dR=min(dR,hypot(X-memX(p),Y-memY(p)));
    end
    dR(~OmegaR)=NaN;

    chiBase=connMatrix*ones(size(X));
    chiBase(cristaLumenMask)=connLumen;
    chiBase(~OmegaR)=NaN;

    Kneck=zeros(size(X));
    for j=1:nCrista
        cr=crista{j};
        dNeck=inf(size(X));
        for p=1:numel(cr.neckX)
            dNeck=min(dNeck,hypot(X-cr.neckX(p),Y-cr.neckY(p)));
        end
        sigmaNeck=0.75*cr.neckRadius;
        thisK=exp(-(dNeck.^2)/(2*sigmaNeck^2));
        Kneck=max(Kneck,thisK);
    end
    Kneck(~OmegaR)=NaN;

    chiG=chiBase-(chiBase-connNeck).*Kneck;
    chiG(~OmegaR)=NaN;
    chiG=min(max(chiG,connNeck),connMatrix);
    chiG(~OmegaR)=NaN;

    G=(1-exp(-dR./ellG)).*chiG;
    G(~OmegaR)=NaN;
    G=min(max(G,0),1);

    S0=zeros(size(X));
    for p=1:numel(memX)
        S0=S0+exp(-((X-memX(p)).^2+(Y-memY(p)).^2)/(2*ellS^2));
    end
    S0(~OmegaR)=NaN;
    S=normalizeInside(S0,OmegaR);

    geom.x=x; geom.y=y; geom.X=X; geom.Y=Y;
    geom.maskOut=maskOut; geom.OmegaR=OmegaR;
    geom.matrixMask=matrixMask;
    geom.cristaLumenMask=cristaLumenMask;
    geom.crista=crista;
    geom.xOuter=xOuter; geom.yOuter=yOuter;
    geom.xInnerPlot=xInnerPlot; geom.yInnerPlot=yInnerPlot;
    geom.G=G; geom.S=S;
    geom.cristaOverlapFraction=cristaOverlapFraction;
end

function [alpha,momentOut] = matchModelToTargetFinal( ...
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
            error('Could not bracket the geometry-specific Full target.');
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

function FPT = runRegionFPTensembleNSM( ...
    c1,c2,outgoing,A,x1_0,x2_0,M_A,Tobs,nRep,seedBase,usePar)

    FPT=inf(nRep,1);

    if usePar
        parfor r=1:nRep
            rng(seedBase+r,'twister');
            FPT(r)=simulateRegionFPTOneNSM( ...
                c1,c2,outgoing,A,x1_0,x2_0,M_A,Tobs);
        end
    else
        for r=1:nRep
            rng(seedBase+r,'twister');
            FPT(r)=simulateRegionFPTOneNSM( ...
                c1,c2,outgoing,A,x1_0,x2_0,M_A,Tobs);
        end
    end
end

function fpt = simulateRegionFPTOneNSM( ...
    c1,c2,out,A,x1_0,x2_0,M_A,Tobs)

    x1=x1_0(:);
    x2=x2_0(:);

    localROS=0;
    fpt=inf;

    voxelRate=x1.*(c1+out.q1sum)+x2.*(c2+out.q2sum);
    tree=fenwickBuild(voxelRate);

    t=0;

    while t<Tobs
        a0=fenwickTotal(tree);
        if ~(a0>0)||~isfinite(a0), break; end

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

            if A(k)
                localROS=localROS+1;
                if localROS>=M_A
                    fpt=tNext;
                    return;
                end
            end

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

function [delta,ci] = bootstrapPairedDifference(x,y,nBoot,alpha,seed)
    x=x(:); y=y(:);
    assert(numel(x)==numel(y),'Paired samples must have equal lengths.');

    d=x-y;
    n=numel(d);
    delta=mean(d);

    boot=zeros(nBoot,1);
    rng(seed,'twister');
    for b=1:nBoot
        idx=randi(n,n,1);
        boot(b)=mean(d(idx));
    end

    ci=[quantileLocal(boot,alpha/2), ...
        quantileLocal(boot,1-alpha/2)];
end

function [med,ci] = bootstrapFamilyBalancedMedianCI( ...
    x,familyIndex,nFamily,nBoot,alpha,seed)

    x=x(:);
    familyIndex=familyIndex(:);
    med=median(x);

    boot=zeros(nBoot,1);
    rng(seed,'twister');

    for b=1:nBoot
        sampledFamilies=randi(nFamily,nFamily,1);
        sampleVals=zeros(nFamily,1);

        for j=1:nFamily
            idx=find(familyIndex==sampledFamilies(j));
            sampleVals(j)=x(idx(randi(numel(idx))));
        end

        boot(b)=median(sampleVals);
    end

    ci=[quantileLocal(boot,alpha/2), ...
        quantileLocal(boot,1-alpha/2)];
end

function H = normalizedP90P10(x)
    x=x(:);
    spread=quantileLocal(x,0.90)-quantileLocal(x,0.10);
    H=spread/max(mean(x),eps);
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
    lo=min(v);
    hi=max(v);
    v=(v-lo)/(hi-lo+eps);
end

function maskTop = topFractionMask(score,topFraction)
    score=score(:);
    K=numel(score);
    nTop=max(1,ceil(topFraction*K));

    [~,ord]=sortrows([-score,(1:K).'],[1 2]);

    maskTop=false(K,1);
    maskTop(ord(1:nTop))=true;
end

function j = deterministicJitter(index,nTotal,amplitude)
    if nTotal<=1
        j=0;
        return;
    end
    z=(index-1)/(nTotal-1);
    j=amplitude*(2*z-1);
end

function plotCompactGeometry(ax,rep,cmap)
    hold(ax,'on');

    h=imagesc(ax,rep.x,rep.y,rep.c1Map);
    set(h,'AlphaData',double(rep.mask));
    set(ax,'YDir','normal');
    colormap(ax,cmap);
    caxis(ax,[0 1]);

    contour(ax,rep.x,rep.y,double(rep.mask),[0.5 0.5], ...
        'Color',[0.16 0.16 0.16],'LineWidth',0.65);
    contour(ax,rep.x,rep.y,double(rep.Amap),[0.5 0.5], ...
        'Color','w','LineStyle','--','LineWidth',0.80);

    axis(ax,'image');
    xlim(ax,[-1.23 1.23]);
    ylim(ax,[-0.70 0.70]);
    set(ax,'XTick',[],'YTick',[],'Box','off');
end

function c = lightenColor(c0,amount)
    amount=min(max(amount,0),1);
    c=c0+(1-c0)*amount;
end

function cmap = cmapSequentialBlue(n)
    anchors=[
        0.96 0.98 0.99
        0.73 0.87 0.91
        0.29 0.63 0.73
        0.04 0.31 0.47
        0.02 0.16 0.28
    ];
    cmap=interpCmap(anchors,n);
end

%% ================= Reused validated Fig. 3 helpers =====================

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

function cmap = interpCmap(anchors,n)
    x0=linspace(0,1,size(anchors,1));
    xq=linspace(0,1,n);
    cmap=zeros(n,3);
    for c=1:3
        cmap(:,c)=interp1(x0,anchors(:,c),xq,'pchip');
    end
    cmap=max(min(cmap,1),0);
end

function Y = normalizeInside(Y,mask)
    vals=Y(mask & ~isnan(Y));
    ymin=min(vals); ymax=max(vals);
    Y=(Y-ymin)/(ymax-ymin+eps);
    Y(~mask)=NaN;
end
