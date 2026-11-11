%% SDRDME_RedoxBufferScan.m
% -------------------------------------------------------------------------
% Fe3+->Fe2+ recycling-capacity x diffusion-heterogeneity condition map.
%
% CENTRAL QUESTION
%   Under one fixed global ROS budget, when does diffusion heterogeneity:
%     (i) act as a retention amplifier of hotspot-driven Fenton risk,
%    (ii) remain neutral, or
%   (iii) become a Fe2+ supply-limiting constraint?
%
% CONDITION AXES
%   beta_H  : reaction-hotspot demand
%             c1_base(r) = c1_0*[1 + beta_H*S_tilde(r)]
%
%   gamma_D : diffusion-heterogeneity strength
%             D_i^(gamma)(r) = mean(D_i) * z_i(r)^gamma /
%                               mean[z_i(r)^gamma],
%             z_i(r) = D_i^Full(r)/mean(D_i^Full).
%
%   gamma_D = 0 gives uniform diffusion.
%   gamma_D = 1 recovers the finalized Full diffusion field.
%   gamma_D > 1 strengthens the same spatial diffusion contrast while
%   preserving the mean diffusion coefficient.
%
% COUNTERFACTUAL REFERENCE
%   For every beta_H, each gamma_D condition is compared with the
%   same-beta_H uniform-diffusion model (gamma_D=0). Both are independently
%   matched to one fixed global target: the finalized Full model
%   (beta_H=3.40, gamma_D=1) at Tobs=60 min.
%
% PRIMARY TRANSPORT EFFECT
%   Pi_T(beta,gamma) = RMST_60(beta,gamma)
%                    - RMST_60(beta,0).
%
%   Pi_T > 0 : diffusion heterogeneity delays FPT -> supply limitation.
%   Pi_T < 0 : diffusion heterogeneity advances FPT -> retention amplification.
%
% MAIN OUTPUTS
%   A. Pi_T signed heat map
%   B. local burden ratio relative to same-beta uniform diffusion
%   C. change in net diffusive Fe2+ supply
%   D. categorical transport-role map
%
% SI OUTPUTS
%   global-matching audit, c1 rescaling audit, supply-demand balance,
%   Monte Carlo uncertainty, U/H/B/F representative FPT curves, and the
%   critical diffusion-heterogeneity boundary for supply limitation.
%
% IMPORTANT
%   This condition analysis varies reaction demand and diffusion contrast.
%   A later hotspot-scale analysis should instead keep hotspot area,
%   mean, variance, and the global burden fixed while changing spatial
%   fragmentation/correlation length only.
% -------------------------------------------------------------------------

clear; close all; clc;
clear functions;

%% --------------------------- Run configuration --------------------------
% The run mode can be selected without editing this file:
%   set SDRDME_RUN_MODE=pilot
%   matlab -batch "Fig3F_reactionTransportRoleSwitch_RESEARCH_v3"
RUN_MODE = strtrim(getenv('SDRDME_REDOX_RUN_MODE'));
if isempty(RUN_MODE)
    RUN_MODE = 'final';
end
RUN_MODE = lower(RUN_MODE);     % 'smoke' | 'pilot' | 'final'
USE_PARALLEL = true;
AUTO_START_POOL = true;
N_WORKERS = 8;
SAVE_FIG = true;
SAVE_DATA = true;

switch lower(RUN_MODE)
    case 'smoke'
        betaList  = [0.25 1.00 4.00];
        gammaList = [0 1.00 2.00];
        nRepMechanism = 20;
        nRepFPT       = 30;
        nBoot         = 200;

    case 'pilot'
        betaList  = [0.25 0.50 1.00 2.00 4.00 8.00];
        gammaList = [0 0.50 1.00 1.50 2.00];
        nRepMechanism = 60;
        nRepFPT       = 100;
        nBoot         = 500;

    case 'focused'
        % Focused confirmation of the demand-dependent transport ceiling.
        % It keeps the pilot grid but increases trajectory counts enough for
        % stable joint endpoint and threshold-sensitivity inference.
        betaList  = [0.25 0.50 1.00 2.00 4.00 8.00];
        gammaList = [0 0.50 1.00 1.50 2.00];
        nRepMechanism = 150;
        nRepFPT       = 200;
        nBoot         = 1000;

    case 'final'
        % Final grid concentrates resolution around the recycling-buffer
        % transition identified by pilot (c2/c2_0 approximately 1--2),
        % while retaining low- and high-capacity anchors.
        betaList  = [0.25 0.50 0.75 1.00 1.25 1.50 2.00 4.00 8.00];
        gammaList = [0 0.50 1.00 1.50 2.00];
        nRepMechanism = 200;
        nRepFPT       = 300;
        nBoot         = 2000;

    otherwise
        error('Unknown RUN_MODE: %s',RUN_MODE);
end

seedBase = 20260801;
alphaCI = 0.05;

% Operational threshold for the critical supply-limitation boundary.
% This threshold is not a biological constant; it prevents very small but
% statistically nonzero changes from being labelled as a regime transition.
limitationEffectThreshold = 2.0;   % min

scriptDir = fileparts(mfilename('fullpath'));
resultRoot = strtrim(getenv('SDRDME_REDOX_OUTPUT_DIR'));
if isempty(resultRoot)
    resultRoot = fullfile(scriptDir,'condition_redox_buffer_results');
end
outDir = fullfile(resultRoot,RUN_MODE);
if ~exist(outDir,'dir')
    mkdir(outDir);
end

outMain = fullfile(outDir,'Fig3H_redoxBuffer_main');
outSI   = fullfile(outDir,'FigS_redoxBuffer_diagnostics');
outCSV  = fullfile(outDir,'Fig3H_redoxBuffer_results.csv');
outThresholdCSV = fullfile(outDir, ...
    'Fig3H_redoxBuffer_thresholdSensitivity.csv');
outMAT  = fullfile(outDir,'Fig3H_redoxBuffer_data.mat');
checkpointFile = fullfile(outDir,'Fig3H_redoxBuffer_checkpoint.mat');

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

colAmplifier = [0.19 0.42 0.72];
colNeutral   = [0.92 0.92 0.90];
colLimiter   = [0.78 0.20 0.14];

%% ---------------------- Locked Fig. 3 parameters ------------------------
% Geometry
nx = 60;
ny = 36;
aOut = 1.12;  bOut = 0.62;
aIn  = 0.94;  bIn  = 0.49;
ellG = 0.085;
ellS = 0.042;
connMatrix = 1.00;
connLumen  = 0.82;
connNeck   = 0.18;

% Reaction and diffusion parameters, min-based units
D0_Fe2   = 3.0e-4;
D2_scale = 0.70;
D0_Fe3   = D2_scale*D0_Fe2;

c1_0     = 6.0e-3;
betaRef  = 3.40;
c2_0     = 1.8e-2;

% Observation and state
Tobs     = 60;
dtRecord = 0.50;
tvec     = 0:dtRecord:Tobs;

FeTot_per_voxel = 6;
fracFe2_0       = 0.65;

% Matching and local region
dtMoment = 0.05;
nBisect  = 14;
windowRadius = 0.10;
qA = 0.85;

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

Svec = normalizeVector(S(mask));

D1Full = D0_Fe2*G;
D2Full = D0_Fe3*G;
D1Full(~mask) = NaN;
D2Full(~mask) = NaN;

D1FullVec = D1Full(mask);
D2FullVec = D2Full(mask);
D1Mean = mean(D1FullVec);
D2Mean = mean(D2FullVec);

c2vec = c2_0*ones(K,1);

N0   = FeTot_per_voxel*ones(K,1);
x1_0 = round(fracFe2_0*N0);
x2_0 = N0-x1_0;

% Fix one region A_0.85 for the entire condition map. It uses the same
% radius-0.10 window-mean c1-hat convention as candidate screening.
coordX = Xg(mask);
coordY = Yg(mask);
win = buildCircularWindowsSparse(coordX,coordY,windowRadius);

c1RefCenter = normalizeVector(c1_0*(1+betaRef*Svec));
c1RefWindow = (win.W*c1RefCenter)./win.size;
c1RefWindow = normalizeVector(c1RefWindow);

A = topFractionMask(c1RefWindow,1-qA);
nA = nnz(A);

Amap = false(size(mask));
Amap(mask) = A;

%% ---------------------- Fixed global target and threshold ---------------
% The reference Full condition is beta_H=3.40 and gamma_D=1.
c1ReferenceBase = c1_0*(1+betaRef*Svec);

[D1Ref,D2Ref] = buildMeanPreservingDiffusionFields( ...
    D1FullVec,D2FullVec,1.0,mask);

[fromRef1,toRef1,qRef1] = buildDirectedEdgesHarmonic( ...
    mask,idxMap,D1Ref,mean(diff(x)),mean(diff(y)));
[fromRef2,toRef2,qRef2] = buildDirectedEdgesHarmonic( ...
    mask,idxMap,D2Ref,mean(diff(x)),mean(diff(y)));

referenceMoment = momentEulerAll( ...
    c1ReferenceBase,c2vec, ...
    fromRef1,toRef1,qRef1,fromRef2,toRef2,qRef2, ...
    x1_0,x2_0,tvec,dtMoment);

targetGlobal = referenceMoment.globalM3(end);
mStar = targetGlobal/K;
M_A = max(1,round(mStar*nA));
thresholdMultipliers = [0.80 1.00 1.20 1.40];
M_A_list = max(1,round(M_A*thresholdMultipliers));
[~,idxPrimaryThreshold] = min(abs(thresholdMultipliers-1));

fprintf('============================================================\n');
fprintf('Redox-recycling x diffusion-heterogeneity map (%s)\n', ...
    upper(RUN_MODE));
fprintf('Grid                            = %d x %d\n',nx,ny);
fprintf('Omega_R voxels K                = %d\n',K);
fprintf('|A_%.2f|/K                      = %.4f (%d voxels)\n', ...
    qA,nA/K,nA);
fprintf('Window radius                   = %.3f\n',windowRadius);
fprintf('Window-size range               = [%d, %d]\n', ...
    min(win.size),max(win.size));
fprintf('Fixed target at 60 min          = %.6f events\n',targetGlobal);
fprintf('Regional threshold M_A          = %d\n',M_A);
fprintf('M_A/|A|                         = %.4f events/voxel\n',M_A/nA);
fprintf('Threshold multipliers           = %s\n', ...
    mat2str(thresholdMultipliers,2));
fprintf('recycling-capacity conditions   = %d\n',numel(betaList));
fprintf('gamma conditions                = %d\n',numel(gammaList));
fprintf('Mechanism trajectories/cell     = %d\n',nRepMechanism);
fprintf('FPT trajectories/cell           = %d\n',nRepFPT);
fprintf('Common-random-number pairing    = enabled\n');
fprintf('Output directory                = %s\n',outDir);
fprintf('============================================================\n');

%% ------------------ Precompute diffusion models by gamma ----------------
nGamma = numel(gammaList);
Diff = repmat(struct(),1,nGamma);

dx = mean(diff(x));
dy = mean(diff(y));

for ig = 1:nGamma
    gammaD = gammaList(ig);

    [D1Map,D2Map] = buildMeanPreservingDiffusionFields( ...
        D1FullVec,D2FullVec,gammaD,mask);

    [from1,to1,q1] = buildDirectedEdgesHarmonic( ...
        mask,idxMap,D1Map,dx,dy);
    [from2,to2,q2] = buildDirectedEdgesHarmonic( ...
        mask,idxMap,D2Map,dx,dy);

    Diff(ig).gamma = gammaD;
    Diff(ig).D1 = D1Map;
    Diff(ig).D2 = D2Map;
    Diff(ig).from1 = from1;
    Diff(ig).to1 = to1;
    Diff(ig).q1 = q1;
    Diff(ig).from2 = from2;
    Diff(ig).to2 = to2;
    Diff(ig).q2 = q2;
    Diff(ig).out = buildOutgoingLists( ...
        K,from1,to1,q1,from2,to2,q2);
end

idxGamma0 = nearestIndex(gammaList,0);
assert(abs(gammaList(idxGamma0))<1e-12, ...
    'gammaList must contain gamma_D=0.');

%% --------------------------- Allocate outputs ---------------------------
nBeta = numel(betaList);
nThreshold = numel(thresholdMultipliers);

RecyclingFraction = nan(nBeta,nGamma);
DiffusiveSupplyFraction = nan(nBeta,nGamma);
Fe2FinalMean = nan(nBeta,nGamma);
LocalBurdenMean = nan(nBeta,nGamma);
RMSTCondition = nan(nBeta,nGamma);

PiTransport = nan(nBeta,nGamma);
Pi_CI_L = nan(nBeta,nGamma);
Pi_CI_U = nan(nBeta,nGamma);

PiTransportThreshold = nan(nBeta,nGamma,nThreshold);
PiThreshold_CI_L = nan(nBeta,nGamma,nThreshold);
PiThreshold_CI_U = nan(nBeta,nGamma,nThreshold);
P60Gamma0Threshold = nan(nBeta,nGamma,nThreshold);
P60ConditionThreshold = nan(nBeta,nGamma,nThreshold);

BurdenRatio = nan(nBeta,nGamma);
BurdenRatio_CI_L = nan(nBeta,nGamma);
BurdenRatio_CI_U = nan(nBeta,nGamma);

DeltaNetSupply = nan(nBeta,nGamma);
DeltaNetSupply_CI_L = nan(nBeta,nGamma);
DeltaNetSupply_CI_U = nan(nBeta,nGamma);

DeltaSupplyDemand = nan(nBeta,nGamma);

AlphaCondition = nan(nBeta,nGamma);
AlphaRatioToGamma0 = nan(nBeta,nGamma);

GlobalMismatch = nan(nBeta,nGamma);
RMST_CI_HalfWidth = nan(nBeta,nGamma);
MassBalanceResidual = nan(nBeta,nGamma);

P60Gamma0 = nan(nBeta,nGamma);
P60Condition = nan(nBeta,nGamma);

% -1 retention amplifier, 0 neutral/mixed, +1 supply limiter.
TransportRole = zeros(nBeta,nGamma);

FPT_gamma0_cell = cell(nBeta,1);
FPT_condition_cell = cell(nBeta,nGamma);
FPT_gamma0_multithreshold_cell = cell(nBeta,1);
FPT_condition_multithreshold_cell = cell(nBeta,nGamma);
Mechanism_condition_cell = cell(nBeta,nGamma);
completedRecyclingRows = false(nBeta,1);

%% --------------------------- Condition scan -----------------------------
for ib = 1:nBeta
    recyclingScale = betaList(ib);
    c1Base = c1_0*(1+betaRef*Svec);
    c2vec = (c2_0*recyclingScale)*ones(K,1);

    fprintf('\n============================================================\n');
    fprintf('recycling capacity c2/c2_0 = %.3f (%d/%d)\n', ...
        recyclingScale,ib,nBeta);
    fprintf('============================================================\n');

    % Same-beta uniform-diffusion reference, independently matched to the
    % one fixed global target.
    [alpha0,moment0] = matchAlphaToTarget( ...
        c1Base,c2vec,Diff(idxGamma0), ...
        x1_0,x2_0,tvec,dtMoment,targetGlobal,nBisect);

    c1Gamma0 = alpha0*c1Base;

    mechanismSeed = seedBase + 100000*ib + 1000;
    fptSeed       = seedBase + 100000*ib + 2000;

    baseMechanism = runRegionalMechanismEnsembleNSM( ...
        c1Gamma0,c2vec,Diff(idxGamma0).out,A, ...
        x1_0,x2_0,Tobs,nRepMechanism,mechanismSeed,usePar);

    FPT0All = runRegionFPTensembleNSM( ...
        c1Gamma0,c2vec,Diff(idxGamma0).out,A, ...
        x1_0,x2_0,M_A_list,Tobs,nRepFPT,fptSeed,usePar);
    FPT0 = FPT0All(:,idxPrimaryThreshold);

    FPT_gamma0_cell{ib} = FPT0;
    FPT_gamma0_multithreshold_cell{ib} = FPT0All;

    for ig = 1:nGamma
        gammaD = gammaList(ig);

        fprintf('  gamma_D = %.3f (%d/%d)\n',gammaD,ig,nGamma);

        [alphaG,momentG] = matchAlphaToTarget( ...
            c1Base,c2vec,Diff(ig), ...
            x1_0,x2_0,tvec,dtMoment,targetGlobal,nBisect);

        AlphaCondition(ib,ig) = alphaG;
        AlphaRatioToGamma0(ib,ig) = alphaG/alpha0;

        if ig==idxGamma0
            conditionMechanism = baseMechanism;
            FPTGAll = FPT0All;
            FPTG = FPT0;

            PiTransport(ib,ig) = 0;
            Pi_CI_L(ib,ig) = 0;
            Pi_CI_U(ib,ig) = 0;

            BurdenRatio(ib,ig) = 1;
            BurdenRatio_CI_L(ib,ig) = 1;
            BurdenRatio_CI_U(ib,ig) = 1;

            DeltaNetSupply(ib,ig) = 0;
            DeltaNetSupply_CI_L(ib,ig) = 0;
            DeltaNetSupply_CI_U(ib,ig) = 0;

        else
            c1Gamma = alphaG*c1Base;

            % Use the same seed sequence as gamma_D=0 to implement common
            % random numbers and reduce pairwise Monte Carlo noise.
            conditionMechanism = runRegionalMechanismEnsembleNSM( ...
                c1Gamma,c2vec,Diff(ig).out,A, ...
                x1_0,x2_0,Tobs,nRepMechanism,mechanismSeed,usePar);

            FPTGAll = runRegionFPTensembleNSM( ...
                c1Gamma,c2vec,Diff(ig).out,A, ...
                x1_0,x2_0,M_A_list,Tobs,nRepFPT,fptSeed,usePar);
            FPTG = FPTGAll(:,idxPrimaryThreshold);

            [BurdenRatio(ib,ig),burdenCI] = bootstrapPairedRatio( ...
                conditionMechanism.localBurdenFinal, ...
                baseMechanism.localBurdenFinal, ...
                nBoot,alphaCI, ...
                seedBase+400000+10000*ib+ig);

            BurdenRatio_CI_L(ib,ig) = burdenCI(1);
            BurdenRatio_CI_U(ib,ig) = burdenCI(2);

            [DeltaNetSupply(ib,ig),supplyCI] = ...
                bootstrapPairedDifference( ...
                conditionMechanism.netDiffFinal, ...
                baseMechanism.netDiffFinal, ...
                nBoot,alphaCI, ...
                seedBase+500000+10000*ib+ig);

            DeltaNetSupply_CI_L(ib,ig) = supplyCI(1);
            DeltaNetSupply_CI_U(ib,ig) = supplyCI(2);
        end

        for ith = 1:nThreshold
            if ig==idxGamma0
                piThreshold = 0;
                piThresholdCI = [0 0];
            else
                [piThreshold,piThresholdCI] = bootstrapPairedDifference( ...
                    min(FPTGAll(:,ith),Tobs), ...
                    min(FPT0All(:,ith),Tobs), ...
                    nBoot,alphaCI, ...
                    seedBase+300000+10000*ib+100*ig+ith);
            end

            PiTransportThreshold(ib,ig,ith) = piThreshold;
            PiThreshold_CI_L(ib,ig,ith) = piThresholdCI(1);
            PiThreshold_CI_U(ib,ig,ith) = piThresholdCI(2);
            P60Gamma0Threshold(ib,ig,ith) = ...
                mean(FPT0All(:,ith)<=Tobs);
            P60ConditionThreshold(ib,ig,ith) = ...
                mean(FPTGAll(:,ith)<=Tobs);
        end

        PiTransport(ib,ig) = ...
            PiTransportThreshold(ib,ig,idxPrimaryThreshold);
        Pi_CI_L(ib,ig) = ...
            PiThreshold_CI_L(ib,ig,idxPrimaryThreshold);
        Pi_CI_U(ib,ig) = ...
            PiThreshold_CI_U(ib,ig,idxPrimaryThreshold);

        FPT_condition_cell{ib,ig} = FPTG;
        FPT_condition_multithreshold_cell{ib,ig} = FPTGAll;
        Mechanism_condition_cell{ib,ig} = conditionMechanism;
        RecyclingFraction(ib,ig) = safeMeanRatio( ...
            conditionMechanism.NredFinal, ...
            conditionMechanism.NfenFinal);
        DiffusiveSupplyFraction(ib,ig) = safeMeanRatio( ...
            conditionMechanism.netDiffFinal, ...
            conditionMechanism.NfenFinal);
        Fe2FinalMean(ib,ig) = mean(conditionMechanism.Fe2Final);
        LocalBurdenMean(ib,ig) = mean( ...
            conditionMechanism.localBurdenFinal);
        RMSTCondition(ib,ig) = mean(min(FPTG,Tobs));

        % Supply-demand ratio:
        % (net diffusive Fe2+ supply + R2 regeneration) / R1 consumption.
        lambda0 = safeMeanRatio( ...
            baseMechanism.netDiffFinal+baseMechanism.NredFinal, ...
            baseMechanism.NfenFinal);

        lambdaG = safeMeanRatio( ...
            conditionMechanism.netDiffFinal+conditionMechanism.NredFinal, ...
            conditionMechanism.NfenFinal);

        DeltaSupplyDemand(ib,ig) = lambdaG-lambda0;

        P60Gamma0(ib,ig) = mean(FPT0<=Tobs);
        P60Condition(ib,ig) = mean(FPTG<=Tobs);

        mismatchBase = abs(mean(baseMechanism.globalBurdenFinal) ...
            - targetGlobal/K)/(targetGlobal/K);

        mismatchG = abs(mean(conditionMechanism.globalBurdenFinal) ...
            - targetGlobal/K)/(targetGlobal/K);

        GlobalMismatch(ib,ig) = max(mismatchBase,mismatchG);

        RMST_CI_HalfWidth(ib,ig) = ...
            0.5*(Pi_CI_U(ib,ig)-Pi_CI_L(ib,ig));

        MassBalanceResidual(ib,ig) = max( ...
            max(abs(baseMechanism.balanceResidual)), ...
            max(abs(conditionMechanism.balanceResidual)));

        % Conservative role classification requiring consistent changes in
        % both the time-to-threshold endpoint and local burden.
        if Pi_CI_L(ib,ig)>0 && BurdenRatio_CI_U(ib,ig)<1
            TransportRole(ib,ig) = 1;   % supply limiter
        elseif Pi_CI_U(ib,ig)<0 && BurdenRatio_CI_L(ib,ig)>1
            TransportRole(ib,ig) = -1;  % retention amplifier
        else
            TransportRole(ib,ig) = 0;   % neutral / mixed / uncertain
        end

        fprintf(['    alpha ratio=%.4f; moment ratio=%.6f; ' ...
                 'SSA mismatch=%.3f%%\n'], ...
            AlphaRatioToGamma0(ib,ig), ...
            momentG.globalM3(end)/targetGlobal, ...
            100*GlobalMismatch(ib,ig));

        fprintf(['    Pi_T=%.3f [%.3f, %.3f] min; ' ...
                 'burden ratio=%.3f [%.3f, %.3f]\n'], ...
            PiTransport(ib,ig),Pi_CI_L(ib,ig),Pi_CI_U(ib,ig), ...
            BurdenRatio(ib,ig), ...
            BurdenRatio_CI_L(ib,ig),BurdenRatio_CI_U(ib,ig));
    end

    completedRecyclingRows(ib) = true;
    if SAVE_DATA
        save(checkpointFile, ...
            'RUN_MODE','betaList','gammaList', ...
            'thresholdMultipliers','M_A_list', ...
            'completedRecyclingRows', ...
            'PiTransport','Pi_CI_L','Pi_CI_U', ...
            'PiTransportThreshold','PiThreshold_CI_L', ...
            'PiThreshold_CI_U','P60Gamma0Threshold', ...
            'P60ConditionThreshold','BurdenRatio', ...
            'BurdenRatio_CI_L','BurdenRatio_CI_U', ...
            'DeltaNetSupply','DeltaNetSupply_CI_L', ...
            'DeltaNetSupply_CI_U','DeltaSupplyDemand', ...
            'AlphaCondition','AlphaRatioToGamma0', ...
            'GlobalMismatch','RMST_CI_HalfWidth', ...
            'MassBalanceResidual','P60Gamma0','P60Condition', ...
            'TransportRole','RecyclingFraction', ...
            'DiffusiveSupplyFraction','Fe2FinalMean', ...
            'LocalBurdenMean','RMSTCondition', ...
            'FPT_gamma0_cell','FPT_condition_cell', ...
            'FPT_gamma0_multithreshold_cell', ...
            'FPT_condition_multithreshold_cell', ...
            'Mechanism_condition_cell','-v7.3');
        fprintf('  checkpoint saved after recycling row %d/%d\n', ...
            ib,nBeta);
    end
end

%% --------------------- Critical limitation boundary ---------------------
criticalGamma = nan(nBeta,1);

for ib = 1:nBeta
    eligible = Pi_CI_L(ib,:)>0 ...
        & PiTransport(ib,:)>=limitationEffectThreshold ...
        & BurdenRatio_CI_U(ib,:)<1;

    idx = find(eligible,1,'first');
    if ~isempty(idx)
        criticalGamma(ib) = gammaList(idx);
    end
end

%% ---------------------------- Main figure -------------------------------
figMain = figure('Units','centimeters','Position',[1.0 1.0 18.8 9.3], ...
    'Color','w');
tl = tiledlayout(figMain,2,2,'TileSpacing','compact','Padding','compact');

% A. Signed transport effect on regional FPT.
ax = nexttile(tl,1);
plotCellHeatmap(ax,gammaList,betaList,PiTransport);
pMax = max(abs(PiTransport(:)),[],'omitnan');
pMax = max(pMax,1);
caxis(ax,[-pMax pMax]);
colormap(ax,cmapBlueWhiteRed(256));
cb = colorbar(ax);
ylabel(cb,'\Pi_T = RMST_{\gamma_D}-RMST_{0} (min)');
hold(ax,'on');
contour(ax,gammaList,betaList,PiTransport,[0 2 5 10], ...
    'LineColor',[0.12 0.12 0.12],'LineWidth',0.60, ...
    'ShowText','on');
xlabel(ax,'Diffusion heterogeneity, \gamma_D');
ylabel(ax,'Recycling capacity, c_2/c_2^0');
title(ax,'Recycling-dependent transport effect on FPT', ...
    'FontSize',8.5,'FontWeight','normal');
text(ax,0.02,0.97,'A','Units','normalized','FontWeight','bold', ...
    'FontSize',12,'VerticalAlignment','top');
formatAxes(ax);

% B. Local burden ratio relative to same-beta uniform diffusion.
ax = nexttile(tl,2);
plotCellHeatmap(ax,gammaList,betaList,BurdenRatio);
bDev = max(abs(BurdenRatio(:)-1),[],'omitnan');
bDev = max(bDev,0.05);
caxis(ax,[1-bDev 1+bDev]);
colormap(ax,cmapBlueWhiteRed(256));
cb = colorbar(ax);
ylabel(cb,'B_A(\gamma_D)/B_A(0)');
hold(ax,'on');
contour(ax,gammaList,betaList,BurdenRatio,[0.8 0.9 1.0 1.1 1.2], ...
    'LineColor',[0.12 0.12 0.12],'LineWidth',0.60, ...
    'ShowText','on');
xlabel(ax,'Diffusion heterogeneity, \gamma_D');
ylabel(ax,'Recycling capacity, c_2/c_2^0');
title(ax,'Recycling-dependent local burden', ...
    'FontSize',8.5,'FontWeight','normal');
text(ax,0.02,0.97,'B','Units','normalized','FontWeight','bold', ...
    'FontSize',12,'VerticalAlignment','top');
formatAxes(ax);

% C. Net diffusive Fe2+ supply change.
ax = nexttile(tl,3);
plotCellHeatmap(ax,gammaList,betaList,DeltaNetSupply);
sMax = max(abs(DeltaNetSupply(:)),[],'omitnan');
sMax = max(sMax,0.02);
caxis(ax,[-sMax sMax]);
colormap(ax,cmapBlueWhiteRed(256));
cb = colorbar(ax);
ylabel(cb,'\Delta(N_{in}-N_{out}) per voxel in A');
hold(ax,'on');
contour(ax,gammaList,betaList,DeltaNetSupply,[0], ...
    'LineColor',[0.12 0.12 0.12],'LineWidth',0.75);
xlabel(ax,'Diffusion heterogeneity, \gamma_D');
ylabel(ax,'Recycling capacity, c_2/c_2^0');
title(ax,'Transport effect on net Fe^{2+} supply', ...
    'FontSize',8.5,'FontWeight','normal');
text(ax,0.02,0.97,'C','Units','normalized','FontWeight','bold', ...
    'FontSize',12,'VerticalAlignment','top');
formatAxes(ax);

% D. Categorical transport role.
ax = nexttile(tl,4);
plotCellHeatmap(ax,gammaList,betaList,TransportRole);
caxis(ax,[-1 1]);
colormap(ax,[colAmplifier;colNeutral;colLimiter]);
cb = colorbar(ax);
cb.Ticks = [-1 0 1];
cb.TickLabels = {'retention amplifier','neutral / mixed','supply limiter'};
hold(ax,'on');
xlabel(ax,'Diffusion heterogeneity, \gamma_D');
ylabel(ax,'Recycling capacity, c_2/c_2^0');
title(ax,'Transport-role classification', ...
    'FontSize',8.5,'FontWeight','normal');
text(ax,0.02,0.97,'D','Units','normalized','FontWeight','bold', ...
    'FontSize',12,'VerticalAlignment','top');
formatAxes(ax);

if SAVE_FIG
    exportgraphics(figMain,[outMain '.pdf'],'ContentType','vector');
    exportgraphics(figMain,[outMain '.png'],'Resolution',600);
    exportgraphics(figMain,[outMain '.tif'],'Resolution',600);
    savefig(figMain,[outMain '.fig']);
end

%% ---------------------- Supplementary diagnostics -----------------------
figSI = figure('Units','centimeters','Position',[1.0 1.0 19.0 11.5], ...
    'Color','w');
tlSI = tiledlayout(figSI,2,3,'TileSpacing','compact','Padding','compact');

% A. SSA global-target audit.
ax = nexttile(tlSI,1);
plotCellHeatmap(ax,gammaList,betaList,100*GlobalMismatch);
colormap(ax,natureBlueYellow(256));
cb = colorbar(ax);
ylabel(cb,'Maximum target mismatch (%)');
xlabel(ax,'\gamma_D');
ylabel(ax,'c_2/c_2^0');
title(ax,'Global-burden matching audit', ...
    'FontSize',8.3,'FontWeight','normal');
text(ax,0.02,0.97,'A','Units','normalized','FontWeight','bold', ...
    'FontSize',11.5,'VerticalAlignment','top');
formatAxes(ax);

% B. Reaction-rate rescaling required by the fixed global budget.
ax = nexttile(tlSI,2);
plotCellHeatmap(ax,gammaList,betaList,AlphaRatioToGamma0);
colormap(ax,natureBlueYellow(256));
cb = colorbar(ax);
ylabel(cb,'\alpha(\gamma_D)/\alpha(0)');
xlabel(ax,'\gamma_D');
ylabel(ax,'c_2/c_2^0');
title(ax,'Global-budget c_1 rescaling', ...
    'FontSize',8.3,'FontWeight','normal');
text(ax,0.02,0.97,'B','Units','normalized','FontWeight','bold', ...
    'FontSize',11.5,'VerticalAlignment','top');
formatAxes(ax);

% C. Supply-demand balance change.
ax = nexttile(tlSI,3);
plotCellHeatmap(ax,gammaList,betaList,DeltaSupplyDemand);
sdMax = max(abs(DeltaSupplyDemand(:)),[],'omitnan');
sdMax = max(sdMax,0.02);
caxis(ax,[-sdMax sdMax]);
colormap(ax,cmapBlueWhiteRed(256));
cb = colorbar(ax);
ylabel(cb,'\Delta[(net diffusion + R2)/R1]');
xlabel(ax,'\gamma_D');
ylabel(ax,'c_2/c_2^0');
title(ax,'Fe^{2+} supply--demand balance', ...
    'FontSize',8.3,'FontWeight','normal');
text(ax,0.02,0.97,'C','Units','normalized','FontWeight','bold', ...
    'FontSize',11.5,'VerticalAlignment','top');
formatAxes(ax);

% D. Monte Carlo uncertainty.
ax = nexttile(tlSI,4);
plotCellHeatmap(ax,gammaList,betaList,RMST_CI_HalfWidth);
colormap(ax,natureBlueYellow(256));
cb = colorbar(ax);
ylabel(cb,'95% CI half-width (min)');
xlabel(ax,'\gamma_D');
ylabel(ax,'c_2/c_2^0');
title(ax,'RMST uncertainty', ...
    'FontSize',8.3,'FontWeight','normal');
text(ax,0.02,0.97,'D','Units','normalized','FontWeight','bold', ...
    'FontSize',11.5,'VerticalAlignment','top');
formatAxes(ax);

% E. Low/high recycling anchor curves.
ax = nexttile(tlSI,5); hold(ax,'on');

idxBeta0 = 1;
idxBetaRef = nBeta;
idxGamma1 = nearestIndex(gammaList,1);

[tu,Fu] = censoredCrossingCDF( ...
    FPT_condition_cell{idxBeta0,idxGamma0},Tobs);
[tb,Fb] = censoredCrossingCDF( ...
    FPT_condition_cell{idxBeta0,idxGamma1},Tobs);
[th,Fh] = censoredCrossingCDF( ...
    FPT_condition_cell{idxBetaRef,idxGamma0},Tobs);
[tf,Ff] = censoredCrossingCDF( ...
    FPT_condition_cell{idxBetaRef,idxGamma1},Tobs);

stairs(ax,tu,Fu,'Color',colU,'LineWidth',1.10);
stairs(ax,th,Fh,'Color',colH,'LineWidth',1.10);
stairs(ax,tb,Fb,'Color',colB,'LineWidth',1.10);
stairs(ax,tf,Ff,'Color',colF,'LineWidth',1.10);

xlabel(ax,'Regional FPT (min)');
ylabel(ax,'Cumulative crossing probability');
title(ax,'Low versus high recycling capacity', ...
    'FontSize',8.3,'FontWeight','normal');
legend(ax,{'low recycling, \gamma_D=0', ...
    'high recycling, \gamma_D=0', ...
    'low recycling, \gamma_D=1', ...
    'high recycling, \gamma_D=1'},'Location','northwest', ...
    'Box','off','FontSize',5.7);
text(ax,0.02,0.97,'E','Units','normalized','FontWeight','bold', ...
    'FontSize',11.5,'VerticalAlignment','top');
xlim(ax,[0 Tobs]); ylim(ax,[0 1]);
formatAxes(ax);

% F. Critical supply-limitation boundary.
ax = nexttile(tlSI,6); hold(ax,'on');
plot(ax,betaList,criticalGamma,'-o','Color',colLimiter, ...
    'MarkerFaceColor','w','MarkerSize',4.0,'LineWidth',1.05);
plot(ax,[min(betaList) max(betaList)],[1 1],':', ...
    'Color',colF,'LineWidth',0.8);
xlabel(ax,'Recycling capacity, c_2/c_2^0');
ylabel(ax,'Minimum limiting \gamma_D');
title(ax,'Critical supply-limitation boundary', ...
    'FontSize',8.3,'FontWeight','normal');
text(ax,0.02,0.97,'F','Units','normalized','FontWeight','bold', ...
    'FontSize',11.5,'VerticalAlignment','top');
text(ax,0.98,0.07,sprintf(['criterion: lower CI > 0, burden upper CI < 1, ' ...
    '\\Pi_T \\geq %.0f min'],limitationEffectThreshold), ...
    'Units','normalized','HorizontalAlignment','right', ...
    'FontSize',5.2);
xlim(ax,[min(betaList) max(betaList)]);
ylim(ax,[min(gammaList)-0.05 max(gammaList)+0.05]);
formatAxes(ax);

if SAVE_FIG
    exportgraphics(figSI,[outSI '.pdf'],'ContentType','vector');
    exportgraphics(figSI,[outSI '.png'],'Resolution',600);
    exportgraphics(figSI,[outSI '.tif'],'Resolution',600);
    savefig(figSI,[outSI '.fig']);
end

%% --------------------------- Export results -----------------------------
[BetaGrid,GammaGrid] = ndgrid(betaList,gammaList);

resultTable = table( ...
    BetaGrid(:),GammaGrid(:), ...
    AlphaCondition(:),AlphaRatioToGamma0(:), ...
    GlobalMismatch(:), ...
    PiTransport(:),Pi_CI_L(:),Pi_CI_U(:), ...
    BurdenRatio(:),BurdenRatio_CI_L(:),BurdenRatio_CI_U(:), ...
    DeltaNetSupply(:), ...
    DeltaNetSupply_CI_L(:),DeltaNetSupply_CI_U(:), ...
    DeltaSupplyDemand(:), ...
    P60Gamma0(:),P60Condition(:), ...
    RMST_CI_HalfWidth(:),MassBalanceResidual(:), ...
    TransportRole(:),RecyclingFraction(:), ...
    DiffusiveSupplyFraction(:),Fe2FinalMean(:), ...
    LocalBurdenMean(:),RMSTCondition(:), ...
    'VariableNames',{ ...
    'RecyclingCapacityScale','DiffusionHeterogeneityGamma', ...
    'AlphaCondition','AlphaRatioToGamma0', ...
    'GlobalMismatch', ...
    'PiTransport_RMSTDifference', ...
    'PiTransport_CI_L','PiTransport_CI_U', ...
    'LocalBurdenRatio','LocalBurdenRatio_CI_L', ...
    'LocalBurdenRatio_CI_U', ...
    'DeltaNetDiffusiveSupply', ...
    'DeltaNetDiffusiveSupply_CI_L', ...
    'DeltaNetDiffusiveSupply_CI_U', ...
    'DeltaSupplyDemandRatio', ...
    'P60_Gamma0','P60_Condition', ...
    'RMST_CI_HalfWidth','MassBalanceResidual', ...
    'TransportRole','RecyclingFraction_R2overR1', ...
    'DiffusiveSupplyFraction_overR1','Fe2FinalMean', ...
    'LocalBurdenMean','RMSTCondition'});

[BetaThresholdGrid,GammaThresholdGrid,ThresholdMultiplierGrid] = ...
    ndgrid(betaList,gammaList,thresholdMultipliers);
ThresholdCountGrid = repmat(reshape(M_A_list,1,1,[]), ...
    nBeta,nGamma,1);

thresholdTable = table( ...
    BetaThresholdGrid(:),GammaThresholdGrid(:), ...
    ThresholdMultiplierGrid(:), ...
    ThresholdCountGrid(:), ...
    PiTransportThreshold(:),PiThreshold_CI_L(:),PiThreshold_CI_U(:), ...
    P60Gamma0Threshold(:),P60ConditionThreshold(:), ...
    'VariableNames',{ ...
    'RecyclingCapacityScale','DiffusionHeterogeneityGamma', ...
    'ThresholdMultiplier','ThresholdCount', ...
    'PiTransport_RMSTDifference', ...
    'PiTransport_CI_L','PiTransport_CI_U', ...
    'P60_Gamma0','P60_Condition'});

if SAVE_DATA
    writetable(resultTable,outCSV);
    writetable(thresholdTable,outThresholdCSV);

    metadata = struct();
    metadata.RUN_MODE = RUN_MODE;
    metadata.seedBase = seedBase;
    metadata.betaList = betaList;
    metadata.recyclingCapacityScales = betaList;
    metadata.gammaList = gammaList;
    metadata.nRepMechanism = nRepMechanism;
    metadata.nRepFPT = nRepFPT;
    metadata.nBoot = nBoot;
    metadata.outputDirectory = outDir;
    metadata.windowRadius = windowRadius;
    metadata.qA = qA;
    metadata.nA = nA;
    metadata.M_A = M_A;
    metadata.Tobs = Tobs;
    metadata.thresholdMultipliers = thresholdMultipliers;
    metadata.M_A_list = M_A_list;
    metadata.targetGlobal = targetGlobal;
    metadata.limitationEffectThreshold = limitationEffectThreshold;
    metadata.commonRandomNumbers = true;
    metadata.betaHotspot = betaRef;

    save(outMAT,'metadata','resultTable','thresholdTable', ...
        'betaList','gammaList','PiTransport','Pi_CI_L','Pi_CI_U', ...
        'thresholdMultipliers','M_A_list', ...
        'PiTransportThreshold','PiThreshold_CI_L','PiThreshold_CI_U', ...
        'P60Gamma0Threshold','P60ConditionThreshold', ...
        'BurdenRatio','BurdenRatio_CI_L','BurdenRatio_CI_U', ...
        'DeltaNetSupply','DeltaNetSupply_CI_L','DeltaNetSupply_CI_U', ...
        'DeltaSupplyDemand','AlphaCondition','AlphaRatioToGamma0', ...
        'GlobalMismatch','RMST_CI_HalfWidth','MassBalanceResidual', ...
        'RecyclingFraction','DiffusiveSupplyFraction', ...
        'Fe2FinalMean','LocalBurdenMean','RMSTCondition', ...
        'P60Gamma0','P60Condition','TransportRole','criticalGamma', ...
        'FPT_gamma0_cell','FPT_condition_cell', ...
        'FPT_gamma0_multithreshold_cell', ...
        'FPT_condition_multithreshold_cell', ...
        'Mechanism_condition_cell', ...
        'A','Amap','c1RefWindow','win','mask','x','y','-v7.3');
end

fprintf('\nCompleted. Outputs:\n');
fprintf('  %s.{pdf,png,tif,fig}\n',outMain);
fprintf('  %s.{pdf,png,tif,fig}\n',outSI);
fprintf('  %s\n',outCSV);
fprintf('  %s\n',outThresholdCSV);
fprintf('  %s\n',outMAT);

%% ======================== Analysis-specific functions ===================

function [D1Map,D2Map] = buildMeanPreservingDiffusionFields( ...
    D1FullVec,D2FullVec,gammaD,mask)

    D1Mean = mean(D1FullVec);
    D2Mean = mean(D2FullVec);

    z1 = D1FullVec/D1Mean;
    z2 = D2FullVec/D2Mean;

    shape1 = z1.^gammaD;
    shape2 = z2.^gammaD;

    D1Vec = D1Mean*shape1/mean(shape1);
    D2Vec = D2Mean*shape2/mean(shape2);

    D1Map = NaN(size(mask));
    D2Map = NaN(size(mask));

    D1Map(mask) = D1Vec;
    D2Map(mask) = D2Vec;
end

function [alpha,momentOut] = matchAlphaToTarget( ...
    c1Base,c2,diffStruct,x1_0,x2_0,tvec,dtMoment, ...
    targetGlobal,nBisect)

    lo = 0.05;
    hi = 8.00;

    momLo = momentEulerAll( ...
        lo*c1Base,c2, ...
        diffStruct.from1,diffStruct.to1,diffStruct.q1, ...
        diffStruct.from2,diffStruct.to2,diffStruct.q2, ...
        x1_0,x2_0,tvec,dtMoment);

    momHi = momentEulerAll( ...
        hi*c1Base,c2, ...
        diffStruct.from1,diffStruct.to1,diffStruct.q1, ...
        diffStruct.from2,diffStruct.to2,diffStruct.q2, ...
        x1_0,x2_0,tvec,dtMoment);

    while momLo.globalM3(end)>targetGlobal
        hi = lo;
        lo = 0.5*lo;

        momLo = momentEulerAll( ...
            lo*c1Base,c2, ...
            diffStruct.from1,diffStruct.to1,diffStruct.q1, ...
            diffStruct.from2,diffStruct.to2,diffStruct.q2, ...
            x1_0,x2_0,tvec,dtMoment);
    end

    while momHi.globalM3(end)<targetGlobal
        lo = hi;
        hi = 2*hi;

        momHi = momentEulerAll( ...
            hi*c1Base,c2, ...
            diffStruct.from1,diffStruct.to1,diffStruct.q1, ...
            diffStruct.from2,diffStruct.to2,diffStruct.q2, ...
            x1_0,x2_0,tvec,dtMoment);

        if hi>100
            error('Could not bracket the fixed global target.');
        end
    end

    for it = 1:nBisect
        mid = 0.5*(lo+hi);

        trial = momentEulerAll( ...
            mid*c1Base,c2, ...
            diffStruct.from1,diffStruct.to1,diffStruct.q1, ...
            diffStruct.from2,diffStruct.to2,diffStruct.q2, ...
            x1_0,x2_0,tvec,dtMoment);

        if trial.globalM3(end)<targetGlobal
            lo = mid;
        else
            hi = mid;
        end
    end

    alpha = 0.5*(lo+hi);

    momentOut = momentEulerAll( ...
        alpha*c1Base,c2, ...
        diffStruct.from1,diffStruct.to1,diffStruct.q1, ...
        diffStruct.from2,diffStruct.to2,diffStruct.q2, ...
        x1_0,x2_0,tvec,dtMoment);
end

function out = runRegionalMechanismEnsembleNSM( ...
    c1,c2,outgoing,A,x1_0,x2_0,Tobs,nRep,seedBase,usePar)

    globalBurdenFinal = zeros(nRep,1);
    localBurdenFinal  = zeros(nRep,1);
    netDiffFinal      = zeros(nRep,1);
    NredFinal         = zeros(nRep,1);
    NfenFinal         = zeros(nRep,1);
    Fe2Final          = zeros(nRep,1);
    balanceResidual   = zeros(nRep,1);

    if usePar
        parfor r = 1:nRep
            rng(seedBase+r,'twister');

            one = simulateRegionalMechanismOneNSM( ...
                c1,c2,outgoing,A,x1_0,x2_0,Tobs);

            globalBurdenFinal(r) = one.globalBurdenFinal;
            localBurdenFinal(r)  = one.localBurdenFinal;
            netDiffFinal(r)      = one.netDiffFinal;
            NredFinal(r)         = one.NredFinal;
            NfenFinal(r)         = one.NfenFinal;
            Fe2Final(r)          = one.Fe2Final;
            balanceResidual(r)   = one.balanceResidual;
        end
    else
        for r = 1:nRep
            rng(seedBase+r,'twister');

            one = simulateRegionalMechanismOneNSM( ...
                c1,c2,outgoing,A,x1_0,x2_0,Tobs);

            globalBurdenFinal(r) = one.globalBurdenFinal;
            localBurdenFinal(r)  = one.localBurdenFinal;
            netDiffFinal(r)      = one.netDiffFinal;
            NredFinal(r)         = one.NredFinal;
            NfenFinal(r)         = one.NfenFinal;
            Fe2Final(r)          = one.Fe2Final;
            balanceResidual(r)   = one.balanceResidual;
        end
    end

    out.globalBurdenFinal = globalBurdenFinal;
    out.localBurdenFinal  = localBurdenFinal;
    out.netDiffFinal      = netDiffFinal;
    out.NredFinal         = NredFinal;
    out.NfenFinal         = NfenFinal;
    out.Fe2Final          = Fe2Final;
    out.balanceResidual   = balanceResidual;
end

function one = simulateRegionalMechanismOneNSM( ...
    c1,c2,out,A,x1_0,x2_0,Tobs)

    K = numel(c1);
    nA = nnz(A);

    x1 = x1_0(:);
    x2 = x2_0(:);

    initialFe2 = sum(x1(A))/nA;

    globalROS = 0;
    localROS  = 0;
    Nin = 0;
    Nout = 0;
    Nred = 0;
    Nfen = 0;

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
            oldRate = voxelRate(k);
            voxelRate(k) = 0;
            tree = fenwickAdd(tree,k,-oldRate);
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

    finalFe2 = sum(x1(A))/nA;
    netDiff = (Nin-Nout)/nA;
    nredDensity = Nred/nA;
    nfenDensity = Nfen/nA;

    one.globalBurdenFinal = globalROS/K;
    one.localBurdenFinal  = localROS/nA;
    one.netDiffFinal      = netDiff;
    one.NredFinal         = nredDensity;
    one.NfenFinal         = nfenDensity;
    one.Fe2Final          = finalFe2;

    one.balanceResidual = ...
        (finalFe2-initialFe2) ...
        - (netDiff+nredDensity-nfenDensity);
end

function FPT = runRegionFPTensembleNSM( ...
    c1,c2,outgoing,A,x1_0,x2_0,M_A,Tobs,nRep,seedBase,usePar)

    M_A = M_A(:).';
    assert(all(diff(M_A)>0), ...
        'Regional FPT thresholds must be strictly increasing.');
    FPT = inf(nRep,numel(M_A));

    if usePar
        parfor r = 1:nRep
            rng(seedBase+r,'twister');

            FPT(r,:) = simulateRegionFPTOneNSM( ...
                c1,c2,outgoing,A,x1_0,x2_0,M_A,Tobs);
        end
    else
        for r = 1:nRep
            rng(seedBase+r,'twister');

            FPT(r,:) = simulateRegionFPTOneNSM( ...
                c1,c2,outgoing,A,x1_0,x2_0,M_A,Tobs);
        end
    end
end

function fpt = simulateRegionFPTOneNSM( ...
    c1,c2,out,A,x1_0,x2_0,M_A,Tobs)

    x1 = x1_0(:);
    x2 = x2_0(:);

    M_A = M_A(:).';
    nThreshold = numel(M_A);
    nextThreshold = 1;
    localROS = 0;
    fpt = inf(1,nThreshold);

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
            oldRate = voxelRate(k);
            voxelRate(k) = 0;
            tree = fenwickAdd(tree,k,-oldRate);
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

                while nextThreshold<=nThreshold ...
                        && localROS>=M_A(nextThreshold)
                    fpt(nextThreshold) = tNext;
                    nextThreshold = nextThreshold+1;
                end

                if nextThreshold>nThreshold
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

function [delta,ci] = bootstrapPairedDifference(x,y,nBoot,alpha,seed)
    x = x(:);
    y = y(:);

    assert(numel(x)==numel(y), ...
        'Paired vectors must have equal lengths.');

    d = x-y;
    n = numel(d);

    delta = mean(d);
    boot = zeros(nBoot,1);

    rng(seed,'twister');

    for b = 1:nBoot
        idx = randi(n,n,1);
        boot(b) = mean(d(idx));
    end

    ci = [quantileLocal(boot,alpha/2), ...
          quantileLocal(boot,1-alpha/2)];
end

function [ratio,ci] = bootstrapPairedRatio(x,y,nBoot,alpha,seed)
    x = x(:);
    y = y(:);

    assert(numel(x)==numel(y), ...
        'Paired vectors must have equal lengths.');

    n = numel(x);
    ratio = mean(x)/max(mean(y),eps);

    boot = zeros(nBoot,1);

    rng(seed,'twister');

    for b = 1:nBoot
        idx = randi(n,n,1);

        xx = x(idx);
        yy = y(idx);

        boot(b) = mean(xx)/max(mean(yy),eps);
    end

    ci = [quantileLocal(boot,alpha/2), ...
          quantileLocal(boot,1-alpha/2)];
end

function r = safeMeanRatio(numerator,denominator)
    numerator = numerator(:);
    denominator = denominator(:);

    r = mean(numerator)/max(mean(denominator),eps);
end

function win = buildCircularWindowsSparse(coordX,coordY,radius)
    coordX = coordX(:);
    coordY = coordY(:);
    K = numel(coordX);

    I = cell(K,1);
    J = cell(K,1);

    radius2 = radius^2+1e-14;

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
        error('At least one fixed-radius window is empty.');
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

function idx = nearestIndex(v,target)
    [~,idx] = min(abs(v-target));
end

function plotUHBFMarkers(ax,betaRef,colU,colH,colB,colF)
    markerSize = 6.5;

    plot(ax,0,0,'o','MarkerSize',markerSize, ...
        'MarkerFaceColor',colU,'MarkerEdgeColor','w');
    text(ax,0,0,' U','Color',[0.05 0.05 0.05], ...
        'FontWeight','bold','FontSize',6.0, ...
        'VerticalAlignment','bottom');

    plot(ax,1,0,'o','MarkerSize',markerSize, ...
        'MarkerFaceColor',colB,'MarkerEdgeColor','w');
    text(ax,1,0,' B','Color',[0.05 0.05 0.05], ...
        'FontWeight','bold','FontSize',6.0, ...
        'VerticalAlignment','bottom');

    plot(ax,0,betaRef,'o','MarkerSize',markerSize, ...
        'MarkerFaceColor',colH,'MarkerEdgeColor','w');
    text(ax,0,betaRef,' H','Color',[0.05 0.05 0.05], ...
        'FontWeight','bold','FontSize',6.0, ...
        'VerticalAlignment','bottom');

    plot(ax,1,betaRef,'o','MarkerSize',markerSize, ...
        'MarkerFaceColor',colF,'MarkerEdgeColor','w');
    text(ax,1,betaRef,' F','Color',[0.05 0.05 0.05], ...
        'FontWeight','bold','FontSize',6.0, ...
        'VerticalAlignment','bottom');
end

function cmap = natureBlueYellow(n)
    anchors = [0.07 0.18 0.40
               0.08 0.43 0.59
               0.28 0.65 0.58
               0.72 0.77 0.38
               0.96 0.80 0.26];

    cmap = interpCmap(anchors,n);
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

function h = plotCellHeatmap(ax,xCenters,yCenters,C)
    % Draw one constant-color cell per sampled condition using the true
    % (possibly nonuniform) coordinate spacing. This keeps heat-map cells,
    % contours, and U/H/B/F markers registered to the same axes.
    xEdges = centersToEdges(xCenters);
    yEdges = centersToEdges(yCenters);

    Cpad = C([1:end end],[1:end end]);
    h = surface(ax,xEdges,yEdges,zeros(size(Cpad)),Cpad, ...
        'FaceColor','flat','EdgeColor','none');
    view(ax,2);
    set(ax,'YDir','normal');
    xlim(ax,[xEdges(1) xEdges(end)]);
    ylim(ax,[yEdges(1) yEdges(end)]);
end

function edges = centersToEdges(centers)
    centers = centers(:).';
    assert(~isempty(centers) && all(isfinite(centers)), ...
        'Heat-map coordinates must be finite and nonempty.');
    assert(numel(centers)==1 || all(diff(centers)>0), ...
        'Heat-map coordinates must be strictly increasing.');

    if numel(centers)==1
        width = max(1,abs(centers(1)));
        edges = centers(1)+[-0.5 0.5]*width;
        return;
    end

    mid = 0.5*(centers(1:end-1)+centers(2:end));
    edges = [centers(1)-0.5*(centers(2)-centers(1)), ...
             mid, ...
             centers(end)+0.5*(centers(end)-centers(end-1))];
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
    % Normalize finite values inside mask to [0,1] and preserve NaN outside.
    valid = mask & isfinite(Y);
    vals = Y(valid);

    if isempty(vals)
        error('normalizeInside:EmptyMask', ...
            'No finite values were found inside the supplied mask.');
    end

    ymin = min(vals);
    ymax = max(vals);

    if ymax-ymin <= eps(max(abs([ymin,ymax,1])))
        Y(valid) = 0;
    else
        Y(valid) = (Y(valid)-ymin)/(ymax-ymin);
    end

    Y(~mask) = NaN;
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
