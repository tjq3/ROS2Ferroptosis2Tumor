function CTMC_SpatialArrangementConditionScan()
% =========================================================================
% CTMC_SpatialArrangementConditionScan
%
% Q2 targeted condition map after the spatial-arrangement phenomenon and
% neighbourhood-access mechanism have both been established at the anchor
% point.
%
% PRIMARY QUESTION
% ----------------
% At fixed mean ferroptotic commitment pressure, when does the spatial
% organization of an otherwise identical susceptibility distribution matter?
%
% PARAMETERIZATION
% ----------------
% x-axis:
%   mean commitment pressure
%       chiBar = mean(lambda_F/r)
%
% y-axis:
%   relative susceptibility contrast
%       delta_F = (chi_H - chi_L)/(chi_H + chi_L),  0 <= delta_F < 1
%
% For the exact initial low/high cell counts used in the Q2 pilot, chi_L and
% chi_H are solved so that the INITIAL composition-weighted mean is EXACTLY
% chiBar:
%
%   q      = (1 + delta_F)/(1 - delta_F)
%   chi_H  = q*chi_L
%   chiBar = w_L*chi_L + w_H*chi_H
%
% Thus every Random-vs-Clustered comparison has:
%   - identical low/high counts
%   - identical chi_L and chi_H
%   - identical mean chiBar
%   - identical susceptibility distribution
%   - identical radial-shell composition
% and differs ONLY in spatial arrangement.
%
% FULL-MAP COMPARISON
% -------------------
%   Random versus Clustered
%
% PRIMARY MAP
% -----------
%   Delta_AUTC_CR(%) =
%       100 * [ mean(AUTC_C) - mean(AUTC_R) ] / mean(AUTC_R)
%
% SECONDARY MAPS
% --------------
%   Delta P_ext(168 h) = Pext_C - Pext_R
%
%   Delta access_CR(%) =
%       100 * [ mean(A_LE_C) - mean(A_LE_R) ] / mean(A_LE_R)
%
%   Delta birth_CR(%) =
%       100 * [ mean(B_L_C) - mean(B_L_R) ] / mean(B_L_R)
%
% A_LE is the EXACT event-time integral of the number of low-susceptibility
% viable -> empty neighbour edges, normalized by the initial low-lineage
% count. It is computed by event-log replay using only local edge updates.
%
% REQUIRED FILES
% --------------
% CTMC_EventDrivenSimulation.m
%
% Macro_Q2_SpatialArrangement_PILOT_OUT/
%   Macro_Q2_SpatialArrangement_PILOT.mat
%
% =========================================================================

clc;
close all;

%% 0. FINAL SETTINGS
% FINAL keeps the PILOT grid and spatial manipulation definition fixed and
% increases only Monte-Carlo replication. All FINAL realizations use new
% independent seed families.
RUN_MODE = 'final';
RESUME_FROM_CHECKPOINT = true;

USE_PARALLEL = true;
AUTO_START_POOL = true;

MEAN_CHI_BASE  = [0.40 0.70 1.00 1.60 1.90 2.20];
DELTA_BASE     = [0.00 0.20 0.40 0.60 0.80];

N_REP = 100;
N_BOOT = 10000;

SEED_PATTERN_BASE = 3010000;
SEED_SIM_BASE     = 3110000;
SEED_BOOT_BASE    = 3210000;

% Descriptive minimum effect magnitude used only to annotate robust cells.
% Statistical confirmation still requires the paired 95% CI to exclude zero.
ROBUST_EFFECT_MIN_PCT = 2.0;

%% 1. FIND SOURCE FILES
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir = pwd; end

coreFile = fullfile(scriptDir,'CTMC_EventDrivenSimulation.m');
assert(isfile(coreFile), ...
    'Missing CTMC_EventDrivenSimulation.m in the script directory.');

finalCandidates = { ...
    fullfile(scriptDir,'Macro_Q2_SpatialArrangement_FINAL_OUT', ...
        'Macro_Q2_SpatialArrangement_FINAL.mat'), ...
    fullfile(scriptDir,'Macro_Q2_SpatialArrangement_FINAL.mat')};

finalFile = '';
for k=1:numel(finalCandidates)
    if isfile(finalCandidates{k})
        finalFile = finalCandidates{k};
        break;
    end
end
assert(~isempty(finalFile), ...
    ['Cannot find Macro_Q2_SpatialArrangement_FINAL.mat. Run ', ...
     'CTMC_SpatialArrangementSimulation.m first.']);

S = load(finalFile,'Q2');
assert(isfield(S,'Q2'),'The FINAL spatial MAT file must contain Q2.');
Q2P = S.Q2;

%% 2. FREEZE PREVIOUS PARAMETERS
L        = double(Q2P.settings.L);
R0       = double(Q2P.settings.R0);
TEND_H   = double(Q2P.settings.Tend_h);
REC_DT_H = double(Q2P.settings.record_dt_h);
r_h      = double(Q2P.settings.r_h);
mu_h     = double(Q2P.settings.mu_h);

assert(isfield(Q2P.settings,'pHigh'), ...
    'Q2 pilot settings are missing pHigh.');
assert(isfield(Q2P.settings,'nSwapAttempts'), ...
    'Q2 pilot settings are missing nSwapAttempts.');

P_HIGH = double(Q2P.settings.pHigh);
N_SWAP_ATTEMPTS = double(Q2P.settings.nSwapAttempts);

anchorChiLow  = double(Q2P.hazard.chiLow);
anchorChiHigh = double(Q2P.hazard.chiHigh);
anchorChiBar  = double(Q2P.hazard.chiBar);
anchorDelta   = (anchorChiHigh-anchorChiLow) / ...
                (anchorChiHigh+anchorChiLow);

chiMu = mu_h/r_h;

meanChiGrid = sort(unique([MEAN_CHI_BASE(:).' anchorChiBar]));
deltaGrid   = sort(unique([DELTA_BASE(:).' anchorDelta]));

nMean  = numel(meanChiGrid);
nDelta = numel(deltaGrid);

[~,iAnchorMean] = min(abs(meanChiGrid-anchorChiBar));
[~,iAnchorDelta] = min(abs(deltaGrid-anchorDelta));

assert(abs(meanChiGrid(iAnchorMean)-anchorChiBar)<1e-10);
assert(abs(deltaGrid(iAnchorDelta)-anchorDelta)<1e-10);

%% 3. OUTPUT DIRECTORY / CHECKPOINT
outDir = fullfile(scriptDir,'Macro_Q2_TargetedConditionMap_FINAL_OUT');
if ~exist(outDir,'dir'), mkdir(outDir); end

checkpointFile = fullfile(outDir, ...
    'Macro_Q2_TargetedConditionMap_FINAL_checkpoint.mat');

outMAT = fullfile(outDir,'Macro_Q2_TargetedConditionMap_FINAL.mat');
outCSV = fullfile(outDir,'Macro_Q2_TargetedConditionMap_FINAL.csv');
outPNG = fullfile(outDir,'Macro_Q2_TargetedConditionMap_FINAL.png');
outPDF = fullfile(outDir,'Macro_Q2_TargetedConditionMap_FINAL.pdf');

%% 4. INITIAL GEOMETRY / EXACT LOW-HIGH COUNTS
G = initialGeometryAndShellCounts(L,R0,P_HIGH);

N0    = G.N0;
nLow  = G.nLow;
nHigh = G.nHigh;

fprintf('\n============================================================\n');
fprintf('MACRO Q2 TARGETED CONDITION MAP: %s\n',upper(RUN_MODE));
fprintf('Solver: exact event-driven direct-method SSA\n');
fprintf('Frozen r=%.8f h^-1; mu=%.8f h^-1; chi_mu=%.4f\n', ...
    r_h,mu_h,chiMu);
fprintf('Frozen L=%d; R0=%g; Tend=%.1f h\n',L,R0,TEND_H);
fprintf('Initial low/high counts = %d/%d of N0=%d\n',nLow,nHigh,N0);
fprintf('Full-map comparison: Random vs Clustered\n');
fprintf('Nrep per grid point = %d\n',N_REP);
fprintf('Mean chi grid :'); fprintf(' %.4g',meanChiGrid); fprintf('\n');
fprintf('delta_F grid  :'); fprintf(' %.4g',deltaGrid); fprintf('\n');
fprintf('Exact Q2 anchor: chiBar=%.6f, delta_F=%.6f\n', ...
    anchorChiBar,anchorDelta);
fprintf('============================================================\n');

%% 5. PREBUILD MATCHED RANDOM / CLUSTERED PATTERN LIBRARY
fprintf('\n[1/4] Building matched pattern library once per replicate...\n');

patternLib = cell(N_REP,1);

poolReady = false;
if USE_PARALLEL && license('test','Distrib_Computing_Toolbox')
    p = gcp('nocreate');
    if isempty(p) && AUTO_START_POOL
        try
            p = parpool('threads');
        catch ME
            warning('Could not start thread pool (%s). Running serially.',ME.message);
            p = [];
        end
    end
    poolReady = ~isempty(p);
end

if poolReady
    parfor rep=1:N_REP
        patternLib{rep} = buildPatternLabels( ...
            L,R0,P_HIGH,SEED_PATTERN_BASE+rep,N_SWAP_ATTEMPTS);
    end
else
    for rep=1:N_REP
        patternLib{rep} = buildPatternLabels( ...
            L,R0,P_HIGH,SEED_PATTERN_BASE+rep,N_SWAP_ATTEMPTS);
    end
end

moranR = nan(N_REP,1);
moranC = nan(N_REP,1);
sameR  = nan(N_REP,1);
sameC  = nan(N_REP,1);

for rep=1:N_REP
    P = patternLib{rep};
    assert(P.N0==N0 && P.nLow==nLow && P.nHigh==nHigh);
    moranR(rep)=P.metrics.random.moranI;
    moranC(rep)=P.metrics.clustered.moranI;
    sameR(rep)=P.metrics.random.sameEdgeFraction;
    sameC(rep)=P.metrics.clustered.sameEdgeFraction;
end

fprintf('  Random    Moran I = %.3f +/- %.3f\n',mean(moranR),std(moranR));
fprintf('  Clustered Moran I = %.3f +/- %.3f\n',mean(moranC),std(moranC));
fprintf('  Random/Clustered same-edge fraction = %.3f / %.3f\n', ...
    mean(sameR),mean(sameC));

assert(mean(moranC)>mean(moranR)+0.25, ...
    'Pattern library failed to produce clear clustering.');

%% 6. RESULT ARRAYS / RESUME
R = initializeResultArrays(nMean,nDelta,N_REP);
doneMask = false(nMean,nDelta);

signature = struct();
signature.version = 'Q2_TARGETED_CONDITION_MAP_FINAL_v1';
signature.runMode = RUN_MODE;
signature.meanChiGrid = meanChiGrid;
signature.deltaGrid = deltaGrid;
signature.N_REP = N_REP;
signature.L = L;
signature.R0 = R0;
signature.TEND_H = TEND_H;
signature.r_h = r_h;
signature.mu_h = mu_h;
signature.P_HIGH = P_HIGH;
signature.N_SWAP_ATTEMPTS = N_SWAP_ATTEMPTS;

if RESUME_FROM_CHECKPOINT && isfile(checkpointFile)
    C = load(checkpointFile,'CHK');
    if isfield(C,'CHK') && checkpointCompatible(C.CHK.signature,signature)
        R = C.CHK.R;
        doneMask = C.CHK.doneMask;
        fprintf('\nCheckpoint loaded: %d/%d grid states already complete.\n', ...
            nnz(doneMask),numel(doneMask));
    else
        warning('Checkpoint exists but is incompatible; starting fresh.');
    end
end

%% 7. RUN TARGETED MAP
fprintf('\n[2/4] Running targeted Random-vs-Clustered condition map...\n');

stateCounter = 0;
nStatesTotal = nMean*nDelta;

for id=1:nDelta
    deltaF = deltaGrid(id);

    for im=1:nMean
        stateCounter = stateCounter+1;

        if doneMask(im,id)
            fprintf('  [%2d/%2d] chiBar=%.4f, delta=%.4f : checkpointed, skip\n', ...
                stateCounter,nStatesTotal,meanChiGrid(im),deltaF);
            continue;
        end

        chiBarTarget = meanChiGrid(im);

        [chiLow,chiHigh,meanExact] = exactLowHighFromMeanContrast( ...
            chiBarTarget,deltaF,nLow,nHigh);

        R.chiLow(im,id)  = chiLow;
        R.chiHigh(im,id) = chiHigh;
        R.meanExact(im,id) = meanExact;

        fprintf('  [%2d/%2d] chiBar=%.4f, delta=%.4f -> chiL=%.4f, chiH=%.4f\n', ...
            stateCounter,nStatesTotal,chiBarTarget,deltaF,chiLow,chiHigh);

        tmp = cell(N_REP,1);
        stateSeedOffset = (id-1)*nMean + im;

        if poolReady
            parfor rep=1:N_REP
                tmp{rep} = runConditionPair( ...
                    patternLib{rep},L,R0,TEND_H,REC_DT_H,r_h,mu_h, ...
                    chiLow,chiHigh,deltaF, ...
                    SEED_SIM_BASE + stateSeedOffset*10000 + rep);
            end
        else
            for rep=1:N_REP
                tmp{rep} = runConditionPair( ...
                    patternLib{rep},L,R0,TEND_H,REC_DT_H,r_h,mu_h, ...
                    chiLow,chiHigh,deltaF, ...
                    SEED_SIM_BASE + stateSeedOffset*10000 + rep);
            end
        end

        for rep=1:N_REP
            Z = tmp{rep};

            R.AUTC_R(im,id,rep) = Z.AUTC_R;
            R.AUTC_C(im,id,rep) = Z.AUTC_C;
            R.finalNV_R(im,id,rep) = Z.finalNV_R;
            R.finalNV_C(im,id,rep) = Z.finalNV_C;
            R.ext_R(im,id,rep) = Z.ext_R;
            R.ext_C(im,id,rep) = Z.ext_C;
            R.access_R(im,id,rep) = Z.access_R;
            R.access_C(im,id,rep) = Z.access_C;
            R.birth_R(im,id,rep) = Z.birth_R;
            R.birth_C(im,id,rep) = Z.birth_C;
            R.highClearAccessible_R(im,id,rep) = Z.highClearAccessible_R;
            R.highClearAccessible_C(im,id,rep) = Z.highClearAccessible_C;
            R.replayFinalNVError_R(im,id,rep) = Z.replayFinalNVError_R;
            R.replayFinalNVError_C(im,id,rep) = Z.replayFinalNVError_C;
        end

        tmpR = reshape(R.replayFinalNVError_R(im,id,:),[],1);
        tmpC = reshape(R.replayFinalNVError_C(im,id,:),[],1);
        vv = [tmpR;tmpC];
        vv = vv(isfinite(vv));
        if isempty(vv)
            maxReplayErr = 0;
        else
            maxReplayErr = max(vv);
        end

        if deltaF>0
            assert(maxReplayErr==0, ...
                'Event-log replay failed final NV identity.');
        end

        doneMask(im,id) = true;

        CHK = struct();
        CHK.signature = signature;
        CHK.R = R;
        CHK.doneMask = doneMask;
        save(checkpointFile,'CHK','-v7.3');
    end
end

%% 8. COMPUTE EFFECT MAPS + BOOTSTRAP CIs
fprintf('\n[3/4] Computing paired effect maps and bootstrap uncertainty...\n');

M = initializeMapSummary(nMean,nDelta);

for id=1:nDelta
    for im=1:nMean

        bootSeed = SEED_BOOT_BASE + id*1000 + im;

        aR = reshape(R.AUTC_R(im,id,:),[],1);
        aC = reshape(R.AUTC_C(im,id,:),[],1);
        eR = reshape(R.ext_R(im,id,:),[],1);
        eC = reshape(R.ext_C(im,id,:),[],1);

        M.meanAUTC_R(im,id) = mean(aR);
        M.meanAUTC_C(im,id) = mean(aC);
        M.Pext_R(im,id) = mean(eR);
        M.Pext_C(im,id) = mean(eC);

        if abs(deltaGrid(id))<1e-14
            M.dAUTC_pct(im,id) = 0;
            M.dAUTC_lo(im,id) = 0;
            M.dAUTC_hi(im,id) = 0;
            M.dPext(im,id) = 0;
            M.dPext_lo(im,id) = 0;
            M.dPext_hi(im,id) = 0;

            M.dAccess_pct(im,id) = NaN;
            M.dAccess_lo(im,id) = NaN;
            M.dAccess_hi(im,id) = NaN;
            M.dBirth_pct(im,id) = NaN;
            M.dBirth_lo(im,id) = NaN;
            M.dBirth_hi(im,id) = NaN;
            M.dHighClearAccess(im,id) = NaN;
            M.dHighClearAccess_lo(im,id) = NaN;
            M.dHighClearAccess_hi(im,id) = NaN;
            M.robustAUTC(im,id) = 0;
            continue;
        end

        A = pairedRelativeRatioEffect(aC,aR,N_BOOT,bootSeed+1);
        M.dAUTC_pct(im,id) = A.point;
        M.dAUTC_lo(im,id) = A.ci(1);
        M.dAUTC_hi(im,id) = A.ci(2);

        P = pairedAbsoluteEffect(eC,eR,N_BOOT,bootSeed+2);
        M.dPext(im,id) = P.point;
        M.dPext_lo(im,id) = P.ci(1);
        M.dPext_hi(im,id) = P.ci(2);

        xR = reshape(R.access_R(im,id,:),[],1);
        xC = reshape(R.access_C(im,id,:),[],1);
        Xeff = pairedRelativeRatioEffect(xC,xR,N_BOOT,bootSeed+3);
        M.dAccess_pct(im,id) = Xeff.point;
        M.dAccess_lo(im,id) = Xeff.ci(1);
        M.dAccess_hi(im,id) = Xeff.ci(2);

        bR = reshape(R.birth_R(im,id,:),[],1);
        bC = reshape(R.birth_C(im,id,:),[],1);
        B = pairedRelativeRatioEffect(bC,bR,N_BOOT,bootSeed+4);
        M.dBirth_pct(im,id) = B.point;
        M.dBirth_lo(im,id) = B.ci(1);
        M.dBirth_hi(im,id) = B.ci(2);

        hR = reshape(R.highClearAccessible_R(im,id,:),[],1);
        hC = reshape(R.highClearAccessible_C(im,id,:),[],1);
        H = pairedAbsoluteEffect(hC,hR,N_BOOT,bootSeed+5);
        M.dHighClearAccess(im,id) = H.point;
        M.dHighClearAccess_lo(im,id) = H.ci(1);
        M.dHighClearAccess_hi(im,id) = H.ci(2);

        if M.dAUTC_hi(im,id)<0 && abs(M.dAUTC_pct(im,id))>=ROBUST_EFFECT_MIN_PCT
            M.robustAUTC(im,id) = -1;
        elseif M.dAUTC_lo(im,id)>0 && abs(M.dAUTC_pct(im,id))>=ROBUST_EFFECT_MIN_PCT
            M.robustAUTC(im,id) = +1;
        else
            M.robustAUTC(im,id) = 0;
        end
    end
end

%% 9. CROSS-GRID MECHANISM ALIGNMENT
maskHet = repmat(deltaGrid>0,nMean,1);

xAccess = M.dAccess_pct(maskHet);
xBirth  = M.dBirth_pct(maskHet);
yAUTC   = M.dAUTC_pct(maskHet);

validA = isfinite(xAccess) & isfinite(yAUTC);
validB = isfinite(xBirth) & isfinite(yAUTC);

[rAccessOutcome,slopeAccessOutcome,~] = ...
    simpleLinearAssociation(xAccess(validA),yAUTC(validA));
[rBirthOutcome,slopeBirthOutcome,~] = ...
    simpleLinearAssociation(xBirth(validB),yAUTC(validB));

%% 10. CONSOLE SUMMARY
fprintf('\n================ Q2 TARGETED MAP SUMMARY =================\n');

nNeg = nnz(M.robustAUTC==-1);
nPos = nnz(M.robustAUTC==+1);
nInc = nnz(maskHet & M.robustAUTC==0);

fprintf('Heterogeneous grid cells         : %d\n',nnz(maskHet));
fprintf('Robust negative clustering cells : %d\n',nNeg);
fprintf('Robust positive clustering cells : %d\n',nPos);
fprintf('Inconclusive / <%.1f%% cells      : %d\n',ROBUST_EFFECT_MIN_PCT,nInc);

Z = M.dAUTC_pct;
Z(~maskHet)=NaN;

[minVal,linMin] = min(Z,[],'all','omitnan','linear');
[maxVal,linMax] = max(Z,[],'all','omitnan','linear');

[imMin,idMin] = ind2sub(size(Z),linMin);
[imMax,idMax] = ind2sub(size(Z),linMax);

fprintf('Strongest negative point effect: %.2f%% at chiBar=%.4f, delta=%.4f\n', ...
    minVal,meanChiGrid(imMin),deltaGrid(idMin));
fprintf('Strongest positive point effect: %.2f%% at chiBar=%.4f, delta=%.4f\n', ...
    maxVal,meanChiGrid(imMax),deltaGrid(idMax));

fprintf('\nExact previously successful Q2 anchor:\n');
fprintf('  chiBar=%.6f; delta=%.6f\n', ...
    meanChiGrid(iAnchorMean),deltaGrid(iAnchorDelta));
fprintf('  Delta AUTC(C-R) = %.2f%% [%.2f, %.2f]\n', ...
    M.dAUTC_pct(iAnchorMean,iAnchorDelta), ...
    M.dAUTC_lo(iAnchorMean,iAnchorDelta), ...
    M.dAUTC_hi(iAnchorMean,iAnchorDelta));
fprintf('  Delta access(C-R) = %.2f%% [%.2f, %.2f]\n', ...
    M.dAccess_pct(iAnchorMean,iAnchorDelta), ...
    M.dAccess_lo(iAnchorMean,iAnchorDelta), ...
    M.dAccess_hi(iAnchorMean,iAnchorDelta));
fprintf('  Delta low births(C-R) = %.2f%% [%.2f, %.2f]\n', ...
    M.dBirth_pct(iAnchorMean,iAnchorDelta), ...
    M.dBirth_lo(iAnchorMean,iAnchorDelta), ...
    M.dBirth_hi(iAnchorMean,iAnchorDelta));
fprintf('  Delta Pext(C-R) = %+.3f [%.3f, %.3f]\n', ...
    M.dPext(iAnchorMean,iAnchorDelta), ...
    M.dPext_lo(iAnchorMean,iAnchorDelta), ...
    M.dPext_hi(iAnchorMean,iAnchorDelta));

fprintf('\nCross-grid mechanism alignment:\n');
fprintf('  Delta access -> Delta AUTC: r=%+.3f; slope=%+.3f\n', ...
    rAccessOutcome,slopeAccessOutcome);
fprintf('  Delta births -> Delta AUTC: r=%+.3f; slope=%+.3f\n', ...
    rBirthOutcome,slopeBirthOutcome);
fprintf('==========================================================\n');

%% 11. BUILD STATE TABLE
nRows = nMean*nDelta;

MeanChi = nan(nRows,1);
DeltaF = nan(nRows,1);
ChiLow = nan(nRows,1);
ChiHigh = nan(nRows,1);
ExactMean = nan(nRows,1);

MeanAUTC_R = nan(nRows,1);
MeanAUTC_C = nan(nRows,1);

DeltaAUTC_pct = nan(nRows,1);
DeltaAUTC_CI_lo = nan(nRows,1);
DeltaAUTC_CI_hi = nan(nRows,1);

Pext_R = nan(nRows,1);
Pext_C = nan(nRows,1);
DeltaPext = nan(nRows,1);
DeltaPext_CI_lo = nan(nRows,1);
DeltaPext_CI_hi = nan(nRows,1);

DeltaAccess_pct = nan(nRows,1);
DeltaAccess_CI_lo = nan(nRows,1);
DeltaAccess_CI_hi = nan(nRows,1);

DeltaLowBirth_pct = nan(nRows,1);
DeltaLowBirth_CI_lo = nan(nRows,1);
DeltaLowBirth_CI_hi = nan(nRows,1);

DeltaHighClearAccessible = nan(nRows,1);
RobustAUTC = nan(nRows,1);

k=0;
for id=1:nDelta
    for im=1:nMean
        k=k+1;

        MeanChi(k)=meanChiGrid(im);
        DeltaF(k)=deltaGrid(id);
        ChiLow(k)=R.chiLow(im,id);
        ChiHigh(k)=R.chiHigh(im,id);
        ExactMean(k)=R.meanExact(im,id);

        MeanAUTC_R(k)=M.meanAUTC_R(im,id);
        MeanAUTC_C(k)=M.meanAUTC_C(im,id);

        DeltaAUTC_pct(k)=M.dAUTC_pct(im,id);
        DeltaAUTC_CI_lo(k)=M.dAUTC_lo(im,id);
        DeltaAUTC_CI_hi(k)=M.dAUTC_hi(im,id);

        Pext_R(k)=M.Pext_R(im,id);
        Pext_C(k)=M.Pext_C(im,id);
        DeltaPext(k)=M.dPext(im,id);
        DeltaPext_CI_lo(k)=M.dPext_lo(im,id);
        DeltaPext_CI_hi(k)=M.dPext_hi(im,id);

        DeltaAccess_pct(k)=M.dAccess_pct(im,id);
        DeltaAccess_CI_lo(k)=M.dAccess_lo(im,id);
        DeltaAccess_CI_hi(k)=M.dAccess_hi(im,id);

        DeltaLowBirth_pct(k)=M.dBirth_pct(im,id);
        DeltaLowBirth_CI_lo(k)=M.dBirth_lo(im,id);
        DeltaLowBirth_CI_hi(k)=M.dBirth_hi(im,id);

        DeltaHighClearAccessible(k)=M.dHighClearAccess(im,id);
        RobustAUTC(k)=M.robustAUTC(im,id);
    end
end

T = table( ...
    MeanChi,DeltaF,ChiLow,ChiHigh,ExactMean, ...
    MeanAUTC_R,MeanAUTC_C, ...
    DeltaAUTC_pct,DeltaAUTC_CI_lo,DeltaAUTC_CI_hi, ...
    Pext_R,Pext_C,DeltaPext,DeltaPext_CI_lo,DeltaPext_CI_hi, ...
    DeltaAccess_pct,DeltaAccess_CI_lo,DeltaAccess_CI_hi, ...
    DeltaLowBirth_pct,DeltaLowBirth_CI_lo,DeltaLowBirth_CI_hi, ...
    DeltaHighClearAccessible,RobustAUTC);

writetable(T,outCSV);

%% 12. SAVE FINAL MAT
MAP = struct();
MAP.version = 'Q2_TARGETED_CONDITION_MAP_FINAL_v1';
MAP.runMode = RUN_MODE;
MAP.sourceFinalSpatial = finalFile;

MAP.settings = struct( ...
    'L',L,'R0',R0,'Tend_h',TEND_H,'recordDt_h',REC_DT_H, ...
    'r_h',r_h,'mu_h',mu_h,'chiMu',chiMu, ...
    'pHigh',P_HIGH,'nRep',N_REP,'nBoot',N_BOOT, ...
    'nSwapAttempts',N_SWAP_ATTEMPTS, ...
    'robustEffectMin_pct',ROBUST_EFFECT_MIN_PCT);

MAP.grid.meanChi = meanChiGrid;
MAP.grid.deltaF = deltaGrid;
MAP.grid.anchorChiBar = anchorChiBar;
MAP.grid.anchorDeltaF = anchorDelta;
MAP.grid.anchorIndex = [iAnchorMean iAnchorDelta];

MAP.initial = struct( ...
    'N0',N0,'nLow',nLow,'nHigh',nHigh, ...
    'meanMoranRandom',mean(moranR), ...
    'meanMoranClustered',mean(moranC), ...
    'meanSameEdgeRandom',mean(sameR), ...
    'meanSameEdgeClustered',mean(sameC));

MAP.raw = R;
MAP.summary = M;
MAP.table = T;

MAP.mechanismAlignment = struct( ...
    'r_access_to_AUTC',rAccessOutcome, ...
    'slope_access_to_AUTC',slopeAccessOutcome, ...
    'r_birth_to_AUTC',rBirthOutcome, ...
    'slope_birth_to_AUTC',slopeBirthOutcome);

MAP.counts = struct( ...
    'nRobustNegative',nNeg, ...
    'nRobustPositive',nPos, ...
    'nInconclusive',nInc);

save(outMAT,'MAP','-v7.3');

%% 13. FIGURE
fprintf('\n[4/4] Drawing targeted Q2 condition-map figure...\n');

fig = figure('Color','w','Position',[60 40 1480 850]);
tl = tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');

cmapDiv = divergingMap(257);
cmapBlue = sequentialBlueMap(257);

% a
ax = nexttile(tl,1);
plotDiscreteMap(ax,meanChiGrid,deltaGrid,M.dAUTC_pct.',cmapDiv);
hold(ax,'on');
overlayRobustDots(ax,meanChiGrid,deltaGrid,M.robustAUTC);
plot(ax,anchorChiBar,anchorDelta,'p', ...
    'MarkerSize',11,'MarkerFaceColor','k','MarkerEdgeColor','w', ...
    'LineWidth',0.8);
xlabel(ax,'Mean commitment pressure, $\bar{\chi}_F$','Interpreter','latex');
ylabel(ax,'Susceptibility contrast, $\delta_F$','Interpreter','latex');
title(ax,'a  Spatial effect on integrated viable burden');
cb=colorbar(ax);
cb.Label.String='$\Delta AUTC_{NV}^{C-R}$ (\%)'; cb.Label.Interpreter='latex';
styleAxis(ax);

% b
ax = nexttile(tl,2);
plotDiscreteMap(ax,meanChiGrid,deltaGrid,M.dPext.',cmapDiv);
hold(ax,'on');
plot(ax,anchorChiBar,anchorDelta,'p', ...
    'MarkerSize',11,'MarkerFaceColor','k','MarkerEdgeColor','w', ...
    'LineWidth',0.8);
xlabel(ax,'Mean commitment pressure, $\bar{\chi}_F$','Interpreter','latex');
ylabel(ax,'Susceptibility contrast, $\delta_F$','Interpreter','latex');
title(ax,'b  Spatial effect on extinction probability');
cb=colorbar(ax);
cb.Label.String='$\Delta P_{ext}^{C-R}(168\,h)$'; cb.Label.Interpreter='latex';
styleAxis(ax);

% c
ax = nexttile(tl,3);
plotDiscreteMap(ax,meanChiGrid,deltaGrid,M.dAccess_pct.',cmapDiv);
hold(ax,'on');
plot(ax,anchorChiBar,anchorDelta,'p', ...
    'MarkerSize',11,'MarkerFaceColor','k','MarkerEdgeColor','w', ...
    'LineWidth',0.8);
xlabel(ax,'Mean commitment pressure, $\bar{\chi}_F$','Interpreter','latex');
ylabel(ax,'Susceptibility contrast, $\delta_F$','Interpreter','latex');
title(ax,'c  Spatial effect on low-lineage vacancy access');
cb=colorbar(ax);
cb.Label.String='$\Delta$ access$^{C-R}$ (\%)'; cb.Label.Interpreter='latex';
styleAxis(ax);

% d
ax = nexttile(tl,4);
plotDiscreteMap(ax,meanChiGrid,deltaGrid,M.dBirth_pct.',cmapDiv);
hold(ax,'on');
plot(ax,anchorChiBar,anchorDelta,'p', ...
    'MarkerSize',11,'MarkerFaceColor','k','MarkerEdgeColor','w', ...
    'LineWidth',0.8);
xlabel(ax,'Mean commitment pressure, $\bar{\chi}_F$','Interpreter','latex');
ylabel(ax,'Susceptibility contrast, $\delta_F$','Interpreter','latex');
title(ax,'d  Spatial effect on low-lineage birth output');
cb=colorbar(ax);
cb.Label.String='$\Delta$ low-lineage births$^{C-R}$ (\%)'; cb.Label.Interpreter='latex';
styleAxis(ax);

% e
ax = nexttile(tl,5);
hold(ax,'on');

deltasForScatter = repmat(deltaGrid,nMean,1);
mask = maskHet & isfinite(M.dAccess_pct) & isfinite(M.dAUTC_pct);

scatter(ax,M.dAccess_pct(mask),M.dAUTC_pct(mask), ...
    50,deltasForScatter(mask),'filled', ...
    'MarkerEdgeColor','w','LineWidth',0.5);

xa = M.dAccess_pct(mask);
ya = M.dAUTC_pct(mask);
if numel(xa)>=2
    pfit = polyfit(xa,ya,1);
    xx = linspace(min(xa),max(xa),100);
    plot(ax,xx,polyval(pfit,xx),'-', ...
        'Color',[0.2 0.2 0.2],'LineWidth',1.4);
end

xlabel(ax,'$\Delta$ low-lineage vacancy access$^{C-R}$ (\%)','Interpreter','latex');
ylabel(ax,'$\Delta AUTC_{NV}^{C-R}$ (\%)','Interpreter','latex');
title(ax,sprintf('e  Mechanism-outcome alignment  (r=%.2f)', ...
    rAccessOutcome));
cb=colorbar(ax);
colormap(ax,cmapBlue);
cb.Label.String='$\delta_F$'; cb.Label.Interpreter='latex';
styleAxis(ax);

% f
ax = nexttile(tl,6);
hold(ax,'on');

id = iAnchorDelta;
x = meanChiGrid;
y = M.dAUTC_pct(:,id).';
lo = M.dAUTC_lo(:,id).';
hi = M.dAUTC_hi(:,id).';

fill(ax,[x fliplr(x)],[lo fliplr(hi)], ...
    hex2rgb('#DCEBE7'), ...
    'EdgeColor','none','FaceAlpha',0.75);

plot(ax,x,y,'-o', ...
    'Color',hex2rgb('#2D7F78'), ...
    'MarkerFaceColor','w','MarkerEdgeColor',hex2rgb('#2D7F78'), ...
    'LineWidth',1.6,'MarkerSize',5);

yline(ax,0,':','Color',[0.35 0.35 0.35],'LineWidth',1);

plot(ax,anchorChiBar,M.dAUTC_pct(iAnchorMean,iAnchorDelta),'p', ...
    'MarkerSize',11,'MarkerFaceColor','k','MarkerEdgeColor','w');

xlabel(ax,'Mean commitment pressure, $\bar{\chi}_F$','Interpreter','latex');
ylabel(ax,'$\Delta AUTC_{NV}^{C-R}$ (\%)','Interpreter','latex');
title(ax,sprintf('f  Slice at anchor contrast, delta_F = %.2f',anchorDelta));
styleAxis(ax);

sgtitle(tl, ...
    'Q2 final condition map: when spatial organization changes macro growth', ...
    'FontWeight','normal','FontSize',12);

try
    axs=findall(fig,'Type','axes');
    for k=1:numel(axs)
        axtoolbar(axs(k),{});
    end
catch
end

exportgraphics(fig,outPNG,'Resolution',500);
exportgraphics(fig,outPDF,'ContentType','vector');

fprintf('\nSaved Q2 targeted-map outputs:\n');
fprintf('  %s\n',outMAT);
fprintf('  %s\n',outCSV);
fprintf('  %s\n',outPNG);
fprintf('  %s\n',outPDF);
fprintf('  %s\n',checkpointFile);

fprintf('\nFINAL STATUS:\n');
fprintf(['This is the confirmatory condition map. Combine it with the independent ', ...
    'FINAL spatial-arrangement and neighbourhood-access mechanism results. ', ...
    'If the sign structure and mechanism co-localization remain consistent, ', ...
    'Q2 can be frozen for Fig.5C/SI assembly.\n']);

end

%% =========================================================================
% RUN ONE MATCHED RANDOM-vs-CLUSTERED PAIR
% =========================================================================
function Z = runConditionPair( ...
    P,L,R0,Tend,recordDt,r,mu,chiLow,chiHigh,deltaF,seed)

Z = struct();

if abs(deltaF)<1e-14
    chiBar = chiLow;
    map = chiBar*ones(L,L,L);

    cfg = baseCfg(L,R0,Tend,recordDt,r,mu,seed);
    cfg.lambda_mode = 'cell_mark';
    cfg.lambdaF_h = chiBar*r;
    cfg.lambdaF_initial_map_h = map*r;
    cfg.daughter_lambda_rule = 'inherit';
    cfg.lambdaF_pool_h = [];
    cfg.lambdaF_site_h = [];
    cfg.store_event_log = false;

    O = CTMC_EventDrivenSimulation(cfg);

    Z.AUTC_R = double(O.AUTC_NV);
    Z.AUTC_C = double(O.AUTC_NV);
    Z.finalNV_R = double(O.NV(end));
    Z.finalNV_C = double(O.NV(end));

    ex = double(O.NV(end)==0);
    Z.ext_R = ex;
    Z.ext_C = ex;

    Z.access_R = NaN;
    Z.access_C = NaN;
    Z.birth_R = NaN;
    Z.birth_C = NaN;
    Z.highClearAccessible_R = NaN;
    Z.highClearAccessible_C = NaN;
    Z.replayFinalNVError_R = NaN;
    Z.replayFinalNVError_C = NaN;
    return;
end

mapR = makeChiMap(P,L,chiLow,chiHigh,'random');
mapC = makeChiMap(P,L,chiLow,chiHigh,'clustered');

cfg = baseCfg(L,R0,Tend,recordDt,r,mu,seed);
cfg.lambda_mode = 'cell_mark';
cfg.lambdaF_h = ((P.nLow*chiLow + P.nHigh*chiHigh)/P.N0)*r;
cfg.daughter_lambda_rule = 'inherit';
cfg.lambdaF_pool_h = [];
cfg.lambdaF_site_h = [];
cfg.store_event_log = true;

cfg.lambdaF_initial_map_h = mapR*r;
OR = CTMC_EventDrivenSimulation(cfg);
DR = replayExactMechanism(OR,mapR,P.initialMask,chiLow,chiHigh);

cfg.lambdaF_initial_map_h = mapC*r;
OC = CTMC_EventDrivenSimulation(cfg);
DC = replayExactMechanism(OC,mapC,P.initialMask,chiLow,chiHigh);

Z.AUTC_R = double(OR.AUTC_NV);
Z.AUTC_C = double(OC.AUTC_NV);
Z.finalNV_R = double(OR.NV(end));
Z.finalNV_C = double(OC.NV(end));
Z.ext_R = double(OR.NV(end)==0);
Z.ext_C = double(OC.NV(end)==0);

Z.access_R = DR.integratedLowEmptyEdgeHours / P.nLow;
Z.access_C = DC.integratedLowEmptyEdgeHours / P.nLow;
Z.birth_R = DR.lowBirthTotal / P.nLow;
Z.birth_C = DC.lowBirthTotal / P.nLow;
Z.highClearAccessible_R = DR.highClearAccessibleFrac;
Z.highClearAccessible_C = DC.highClearAccessibleFrac;
Z.replayFinalNVError_R = DR.finalNVError;
Z.replayFinalNVError_C = DC.finalNVError;

end

%% =========================================================================
% BASE SSA CFG
% =========================================================================
function cfg = baseCfg(L,R0,Tend,recordDt,r,mu,seed)

cfg = struct();
cfg.L = L;
cfg.R0 = R0;
cfg.Tend_h = Tend;
cfg.record_dt_h = recordDt;
cfg.r_h = r;
cfg.mu_h = mu;
cfg.seed = seed;

cfg.snapshot_times_h = [];
cfg.return_final_state = false;
cfg.compute_cavity = false;
cfg.verbose = false;

end

%% =========================================================================
% EXACT EVENT-TIME MECHANISM REPLAY
% =========================================================================
function D = replayExactMechanism(O,initialChiMap,initialMask,chiLow,chiHigh)

L = O.cfg.L;
nSite = L^3;

state = zeros(nSite,1,'uint8');
initialIdx = find(initialMask(:));
state(initialIdx) = uint8(1);

mark = nan(nSite,1);
chiVec = double(initialChiMap(:));
mark(initialIdx) = chiVec(initialIdx);

tol = 1e-9*max(1,abs(chiHigh));

EL = O.event_log;
assert(~isempty(EL),'Event log is empty although store_event_log=true.');

eventTime   = double(EL.time_h);
eventType   = double(EL.eventType);
eventSite   = double(EL.site);
eventTarget = double(EL.target);

currentLE = totalLowEmptyEdges(state,mark,chiLow,tol,L);

integratedLE = 0;
tPrev = 0;

nLowBirth = 0;
nBirth = 0;
nFerro = 0;

nHighClear = 0;
nHighClearAccessible = 0;

for ev=1:height(EL)

    te = min(eventTime(ev),O.cfg.Tend_h);

    if te>tPrev
        integratedLE = integratedLE + currentLE*(te-tPrev);
        tPrev = te;
    end

    typ = eventType(ev);
    s   = eventSite(ev);
    tg  = eventTarget(ev);

    if typ==1
        oldLocal = localLowEmptyContribution(tg,state,mark,chiLow,tol,L);

        assert(state(s)==1,'Replay birth parent is not viable.');
        assert(state(tg)==0,'Replay birth target is not empty.');
        assert(isfinite(mark(s)),'Replay birth parent has no mark.');

        parentMark = mark(s);

        state(tg) = uint8(1);
        mark(tg) = parentMark;

        newLocal = localLowEmptyContribution(tg,state,mark,chiLow,tol,L);
        currentLE = currentLE + (newLocal-oldLocal);

        nBirth = nBirth+1;
        if abs(parentMark-chiLow)<=tol
            nLowBirth = nLowBirth+1;
        end

    elseif typ==2
        oldLocal = localLowEmptyContribution(s,state,mark,chiLow,tol,L);

        assert(state(s)==1,'Replay ferroptosis site is not viable.');
        assert(isfinite(mark(s)),'Replay ferroptosis site has no mark.');

        state(s)=uint8(2);

        newLocal = localLowEmptyContribution(s,state,mark,chiLow,tol,L);
        currentLE = currentLE + (newLocal-oldLocal);
        nFerro = nFerro+1;

    elseif typ==3
        oldLocal = localLowEmptyContribution(s,state,mark,chiLow,tol,L);

        assert(state(s)==2,'Replay clearance site is not F.');
        assert(isfinite(mark(s)),'Replay F site has no retained mark.');

        if abs(mark(s)-chiHigh)<=tol
            nHighClear=nHighClear+1;

            nbr=neighbourIndicesLocal(s,L);
            isLowNbr = state(nbr)==1 & ...
                isfinite(mark(nbr)) & ...
                abs(mark(nbr)-chiLow)<=tol;

            if any(isLowNbr)
                nHighClearAccessible=nHighClearAccessible+1;
            end
        end

        state(s)=uint8(0);
        mark(s)=NaN;

        newLocal = localLowEmptyContribution(s,state,mark,chiLow,tol,L);
        currentLE = currentLE + (newLocal-oldLocal);

    else
        error('Unknown event type in event log.');
    end
end

if tPrev<O.cfg.Tend_h
    integratedLE = integratedLE + currentLE*(O.cfg.Tend_h-tPrev);
end

finalNVReplay = nnz(state==1);
finalNVError = abs(finalNVReplay-double(O.NV(end)));

NVidentity = numel(initialIdx)+nBirth-nFerro;
assert(NVidentity==finalNVReplay, ...
    'Replay birth/ferro count identity failed.');

if nHighClear>0
    highClearAccessibleFrac = nHighClearAccessible/nHighClear;
else
    highClearAccessibleFrac = NaN;
end

D = struct();
D.integratedLowEmptyEdgeHours = integratedLE;
D.lowBirthTotal = nLowBirth;
D.highClearAccessibleFrac = highClearAccessibleFrac;
D.finalNVError = finalNVError;

end

%% =========================================================================
% LOW-EMPTY EDGE ACCOUNTING
% =========================================================================
function total = totalLowEmptyEdges(state,mark,chiLow,tol,L)

lowIdx = find(state==1 & isfinite(mark) & abs(mark-chiLow)<=tol);

total=0;
for k=1:numel(lowIdx)
    nbr=neighbourIndicesLocal(lowIdx(k),L);
    total=total+nnz(state(nbr)==0);
end

end

function c = localLowEmptyContribution(idx,state,mark,chiLow,tol,L)

nbr = neighbourIndicesLocal(idx,L);
c=0;

isLowHere = state(idx)==1 && isfinite(mark(idx)) && ...
    abs(mark(idx)-chiLow)<=tol;
isEmptyHere = state(idx)==0;

for q=1:numel(nbr)
    j=nbr(q);

    isLowNbr = state(j)==1 && isfinite(mark(j)) && ...
        abs(mark(j)-chiLow)<=tol;
    isEmptyNbr = state(j)==0;

    if (isLowHere && isEmptyNbr) || (isEmptyHere && isLowNbr)
        c=c+1;
    end
end

end

function nbr=neighbourIndicesLocal(idx,L)

[i,j,k]=ind2sub([L L L],idx);
nbr=zeros(6,1);
n=0;

if i>1, n=n+1; nbr(n)=sub2ind([L L L],i-1,j,k); end
if i<L, n=n+1; nbr(n)=sub2ind([L L L],i+1,j,k); end
if j>1, n=n+1; nbr(n)=sub2ind([L L L],i,j-1,k); end
if j<L, n=n+1; nbr(n)=sub2ind([L L L],i,j+1,k); end
if k>1, n=n+1; nbr(n)=sub2ind([L L L],i,j,k-1); end
if k<L, n=n+1; nbr(n)=sub2ind([L L L],i,j,k+1); end

nbr=nbr(1:n);

end

%% =========================================================================
% EXACT LOW/HIGH VALUES AT FIXED COMPOSITION-WEIGHTED MEAN
% =========================================================================
function [chiLow,chiHigh,meanExact] = ...
    exactLowHighFromMeanContrast(chiBar,deltaF,nLow,nHigh)

assert(deltaF>=0 && deltaF<1, ...
    'delta_F must satisfy 0 <= delta_F < 1.');

N=nLow+nHigh;
wL=nLow/N;
wH=nHigh/N;

if abs(deltaF)<1e-14
    chiLow=chiBar;
    chiHigh=chiBar;
    meanExact=chiBar;
    return;
end

q=(1+deltaF)/(1-deltaF);

chiLow=chiBar/(wL+wH*q);
chiHigh=q*chiLow;

meanExact=wL*chiLow+wH*chiHigh;

assert(chiLow>=0 && chiHigh>=chiLow);
assert(abs(meanExact-chiBar)<=1e-12*max(1,chiBar), ...
    'Exact mean-locking failed.');

end

%% =========================================================================
% INITIAL GEOMETRY / PATTERN LIBRARY
% =========================================================================
function G = initialGeometryAndShellCounts(L,R0,pHigh)

ctr=[(L+1)/2 (L+1)/2 (L+1)/2];
[xg,yg,zg]=ndgrid(1:L,1:L,1:L);

r2=(xg-ctr(1)).^2+(yg-ctr(2)).^2+(zg-ctr(3)).^2;
mask=r2<=R0^2;
idx=find(mask);

N=numel(idx);
nHigh=round(pHigh*N);
nLow=N-nHigh;

G=struct();
G.initialMask=mask;
G.initialIdx=idx;
G.N0=N;
G.nHigh=nHigh;
G.nLow=nLow;

end

function P = buildPatternLabels(L,R0,pHigh,seed,nSwapAttempts)

rng(seed,'twister');

ctr=[(L+1)/2 (L+1)/2 (L+1)/2];
[xg,yg,zg]=ndgrid(1:L,1:L,1:L);

r2=(xg-ctr(1)).^2+(yg-ctr(2)).^2+(zg-ctr(3)).^2;
mask=r2<=R0^2;
idx=find(mask);

N=numel(idx);

[ii,jj,kk]=ind2sub([L L L],idx);
rad=sqrt((ii-ctr(1)).^2+(jj-ctr(2)).^2+(kk-ctr(3)).^2);
shell=floor(rad+1e-12);
uShell=unique(shell(:)).';

nHighTotal=round(pHigh*N);

nByShell=zeros(size(uShell));
rawHigh=zeros(size(uShell));

for s=1:numel(uShell)
    nByShell(s)=nnz(shell==uShell(s));
    rawHigh(s)=pHigh*nByShell(s);
end

nHighShell=floor(rawHigh);
remain=nHighTotal-sum(nHighShell);

[~,ordFrac]=sort(rawHigh-floor(rawHigh),'descend');
for k=1:remain
    nHighShell(ordFrac(k))=nHighShell(ordFrac(k))+1;
end

labelRandom=-ones(N,1);

for s=1:numel(uShell)
    loc=find(shell==uShell(s));
    perm=loc(randperm(numel(loc)));
    nh=nHighShell(s);

    if nh>0
        labelRandom(perm(1:nh))=+1;
    end
end

assert(nnz(labelRandom==1)==nHighTotal);

[e1,e2,incEdges]=initialGraph(idx,L);

labelClustered=optimizeArrangement( ...
    labelRandom,shell,e1,e2,incEdges,+1,nSwapAttempts);

P=struct();
P.N0=N;
P.nHigh=nHighTotal;
P.nLow=N-nHighTotal;

P.initialMask=mask;
P.initialIdx=idx;

P.labelRandom=labelRandom;
P.labelClustered=labelClustered;

P.metrics.random=patternMetrics(labelRandom,e1,e2);
P.metrics.clustered=patternMetrics(labelClustered,e1,e2);

end

function map = makeChiMap(P,L,chiLow,chiHigh,whichPattern)

chiBar=(P.nLow*chiLow+P.nHigh*chiHigh)/P.N0;
map=chiBar*ones(L,L,L);

switch lower(whichPattern)
    case 'random'
        label=P.labelRandom;
    case 'clustered'
        label=P.labelClustered;
    otherwise
        error('Unknown pattern.');
end

idx=P.initialIdx;

map(idx(label<0))=chiLow;
map(idx(label>0))=chiHigh;

vals=map(idx);
expected=(P.nLow*chiLow+P.nHigh*chiHigh)/P.N0;
assert(abs(mean(vals)-expected)<1e-12*max(1,expected));

end

function label=optimizeArrangement(label0,shell,e1,e2,incEdges,direction,nAttempt)

label=label0(:);

uShell=unique(shell(:)).';
eligible=false(size(uShell));

for s=1:numel(uShell)
    loc=find(shell==uShell(s));
    eligible(s)=any(label(loc)>0) && any(label(loc)<0);
end

eligibleShell=uShell(eligible);
if isempty(eligibleShell), return; end

for it=1:nAttempt
    sh=eligibleShell(randi(numel(eligibleShell)));
    loc=find(shell==sh);

    pos=loc(label(loc)>0);
    neg=loc(label(loc)<0);

    if isempty(pos) || isempty(neg), continue; end

    a=pos(randi(numel(pos)));
    b=neg(randi(numel(neg)));

    affected=unique([incEdges{a};incEdges{b}]);
    oldScore=sum(label(e1(affected)).*label(e2(affected)));

    label(a)=-label(a);
    label(b)=-label(b);

    newScore=sum(label(e1(affected)).*label(e2(affected)));
    delta=newScore-oldScore;

    accept=(direction*delta>0);
    if delta==0
        accept=(rand<0.10);
    end

    if ~accept
        label(a)=-label(a);
        label(b)=-label(b);
    end
end

end

function [e1,e2,incEdges]=initialGraph(globalIdx,L)

N=numel(globalIdx);

localOfGlobal=zeros(L^3,1,'uint32');
localOfGlobal(globalIdx)=uint32(1:N);

e1=zeros(3*N,1,'uint32');
e2=zeros(3*N,1,'uint32');
ne=0;

for a=1:N
    g=globalIdx(a);
    [i,j,k]=ind2sub([L L L],g);

    if i<L
        b=localOfGlobal(sub2ind([L L L],i+1,j,k));
        if b>0
            ne=ne+1; e1(ne)=uint32(a); e2(ne)=b;
        end
    end

    if j<L
        b=localOfGlobal(sub2ind([L L L],i,j+1,k));
        if b>0
            ne=ne+1; e1(ne)=uint32(a); e2(ne)=b;
        end
    end

    if k<L
        b=localOfGlobal(sub2ind([L L L],i,j,k+1));
        if b>0
            ne=ne+1; e1(ne)=uint32(a); e2(ne)=b;
        end
    end
end

e1=double(e1(1:ne));
e2=double(e2(1:ne));

incEdges=cell(N,1);

for e=1:ne
    incEdges{e1(e)}(end+1,1)=e; %#ok<AGROW>
    incEdges{e2(e)}(end+1,1)=e; %#ok<AGROW>
end

end

function M=patternMetrics(label,e1,e2)

label=double(label(:));
N=numel(label);
E=numel(e1);

z=label-mean(label);
den=sum(z.^2);

if den<=eps || E==0
    I=NaN;
else
    I=(N/E)*sum(z(e1).*z(e2))/den;
end

M=struct();
M.moranI=I;
M.sameEdgeFraction=mean(label(e1)==label(e2));
M.nEdges=E;

end

%% =========================================================================
% ARRAYS
% =========================================================================
function R=initializeResultArrays(nMean,nDelta,nRep)

sz=[nMean nDelta nRep];

R=struct();

R.chiLow=nan(nMean,nDelta);
R.chiHigh=nan(nMean,nDelta);
R.meanExact=nan(nMean,nDelta);

R.AUTC_R=nan(sz);
R.AUTC_C=nan(sz);
R.finalNV_R=nan(sz);
R.finalNV_C=nan(sz);
R.ext_R=nan(sz);
R.ext_C=nan(sz);
R.access_R=nan(sz);
R.access_C=nan(sz);
R.birth_R=nan(sz);
R.birth_C=nan(sz);
R.highClearAccessible_R=nan(sz);
R.highClearAccessible_C=nan(sz);
R.replayFinalNVError_R=nan(sz);
R.replayFinalNVError_C=nan(sz);

end

function M=initializeMapSummary(nMean,nDelta)

M=struct();

fields={ ...
    'meanAUTC_R','meanAUTC_C', ...
    'Pext_R','Pext_C', ...
    'dAUTC_pct','dAUTC_lo','dAUTC_hi', ...
    'dPext','dPext_lo','dPext_hi', ...
    'dAccess_pct','dAccess_lo','dAccess_hi', ...
    'dBirth_pct','dBirth_lo','dBirth_hi', ...
    'dHighClearAccess','dHighClearAccess_lo','dHighClearAccess_hi'};

for k=1:numel(fields)
    M.(fields{k})=nan(nMean,nDelta);
end

M.robustAUTC=zeros(nMean,nDelta);

end

%% =========================================================================
% CHECKPOINT
% =========================================================================
function tf=checkpointCompatible(A,B)

tf=true;

must={'version','runMode','meanChiGrid','deltaGrid','N_REP','L','R0', ...
    'TEND_H','r_h','mu_h','P_HIGH','N_SWAP_ATTEMPTS'};

for k=1:numel(must)
    f=must{k};

    if ~isfield(A,f) || ~isfield(B,f)
        tf=false; return;
    end

    a=A.(f);
    b=B.(f);

    if isnumeric(a)
        if ~isequal(size(a),size(b))
            tf=false; return;
        end
        scale=max(1,max(abs(double(b(:)))));
        if any(abs(double(a(:))-double(b(:)))>1e-12*scale)
            tf=false; return;
        end
    else
        if ~isequal(a,b)
            tf=false; return;
        end
    end
end

end

%% =========================================================================
% BOOTSTRAP EFFECTS
% =========================================================================
function A=pairedRelativeRatioEffect(xC,xR,nBoot,seed)

xC=double(xC(:));
xR=double(xR(:));

mask=isfinite(xC)&isfinite(xR);
xC=xC(mask);
xR=xR(mask);

A=struct();
A.n=numel(xC);

if isempty(xC) || mean(xR)==0
    A.point=NaN;
    A.ci=[NaN NaN];
    return;
end

A.point=100*(mean(xC)-mean(xR))/mean(xR);

rng(seed,'twister');

n=numel(xC);
boot=nan(nBoot,1);

for b=1:nBoot
    ii=randi(n,n,1);

    mR=mean(xR(ii));
    mC=mean(xC(ii));

    if mR~=0
        boot(b)=100*(mC-mR)/mR;
    end
end

A.ci=percentileCI(boot);

end

function A=pairedAbsoluteEffect(xC,xR,nBoot,seed)

xC=double(xC(:));
xR=double(xR(:));

mask=isfinite(xC)&isfinite(xR);
d=xC(mask)-xR(mask);

A=struct();
A.n=numel(d);

if isempty(d)
    A.point=NaN;
    A.ci=[NaN NaN];
    return;
end

A.point=mean(d);

rng(seed,'twister');

n=numel(d);
boot=zeros(nBoot,1);

for b=1:nBoot
    ii=randi(n,n,1);
    boot(b)=mean(d(ii));
end

A.ci=percentileCI(boot);

end

function ci=percentileCI(x)

x=sort(x(isfinite(x)));

if isempty(x)
    ci=[NaN NaN];
    return;
end

n=numel(x);
ilo=max(1,round(0.025*n));
ihi=min(n,round(0.975*n));

ci=[x(ilo) x(ihi)];

end

function [r,slope,intercept]=simpleLinearAssociation(x,y)

x=double(x(:));
y=double(y(:));

mask=isfinite(x)&isfinite(y);
x=x(mask);
y=y(mask);

if numel(x)<2
    r=NaN; slope=NaN; intercept=NaN;
    return;
end

C=corrcoef(x,y);
r=C(1,2);

p=polyfit(x,y,1);
slope=p(1);
intercept=p(2);

end

%% =========================================================================
% PLOTTING
% =========================================================================
function plotDiscreteMap(ax,xCenters,yCenters,Z,cmap)
% Draw each simulated grid state as a discrete cell. No interpolation is
% performed; this avoids visually inventing continuous islands between Monte-
% Carlo design points.

xEdges = centersToEdges(xCenters);
yEdges = centersToEdges(yCenters);

% pcolor requires a (ny+1)-by-(nx+1) C matrix. The final row/column are
% padding; with flat shading, each rectangle uses the lower-left C value.
Cpad = nan(numel(yCenters)+1,numel(xCenters)+1);
Cpad(1:numel(yCenters),1:numel(xCenters)) = Z;

[Xe,Ye] = meshgrid(xEdges,yEdges);
h = pcolor(ax,Xe,Ye,Cpad);
set(h,'EdgeColor',[1 1 1]*0.88,'LineWidth',0.35);
shading(ax,'flat');

finiteZ=Z(isfinite(Z));
if isempty(finiteZ)
    mx=1;
else
    mx=max(abs(finiteZ));
    if mx<1e-12, mx=1; end
end

caxis(ax,[-mx mx]);
colormap(ax,cmap);
xlim(ax,[xEdges(1) xEdges(end)]);
ylim(ax,[yEdges(1) yEdges(end)]);

end

function e=centersToEdges(c)
c=double(c(:).');
if numel(c)==1
    e=[c-0.5 c+0.5];
    return;
end
mid=(c(1:end-1)+c(2:end))/2;
e=[c(1)-(mid(1)-c(1)), mid, c(end)+(c(end)-mid(end))];
end

function overlayRobustDots(ax,meanGrid,deltaGrid,robustMap)

for id=1:numel(deltaGrid)
    for im=1:numel(meanGrid)
        if robustMap(im,id)==-1
            plot(ax,meanGrid(im),deltaGrid(id),'o', ...
                'MarkerSize',4,'MarkerFaceColor',[0.1 0.1 0.1], ...
                'MarkerEdgeColor','w','LineWidth',0.4);
        elseif robustMap(im,id)==+1
            plot(ax,meanGrid(im),deltaGrid(id),'^', ...
                'MarkerSize',4.5,'MarkerFaceColor',[0.1 0.1 0.1], ...
                'MarkerEdgeColor','w','LineWidth',0.4);
        end
    end
end

end

function cmap=divergingMap(n)

if nargin<1, n=257; end

blue=hex2rgb('#2F7EA8');
mid=hex2rgb('#F7F5EF');
red=hex2rgb('#C85C4A');

n1=ceil(n/2);
n2=n-n1+1;

a=zeros(n1,3);
b=zeros(n2,3);

for k=1:3
    a(:,k)=linspace(blue(k),mid(k),n1);
    b(:,k)=linspace(mid(k),red(k),n2);
end

cmap=[a;b(2:end,:)];

end

function cmap=sequentialBlueMap(n)

if nargin<1, n=257; end

c1=hex2rgb('#EAF4EF');
c2=hex2rgb('#7FC8C6');
c3=hex2rgb('#1F6FA8');

n1=ceil(n/2);
n2=n-n1+1;

a=zeros(n1,3);
b=zeros(n2,3);

for k=1:3
    a(:,k)=linspace(c1(k),c2(k),n1);
    b(:,k)=linspace(c2(k),c3(k),n2);
end

cmap=[a;b(2:end,:)];

end

function styleAxis(ax)

ax.FontName='Arial';
ax.FontSize=8.5;
ax.LineWidth=0.8;
ax.TickDir='out';
ax.TickLength=[0.015 0.015];
ax.Box='off';
ax.Layer='top';
grid(ax,'off');

end

function rgb=hex2rgb(h)

h=char(h);

if h(1)=='#'
    h=h(2:end);
end

rgb=[hex2dec(h(1:2)) hex2dec(h(3:4)) hex2dec(h(5:6))]/255;

end
