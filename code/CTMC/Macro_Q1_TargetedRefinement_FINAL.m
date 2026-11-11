function Macro_Q1_TargetedRefinement_FINAL()
% =========================================================================
% Macro_Q1_TargetedRefinement_FINAL
%
% Targeted refinement of Macro Q1 after the first Q1 pilot.
%
% SCIENTIFIC AIM
% ----------------
% Resolve the sequence:
% sustained growth -> slow/near-arrest growth -> persistent regression
% -> stochastic extinction transition -> extinction-dominated outcome,
% and quantify how post-commitment clearance shifts these boundaries.
%
% FROZEN MACRO MODEL
% ------------------
% V --(r*q_i)--------> V + V
% V --(lambda_F)-----> F
% F --(mu)-----------> empty
%
% Solver: macro_sim_core_hazard_SSA_FINAL.m
%         event-driven direct-method SSA.
%
% KEY IMPROVEMENTS OVER Macro_Q1_RegimeMap_FINAL
% -----------------------------------------------
% 1) Dense chi_F scan only in the true transition region.
% 2) Additional chi_mu interpolation points without recalibrating mu.
% 3) Survivor-aware late-growth metric.
% 4) Late-growth masking in extinction-dominated regions.
% 5) Isotonic extinction curves.
% 6) Bootstrap 95% CIs for transition boundaries.
% 7) Adaptive Monte-Carlo allocation.
%
% IMPORTANT CROSS-SCALE INTERPRETATION
% ------------------------------------
% The microscopic setting was intentionally pro-ferroptotic. Therefore the
% 48 mesoscopic states are expected to project into the high-pressure
% sector. No rescaling is introduced to force them to span macro regimes.
%
% OUTPUT DIRECTORY
% ----------------
% Macro_Q1_TargetedRefinement_FINAL_OUT/
% =========================================================================

clc;
close all;

%% 0. RUN SETTINGS
RUN_MODE = 'final';  % 'smoke' | 'pilot' | 'final'

TEND_H = 168;
RECORD_DT_H = 1;       % output sampling only, not a numerical time step
LATE_WINDOW_H = [120 168];

CHI_F_GRID = [ ...
    0 ...
    0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 1.00 1.10 1.20 ...
    1.30 1.40 1.50 1.60 1.70 1.80 1.90 2.00 ...
    2.10 2.20 2.30 2.40 2.50 2.60];

CHI_F_HIGH_ANCHOR = 8;

TAU_F_BENCHMARK_H = [48 24 12 6 3];
CHI_MU_INTERP = [1 2 4 8];
TAU_F_REFERENCE_H = 12;

MCSE_TARGET = 0.02;
PEXT_WILSON_HALF_TARGET = 0.09;

N_BOOT_BOUNDARY = 1500;
PEXT_LEVELS = [0.10 0.50 0.90];

PEXT_MASK_LATE_GROWTH = 0.90;
MIN_SURVIVAL_AT_LATE_START = 0.10;

CHI_REP = [0 0.4 0.7 1.2 1.6 1.9 2.2 2.6];

GROWTH_RATIO_AGGRESSIVE = 0.50;
GROWTH_RATIO_NEAR_ARREST = 0.10;

SEED_BASE = 810000;
SEED_BOOT = 820000;

USE_PARALLEL = true;
AUTO_START_POOL = true;

switch lower(RUN_MODE)
    case 'smoke'
        N_BASE = 12;
        N_MAX = 24;
        MC_BATCH = 12;
        N_ANCHOR = 12;
        N_BOOT_BOUNDARY = 250;

    case 'pilot'
        N_BASE = 50;
        N_MAX = 140;
        MC_BATCH = 30;
        N_ANCHOR = 40;

    case 'final'
        N_BASE = 80;
        N_MAX = 240;
        MC_BATCH = 40;
        N_ANCHOR = 80;
        N_BOOT_BOUNDARY = 3000;

    otherwise
        error('RUN_MODE must be smoke, pilot, or final.');
end

%% 1. LOAD FROZEN FRAMEWORK
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir = pwd; end

coreFile = fullfile(scriptDir,'macro_sim_core_hazard_SSA_FINAL.m');
rFile = fullfile(scriptDir,'Macro_r_Calibration_OUT','Macro_r_LOCKED.mat');
muFile = fullfile(scriptDir,'Macro_mu_Calibration_OUT','Macro_mu_LOCKED.mat');
ssaLockFile = fullfile(scriptDir,'Macro_Hazard_SSA_Domain_Confirm_OUT','Macro_SSA_LOCKED.mat');
interfaceFile = fullfile(scriptDir,'Meso_to_Macro_Risk_Interface_FINAL_OUT', ...
    'Meso_to_Macro_Risk_Interface_FINAL.mat');

assert(isfile(coreFile),'Missing macro_sim_core_hazard_SSA_FINAL.m.');
assert(isfile(rFile),'Missing Macro_r_LOCKED.mat.');
assert(isfile(muFile),'Missing Macro_mu_LOCKED.mat.');
assert(isfile(ssaLockFile),'Missing Macro_SSA_LOCKED.mat.');
assert(isfile(interfaceFile),'Missing meso->macro interface MAT.');

R = load(rFile,'LOCK');
M = load(muFile,'LOCK');
S = load(ssaLockFile);
I = load(interfaceFile,'INTERFACE');

r_h = double(R.LOCK.r_star_per_h);
mu_ref_h = double(M.LOCK.mu_star_per_h);

if isfield(S,'SSA_LOCK')
    SSA_LOCK = S.SSA_LOCK;
elseif isfield(S,'LOCK')
    SSA_LOCK = S.LOCK;
else
    error('Macro_SSA_LOCKED.mat contains neither SSA_LOCK nor LOCK.');
end

assert(isfield(SSA_LOCK,'readyForQ1Pilot') && SSA_LOCK.readyForQ1Pilot, ...
    'Frozen SSA framework is not marked ready for Q1.');

L = double(SSA_LOCK.L);
R0 = double(SSA_LOCK.R0);

assert(L==61,'Targeted Q1 requires frozen L=61.');
assert(abs(double(SSA_LOCK.Tend_h)-TEND_H)<1e-12,'Tend mismatch.');
assert(abs(double(SSA_LOCK.r_h)-r_h)<1e-12,'r mismatch.');
assert(abs(double(SSA_LOCK.mu_reference_h)-mu_ref_h)<1e-12,'mu mismatch.');

%% 2. BUILD chi_mu GRID
MU_BENCHMARK_H = 1./TAU_F_BENCHMARK_H;
CHI_MU_BENCHMARK = MU_BENCHMARK_H/r_h;

CHI_MU_GRID = unique([CHI_MU_BENCHMARK CHI_MU_INTERP],'sorted');
MU_GRID_H = CHI_MU_GRID*r_h;
TAU_F_GRID_H = 1./MU_GRID_H;

isBenchmarkMu = false(size(CHI_MU_GRID));
for k=1:numel(CHI_MU_BENCHMARK)
    [d,idx] = min(abs(CHI_MU_GRID-CHI_MU_BENCHMARK(k)));
    assert(d<1e-10);
    isBenchmarkMu(idx)=true;
end

[~,iMuRef] = min(abs(TAU_F_GRID_H-TAU_F_REFERENCE_H));
assert(abs(MU_GRID_H(iMuRef)-mu_ref_h)<1e-10, ...
    'Reference 12-h clearance is not reproduced.');

%% 3. LOAD 48 MESOSCOPIC STATES
INTERFACE = I.INTERFACE;

tauL = double(INTERFACE.tauL_min(:));
condNames = cellstr(string(INTERFACE.condition.names(:)));
lambdaReal = double(INTERFACE.macro.lambdaF_eq_h);

nTau = numel(tauL);
nCond = numel(condNames);

if isequal(size(lambdaReal),[nTau nCond])
elseif isequal(size(lambdaReal),[nCond nTau])
    lambdaReal=lambdaReal.';
else
    error('Unexpected lambdaF_eq_h size in interface.');
end

realChiF = lambdaReal/r_h;

assert(min(realChiF(:))>CHI_F_GRID(end), ...
    'The targeted generalized regime grid unexpectedly overlaps the current 48-state range.');

%% 4. OUTPUT / CHECKPOINT
outDir = fullfile(scriptDir,'Macro_Q1_TargetedRefinement_FINAL_OUT');
if ~exist(outDir,'dir'), mkdir(outDir); end

outPNG = fullfile(outDir,'Macro_Q1_TargetedRefinement_FINAL.png');
outPDF = fullfile(outDir,'Macro_Q1_TargetedRefinement_FINAL.pdf');
outFIG = fullfile(outDir,'Macro_Q1_TargetedRefinement_FINAL.fig');
outMAT = fullfile(outDir,'Macro_Q1_TargetedRefinement_FINAL.mat');
outMetricsCSV = fullfile(outDir,'Macro_Q1_Targeted_StateMetrics.csv');
outBoundaryCSV = fullfile(outDir,'Macro_Q1_Targeted_Boundaries.csv');
outProjectionCSV = fullfile(outDir,'Macro_Q1_Targeted_MesoscopicProjection.csv');
checkpointFile = fullfile(outDir,'Macro_Q1_Targeted_checkpoint.mat');

%% 5. PARALLEL POOL
poolReady = false;
if USE_PARALLEL && license('test','Distrib_Computing_Toolbox')
    p = gcp('nocreate');
    if isempty(p) && AUTO_START_POOL
        try
            p = parpool('threads');
        catch ME
            warning('Could not start thread pool (%s). Using serial.',ME.message);
            p=[];
        end
    end
    poolReady = ~isempty(p);
end

%% 6. SETTINGS SIGNATURE / RESUME
SETTINGS = struct();
SETTINGS.version='MACRO_Q1_TARGETED_REFINEMENT_v1';
SETTINGS.runMode=RUN_MODE;
SETTINGS.L=L;
SETTINGS.R0=R0;
SETTINGS.Tend_h=TEND_H;
SETTINGS.record_dt_h=RECORD_DT_H;
SETTINGS.lateWindow_h=LATE_WINDOW_H;
SETTINGS.r_h=r_h;
SETTINGS.chiF=CHI_F_GRID;
SETTINGS.chiMu=CHI_MU_GRID;
SETTINGS.mu_h=MU_GRID_H;
SETTINGS.nBase=N_BASE;
SETTINGS.nMax=N_MAX;
SETTINGS.mcBatch=MC_BATCH;
SETTINGS.mcseTarget=MCSE_TARGET;
SETTINGS.pExtHalfTarget=PEXT_WILSON_HALF_TARGET;
SETTINGS.pExtLevels=PEXT_LEVELS;
SETTINGS.pExtLateGrowthMask=PEXT_MASK_LATE_GROWTH;
SETTINGS.minLateStartSurvival=MIN_SURVIVAL_AT_LATE_START;

nChi = numel(CHI_F_GRID);
nMu = numel(CHI_MU_GRID);

RES = cell(nMu,nChi);
COMPLETE = false(nMu,nChi);

if isfile(checkpointFile)
    C = load(checkpointFile,'SETTINGS_CK','RES','COMPLETE');
    if isfield(C,'SETTINGS_CK') && isequaln(C.SETTINGS_CK,SETTINGS)
        RES=C.RES;
        COMPLETE=C.COMPLETE;
        fprintf('Compatible checkpoint found: %d/%d states complete.\n', ...
            nnz(COMPLETE),numel(COMPLETE));
    else
        warning('Existing checkpoint is incompatible; starting a new run.');
    end
end

%% 7. CONSOLE HEADER
fprintf('\n============================================================\n');
fprintf('MACRO Q1 TARGETED REFINEMENT: %s\n',upper(RUN_MODE));
fprintf('Solver: event-driven direct SSA; no numerical dt\n');
fprintf('Locked r = %.8f h^-1; L=%d; Tend=%.1f h\n',r_h,L,TEND_H);
fprintf('chi_F grid (%d): ',nChi); fprintf('%g ',CHI_F_GRID); fprintf('\n');
fprintf('chi_mu grid (%d): ',nMu); fprintf('%.3f ',CHI_MU_GRID); fprintf('\n');
fprintf('mean F residence (h): '); fprintf('%.2f ',TAU_F_GRID_H); fprintf('\n');
fprintf('Real 48-state chi_F range = %.3f--%.3f\n',min(realChiF(:)),max(realChiF(:)));
fprintf('N_base/N_max = %d/%d; Pext CI-half target=%.3f\n', ...
    N_BASE,N_MAX,PEXT_WILSON_HALF_TARGET);
fprintf('============================================================\n');

%% 8. RUN TARGETED GRID
for j=1:nChi
    chiF = CHI_F_GRID(j);
    lambda_h = chiF*r_h;

    if chiF==0
        if ~COMPLETE(1,j)
            cfg = makeCfg(L,R0,TEND_H,RECORD_DT_H,r_h,MU_GRID_H(1),lambda_h);
            fprintf('\nchi_F=0: run once and reuse across all clearance points.\n');
            RES{1,j} = runAdaptiveState(cfg,N_BASE,N_MAX,MC_BATCH, ...
                MCSE_TARGET,PEXT_WILSON_HALF_TARGET, ...
                SEED_BASE+j*10000,poolReady,LATE_WINDOW_H);
            COMPLETE(1,j)=true;
            saveCheckpoint(checkpointFile,SETTINGS,RES,COMPLETE);
        end

        for i=2:nMu
            if ~COMPLETE(i,j)
                RES{i,j}=RES{1,j};
                COMPLETE(i,j)=true;
            end
        end

        saveCheckpoint(checkpointFile,SETTINGS,RES,COMPLETE);
        continue;
    end

    for i=1:nMu
        if COMPLETE(i,j), continue; end

        mu_h = MU_GRID_H(i);
        chiMu = CHI_MU_GRID(i);
        cfg = makeCfg(L,R0,TEND_H,RECORD_DT_H,r_h,mu_h,lambda_h);
        seedState = SEED_BASE + j*10000 + i*500;

        fprintf('\nchi_F=%5.2f | chi_mu=%6.3f | lambda=%.6g | mu=%.6g\n', ...
            chiF,chiMu,lambda_h,mu_h);

        RES{i,j} = runAdaptiveState(cfg,N_BASE,N_MAX,MC_BATCH, ...
            MCSE_TARGET,PEXT_WILSON_HALF_TARGET, ...
            seedState,poolReady,LATE_WINDOW_H);

        COMPLETE(i,j)=true;
        saveCheckpoint(checkpointFile,SETTINGS,RES,COMPLETE);
    end
end

assert(all(COMPLETE(:)),'Targeted Q1 grid is incomplete.');

%% 9. HIGH-PRESSURE ANCHOR
ANCHOR = cell(nMu,1);

fprintf('\n================ HIGH-PRESSURE ANCHOR chi_F=%g ================\n', ...
    CHI_F_HIGH_ANCHOR);

for i=1:nMu
    cfg = makeCfg(L,R0,TEND_H,RECORD_DT_H,r_h,MU_GRID_H(i), ...
        CHI_F_HIGH_ANCHOR*r_h);

    ANCHOR{i}=runFixedState(cfg,N_ANCHOR, ...
        SEED_BASE+900000+i*1000,poolReady,LATE_WINDOW_H);

    fprintf('chi_mu=%6.3f | Pext=%.3f [%0.3f,%0.3f]\n', ...
        CHI_MU_GRID(i),ANCHOR{i}.Pext,ANCHOR{i}.PextCI(1),ANCHOR{i}.PextCI(2));
end

%% 10. ASSEMBLE STATE METRICS
Pext = nan(nMu,nChi);
PextLo = nan(nMu,nChi);
PextHi = nan(nMu,nChi);
AUTC_NV = nan(nMu,nChi);
AUTC_NT = nan(nMu,nChi);
FinalNV = nan(nMu,nChi);
FinalNT = nan(nMu,nChi);
LateSurvGrowth = nan(nMu,nChi);
LateStartSurvival = nan(nMu,nChi);
MCSE_NV = nan(nMu,nChi);
MCSE_NT = nan(nMu,nChi);
NrepUsed = nan(nMu,nChi);

for i=1:nMu
    for j=1:nChi
        Q=RES{i,j};
        Pext(i,j)=Q.Pext;
        PextLo(i,j)=Q.PextCI(1);
        PextHi(i,j)=Q.PextCI(2);
        AUTC_NV(i,j)=Q.meanAUTC_NV;
        AUTC_NT(i,j)=Q.meanAUTC_NT;
        FinalNV(i,j)=Q.meanFinalNV;
        FinalNT(i,j)=Q.meanFinalNT;
        LateSurvGrowth(i,j)=Q.survivorAwareLateGrowth_per_h;
        LateStartSurvival(i,j)=Q.lateStartSurvivalFraction;
        MCSE_NV(i,j)=Q.relMCSE_AUTC_NV;
        MCSE_NT(i,j)=Q.relMCSE_AUTC_NT;
        NrepUsed(i,j)=Q.nRep;
    end
end

iChi0=find(CHI_F_GRID==0,1);
baselineAUTC_NV=AUTC_NV(1,iChi0);
baselineAUTC_NT=AUTC_NT(1,iChi0);
baselineFinalNV=FinalNV(1,iChi0);
baselineFinalNT=FinalNT(1,iChi0);
baselineLateGrowth=LateSurvGrowth(1,iChi0);

assert(baselineAUTC_NV>0);
assert(baselineAUTC_NT>0);
assert(isfinite(baselineLateGrowth) && baselineLateGrowth>0);

AUTC_NV_norm=AUTC_NV/baselineAUTC_NV;
AUTC_NT_norm=AUTC_NT/baselineAUTC_NT;
FinalNV_norm=FinalNV/baselineFinalNV;
FinalNT_norm=FinalNT/baselineFinalNT;
LateSurvGrowthRatio=LateSurvGrowth/baselineLateGrowth;

LateGrowthInterpretationMask = ...
    Pext<PEXT_MASK_LATE_GROWTH & ...
    LateStartSurvival>=MIN_SURVIVAL_AT_LATE_START & ...
    isfinite(LateSurvGrowthRatio);

LateGrowthPlot = LateSurvGrowthRatio;
LateGrowthPlot(~LateGrowthInterpretationMask)=NaN;

%% 11. ISOTONIC EXTINCTION CURVES + BOOTSTRAP BOUNDARIES
PextIso = nan(nMu,nChi);

nLevel=numel(PEXT_LEVELS);
Boundary = nan(nMu,nLevel);
BoundaryLo = nan(nMu,nLevel);
BoundaryHi = nan(nMu,nLevel);

GrowthArrest = nan(nMu,1);
GrowthArrestLo = nan(nMu,1);
GrowthArrestHi = nan(nMu,1);

for i=1:nMu
    weights=zeros(1,nChi);
    for j=1:nChi
        weights(j)=RES{i,j}.nRep;
    end

    PextIso(i,:)=isotonicIncreasing(Pext(i,:),weights);

    for q=1:nLevel
        Boundary(i,q)=thresholdCrossing(CHI_F_GRID,PextIso(i,:),PEXT_LEVELS(q));
    end

    growthForBoundary=LateSurvGrowth(i,:);
    growthForBoundary(~LateGrowthInterpretationMask(i,:))=NaN;
    GrowthArrest(i)=firstZeroCrossingWithNaN(CHI_F_GRID,growthForBoundary);

    [bootP,bootG] = bootstrapBoundariesForMu(RES(i,:),CHI_F_GRID,PEXT_LEVELS, ...
        PEXT_MASK_LATE_GROWTH,MIN_SURVIVAL_AT_LATE_START, ...
        LATE_WINDOW_H,N_BOOT_BOUNDARY,SEED_BOOT+i*10000);

    for q=1:nLevel
        BoundaryLo(i,q)=percentileFinite(bootP(:,q),2.5);
        BoundaryHi(i,q)=percentileFinite(bootP(:,q),97.5);
    end

    GrowthArrestLo(i)=percentileFinite(bootG,2.5);
    GrowthArrestHi(i)=percentileFinite(bootG,97.5);
end

%% 12. OPERATIONAL DYNAMIC LABELS
DYN = strings(nMu,nChi);

for i=1:nMu
    for j=1:nChi
        if Pext(i,j)>=0.90
            DYN(i,j)="extinction-dominated";
        elseif ~LateGrowthInterpretationMask(i,j)
            DYN(i,j)="transition / unresolved late growth";
        elseif LateSurvGrowthRatio(i,j)>=GROWTH_RATIO_AGGRESSIVE
            DYN(i,j)="sustained-growth-like";
        elseif LateSurvGrowthRatio(i,j)>GROWTH_RATIO_NEAR_ARREST
            DYN(i,j)="slow-growth-like";
        elseif abs(LateSurvGrowthRatio(i,j))<=GROWTH_RATIO_NEAR_ARREST
            DYN(i,j)="near-arrest-like";
        else
            DYN(i,j)="persistent-regression-like";
        end
    end
end

%% 13. REPRESENTATIVE TRAJECTORIES
repIdx=nan(size(CHI_REP));
for k=1:numel(CHI_REP)
    [d,idx]=min(abs(CHI_F_GRID-CHI_REP(k)));
    assert(d<1e-12,'Representative chi_F=%g is absent.',CHI_REP(k));
    repIdx(k)=idx;
end

tCurve=RES{iMuRef,repIdx(1)}.t_h;
repMeanNV=nan(numel(tCurve),numel(repIdx));

for k=1:numel(repIdx)
    repMeanNV(:,k)=RES{iMuRef,repIdx(k)}.meanNVcurve;
end

%% 14. MESOSCOPIC PROJECTION
iP50=find(abs(PEXT_LEVELS-0.50)<1e-12,1);
iP90=find(abs(PEXT_LEVELS-0.90)<1e-12,1);

chi50Ref=Boundary(iMuRef,iP50);
chi90Ref=Boundary(iMuRef,iP90);
chiArrestRef=GrowthArrest(iMuRef);

ratioToP50=realChiF/chi50Ref;
ratioToP90=realChiF/chi90Ref;

anchorPext=cellfun(@(q) q.Pext,ANCHOR);
anchorPextLo=cellfun(@(q) q.PextCI(1),ANCHOR);
anchorPextHi=cellfun(@(q) q.PextCI(2),ANCHOR);

%% 15. TABLES
nRows=nMu*nChi;

chiF_col=zeros(nRows,1);
chiMu_col=zeros(nRows,1);
mu_h_col=zeros(nRows,1);
tauF_h_col=zeros(nRows,1);
isBenchmarkMu_col=false(nRows,1);
nRep_col=zeros(nRows,1);
Pext_col=zeros(nRows,1);
Pext_lo_col=zeros(nRows,1);
Pext_hi_col=zeros(nRows,1);
Pext_iso_col=zeros(nRows,1);
AUTC_NV_col=zeros(nRows,1);
AUTC_NV_norm_col=zeros(nRows,1);
AUTC_NT_col=zeros(nRows,1);
AUTC_NT_norm_col=zeros(nRows,1);
FinalNV_col=zeros(nRows,1);
FinalNV_norm_col=zeros(nRows,1);
FinalNT_col=zeros(nRows,1);
FinalNT_norm_col=zeros(nRows,1);
LateSurvGrowth_col=zeros(nRows,1);
LateSurvGrowthRatio_col=zeros(nRows,1);
LateStartSurvival_col=zeros(nRows,1);
LateGrowthInterpretable_col=false(nRows,1);
MCSE_NV_col=zeros(nRows,1);
MCSE_NT_col=zeros(nRows,1);
DynamicLabel_col=strings(nRows,1);

rr=0;
for i=1:nMu
    for j=1:nChi
        rr=rr+1;
        chiF_col(rr)=CHI_F_GRID(j);
        chiMu_col(rr)=CHI_MU_GRID(i);
        mu_h_col(rr)=MU_GRID_H(i);
        tauF_h_col(rr)=TAU_F_GRID_H(i);
        isBenchmarkMu_col(rr)=isBenchmarkMu(i);
        nRep_col(rr)=NrepUsed(i,j);
        Pext_col(rr)=Pext(i,j);
        Pext_lo_col(rr)=PextLo(i,j);
        Pext_hi_col(rr)=PextHi(i,j);
        Pext_iso_col(rr)=PextIso(i,j);
        AUTC_NV_col(rr)=AUTC_NV(i,j);
        AUTC_NV_norm_col(rr)=AUTC_NV_norm(i,j);
        AUTC_NT_col(rr)=AUTC_NT(i,j);
        AUTC_NT_norm_col(rr)=AUTC_NT_norm(i,j);
        FinalNV_col(rr)=FinalNV(i,j);
        FinalNV_norm_col(rr)=FinalNV_norm(i,j);
        FinalNT_col(rr)=FinalNT(i,j);
        FinalNT_norm_col(rr)=FinalNT_norm(i,j);
        LateSurvGrowth_col(rr)=LateSurvGrowth(i,j);
        LateSurvGrowthRatio_col(rr)=LateSurvGrowthRatio(i,j);
        LateStartSurvival_col(rr)=LateStartSurvival(i,j);
        LateGrowthInterpretable_col(rr)=LateGrowthInterpretationMask(i,j);
        MCSE_NV_col(rr)=MCSE_NV(i,j);
        MCSE_NT_col(rr)=MCSE_NT(i,j);
        DynamicLabel_col(rr)=DYN(i,j);
    end
end

Tmetrics=table( ...
    chiF_col,chiMu_col,mu_h_col,tauF_h_col,isBenchmarkMu_col,nRep_col, ...
    Pext_col,Pext_lo_col,Pext_hi_col,Pext_iso_col, ...
    AUTC_NV_col,AUTC_NV_norm_col,AUTC_NT_col,AUTC_NT_norm_col, ...
    FinalNV_col,FinalNV_norm_col,FinalNT_col,FinalNT_norm_col, ...
    LateSurvGrowth_col,LateSurvGrowthRatio_col,LateStartSurvival_col, ...
    LateGrowthInterpretable_col,MCSE_NV_col,MCSE_NT_col,DynamicLabel_col, ...
    'VariableNames',{ ...
    'chiF_lambda_over_r','chiMu_mu_over_r','mu_h','meanFResidence_h', ...
    'isOriginalMuBenchmark','nRep','PextV_168h','Pext_lo','Pext_hi', ...
    'Pext_isotonic','mean_AUTC_NV','AUTC_NV_vs_noF','mean_AUTC_NT', ...
    'AUTC_NT_vs_noF','mean_finalNV','finalNV_vs_noF','mean_finalNT', ...
    'finalNT_vs_noF','survivorAwareLateGrowth_per_h', ...
    'survivorLateGrowth_vs_noF','lateStartSurvivalFraction', ...
    'lateGrowthInterpretable','relMCSE_AUTC_NV','relMCSE_AUTC_NT','dynamicLabel'});

chiMu_b=CHI_MU_GRID(:);
mu_h_b=MU_GRID_H(:);
tauF_h_b=TAU_F_GRID_H(:);
isBenchmark_b=isBenchmarkMu(:);

Tboundary=table( ...
    chiMu_b,mu_h_b,tauF_h_b,isBenchmark_b, ...
    Boundary(:,1),BoundaryLo(:,1),BoundaryHi(:,1), ...
    Boundary(:,2),BoundaryLo(:,2),BoundaryHi(:,2), ...
    Boundary(:,3),BoundaryLo(:,3),BoundaryHi(:,3), ...
    GrowthArrest,GrowthArrestLo,GrowthArrestHi, ...
    anchorPext(:),anchorPextLo(:),anchorPextHi(:), ...
    'VariableNames',{ ...
    'chiMu_mu_over_r','mu_h','meanFResidence_h','isOriginalMuBenchmark', ...
    'chiF_Pext10','chiF_Pext10_lo','chiF_Pext10_hi', ...
    'chiF_Pext50','chiF_Pext50_lo','chiF_Pext50_hi', ...
    'chiF_Pext90','chiF_Pext90_lo','chiF_Pext90_hi', ...
    'chiF_survivorLateGrowthArrest','chiF_growthArrest_lo','chiF_growthArrest_hi', ...
    'Pext_at_chiF8','Pext_chiF8_lo','Pext_chiF8_hi'});

nProj=nTau*nCond;

Condition=strings(nProj,1);
tauL_min_col=zeros(nProj,1);
lambdaReal_col=zeros(nProj,1);
realChiF_col=zeros(nProj,1);
ratioP50_col=zeros(nProj,1);
ratioP90_col=zeros(nProj,1);

rr=0;
for c=1:nCond
    for j=1:nTau
        rr=rr+1;
        Condition(rr)=string(condNames{c});
        tauL_min_col(rr)=tauL(j);
        lambdaReal_col(rr)=lambdaReal(j,c);
        realChiF_col(rr)=realChiF(j,c);
        ratioP50_col(rr)=ratioToP50(j,c);
        ratioP90_col(rr)=ratioToP90(j,c);
    end
end

Tprojection=table( ...
    Condition,tauL_min_col,lambdaReal_col,realChiF_col,ratioP50_col,ratioP90_col, ...
    'VariableNames',{ ...
    'Condition','tauL_min','lambdaF_eq_h','chiF_lambda_over_r', ...
    'chiF_over_referenceMu_Pext50Boundary', ...
    'chiF_over_referenceMu_Pext90Boundary'});

writetable(Tmetrics,outMetricsCSV);
writetable(Tboundary,outBoundaryCSV);
writetable(Tprojection,outProjectionCSV);

%% 16. CONSOLE SUMMARY
fprintf('\n================ TARGETED Q1 SUMMARY ================\n');
fprintf('Reference clearance: mean F residence = %.2f h, chi_mu=%.3f\n', ...
    TAU_F_GRID_H(iMuRef),CHI_MU_GRID(iMuRef));

fprintf('Survivor-aware late-growth arrest chi_F = %s [%s,%s]\n', ...
    numOrNA(chiArrestRef), ...
    numOrNA(GrowthArrestLo(iMuRef)),numOrNA(GrowthArrestHi(iMuRef)));

fprintf('Pext=0.10 chi_F = %.3f [%.3f,%.3f]\n', ...
    Boundary(iMuRef,1),BoundaryLo(iMuRef,1),BoundaryHi(iMuRef,1));
fprintf('Pext=0.50 chi_F = %.3f [%.3f,%.3f]\n', ...
    Boundary(iMuRef,2),BoundaryLo(iMuRef,2),BoundaryHi(iMuRef,2));
fprintf('Pext=0.90 chi_F = %.3f [%.3f,%.3f]\n', ...
    Boundary(iMuRef,3),BoundaryLo(iMuRef,3),BoundaryHi(iMuRef,3));

fprintf('Real 48-state chi_F range = %.3f--%.3f\n', ...
    min(realChiF(:)),max(realChiF(:)));
fprintf('Min real chi_F / reference Pext50 boundary = %.3f\n', ...
    min(ratioToP50(:)));
fprintf('Min real chi_F / reference Pext90 boundary = %.3f\n', ...
    min(ratioToP90(:)));

fprintf('High-pressure anchor chi_F=8: Pext range = %.3f--%.3f\n', ...
    min(anchorPext),max(anchorPext));

fprintf('Maximum MCSE AUTC_NV = %.2f%%\n',100*max(MCSE_NV(:)));
fprintf('Maximum MCSE AUTC_NT = %.2f%%\n',100*max(MCSE_NT(:)));
fprintf('=====================================================\n');

%% 17. FIGURE
fig=figure('Color','w','Position',[60 50 1510 860]);
tl=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');

% A. Refined extinction landscape
ax=nexttile(tl,1);
contourf(ax,CHI_F_GRID,CHI_MU_GRID,PextIso,0:0.1:1,'LineStyle','none');
hold(ax,'on');
plot(ax,Boundary(:,2),CHI_MU_GRID,'k-','LineWidth',1.8);
plot(ax,Boundary(:,1),CHI_MU_GRID,'k:','LineWidth',1.2);
plot(ax,Boundary(:,3),CHI_MU_GRID,'k--','LineWidth',1.2);
scatter(ax,Boundary(isBenchmarkMu,2),CHI_MU_GRID(isBenchmarkMu),28,'k','filled');
xlabel(ax,'\chi_F=\lambda_F/r');
ylabel(ax,'\chi_\mu=\mu/r');
title(ax,'a  Refined viable-cell extinction landscape');
cb=colorbar(ax); cb.Label.String='isotonic P(V extinction by 168 h)';
xlim(ax,[CHI_F_GRID(1) CHI_F_GRID(end)]);
ylim(ax,[CHI_MU_GRID(1) CHI_MU_GRID(end)]);
box(ax,'off');

% B. Survivor-aware late growth
ax=nexttile(tl,2);
contourf(ax,CHI_F_GRID,CHI_MU_GRID,LateGrowthPlot,14,'LineStyle','none');
hold(ax,'on');
plot(ax,GrowthArrest,CHI_MU_GRID,'k-','LineWidth',1.7);
xlabel(ax,'\chi_F=\lambda_F/r');
ylabel(ax,'\chi_\mu=\mu/r');
title(ax,'b  Survivor-aware late viable growth');
cb=colorbar(ax); cb.Label.String='late growth / no-ferroptosis baseline';
xlim(ax,[CHI_F_GRID(1) CHI_F_GRID(end)]);
ylim(ax,[CHI_MU_GRID(1) CHI_MU_GRID(end)]);
box(ax,'off');

% C. Integrated viable burden
ax=nexttile(tl,3);
contourf(ax,CHI_F_GRID,CHI_MU_GRID,AUTC_NV_norm,14,'LineStyle','none');
hold(ax,'on');
plot(ax,Boundary(:,2),CHI_MU_GRID,'k-','LineWidth',1.5);
xlabel(ax,'\chi_F=\lambda_F/r');
ylabel(ax,'\chi_\mu=\mu/r');
title(ax,'c  Integrated viable tumour burden');
cb=colorbar(ax); cb.Label.String='AUTC_{NV} / no-ferroptosis baseline';
xlim(ax,[CHI_F_GRID(1) CHI_F_GRID(end)]);
ylim(ax,[CHI_MU_GRID(1) CHI_MU_GRID(end)]);
box(ax,'off');

% D. Transition boundaries with CI
ax=nexttile(tl,4);
hold(ax,'on');
plotBoundaryCI(ax,CHI_MU_GRID,Boundary(:,1),BoundaryLo(:,1),BoundaryHi(:,1), ...
    '-o','P_{ext}=0.10');
plotBoundaryCI(ax,CHI_MU_GRID,Boundary(:,2),BoundaryLo(:,2),BoundaryHi(:,2), ...
    '-s','P_{ext}=0.50');
plotBoundaryCI(ax,CHI_MU_GRID,Boundary(:,3),BoundaryLo(:,3),BoundaryHi(:,3), ...
    '-^','P_{ext}=0.90');
plotBoundaryCI(ax,CHI_MU_GRID,GrowthArrest,GrowthArrestLo,GrowthArrestHi, ...
    '--d','late-growth arrest');
scatter(ax,CHI_MU_GRID(isBenchmarkMu),Boundary(isBenchmarkMu,2), ...
    25,'k','filled','HandleVisibility','off');
xlabel(ax,'\chi_\mu=\mu/r');
ylabel(ax,'Boundary \chi_F');
title(ax,'d  Refined transition boundaries');
legend(ax,'Location','best','Box','off');
grid(ax,'on');
box(ax,'off');

% E. Representative trajectories
ax=nexttile(tl,5);
hold(ax,'on');
for k=1:numel(CHI_REP)
    plot(ax,tCurve,repMeanNV(:,k),'LineWidth',1.45, ...
        'DisplayName',sprintf('\\chi_F=%g',CHI_REP(k)));
end
xlabel(ax,'Time (h)');
ylabel(ax,'Mean viable cells N_V');
title(ax,sprintf('e  Dynamics at reference \\chi_\\mu=%.3f',CHI_MU_GRID(iMuRef)));
legend(ax,'Location','best','Box','off');
grid(ax,'on');
box(ax,'off');

% F. Cross-scale projection
ax=nexttile(tl,6);
hold(ax,'on');
markers={'o','s','^','d'};
for c=1:nCond
    semilogy(ax,tauL,realChiF(:,c),['-' markers{min(c,numel(markers))}], ...
        'LineWidth',1.45,'MarkerFaceColor','w','DisplayName',condNames{c});
end
yline(ax,chi50Ref,'k-','P_{ext}=0.50','LineWidth',1.4);
yline(ax,chi90Ref,'k--','P_{ext}=0.90','LineWidth',1.2);
yline(ax,CHI_F_HIGH_ANCHOR,':','\chi_F=8 anchor','LineWidth',1.0);
xlabel(ax,'Lipid-memory time \tau_L (min)');
ylabel(ax,'Mesoscopic \chi_F=\lambda_F^{eq}/r');
title(ax,'f  Cross-scale projection of pro-ferroptotic states');
legend(ax,'Location','northwest','Box','off');
grid(ax,'on');
box(ax,'off');

sgtitle(tl, ...
    'Q1 targeted refinement: commitment pressure sets growth and extinction transitions');

%% 18. SAVE
try
    axs=findall(fig,'Type','axes');
    for k=1:numel(axs), axtoolbar(axs(k),{}); end
catch
end

exportgraphics(fig,outPNG,'Resolution',600);
exportgraphics(fig,outPDF,'ContentType','vector');
savefig(fig,outFIG);

Q1R=struct();
Q1R.version='MACRO_Q1_TARGETED_REFINEMENT_v1';
Q1R.created=datestr(now,30);
Q1R.runMode=RUN_MODE;
Q1R.settings=SETTINGS;
Q1R.chiF=CHI_F_GRID;
Q1R.chiMu=CHI_MU_GRID;
Q1R.mu_h=MU_GRID_H;
Q1R.meanFResidence_h=TAU_F_GRID_H;
Q1R.isOriginalMuBenchmark=isBenchmarkMu;
Q1R.Pext=Pext;
Q1R.PextLo=PextLo;
Q1R.PextHi=PextHi;
Q1R.PextIsotonic=PextIso;
Q1R.AUTC_NV=AUTC_NV;
Q1R.AUTC_NT=AUTC_NT;
Q1R.AUTC_NV_norm=AUTC_NV_norm;
Q1R.AUTC_NT_norm=AUTC_NT_norm;
Q1R.FinalNV=FinalNV;
Q1R.FinalNT=FinalNT;
Q1R.survivorAwareLateGrowth=LateSurvGrowth;
Q1R.survivorLateGrowthRatio=LateSurvGrowthRatio;
Q1R.lateStartSurvivalFraction=LateStartSurvival;
Q1R.lateGrowthInterpretationMask=LateGrowthInterpretationMask;
Q1R.boundary.pLevels=PEXT_LEVELS;
Q1R.boundary.pExt=Boundary;
Q1R.boundary.pExtLo=BoundaryLo;
Q1R.boundary.pExtHi=BoundaryHi;
Q1R.boundary.growthArrest=GrowthArrest;
Q1R.boundary.growthArrestLo=GrowthArrestLo;
Q1R.boundary.growthArrestHi=GrowthArrestHi;
Q1R.highPressureAnchor.chiF=CHI_F_HIGH_ANCHOR;
Q1R.highPressureAnchor.Pext=anchorPext;
Q1R.highPressureAnchor.PextLo=anchorPextLo;
Q1R.highPressureAnchor.PextHi=anchorPextHi;
Q1R.dynamicLabel=DYN;
Q1R.realMesoscopic.tauL_min=tauL;
Q1R.realMesoscopic.conditionNames=condNames;
Q1R.realMesoscopic.lambdaF_eq_h=lambdaReal;
Q1R.realMesoscopic.chiF=realChiF;
Q1R.referenceMu.index=iMuRef;
Q1R.referenceMu.chiMu=CHI_MU_GRID(iMuRef);
Q1R.referenceMu.mu_h=MU_GRID_H(iMuRef);
Q1R.referenceMu.meanFResidence_h=TAU_F_GRID_H(iMuRef);
Q1R.tables.metrics=Tmetrics;
Q1R.tables.boundaries=Tboundary;
Q1R.tables.projection=Tprojection;

save(outMAT,'Q1R','RES','ANCHOR','Tmetrics','Tboundary','Tprojection','-v7.3');

fprintf('\nSaved targeted-Q1 outputs:\n');
fprintf('  %s\n',outPNG);
fprintf('  %s\n',outPDF);
fprintf('  %s\n',outFIG);
fprintf('  %s\n',outMAT);
fprintf('  %s\n',outMetricsCSV);
fprintf('  %s\n',outBoundaryCSV);
fprintf('  %s\n',outProjectionCSV);

fprintf('\nNEXT STEP:\n');
fprintf(['Inspect whether sustained-growth, near-arrest, persistent-regression, ', ...
    'stochastic-extinction, and extinction-dominated sectors emerge naturally. ', ...
    'Do not rename them aggressive/bounded/indolent until the refined dynamics ', ...
    'support those labels.\n']);

end

%% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function cfg=makeCfg(L,R0,Tend,recordDt,r,mu,lambda)
cfg=struct();
cfg.L=L;
cfg.R0=R0;
cfg.Tend_h=Tend;
cfg.record_dt_h=recordDt;
cfg.r_h=r;
cfg.mu_h=mu;
cfg.seed=1;
cfg.lambda_mode='scalar';
cfg.lambdaF_h=lambda;
cfg.snapshot_times_h=[];
cfg.return_final_state=false;
cfg.compute_cavity=false;
cfg.store_event_log=false;
cfg.verbose=false;
end

function Q=runAdaptiveState(cfg,nBase,nMax,batchSize,mcseTarget,pExtHalfTarget, ...
    seedBase,poolReady,lateWindow)

E=runBatch(cfg,nBase,seedBase,poolReady);
nCurrent=nBase;

while true
    Q=summarizeEnsemble(E,lateWindow);

    precisionOK = ...
        Q.relMCSE_AUTC_NV<=mcseTarget && ...
        Q.relMCSE_AUTC_NT<=mcseTarget && ...
        Q.PextCIHalf<=pExtHalfTarget && ...
        all(E.bookkeepingPASS);

    if precisionOK || nCurrent>=nMax, break; end

    nAdd=min(batchSize,nMax-nCurrent);
    Eadd=runBatch(cfg,nAdd,seedBase+nCurrent,poolReady);
    E=appendBatch(E,Eadd);
    nCurrent=nCurrent+nAdd;
end

Q=summarizeEnsemble(E,lateWindow);
Q.nRep=nCurrent;

fprintf(['  n=%3d | Pext=%.3f [%0.3f,%0.3f] | ', ...
    'MCSE NV/NT=%.2f/%.2f%% | late-start surv=%.3f | g_late,surv=%s\n'], ...
    Q.nRep,Q.Pext,Q.PextCI(1),Q.PextCI(2), ...
    100*Q.relMCSE_AUTC_NV,100*Q.relMCSE_AUTC_NT, ...
    Q.lateStartSurvivalFraction,numOrNA(Q.survivorAwareLateGrowth_per_h));
end

function Q=runFixedState(cfg,nRep,seedBase,poolReady,lateWindow)
E=runBatch(cfg,nRep,seedBase,poolReady);
Q=summarizeEnsemble(E,lateWindow);
Q.nRep=nRep;
end

function E=runBatch(cfg,nRep,seedBase,poolReady)
tmp=cell(nRep,1);

if poolReady
    parfor k=1:nRep
        c=cfg;
        c.seed=seedBase+k;
        O=macro_sim_core_hazard_SSA_FINAL(c);
        tmp{k}=reduceSSA(O);
    end
else
    for k=1:nRep
        c=cfg;
        c.seed=seedBase+k;
        O=macro_sim_core_hazard_SSA_FINAL(c);
        tmp{k}=reduceSSA(O);
    end
end

E=stackReduced(tmp);
end

function z=reduceSSA(O)
z=struct();
z.t_h=double(O.t_h(:));
z.NV=double(O.NV(:));
z.NT=double(O.NT(:));
z.finalNV=double(O.NV(end));
z.finalNT=double(O.NT(end));
z.AUTC_NV=double(O.AUTC_NV);
z.AUTC_NT=double(O.AUTC_NT);
z.extinctV=isfinite(O.viableExtinctionTime_h) && O.viableExtinctionTime_h<=O.cfg.Tend_h;
z.edgeTouch=any(O.edgeOcc);
z.bookkeepingPASS=logical(O.audit.bookkeepingPASS);
end

function E=stackReduced(tmp)
n=numel(tmp);
E=struct();
E.t_h=tmp{1}.t_h;
nT=numel(E.t_h);
E.NVcurve=zeros(nT,n);
E.NTcurve=zeros(nT,n);
E.finalNV=zeros(n,1);
E.finalNT=zeros(n,1);
E.AUTC_NV=zeros(n,1);
E.AUTC_NT=zeros(n,1);
E.extinctV=false(n,1);
E.edgeTouch=false(n,1);
E.bookkeepingPASS=false(n,1);

for k=1:n
    z=tmp{k};
    assert(isequal(z.t_h,E.t_h));
    E.NVcurve(:,k)=z.NV;
    E.NTcurve(:,k)=z.NT;
    E.finalNV(k)=z.finalNV;
    E.finalNT(k)=z.finalNT;
    E.AUTC_NV(k)=z.AUTC_NV;
    E.AUTC_NT(k)=z.AUTC_NT;
    E.extinctV(k)=z.extinctV;
    E.edgeTouch(k)=z.edgeTouch;
    E.bookkeepingPASS(k)=z.bookkeepingPASS;
end
end

function E=appendBatch(A,B)
assert(isequal(A.t_h,B.t_h));
E=struct();
E.t_h=A.t_h;
E.NVcurve=[A.NVcurve B.NVcurve];
E.NTcurve=[A.NTcurve B.NTcurve];
E.finalNV=[A.finalNV;B.finalNV];
E.finalNT=[A.finalNT;B.finalNT];
E.AUTC_NV=[A.AUTC_NV;B.AUTC_NV];
E.AUTC_NT=[A.AUTC_NT;B.AUTC_NT];
E.extinctV=[A.extinctV;B.extinctV];
E.edgeTouch=[A.edgeTouch;B.edgeTouch];
E.bookkeepingPASS=[A.bookkeepingPASS;B.bookkeepingPASS];
end

function Q=summarizeEnsemble(E,lateWindow)
Q=struct();
Q.t_h=E.t_h;
Q.meanNVcurve=mean(E.NVcurve,2);
Q.meanNTcurve=mean(E.NTcurve,2);
Q.meanFinalNV=mean(E.finalNV);
Q.meanFinalNT=mean(E.finalNT);
Q.meanAUTC_NV=mean(E.AUTC_NV);
Q.meanAUTC_NT=mean(E.AUTC_NT);
Q.relMCSE_AUTC_NV=relativeMCSE(E.AUTC_NV);
Q.relMCSE_AUTC_NT=relativeMCSE(E.AUTC_NT);
Q.Pext=mean(E.extinctV);
[lo,hi]=wilsonCI(sum(E.extinctV),numel(E.extinctV),0.05);
Q.PextCI=[lo hi];
Q.PextCIHalf=(hi-lo)/2;
Q.edgeTouchProb=mean(E.edgeTouch);
Q.bookkeepingPASS=all(E.bookkeepingPASS);
Q.rawExtinctV=logical(E.extinctV(:));

[~,i1]=min(abs(E.t_h-lateWindow(1)));
[~,i2]=min(abs(E.t_h-lateWindow(2)));

dt=E.t_h(i2)-E.t_h(i1);

NV1=E.NVcurve(i1,:);
NV2=E.NVcurve(i2,:);

Q.rawLateNV1=double(NV1(:));
Q.rawLateNV2=double(NV2(:));

aliveAtStart=NV1>0;
Q.lateStartSurvivalFraction=mean(aliveAtStart);

if nnz(aliveAtStart)>=2
    meanN1=mean(NV1(aliveAtStart));
    meanN2=mean(NV2(aliveAtStart));

    if meanN1>0 && meanN2>0
        Q.survivorAwareLateGrowth_per_h=log(meanN2/meanN1)/dt;
    else
        Q.survivorAwareLateGrowth_per_h=NaN;
    end
else
    Q.survivorAwareLateGrowth_per_h=NaN;
end
end

function m=relativeMCSE(x)
x=double(x(:));
mu=mean(x);
if abs(mu)<=eps
    m=0;
else
    m=std(x,0)/sqrt(numel(x))/abs(mu);
end
end

function [lo,hi]=wilsonCI(k,n,alpha)
if abs(alpha-0.05)<1e-12
    z=1.959963984540054;
else
    z=sqrt(2)*erfcinv(alpha);
end

p=k/n;
den=1+z^2/n;
ctr=(p+z^2/(2*n))/den;
half=z*sqrt(p*(1-p)/n+z^2/(4*n^2))/den;
lo=max(0,ctr-half);
hi=min(1,ctr+half);
end

function yIso=isotonicIncreasing(y,w)
y=double(y(:).');
w=double(w(:).');

n=numel(y);
blockValue=zeros(1,n);
blockWeight=zeros(1,n);
blockStart=zeros(1,n);
blockEnd=zeros(1,n);
m=0;

for i=1:n
    m=m+1;
    blockValue(m)=y(i);
    blockWeight(m)=w(i);
    blockStart(m)=i;
    blockEnd(m)=i;

    while m>=2 && blockValue(m-1)>blockValue(m)
        newW=blockWeight(m-1)+blockWeight(m);
        newV=(blockValue(m-1)*blockWeight(m-1)+blockValue(m)*blockWeight(m))/newW;
        blockValue(m-1)=newV;
        blockWeight(m-1)=newW;
        blockEnd(m-1)=blockEnd(m);
        m=m-1;
    end
end

yIso=zeros(1,n);
for b=1:m
    yIso(blockStart(b):blockEnd(b))=blockValue(b);
end
end

function xCross=thresholdCrossing(x,y,target)
x=double(x(:));
y=double(y(:));
xCross=NaN;

if y(1)>=target
    xCross=x(1);
    return;
end

for k=1:numel(x)-1
    if y(k)<=target && y(k+1)>=target
        if abs(y(k+1)-y(k))<eps
            xCross=(x(k)+x(k+1))/2;
        else
            a=(target-y(k))/(y(k+1)-y(k));
            xCross=x(k)+a*(x(k+1)-x(k));
        end
        return;
    end
end
end

function xCross=firstZeroCrossingWithNaN(x,y)
x=double(x(:));
y=double(y(:));
xCross=NaN;

for k=1:numel(x)-1
    if ~isfinite(y(k)) || ~isfinite(y(k+1)), continue; end

    if y(k)==0
        xCross=x(k);
        return;
    end

    if y(k)>0 && y(k+1)<=0
        if abs(y(k+1)-y(k))<eps
            xCross=(x(k)+x(k+1))/2;
        else
            a=(0-y(k))/(y(k+1)-y(k));
            xCross=x(k)+a*(x(k+1)-x(k));
        end
        return;
    end
end
end

function [bootP,bootG]=bootstrapBoundariesForMu( ...
    rowRES,chiF,pLevels,pMask,minLateSurvival,lateWindow,nBoot,seed)

rng(seed,'twister');

nChi=numel(chiF);
nLevel=numel(pLevels);
bootP=nan(nBoot,nLevel);
bootG=nan(nBoot,1);

for b=1:nBoot
    pB=nan(1,nChi);
    gB=nan(1,nChi);
    wB=nan(1,nChi);

    for j=1:nChi
        Q=rowRES{j};
        ext=double(Q.rawExtinctV(:));
        n=numel(ext);

        idx=randi(n,n,1);
        extB=ext(idx);

        pB(j)=mean(extB);
        wB(j)=n;

        NV1=double(Q.rawLateNV1(:));
        NV2=double(Q.rawLateNV2(:));

        NV1B=NV1(idx);
        NV2B=NV2(idx);

        alive=NV1B>0;
        lateSurv=mean(alive);

        if pB(j)<pMask && lateSurv>=minLateSurvival && nnz(alive)>=2
            m1=mean(NV1B(alive));
            m2=mean(NV2B(alive));

            if m1>0 && m2>0
                dt=lateWindow(2)-lateWindow(1);
                gB(j)=log(m2/m1)/dt;
            end
        end
    end

    pIso=isotonicIncreasing(pB,wB);

    for q=1:nLevel
        bootP(b,q)=thresholdCrossing(chiF,pIso,pLevels(q));
    end

    bootG(b)=firstZeroCrossingWithNaN(chiF,gB);
end
end

function q=percentileFinite(v,p)
v=v(isfinite(v));

if isempty(v)
    q=NaN;
    return;
end

v=sort(v(:));
pos=1+(numel(v)-1)*p/100;
a=floor(pos);
b=ceil(pos);

if a==b
    q=v(a);
else
    q=v(a)+(pos-a)*(v(b)-v(a));
end
end

function plotBoundaryCI(ax,x,y,lo,hi,lineSpec,labelText)
x=double(x(:));
y=double(y(:));
lo=double(lo(:));
hi=double(hi(:));

mask=isfinite(x) & isfinite(y) & isfinite(lo) & isfinite(hi);

if nnz(mask)>=2
    xx=x(mask);
    yy=y(mask);
    ll=lo(mask);
    hh=hi(mask);

    patch(ax,[xx;flipud(xx)],[ll;flipud(hh)], ...
        [0.75 0.75 0.75],'FaceAlpha',0.15,'EdgeColor','none', ...
        'HandleVisibility','off');

    plot(ax,xx,yy,lineSpec,'LineWidth',1.5, ...
        'MarkerFaceColor','w','DisplayName',labelText);

elseif nnz(mask)==1
    idx=find(mask,1);
    plot(ax,x(idx),y(idx),lineSpec,'LineWidth',1.5, ...
        'MarkerFaceColor','w','DisplayName',labelText);
end
end

function saveCheckpoint(file,SETTINGS,RES,COMPLETE)
SETTINGS_CK=SETTINGS; %#ok<NASGU>
save(file,'SETTINGS_CK','RES','COMPLETE','-v7.3');
end

function s=numOrNA(x)
if isfinite(x)
    s=sprintf('%.4g',x);
else
    s='NA';
end
end
