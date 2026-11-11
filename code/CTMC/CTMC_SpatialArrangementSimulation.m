function CTMC_SpatialArrangementSimulation()
% =========================================================================
% CTMC_SpatialArrangementSimulation
%
% Q2 primary causal test:
% Does the SPATIAL ARRANGEMENT of cell-level ferroptotic susceptibility
% alter macroscopic tumour growth and morphology when the hazard histogram,
% mean hazard, initial tumour geometry, r, mu, L and observation time are
% all held fixed?
%
% FROZEN MACRO MODEL
% ------------------
% V --(r*q_i)----------> V + V
% V --(lambda_F,i)-----> F
% F --(mu)-------------> empty
%
% Solver:
%   CTMC_EventDrivenSimulation.m
%   exact event-driven direct-method SSA.
%
% Q2 PRIMARY HETEROGENEITY MODEL
% ------------------------------
% We use lambda_mode = 'cell_mark' with daughter_lambda_rule = 'inherit'.
% Thus lambda_F is a persistent CELL-LEVEL susceptibility mark.
%
% Four matched counterfactuals are compared:
%   1) Homogeneous: every initial cell has lambda_bar.
%   2) Dispersed: same low/high lambda_F counts, spatially anti-assorted.
%   3) Random: same low/high lambda_F counts, randomly arranged.
%   4) Clustered: same low/high lambda_F counts, spatially assorted.
%
% Low/high counts are also matched WITHIN RADIAL SHELLS before spatial
% optimization. Hence clustered/dispersed/random conditions have the same
% global hazard histogram and nearly identical radial profile; the intended
% manipulated variable is nearest-neighbour arrangement.
%
% Q1-ANCHORED HAZARDS
% -------------------
% This is a GENERALIZED macro mechanism experiment anchored to Q1
% transition boundaries. It does NOT rescale the actual 48 pro-ferroptotic
% mesoscopic states.
%
% Default:
%   chi_low  = (2/3)*chi_F^arrest
%   chi_high = chi_F^90
%
% With an approximately 50:50 mixture, lambda_bar is fixed exactly from the
% realized initial counts. A homogeneous lambda_bar control is included.
%
% OUTPUTS
% -------
% Macro_Q2_SpatialArrangement_FINAL_OUT/
%   Macro_Q2_SpatialArrangement_FINAL.mat
%   Macro_Q2_SpatialArrangement_Summary.csv
%   Macro_Q2_SpatialArrangement_PairedEffects.csv
%   Macro_Q2_SpatialArrangement_FINAL.png
%   Macro_Q2_SpatialArrangement_PatternExamples_FINAL.png
% =========================================================================

clc;
close all;

%% 0. FINAL SETTINGS
% IMPORTANT:
% The spatial manipulation itself is frozen to the PILOT definition.
% Therefore N_SWAP_ATTEMPTS is kept at 18000.  Increasing it in FINAL would
% silently strengthen the clustering/dispersal intervention and would no
% longer be a confirmatory replication of the same Q2 experiment.
RUN_MODE = 'final';

TEND_H = 168;
RECORD_DT_H = 1;
SNAPSHOT_TIMES_H = [0 72 120 168];

P_HIGH = 0.50;
LOW_RULE  = 'two_thirds_arrest';
HIGH_RULE = 'P90';

N_REP = 100;
N_SWAP_ATTEMPTS = 18000;     % LOCKED manipulation intensity from PILOT
N_BOOT = 10000;

% Independent FINAL seed families: do not recycle the exploratory PILOT
% realizations into the confirmatory ensemble.
SEED_PATTERN_BASE = 1910000;
SEED_SIM_BASE     = 1920000;
SEED_BOOT         = 1930000;

USE_PARALLEL = true;
AUTO_START_POOL = true;

arrNames = {'Homogeneous','Dispersed','Random','Clustered'};
nArr = numel(arrNames);

%% 1. LOAD FROZEN MACRO FRAMEWORK + FINAL Q1 BOUNDARIES
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir = pwd; end

coreFile = fullfile(scriptDir,'CTMC_EventDrivenSimulation.m');
rFile = fullfile(scriptDir,'Macro_r_Calibration_OUT','Macro_r_LOCKED.mat');
muFile = fullfile(scriptDir,'Macro_mu_Calibration_OUT','Macro_mu_LOCKED.mat');
ssaLockFile = fullfile(scriptDir,'Macro_Hazard_SSA_Domain_Confirm_OUT','Macro_SSA_LOCKED.mat');

q1Candidates = { ...
    fullfile(scriptDir,'Macro_Q1_TargetedRefinement_FINAL_OUT', ...
        'Macro_Q1_TargetedRefinement_FINAL.mat'), ...
    fullfile(scriptDir,'Macro_Q1_TargetedRefinement_FINAL.mat')};

assert(isfile(coreFile),'Missing CTMC_EventDrivenSimulation.m.');
assert(isfile(rFile),'Missing Macro_r_LOCKED.mat.');
assert(isfile(muFile),'Missing Macro_mu_LOCKED.mat.');
assert(isfile(ssaLockFile),'Missing Macro_SSA_LOCKED.mat.');

q1File = '';
for k=1:numel(q1Candidates)
    if isfile(q1Candidates{k})
        q1File = q1Candidates{k};
        break;
    end
end
assert(~isempty(q1File),'Missing Macro_Q1_TargetedRefinement_FINAL.mat.');

R = load(rFile,'LOCK');
M = load(muFile,'LOCK');
SLOCK = load(ssaLockFile);
QF = load(q1File,'Q1R');

r_h = double(R.LOCK.r_star_per_h);
mu_h = double(M.LOCK.mu_star_per_h);
Q = QF.Q1R;

if isfield(SLOCK,'SSA_LOCK')
    SSA_LOCK = SLOCK.SSA_LOCK;
elseif isfield(SLOCK,'LOCK')
    SSA_LOCK = SLOCK.LOCK;
else
    error('Macro_SSA_LOCKED.mat contains neither SSA_LOCK nor LOCK.');
end

L = double(SSA_LOCK.L);
R0 = double(SSA_LOCK.R0);

assert(L==61,'Q2 requires the frozen L=61 domain.');
assert(abs(double(SSA_LOCK.r_h)-r_h)<1e-12,'Locked r mismatch.');
assert(abs(double(SSA_LOCK.mu_reference_h)-mu_h)<1e-12,'Locked mu mismatch.');

%% 2. EXTRACT REFERENCE Q1 BOUNDARIES AND DEFINE Q2 HAZARDS
chiMu_Q1 = double(Q.chiMu(:)).';
nMuQ = numel(chiMu_Q1);

pLevels = double(Q.boundary.pLevels(:)).';
B = orientBoundaryMatrix(double(Q.boundary.pExt),nMuQ,numel(pLevels), ...
    'Q1R.boundary.pExt');
gArrest = double(Q.boundary.growthArrest(:));

iP50 = find(abs(pLevels-0.50)<1e-12,1);
iP90 = find(abs(pLevels-0.90)<1e-12,1);
iMuRef = round(double(Q.referenceMu.index));

chiArrest = gArrest(iMuRef);
chiP50 = B(iMuRef,iP50);
chiP90 = B(iMuRef,iP90);
chiMuRef = double(Q.referenceMu.chiMu);

switch lower(LOW_RULE)
    case 'two_thirds_arrest'
        chiLow = (2/3)*chiArrest;
    otherwise
        error('Unknown LOW_RULE.');
end

switch lower(HIGH_RULE)
    case 'p50'
        chiHigh = chiP50;
    case 'p90'
        chiHigh = chiP90;
    otherwise
        error('HIGH_RULE must be P50 or P90.');
end

lambdaLow_h = chiLow*r_h;
lambdaHigh_h = chiHigh*r_h;

%% 3. INITIAL GEOMETRY / EXACT HAZARD COUNTS
ctr = [(L+1)/2 (L+1)/2 (L+1)/2];
[xg,yg,zg] = ndgrid(1:L,1:L,1:L);
r2 = (xg-ctr(1)).^2 + (yg-ctr(2)).^2 + (zg-ctr(3)).^2;
initialMask = r2 <= R0^2;
N0 = nnz(initialMask);

nHigh = round(P_HIGH*N0);
nLow = N0-nHigh;
chiBar = (nLow*chiLow + nHigh*chiHigh)/N0;
lambdaBar_h = chiBar*r_h;

fprintf('\n============================================================\n');
fprintf('MACRO Q2 SPATIAL-ARRANGEMENT TEST: %s\n',upper(RUN_MODE));
fprintf('Solver: exact event-driven direct-method SSA\n');
fprintf('Frozen r=%.8f h^-1; mu=%.8f h^-1; L=%d; R0=%g; Tend=%.1f h\n', ...
    r_h,mu_h,L,R0,TEND_H);
fprintf('Reference chi_mu = %.4f\n',chiMuRef);
fprintf('Q1 boundaries: arrest=%.4f, P50=%.4f, P90=%.4f\n', ...
    chiArrest,chiP50,chiP90);
fprintf('Q2 marks: chi_low=%.4f, chi_high=%.4f, exact mean chi=%.4f\n', ...
    chiLow,chiHigh,chiBar);
fprintf('Initial counts low/high = %d/%d of N0=%d\n',nLow,nHigh,N0);
fprintf('Primary manipulation: spatial arrangement only.\n');
fprintf('============================================================\n');

%% 4. PARALLEL POOL
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

%% 5. ONE PATTERN EXAMPLE FOR VISUAL AUDIT
PATTERN_EXAMPLE = buildMatchedPatterns( ...
    L,R0,chiLow,chiHigh,P_HIGH,SEED_PATTERN_BASE+1,N_SWAP_ATTEMPTS);

fprintf('\nPattern example audit:\n');
fprintf('  Dispersed Moran I = %+0.4f; same-edge fraction = %.3f\n', ...
    PATTERN_EXAMPLE.metrics.dispersed.moranI, ...
    PATTERN_EXAMPLE.metrics.dispersed.sameEdgeFraction);
fprintf('  Random    Moran I = %+0.4f; same-edge fraction = %.3f\n', ...
    PATTERN_EXAMPLE.metrics.random.moranI, ...
    PATTERN_EXAMPLE.metrics.random.sameEdgeFraction);
fprintf('  Clustered Moran I = %+0.4f; same-edge fraction = %.3f\n', ...
    PATTERN_EXAMPLE.metrics.clustered.moranI, ...
    PATTERN_EXAMPLE.metrics.clustered.sameEdgeFraction);

%% 6. RUN PAIRED Q2 ENSEMBLES
tmp = cell(N_REP,1);

if poolReady
    parfor rep=1:N_REP
        tmp{rep} = runOneRep( ...
            rep,L,R0,TEND_H,RECORD_DT_H,SNAPSHOT_TIMES_H, ...
            r_h,mu_h,chiLow,chiHigh,P_HIGH, ...
            N_SWAP_ATTEMPTS,SEED_PATTERN_BASE,SEED_SIM_BASE);
    end
else
    for rep=1:N_REP
        tmp{rep} = runOneRep( ...
            rep,L,R0,TEND_H,RECORD_DT_H,SNAPSHOT_TIMES_H, ...
            r_h,mu_h,chiLow,chiHigh,P_HIGH, ...
            N_SWAP_ATTEMPTS,SEED_PATTERN_BASE,SEED_SIM_BASE);
        fprintf('  completed replicate %d/%d\n',rep,N_REP);
    end
end

E = stackReplicates(tmp,nArr);
assert(all(E.bookkeepingPASS(:)), ...
    'At least one Q2 SSA realization failed bookkeeping.');
assert(~any(E.edgeTouch(:)), ...
    'Boundary touch detected; L=61 cannot be used for this Q2 setting.');

% Confirm that the FINAL independent ensemble reproduces the same spatial
% manipulation class used in PILOT, without strengthening the intervention.
meanMoranDisp = mean(E.initialMoran(:,2),'omitnan');
meanMoranRand = mean(E.initialMoran(:,3),'omitnan');
meanMoranClus = mean(E.initialMoran(:,4),'omitnan');
PATTERN_PASS = meanMoranDisp < -0.60 && abs(meanMoranRand) < 0.10 && meanMoranClus > 0.60;

fprintf('\nFINAL spatial-manipulation audit:\n');
fprintf('  mean Moran I: Dispersed=%+.3f, Random=%+.3f, Clustered=%+.3f\n', ...
    meanMoranDisp,meanMoranRand,meanMoranClus);
fprintf('  manipulation class preserved: %s\n',passFail(PATTERN_PASS));
assert(PATTERN_PASS, ...
    'FINAL spatial manipulation no longer matches the frozen PILOT design.');

%% 7. SUMMARIZE ARRANGEMENT CONDITIONS
Summary = struct();
for a=1:nArr
    Summary(a).name = arrNames{a};
    Summary(a).meanInitialMoran = mean(E.initialMoran(:,a),'omitnan');
    Summary(a).sdInitialMoran = std(E.initialMoran(:,a),'omitnan');
    Summary(a).meanAUTC_NV = mean(E.AUTC_NV(:,a));
    Summary(a).seAUTC_NV = std(E.AUTC_NV(:,a))/sqrt(N_REP);
    Summary(a).meanFinalNV = mean(E.finalNV(:,a));
    Summary(a).Pext = mean(E.extinctV(:,a));
    [Summary(a).PextLo,Summary(a).PextHi] = wilsonCI(sum(E.extinctV(:,a)),N_REP,1.96);
    Summary(a).meanRough120 = mean(E.rough120(:,a));
    Summary(a).meanPeakRough = mean(E.peakRough(:,a));
    Summary(a).meanCavity120 = mean(E.cavity120(:,a));
    Summary(a).meanMaxCavity = mean(E.maxCavity(:,a));
    Summary(a).meanLowFracFinal = mean(E.lowFracFinal(:,a),'omitnan');
    Summary(a).meanFinalMoran = mean(E.finalMoran(:,a),'omitnan');
end

AUTCrelHom = E.AUTC_NV ./ E.AUTC_NV(:,1);

%% 8. PAIRED EFFECTS
A_cluster_vs_random = pairedMeanDiff(E.AUTC_NV(:,4),E.AUTC_NV(:,3),N_BOOT,SEED_BOOT+1);
A_disperse_vs_random = pairedMeanDiff(E.AUTC_NV(:,2),E.AUTC_NV(:,3),N_BOOT,SEED_BOOT+2);
A_random_vs_hom = pairedMeanDiff(E.AUTC_NV(:,3),E.AUTC_NV(:,1),N_BOOT,SEED_BOOT+3);

R_cluster_vs_random = pairedMeanDiff(E.rough120(:,4),E.rough120(:,3),N_BOOT,SEED_BOOT+11);
R_disperse_vs_random = pairedMeanDiff(E.rough120(:,2),E.rough120(:,3),N_BOOT,SEED_BOOT+12);
C_cluster_vs_random = pairedMeanDiff(E.cavity120(:,4),E.cavity120(:,3),N_BOOT,SEED_BOOT+21);
C_disperse_vs_random = pairedMeanDiff(E.cavity120(:,2),E.cavity120(:,3),N_BOOT,SEED_BOOT+22);

xI = E.initialMoran(:,2:4); xI=xI(:);
yA = AUTCrelHom(:,2:4); yA=yA(:);
valid=isfinite(xI)&isfinite(yA);
[rPearson,slopeOLS,interceptOLS] = simpleLinearAssociation(xI(valid),yA(valid));
slopeCI = bootstrapSlopeCI(xI(valid),yA(valid),N_BOOT,SEED_BOOT+40);

%% 9. CONSOLE SUMMARY
fprintf('\n================ Q2 FINAL SUMMARY =================\n');
for a=1:nArr
    fprintf(['%-11s | Moran=%+.3f | AUTC_NV=%.1f | finalNV=%.1f | ', ...
        'Pext=%.3f | rough120=%.3f | cavity120=%.4f | lowFracFinal=%.3f\n'], ...
        Summary(a).name,Summary(a).meanInitialMoran,Summary(a).meanAUTC_NV, ...
        Summary(a).meanFinalNV,Summary(a).Pext,Summary(a).meanRough120, ...
        Summary(a).meanCavity120,Summary(a).meanLowFracFinal);
end

fprintf('\nPrimary paired effects (first-minus-second):\n');
fprintf('AUTC clustered-random  = %+g [%.3g, %.3g]\n', ...
    A_cluster_vs_random.point,A_cluster_vs_random.ci(1),A_cluster_vs_random.ci(2));
fprintf('AUTC dispersed-random  = %+g [%.3g, %.3g]\n', ...
    A_disperse_vs_random.point,A_disperse_vs_random.ci(1),A_disperse_vs_random.ci(2));
fprintf('AUTC random-homogeneous= %+g [%.3g, %.3g]\n', ...
    A_random_vs_hom.point,A_random_vs_hom.ci(1),A_random_vs_hom.ci(2));
fprintf('rough120 clustered-random = %+g [%.3g, %.3g]\n', ...
    R_cluster_vs_random.point,R_cluster_vs_random.ci(1),R_cluster_vs_random.ci(2));
fprintf('rough120 dispersed-random = %+g [%.3g, %.3g]\n', ...
    R_disperse_vs_random.point,R_disperse_vs_random.ci(1),R_disperse_vs_random.ci(2));
fprintf('cavity120 clustered-random= %+g [%.3g, %.3g]\n', ...
    C_cluster_vs_random.point,C_cluster_vs_random.ci(1),C_cluster_vs_random.ci(2));
fprintf('cavity120 dispersed-random= %+g [%.3g, %.3g]\n', ...
    C_disperse_vs_random.point,C_disperse_vs_random.ci(1),C_disperse_vs_random.ci(2));
fprintf('Moran I -> AUTC/homogeneous: r=%+.3f, slope=%+.4f [%.4f, %.4f]\n', ...
    rPearson,slopeOLS,slopeCI(1),slopeCI(2));

PHENOMENON_PASS = PATTERN_PASS && A_cluster_vs_random.ci(2) < 0;
fprintf('FINAL spatial-arrangement phenomenon: %s\n',passFail(PHENOMENON_PASS));
fprintf('===================================================\n');

%% 10. TABLES
Name = string(arrNames(:));
MeanInitialMoran = arrayfun(@(s)s.meanInitialMoran,Summary).';
SDInitialMoran = arrayfun(@(s)s.sdInitialMoran,Summary).';
MeanAUTC_NV = arrayfun(@(s)s.meanAUTC_NV,Summary).';
MeanFinalNV = arrayfun(@(s)s.meanFinalNV,Summary).';
Pext = arrayfun(@(s)s.Pext,Summary).';
PextLo = arrayfun(@(s)s.PextLo,Summary).';
PextHi = arrayfun(@(s)s.PextHi,Summary).';
MeanRough120 = arrayfun(@(s)s.meanRough120,Summary).';
MeanPeakRough = arrayfun(@(s)s.meanPeakRough,Summary).';
MeanCavity120 = arrayfun(@(s)s.meanCavity120,Summary).';
MeanMaxCavity = arrayfun(@(s)s.meanMaxCavity,Summary).';
MeanLowFracFinal = arrayfun(@(s)s.meanLowFracFinal,Summary).';
MeanFinalMoran = arrayfun(@(s)s.meanFinalMoran,Summary).';

Tsummary = table(Name,MeanInitialMoran,SDInitialMoran,MeanAUTC_NV,MeanFinalNV, ...
    Pext,PextLo,PextHi,MeanRough120,MeanPeakRough,MeanCavity120,MeanMaxCavity, ...
    MeanLowFracFinal,MeanFinalMoran);

Comparison = ["Clustered - Random";"Dispersed - Random";"Random - Homogeneous"; ...
              "Clustered - Random";"Dispersed - Random"; ...
              "Clustered - Random";"Dispersed - Random"];
Metric = ["AUTC_NV";"AUTC_NV";"AUTC_NV"; ...
          "Roughness_120h";"Roughness_120h";"Cavity_120h";"Cavity_120h"];
Effect = [A_cluster_vs_random.point;A_disperse_vs_random.point;A_random_vs_hom.point; ...
          R_cluster_vs_random.point;R_disperse_vs_random.point; ...
          C_cluster_vs_random.point;C_disperse_vs_random.point];
CI95_lo = [A_cluster_vs_random.ci(1);A_disperse_vs_random.ci(1);A_random_vs_hom.ci(1); ...
           R_cluster_vs_random.ci(1);R_disperse_vs_random.ci(1); ...
           C_cluster_vs_random.ci(1);C_disperse_vs_random.ci(1)];
CI95_hi = [A_cluster_vs_random.ci(2);A_disperse_vs_random.ci(2);A_random_vs_hom.ci(2); ...
           R_cluster_vs_random.ci(2);R_disperse_vs_random.ci(2); ...
           C_cluster_vs_random.ci(2);C_disperse_vs_random.ci(2)];
Tpaired = table(Comparison,Metric,Effect,CI95_lo,CI95_hi);

%% 11. SAVE DATA
outDir = fullfile(scriptDir,'Macro_Q2_SpatialArrangement_FINAL_OUT');
if ~exist(outDir,'dir'), mkdir(outDir); end
outMAT = fullfile(outDir,'Macro_Q2_SpatialArrangement_FINAL.mat');
outCSV = fullfile(outDir,'Macro_Q2_SpatialArrangement_Summary.csv');
outPairCSV = fullfile(outDir,'Macro_Q2_SpatialArrangement_PairedEffects.csv');
outPNG = fullfile(outDir,'Macro_Q2_SpatialArrangement_FINAL.png');
outPDF = fullfile(outDir,'Macro_Q2_SpatialArrangement_FINAL.pdf');
outPatterns = fullfile(outDir,'Macro_Q2_SpatialArrangement_PatternExamples_FINAL.png');

writetable(Tsummary,outCSV);
writetable(Tpaired,outPairCSV);

Q2=struct();
Q2.version='MACRO_Q2_SPATIAL_ARRANGEMENT_FINAL_v1';
Q2.runMode=RUN_MODE;
Q2.source.q1File=q1File;
Q2.settings=struct('L',L,'R0',R0,'Tend_h',TEND_H,'record_dt_h',RECORD_DT_H, ...
    'r_h',r_h,'mu_h',mu_h,'chiMu_ref',chiMuRef,'pHigh',P_HIGH, ...
    'nRep',N_REP,'nSwapAttempts',N_SWAP_ATTEMPTS);
Q2.q1Boundary=struct('chiArrest',chiArrest,'chiP50',chiP50,'chiP90',chiP90);
Q2.hazard=struct('chiLow',chiLow,'chiHigh',chiHigh,'chiBar',chiBar, ...
    'lambdaLow_h',lambdaLow_h,'lambdaHigh_h',lambdaHigh_h,'lambdaBar_h',lambdaBar_h);
Q2.initialCounts=struct('N0',N0,'nLow',nLow,'nHigh',nHigh);
Q2.arrangementNames=arrNames;
Q2.ensemble=E;
Q2.summary=Tsummary;
Q2.pairedEffects=Tpaired;
Q2.moranAUTCAssociation=struct('r',rPearson,'slope',slopeOLS,'slopeCI95',slopeCI);
Q2.patternExample=PATTERN_EXAMPLE;
Q2.finalDecision=struct('patternPASS',PATTERN_PASS,'phenomenonPASS',PHENOMENON_PASS);
save(outMAT,'Q2','-v7.3');

%% 12. PILOT SUMMARY FIGURE
cHom=hex2rgb('#4A4A4A');
cDisp=hex2rgb('#6AA6A1');
cRand=hex2rgb('#D3A536');
cClust=hex2rgb('#B4534B');
cc=[cHom;cDisp;cRand;cClust];

fig=figure('Color','w','Position',[70 70 1350 800]);
tl=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');

ax=nexttile(tl,1); hold(ax,'on');
x=1:nArr;
m=arrayfun(@(s)s.meanInitialMoran,Summary);
sd=arrayfun(@(s)s.sdInitialMoran,Summary);
errorbar(ax,x(2:4),m(2:4),sd(2:4),'o','LineWidth',1.3,'MarkerSize',6, ...
    'CapSize',6,'Color',[0.25 0.25 0.25]);
for a=2:4, scatter(ax,a,m(a),55,cc(a,:),'filled','MarkerEdgeColor','w'); end
yline(ax,0,':','Color',[0.6 0.6 0.6]);
xlim(ax,[1.5 4.5]); xticks(ax,2:4); xticklabels(ax,arrNames(2:4));
ylabel(ax,'Initial Moran''s I of \chi_F marks');
title(ax,'a  Spatial arrangement at fixed composition'); styleAxis(ax);

ax=nexttile(tl,2); hold(ax,'on');
for a=1:nArr
    plot(ax,E.t_h,mean(E.NVcurve(:,:,a),2),'LineWidth',1.7,'Color',cc(a,:), ...
        'DisplayName',arrNames{a});
end
xlabel(ax,'Time (h)'); ylabel(ax,'Mean viable cells, N_V');
title(ax,'b  Viable-burden dynamics');
legend(ax,'Location','best','Box','off','FontSize',7.5); styleAxis(ax);

ax=nexttile(tl,3); hold(ax,'on');
for a=1:nArr
    yy=AUTCrelHom(:,a); mm=mean(yy); ci=bootstrapMeanCI(yy,N_BOOT,SEED_BOOT+100+a);
    errorbar(ax,a,mm,mm-ci(1),ci(2)-mm,'o','Color',cc(a,:), ...
        'MarkerFaceColor',cc(a,:),'MarkerEdgeColor','w','LineWidth',1.2,'CapSize',6,'MarkerSize',7);
end
yline(ax,1,':','Color',[0.6 0.6 0.6]);
xlim(ax,[0.5 4.5]); xticks(ax,1:4); xticklabels(ax,arrNames);
ylabel(ax,'AUTC_{NV} / paired homogeneous'); title(ax,'c  Integrated viable burden'); styleAxis(ax);

ax=nexttile(tl,4); hold(ax,'on');
for a=1:nArr
    yy=E.rough120(:,a); mm=mean(yy); ci=bootstrapMeanCI(yy,N_BOOT,SEED_BOOT+200+a);
    errorbar(ax,a,mm,mm-ci(1),ci(2)-mm,'o','Color',cc(a,:), ...
        'MarkerFaceColor',cc(a,:),'MarkerEdgeColor','w','LineWidth',1.2,'CapSize',6,'MarkerSize',7);
end
xlim(ax,[0.5 4.5]); xticks(ax,1:4); xticklabels(ax,arrNames);
ylabel(ax,'Boundary roughness at 120 h (sites)'); title(ax,'d  Tumour-boundary roughness'); styleAxis(ax);

ax=nexttile(tl,5); hold(ax,'on');
for a=1:nArr
    yy=E.cavity120(:,a); mm=mean(yy); ci=bootstrapMeanCI(yy,N_BOOT,SEED_BOOT+300+a);
    errorbar(ax,a,mm,mm-ci(1),ci(2)-mm,'o','Color',cc(a,:), ...
        'MarkerFaceColor',cc(a,:),'MarkerEdgeColor','w','LineWidth',1.2,'CapSize',6,'MarkerSize',7);
end
xlim(ax,[0.5 4.5]); xticks(ax,1:4); xticklabels(ax,arrNames);
ylabel(ax,'Cavity fraction at 120 h'); title(ax,'e  Internal spatial disruption'); styleAxis(ax);

ax=nexttile(tl,6); hold(ax,'on');
for a=2:4
    scatter(ax,E.initialMoran(:,a),AUTCrelHom(:,a),24,cc(a,:),'filled', ...
        'MarkerFaceAlpha',0.65,'DisplayName',arrNames{a});
end
xx=linspace(min(xI(valid)),max(xI(valid)),100);
plot(ax,xx,interceptOLS+slopeOLS*xx,'-','Color',[0.15 0.15 0.15], ...
    'LineWidth',1.4,'HandleVisibility','off');
xlabel(ax,'Initial Moran''s I of \chi_F marks');
ylabel(ax,'AUTC_{NV} / paired homogeneous');
title(ax,sprintf('f  Arrangement-response relation  (r=%.2f)',rPearson));
legend(ax,'Location','best','Box','off','FontSize',7.5); styleAxis(ax);

sgtitle(tl,'Q2 final: spatial arrangement of ferroptotic susceptibility at fixed hazard distribution', ...
    'FontWeight','normal','FontSize',12);
exportgraphics(fig,outPNG,'Resolution',500);
exportgraphics(fig,outPDF,'ContentType','vector');

%% 13. PATTERN EXAMPLES
fig2=figure('Color','w','Position',[120 120 1050 340]);
tl2=tiledlayout(fig2,1,3,'TileSpacing','compact','Padding','compact');
patternFields={'dispersed','random','clustered'};
patternTitles={'Dispersed','Random','Clustered'};
for k=1:3
    ax=nexttile(tl2,k);
    showInitialCentralSlice(ax,PATTERN_EXAMPLE,patternFields{k},L,R0,ctr,chiLow,chiHigh);
    title(ax,sprintf('%s  (Moran I=%+.2f)',patternTitles{k}, ...
        PATTERN_EXAMPLE.metrics.(patternFields{k}).moranI),'FontWeight','normal');
end
sgtitle(tl2,'Matched initial ferroptotic-susceptibility patterns', ...
    'FontWeight','normal','FontSize',11);
exportgraphics(fig2,outPatterns,'Resolution',500);

fprintf('\nSaved Q2 FINAL outputs:\n  %s\n  %s\n  %s\n  %s\n  %s\n  %s\n', ...
    outMAT,outCSV,outPairCSV,outPNG,outPDF,outPatterns);
fprintf(['\nFINAL STATUS: spatial-arrangement phenomenon confirmed on an independent ', ...
    'N=100 ensemble. Use this file as the source-of-truth anchor for the FINAL ', ...
    'mechanism replay and Q2 condition-map confirmation.\n']);

end

%% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================
function R=runOneRep(rep,L,R0,Tend,recordDt,snapTimes,r,mu, ...
    chiLow,chiHigh,pHigh,nSwapAttempts,patternSeedBase,simSeedBase)

P=buildMatchedPatterns(L,R0,chiLow,chiHigh,pHigh, ...
    patternSeedBase+rep,nSwapAttempts);

lambdaLow=chiLow*r; lambdaHigh=chiHigh*r;
N0=P.N0; nHigh=P.nHigh; nLow=P.nLow;
lambdaBar=(nLow*lambdaLow+nHigh*lambdaHigh)/N0;

mapsChi={P.mapHomogeneous,P.mapDispersed,P.mapRandom,P.mapClustered};
R=struct();
R.initialMoran=[NaN,P.metrics.dispersed.moranI,P.metrics.random.moranI,P.metrics.clustered.moranI];
R.sameEdgeFraction=[NaN,P.metrics.dispersed.sameEdgeFraction, ...
    P.metrics.random.sameEdgeFraction,P.metrics.clustered.sameEdgeFraction];

for a=1:4
    cfg=struct();
    cfg.L=L; cfg.R0=R0; cfg.Tend_h=Tend; cfg.record_dt_h=recordDt;
    cfg.r_h=r; cfg.mu_h=mu; cfg.seed=simSeedBase+rep;
    cfg.lambda_mode='cell_mark'; cfg.lambdaF_h=lambdaBar;
    cfg.lambdaF_initial_map_h=mapsChi{a}*r;
    cfg.daughter_lambda_rule='inherit';
    cfg.lambdaF_pool_h=[]; cfg.lambdaF_site_h=[];
    cfg.snapshot_times_h=snapTimes; cfg.return_final_state=true;
    cfg.compute_cavity=true; cfg.store_event_log=false; cfg.verbose=false;

    O=CTMC_EventDrivenSimulation(cfg);
    R.out(a)=reduceQ2Output(O,chiLow,chiHigh); %#ok<AGROW>
end
end

function Z=reduceQ2Output(O,chiLow,chiHigh)
Z=struct();
Z.t_h=double(O.t_h(:)); Z.NV=double(O.NV(:)); Z.NT=double(O.NT(:));
Z.finalNV=double(O.NV(end)); Z.finalNF=double(O.NF(end)); Z.finalNT=double(O.NT(end));
Z.AUTC_NV=double(O.AUTC_NV); Z.AUTC_NT=double(O.AUTC_NT);
Z.extinctV=isfinite(O.viableExtinctionTime_h) && O.viableExtinctionTime_h<=O.cfg.Tend_h;
Z.peakRough=max(double(O.roughness_sites)); Z.maxCavity=max(double(O.cavity_fraction));
Z.finalRmax=double(O.Rmax_sites(end));
[~,i120]=min(abs(double(O.t_h(:))-120));
Z.rough120=double(O.roughness_sites(i120)); Z.cavity120=double(O.cavity_fraction(i120));
Z.edgeTouch=any(O.edgeOcc); Z.bookkeepingPASS=logical(O.audit.bookkeepingPASS);
Z.lowFracFinal=NaN; Z.finalMoran=NaN;
if ~isempty(O.final_state) && ~isempty(O.final_lambda_cell_h)
    Vmask=(O.final_state==1);
    if any(Vmask(:))
        chiMap=double(O.final_lambda_cell_h)/double(O.cfg.r_h);
        vals=chiMap(Vmask); mid=(chiLow+chiHigh)/2;
        Z.lowFracFinal=mean(vals<mid);
        Z.finalMoran=moranIGrid(Vmask,chiMap);
    end
end
end

function E=stackReplicates(tmp,nArr)
nRep=numel(tmp); t=tmp{1}.out(1).t_h; nT=numel(t);
E=struct(); E.t_h=t;
E.initialMoran=nan(nRep,nArr); E.sameEdgeFraction=nan(nRep,nArr);
fields={'finalNV','finalNF','finalNT','AUTC_NV','AUTC_NT','rough120','peakRough', ...
    'cavity120','maxCavity','finalRmax'};
for f=1:numel(fields), E.(fields{f})=zeros(nRep,nArr); end
E.extinctV=false(nRep,nArr); E.lowFracFinal=nan(nRep,nArr); E.finalMoran=nan(nRep,nArr);
E.edgeTouch=false(nRep,nArr); E.bookkeepingPASS=false(nRep,nArr);
E.NVcurve=zeros(nT,nRep,nArr); E.NTcurve=zeros(nT,nRep,nArr);
for rep=1:nRep
    E.initialMoran(rep,:)=tmp{rep}.initialMoran;
    E.sameEdgeFraction(rep,:)=tmp{rep}.sameEdgeFraction;
    for a=1:nArr
        z=tmp{rep}.out(a);
        for f=1:numel(fields), E.(fields{f})(rep,a)=z.(fields{f}); end
        E.extinctV(rep,a)=z.extinctV; E.lowFracFinal(rep,a)=z.lowFracFinal;
        E.finalMoran(rep,a)=z.finalMoran; E.edgeTouch(rep,a)=z.edgeTouch;
        E.bookkeepingPASS(rep,a)=z.bookkeepingPASS;
        E.NVcurve(:,rep,a)=z.NV; E.NTcurve(:,rep,a)=z.NT;
    end
end
end

function P=buildMatchedPatterns(L,R0,chiLow,chiHigh,pHigh,seed,nSwapAttempts)
rng(seed,'twister');
ctr=[(L+1)/2 (L+1)/2 (L+1)/2];
[xg,yg,zg]=ndgrid(1:L,1:L,1:L);
r2=(xg-ctr(1)).^2+(yg-ctr(2)).^2+(zg-ctr(3)).^2;
mask=r2<=R0^2; idx=find(mask); N=numel(idx);
[ii,jj,kk]=ind2sub([L L L],idx);
rad=sqrt((ii-ctr(1)).^2+(jj-ctr(2)).^2+(kk-ctr(3)).^2);
shell=floor(rad+1e-12); uShell=unique(shell(:)).';

nHighTotal=round(pHigh*N);
nByShell=zeros(size(uShell)); rawHigh=zeros(size(uShell));
for s=1:numel(uShell)
    nByShell(s)=nnz(shell==uShell(s)); rawHigh(s)=pHigh*nByShell(s);
end
nHighShell=floor(rawHigh); remain=nHighTotal-sum(nHighShell);
[~,ordFrac]=sort(rawHigh-floor(rawHigh),'descend');
for k=1:remain, nHighShell(ordFrac(k))=nHighShell(ordFrac(k))+1; end

labelRandom=-ones(N,1);
for s=1:numel(uShell)
    loc=find(shell==uShell(s)); perm=loc(randperm(numel(loc))); nh=nHighShell(s);
    if nh>0, labelRandom(perm(1:nh))=+1; end
end

[e1,e2,incEdges]=initialGraph(idx,L);
labelClustered=optimizeArrangement(labelRandom,shell,e1,e2,incEdges,+1,nSwapAttempts);
labelDispersed=optimizeArrangement(labelRandom,shell,e1,e2,incEdges,-1,nSwapAttempts);

nLow=N-nHighTotal; chiBar=(nLow*chiLow+nHighTotal*chiHigh)/N;
mapHom=chiBar*ones(L,L,L); mapRand=mapHom; mapDisp=mapHom; mapClust=mapHom;
mapRand(idx(labelRandom<0))=chiLow; mapRand(idx(labelRandom>0))=chiHigh;
mapDisp(idx(labelDispersed<0))=chiLow; mapDisp(idx(labelDispersed>0))=chiHigh;
mapClust(idx(labelClustered<0))=chiLow; mapClust(idx(labelClustered>0))=chiHigh;

P=struct('N0',N,'nHigh',nHighTotal,'nLow',nLow,'initialMask',mask,'initialIdx',idx, ...
    'radialShell',shell,'mapHomogeneous',mapHom,'mapRandom',mapRand, ...
    'mapDispersed',mapDisp,'mapClustered',mapClust,'labelRandom',labelRandom, ...
    'labelDispersed',labelDispersed,'labelClustered',labelClustered);
P.metrics.random=patternMetrics(labelRandom,e1,e2);
P.metrics.dispersed=patternMetrics(labelDispersed,e1,e2);
P.metrics.clustered=patternMetrics(labelClustered,e1,e2);

assert(sum(labelRandom>0)==sum(labelDispersed>0) && sum(labelRandom>0)==sum(labelClustered>0));
for s=1:numel(uShell)
    loc=(shell==uShell(s));
    assert(nnz(labelRandom(loc)>0)==nnz(labelDispersed(loc)>0));
    assert(nnz(labelRandom(loc)>0)==nnz(labelClustered(loc)>0));
end
end

function label=optimizeArrangement(label0,shell,e1,e2,incEdges,direction,nAttempt)
label=label0(:); uShell=unique(shell(:)).'; eligible=false(size(uShell));
for s=1:numel(uShell)
    loc=find(shell==uShell(s)); eligible(s)=any(label(loc)>0)&&any(label(loc)<0);
end
eligibleShell=uShell(eligible); if isempty(eligibleShell), return; end
for it=1:nAttempt
    sh=eligibleShell(randi(numel(eligibleShell))); loc=find(shell==sh);
    pos=loc(label(loc)>0); neg=loc(label(loc)<0);
    if isempty(pos)||isempty(neg), continue; end
    a=pos(randi(numel(pos))); b=neg(randi(numel(neg)));
    affected=unique([incEdges{a};incEdges{b}]);
    oldScore=sum(label(e1(affected)).*label(e2(affected)));
    label(a)=-label(a); label(b)=-label(b);
    newScore=sum(label(e1(affected)).*label(e2(affected)));
    delta=newScore-oldScore; accept=(direction*delta>0);
    if delta==0, accept=(rand<0.10); end
    if ~accept, label(a)=-label(a); label(b)=-label(b); end
end
end

function [e1,e2,incEdges]=initialGraph(globalIdx,L)
N=numel(globalIdx); localOfGlobal=zeros(L^3,1,'uint32'); localOfGlobal(globalIdx)=uint32(1:N);
e1=zeros(3*N,1,'uint32'); e2=zeros(3*N,1,'uint32'); ne=0;
for a=1:N
    g=globalIdx(a); [i,j,k]=ind2sub([L L L],g);
    if i<L
        b=localOfGlobal(sub2ind([L L L],i+1,j,k)); if b>0, ne=ne+1; e1(ne)=a; e2(ne)=b; end
    end
    if j<L
        b=localOfGlobal(sub2ind([L L L],i,j+1,k)); if b>0, ne=ne+1; e1(ne)=a; e2(ne)=b; end
    end
    if k<L
        b=localOfGlobal(sub2ind([L L L],i,j,k+1)); if b>0, ne=ne+1; e1(ne)=a; e2(ne)=b; end
    end
end
e1=double(e1(1:ne)); e2=double(e2(1:ne)); incEdges=cell(N,1);
for e=1:ne
    incEdges{e1(e)}(end+1,1)=e; %#ok<AGROW>
    incEdges{e2(e)}(end+1,1)=e; %#ok<AGROW>
end
end

function M=patternMetrics(label,e1,e2)
label=double(label(:)); N=numel(label); E=numel(e1); z=label-mean(label); den=sum(z.^2);
if den<=eps||E==0, I=NaN; else, I=(N/E)*sum(z(e1).*z(e2))/den; end
M=struct('moranI',I,'sameEdgeFraction',mean(label(e1)==label(e2)),'nEdges',E);
end

function I=moranIGrid(mask,valueMap)
L=size(mask,1); idx=find(mask); N=numel(idx);
if N<2, I=NaN; return; end
vals=double(valueMap(idx)); if all(abs(vals-vals(1))<1e-14), I=NaN; return; end
localOfGlobal=zeros(L^3,1,'uint32'); localOfGlobal(idx)=uint32(1:N);
e1=[]; e2=[];
for a=1:N
    g=idx(a); [i,j,k]=ind2sub([L L L],g);
    if i<L
        b=localOfGlobal(sub2ind([L L L],i+1,j,k)); if b>0, e1(end+1,1)=a; e2(end+1,1)=double(b); end %#ok<AGROW>
    end
    if j<L
        b=localOfGlobal(sub2ind([L L L],i,j+1,k)); if b>0, e1(end+1,1)=a; e2(end+1,1)=double(b); end %#ok<AGROW>
    end
    if k<L
        b=localOfGlobal(sub2ind([L L L],i,j,k+1)); if b>0, e1(end+1,1)=a; e2(end+1,1)=double(b); end %#ok<AGROW>
    end
end
if isempty(e1), I=NaN; return; end
z=vals-mean(vals); den=sum(z.^2);
if den<=eps, I=NaN; else, I=(N/numel(e1))*sum(z(e1).*z(e2))/den; end
end

function A=pairedMeanDiff(x,y,nBoot,seed)
d=double(x(:))-double(y(:)); A.point=mean(d); A.ci=bootstrapMeanCI(d,nBoot,seed);
end

function ci=bootstrapMeanCI(x,nBoot,seed)
x=double(x(:)); x=x(isfinite(x)); if isempty(x), ci=[NaN NaN]; return; end
rng(seed,'twister'); n=numel(x); boot=zeros(nBoot,1);
for b=1:nBoot, ii=randi(n,n,1); boot(b)=mean(x(ii)); end
boot=sort(boot); lo=max(1,round(0.025*nBoot)); hi=min(nBoot,round(0.975*nBoot)); ci=[boot(lo) boot(hi)];
end

function ci=bootstrapSlopeCI(x,y,nBoot,seed)
x=double(x(:)); y=double(y(:)); mask=isfinite(x)&isfinite(y); x=x(mask); y=y(mask);
rng(seed,'twister'); n=numel(x); slopes=zeros(nBoot,1);
for b=1:nBoot, ii=randi(n,n,1); p=polyfit(x(ii),y(ii),1); slopes(b)=p(1); end
slopes=sort(slopes); lo=max(1,round(0.025*nBoot)); hi=min(nBoot,round(0.975*nBoot)); ci=[slopes(lo) slopes(hi)];
end

function [r,slope,intercept]=simpleLinearAssociation(x,y)
C=corrcoef(double(x(:)),double(y(:))); if numel(C)>=4, r=C(1,2); else, r=NaN; end
p=polyfit(double(x(:)),double(y(:)),1); slope=p(1); intercept=p(2);
end

function [lo,hi]=wilsonCI(k,n,z)
if n==0, lo=NaN; hi=NaN; return; end
ph=k/n; den=1+z^2/n; cen=(ph+z^2/(2*n))/den;
half=z*sqrt(ph*(1-ph)/n+z^2/(4*n^2))/den;
lo=max(0,cen-half); hi=min(1,cen+half);
end

function B=orientBoundaryMatrix(B,nMu,nLevel,label)
if isequal(size(B),[nMu nLevel]), return; elseif isequal(size(B),[nLevel nMu]), B=B.'; else, error('%s has unexpected size.',label); end
end

function showInitialCentralSlice(ax,P,fieldName,L,R0,ctr,chiLow,chiHigh)
switch fieldName
    case 'dispersed', M=P.mapDispersed;
    case 'random', M=P.mapRandom;
    case 'clustered', M=P.mapClustered;
    otherwise, error('Unknown fieldName.');
end
mask=P.initialMask; z0=round(ctr(3)); slice=M(:,:,z0); ms=mask(:,:,z0); slice(~ms)=NaN;
half=R0+1; ix=(round(ctr(1))-half):(round(ctr(1))+half); iy=(round(ctr(2))-half):(round(ctr(2))+half);
imagesc(ax,iy-ctr(2),ix-ctr(1),slice(ix,iy)); axis(ax,'image'); set(ax,'YDir','normal');
clim(ax,[chiLow chiHigh]); colormap(ax,[hex2rgb('#79A9A4');hex2rgb('#C75B50')]);
xlabel(ax,'y (sites)'); ylabel(ax,'x (sites)'); box(ax,'off');
end

function styleAxis(ax)
ax.FontName='Arial'; ax.FontSize=8.5; ax.LineWidth=0.8; ax.TickDir='out';
ax.TickLength=[0.015 0.015]; ax.Box='off'; ax.Layer='top'; grid(ax,'off');
end

function s=passFail(tf)
if tf
    s='PASS';
else
    s='FAIL';
end
end

function rgb=hex2rgb(h)
h=char(h); if h(1)=='#', h=h(2:end); end
rgb=[hex2dec(h(1:2)) hex2dec(h(3:4)) hex2dec(h(5:6))]/255;
end
