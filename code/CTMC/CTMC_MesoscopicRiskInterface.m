function CTMC_MesoscopicRiskInterface()
% =========================================================================
% CTMC_MesoscopicRiskInterface
%
% FINAL mesoscopic -> macroscopic ferroptotic-risk interface
%
% INPUT (actual files from the completed mesoscopic study)
%   1) Fig4D_CellStateDependentFateRewiring_FINAL_LAYOUTFIX_data(1).mat
%      - primary interface source
%      - 12 tau_L states x 4 microscopic backgrounds
%      - contains 0--60 min FPT samples and summary statistics
%
%   2) Fig4B_fateRoutes_data(3).mat
%      - reference anchor at tau_L = 10 min
%
% SCIENTIFIC INTERFACE
%   Mesoscopic 0--T_I commitment statistics are NOT used as lifetime FPT
%   deadlines in the macro model. Instead, for each mesoscopic state z,
%
%       RMST_z(T_I)
%          -> unique RMST-equivalent hazard lambda_F,z^eq
%
%   defined by
%
%       RMST_z(T_I) = [1-exp(-lambda*T_I)] / lambda .
%
%   T_I = 60 min is a coarse-graining window, not cell lifetime and not
%   macro experiment duration.  lambda_F^eq is stored in h^-1.
%
% TIME-SCALE HANDLING
%   A) Primary interface window is locked at 60 min, consistent with the
%      completed micro/meso analyses.
%   B) 30/45/60 min window audit is computed from the SAME FPT samples.
%      This is a robustness audit of the coarse-graining choice; it does not
%      redefine the primary 60-min interface.
%   C) If a locked macro r file is available, the script additionally audits
%      scale separation:
%
%           epsilon_r = T_I * r
%
%      and reports lambda_F^eq/r.  If r is not yet locked, this audit is
%      explicitly deferred rather than filled with a placeholder.
%   D) Continuous-time hazards are converted to any macro time step dt by
%
%           p_F(dt) = 1-exp(-lambda_F^eq*dt)
%
%      so the interface itself does not depend on the macro numerical dt.
%
% IMPORTANT
%   - The 12,000 FPT realizations within a single mesoscopic state are
%     stochastic outcomes, NOT 12,000 persistent cell-specific hazard marks.
%   - Stable macro heterogeneity should be constructed from mixtures of
%     mesoscopic states (e.g. different tau_L / repair states), each carrying
%     its corresponding lambda_F^eq.
%   - No 28/36/60-h FPT-tail extension files are used.
%
% OUTPUT
%   Meso_to_Macro_Risk_Interface_FINAL_OUT/
%       Meso_to_Macro_Risk_Interface_FINAL.mat
%       Meso_to_Macro_Risk_Interface_48states.csv
%       Meso_to_Macro_WindowAudit.csv
%       Meso_to_Macro_Risk_Interface_FINAL.png
%       Meso_to_Macro_Risk_Interface_console.txt
% =========================================================================

clc;
close all;

%% ========================================================================
% 0. USER SETTINGS
% =========================================================================
PRIMARY_T_MIN = 60;
AUDIT_T_MIN   = [30 45 60];

% Optional macro numerical steps for exact hazard -> step-probability table.
MACRO_DT_H = [0.125 0.25 0.5];

% Optional time-scale separation threshold.
% This is an audit criterion, not a fitted biological parameter.
EPS_SCALE_TARGET = 0.10;

% If Macro_r_LOCKED.mat exists, the code reads r automatically.
% Otherwise r-dependent diagnostics are deferred.
R_LOCK_CANDIDATES = { ...
    fullfile('Macro_r_Calibration_OUT','Macro_r_LOCKED.mat'), ...
    'Macro_r_LOCKED.mat'};

% Numerical tolerances.
ANCHOR_TOL_RMST = 1e-10;
ANCHOR_TOL_P    = 1e-12;
SOURCE_TOL_RMST = 1e-9;
SOURCE_TOL_P    = 1e-12;

scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir = pwd; end

mesoDir = fullfile('..','LPO-FPT');
fig4DFile = findFirstExisting(scriptDir,{ ...
    fullfile(mesoDir,'Fig4D_FINAL_outputs', ...
        'Fig4D_CellStateDependentFateRewiring_FINAL_LAYOUTFIX_data.mat'), ...
    'Fig4D_CellStateDependentFateRewiring_FINAL_LAYOUTFIX_data(1).mat', ...
    'Fig4D_CellStateDependentFateRewiring_FINAL_LAYOUTFIX_data.mat'});

fig4BFile = findFirstExisting(scriptDir,{ ...
    fullfile(mesoDir,'Fig4B_fate_routes_final','Fig4B_fateRoutes_data.mat'), ...
    'Fig4B_fateRoutes_data(3).mat', ...
    'Fig4B_fateRoutes_data.mat'});

assert(~isempty(fig4DFile),'Cannot find Fig4D mesoscopic data file.');
assert(~isempty(fig4BFile),'Cannot find Fig4B anchor data file.');

outDir = fullfile(scriptDir,'Meso_to_Macro_Risk_Interface_FINAL_OUT');
if ~exist(outDir,'dir'), mkdir(outDir); end

outMAT = fullfile(outDir,'Meso_to_Macro_Risk_Interface_FINAL.mat');
outCSV = fullfile(outDir,'Meso_to_Macro_Risk_Interface_48states.csv');
outAuditCSV = fullfile(outDir,'Meso_to_Macro_WindowAudit.csv');
outPNG = fullfile(outDir,'Meso_to_Macro_Risk_Interface_FINAL.png');
outTXT = fullfile(outDir,'Meso_to_Macro_Risk_Interface_console.txt');

diary(outTXT);
diaryCleanup = onCleanup(@() diary('off')); %#ok<NASGU>

fprintf('\n============================================================\n');
fprintf('MESOSCOPIC -> MACROSCOPIC FERROPTOTIC-RISK INTERFACE\n');
fprintf('Primary coarse-graining window T_I = %.1f min = %.3f h\n', ...
    PRIMARY_T_MIN,PRIMARY_T_MIN/60);
fprintf('Fig4D source : %s\n',fig4DFile);
fprintf('Fig4B anchor : %s\n',fig4BFile);
fprintf('============================================================\n');

%% ========================================================================
% 1. LOAD AND VERIFY ACTUAL MESOSCOPIC FILES
% =========================================================================
D = load(fig4DFile);
B = load(fig4BFile);

requiredD = {'R','RMST','Pcommit','condNames','condLabels','meta'};
for k=1:numel(requiredD)
    assert(isfield(D,requiredD{k}),'Fig4D lacks variable %s.',requiredD{k});
end
requiredB = {'R','condNames','condLabels','meta'};
for k=1:numel(requiredB)
    assert(isfield(B,requiredB{k}),'Fig4B lacks variable %s.',requiredB{k});
end

tau = double(D.meta.tauScan_min(:));
condNames = cellstr(string(D.condNames(:)));
condLabels = cellstr(string(D.condLabels(:)));

nTau = numel(tau);
nCond = numel(condNames);

assert(nTau==12,'Expected 12 tau_L scan states.');
assert(nCond==4,'Expected 4 mesoscopic background conditions.');
assert(abs(double(D.meta.Tobs_min)-PRIMARY_T_MIN)<1e-12, ...
    'Fig4D Tobs differs from the locked 60-min interface window.');

% ------------------------------------------------------------------------
% IMPORTANT ORIENTATION NORMALIZATION
% ------------------------------------------------------------------------
% The current Fig4D FINAL file is stored as condition x tau in MATLAB
% (4 x 12), although HDF5 viewers often display the reversed 12 x 4 shape.
% Normalize all state arrays internally to tau x condition (12 x 4).
RMST_SRC = orientTauByCondition(double(D.RMST),nTau,nCond,'Fig4D RMST');
PCOMMIT_SRC = orientTauByCondition(double(D.Pcommit),nTau,nCond,'Fig4D Pcommit');
R_SRC = orientTauByConditionCell(D.R,nTau,nCond,'Fig4D R');

expectedCond = {'CU','GU','CF','GF'};
assert(isequal(condNames(:)',expectedCond), ...
    'Unexpected condition ordering in Fig4D.');

fprintf('[1] Fig4D verified: %d tau states x %d conditions = %d mesoscopic states.\n', ...
    nTau,nCond,nTau*nCond);
fprintf('    Original MATLAB sizes: RMST=%s, Pcommit=%s, R=%s.\n', ...
    mat2str(size(D.RMST)),mat2str(size(D.Pcommit)),mat2str(size(D.R)));
fprintf('    Internal normalized size: %d x %d (tau x condition).\n',nTau,nCond);

%% ========================================================================
% 2. RECOMPUTE 60-min SUMMARY FROM RAW FPT_TS AND VERIFY SOURCE TABLE
% =========================================================================
RMST60_fromFPT = nan(nTau,nCond);
P60_fromFPT = nan(nTau,nCond);
Psi60 = nan(nTau,nCond);
lambdaRMST_h = nan(nTau,nCond);

nFPT = nan(nTau,nCond);

for j=1:nTau
    for c=1:nCond
        Rij = R_SRC{j,c};
        assert(isstruct(Rij) && isfield(Rij,'FPT_TS'), ...
            'Fig4D R{%d,%d} lacks FPT_TS.',j,c);

        tf = double(Rij.FPT_TS(:));
        assert(~any(isnan(tf)),'NaN FPT marks in tau=%g, %s.',tau(j),condNames{c});
        assert(all(tf(isfinite(tf))>=0), ...
            'Negative FPT mark in tau=%g, %s.',tau(j),condNames{c});

        nFPT(j,c) = numel(tf);
        RMST60_fromFPT(j,c) = empiricalRMST(tf,PRIMARY_T_MIN);
        P60_fromFPT(j,c) = mean(tf<=PRIMARY_T_MIN);
        Psi60(j,c) = 1 - RMST60_fromFPT(j,c)/PRIMARY_T_MIN;
        lambdaRMST_h(j,c) = rmstEquivalentHazard( ...
            RMST60_fromFPT(j,c)/60, PRIMARY_T_MIN/60);
    end
end

maxRMSTSourceErr = max(abs(RMST60_fromFPT(:)-RMST_SRC(:)));
maxPSourceErr = max(abs(P60_fromFPT(:)-PCOMMIT_SRC(:)));

assert(maxRMSTSourceErr<SOURCE_TOL_RMST, ...
    'Raw FPT-derived RMST does not reproduce Fig4D RMST.');
assert(maxPSourceErr<SOURCE_TOL_P, ...
    'Raw FPT-derived Pcommit does not reproduce Fig4D Pcommit.');

fprintf('[2] Raw FPT -> Fig4D summary reproduction: PASS.\n');
fprintf('    max |dRMST| = %.3e min; max |dPcommit| = %.3e.\n', ...
    maxRMSTSourceErr,maxPSourceErr);
fprintf('    FPT sample count range = %d--%d per state.\n', ...
    min(nFPT(:)),max(nFPT(:)));

%% ========================================================================
% 3. FIG4B REFERENCE-ANCHOR AUDIT AT tau_L = 10 min
% =========================================================================
[~,iRef] = min(abs(tau-double(D.meta.tauReference_min)));
assert(abs(tau(iRef)-10)<1e-12,'Expected tau_L reference = 10 min.');

BcondNames = cellstr(string(B.condNames(:)));
assert(isequal(BcondNames(:)',expectedCond), ...
    'Fig4B and Fig4D condition ordering differs.');

B_RMST = nan(1,nCond);
B_P = nan(1,nCond);
fptAnchorExact = true(1,nCond);

for c=1:nCond
    [B_RMST(c),B_P(c),tfB] = getFig4BCondition(B.R,c,nCond);
    tfD = double(R_SRC{iRef,c}.FPT_TS(:));

    if ~isequaln(tfB(:),tfD(:))
        % Numerical fallback; Inf masks must still match.
        sameSize = numel(tfB)==numel(tfD);
        sameInf = sameSize && isequal(isinf(tfB(:)),isinf(tfD(:)));
        if sameInf
            finiteMask = isfinite(tfB(:)) & isfinite(tfD(:));
            fptAnchorExact(c)=all(abs(tfB(finiteMask)-tfD(finiteMask))<1e-12);
        else
            fptAnchorExact(c)=false;
        end
    end
end

anchorRMSTErr = max(abs(B_RMST - RMST60_fromFPT(iRef,:)));
anchorPErr = max(abs(B_P - P60_fromFPT(iRef,:)));

anchorPASS = anchorRMSTErr<ANCHOR_TOL_RMST && ...
             anchorPErr<ANCHOR_TOL_P && all(fptAnchorExact);

assert(anchorPASS,'Fig4D tau=10 does not exactly reproduce Fig4B anchor.');

fprintf('[3] Fig4B -> Fig4D tau=10 anchor: PASS.\n');
fprintf('    max |dRMST| = %.3e min; max |dP| = %.3e; raw FPT match = %d.\n', ...
    anchorRMSTErr,anchorPErr,all(fptAnchorExact));

%% ========================================================================
% 4. TIME-WINDOW COARSE-GRAINING AUDIT: 30 / 45 / 60 min
%
% The primary interface stays at 60 min. Shorter windows only test whether
% the state ordering / relative pressure structure is grossly dependent on
% a single endpoint choice.
% =========================================================================
nWin = numel(AUDIT_T_MIN);
RMST_AUDIT_min = nan(nTau,nCond,nWin);
P_AUDIT = nan(nTau,nCond,nWin);
Psi_AUDIT = nan(nTau,nCond,nWin);
Lambda_AUDIT_h = nan(nTau,nCond,nWin);

for w=1:nWin
    Tw = AUDIT_T_MIN(w);

    for j=1:nTau
        for c=1:nCond
            tf = double(R_SRC{j,c}.FPT_TS(:));

            rmst = empiricalRMST(tf,Tw);
            pc = mean(tf<=Tw);

            RMST_AUDIT_min(j,c,w)=rmst;
            P_AUDIT(j,c,w)=pc;
            Psi_AUDIT(j,c,w)=1-rmst/Tw;
            Lambda_AUDIT_h(j,c,w)=rmstEquivalentHazard(rmst/60,Tw/60);
        end
    end
end

iPrimaryWin = find(abs(AUDIT_T_MIN-PRIMARY_T_MIN)<1e-12,1);
assert(~isempty(iPrimaryWin),'AUDIT_T_MIN must include PRIMARY_T_MIN.');
assert(max(abs(Lambda_AUDIT_h(:,:,iPrimaryWin)-lambdaRMST_h),[],'all')<1e-10);

% Rank concordance versus the primary 60-min interface.
rankCorrVs60 = nan(nWin,1);
orderAgreementVs60 = nan(nWin,1);

primaryFlat = lambdaRMST_h(:);
primaryRank = tiedRankLocal(primaryFlat);

for w=1:nWin
    x = Lambda_AUDIT_h(:,:,w);
    x = x(:);

    rw = tiedRankLocal(x);
    rankCorrVs60(w)=pearsonLocal(rw,primaryRank);

    orderAgreementVs60(w)=pairwiseOrderAgreement(x,primaryFlat);
end

fprintf('[4] Coarse-graining window audit completed.\n');
for w=1:nWin
    fprintf('    T=%4.0f min: rank corr vs 60 = %.4f; pairwise order agreement = %.4f.\n', ...
        AUDIT_T_MIN(w),rankCorrVs60(w),orderAgreementVs60(w));
end
fprintf('    NOTE: 60 min remains the predefined primary interface window.\n');

%% ========================================================================
% 5. OPTIONAL MACRO TIME-SCALE-SEPARATION AUDIT USING A LOCKED r
%
% No placeholder r is invented here. If no locked r exists, this part is
% explicitly left pending and can be rerun later without recomputing meso.
% =========================================================================
rMacro_h = NaN;
rSource = '';

for k=1:numel(R_LOCK_CANDIDATES)
    f = fullfile(scriptDir,R_LOCK_CANDIDATES{k});
    if isfile(f)
        Rlock = load(f);
        if isfield(Rlock,'LOCK') && isfield(Rlock.LOCK,'r_star_per_h')
            rMacro_h = double(Rlock.LOCK.r_star_per_h);
            rSource=f;
            break;
        end
    end
end

TI_h = PRIMARY_T_MIN/60;

if isfinite(rMacro_h) && rMacro_h>0
    epsilon_r = TI_h*rMacro_h;
    scaleSeparationPASS = epsilon_r <= EPS_SCALE_TARGET;
    lambdaOverR = lambdaRMST_h ./ rMacro_h;

    fprintf('[5] Locked-r time-scale audit available.\n');
    fprintf('    r = %.6g h^-1 from %s\n',rMacro_h,rSource);
    fprintf('    epsilon_r = T_I*r = %.4g -> %s (target <= %.3g).\n', ...
        epsilon_r,passfail(scaleSeparationPASS),EPS_SCALE_TARGET);
    fprintf('    lambda_F^eq/r range across 48 meso states = %.3g--%.3g.\n', ...
        min(lambdaOverR(:)),max(lambdaOverR(:)));
else
    epsilon_r = NaN;
    scaleSeparationPASS = NaN;
    lambdaOverR = nan(size(lambdaRMST_h));

    fprintf('[5] Locked-r time-scale audit DEFERRED.\n');
    fprintf('    No Macro_r_LOCKED.mat was found. No placeholder r was used.\n');
    fprintf(['    Re-run this interface after r calibration to obtain T_I*r ', ...
        'and lambda_F^eq/r diagnostics.\n']);
end

%% ========================================================================
% 6. EXACT CONTINUOUS-HAZARD -> MACRO STEP PROBABILITIES
%
% This resolves numerical time-step conversion without redefining lambda.
% =========================================================================
nDt = numel(MACRO_DT_H);
Pstep = nan(nTau,nCond,nDt);
for q=1:nDt
    Pstep(:,:,q) = 1-exp(-lambdaRMST_h*MACRO_DT_H(q));
end

fprintf('[6] Exact hazard -> macro-step probability table generated.\n');
for q=1:nDt
    fprintf('    dt=%.3f h: p_F range = %.4f--%.4f.\n', ...
        MACRO_DT_H(q),min(Pstep(:,:,q),[],'all'),max(Pstep(:,:,q),[],'all'));
end

%% ========================================================================
% 7. BUILD THE 48-STATE INTERFACE TABLE
% =========================================================================
nState = nTau*nCond;

StateID = strings(nState,1);
Condition = strings(nState,1);
ConditionLabel = strings(nState,1);
TauL_min = nan(nState,1);
RMST60_min = nan(nState,1);
Pcommit60 = nan(nState,1);
PsiF60 = nan(nState,1);
LambdaF_eq_h = nan(nState,1);
LambdaOverR = nan(nState,1);

ip=0;
for c=1:nCond
    for j=1:nTau
        ip=ip+1;

        StateID(ip)=sprintf('%s_tau%g',condNames{c},tau(j));
        Condition(ip)=string(condNames{c});
        ConditionLabel(ip)=string(condLabels{c});
        TauL_min(ip)=tau(j);
        RMST60_min(ip)=RMST60_fromFPT(j,c);
        Pcommit60(ip)=P60_fromFPT(j,c);
        PsiF60(ip)=Psi60(j,c);
        LambdaF_eq_h(ip)=lambdaRMST_h(j,c);
        LambdaOverR(ip)=lambdaOverR(j,c);
    end
end

stateTable = table(StateID,Condition,ConditionLabel,TauL_min, ...
    RMST60_min,Pcommit60,PsiF60,LambdaF_eq_h,LambdaOverR);

writetable(stateTable,outCSV);

%% ========================================================================
% 8. BUILD WINDOW-AUDIT TABLE
% =========================================================================
Window_min = AUDIT_T_MIN(:);
RankCorr_vs60 = rankCorrVs60;
PairwiseOrderAgreement_vs60 = orderAgreementVs60;

windowAuditTable = table(Window_min,RankCorr_vs60,PairwiseOrderAgreement_vs60);
writetable(windowAuditTable,outAuditCSV);

%% ========================================================================
% 9. SAVE A CLEAN MACRO INTERFACE STRUCT
% =========================================================================
INTERFACE = struct();

INTERFACE.version = 'MESO_TO_MACRO_RISK_INTERFACE_v1';
INTERFACE.created = datestr(now,30);

INTERFACE.source.fig4D = fig4DFile;
INTERFACE.source.fig4B = fig4BFile;

INTERFACE.condition.names = condNames;
INTERFACE.condition.labels = condLabels;
INTERFACE.tauL_min = tau;

INTERFACE.coarseGraining.primaryWindow_min = PRIMARY_T_MIN;
INTERFACE.coarseGraining.primaryWindow_h = TI_h;
INTERFACE.coarseGraining.auditWindows_min = AUDIT_T_MIN;
INTERFACE.coarseGraining.interpretation = [ ...
    'T_I is a fast-process coarse-graining window, not cell lifetime and ', ...
    'not macro experiment duration.'];

INTERFACE.meso.RMST60_min = RMST60_fromFPT;
INTERFACE.meso.Pcommit60 = P60_fromFPT;
INTERFACE.meso.PsiF60 = Psi60;
INTERFACE.meso.sourceOrientation.originalRMSTSize = size(D.RMST);
INTERFACE.meso.sourceOrientation.originalPcommitSize = size(D.Pcommit);
INTERFACE.meso.sourceOrientation.originalRSize = size(D.R);
INTERFACE.meso.sourceOrientation.internalConvention = 'tau x condition';

INTERFACE.macro.lambdaF_eq_h = lambdaRMST_h;
INTERFACE.macro.definition = [ ...
    'RMST-equivalent constant commitment hazard over the locked 60-min ', ...
    'coarse-graining window.'];
INTERFACE.macro.transition = 'V -> F';
INTERFACE.macro.stepDt_h = MACRO_DT_H;
INTERFACE.macro.Pcommit_per_step = Pstep;

% Stable heterogeneity library:
% one hazard for each mesoscopic state, not one hazard per stochastic FPT draw.
INTERFACE.heterogeneity.stateTable = stateTable;
INTERFACE.heterogeneity.note = [ ...
    'For Q2, construct cell-specific hazards by assigning mesoscopic states ', ...
    '(e.g. repair/tau_L states) and use the corresponding lambdaF_eq_h. ', ...
    'Do not convert individual stochastic FPT realizations into persistent ', ...
    'cell-specific hazards.'];

INTERFACE.audit.anchorPASS = anchorPASS;
INTERFACE.audit.maxRMSTSourceError_min = maxRMSTSourceErr;
INTERFACE.audit.maxPcommitSourceError = maxPSourceErr;
INTERFACE.audit.window.RMST_min = RMST_AUDIT_min;
INTERFACE.audit.window.Pcommit = P_AUDIT;
INTERFACE.audit.window.PsiF = Psi_AUDIT;
INTERFACE.audit.window.lambdaF_eq_h = Lambda_AUDIT_h;
INTERFACE.audit.window.rankCorrVs60 = rankCorrVs60;
INTERFACE.audit.window.orderAgreementVs60 = orderAgreementVs60;

INTERFACE.audit.r.available = isfinite(rMacro_h);
INTERFACE.audit.r.r_h = rMacro_h;
INTERFACE.audit.r.source = rSource;
INTERFACE.audit.r.epsilon_TIr = epsilon_r;
INTERFACE.audit.r.target = EPS_SCALE_TARGET;
INTERFACE.audit.r.scaleSeparationPASS = scaleSeparationPASS;
INTERFACE.audit.r.lambdaOverR = lambdaOverR;

INTERFACE.units.RMST = 'min';
INTERFACE.units.lambdaF = 'h^-1';
INTERFACE.units.macro_dt = 'h';

save(outMAT,'INTERFACE','stateTable','windowAuditTable','-v7.3');

%% ========================================================================
% 10. FIGURE
% =========================================================================
fig=figure('Color','w','Position',[80 80 1320 820]);
tl=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');

cc=lines(nCond);

% A: RMST60
ax=nexttile(tl,1); hold(ax,'on');
for c=1:nCond
    plot(ax,tau,RMST60_fromFPT(:,c),'-o','LineWidth',1.4, ...
        'Color',cc(c,:),'DisplayName',condNames{c});
end
xlabel(ax,'\tau_L (min)'); ylabel(ax,'RMST_{60} (min)');
title(ax,'a  Mesoscopic commitment timing');
legend(ax,'Location','best'); grid(ax,'on'); box(ax,'off');

% B: Psi
ax=nexttile(tl,2); hold(ax,'on');
for c=1:nCond
    plot(ax,tau,Psi60(:,c),'-o','LineWidth',1.4,'Color',cc(c,:));
end
xlabel(ax,'\tau_L (min)'); ylabel(ax,'\Psi_F = 1-RMST_{60}/60');
title(ax,'b  Mesoscopic ferroptotic pressure');
grid(ax,'on'); box(ax,'off');

% C: lambda eq
ax=nexttile(tl,3); hold(ax,'on');
for c=1:nCond
    plot(ax,tau,lambdaRMST_h(:,c),'-o','LineWidth',1.4,'Color',cc(c,:));
end
xlabel(ax,'\tau_L (min)'); ylabel(ax,'\lambda_F^{eq} (h^{-1})');
title(ax,'c  RMST-equivalent macro hazard');
grid(ax,'on'); box(ax,'off');

% D: window audit, all states
ax=nexttile(tl,4); hold(ax,'on');
for w=1:nWin
    vals=Lambda_AUDIT_h(:,:,w);
    plot(ax,1:nState,sort(vals(:)),'LineWidth',1.3, ...
        'DisplayName',sprintf('%g min',AUDIT_T_MIN(w)));
end
xlabel(ax,'State rank'); ylabel(ax,'\lambda_F^{eq} (h^{-1})');
title(ax,'d  Coarse-graining-window audit');
legend(ax,'Location','best'); grid(ax,'on'); box(ax,'off');

% E: rank/order robustness
ax=nexttile(tl,5); hold(ax,'on');
plot(ax,AUDIT_T_MIN,rankCorrVs60,'-o','LineWidth',1.4, ...
    'DisplayName','rank correlation');
plot(ax,AUDIT_T_MIN,orderAgreementVs60,'-s','LineWidth',1.4, ...
    'DisplayName','pairwise order agreement');
ylim(ax,[0 1.02]);
xlabel(ax,'Interface window (min)'); ylabel(ax,'Agreement vs 60 min');
title(ax,'e  Time-window robustness');
legend(ax,'Location','best'); grid(ax,'on'); box(ax,'off');

% F: exact step probabilities or lambda/r if r available
ax=nexttile(tl,6); hold(ax,'on');
if isfinite(rMacro_h)
    for c=1:nCond
        plot(ax,tau,lambdaOverR(:,c),'-o','LineWidth',1.4,'Color',cc(c,:));
    end
    yline(ax,1,'k:','\lambda_F=r');
    ylabel(ax,'\lambda_F^{eq}/r');
    title(ax,'f  Proliferation-commitment competition');
else
    for q=1:nDt
        vals=Pstep(:,:,q);
        plot(ax,1:nState,sort(vals(:)),'LineWidth',1.3, ...
            'DisplayName',sprintf('\\Deltat=%g h',MACRO_DT_H(q)));
    end
    ylabel(ax,'p_F = 1-e^{-\lambda_F\Deltat}');
    title(ax,'f  Exact hazard-to-step conversion');
    legend(ax,'Location','best');
end
xlabel(ax,'\tau_L / state rank');
grid(ax,'on'); box(ax,'off');

sgtitle(tl,'Mesoscopic-to-macroscopic ferroptotic-risk interface');
exportgraphics(fig,outPNG,'Resolution',300);

%% ========================================================================
% 11. FINAL CONSOLE SUMMARY
% =========================================================================
fprintf('\n==================== INTERFACE SUMMARY ====================\n');
fprintf('Fig4D raw-FPT reproduction : PASS\n');
fprintf('Fig4B tau=10 anchor        : PASS\n');
fprintf('Primary window             : %.0f min\n',PRIMARY_T_MIN);
fprintf('RMST60 range               : %.3f--%.3f min\n', ...
    min(RMST60_fromFPT(:)),max(RMST60_fromFPT(:)));
fprintf('PsiF60 range               : %.4f--%.4f\n', ...
    min(Psi60(:)),max(Psi60(:)));
fprintf('lambda_F^eq range          : %.4f--%.4f h^-1\n', ...
    min(lambdaRMST_h(:)),max(lambdaRMST_h(:)));
fprintf('30/45/60-min rank corr     : ');
fprintf('%.4f ',rankCorrVs60);
fprintf('\n');
fprintf('30/45/60-min order agree   : ');
fprintf('%.4f ',orderAgreementVs60);
fprintf('\n');

if isfinite(rMacro_h)
    fprintf('Locked r                   : %.6g h^-1\n',rMacro_h);
    fprintf('T_I*r                      : %.4g (%s)\n', ...
        epsilon_r,passfail(scaleSeparationPASS));
else
    fprintf('Locked r                   : not yet available; r-dependent audit deferred.\n');
end

fprintf('\nSaved:\n  %s\n  %s\n  %s\n  %s\n', ...
    outMAT,outCSV,outAuditCSV,outPNG);

fprintf('\nNEXT STEP:\n');
if isfinite(rMacro_h)
    fprintf(['Use INTERFACE.macro.lambdaF_eq_h directly in the updated macro ', ...
        'V->F transition and proceed to the light Pre-Q1 numerical audit.\n']);
else
    fprintf(['Calibrate/lock macro r independently (ferroptosis OFF), then rerun ', ...
        'this interface once to complete the T_I*r and lambda_F/r audit. ', ...
        'No mesoscopic simulation needs to be rerun.\n']);
end
fprintf('===========================================================\n');

end

%% ========================================================================
% LOCAL FUNCTIONS
% =========================================================================
function f = findFirstExisting(baseDir,candidates)
f='';
for k=1:numel(candidates)
    q=fullfile(baseDir,candidates{k});
    if isfile(q)
        f=q;
        return;
    end
end
end

function rmst = empiricalRMST(tf,T)
tf=double(tf(:));
assert(T>0);
% min(Inf,T)=T, which is exactly the restricted-survival contribution of a
% right-censored observation at the interface horizon.
rmst=mean(min(tf,T));
end

function lambda_h = rmstEquivalentHazard(rmst_h,T_h)
% Solve rmst/T = (1-exp(-x))/x, x=lambda*T.
assert(T_h>0);
ratio=rmst_h/T_h;

if ratio>=1-1e-14
    lambda_h=0;
    return;
end
if ratio<=0
    lambda_h=Inf;
    return;
end

fun=@(x) (-expm1(-x))/x-ratio;

xLo=1e-12;
xHi=max(10,2/ratio);
while fun(xHi)>0
    xHi=xHi*2;
    if xHi>1e8
        error('Could not bracket RMST-equivalent hazard.');
    end
end

x=fzero(fun,[xLo xHi]);
lambda_h=x/T_h;
end

function A = orientTauByCondition(A,nTau,nCond,label)
% Normalize a numeric state matrix to tau x condition.
sz=size(A);
if isequal(sz,[nTau nCond])
    return;
elseif isequal(sz,[nCond nTau])
    A=A.';
else
    error('%s has unexpected size %s; expected %dx%d or %dx%d.', ...
        label,mat2str(sz),nTau,nCond,nCond,nTau);
end
end

function C = orientTauByConditionCell(C,nTau,nCond,label)
% Normalize a cell state array to tau x condition.
assert(iscell(C),'%s must be a cell array.',label);
sz=size(C);
if isequal(sz,[nTau nCond])
    return;
elseif isequal(sz,[nCond nTau])
    C=C.';
else
    error('%s has unexpected size %s; expected %dx%d or %dx%d.', ...
        label,mat2str(sz),nTau,nCond,nCond,nTau);
end
end

function [rmst,p,tf] = getFig4BCondition(R,c,nCond)
% Robust extraction from either:
%   (i) a 1x4 / 4x1 struct array, or
%   (ii) a scalar struct whose fields contain 4 cells/elements.
if isstruct(R) && numel(R)==nCond
    rc=R(c);
    assert(isfield(rc,'RMST60') && isfield(rc,'Pcross60') && isfield(rc,'FPT_TS'), ...
        'Fig4B condition struct lacks required fields.');
    rmst=double(rc.RMST60);
    p=double(rc.Pcross60);
    tf=double(rc.FPT_TS(:));
    rmst=rmst(1);
    p=p(1);
elseif isstruct(R) && isscalar(R)
    assert(isfield(R,'RMST60') && isfield(R,'Pcross60') && isfield(R,'FPT_TS'), ...
        'Fig4B scalar struct lacks required fields.');
    rmst=scalarFromMaybeCell(R.RMST60,c);
    p=scalarFromMaybeCell(R.Pcross60,c);
    tf=vectorFromMaybeCell(R.FPT_TS,c);
else
    error('Unexpected Fig4B R structure: class=%s, size=%s.', ...
        class(R),mat2str(size(R)));
end
end

function x = scalarFromMaybeCell(v,idx)
if iscell(v)
    x=double(v{idx});
else
    x=double(v(idx));
end
x=x(1);
end

function x = vectorFromMaybeCell(v,idx)
if iscell(v)
    x=double(v{idx}(:));
else
    if isvector(v)
        x=double(v(:));
    else
        x=double(v(:,idx));
    end
end
end

function r = tiedRankLocal(x)
x=double(x(:));
[sv,ord]=sort(x);
r=zeros(size(x));
n=numel(x);
i=1;
while i<=n
    j=i;
    while j<n && sv(j+1)==sv(i)
        j=j+1;
    end
    rr=(i+j)/2;
    r(ord(i:j))=rr;
    i=j+1;
end
end

function rho = pearsonLocal(x,y)
x=double(x(:)); y=double(y(:));
x=x-mean(x); y=y-mean(y);
den=sqrt(sum(x.^2)*sum(y.^2));
if den==0
    rho=NaN;
else
    rho=sum(x.*y)/den;
end
end

function a = pairwiseOrderAgreement(x,ref)
% Fraction of unordered pairs whose ordering relative to the primary
% interface is preserved. Ties in either vector are counted as agreement
% only when both vectors tie.
x=double(x(:));
ref=double(ref(:));
n=numel(x);
agree=0;
total=0;

for i=1:n-1
    for j=i+1:n
        sx=sign(x(i)-x(j));
        sr=sign(ref(i)-ref(j));
        agree=agree+(sx==sr);
        total=total+1;
    end
end
a=agree/total;
end

function s=passfail(tf)
if tf
    s='PASS';
else
    s='FAIL';
end
end
