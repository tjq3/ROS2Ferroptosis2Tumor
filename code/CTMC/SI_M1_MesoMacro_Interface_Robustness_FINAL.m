function SI_M1_MesoMacro_Interface_Robustness_FINAL()
% =========================================================================
% SI_M1_MesoMacro_Interface_Robustness_FINAL
%
% FINAL SI figure for robustness of mesoscopic -> macroscopic
% ferroptotic-risk coarse-graining.
%
% A. 30/45/60-min RMST-equivalent hazard profiles across 48 states.
% B. Rank correlation + pairwise-order agreement relative to 60 min.
% C. Meso-macro time-scale separation epsilon_r = T_I*r over the locked
%    literature-informed r sensitivity range.
%
% IMPORTANT:
% - NO macro numerical dt is used here.
% - NO hazard-to-step conversion is computed.
% - T_I=60 min is a coarse-graining window, not a cell lifetime.
% =========================================================================

clc; close all;

%% 0. SETTINGS
PRIMARY_T_MIN = 60;
AUDIT_T_MIN = [30 45 60];
EPS_SCALE_TARGET = 0.10;
ANCHOR_TOL_RMST = 1e-10;
ANCHOR_TOL_P = 1e-12;
SOURCE_TOL_RMST = 1e-9;
SOURCE_TOL_P = 1e-12;

%% 1. FILES
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

rFile = fullfile(scriptDir,'Macro_r_Calibration_OUT','Macro_r_LOCKED.mat');

assert(~isempty(fig4DFile),'Cannot find Fig4D mesoscopic data file.');
assert(~isempty(fig4BFile),'Cannot find Fig4B anchor data file.');
assert(isfile(rFile),'Cannot find Macro_r_LOCKED.mat.');

D = load(fig4DFile);
B = load(fig4BFile);
RLOCK = load(rFile,'LOCK');

rStar = double(RLOCK.LOCK.r_star_per_h);

if isfield(RLOCK.LOCK,'sensitivityDoublingReference_h')
    TdSens = double(RLOCK.LOCK.sensitivityDoublingReference_h(:).');
else
    TdSens = [16 24 36 48];
end

if isfield(RLOCK.LOCK,'sensitivity_r_per_h')
    rSens = double(RLOCK.LOCK.sensitivity_r_per_h(:).');
else
    rSens = log(2)./TdSens;
end

%% 2. VERIFY MESOSCOPIC SOURCE
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
nTau = numel(tau);
nCond = numel(condNames);

assert(nTau==12,'Expected 12 tau_L states.');
assert(nCond==4,'Expected 4 conditions.');
assert(abs(double(D.meta.Tobs_min)-PRIMARY_T_MIN)<1e-12, ...
    'Fig4D Tobs differs from the locked 60-min interface window.');

expectedCond = {'CU','GU','CF','GF'};
assert(isequal(condNames(:)',expectedCond),'Unexpected condition ordering.');

RMST_SRC = orientTauByCondition(double(D.RMST),nTau,nCond,'Fig4D RMST');
PCOMMIT_SRC = orientTauByCondition(double(D.Pcommit),nTau,nCond,'Fig4D Pcommit');
R_SRC = orientTauByConditionCell(D.R,nTau,nCond,'Fig4D R');

%% 3. REPRODUCE 60-min SOURCE SUMMARIES
RMST60 = nan(nTau,nCond);
P60 = nan(nTau,nCond);
Lambda60_h = nan(nTau,nCond);

for j=1:nTau
    for c=1:nCond
        tf = double(R_SRC{j,c}.FPT_TS(:));
        assert(~any(isnan(tf)),'NaN FPT in tau=%g, %s.',tau(j),condNames{c});
        assert(all(tf(isfinite(tf))>=0),'Negative FPT in tau=%g, %s.',tau(j),condNames{c});
        RMST60(j,c) = empiricalRMST(tf,PRIMARY_T_MIN);
        P60(j,c) = mean(tf<=PRIMARY_T_MIN);
        Lambda60_h(j,c) = rmstEquivalentHazard(RMST60(j,c)/60,PRIMARY_T_MIN/60);
    end
end

maxRMSTErr = max(abs(RMST60(:)-RMST_SRC(:)));
maxPErr = max(abs(P60(:)-PCOMMIT_SRC(:)));
assert(maxRMSTErr<SOURCE_TOL_RMST,'Raw FPT does not reproduce Fig4D RMST.');
assert(maxPErr<SOURCE_TOL_P,'Raw FPT does not reproduce Fig4D Pcommit.');

%% 4. FIG4B ANCHOR AT tau_L=10 min
[~,iRef] = min(abs(tau-double(D.meta.tauReference_min)));
assert(abs(tau(iRef)-10)<1e-12,'Expected tau_L reference = 10 min.');

BcondNames = cellstr(string(B.condNames(:)));
assert(isequal(BcondNames(:)',expectedCond),'Fig4B and Fig4D condition ordering differs.');

B_RMST = nan(1,nCond); B_P = nan(1,nCond); rawMatch = true(1,nCond);
for c=1:nCond
    [B_RMST(c),B_P(c),tfB] = getFig4BCondition(B.R,c,nCond);
    tfD = double(R_SRC{iRef,c}.FPT_TS(:));
    if ~isequaln(tfB(:),tfD(:))
        sameSize = numel(tfB)==numel(tfD);
        sameInf = sameSize && isequal(isinf(tfB(:)),isinf(tfD(:)));
        if sameInf
            m = isfinite(tfB(:)) & isfinite(tfD(:));
            rawMatch(c) = all(abs(tfB(m)-tfD(m))<1e-12);
        else
            rawMatch(c) = false;
        end
    end
end

anchorRMSTErr = max(abs(B_RMST-RMST60(iRef,:)));
anchorPErr = max(abs(B_P-P60(iRef,:)));
anchorPASS = anchorRMSTErr<ANCHOR_TOL_RMST && anchorPErr<ANCHOR_TOL_P && all(rawMatch);
assert(anchorPASS,'Fig4B -> Fig4D anchor failed.');

%% 5. 30/45/60-min COARSE-GRAINING AUDIT
nWin = numel(AUDIT_T_MIN);
LambdaAudit_h = nan(nTau,nCond,nWin);

for w=1:nWin
    Tw = AUDIT_T_MIN(w);
    for j=1:nTau
        for c=1:nCond
            tf = double(R_SRC{j,c}.FPT_TS(:));
            rmst = empiricalRMST(tf,Tw);
            LambdaAudit_h(j,c,w) = rmstEquivalentHazard(rmst/60,Tw/60);
        end
    end
end

iPrimary = find(abs(AUDIT_T_MIN-PRIMARY_T_MIN)<1e-12,1);
assert(~isempty(iPrimary));
assert(max(abs(LambdaAudit_h(:,:,iPrimary)-Lambda60_h),[],'all')<1e-10);

primaryFlat = Lambda60_h(:);
primaryRank = tiedRankLocal(primaryFlat);
rankCorrVs60 = nan(nWin,1);
orderAgreementVs60 = nan(nWin,1);

for w=1:nWin
    x = LambdaAudit_h(:,:,w); x = x(:);
    rankCorrVs60(w) = pearsonLocal(tiedRankLocal(x),primaryRank);
    orderAgreementVs60(w) = pairwiseOrderAgreement(x,primaryFlat);
end

%% 6. MESO-MACRO TIME-SCALE SEPARATION
TI_h = PRIMARY_T_MIN/60;
epsilonSens = TI_h*rSens;
epsilonStar = TI_h*rStar;
scalePASS = epsilonStar<=EPS_SCALE_TARGET && all(epsilonSens<=EPS_SCALE_TARGET);

%% 7. TABLES
Window_min = AUDIT_T_MIN(:);
RankCorr_vs60 = rankCorrVs60;
PairwiseOrderAgreement_vs60 = orderAgreementVs60;
Twindow = table(Window_min,RankCorr_vs60,PairwiseOrderAgreement_vs60);

Td_ref_h = TdSens(:);
r_per_h = rSens(:);
epsilon_TIr = epsilonSens(:);
Tscale = table(Td_ref_h,r_per_h,epsilon_TIr);

%% 8. CONSOLE
fprintf('\n============================================================\n');
fprintf('SI-M1 MESO->MACRO INTERFACE ROBUSTNESS\n');
fprintf('Fig4D raw-FPT reproduction: PASS\n');
fprintf('Fig4B tau=10 anchor       : PASS\n');
fprintf('Primary T_I              : %.0f min\n',PRIMARY_T_MIN);
for w=1:nWin
    fprintf('T=%2.0f min | rank corr=%.4f | order agreement=%.4f\n', ...
        AUDIT_T_MIN(w),rankCorrVs60(w),orderAgreementVs60(w));
end
fprintf('Locked r = %.8f h^-1\n',rStar);
fprintf('T_I*r = %.5f\n',epsilonStar);
fprintf('All r-sensitivity T_I*r <= %.2f: %s\n',EPS_SCALE_TARGET,passfail(scalePASS));
fprintf('============================================================\n');

%% 9. FIGURE
fig = figure('Color','w','Position',[100 100 1450 430]);
tl = tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');

ax = nexttile(tl,1); hold(ax,'on');
for w=1:nWin
    vals = LambdaAudit_h(:,:,w);
    plot(ax,1:numel(vals),sort(vals(:)),'LineWidth',1.5, ...
        'DisplayName',sprintf('%g min',AUDIT_T_MIN(w)));
end
xlabel(ax,'Mesoscopic state rank'); ylabel(ax,'\lambda_F^{eq} (h^{-1})');
title(ax,'a  Coarse-graining-window robustness');
legend(ax,'Location','northwest','Box','off'); grid(ax,'on'); box(ax,'off');

ax = nexttile(tl,2); hold(ax,'on');
plot(ax,AUDIT_T_MIN,rankCorrVs60,'-o','LineWidth',1.5,'DisplayName','rank correlation');
plot(ax,AUDIT_T_MIN,orderAgreementVs60,'-s','LineWidth',1.5,'DisplayName','pairwise order agreement');
ylim(ax,[0.95 1.005]);
xlabel(ax,'Interface window T_I (min)'); ylabel(ax,'Agreement with 60-min interface');
title(ax,'b  Preservation of state ordering');
legend(ax,'Location','southeast','Box','off'); grid(ax,'on'); box(ax,'off');

ax = nexttile(tl,3);
plot(ax,TdSens,epsilonSens,'-o','LineWidth',1.5); hold(ax,'on');
plot(ax,24,epsilonStar,'ks','MarkerFaceColor','w','MarkerSize',8,'LineWidth',1.2);
yline(ax,EPS_SCALE_TARGET,'k--','diagnostic target');
xlabel(ax,'Reference free-proliferation doubling time (h)'); ylabel(ax,'\epsilon_r = T_I r');
title(ax,'c  Meso-macro time-scale separation'); grid(ax,'on'); box(ax,'off');

sgtitle(tl,'Robustness of mesoscopic-to-macroscopic coarse-graining');

%% 10. SAVE
outDir = fullfile(scriptDir,'SI_M1_MesoMacro_Interface_Robustness_FINAL_OUT');
if ~exist(outDir,'dir'), mkdir(outDir); end
outPNG = fullfile(outDir,'SI_M1_MesoMacro_Interface_Robustness_FINAL.png');
outPDF = fullfile(outDir,'SI_M1_MesoMacro_Interface_Robustness_FINAL.pdf');
outFIG = fullfile(outDir,'SI_M1_MesoMacro_Interface_Robustness_FINAL.fig');
outMAT = fullfile(outDir,'SI_M1_MesoMacro_Interface_Robustness_FINAL.mat');
outCSV1 = fullfile(outDir,'SI_M1_MesoMacro_WindowAudit.csv');
outCSV2 = fullfile(outDir,'SI_M1_MesoMacro_TimeScaleAudit.csv');

try
    axs=findall(fig,'Type','axes');
    for k=1:numel(axs), axtoolbar(axs(k),{}); end
catch
end

exportgraphics(fig,outPNG,'Resolution',600);
exportgraphics(fig,outPDF,'ContentType','vector');
savefig(fig,outFIG);
writetable(Twindow,outCSV1); writetable(Tscale,outCSV2);

SI_M1 = struct();
SI_M1.primaryWindow_min = PRIMARY_T_MIN;
SI_M1.auditWindows_min = AUDIT_T_MIN;
SI_M1.lambdaAudit_h = LambdaAudit_h;
SI_M1.rankCorrVs60 = rankCorrVs60;
SI_M1.orderAgreementVs60 = orderAgreementVs60;
SI_M1.rStar_h = rStar;
SI_M1.rSensitivity_h = rSens;
SI_M1.TdSensitivity_h = TdSens;
SI_M1.epsilonStar = epsilonStar;
SI_M1.epsilonSensitivity = epsilonSens;
SI_M1.scalePASS = scalePASS;
SI_M1.anchorPASS = anchorPASS;
save(outMAT,'SI_M1','Twindow','Tscale','-v7.3');

fprintf('\nSaved SI-M1 outputs:\n  %s\n  %s\n  %s\n',outPNG,outPDF,outMAT);
end

%% LOCAL FUNCTIONS
function f = findFirstExisting(baseDir,candidates)
f='';
for k=1:numel(candidates)
    q=fullfile(baseDir,candidates{k});
    if isfile(q), f=q; return; end
end
end

function rmst = empiricalRMST(tf,T)
tf=double(tf(:)); rmst=mean(min(tf,T));
end

function lambda_h = rmstEquivalentHazard(rmst_h,T_h)
ratio=rmst_h/T_h;
if ratio>=1-1e-14, lambda_h=0; return; end
if ratio<=0, lambda_h=Inf; return; end
fun=@(x) (-expm1(-x))/x-ratio;
xLo=1e-12; xHi=max(10,2/ratio);
while fun(xHi)>0
    xHi=xHi*2;
    if xHi>1e8, error('Could not bracket RMST-equivalent hazard.'); end
end
x=fzero(fun,[xLo xHi]); lambda_h=x/T_h;
end

function A = orientTauByCondition(A,nTau,nCond,label)
sz=size(A);
if isequal(sz,[nTau nCond]), return;
elseif isequal(sz,[nCond nTau]), A=A.';
else, error('%s has unexpected size %s.',label,mat2str(sz)); end
end

function C = orientTauByConditionCell(C,nTau,nCond,label)
assert(iscell(C),'%s must be a cell array.',label);
sz=size(C);
if isequal(sz,[nTau nCond]), return;
elseif isequal(sz,[nCond nTau]), C=C.';
else, error('%s has unexpected size %s.',label,mat2str(sz)); end
end

function [rmst,p,tf] = getFig4BCondition(R,c,nCond)
if isstruct(R) && numel(R)==nCond
    rc=R(c); rmst=double(rc.RMST60); p=double(rc.Pcross60); tf=double(rc.FPT_TS(:));
    rmst=rmst(1); p=p(1);
elseif isstruct(R) && isscalar(R)
    rmst=scalarFromMaybeCell(R.RMST60,c); p=scalarFromMaybeCell(R.Pcross60,c); tf=vectorFromMaybeCell(R.FPT_TS,c);
else
    error('Unexpected Fig4B R structure.');
end
end

function x = scalarFromMaybeCell(v,idx)
if iscell(v), x=double(v{idx}); else, x=double(v(idx)); end
x=x(1);
end

function x = vectorFromMaybeCell(v,idx)
if iscell(v), x=double(v{idx}(:));
else
    if isvector(v), x=double(v(:)); else, x=double(v(:,idx)); end
end
end

function r = tiedRankLocal(x)
x=double(x(:)); [sv,ord]=sort(x); r=zeros(size(x)); n=numel(x); i=1;
while i<=n
    j=i; while j<n && sv(j+1)==sv(i), j=j+1; end
    rr=(i+j)/2; r(ord(i:j))=rr; i=j+1;
end
end

function rho = pearsonLocal(x,y)
x=double(x(:)); y=double(y(:)); x=x-mean(x); y=y-mean(y);
den=sqrt(sum(x.^2)*sum(y.^2));
if den==0, rho=NaN; else, rho=sum(x.*y)/den; end
end

function a = pairwiseOrderAgreement(x,ref)
x=double(x(:)); ref=double(ref(:)); n=numel(x); agree=0; total=0;
for i=1:n-1
    for j=i+1:n
        sx=sign(x(i)-x(j)); sr=sign(ref(i)-ref(j));
        agree=agree+(sx==sr); total=total+1;
    end
end
a=agree/total;
end

function s=passfail(tf)
if tf, s='PASS'; else, s='FAIL'; end
end
