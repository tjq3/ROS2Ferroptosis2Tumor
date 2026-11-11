function SI_M2_Macro_Parameter_Benchmarks_FINAL()
% =========================================================================
% SI_M2_Macro_Parameter_Benchmarks_FINAL
%
% FINAL SI figure for literature-informed macro parameter benchmarks.
%
% A. Free-proliferation reference doubling time -> intrinsic r.
% B. Ferroptosis-OFF spatial growth generated with the FINAL event-driven
%    direct SSA over the locked r sensitivity range.
% C. Mean F-state residence time -> dimensionless clearance coordinate mu/r.
%
% IMPORTANT
% ---------
% - r and mu are NOT re-fitted here. Frozen values are read and verified.
% - No fixed numerical dt is used in the spatial validation.
% - record_dt_h controls output sampling only.
% =========================================================================

clc; close all;

%% 0. SETTINGS
TD_SENS_H = [16 24 36 48];
TD_REFERENCE_H = 24;
RESIDENCE_SENS_H = [3 6 12 24 48];
RESIDENCE_REFERENCE_H = 12;

L = 61;
R0 = 4;
TEND_H = 168;
RECORD_DT_H = 1;
FIT_WINDOW_H = [12 72];

NREP_SCAN = 40;
NREP_REFERENCE = 80;
SEED_SCAN = 610000;
SEED_REFERENCE = 620000;

USE_PARALLEL = true;
AUTO_START_POOL = true;

%% 1. FILES / FROZEN PARAMETERS
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir = pwd; end

coreFile = fullfile(scriptDir,'macro_sim_core_hazard_SSA_FINAL.m');
rFile = fullfile(scriptDir,'Macro_r_Calibration_OUT','Macro_r_LOCKED.mat');
muFile = fullfile(scriptDir,'Macro_mu_Calibration_OUT','Macro_mu_LOCKED.mat');
ssaLockFile = fullfile(scriptDir,'Macro_Hazard_SSA_Domain_Confirm_OUT','Macro_SSA_LOCKED.mat');

assert(isfile(coreFile),'Missing macro_sim_core_hazard_SSA_FINAL.m.');
assert(isfile(rFile),'Missing Macro_r_LOCKED.mat.');
assert(isfile(muFile),'Missing Macro_mu_LOCKED.mat.');

RLOCK = load(rFile,'LOCK');
MLOCK = load(muFile,'LOCK');

rStar = double(RLOCK.LOCK.r_star_per_h);
muStar = double(MLOCK.LOCK.mu_star_per_h);

rAnalytic = log(2)./TD_SENS_H;
rReferenceExpected = log(2)/TD_REFERENCE_H;
muGrid = 1./RESIDENCE_SENS_H;
muReferenceExpected = 1/RESIDENCE_REFERENCE_H;

assert(abs(rStar-rReferenceExpected)<1e-12, ...
    'Frozen r does not match the locked 24-h reference.');
assert(abs(muStar-muReferenceExpected)<1e-12, ...
    'Frozen mu does not match the locked 12-h residence reference.');

if isfile(ssaLockFile)
    S=load(ssaLockFile);
    if isfield(S,'SSA_LOCK')
        assert(S.SSA_LOCK.readyForQ1Pilot,'SSA lock is not ready for Q1.');
        assert(S.SSA_LOCK.L==L,'Frozen SSA domain differs from L=61.');
    elseif isfield(S,'LOCK')
        assert(S.LOCK.readyForQ1Pilot,'SSA lock is not ready for Q1.');
    end
end

%% 2. PARALLEL POOL
poolReady=false;
if USE_PARALLEL && license('test','Distrib_Computing_Toolbox')
    p=gcp('nocreate');
    if isempty(p) && AUTO_START_POOL
        try
            p=parpool('threads');
        catch ME
            warning('Could not start thread pool (%s). Using serial.',ME.message);
            p=[];
        end
    end
    poolReady=~isempty(p);
end

%% 3. EVENT-DRIVEN FERROPTOSIS-OFF r SENSITIVITY
nR=numel(rAnalytic);
SCAN=cell(nR,1);

fprintf('\n============================================================\n');
fprintf('SI-M2 MACRO PARAMETER BENCHMARKS\n');
fprintf('Frozen r  = %.8f h^-1 (T_D,ref=24 h)\n',rStar);
fprintf('Frozen mu = %.8f h^-1 (mean F residence=12 h)\n',muStar);
fprintf('Spatial solver: event-driven direct SSA; lambda_F=0\n');
fprintf('============================================================\n');

for k=1:nR
    cfg=makeCfg(L,R0,TEND_H,RECORD_DT_H,rAnalytic(k));
    fprintf('Running r benchmark T_D=%g h, r=%.8f h^-1 ...\n',TD_SENS_H(k),rAnalytic(k));
    SCAN{k}=runSSAEnsemble(cfg,NREP_SCAN,SEED_SCAN+10000*k,poolReady,true);
end

%% 4. INDEPENDENT REFERENCE-r ENSEMBLE
cfgRef=makeCfg(L,R0,TEND_H,RECORD_DT_H,rStar);
REF=runSSAEnsemble(cfgRef,NREP_REFERENCE,SEED_REFERENCE,poolReady,true);

[tRef,meanRef,loRef,hiRef]=ensembleCurve(REF.t_h,REF.NVcurve);
gRef=fitLogSlope(tRef,meanRef,FIT_WINDOW_H);
finalRef=REF.finalNV;
relMCSERef=std(finalRef,0)/sqrt(numel(finalRef))/mean(finalRef);
edgeRef=mean(REF.edgeTouch);

%% 5. SUMMARY TABLES
Td_ref_h = TD_SENS_H(:);
r_per_h = rAnalytic(:);
earlyGrowth_h = nan(nR,1);
finalMeanNV = nan(nR,1);
edgeTouchProb = nan(nR,1);

for k=1:nR
    [t,mu,~,~]=ensembleCurve(SCAN{k}.t_h,SCAN{k}.NVcurve);
    earlyGrowth_h(k)=fitLogSlope(t,mu,FIT_WINDOW_H);
    finalMeanNV(k)=mean(SCAN{k}.finalNV);
    edgeTouchProb(k)=mean(SCAN{k}.edgeTouch);
end

Tr = table(Td_ref_h,r_per_h,earlyGrowth_h,finalMeanNV,edgeTouchProb);

Residence_h = RESIDENCE_SENS_H(:);
mu_per_h = muGrid(:);
mu_over_r = muGrid(:)/rStar;
Tmu = table(Residence_h,mu_per_h,mu_over_r);

fprintf('\nReference-r event-driven SSA validation:\n');
fprintf('  mean final NV = %.3f\n',mean(finalRef));
fprintf('  early g       = %.7f h^-1\n',gRef);
fprintf('  edge touch    = %.4f\n',edgeRef);
fprintf('  final-NV MCSE = %.3f%%\n',100*relMCSERef);
fprintf('mu/r reference  = %.3f\n',muStar/rStar);

%% 6. FIGURE
fig=figure('Color','w','Position',[100 100 1460 440]);
tl=tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');

% A. Analytic r scale.
ax=nexttile(tl,1);
plot(ax,TD_SENS_H,rAnalytic,'-o','LineWidth',1.5); hold(ax,'on');
plot(ax,TD_REFERENCE_H,rStar,'ks','MarkerFaceColor','w','MarkerSize',8,'LineWidth',1.2);
xlabel(ax,'Reference free-proliferation doubling time (h)');
ylabel(ax,'r = ln(2)/T_D (h^{-1})');
title(ax,'a  Intrinsic proliferation benchmark');
grid(ax,'on'); box(ax,'off');

% B. FINAL event-driven SSA no-ferroptosis baseline.
ax=nexttile(tl,2); hold(ax,'on');
for k=1:nR
    [t,mu,~,~]=ensembleCurve(SCAN{k}.t_h,SCAN{k}.NVcurve);
    plot(ax,t,mu,'LineWidth',1.45,'DisplayName',sprintf('T_D=%g h',TD_SENS_H(k)));
end
% Show the independent reference validation band lightly.
fill(ax,[tRef;flipud(tRef)],[loRef;flipud(hiRef)],[0.8 0.8 0.8], ...
    'FaceAlpha',0.18,'EdgeColor','none','HandleVisibility','off');
plot(ax,tRef,meanRef,'k--','LineWidth',1.1,'DisplayName','24-h independent validation');
xlabel(ax,'Time (h)'); ylabel(ax,'Mean viable cells N_V');
title(ax,'b  Ferroptosis-OFF event-driven SSA baseline');
legend(ax,'Location','northwest','Box','off'); grid(ax,'on'); box(ax,'off');

% C. Clearance scale.
ax=nexttile(tl,3);
plot(ax,RESIDENCE_SENS_H,muGrid/rStar,'-o','LineWidth',1.5); hold(ax,'on');
plot(ax,RESIDENCE_REFERENCE_H,muStar/rStar,'ks','MarkerFaceColor','w','MarkerSize',8,'LineWidth',1.2);
yline(ax,1,'k:','\mu=r');
xlabel(ax,'Mean F-state residence 1/\mu (h)');
ylabel(ax,'\chi_\mu = \mu/r');
title(ax,'c  Clearance benchmark and sensitivity');
grid(ax,'on'); box(ax,'off');

sgtitle(tl,'Literature-informed macro parameter benchmarks');

%% 7. SAVE
outDir=fullfile(scriptDir,'SI_M2_Macro_Parameter_Benchmarks_FINAL_OUT');
if ~exist(outDir,'dir'), mkdir(outDir); end
outPNG=fullfile(outDir,'SI_M2_Macro_Parameter_Benchmarks_FINAL.png');
outPDF=fullfile(outDir,'SI_M2_Macro_Parameter_Benchmarks_FINAL.pdf');
outFIG=fullfile(outDir,'SI_M2_Macro_Parameter_Benchmarks_FINAL.fig');
outMAT=fullfile(outDir,'SI_M2_Macro_Parameter_Benchmarks_FINAL.mat');
outCSVr=fullfile(outDir,'SI_M2_r_Benchmark.csv');
outCSVmu=fullfile(outDir,'SI_M2_mu_Benchmark.csv');

try
    axs=findall(fig,'Type','axes');
    for k=1:numel(axs), axtoolbar(axs(k),{}); end
catch
end

exportgraphics(fig,outPNG,'Resolution',600);
exportgraphics(fig,outPDF,'ContentType','vector');
savefig(fig,outFIG);
writetable(Tr,outCSVr); writetable(Tmu,outCSVmu);

SI_M2=struct();
SI_M2.rStar_h=rStar;
SI_M2.muStar_h=muStar;
SI_M2.TdSensitivity_h=TD_SENS_H;
SI_M2.rSensitivity_h=rAnalytic;
SI_M2.residenceSensitivity_h=RESIDENCE_SENS_H;
SI_M2.muSensitivity_h=muGrid;
SI_M2.muOverR=muGrid/rStar;
SI_M2.referenceSSA.meanFinalNV=mean(finalRef);
SI_M2.referenceSSA.earlyGrowth_h=gRef;
SI_M2.referenceSSA.edgeTouchProb=edgeRef;
SI_M2.referenceSSA.relMCSEFinal=relMCSERef;
SI_M2.note='Frozen r and mu are verified, not refitted.';
save(outMAT,'SI_M2','SCAN','REF','Tr','Tmu','-v7.3');

fprintf('\nSaved SI-M2 outputs:\n  %s\n  %s\n  %s\n',outPNG,outPDF,outMAT);
end

%% LOCAL FUNCTIONS
function cfg=makeCfg(L,R0,Tend,recordDt,r)
cfg=struct();
cfg.L=L; cfg.R0=R0; cfg.Tend_h=Tend; cfg.record_dt_h=recordDt;
cfg.r_h=r; cfg.mu_h=0; cfg.seed=1;
cfg.lambda_mode='scalar'; cfg.lambdaF_h=0;
cfg.snapshot_times_h=[]; cfg.return_final_state=false;
cfg.compute_cavity=false; cfg.store_event_log=false; cfg.verbose=false;
end

function E=runSSAEnsemble(cfg,nRep,seedBase,poolReady,storeCurves)
tmp=cell(nRep,1);
if poolReady
    parfor k=1:nRep
        c=cfg; c.seed=seedBase+k;
        O=macro_sim_core_hazard_SSA_FINAL(c);
        tmp{k}=reduceSSA(O,storeCurves);
    end
else
    for k=1:nRep
        c=cfg; c.seed=seedBase+k;
        O=macro_sim_core_hazard_SSA_FINAL(c);
        tmp{k}=reduceSSA(O,storeCurves);
    end
end
E=stackReduced(tmp,storeCurves);
end

function z=reduceSSA(O,storeCurves)
z.finalNV=double(O.NV(end));
z.edgeTouch=any(O.edgeOcc);
z.bookkeepingPASS=logical(O.audit.bookkeepingPASS);
if storeCurves
    z.t_h=double(O.t_h(:)); z.NV=double(O.NV(:));
else
    z.t_h=[]; z.NV=[];
end
end

function E=stackReduced(tmp,storeCurves)
n=numel(tmp); E.finalNV=zeros(n,1); E.edgeTouch=false(n,1); E.bookkeepingPASS=false(n,1);
if storeCurves
    E.t_h=tmp{1}.t_h; E.NVcurve=zeros(numel(E.t_h),n);
else
    E.t_h=[]; E.NVcurve=[];
end
for k=1:n
    E.finalNV(k)=tmp{k}.finalNV; E.edgeTouch(k)=tmp{k}.edgeTouch; E.bookkeepingPASS(k)=tmp{k}.bookkeepingPASS;
    if storeCurves, E.NVcurve(:,k)=tmp{k}.NV; end
end
assert(all(E.bookkeepingPASS),'SSA bookkeeping failure detected.');
end

function [t,mu,lo,hi]=ensembleCurve(t,M)
mu=mean(M,2); lo=zeros(size(mu)); hi=zeros(size(mu));
for k=1:numel(mu)
    lo(k)=pctLocal(M(k,:),2.5); hi(k)=pctLocal(M(k,:),97.5);
end
end

function g=fitLogSlope(t,mu,win)
mask=t>=win(1) & t<=win(2) & mu>0;
if nnz(mask)<3, g=NaN; return; end
p=polyfit(t(mask),log(mu(mask)),1); g=p(1);
end

function q=pctLocal(v,p)
v=sort(v(isfinite(v)));
if isempty(v), q=NaN; return; end
pos=1+(numel(v)-1)*p/100; a=floor(pos); b=ceil(pos);
if a==b, q=v(a); else, q=v(a)+(pos-a)*(v(b)-v(a)); end
end
