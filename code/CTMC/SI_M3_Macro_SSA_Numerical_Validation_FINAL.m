function SI_M3_Macro_SSA_Numerical_Validation_FINAL()
% =========================================================================
% SI_M3_Macro_SSA_Numerical_Validation_FINAL
%
% FINAL SI numerical validation of the event-driven hazard-based macro SSA.
%
% A. No-ferroptosis ensemble regression against the locked r benchmark.
% B. Supportive solver-transition consistency versus the previous finest
%    fixed-step implementation (NOT treated as ground truth).
% C. Monte-Carlo precision across representative chi=lambda_F/r states.
% D. Final common-seed domain confirmation: L=61 vs L=71 at chi=[0,2].
%
% PRIMARY VALIDATION GATES
% ------------------------
% Regression:
%   <=10% relative difference in final N_V and early growth rate,
%   no boundary touch, exact bookkeeping.
%
% Monte-Carlo:
%   relative MCSE <=2% for exact AUTC_NV and AUTC_NT;
%   Wilson 95% half-width <=0.12 for extinction probability.
%
% Domain:
%   paired AUTC effects <=2%, bootstrap 95% CI within +/-5%,
%   zero boundary touch in L=61 and L=71.
%
% The old fixed-step comparison is supportive only and cannot veto an
% otherwise clean event-driven SSA validation.
% =========================================================================

clc; close all;

%% 0. SETTINGS
L_MAIN = 61;
L_REF = 71;
R0 = 4;
TEND_H = 168;
RECORD_DT_H = 1;  % output sampling only; no numerical solver dt
FIT_WINDOW_H = [12 72];

CHI_XVAL = [0 1 2 8 135];
CHI_MC = [0 1 1.5 2 3 8 135];
CHI_DOMAIN = [0 2];

N_REG = 80;
N_XVAL = 60;
N_MC_BASE = 80;
N_MC_MAX = 160;
MC_BATCH = 20;
N_DOMAIN = 50;
N_BOOT = 5000;

REG_REL_TOL = 0.10;
MCSE_TOL = 0.02;
EXT_WILSON_HALF_TOL = 0.12;
DOMAIN_POINT_TOL = 0.02;
DOMAIN_CI_TOL = 0.05;

SEED_REG = 710000;
SEED_XVAL = 720000;
SEED_MC = 730000;
SEED_DOMAIN = 740000;
SEED_BOOT = 750000;

USE_PARALLEL = true;
AUTO_START_POOL = true;

%% 1. FILES / LOCKED PARAMETERS
scriptDir=fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir=pwd; end

coreFile=fullfile(scriptDir,'macro_sim_core_hazard_SSA_FINAL.m');
rFile=fullfile(scriptDir,'Macro_r_Calibration_OUT','Macro_r_LOCKED.mat');
muFile=fullfile(scriptDir,'Macro_mu_Calibration_OUT','Macro_mu_LOCKED.mat');
oldDtFile=fullfile(scriptDir,'Macro_Hazard_PreQ1_DT_Refinement_OUT', ...
    'Macro_Hazard_PreQ1_DT_Refinement_FINAL.mat');

assert(isfile(coreFile),'Missing macro_sim_core_hazard_SSA_FINAL.m.');
assert(isfile(rFile),'Missing Macro_r_LOCKED.mat.');
assert(isfile(muFile),'Missing Macro_mu_LOCKED.mat.');

RLOCK=load(rFile);
MLOCK=load(muFile,'LOCK');
r_h=double(RLOCK.LOCK.r_star_per_h);
mu_h=double(MLOCK.LOCK.mu_star_per_h);

targetFinalNV=double(RLOCK.LOCK.spatialValidation.finalMeanNV);
targetEarlyG=double(RLOCK.LOCK.spatialValidation.earlyGrowth_per_h);

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

fprintf('\n============================================================\n');
fprintf('SI-M3 FINAL EVENT-DRIVEN MACRO SSA NUMERICAL VALIDATION\n');
fprintf('r=%.8f h^-1; mu=%.8f h^-1; L=%d; Tend=%.0f h\n',r_h,mu_h,L_MAIN,TEND_H);
fprintf('Solver: direct-method event-driven SSA; no numerical dt\n');
fprintf('============================================================\n');

%% 3. A: NO-FERROPTOSIS ENSEMBLE REGRESSION
cfgReg=makeCfg(L_MAIN,R0,TEND_H,RECORD_DT_H,r_h,mu_h,0);
REG=runSSAEnsemble(cfgReg,N_REG,SEED_REG,poolReady,true);

[tReg,meanReg,loReg,hiReg]=ensembleCurve(REG.t_h,REG.NVcurve);
earlyG=fitLogSlope(tReg,meanReg,FIT_WINDOW_H);
meanFinalNV=mean(REG.finalNV);
relFinal=(meanFinalNV-targetFinalNV)/targetFinalNV;
relG=(earlyG-targetEarlyG)/targetEarlyG;
edgeReg=mean(REG.edgeTouch);
bookReg=all(REG.bookkeepingPASS);
regPASS=abs(relFinal)<=REG_REL_TOL && abs(relG)<=REG_REL_TOL && edgeReg==0 && bookReg;

fprintf('\n[A] No-F regression\n');
fprintf('locked final NV = %.3f; SSA = %.3f (%+.2f%%)\n',targetFinalNV,meanFinalNV,100*relFinal);
fprintf('locked early g  = %.7f; SSA = %.7f (%+.2f%%)\n',targetEarlyG,earlyG,100*relG);
fprintf('edge touch = %.4f; bookkeeping=%d; %s\n',edgeReg,bookReg,passfail(regPASS));

%% 4. B: SUPPORTIVE SOLVER-TRANSITION CONSISTENCY
xvalAvailable=isfile(oldDtFile);
xRows=[];

if xvalAvailable
    OLD=load(oldDtFile,'AUDIT','RES');
    oldChi=double(OLD.AUDIT.settings.chi(:).');
    oldDtGrid=double(OLD.AUDIT.settings.dtGrid_h(:).');
    oldRefDt=double(OLD.AUDIT.settings.dtReference_h);
    iOldRef=find(abs(oldDtGrid-oldRefDt)<1e-14,1);
    assert(~isempty(iOldRef),'Cannot identify old finest dt reference.');

    for ix=1:numel(CHI_XVAL)
        chi=CHI_XVAL(ix);
        iOldChi=find(abs(oldChi-chi)<1e-12,1);
        if isempty(iOldChi), continue; end

        cfg=makeCfg(L_MAIN,R0,TEND_H,RECORD_DT_H,r_h,mu_h,chi*r_h);
        E=runSSAEnsemble(cfg,N_XVAL,SEED_XVAL+ix*10000,poolReady,false);
        O=OLD.RES{iOldChi,iOldRef};

        ssaPext=mean(E.extinctV);
        oldPext=mean(double(O.extinctV));
        ssaTau=mean(E.tauV_R_h);
        oldTau=mean(double(O.tauV_R_h));
        ssaAUTC=mean(E.AUTC_NT);
        oldAUTC=mean(double(O.AUTC_NT));
        dAUTC=symRelScalar(ssaAUTC,oldAUTC);

        xRows=[xRows;chi,ssaPext,oldPext,ssaPext-oldPext,ssaTau,oldTau,dAUTC]; %#ok<AGROW>
    end
end

if ~isempty(xRows)
    TX=array2table(xRows,'VariableNames',{ ...
        'chi','SSA_Pext','oldFine_Pext','dPext','SSA_tauV_h','oldFine_tauV_h','context_dAUTC_NT'});
else
    TX=table();
    xvalAvailable=false;
end

%% 5. C: MONTE-CARLO PRECISION
MC=cell(numel(CHI_MC),1);
mcRows=zeros(numel(CHI_MC),11);
mcPASS=true;

for i=1:numel(CHI_MC)
    chi=CHI_MC(i);
    cfg=makeCfg(L_MAIN,R0,TEND_H,RECORD_DT_H,r_h,mu_h,chi*r_h);
    seedBase=SEED_MC+i*10000;

    E=runSSAEnsemble(cfg,N_MC_BASE,seedBase,poolReady,false);
    nCurrent=N_MC_BASE;

    while true
        relNV=relativeMCSE(E.AUTC_NV);
        relNT=relativeMCSE(E.AUTC_NT);
        [loP,hiP]=wilsonCI(sum(E.extinctV),numel(E.extinctV),0.05);
        halfP=(hiP-loP)/2;
        passNow=relNV<=MCSE_TOL && relNT<=MCSE_TOL && halfP<=EXT_WILSON_HALF_TOL && all(E.bookkeepingPASS);
        if passNow || nCurrent>=N_MC_MAX, break; end

        nAdd=min(MC_BATCH,N_MC_MAX-nCurrent);
        Eadd=runSSAEnsemble(cfg,nAdd,seedBase+nCurrent,poolReady,false);
        E=appendEnsembles(E,Eadd);
        nCurrent=nCurrent+nAdd;
    end

    relNV=relativeMCSE(E.AUTC_NV);
    relNT=relativeMCSE(E.AUTC_NT);
    pExt=mean(E.extinctV);
    [loP,hiP]=wilsonCI(sum(E.extinctV),numel(E.extinctV),0.05);
    halfP=(hiP-loP)/2;
    passRow=relNV<=MCSE_TOL && relNT<=MCSE_TOL && halfP<=EXT_WILSON_HALF_TOL && all(E.bookkeepingPASS);
    mcPASS=mcPASS && passRow;
    MC{i}=E;

    mcRows(i,:)=[chi,numel(E.AUTC_NV),mean(E.AUTC_NV),relNV,mean(E.AUTC_NT),relNT,pExt,loP,hiP,halfP,double(passRow)];
end

TMC=array2table(mcRows,'VariableNames',{ ...
    'chi','nRep','mean_AUTC_NV','relMCSE_AUTC_NV','mean_AUTC_NT','relMCSE_AUTC_NT', ...
    'PextV','PextV_lo','PextV_hi','PextV_halfWidth','PASS'});

%% 6. D: FINAL COMMON-SEED DOMAIN CONFIRMATION
DOMC=cell(numel(CHI_DOMAIN),1);
DOMR=cell(numel(CHI_DOMAIN),1);
domRows=zeros(numel(CHI_DOMAIN),10);
domainPASS=true;

for i=1:numel(CHI_DOMAIN)
    chi=CHI_DOMAIN(i);
    seedFamily=SEED_DOMAIN+i*10000;

    cfgC=makeCfg(L_MAIN,R0,TEND_H,RECORD_DT_H,r_h,mu_h,chi*r_h);
    cfgR=makeCfg(L_REF,R0,TEND_H,RECORD_DT_H,r_h,mu_h,chi*r_h);

    C=runSSAEnsemble(cfgC,N_DOMAIN,seedFamily,poolReady,false);
    F=runSSAEnsemble(cfgR,N_DOMAIN,seedFamily,poolReady,false);
    DOMC{i}=C; DOMR{i}=F;

    aNV=pairedSymRelAudit(C.AUTC_NV,F.AUTC_NV,N_BOOT,SEED_BOOT+i*1000+1);
    aNT=pairedSymRelAudit(C.AUTC_NT,F.AUTC_NT,N_BOOT,SEED_BOOT+i*1000+2);
    edgeC=mean(C.edgeTouch); edgeF=mean(F.edgeTouch);

    passRow=effectPass(aNV,DOMAIN_POINT_TOL,DOMAIN_CI_TOL) && ...
        effectPass(aNT,DOMAIN_POINT_TOL,DOMAIN_CI_TOL) && edgeC==0 && edgeF==0 && ...
        all(C.bookkeepingPASS) && all(F.bookkeepingPASS);
    domainPASS=domainPASS && passRow;

    domRows(i,:)=[chi,aNV.point,aNV.ci(1),aNV.ci(2),aNT.point,aNT.ci(1),aNT.ci(2),edgeC,edgeF,double(passRow)];
end

TDOM=array2table(domRows,'VariableNames',{ ...
    'chi','dAUTC_NV','dAUTC_NV_lo','dAUTC_NV_hi','dAUTC_NT','dAUTC_NT_lo','dAUTC_NT_hi', ...
    'edge_L61','edge_L71','PASS'});

%% 7. DECISION
primaryPASS=regPASS && mcPASS && domainPASS;
if primaryPASS
    decision='FINAL_EVENT_DRIVEN_SSA_VALIDATED';
else
    decision='REVIEW_FAILED_PRIMARY_VALIDATION_COMPONENT';
end

fprintf('\n================ SI-M3 SUMMARY ================\n');
fprintf('No-F regression : %s\n',passfail(regPASS));
fprintf('MC precision    : %s\n',passfail(mcPASS));
fprintf('Domain          : %s\n',passfail(domainPASS));
fprintf('Old solver check: supportive only; available=%d\n',xvalAvailable);
fprintf('OVERALL         : %s\n',passfail(primaryPASS));
fprintf('Decision        : %s\n',decision);
fprintf('================================================\n');

%% 8. FIGURE
fig=figure('Color','w','Position',[80 80 1380 820]);
tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

% A. no-F regression
ax=nexttile(tl,1);
plot(ax,tReg,meanReg,'LineWidth',1.6,'DisplayName','event-driven SSA'); hold(ax,'on');
if isfield(RLOCK,'tStar') && isfield(RLOCK,'muStar')
    plot(ax,double(RLOCK.tStar(:)),double(RLOCK.muStar(:)),'--','LineWidth',1.4,'DisplayName','locked reference');
else
    plot(ax,TEND_H,targetFinalNV,'ks','MarkerFaceColor','w','DisplayName','locked final-N_V target');
end
fill(ax,[tReg;flipud(tReg)],[loReg;flipud(hiReg)],[0.8 0.8 0.8], ...
    'FaceAlpha',0.18,'EdgeColor','none','HandleVisibility','off');
xlabel(ax,'Time (h)'); ylabel(ax,'Mean N_V');
title(ax,'a  No-ferroptosis ensemble regression');
legend(ax,'Location','northwest','Box','off'); grid(ax,'on'); box(ax,'off');

% B. solver transition
ax=nexttile(tl,2);
if xvalAvailable
    plot(ax,TX.chi,TX.SSA_Pext,'-o','LineWidth',1.4,'DisplayName','event-driven SSA'); hold(ax,'on');
    plot(ax,TX.chi,TX.oldFine_Pext,'--s','LineWidth',1.4,'DisplayName','previous fine-step');
    xlabel(ax,'\chi_F=\lambda_F/r'); ylabel(ax,'P(V extinction by 168 h)'); ylim(ax,[-0.03 1.03]);
    title(ax,'b  Supportive solver-transition consistency');
    legend(ax,'Location','best','Box','off'); grid(ax,'on'); box(ax,'off');
else
    axis(ax,'off'); text(ax,0.5,0.5,'Previous fine-step audit unavailable','HorizontalAlignment','center');
end

% C. MC precision
ax=nexttile(tl,3);
yyaxis(ax,'left');
plot(ax,CHI_MC,100*TMC.relMCSE_AUTC_NT,'-o','LineWidth',1.4); hold(ax,'on');
yline(ax,100*MCSE_TOL,'--','MCSE target');
ylabel(ax,'Relative MCSE of AUTC_{NT} (%)');
yyaxis(ax,'right');
errorbar(ax,CHI_MC,TMC.PextV,TMC.PextV-TMC.PextV_lo,TMC.PextV_hi-TMC.PextV,'-s','LineWidth',1.2);
ylabel(ax,'P(V extinction by 168 h)'); ylim(ax,[-0.03 1.03]);
xlabel(ax,'\chi_F=\lambda_F/r'); title(ax,'c  Monte-Carlo precision'); grid(ax,'on'); box(ax,'off');

% D. final domain confirmation
ax=nexttile(tl,4); hold(ax,'on');
e1lo=100*(TDOM.dAUTC_NV-TDOM.dAUTC_NV_lo); e1hi=100*(TDOM.dAUTC_NV_hi-TDOM.dAUTC_NV);
errorbar(ax,CHI_DOMAIN,100*TDOM.dAUTC_NV,e1lo,e1hi,'-o','LineWidth',1.3,'DisplayName','AUTC_{NV}');
e2lo=100*(TDOM.dAUTC_NT-TDOM.dAUTC_NT_lo); e2hi=100*(TDOM.dAUTC_NT_hi-TDOM.dAUTC_NT);
errorbar(ax,CHI_DOMAIN,100*TDOM.dAUTC_NT,e2lo,e2hi,'-s','LineWidth',1.3,'DisplayName','AUTC_{NT}');
yline(ax,100*DOMAIN_POINT_TOL,'k:','HandleVisibility','off'); yline(ax,-100*DOMAIN_POINT_TOL,'k:','HandleVisibility','off');
yline(ax,100*DOMAIN_CI_TOL,'k--','HandleVisibility','off'); yline(ax,-100*DOMAIN_CI_TOL,'k--','HandleVisibility','off');
xlabel(ax,'\chi_F=\lambda_F/r'); ylabel(ax,'L=61 vs L=71 effect (%)');
title(ax,'d  Final domain confirmation'); legend(ax,'Location','best','Box','off'); grid(ax,'on'); box(ax,'off');

sgtitle(tl,sprintf('Numerical validation of the event-driven macro SSA: %s',decision));

%% 9. SAVE
outDir=fullfile(scriptDir,'SI_M3_Macro_SSA_Numerical_Validation_FINAL_OUT');
if ~exist(outDir,'dir'), mkdir(outDir); end
outPNG=fullfile(outDir,'SI_M3_Macro_SSA_Numerical_Validation_FINAL.png');
outPDF=fullfile(outDir,'SI_M3_Macro_SSA_Numerical_Validation_FINAL.pdf');
outFIG=fullfile(outDir,'SI_M3_Macro_SSA_Numerical_Validation_FINAL.fig');
outMAT=fullfile(outDir,'SI_M3_Macro_SSA_Numerical_Validation_FINAL.mat');
outCSVreg=fullfile(outDir,'SI_M3_Regression.csv');
outCSVx=fullfile(outDir,'SI_M3_SolverTransition.csv');
outCSVmc=fullfile(outDir,'SI_M3_MC.csv');
outCSVdom=fullfile(outDir,'SI_M3_Domain.csv');

try
    axs=findall(fig,'Type','axes');
    for k=1:numel(axs), axtoolbar(axs(k),{}); end
catch
end

exportgraphics(fig,outPNG,'Resolution',600);
exportgraphics(fig,outPDF,'ContentType','vector');
savefig(fig,outFIG);

TREG=table(targetFinalNV,meanFinalNV,relFinal,targetEarlyG,earlyG,relG,edgeReg,double(regPASS), ...
    'VariableNames',{'locked_finalNV','SSA_finalNV','rel_finalNV','locked_early_g','SSA_early_g','rel_early_g','edgeTouchProb','PASS'});
writetable(TREG,outCSVreg);
if xvalAvailable, writetable(TX,outCSVx); end
writetable(TMC,outCSVmc); writetable(TDOM,outCSVdom);

SI_M3=struct();
SI_M3.r_h=r_h; SI_M3.mu_h=mu_h; SI_M3.L=L_MAIN; SI_M3.Tend_h=TEND_H;
SI_M3.regressionPASS=regPASS; SI_M3.mcPASS=mcPASS; SI_M3.domainPASS=domainPASS;
SI_M3.crossValidationSupportive=xvalAvailable; SI_M3.overallPASS=primaryPASS; SI_M3.decision=decision;
save(outMAT,'SI_M3','REG','MC','DOMC','DOMR','TREG','TX','TMC','TDOM','-v7.3');

fprintf('\nSaved SI-M3 outputs:\n  %s\n  %s\n  %s\n',outPNG,outPDF,outMAT);
end

%% LOCAL FUNCTIONS
function cfg=makeCfg(L,R0,Tend,recordDt,r,mu,lambda)
cfg=struct(); cfg.L=L; cfg.R0=R0; cfg.Tend_h=Tend; cfg.record_dt_h=recordDt;
cfg.r_h=r; cfg.mu_h=mu; cfg.seed=1; cfg.lambda_mode='scalar'; cfg.lambdaF_h=lambda;
cfg.snapshot_times_h=[]; cfg.return_final_state=false; cfg.compute_cavity=false;
cfg.store_event_log=false; cfg.verbose=false;
end

function E=runSSAEnsemble(cfg,nRep,seedBase,poolReady,storeCurves)
tmp=cell(nRep,1);
if poolReady
    parfor k=1:nRep
        c=cfg; c.seed=seedBase+k; O=macro_sim_core_hazard_SSA_FINAL(c); tmp{k}=reduceSSA(O,storeCurves);
    end
else
    for k=1:nRep
        c=cfg; c.seed=seedBase+k; O=macro_sim_core_hazard_SSA_FINAL(c); tmp{k}=reduceSSA(O,storeCurves);
    end
end
E=stackReduced(tmp,storeCurves);
end

function z=reduceSSA(O,storeCurves)
z.finalNV=double(O.NV(end)); z.finalNT=double(O.NT(end));
z.AUTC_NV=double(O.AUTC_NV); z.AUTC_NT=double(O.AUTC_NT);
z.extinctV=isfinite(O.viableExtinctionTime_h) && O.viableExtinctionTime_h<=O.cfg.Tend_h;
z.tauV_R_h=min(double(O.viableExtinctionTime_h),double(O.cfg.Tend_h));
z.edgeTouch=any(O.edgeOcc); z.bookkeepingPASS=logical(O.audit.bookkeepingPASS);
if storeCurves, z.t_h=double(O.t_h(:)); z.NV=double(O.NV(:)); else, z.t_h=[]; z.NV=[]; end
end

function E=stackReduced(tmp,storeCurves)
n=numel(tmp); f={'finalNV','finalNT','AUTC_NV','AUTC_NT','tauV_R_h'};
for j=1:numel(f), E.(f{j})=zeros(n,1); end
E.extinctV=false(n,1); E.edgeTouch=false(n,1); E.bookkeepingPASS=false(n,1);
if storeCurves, E.t_h=tmp{1}.t_h; E.NVcurve=zeros(numel(E.t_h),n); else, E.t_h=[]; E.NVcurve=[]; end
for k=1:n
    E.finalNV(k)=tmp{k}.finalNV; E.finalNT(k)=tmp{k}.finalNT; E.AUTC_NV(k)=tmp{k}.AUTC_NV; E.AUTC_NT(k)=tmp{k}.AUTC_NT;
    E.tauV_R_h(k)=tmp{k}.tauV_R_h; E.extinctV(k)=tmp{k}.extinctV; E.edgeTouch(k)=tmp{k}.edgeTouch; E.bookkeepingPASS(k)=tmp{k}.bookkeepingPASS;
    if storeCurves, E.NVcurve(:,k)=tmp{k}.NV; end
end
end

function E=appendEnsembles(A,B)
E=struct(); fields={'finalNV','finalNT','AUTC_NV','AUTC_NT','tauV_R_h','extinctV','edgeTouch','bookkeepingPASS'};
for k=1:numel(fields), f=fields{k}; E.(f)=[A.(f);B.(f)]; end
E.t_h=[]; E.NVcurve=[];
end

function [t,mu,lo,hi]=ensembleCurve(t,M)
mu=mean(M,2); lo=zeros(size(mu)); hi=zeros(size(mu));
for k=1:numel(mu), lo(k)=pctLocal(M(k,:),2.5); hi(k)=pctLocal(M(k,:),97.5); end
end

function g=fitLogSlope(t,mu,win)
mask=t>=win(1) & t<=win(2) & mu>0;
if nnz(mask)<3, g=NaN; return; end
p=polyfit(t(mask),log(mu(mask)),1); g=p(1);
end

function m=relativeMCSE(x)
x=double(x(:)); mu=mean(x); if abs(mu)<=eps, m=0; else, m=std(x,0)/sqrt(numel(x))/abs(mu); end
end

function [lo,hi]=wilsonCI(k,n,alpha)
if nargin<3, alpha=0.05; end
if abs(alpha-0.05)<1e-12, z=1.959963984540054; else, z=sqrt(2)*erfcinv(alpha); end
p=k/n; den=1+z^2/n; ctr=(p+z^2/(2*n))/den; half=z*sqrt(p*(1-p)/n+z^2/(4*n^2))/den;
lo=max(0,ctr-half); hi=min(1,ctr+half);
end

function A=pairedSymRelAudit(x,y,nBoot,seed)
x=double(x(:)); y=double(y(:)); assert(numel(x)==numel(y)); den=abs(x)+abs(y); d=zeros(size(x)); m=den>eps; d(m)=2*(x(m)-y(m))./den(m);
A.point=mean(d); A.ci=bootstrapMeanCI(d,nBoot,seed);
end

function ci=bootstrapMeanCI(d,nBoot,seed)
d=double(d(:)); rng(seed,'twister'); n=numel(d); b=zeros(nBoot,1);
for k=1:nBoot, idx=randi(n,n,1); b(k)=mean(d(idx)); end
ci=[pctLocal(b,2.5) pctLocal(b,97.5)];
end

function tf=effectPass(A,pointTol,ciTol)
tf=abs(A.point)<=pointTol && A.ci(1)>=-ciTol && A.ci(2)<=ciTol;
end

function d=symRelScalar(x,y)
den=abs(x)+abs(y); if den<=eps, d=0; else, d=2*(x-y)/den; end
end

function q=pctLocal(v,p)
v=sort(v(isfinite(v))); if isempty(v), q=NaN; return; end
pos=1+(numel(v)-1)*p/100; a=floor(pos); b=ceil(pos);
if a==b, q=v(a); else, q=v(a)+(pos-a)*(v(b)-v(a)); end
end

function s=passfail(tf)
if tf, s='PASS'; else, s='FAIL'; end
end
