function LPO_RepairTimeAndFateRoutes(DATA_FILE)
% LPO_RepairTimeAndFateRoutes
% FINAL main-text implementation for Fig.4D.
% Identical microscopic histories; only cellular repair/memory state varies.

if nargin < 1 || isempty(DATA_FILE)
    scriptDir = fileparts(mfilename('fullpath'));
    DATA_FILE = fullfile(scriptDir,'Fig4B_fate_routes_final','Fig4B_fateRoutes_data.mat');
end
close all; clc;

%% FINAL SETTINGS
tauScan = [3 5 7.5 10 12.5 15 17.5 20 21 22.5 25 30];
N_BOOT = 1000;
ALPHA_CI = 0.05;
SEED_BOOT = 20261021;
BOOT_CHUNK = 50;
TAU_HIGH_REPAIR = 5;
TAU_REFERENCE = 10;
TAU_LOW_REPAIR = 20;
N_MARGIN_BINS = 10;
MIN_MARGIN_COUNT = 100;
SAVE_DPI = 600;

%% STYLE
set(groot,'defaultFigureColor','w');
set(groot,'defaultAxesFontName','Arial');
set(groot,'defaultTextFontName','Arial');
set(groot,'defaultAxesFontSize',7.5);
set(groot,'defaultTextFontSize',7.5);
set(groot,'defaultAxesLineWidth',0.75);
set(groot,'defaultAxesTickDir','out');
set(groot,'defaultAxesTickLength',[0.014 0.014]);
set(groot,'defaultAxesBox','off');
set(groot,'defaultAxesLayer','top');

% Restrained, color-blind-aware palette for publication figures.
colCF=[0.20 0.44 0.69];     % muted blue
colGF=[0.80 0.32 0.24];     % muted vermilion
routeColors=[ ...
    0.26 0.50 0.72; ...     % T-sufficient
    0.32 0.66 0.52; ...     % S-sufficient
    0.48 0.36 0.66; ...     % Redundant
    0.93 0.66 0.18; ...     % Joint-required
    0.79 0.80 0.81];        % Resistant
routeNames={'T-sufficient','S-sufficient','Redundant','Joint-required','Resistant'};

%% LOAD
if ~isfile(DATA_FILE), error('Cannot find input MAT file:\n%s',DATA_FILE); end
D=load(DATA_FILE);
if ~isfield(D,'meta') || ~isfield(D,'inputs'), error('Input MAT lacks meta/inputs.'); end
req={'IT_C','IT_G','IS_U','IS_F','idxT_pair','idxS_pair'};
for k=1:numel(req)
    if ~isfield(D.inputs,req{k}), error('Missing D.inputs.%s',req{k}); end
end
IT_C=D.inputs.IT_C; IT_G=D.inputs.IT_G; IS_U=D.inputs.IS_U; IS_F=D.inputs.IS_F;
idxT=D.inputs.idxT_pair(:); idxS=D.inputs.idxS_pair(:);
meso0=D.meta.meso; dt=D.meta.dtBin_min; Tobs=D.meta.Tobs_min;
if isfield(meso0,'tauL_min'), tauRef=meso0.tauL_min; else, tauRef=1/meso0.rhoL; end
nTD=size(IT_C,1); nSD=size(IS_F,1); nPairs=numel(idxT);

fprintf('Loading: %s\n',DATA_FILE);
fprintf('nTD / nSD / nPairs = %d / %d / %d\n',nTD,nSD,nPairs);
fprintf('dt / Tobs = %.2f / %.1f min\n',dt,Tobs);
fprintf('M_L / omegaT / omegaS = %.4f / %.6f / %.6f\n',meso0.ML,meso0.omegaT,meso0.omegaS);
fprintf('Locked reference tau_L = %.2f min\n',tauRef);

%% FIXED CELL PAIRING
ITCp=IT_C(idxT,:); ITGp=IT_G(idxT,:); ISUp=IS_U(idxS,:); ISFp=IS_F(idxS,:);
Tbank={ITCp,ITGp,ITCp,ITGp};
Sbank={ISUp,ISUp,ISFp,ISFp};
condNames={'CU','GU','CF','GF'};
condLabels={'CIR + Uniform','GammaOU + Uniform','CIR + Full','GammaOU + Full'};
nCond=4; nTau=numel(tauScan);
[~,iHigh]=min(abs(tauScan-TAU_HIGH_REPAIR));
[~,iRef]=min(abs(tauScan-TAU_REFERENCE));
[~,iLow]=min(abs(tauScan-TAU_LOW_REPAIR));

%% DECODER SCAN
R=cell(nCond,nTau);
fprintf('\nEvaluating final decoder-state scan...\n');
for j=1:nTau
    meso=meso0;
    meso.tauL_min=tauScan(j);
    meso.rhoL=1/tauScan(j);   % ONLY changing cellular decoder state
    fprintf('  tau_L = %5.1f min | rho_L = %.5f min^-1\n',tauScan(j),meso.rhoL);
    for c=1:nCond
        R{c,j}=evaluateFateRoutesWithPeak(Tbank{c},Sbank{c},meso,dt,Tobs);
        R{c,j}.name=condNames{c}; R{c,j}.label=condLabels{c};
    end
end

Pcommit=zeros(nCond,nTau); RMST=zeros(nCond,nTau); routeFrac=zeros(nCond,nTau,5);
for c=1:nCond
    for j=1:nTau
        Pcommit(c,j)=R{c,j}.Pcross;
        RMST(c,j)=R{c,j}.RMST;
        routeFrac(c,j,:)=reshape(R{c,j}.routeFrac,1,1,5);
    end
end

%% BASELINE AUDIT
fprintf('\nBaseline reproduction audit @ tau_L = %.1f min\n',tauScan(iRef));
if isfield(D,'R') && numel(D.R)>=4
    for c=1:4
        dP=Pcommit(c,iRef)-D.R(c).Pcross60;
        dR=RMST(c,iRef)-D.R(c).RMST60;
        dRoute=max(abs(squeeze(routeFrac(c,iRef,:)).'-D.R(c).routeFrac));
        fprintf('%s: dP=%+.2e, dRMST=%+.2e, max dRoute=%.2e\n',condNames{c},dP,dR,dRoute);
    end
end

%% MONOTONICITY AUDIT
maxViolation=0;
for c=1:nCond
    for j=1:nTau-1
        f1=capInf(R{c,j}.FPT_TS,Tobs+dt); f2=capInf(R{c,j+1}.FPT_TS,Tobs+dt);
        maxViolation=max(maxViolation,max(f2-f1));
    end
end
fprintf('Maximum pathwise monotonicity violation = %.3e min\n',maxViolation);
if maxViolation>1e-8, error('Repair-state monotonicity audit failed.'); end

%% APPROXIMATE S-ONLY CRITICAL MEMORY
meanISF=mean(IS_F(:),'omitnan');
meanUS=meso0.omegaS*meanISF;
tauScrit=meso0.ML/((1-meso0.ML)*meanUS);
fprintf('Approximate S-only critical tau_L = %.2f min\n',tauScrit);

%% BOOTSTRAP
fprintf('\nTwo-source cluster bootstrap: %d replicates...\n',N_BOOT);
BOOT=bootstrapDecoderScanCluster(R,idxT,idxS,nTD,nSD,Tobs,N_BOOT,ALPHA_CI,SEED_BOOT,BOOT_CHUNK);

%% TRANSITIONS
TRANS=struct();
for c=[3 4]
    rcHigh=R{c,iHigh}.routeCode; rcRef=R{c,iRef}.routeCode; rcLow=R{c,iLow}.routeCode;
    TRANS(c).high_to_ref=routeTransitionMatrix(rcHigh,rcRef,5);
    TRANS(c).ref_to_low=routeTransitionMatrix(rcRef,rcLow,5);
    TRANS(c).fracSwitch_high_to_ref=mean(rcHigh~=rcRef);
    TRANS(c).fracSwitch_ref_to_low=mean(rcRef~=rcLow);
end

%% FATE-MARGIN MAP
marginEdges=linspace(0,meso0.ML,N_MARGIN_BINS+1);
marginCenters=0.5*(marginEdges(1:end-1)+marginEdges(2:end));
MARGIN=struct();
for c=[3 4]
    MARGIN(c).prob=nan(nTau,N_MARGIN_BINS);
    MARGIN(c).count=zeros(nTau,N_MARGIN_BINS);
    for j=1:nTau
        Rc=R{c,j}; neither=~Rc.crossT & ~Rc.crossS;
        margin=max(meso0.ML-Rc.peakL_T,0);
        for b=1:N_MARGIN_BINS
            if b<N_MARGIN_BINS
                inBin=margin>=marginEdges(b) & margin<marginEdges(b+1);
            else
                inBin=margin>=marginEdges(b) & margin<=marginEdges(b+1);
            end
            use=neither & inBin;
            MARGIN(c).count(j,b)=sum(use);
            if any(use), MARGIN(c).prob(j,b)=mean(Rc.routeCode(use)==4); end
        end
    end
    MARGIN(c).probMasked=MARGIN(c).prob;
    MARGIN(c).probMasked(MARGIN(c).count<MIN_MARGIN_COUNT)=NaN;
end

%% SUMMARY
fprintf('\n================ FIG.4D FINAL SUMMARY ================\n');
for c=[3 4]
    fprintf('\n%s | %s\n',condNames{c},condLabels{c});
    for jj=[iHigh iRef iLow]
        ff=squeeze(routeFrac(c,jj,:));
        fprintf('tau=%.1f: T %.3f | S %.3f | Red %.3f | Joint %.3f | Res %.3f\n',tauScan(jj),ff(1),ff(2),ff(3),ff(4),ff(5));
    end
    fprintf('Same-cell switch 10->20 min = %.3f\n',TRANS(c).fracSwitch_ref_to_low);
    [mxJ,ij]=max(squeeze(routeFrac(c,:,4)));
    fprintf('Max Joint-required = %.3f at tau_L = %.1f min\n',mxJ,tauScan(ij));
end
fprintf('=======================================================\n');

%% OUTPUT DIR
scriptDir=fileparts(mfilename('fullpath')); if isempty(scriptDir), scriptDir=pwd; end
outDir=fullfile(scriptDir,'Fig4D_FINAL_outputs'); if ~exist(outDir,'dir'), mkdir(outDir); end

%% MAIN FIGURE
% -------------------------------------------------------------------------
% Publication layout note:
% All top/bottom axes within a column are assigned EXACTLY the same Position.
% The legend and colorbar are positioned manually so that MATLAB cannot
% shrink one of the data axes (the previous 'southoutside' legend caused
% the lower-left axis to become shorter than the upper-left axis).
% -------------------------------------------------------------------------
fig=figure('Units','centimeters','Position',[1 1 18.3 9.2],'Color','w');
set(fig,'Renderer','painters');

% ---- locked axes geometry [left bottom width height] ----
P11=[0.060 0.575 0.250 0.305];
P12=[0.060 0.165 0.250 0.305];
P21=[0.385 0.575 0.225 0.305];
P22=[0.385 0.165 0.225 0.305];
P31=[0.695 0.575 0.225 0.305];
P32=[0.695 0.165 0.225 0.305];

ax11=axes(fig,'Units','normalized','Position',P11); hold(ax11,'on');
ax12=axes(fig,'Units','normalized','Position',P12); hold(ax12,'on');
ax21=axes(fig,'Units','normalized','Position',P21);
ax22=axes(fig,'Units','normalized','Position',P22);
ax31=axes(fig,'Units','normalized','Position',P31);
ax32=axes(fig,'Units','normalized','Position',P32);
tauRatio=tauScan/tauRef;

% ======================== D1: fate-route rewiring ========================
plotRouteComposition(ax11,tauRatio,squeeze(routeFrac(3,:,:)),routeColors);
title(ax11,'CIR + Full','FontWeight','normal','FontSize',8.3);
ylabel(ax11,'Population fraction');
xline(ax11,1,':','Color',[0.25 0.25 0.25],'LineWidth',0.8);
ylim(ax11,[0 1]); xlim(ax11,[min(tauRatio) max(tauRatio)]);
set(ax11,'XTickLabel',[]); styleAxes(ax11);

plotRouteComposition(ax12,tauRatio,squeeze(routeFrac(4,:,:)),routeColors);
title(ax12,'GammaOU + Full','FontWeight','normal','FontSize',8.3);
xlabel(ax12,'Lipid-damage memory, \tau_L/\tau_{L,0}');
ylabel(ax12,'Population fraction');
xline(ax12,1,':','Color',[0.25 0.25 0.25],'LineWidth',0.8);
ylim(ax12,[0 1]); xlim(ax12,[min(tauRatio) max(tauRatio)]); styleAxes(ax12);

annotation(fig,'textbox',[P11(1) 0.915 P11(3) 0.035], ...
    'String','Cell-state-dependent fate-route rewiring', ...
    'EdgeColor','none','HorizontalAlignment','center', ...
    'FontSize',8.8,'FontWeight','normal');

% Manual legend placement: do NOT use 'southoutside', because it resizes ax12.
dummy=gobjects(5,1);
for rr=1:5
    dummy(rr)=plot(ax11,nan,nan,'s','MarkerSize',5, ...
        'MarkerFaceColor',routeColors(rr,:),'MarkerEdgeColor','none');
end
lgd=legend(ax12,dummy,routeNames,'NumColumns',2,'Box','off','FontSize',5.8);
set(lgd,'Units','normalized','Position',[0.067 0.020 0.235 0.075]);

% CRITICAL: restore locked axes positions after creating the legend.
set(ax11,'Units','normalized','Position',P11);
set(ax12,'Units','normalized','Position',P12);

% ====================== D2: same-cell switching ==========================
plotTransitionMatrix(ax21,TRANS(3).ref_to_low);
title(ax21,'CIR + Full: 10 \rightarrow 20 min','FontWeight','normal','FontSize',8.0);
plotTransitionMatrix(ax22,TRANS(4).ref_to_low);
title(ax22,'GammaOU + Full: 10 \rightarrow 20 min','FontWeight','normal','FontSize',8.0);
annotation(fig,'textbox',[P21(1) 0.915 P21(3) 0.035], ...
    'String','Same-cell fate-route switching', ...
    'EdgeColor','none','HorizontalAlignment','center', ...
    'FontSize',8.8,'FontWeight','normal');

% Restore exact container sizes (axis square only changes the inner plot box).
set(ax21,'Units','normalized','Position',P21);
set(ax22,'Units','normalized','Position',P22);

% ======================== D3: fate boundary ==============================
plotMarginMap(ax31,marginCenters,tauRatio,MARGIN(3).probMasked);
title(ax31,'CIR + Full','FontWeight','normal','FontSize',8.0);
ylabel(ax31,'\tau_L/\tau_{L,0}'); set(ax31,'XTickLabel',[]);
yline(ax31,1,':','Color',[0.2 0.2 0.2],'LineWidth',0.8);

plotMarginMap(ax32,marginCenters,tauRatio,MARGIN(4).probMasked);
title(ax32,'GammaOU + Full','FontWeight','normal','FontSize',8.0);
xlabel(ax32,'Temporal fate margin, M_L - max L_T');
ylabel(ax32,'\tau_L/\tau_{L,0}');
yline(ax32,1,':','Color',[0.2 0.2 0.2],'LineWidth',0.8);

annotation(fig,'textbox',[P31(1) 0.915 P31(3) 0.035], ...
    'String','Repair state reshapes the fate boundary', ...
    'EdgeColor','none','HorizontalAlignment','center', ...
    'FontSize',8.8,'FontWeight','normal');

% Manual colorbar placement; then restore both right-column axes positions.
cb=colorbar(ax31);
set(cb,'Units','normalized','Position',[0.934 0.590 0.010 0.275]);
cb.Label.String='Joint-route probability'; cb.FontSize=6.5;
set(ax31,'Units','normalized','Position',P31);
set(ax32,'Units','normalized','Position',P32);

% Final geometry lock: guarantees identical width/height for every vertical pair.
drawnow;
set(ax11,'Position',P11); set(ax12,'Position',P12);
set(ax21,'Position',P21); set(ax22,'Position',P22);
set(ax31,'Position',P31); set(ax32,'Position',P32);

annotation(fig,'textbox',[0.013 0.920 0.025 0.035], ...
    'String','D','EdgeColor','none','FontSize',12,'FontWeight','bold');

stem=fullfile(outDir,'Fig4D_CellStateDependentFateRewiring_FINAL_LAYOUTFIX');
exportgraphics(fig,[stem '.png'],'Resolution',SAVE_DPI);
exportgraphics(fig,[stem '.tif'],'Resolution',SAVE_DPI);
exportgraphics(fig,[stem '.pdf'],'ContentType','vector');
savefig(fig,[stem '.fig']);

%% CSV
rows={};
for c=1:nCond
    for j=1:nTau
        pCI=squeeze(BOOT.P_CI(c,j,:)); rCI=squeeze(BOOT.RMST_CI(c,j,:)); routeCI=squeeze(BOOT.route_CI(c,j,:,:));
        rows(end+1,:)={string(condNames{c}),string(condLabels{c}),tauScan(j),1/tauScan(j),tauRatio(j),Pcommit(c,j),pCI(1),pCI(2),RMST(c,j),rCI(1),rCI(2),routeFrac(c,j,1),routeCI(1,1),routeCI(1,2),routeFrac(c,j,2),routeCI(2,1),routeCI(2,2),routeFrac(c,j,3),routeCI(3,1),routeCI(3,2),routeFrac(c,j,4),routeCI(4,1),routeCI(4,2),routeFrac(c,j,5),routeCI(5,1),routeCI(5,2)};
    end
end
summaryTable=cell2table(rows,'VariableNames',{'Condition','ConditionLabel','tauL_min','rhoL_per_min','tauRatio','Pcommit','Pcommit_CI_L','Pcommit_CI_U','RMST60_min','RMST_CI_L','RMST_CI_U','T_sufficient','T_CI_L','T_CI_U','S_sufficient','S_CI_L','S_CI_U','Redundant','Red_CI_L','Red_CI_U','Joint_required','Joint_CI_L','Joint_CI_U','Resistant','Res_CI_L','Res_CI_U'});
writetable(summaryTable,fullfile(outDir,'Fig4D_CellStateDependentFateRewiring_FINAL_LAYOUTFIX_summary.csv'));

%% SAVE DATA
meta=struct(); meta.sourceFile=DATA_FILE; meta.tauScan_min=tauScan; meta.tauReference_min=tauRef; meta.tauScrit_approx_min=tauScrit; meta.MIN_MARGIN_COUNT=MIN_MARGIN_COUNT; meta.ML=meso0.ML; meta.omegaT=meso0.omegaT; meta.omegaS=meso0.omegaS; meta.dt_min=dt; meta.Tobs_min=Tobs;
save(fullfile(outDir,'Fig4D_CellStateDependentFateRewiring_FINAL_LAYOUTFIX_data.mat'),'meta','R','Pcommit','RMST','routeFrac','BOOT','TRANS','MARGIN','summaryTable','routeNames','routeColors','condNames','condLabels','-v7.3');
fprintf('\nSaved FINAL outputs to:\n%s\n',outDir);
end

function R=evaluateFateRoutesWithPeak(IT,IS,meso,dt,Tobs)
UT=meso.omegaT*IT; US=meso.omegaS*IS; UTS=UT+US;
[FPT_T,peakL_T]=propagateLPOExactPeak(UT,meso.rhoL,meso.ML,meso.L0,dt,Tobs);
[FPT_S,peakL_S]=propagateLPOExactPeak(US,meso.rhoL,meso.ML,meso.L0,dt,Tobs);
[FPT_TS,peakL_TS]=propagateLPOExactPeak(UTS,meso.rhoL,meso.ML,meso.L0,dt,Tobs);
crossT=FPT_T<=Tobs; crossS=FPT_S<=Tobs; crossTS=FPT_TS<=Tobs;
if any((crossT|crossS)&~crossTS), error('Pathwise ordering violation.'); end
routeCode=zeros(size(FPT_TS));
routeCode(crossT & ~crossS)=1; routeCode(~crossT & crossS)=2; routeCode(crossT & crossS)=3; routeCode(~crossT & ~crossS & crossTS)=4; routeCode(~crossTS)=5;
if any(routeCode==0), error('Unclassified route detected.'); end
R=struct(); R.FPT_T=FPT_T; R.FPT_S=FPT_S; R.FPT_TS=FPT_TS; R.peakL_T=peakL_T; R.peakL_S=peakL_S; R.peakL_TS=peakL_TS; R.crossT=crossT; R.crossS=crossS; R.crossTS=crossTS; R.routeCode=routeCode; R.routeFrac=arrayfun(@(r)mean(routeCode==r),1:5); R.Pcross=mean(crossTS); R.RMST=mean(min(FPT_TS,Tobs));
end

function [FPT,peakL]=propagateLPOExactPeak(U,rho,ML,L0,dt,Tobs)
[n,nBins]=size(U); L=L0*ones(n,1); FPT=inf(n,1); peakL=L;
for b=1:nBins
    Ub=max(U(:,b),0); a=Ub+rho; Lnext=L; active=a>eps;
    if any(active)
        aa=a(active); uu=Ub(active); Lss=uu./aa; E=exp(-aa*dt); Lnext(active)=Lss+(L(active)-Lss).*E;
    end
    Lnext=min(max(Lnext,0),1); cross=isinf(FPT)&(L<ML)&(Lnext>=ML);
    if any(cross)
        idx=find(cross); aa=a(idx); uu=Ub(idx); Lss=uu./aa; ratio=(ML-Lss)./(L(idx)-Lss); ratio=max(min(ratio,1),realmin); tau=-log(ratio)./aa; tau=max(min(tau,dt),0); FPT(idx)=(b-1)*dt+tau;
    end
    peakL=max(peakL,Lnext); L=Lnext;
end
FPT(FPT>Tobs)=inf;
end

function y=capInf(x,val), y=x; y(~isfinite(y))=val; end

function T=routeTransitionMatrix(routeFrom,routeTo,nRoute)
T=nan(nRoute,nRoute);
for i=1:nRoute
    from=routeFrom==i; denom=sum(from); if denom==0, continue; end
    for j=1:nRoute, T(i,j)=sum(from & routeTo==j)/denom; end
end
end

function BOOT=bootstrapDecoderScanCluster(R,idxT,idxS,nTD,nSD,Tobs,nBoot,alpha,seed,chunkSize)
rng(seed,'twister'); [nCond,nTau]=size(R); nPairs=numel(idxT); nMetric=7; X=cell(nCond,1);
for c=1:nCond
    M=zeros(nPairs,nTau*nMetric);
    for j=1:nTau
        col=(j-1)*nMetric; f=R{c,j}.FPT_TS; rc=R{c,j}.routeCode; M(:,col+1)=double(f<=Tobs); M(:,col+2)=min(f,Tobs);
        for rr=1:5, M(:,col+2+rr)=double(rc==rr); end
    end
    X{c}=M;
end
bootVals=nan(nBoot,nCond,nTau*nMetric);
for first=1:chunkSize:nBoot
    last=min(nBoot,first+chunkSize-1); nb=last-first+1; WT=zeros(nTD,nb); WS=zeros(nSD,nb);
    for bb=1:nb
        drawT=randi(nTD,nTD,1); drawS=randi(nSD,nSD,1); WT(:,bb)=accumarray(drawT,1,[nTD 1]); WS(:,bb)=accumarray(drawS,1,[nSD 1]);
    end
    W=WT(idxT,:).*WS(idxS,:); den=sum(W,1).'; good=den>0;
    for c=1:nCond
        V=W.'*X{c}; V(good,:)=V(good,:)./den(good); bootVals(first:last,c,:)=reshape(V,[nb 1 size(V,2)]);
    end
    fprintf('bootstrap %d / %d\n',last,nBoot);
end
BOOT.P_CI=zeros(nCond,nTau,2); BOOT.RMST_CI=zeros(nCond,nTau,2); BOOT.route_CI=zeros(nCond,nTau,5,2);
for c=1:nCond
    for j=1:nTau
        col=(j-1)*nMetric; BOOT.P_CI(c,j,:)=quantileLocal(squeeze(bootVals(:,c,col+1)),[alpha/2 1-alpha/2]); BOOT.RMST_CI(c,j,:)=quantileLocal(squeeze(bootVals(:,c,col+2)),[alpha/2 1-alpha/2]);
        for rr=1:5, BOOT.route_CI(c,j,rr,:)=quantileLocal(squeeze(bootVals(:,c,col+2+rr)),[alpha/2 1-alpha/2]); end
    end
end
end

function q=quantileLocal(x,p)
x=x(isfinite(x)); x=sort(x(:)); q=nan(size(p)); if isempty(x), return; end; n=numel(x);
for k=1:numel(p)
    z=1+(n-1)*p(k); lo=floor(z); hi=ceil(z); if lo==hi, q(k)=x(lo); else, q(k)=x(lo)+(z-lo)*(x(hi)-x(lo)); end
end
end

function plotRouteComposition(ax,x,F,colors)
hold(ax,'on'); if size(F,2)~=5, F=F.'; end; h=area(ax,x,F,'LineStyle','none');
for k=1:5, h(k).FaceColor=colors(k,:); h(k).FaceAlpha=0.96; end; styleAxes(ax);
end

function plotTransitionMatrix(ax,M)
imagesc(ax,M,[0 1]); axis(ax,'square'); short={'T','S','Red','Joint','Res'}; set(ax,'XTick',1:5,'YTick',1:5,'XTickLabel',short,'YTickLabel',short,'XTickLabelRotation',35); xlabel(ax,'Route after repair change'); ylabel(ax,'Route before repair change');
cmap=[linspace(0.95,0.18,256)' linspace(0.97,0.43,256)' linspace(1.00,0.70,256)']; colormap(ax,cmap);
for i=1:5
    for j=1:5
        if isfinite(M(i,j)) && M(i,j)>=0.045
            tc='k'; if M(i,j)>0.55, tc='w'; end; text(ax,j,i,sprintf('%.2f',M(i,j)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',5.7,'Color',tc);
        end
    end
end
styleAxes(ax);
end

function plotMarginMap(ax,x,y,Z)
h=imagesc(ax,x,y,Z);
set(ax,'YDir','normal');
h.AlphaData=isfinite(Z);
colormap(ax,publicationSequentialMap(256));
caxis(ax,[0 1]);
xlim(ax,[min(x) max(x)]);
ylim(ax,[min(y) max(y)]);
styleAxes(ax);
end

function cmap=publicationSequentialMap(n)
% Low-saturation blue -> cyan/green -> warm yellow.
anchors=[ ...
    0.15 0.22 0.55; ...
    0.18 0.50 0.72; ...
    0.24 0.68 0.62; ...
    0.73 0.78 0.36; ...
    0.96 0.82 0.25];
x=linspace(0,1,size(anchors,1));
xi=linspace(0,1,n);
cmap=zeros(n,3);
for k=1:3
    cmap(:,k)=interp1(x,anchors(:,k),xi,'pchip');
end
cmap=min(max(cmap,0),1);
end

function styleAxes(ax)
set(ax,'FontName','Arial','FontSize',7.5,'LineWidth',0.75, ...
    'TickDir','out','TickLength',[0.014 0.014], ...
    'Box','off','Layer','top');
end
