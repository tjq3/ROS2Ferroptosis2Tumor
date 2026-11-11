%% Fig4C_SI_FINAL_COLORFIX.m
% =========================================================================
% Fig. 4C SUPPLEMENTARY-INFORMATION FINAL
%
% This script generates ONLY the supporting analyses needed to validate the
% Fig. 4C main-text mechanism without over-expanding the SI.
%
% Run Fig4C_MAIN_FINAL.m first.
%
% SI FIGURES
% ----------
% SFig. C1 | Relative-timing reciprocity and finite-horizon boundary audit
%   - S-shift, T-shift and symmetric counter-shift
%   - source-overlap manipulation
%   - Delta P(commitment), Delta RMST
%   - T-only boundary audit
% Purpose: show that strict T-S coincidence is only a secondary mechanism,
% while shifting T can introduce absolute-time / finite-horizon effects.
%
% SFig. C2 | Priming robustness and route architecture
%   - RMST across T-only / Original / Flat / Random / Aligned
%   - Joint-required fraction across the same counterfactuals
%   - route composition
%   - rank-based overlap association audit
% Purpose: show that removing or randomizing S timing barely changes the
% population fate architecture even though realized overlap predicts fate.
%
% SFig. C3 | Critical-priming validation
%   - lambda*_orig vs lambda*_flat
%   - distribution of Delta lambda*
%   - population distribution of lambda*
%   - route-boundary consistency audit
%   - constant-input S-only asymptotic benchmark
%
% INPUTS
% ------
% 1) Fig4C_MAIN_FINAL_data.mat produced by Fig4C_MAIN_FINAL.m
% 2) original Fig4B_fateRoutes_data.mat source bank
% =========================================================================

clear; close all; clc;

%% ============================== Configuration =============================
MAIN_DATA_FILE = ''; % optional full path
MAIN_DATA_NAME = 'Fig4C_MAIN_FINAL_data.mat';

SOURCE_DATA_FILE = ''; % optional override for Fig4B source MAT
SOURCE_DATA_NAME = 'Fig4B_fateRoutes_data.mat';

N_BOOT = 1000;
ALPHA_CI = 0.05;

lagNorm = [-2 -1.5 -1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1 1.5 2];

MAX_SCATTER_POINTS = 2500;
SEED_BOOT = 20261101;
SEED_SCATTER = 20261102;

SAVE_FIG=true;
SAVE_CSV=true;
SAVE_DATA=true;
EXPORT_DPI=600;

%% ================================= Style =================================
set(groot,'defaultFigureColor','w');
set(groot,'defaultAxesFontName','Arial');
set(groot,'defaultTextFontName','Arial');
set(groot,'defaultAxesFontSize',7);
set(groot,'defaultTextFontSize',7);
set(groot,'defaultAxesLineWidth',0.70);
set(groot,'defaultAxesTickDir','out');

colCF=[0.18 0.43 0.70];
colGF=[0.82 0.31 0.24];

% Reciprocity-panel palettes:
% keep the biological condition encoded by hue
%   CIR + Full      -> blue family
%   GammaOU + Full  -> red/orange family
% and distinguish S-shift / T-shift / Symmetric by line style + shade.
methodColsCF=[ ...
    0.18 0.43 0.70; ...  % S-shift: base blue
    0.43 0.63 0.80; ...  % T-shift: lighter blue
    0.08 0.28 0.50];     % Symmetric: darker blue

methodColsGF=[ ...
    0.82 0.31 0.24; ...  % S-shift: base red
    0.90 0.55 0.47; ...  % T-shift: lighter red/orange
    0.62 0.18 0.13];     % Symmetric: darker red

colJoint=[0.93 0.66 0.20];
colRes=[0.62 0.62 0.62];

routeColors=[ ...
    0.24 0.50 0.74; ...
    0.28 0.65 0.50; ...
    0.50 0.34 0.67; ...
    0.93 0.66 0.20; ...
    0.78 0.79 0.80];

%% ================================= Load ==================================
scriptDir=fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir=pwd; end

if isempty(MAIN_DATA_FILE)
    mainFile=resolveMainFile(scriptDir,MAIN_DATA_NAME);
else
    mainFile=MAIN_DATA_FILE;
end
if isempty(mainFile) || ~isfile(mainFile)
    error(['Cannot find %s. Run Fig4C_MAIN_FINAL.m first, then set ' ...
           'MAIN_DATA_FILE if needed.'],MAIN_DATA_NAME);
end

M=load(mainFile);
required={'MAIN','PRIM_CF','PRIM_GF','LAMBDA_CF','LAMBDA_GF'};
for k=1:numel(required)
    if ~isfield(M,required{k}), error('Main data missing %s.',required{k}); end
end

MAIN=M.MAIN;
PRIM_CF=M.PRIM_CF;
PRIM_GF=M.PRIM_GF;
LAMBDA_CF=M.LAMBDA_CF;
LAMBDA_GF=M.LAMBDA_GF;

meso=MAIN.meso;
dt=MAIN.dt;
Tobs=MAIN.Tobs;
tauL=MAIN.tauL;

% Locate locked Fig.4B source bank.
if ~isempty(SOURCE_DATA_FILE)
    sourceFile=SOURCE_DATA_FILE;
elseif isfield(MAIN,'sourceFile') && isfile(MAIN.sourceFile)
    sourceFile=MAIN.sourceFile;
else
    sourceFile=resolveInputFile(scriptDir,SOURCE_DATA_NAME);
end
if isempty(sourceFile) || ~isfile(sourceFile)
    error('Cannot find the locked Fig4B source MAT.');
end

D=load(sourceFile);
IT_C=D.inputs.IT_C;
IT_G=D.inputs.IT_G;
IS_F=D.inputs.IS_F;

idxT=MAIN.idxT(:);
idxS=MAIN.idxS(:);
nTD=MAIN.nTD;
nSD=MAIN.nSD;
nPairs=MAIN.nPairs;
nBins=size(IT_C,2);

TC=IT_C(idxT,:);
TG=IT_G(idxT,:);
SF=IS_F(idxS,:);

SFlat=repmat(mean(SF,2),1,nBins);
SRandom=rowCircularShift(SF,MAIN.randomShift);
SAlignCF=rowCircularShift(SF,MAIN.bestShiftCF);
SAlignGF=rowCircularShift(SF,MAIN.bestShiftGF);

outDir=fullfile(scriptDir,'Fig4C_SI_FINAL_output');
if ~exist(outDir,'dir'), mkdir(outDir); end

fprintf('===============================================================\n');
fprintf('Fig. 4C SI FINAL\n');
fprintf('===============================================================\n');
fprintf('Main data : %s\n',mainFile);
fprintf('Source MAT: %s\n',sourceFile);
fprintf('nPairs=%d | nTD=%d | nSD=%d | bootstrap=%d\n', ...
    nPairs,nTD,nSD,N_BOOT);

%% ========================================================================
% SFIG C1. RELATIVE-TIMING RECIPROCITY + BOUNDARY AUDIT
% ========================================================================
lagMin=lagNorm*tauL;
lagBins=round(lagMin/dt);
lagMin=lagBins*dt;
lagNorm=lagMin/tauL;
[~,iZero]=min(abs(lagBins));

methodNames={'S-shift','T-shift','Symmetric'};
nMethod=3;
nLag=numel(lagBins);

REC_CF=initializeReciprocity(nPairs,nLag,nMethod);
REC_GF=initializeReciprocity(nPairs,nLag,nMethod);

fprintf('\nSFig C1 | reciprocity scan...\n');
for j=1:nLag
    d=lagBins(j);
    tHalf=-floor(d/2);
    sHalf= ceil(d/2);

    % CF
    REC_CF=fillReciprocity(REC_CF,j,1,TC,rowCircularShift(SAlignCF,d),meso,dt,Tobs);
    REC_CF=fillReciprocity(REC_CF,j,2,rowCircularShift(TC,-d),SAlignCF,meso,dt,Tobs);
    REC_CF=fillReciprocity(REC_CF,j,3,rowCircularShift(TC,tHalf), ...
        rowCircularShift(SAlignCF,sHalf),meso,dt,Tobs);

    % GF
    REC_GF=fillReciprocity(REC_GF,j,1,TG,rowCircularShift(SAlignGF,d),meso,dt,Tobs);
    REC_GF=fillReciprocity(REC_GF,j,2,rowCircularShift(TG,-d),SAlignGF,meso,dt,Tobs);
    REC_GF=fillReciprocity(REC_GF,j,3,rowCircularShift(TG,tHalf), ...
        rowCircularShift(SAlignGF,sHalf),meso,dt,Tobs);
end

BOOT_REC_CF=bootstrapReciprocity(REC_CF,idxT,idxS,nTD,nSD,Tobs, ...
    N_BOOT,ALPHA_CI,SEED_BOOT,iZero);
BOOT_REC_GF=bootstrapReciprocity(REC_GF,idxT,idxS,nTD,nSD,Tobs, ...
    N_BOOT,ALPHA_CI,SEED_BOOT+1,iZero);

methodLS={'-','--','-.'};
methodMk={'o','s','d'};

fig1=figure('Units','centimeters','Position',[1 1 18.3 10.6],'Color','w');
tl=tiledlayout(fig1,2,4,'TileSpacing','compact','Padding','compact');

for row=1:2
    if row==1
        X=REC_CF; B=BOOT_REC_CF; ttl='CIR + Full';
        methodCols=methodColsCF;
    else
        X=REC_GF; B=BOOT_REC_GF; ttl='GammaOU + Full';
        methodCols=methodColsGF;
    end

    ax=nexttile(tl,(row-1)*4+1); hold(ax,'on');
    for m=1:nMethod
        plot(ax,lagNorm,X.overlap(:,m),methodLS{m},'Color',methodCols(m,:), ...
            'Marker',methodMk{m},'MarkerFaceColor','w','LineWidth',1.05,'MarkerSize',3.8);
    end
    xline(ax,0,':','Color',[0.65 0.65 0.65]);
    xlabel(ax,'\delta/\tau_L'); ylabel(ax,'Mean source overlap');
    title(ax,[ttl ' | overlap'],'FontWeight','normal'); formatAxes(ax);
    if row==1
        legend(ax,methodNames,'Box','off','Location','best','FontSize',5.7);
    end

    ax=nexttile(tl,(row-1)*4+2); hold(ax,'on');
    for m=1:nMethod
        y=X.Pcommit(:,m)-X.Pcommit(iZero,m);
        ci=squeeze(B.P_delta_CI(:,m,:));
        plotWithBand(ax,lagNorm,y,ci,methodCols(m,:),methodLS{m},methodMk{m});
    end
    yline(ax,0,':','Color',[0.65 0.65 0.65]); xline(ax,0,':','Color',[0.75 0.75 0.75]);
    xlabel(ax,'\delta/\tau_L'); ylabel(ax,'\DeltaP(commitment)');
    title(ax,'Combined-fate effect','FontWeight','normal'); formatAxes(ax);

    ax=nexttile(tl,(row-1)*4+3); hold(ax,'on');
    for m=1:nMethod
        y=X.RMST(:,m)-X.RMST(iZero,m);
        ci=squeeze(B.RMST_delta_CI(:,m,:));
        plotWithBand(ax,lagNorm,y,ci,methodCols(m,:),methodLS{m},methodMk{m});
    end
    yline(ax,0,':','Color',[0.65 0.65 0.65]); xline(ax,0,':','Color',[0.75 0.75 0.75]);
    xlabel(ax,'\delta/\tau_L'); ylabel(ax,'RMST(\delta)-RMST(0) (min)');
    title(ax,'Timing effect','FontWeight','normal'); formatAxes(ax);

    ax=nexttile(tl,(row-1)*4+4); hold(ax,'on');
    for m=1:nMethod
        plot(ax,lagNorm,X.P_Tonly(:,m),methodLS{m},'Color',methodCols(m,:), ...
            'Marker',methodMk{m},'MarkerFaceColor','w','LineWidth',1.05,'MarkerSize',3.8);
    end
    xline(ax,0,':','Color',[0.75 0.75 0.75]);
    xlabel(ax,'\delta/\tau_L'); ylabel(ax,'P(T-only commitment)');
    title(ax,'Finite-horizon boundary audit','FontWeight','normal'); formatAxes(ax);
end

saveFigureSet(fig1,fullfile(outDir,'SFig4C1_reciprocity_boundary_FINAL'),SAVE_FIG,EXPORT_DPI);

%% ========================================================================
% SFIG C2. PRIMING ROBUSTNESS + ROUTE ARCHITECTURE
% ========================================================================
variantNames=MAIN.variantNames;
nVariant=numel(variantNames);

BOOT_PRIM_CF=bootstrapVariants(PRIM_CF,idxT,idxS,nTD,nSD,Tobs, ...
    N_BOOT,ALPHA_CI,SEED_BOOT+10);
BOOT_PRIM_GF=bootstrapVariants(PRIM_GF,idxT,idxS,nTD,nSD,Tobs, ...
    N_BOOT,ALPHA_CI,SEED_BOOT+11);

% Rank-based overlap association audit.
overlapRefCF=normalizedOverlapRows(positiveExcess(TC),positiveExcess(SF));
overlapRefGF=normalizedOverlapRows(positiveExcess(TG),positiveExcess(SF));

REF_CF=evaluateCondition(TC,SF,meso,dt,Tobs,false);
REF_GF=evaluateCondition(TG,SF,meso,dt,Tobs,false);

eligCF=(REF_CF.FPT_T>Tobs)&(REF_CF.FPT_S>Tobs);
eligGF=(REF_GF.FPT_T>Tobs)&(REF_GF.FPT_S>Tobs);

qCF=equalCountBins(overlapRefCF,eligCF,4,'ascend');
qGF=equalCountBins(overlapRefGF,eligGF,4,'ascend');

pOverlapCF=nan(1,4); pOverlapGF=nan(1,4);
for k=1:4
    ii=eligCF&qCF==k; pOverlapCF(k)=mean(REF_CF.routeCode(ii)==4);
    ii=eligGF&qGF==k; pOverlapGF(k)=mean(REF_GF.routeCode(ii)==4);
end

fig2=figure('Units','centimeters','Position',[1 1 18.3 10.5],'Color','w');
tl=tiledlayout(fig2,2,3,'TileSpacing','compact','Padding','compact');

x=1:nVariant;
ax=nexttile(tl,1); hold(ax,'on');
errorbar(ax,x-0.08,PRIM_CF.RMST, ...
    PRIM_CF.RMST-BOOT_PRIM_CF.RMST_CI(:,1)',BOOT_PRIM_CF.RMST_CI(:,2)'-PRIM_CF.RMST, ...
    'o','Color',colCF,'MarkerFaceColor','w','LineWidth',1,'CapSize',3);
errorbar(ax,x+0.08,PRIM_GF.RMST, ...
    PRIM_GF.RMST-BOOT_PRIM_GF.RMST_CI(:,1)',BOOT_PRIM_GF.RMST_CI(:,2)'-PRIM_GF.RMST, ...
    's','Color',colGF,'MarkerFaceColor','w','LineWidth',1,'CapSize',3);
setVariantAxis(ax,variantNames); ylabel(ax,'RMST_{60} (min)');
title(ax,'Priming counterfactual: timing','FontWeight','normal');
legend(ax,{'CIR + Full','GammaOU + Full'},'Box','off','Location','best','FontSize',5.8);
formatAxes(ax);

ax=nexttile(tl,2); hold(ax,'on');
errorbar(ax,x-0.08,PRIM_CF.jointFrac, ...
    PRIM_CF.jointFrac-BOOT_PRIM_CF.Joint_CI(:,1)',BOOT_PRIM_CF.Joint_CI(:,2)'-PRIM_CF.jointFrac, ...
    'o','Color',colCF,'MarkerFaceColor','w','LineWidth',1,'CapSize',3);
errorbar(ax,x+0.08,PRIM_GF.jointFrac, ...
    PRIM_GF.jointFrac-BOOT_PRIM_GF.Joint_CI(:,1)',BOOT_PRIM_GF.Joint_CI(:,2)'-PRIM_GF.jointFrac, ...
    's','Color',colGF,'MarkerFaceColor','w','LineWidth',1,'CapSize',3);
setVariantAxis(ax,variantNames); ylabel(ax,'Joint-required fraction');
title(ax,'Priming counterfactual: joint route','FontWeight','normal'); formatAxes(ax);

ax=nexttile(tl,3); hold(ax,'on');
plot(ax,1:4,pOverlapCF,'o-','Color',colCF,'MarkerFaceColor','w','LineWidth',1.2);
plot(ax,1:4,pOverlapGF,'s--','Color',colGF,'MarkerFaceColor','w','LineWidth',1.2);
set(ax,'XTick',1:4,'XTickLabel',{'Q1','Q2','Q3','Q4'});
xlabel(ax,'Realized overlap rank quartile');
ylabel(ax,'P(Joint-required | neither alone)');
ylim(ax,[0 1]);
title(ax,'Overlap predicts fate but is not the dominant cause','FontWeight','normal');
formatAxes(ax);

% route composition CF/GF
for row=1:2
    if row==1, X=PRIM_CF; ttl='CIR + Full'; tile=4;
    else, X=PRIM_GF; ttl='GammaOU + Full'; tile=5; end
    frac=zeros(nVariant,5);
    for v=1:nVariant
        rc=X.route(:,v);
        frac(v,:)=arrayfun(@(z)mean(rc==z),1:5);
    end
    ax=nexttile(tl,tile); hold(ax,'on');
    h=bar(ax,1:nVariant,frac,'stacked','BarWidth',0.72);
    for z=1:5
        h(z).FaceColor=routeColors(z,:); h(z).EdgeColor='none';
    end
    setVariantAxis(ax,variantNames);
    ylabel(ax,'Population fraction'); ylim(ax,[0 1]);
    title(ax,[ttl ' | fate-route composition'],'FontWeight','normal');
    formatAxes(ax);
end
legend(nexttile(tl,5),{'T-sufficient','S-sufficient','Redundant','Joint-required','Resistant'}, ...
    'Box','off','Location','southoutside','NumColumns',2,'FontSize',5.5);

% Original vs Flat S effective support at T-only peak.
SUP_CF=supportAtTemporalPeak(TC,SF,SFlat,meso,dt,Tobs);
SUP_GF=supportAtTemporalPeak(TG,SF,SFlat,meso,dt,Tobs);

ax=nexttile(tl,6); hold(ax,'on');
rng(SEED_SCATTER,'twister');
ii=find(SUP_CF.eligible);
if numel(ii)>MAX_SCATTER_POINTS, ii=ii(randperm(numel(ii),MAX_SCATTER_POINTS)); end
scatter(ax,SUP_CF.original(ii),SUP_CF.flat(ii),7,colCF,'filled','MarkerFaceAlpha',0.20);
jj=find(SUP_GF.eligible);
if numel(jj)>MAX_SCATTER_POINTS, jj=jj(randperm(numel(jj),MAX_SCATTER_POINTS)); end
scatter(ax,SUP_GF.original(jj),SUP_GF.flat(jj),7,colGF,'filled','MarkerFaceAlpha',0.18);
mx=max([SUP_CF.original(ii);SUP_CF.flat(ii);SUP_GF.original(jj);SUP_GF.flat(jj)]);
plot(ax,[0 mx],[0 mx],':','Color',[0.20 0.20 0.20]);
xlabel(ax,'Original-S support at T peak');
ylabel(ax,'Flat-S support at T peak');
title(ax,'Removing S timing preserves effective support','FontWeight','normal');
formatAxes(ax);

saveFigureSet(fig2,fullfile(outDir,'SFig4C2_priming_robustness_FINAL'),SAVE_FIG,EXPORT_DPI);

%% ========================================================================
% SFIG C3. CRITICAL-PRIMING VALIDATION
% ========================================================================
lambdaScrit = meso.rhoL*meso.ML / ...
    ((1-meso.ML)*meso.omegaS*mean(SF(:)));

fig3=figure('Units','centimeters','Position',[1 1 18.3 10.5],'Color','w');
tl=tiledlayout(fig3,2,3,'TileSpacing','compact','Padding','compact');

rng(SEED_SCATTER+10,'twister');

% orig vs flat, CF
ax=nexttile(tl,1); hold(ax,'on');
plotLambdaOrigFlat(ax,LAMBDA_CF,colCF,MAX_SCATTER_POINTS);
title(ax,sprintf('CIR + Full | r = %.3f',LAMBDA_CF.corrOrigFlat),'FontWeight','normal');
formatAxes(ax);

% orig vs flat, GF
ax=nexttile(tl,2); hold(ax,'on');
plotLambdaOrigFlat(ax,LAMBDA_GF,colGF,MAX_SCATTER_POINTS);
title(ax,sprintf('GammaOU + Full | r = %.3f',LAMBDA_GF.corrOrigFlat),'FontWeight','normal');
formatAxes(ax);

% delta lambda CDF
ax=nexttile(tl,3); hold(ax,'on');
dCF=LAMBDA_CF.lambdaFlat-LAMBDA_CF.lambdaOriginal;
dGF=LAMBDA_GF.lambdaFlat-LAMBDA_GF.lambdaOriginal;
plotEmpiricalCDF(ax,dCF(LAMBDA_CF.Tinsufficient & isfinite(dCF)),colCF,'-');
plotEmpiricalCDF(ax,dGF(LAMBDA_GF.Tinsufficient & isfinite(dGF)),colGF,'--');
xline(ax,0,':','Color',[0.15 0.15 0.15]);
xlabel(ax,'\lambda_{i,flat}^*-\lambda_{i,orig}^*');
ylabel(ax,'Cumulative probability');
title(ax,'Residual effect of S temporal structure','FontWeight','normal');
legend(ax,{'CIR + Full','GammaOU + Full'},'Box','off','Location','best','FontSize',5.8);
formatAxes(ax);

% lambda population CDF + theoretical S-only asymptote
ax=nexttile(tl,4); hold(ax,'on');
plotEmpiricalCDF(ax,LAMBDA_CF.lambdaOriginal( ...
    LAMBDA_CF.Tinsufficient & isfinite(LAMBDA_CF.lambdaOriginal)),colCF,'-');
plotEmpiricalCDF(ax,LAMBDA_GF.lambdaOriginal( ...
    LAMBDA_GF.Tinsufficient & isfinite(LAMBDA_GF.lambdaOriginal)),colGF,'--');
xline(ax,1,':','Color',[0.15 0.15 0.15]);
xline(ax,lambdaScrit,'-.','Color',[0.45 0.45 0.45]);
xlabel(ax,'Critical priming requirement, \lambda_i^*');
ylabel(ax,'Cumulative probability');
title(ax,sprintf('Population demand; S-only limit \\approx %.2f',lambdaScrit), ...
    'FontWeight','normal');
legend(ax,{'CIR + Full','GammaOU + Full','Observed \lambda=1','S-only constant-input limit'}, ...
    'Box','off','Location','southeast','FontSize',5.3);
formatAxes(ax);

% route-boundary audit CF/GF
ax=nexttile(tl,5); hold(ax,'on');
plotRouteAudit(ax,LAMBDA_CF);
title(ax,'CIR + Full | route-boundary audit','FontWeight','normal'); formatAxes(ax);

ax=nexttile(tl,6); hold(ax,'on');
plotRouteAudit(ax,LAMBDA_GF);
title(ax,'GammaOU + Full | route-boundary audit','FontWeight','normal'); formatAxes(ax);

saveFigureSet(fig3,fullfile(outDir,'SFig4C3_lambda_validation_FINAL'),SAVE_FIG,EXPORT_DPI);

%% =============================== Export ===================================
recRows=[];
for c=1:2
    if c==1, X=REC_CF; nm='CF'; else, X=REC_GF; nm='GF'; end
    for m=1:nMethod
        for j=1:nLag
            recRows=[recRows;table(string(nm),string(methodNames{m}), ...
                lagNorm(j),lagMin(j),X.overlap(j,m),X.Pcommit(j,m), ...
                X.RMST(j,m),X.jointFrac(j,m),X.P_Tonly(j,m),X.P_Sonly(j,m), ...
                'VariableNames',{'Condition','Manipulation','LagOverTauL','LagMin', ...
                'MeanOverlap','Pcommit','RMST','JointFrac','P_Tonly','P_Sonly'})]; %#ok<AGROW>
        end
    end
end

siSummary=table( ...
    ["CF";"GF"], ...
    [LAMBDA_CF.corrOrigFlat;LAMBDA_GF.corrOrigFlat], ...
    [LAMBDA_CF.spearmanMarginLambda;LAMBDA_GF.spearmanMarginLambda], ...
    [LAMBDA_CF.auditOriginal;LAMBDA_GF.auditOriginal], ...
    [LAMBDA_CF.auditFlat;LAMBDA_GF.auditFlat], ...
    repmat(lambdaScrit,2,1), ...
    'VariableNames',{'Condition','CorrLambdaOrigFlat','SpearmanMarginLambda', ...
    'RouteAuditOriginal','RouteAuditFlat','ConstantInputSOnlyLambdaCrit'});

if SAVE_CSV
    writetable(recRows,fullfile(outDir,'SFig4C1_reciprocity_data.csv'));
    writetable(siSummary,fullfile(outDir,'SFig4C3_lambda_validation_summary.csv'));
end

if SAVE_DATA
    save(fullfile(outDir,'Fig4C_SI_FINAL_data.mat'), ...
        'REC_CF','REC_GF','BOOT_REC_CF','BOOT_REC_GF', ...
        'BOOT_PRIM_CF','BOOT_PRIM_GF','SUP_CF','SUP_GF', ...
        'pOverlapCF','pOverlapGF','lambdaScrit','-v7.3');
end

fprintf('\n====================== FIG.4C SI FINAL ==========================\n');
fprintf('S-only commitment max across reciprocity scan:\n');
fprintf('  CF %.4g | GF %.4g\n',max(REC_CF.P_Sonly(:)),max(REC_GF.P_Sonly(:)));
fprintf('Constant-input S-only lambda benchmark = %.4f\n',lambdaScrit);
fprintf('Outputs: %s\n',outDir);
fprintf('=================================================================\n');

%% =========================================================================
% Local functions
% =========================================================================
function path=resolveMainFile(scriptDir,name)
    candidates={ ...
        fullfile(pwd,name), ...
        fullfile(scriptDir,name), ...
        fullfile(pwd,'Fig4C_MAIN_FINAL_output',name), ...
        fullfile(scriptDir,'Fig4C_MAIN_FINAL_output',name)};
    path='';
    for k=1:numel(candidates)
        if isfile(candidates{k}), path=candidates{k}; return; end
    end
end

function path=resolveInputFile(scriptDir,name)
    candidates = {fullfile(pwd,name),fullfile(scriptDir,name), ...
        fullfile(scriptDir,'Fig4B_fate_routes_final',name), ...
        fullfile(fileparts(scriptDir),name)};
    path='';
    for k=1:numel(candidates)
        if isfile(candidates{k}), path=candidates{k}; return; end
    end
end

function E=positiveExcess(I)
    E=max(I-mean(I,2),0);
end

function s=normalizedOverlapRows(A,B)
    den=sqrt(sum(A.^2,2).*sum(B.^2,2));
    s=sum(A.*B,2)./max(den,eps);
    s(den<=eps)=0;
end

function Y=rowCircularShift(X,shiftBins)
    [n,m]=size(X);
    s=mod(round(shiftBins(:)),m);
    if numel(s)==1, s=repmat(s,n,1); end
    rows=repmat((1:n)',1,m);
    cols=mod(repmat(0:m-1,n,1)-s,m)+1;
    Y=X(sub2ind([n,m],rows,cols));
end

function X=initializeReciprocity(nPairs,nLag,nMethod)
    X.FPT_TS=inf(nPairs,nLag,nMethod);
    X.route=zeros(nPairs,nLag,nMethod,'uint8');
    X.Pcommit=zeros(nLag,nMethod);
    X.RMST=zeros(nLag,nMethod);
    X.jointFrac=zeros(nLag,nMethod);
    X.P_Tonly=zeros(nLag,nMethod);
    X.P_Sonly=zeros(nLag,nMethod);
    X.overlap=zeros(nLag,nMethod);
end

function X=fillReciprocity(X,j,m,T,S,meso,dt,Tobs)
    R=evaluateCondition(T,S,meso,dt,Tobs,false);
    X.FPT_TS(:,j,m)=R.FPT_TS;
    X.route(:,j,m)=uint8(R.routeCode);
    X.Pcommit(j,m)=R.Pcommit;
    X.RMST(j,m)=R.RMST;
    X.jointFrac(j,m)=R.routeFrac(4);
    X.P_Tonly(j,m)=mean(R.FPT_T<=Tobs);
    X.P_Sonly(j,m)=mean(R.FPT_S<=Tobs);
    X.overlap(j,m)=mean(normalizedOverlapRows(positiveExcess(T),positiveExcess(S)));
end

function B=bootstrapReciprocity(X,idxT,idxS,nTD,nSD,Tobs,nBoot,alpha,seed,iZero)
    rng(seed,'twister');
    nLag=size(X.FPT_TS,2);
    nMethod=size(X.FPT_TS,3);
    bP=nan(nBoot,nLag,nMethod);
    bR=nan(nBoot,nLag,nMethod);
    bJ=nan(nBoot,nLag,nMethod);

    for b=1:nBoot
        w=twoWayWeights(idxT,idxS,nTD,nSD);
        for m=1:nMethod
            for j=1:nLag
                f=X.FPT_TS(:,j,m); rc=X.route(:,j,m);
                bP(b,j,m)=weightedMean(f<=Tobs,w);
                bR(b,j,m)=weightedMean(min(f,Tobs),w);
                bJ(b,j,m)=weightedMean(rc==4,w);
            end
        end
    end

    B.P_delta_CI=ci3D(bP-bP(:,iZero,:),alpha);
    B.RMST_delta_CI=ci3D(bR-bR(:,iZero,:),alpha);
    B.Joint_delta_CI=ci3D(bJ-bJ(:,iZero,:),alpha);
end

function B=bootstrapVariants(X,idxT,idxS,nTD,nSD,Tobs,nBoot,alpha,seed)
    rng(seed,'twister');
    nVariant=size(X.FPT_TS,2);
    bP=nan(nBoot,nVariant); bR=nan(nBoot,nVariant); bJ=nan(nBoot,nVariant);
    for b=1:nBoot
        w=twoWayWeights(idxT,idxS,nTD,nSD);
        for v=1:nVariant
            f=X.FPT_TS(:,v); rc=X.route(:,v);
            bP(b,v)=weightedMean(f<=Tobs,w);
            bR(b,v)=weightedMean(min(f,Tobs),w);
            bJ(b,v)=weightedMean(rc==4,w);
        end
    end
    B.P_CI=ci2D(bP,alpha);
    B.RMST_CI=ci2D(bR,alpha);
    B.Joint_CI=ci2D(bJ,alpha);
end

function S=supportAtTemporalPeak(T,Sorig,Sflat,meso,dt,Tobs)
    Rorig=evaluateCondition(T,Sorig,meso,dt,Tobs,true);
    Rflat=evaluateCondition(T,Sflat,meso,dt,Tobs,true);
    LT=Rorig.L_T;
    [~,peak]=max(LT,[],2);
    n=size(T,1);
    lin=sub2ind(size(LT),(1:n)',peak);
    S.original=Rorig.L_TS(lin)-LT(lin);
    S.flat=Rflat.L_TS(lin)-LT(lin);
    S.eligible=(Rorig.FPT_T>Tobs)&(Rorig.FPT_S>Tobs);
end

function R=evaluateCondition(T,S,meso,dt,Tobs,returnHist)
    if nargin<6, returnHist=false; end
    [FT,LT]=propagateLPOExact(meso.omegaT*T,meso.rhoL,meso.ML,meso.L0,dt,Tobs,returnHist);
    [FS,LS]=propagateLPOExact(meso.omegaS*S,meso.rhoL,meso.ML,meso.L0,dt,Tobs,returnHist);
    [FTS,LTS]=propagateLPOExact(meso.omegaT*T+meso.omegaS*S, ...
        meso.rhoL,meso.ML,meso.L0,dt,Tobs,returnHist);

    cT=FT<=Tobs; cS=FS<=Tobs; cTS=FTS<=Tobs;
    rc=zeros(size(FTS),'uint8');
    rc(cT&~cS)=1; rc(~cT&cS)=2; rc(cT&cS)=3;
    rc(~cT&~cS&cTS)=4; rc(~cTS)=5;

    R.FPT_T=FT; R.FPT_S=FS; R.FPT_TS=FTS;
    R.routeCode=rc;
    R.routeFrac=arrayfun(@(z)mean(rc==z),1:5);
    R.Pcommit=mean(cTS);
    R.RMST=mean(min(FTS,Tobs));
    if returnHist
        R.L_T=LT; R.L_S=LS; R.L_TS=LTS;
    end
end

function [FPT,Lhist]=propagateLPOExact(U,rho,ML,L0,dt,Tobs,returnHistory)
    if nargin<7, returnHistory=false; end
    [n,nBins]=size(U);
    L=L0*ones(n,1); FPT=inf(n,1);
    if returnHistory, Lhist=zeros(n,nBins); else, Lhist=[]; end

    for b=1:nBins
        Ub=max(U(:,b),0); a=Ub+rho; Lnext=L;
        ok=a>eps;
        if any(ok)
            av=a(ok); uv=Ub(ok); Lss=uv./av;
            Lnext(ok)=Lss+(L(ok)-Lss).*exp(-av*dt);
        end
        Lnext=min(max(Lnext,0),1);

        cross=isinf(FPT)&(L<ML)&(Lnext>=ML);
        if any(cross)
            ii=find(cross); av=a(ii); uv=Ub(ii); Lss=uv./av;
            ratio=(ML-Lss)./(L(ii)-Lss);
            ratio=max(min(ratio,1),realmin);
            tc=-log(ratio)./av; tc=max(min(tc,dt),0);
            FPT(ii)=(b-1)*dt+tc;
        end
        L=Lnext;
        if returnHistory, Lhist(:,b)=L; end
    end
    FPT(FPT>Tobs)=inf;
end

function plotLambdaOrigFlat(ax,X,cMain,maxPoints)
    valid=X.Tinsufficient & X.resolvedOriginal & X.resolvedFlat & ...
        isfinite(X.lambdaOriginal)&isfinite(X.lambdaFlat);
    ii=find(valid);
    if numel(ii)>maxPoints, ii=ii(randperm(numel(ii),maxPoints)); end
    scatter(ax,X.lambdaOriginal(ii),X.lambdaFlat(ii),7,cMain, ...
        'filled','MarkerFaceAlpha',0.18);
    mx=max([X.lambdaOriginal(ii);X.lambdaFlat(ii);1]);
    plot(ax,[0 mx],[0 mx],':','Color',[0.10 0.10 0.10]);
    xline(ax,1,':','Color',[0.65 0.65 0.65]);
    yline(ax,1,':','Color',[0.65 0.65 0.65]);
    xlabel(ax,'\lambda_{i,orig}^*'); ylabel(ax,'\lambda_{i,flat}^*');
    xlim(ax,[0 mx]); ylim(ax,[0 mx]); axis(ax,'square');
end

function plotRouteAudit(ax,X)
    comparable=X.Tinsufficient & X.SinsufficientOriginal;
    joint=comparable & X.routeOriginal==4 & isfinite(X.lambdaOriginal);
    resist=comparable & X.routeOriginal==5 & isfinite(X.lambdaOriginal);
    plotEmpiricalCDF(ax,X.lambdaOriginal(joint),colLocalJoint(),'-');
    plotEmpiricalCDF(ax,X.lambdaOriginal(resist),[0.62 0.62 0.62],'--');
    xline(ax,1,':','Color',[0.12 0.12 0.12]);
    xlabel(ax,'Critical priming requirement, \lambda_i^*');
    ylabel(ax,'Cumulative probability');
    legend(ax,{'Joint-required','Resistant','Observed \lambda=1'}, ...
        'Box','off','Location','best','FontSize',5.5);
end

function c=colLocalJoint()
    c=[0.93 0.66 0.20];
end

function plotEmpiricalCDF(ax,x,col,ls)
    x=x(isfinite(x)); x=sort(x(:));
    if isempty(x), return; end
    y=(1:numel(x))'/numel(x);
    plot(ax,x,y,ls,'Color',col,'LineWidth',1.2);
end

function bins=equalCountBins(score,eligible,nBin,order)
    bins=zeros(size(score));
    idx=find(eligible&isfinite(score));
    if isempty(idx), return; end
    [~,ord]=sort(score(idx),order); idx=idx(ord);
    N=numel(idx); edges=round(linspace(0,N,nBin+1));
    for q=1:nBin
        lo=edges(q)+1; hi=edges(q+1);
        if hi>=lo, bins(idx(lo:hi))=q; end
    end
end

function w=twoWayWeights(idxT,idxS,nTD,nSD)
    cT=accumarray(randi(nTD,nTD,1),1,[nTD 1]);
    cS=accumarray(randi(nSD,nSD,1),1,[nSD 1]);
    w=double(cT(idxT)).*double(cS(idxS));
end

function m=weightedMean(x,w)
    x=double(x); ok=isfinite(x)&isfinite(w)&w>0;
    if ~any(ok), m=NaN; else, m=sum(w(ok).*x(ok))/sum(w(ok)); end
end

function CI=ci2D(X,alpha)
    CI=nan(size(X,2),2);
    for j=1:size(X,2)
        x=X(:,j); x=x(isfinite(x));
        if ~isempty(x), CI(j,:)=quantileLocal(x,[alpha/2 1-alpha/2]); end
    end
end

function CI=ci3D(X,alpha)
    nLag=size(X,2); nMethod=size(X,3);
    CI=nan(nLag,nMethod,2);
    for j=1:nLag
        for m=1:nMethod
            x=X(:,j,m); x=x(isfinite(x));
            if ~isempty(x)
                q=quantileLocal(x,[alpha/2 1-alpha/2]);
                CI(j,m,1)=q(1); CI(j,m,2)=q(2);
            end
        end
    end
end

function q=quantileLocal(x,p)
    x=x(isfinite(x)); x=sort(x(:));
    n=numel(x); p=p(:).'; q=nan(size(p));
    if n==0, return; end
    for j=1:numel(p)
        pp=min(max(p(j),0),1); pos=1+pp*(n-1);
        lo=floor(pos); hi=ceil(pos);
        if lo==hi, q(j)=x(lo);
        else, q(j)=x(lo)+(pos-lo)*(x(hi)-x(lo)); end
    end
end

function plotWithBand(ax,x,y,CI,col,ls,mk)
    x=x(:).'; y=y(:).';
    lo=CI(:,1).'; hi=CI(:,2).';
    valid=isfinite(x)&isfinite(y)&isfinite(lo)&isfinite(hi);
    if any(valid)
        patch(ax,[x(valid) fliplr(x(valid))],[lo(valid) fliplr(hi(valid))], ...
            col,'FaceAlpha',0.08,'EdgeColor','none','HandleVisibility','off');
    end
    plot(ax,x,y,ls,'Color',col,'LineWidth',1.05,'Marker',mk, ...
        'MarkerFaceColor','w','MarkerSize',3.8);
end

function setVariantAxis(ax,names)
    set(ax,'XTick',1:numel(names),'XTickLabel', ...
        {'T only','Original','Flat S','Random','Aligned'}, ...
        'XTickLabelRotation',25);
    xlim(ax,[0.5 numel(names)+0.5]);
end

function saveFigureSet(fig,stem,saveFlag,dpi)
    if ~saveFlag, return; end
    exportgraphics(fig,[stem '.png'],'Resolution',dpi);
    exportgraphics(fig,[stem '.tif'],'Resolution',dpi);
    exportgraphics(fig,[stem '.pdf'],'ContentType','vector');
    savefig(fig,[stem '.fig']);
end

function formatAxes(ax)
    set(ax,'Box','on','LineWidth',0.70,'TickDir','out','Layer','top', ...
        'FontName','Arial','FontSize',7, ...
        'XColor',[0.08 0.08 0.08],'YColor',[0.08 0.08 0.08]);
end
