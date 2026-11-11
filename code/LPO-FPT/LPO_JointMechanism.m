%% Fig4C_MAIN_FINAL.m
% =========================================================================
% Fig. 4C MAIN-TEXT FINAL
% Persistent mitochondrial priming closes the residual temporal fate margin
%
% This script is the FINAL main-text implementation for Fig. 4C.
%
% Scientific logic
% ----------------
% Fig. 4B established a Joint-required population:
%   T-only does not cross the ferroptotic boundary,
%   S-only does not cross the boundary,
%   but T+S does.
%
% Mechanism-screening experiments then showed that:
%   1) removing the temporal structure of S (Flat-S) barely changes fate;
%   2) randomizing or maximizing T-S phase produces only a secondary effect;
%   3) Joint-required cells are concentrated near the T-only fate boundary.
%
% Fig. 4C therefore focuses on the mechanism supported by those tests:
%
%   temporal history -> residual fate margin -> critical mitochondrial
%   priming requirement -> ferroptotic commitment.
%
% The critical priming coordinate is
%
%   lambda_i* = inf{lambda >= 0 :
%                   T_F,i[omega_T I_T + lambda omega_S I_S] <= Tobs }.
%
% MAIN FIGURE CONTENT
% -------------------
% C1. Representative Joint-required cell:
%       T is burst-like; S is persistent; Original-S and Flat-S give nearly
%       identical LPO trajectories and both close the T-only fate gap.
%
% C2. Population counterfactual:
%       T-only vs Original-S vs Flat-S vs Random-phase-S vs Max-aligned-S.
%       Main outcome: P(commitment by 60 min), with two-source bootstrap CI.
%
% C3. Critical priming map:
%       temporal fate margin Delta_T = M_L - max L_T
%       versus lambda_i*, separately for CIR+Full and GammaOU+Full.
%
% INPUT
% -----
% Locked FINAL Fig. 4B source bank:
%       Fig4B_fateRoutes_data.mat
%
% IMPORTANT
% ---------
% This script does NOT regenerate TD-CME or SD-RDME trajectories.
% It uses the locked Fig. 4B source banks and the same mesoscopic parameters.
% =========================================================================

clear; close all; clc;

%% ============================== Configuration =============================
DATA_FILE = '';  % optional full path; leave empty for auto-detection
DATA_NAME = 'Fig4B_fateRoutes_data.mat';

N_BOOT = 1000;
ALPHA_CI = 0.05;

LAMBDA_MAX  = 8.0;
LAMBDA_TOL  = 1e-4;
LAMBDA_ITER = 20;

MAX_SCATTER_POINTS = 1800;

SEED_RANDOM_PHASE = 20261001;
SEED_BOOTSTRAP    = 20261002;
SEED_SCATTER      = 20261003;

SAVE_FIG  = true;
SAVE_CSV  = true;
SAVE_DATA = true;
EXPORT_DPI = 600;

%% ================================ Style ==================================
set(groot,'defaultFigureColor','w');
set(groot,'defaultAxesFontName','Arial');
set(groot,'defaultTextFontName','Arial');
set(groot,'defaultAxesFontSize',7);
set(groot,'defaultTextFontSize',7);
set(groot,'defaultAxesLineWidth',0.70);
set(groot,'defaultAxesTickDir','out');

colCF = [0.18 0.43 0.70];
colGF = [0.82 0.31 0.24];
colS  = [0.93 0.66 0.20];
colGray = [0.42 0.42 0.42];
colLight = [0.72 0.72 0.72];

%% ================================ Load ===================================
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir = pwd; end

if isempty(DATA_FILE)
    sourceFile = resolveInputFile(scriptDir,DATA_NAME);
else
    sourceFile = DATA_FILE;
end
if isempty(sourceFile) || ~isfile(sourceFile)
    error(['Cannot find %s. Set DATA_FILE or place the MAT file beside ' ...
           'this script.'],DATA_NAME);
end

D=load(sourceFile);

if ~isfield(D,'meta') || ~isfield(D,'inputs')
    error('The source MAT must contain meta and inputs.');
end
if ~isfield(D.inputs,'idxT_pair') || ~isfield(D.inputs,'idxS_pair')
    error('The source MAT must contain idxT_pair and idxS_pair.');
end

IT_C=D.inputs.IT_C;
IT_G=D.inputs.IT_G;
IS_F=D.inputs.IS_F;

idxT=D.inputs.idxT_pair(:);
idxS=D.inputs.idxS_pair(:);

meso=D.meta.meso;
dt=D.meta.dtBin_min;
Tobs=D.meta.Tobs_min;

if isfield(meso,'tauL_min')
    tauL=meso.tauL_min;
else
    tauL=1/meso.rhoL;
end

nTD=size(IT_C,1);
nSD=size(IS_F,1);
nPairs=numel(idxT);
nBins=size(IT_C,2);

assert(size(IT_G,2)==nBins && size(IS_F,2)==nBins, ...
    'TD/SD source banks do not share a common time grid.');
assert(abs(nBins*dt-Tobs)<1e-8,'nBins*dt and Tobs are inconsistent.');

TC=IT_C(idxT,:);
TG=IT_G(idxT,:);
SF=IS_F(idxS,:);

outDir=fullfile(scriptDir,'Fig4C_MAIN_FINAL_output');
if ~exist(outDir,'dir'), mkdir(outDir); end

fprintf('===============================================================\n');
fprintf('Fig. 4C MAIN FINAL\n');
fprintf('===============================================================\n');
fprintf('Source MAT            : %s\n',sourceFile);
fprintf('nTD / nSD / nPairs    : %d / %d / %d\n',nTD,nSD,nPairs);
fprintf('dt / Tobs / tau_L     : %.2f / %.1f / %.1f min\n',dt,Tobs,tauL);
fprintf('M_L                   : %.3f\n',meso.ML);
fprintf('omega_T / omega_S     : %.5f / %.5f min^-1\n',meso.omegaT,meso.omegaS);
fprintf('Bootstrap             : %d two-source replicates\n',N_BOOT);

%% ===================== Source counterfactual construction =================
% Flat-S: remove all temporal structure while preserving each paired cell's
% exact integrated spatial dose.
SFlat=repmat(mean(SF,2),1,nBins);

% Random-phase S: preserve the exact S path and dose, change only phase.
rng(SEED_RANDOM_PHASE,'twister');
randomShift=randi([0 nBins-1],nPairs,1);
SRandom=rowCircularShift(SF,randomShift);

% Max-aligned S: source-only alignment, independent of fate.
ETC=positiveExcess(TC);
ETG=positiveExcess(TG);
ESF=positiveExcess(SF);

[bestShiftCF,overlapRefCF,overlapBestCF]=maxOverlapShift(ETC,ESF);
[bestShiftGF,overlapRefGF,overlapBestGF]=maxOverlapShift(ETG,ESF);

SAlignCF=rowCircularShift(SF,bestShiftCF);
SAlignGF=rowCircularShift(SF,bestShiftGF);

fprintf('\nSpatial-dose preservation audit:\n');
fprintf('  Flat-S max error      = %.3e\n', ...
    max(abs(sum(SFlat,2)*dt-sum(SF,2)*dt)));
fprintf('  Random-phase max error= %.3e\n', ...
    max(abs(sum(SRandom,2)*dt-sum(SF,2)*dt)));
fprintf('Source overlap reference -> maximal alignment:\n');
fprintf('  CF %.4f -> %.4f\n',mean(overlapRefCF),mean(overlapBestCF));
fprintf('  GF %.4f -> %.4f\n',mean(overlapRefGF),mean(overlapBestGF));

%% ======================= Main counterfactual outcomes =====================
variantNames={'T-only','Original S','Flat-S','Random-phase S','Max-aligned S'};
nVariant=numel(variantNames);

S_C={zeros(size(SF)),SF,SFlat,SRandom,SAlignCF};
S_G={zeros(size(SF)),SF,SFlat,SRandom,SAlignGF};

PRIM_CF=initializeVariantResult(nPairs,nVariant);
PRIM_GF=initializeVariantResult(nPairs,nVariant);

for v=1:nVariant
    PRIM_CF=fillVariant(PRIM_CF,v,TC,S_C{v},meso,dt,Tobs);
    PRIM_GF=fillVariant(PRIM_GF,v,TG,S_G{v},meso,dt,Tobs);
end

fprintf('\nTwo-source bootstrap for main counterfactuals...\n');
BOOT_CF=bootstrapVariants(PRIM_CF,idxT,idxS,nTD,nSD,Tobs, ...
    N_BOOT,ALPHA_CI,SEED_BOOTSTRAP);
BOOT_GF=bootstrapVariants(PRIM_GF,idxT,idxS,nTD,nSD,Tobs, ...
    N_BOOT,ALPHA_CI,SEED_BOOTSTRAP+1);

%% ===================== Temporal fate margin + lambda_i* ==================
fprintf('\nComputing temporal fate margins...\n');
MARGIN_CF=temporalMarginOnly(TC,meso,dt,Tobs);
MARGIN_GF=temporalMarginOnly(TG,meso,dt,Tobs);

fprintf('Computing critical priming lambda_i* (Original and Flat-S)...\n');
LAMBDA_CF=criticalLambdaAnalysis(TC,SF,SFlat,MARGIN_CF,meso,dt,Tobs, ...
    LAMBDA_MAX,LAMBDA_TOL,LAMBDA_ITER);
LAMBDA_GF=criticalLambdaAnalysis(TG,SF,SFlat,MARGIN_GF,meso,dt,Tobs, ...
    LAMBDA_MAX,LAMBDA_TOL,LAMBDA_ITER);

fprintf('\nCritical-priming audits:\n');
fprintf('  CF resolved original/flat = %.4f / %.4f\n', ...
    mean(LAMBDA_CF.resolvedOriginal(LAMBDA_CF.Tinsufficient)), ...
    mean(LAMBDA_CF.resolvedFlat(LAMBDA_CF.Tinsufficient)));
fprintf('  GF resolved original/flat = %.4f / %.4f\n', ...
    mean(LAMBDA_GF.resolvedOriginal(LAMBDA_GF.Tinsufficient)), ...
    mean(LAMBDA_GF.resolvedFlat(LAMBDA_GF.Tinsufficient)));
fprintf('  CF route audit original/flat = %.4f / %.4f\n', ...
    LAMBDA_CF.auditOriginal,LAMBDA_CF.auditFlat);
fprintf('  GF route audit original/flat = %.4f / %.4f\n', ...
    LAMBDA_GF.auditOriginal,LAMBDA_GF.auditFlat);
fprintf('  CF corr(lambda*_orig,lambda*_flat)=%.4f | r_s(margin,lambda*)=%.4f\n', ...
    LAMBDA_CF.corrOrigFlat,LAMBDA_CF.spearmanMarginLambda);
fprintf('  GF corr(lambda*_orig,lambda*_flat)=%.4f | r_s(margin,lambda*)=%.4f\n', ...
    LAMBDA_GF.corrOrigFlat,LAMBDA_GF.spearmanMarginLambda);

fprintf('Bootstrap summary for lambda_i*...\n');
BOOT_LAM_CF=bootstrapLambdaSummary(LAMBDA_CF,idxT,idxS,nTD,nSD, ...
    N_BOOT,ALPHA_CI,SEED_BOOTSTRAP+10);
BOOT_LAM_GF=bootstrapLambdaSummary(LAMBDA_GF,idxT,idxS,nTD,nSD, ...
    N_BOOT,ALPHA_CI,SEED_BOOTSTRAP+11);

%% ======================= Representative priming cell =====================
% Prefer a cell that is Joint-required under Original-S AND remains
% Joint-required under Flat-S. Choose one near the median temporal margin.
origRouteCF=PRIM_CF.route(:,2);
flatRouteCF=PRIM_CF.route(:,3);
candidate=find(origRouteCF==4 & flatRouteCF==4);

if isempty(candidate)
    candidate=find(origRouteCF==4);
end
if isempty(candidate)
    error('No CF Joint-required cell found for the representative panel.');
end

mvals=MARGIN_CF.margin(candidate);
medm=median(mvals,'omitnan');
[~,kk]=min(abs(mvals-medm));
repIndex=candidate(kk);

[~,LT_rep]=propagateLPOExact(meso.omegaT*TC(repIndex,:), ...
    meso.rhoL,meso.ML,meso.L0,dt,Tobs,true);
[~,LOrig_rep]=propagateLPOExact( ...
    meso.omegaT*TC(repIndex,:)+meso.omegaS*SF(repIndex,:), ...
    meso.rhoL,meso.ML,meso.L0,dt,Tobs,true);
[~,LFlat_rep]=propagateLPOExact( ...
    meso.omegaT*TC(repIndex,:)+meso.omegaS*SFlat(repIndex,:), ...
    meso.rhoL,meso.ML,meso.L0,dt,Tobs,true);

t=(1:nBins)*dt;

%% =============================== MAIN FIGURE ==============================
% Layout:
% left block  : representative source histories + LPO trajectories
% middle      : population counterfactual commitment
% right block : critical priming maps for CF and GF

fig=figure('Units','centimeters','Position',[1 1 19.0 7.7],'Color','w');

ax1=axes(fig,'Position',[0.055 0.57 0.295 0.35]); hold(ax1,'on');
plot(ax1,t,TC(repIndex,:),'Color',colCF,'LineWidth',1.05);
plot(ax1,t,SF(repIndex,:),'Color',colS,'LineWidth',1.00);
plot(ax1,t,SFlat(repIndex,:),'--','Color',colGray,'LineWidth',1.00);
ylabel(ax1,'Normalized input');
title(ax1,'Persistent mitochondrial priming','FontWeight','normal','FontSize',8);
legend(ax1,{'Temporal input','Original spatial input','Flat spatial input'}, ...
    'Box','off','Location','northwest','FontSize',5.7);
formatAxes(ax1);
text(ax1,-0.18,1.08,'C','Units','normalized','FontWeight','bold','FontSize',12);

ax2=axes(fig,'Position',[0.055 0.13 0.295 0.32]); hold(ax2,'on');
plot(ax2,t,LT_rep,'Color',colCF,'LineWidth',1.15);
plot(ax2,t,LOrig_rep,'Color',colS,'LineWidth',1.25);
plot(ax2,t,LFlat_rep,'--','Color',colGray,'LineWidth',1.20);
yline(ax2,meso.ML,':','Color',[0.10 0.10 0.10],'LineWidth',0.9);
xlabel(ax2,'Time (min)');
ylabel(ax2,'Lipid-peroxidation burden, L');
legend(ax2,{'T-only','Original T+S','Flat-S T+S','M_L'}, ...
    'Box','off','Location','northwest','FontSize',5.7);
formatAxes(ax2);

% Population counterfactual.
ax3=axes(fig,'Position',[0.405 0.18 0.205 0.70]); hold(ax3,'on');
x=1:nVariant;
errorbar(ax3,x-0.08,PRIM_CF.Pcommit, ...
    PRIM_CF.Pcommit-BOOT_CF.P_CI(:,1)', ...
    BOOT_CF.P_CI(:,2)'-PRIM_CF.Pcommit, ...
    'o','Color',colCF,'MarkerFaceColor','w','MarkerSize',4.5, ...
    'LineWidth',1.0,'CapSize',3);
errorbar(ax3,x+0.08,PRIM_GF.Pcommit, ...
    PRIM_GF.Pcommit-BOOT_GF.P_CI(:,1)', ...
    BOOT_GF.P_CI(:,2)'-PRIM_GF.Pcommit, ...
    's','Color',colGF,'MarkerFaceColor','w','MarkerSize',4.5, ...
    'LineWidth',1.0,'CapSize',3);
set(ax3,'XTick',x,'XTickLabel',{'T only','Original','Flat S','Random','Aligned'}, ...
    'XTickLabelRotation',28);
xlim(ax3,[0.5 nVariant+0.5]);
ylim(ax3,[0 0.60]);
ylabel(ax3,'P(commitment by 60 min)');
title(ax3,'Persistent spatial priming dominates over precise timing', ...
    'FontWeight','normal','FontSize',8);
legend(ax3,{'CIR + Full','GammaOU + Full'}, ...
    'Box','off','Location','northwest','FontSize',5.8);
formatAxes(ax3);

% Critical priming maps.
rng(SEED_SCATTER,'twister');
ax4=axes(fig,'Position',[0.665 0.18 0.145 0.70]); hold(ax4,'on');
plotCriticalMap(ax4,LAMBDA_CF,colCF,MAX_SCATTER_POINTS,LAMBDA_MAX);
xlabel(ax4,'Temporal fate margin');
ylabel(ax4,'Critical priming, \lambda_i^*');
title(ax4,sprintf('CIR + Full  r_s = %.2f',LAMBDA_CF.spearmanMarginLambda), ...
    'FontWeight','normal','FontSize',7.5);
formatAxes(ax4);

ax5=axes(fig,'Position',[0.835 0.18 0.145 0.70]); hold(ax5,'on');
plotCriticalMap(ax5,LAMBDA_GF,colGF,MAX_SCATTER_POINTS,LAMBDA_MAX);
xlabel(ax5,'Temporal fate margin');
ylabel(ax5,'');
title(ax5,sprintf('GammaOU + Full  r_s = %.2f',LAMBDA_GF.spearmanMarginLambda), ...
    'FontWeight','normal','FontSize',7.5);
formatAxes(ax5);

mainStem=fullfile(outDir,'Fig4C_MAIN_FINAL');
saveFigureSet(fig,mainStem,SAVE_FIG,EXPORT_DPI);

%% ================================ Export ==================================
summaryTable=table( ...
    ["CF";"GF"], ...
    [PRIM_CF.Pcommit(1);PRIM_GF.Pcommit(1)], ...
    [PRIM_CF.Pcommit(2);PRIM_GF.Pcommit(2)], ...
    [PRIM_CF.Pcommit(3);PRIM_GF.Pcommit(3)], ...
    [PRIM_CF.Pcommit(4);PRIM_GF.Pcommit(4)], ...
    [PRIM_CF.Pcommit(5);PRIM_GF.Pcommit(5)], ...
    [LAMBDA_CF.medianOriginal;LAMBDA_GF.medianOriginal], ...
    [LAMBDA_CF.medianFlat;LAMBDA_GF.medianFlat], ...
    [LAMBDA_CF.corrOrigFlat;LAMBDA_GF.corrOrigFlat], ...
    [LAMBDA_CF.spearmanMarginLambda;LAMBDA_GF.spearmanMarginLambda], ...
    [LAMBDA_CF.fracAtOneOriginal;LAMBDA_GF.fracAtOneOriginal], ...
    'VariableNames',{'Condition','P_Tonly','P_Original','P_Flat', ...
    'P_Random','P_Aligned','MedianLambdaOriginal','MedianLambdaFlat', ...
    'CorrLambdaOrigFlat','SpearmanMarginLambda','FracLambdaLE1'});

variantRows=[];
for c=1:2
    if c==1
        X=PRIM_CF; B=BOOT_CF; nm='CF';
    else
        X=PRIM_GF; B=BOOT_GF; nm='GF';
    end
    for v=1:nVariant
        variantRows=[variantRows; table( ...
            string(nm),string(variantNames{v}), ...
            X.Pcommit(v),B.P_CI(v,1),B.P_CI(v,2), ...
            X.RMST(v),B.RMST_CI(v,1),B.RMST_CI(v,2), ...
            X.jointFrac(v),B.Joint_CI(v,1),B.Joint_CI(v,2), ...
            'VariableNames',{'Condition','Variant','Pcommit','P_L95','P_U95', ...
            'RMST','RMST_L95','RMST_U95','JointFrac','Joint_L95','Joint_U95'})]; %#ok<AGROW>
    end
end

lambdaRows=[buildLambdaTable(LAMBDA_CF,'CF');buildLambdaTable(LAMBDA_GF,'GF')];

if SAVE_CSV
    writetable(summaryTable,fullfile(outDir,'Fig4C_MAIN_FINAL_summary.csv'));
    writetable(variantRows,fullfile(outDir,'Fig4C_MAIN_FINAL_counterfactuals.csv'));
    writetable(lambdaRows,fullfile(outDir,'Fig4C_MAIN_FINAL_lambda_cells.csv'));
end

if SAVE_DATA
    MAIN=struct();
    MAIN.sourceFile=sourceFile;
    MAIN.dt=dt;
    MAIN.Tobs=Tobs;
    MAIN.tauL=tauL;
    MAIN.meso=meso;
    MAIN.nTD=nTD;
    MAIN.nSD=nSD;
    MAIN.nPairs=nPairs;
    MAIN.idxT=idxT;
    MAIN.idxS=idxS;
    MAIN.variantNames=variantNames;
    MAIN.randomShift=randomShift;
    MAIN.bestShiftCF=bestShiftCF;
    MAIN.bestShiftGF=bestShiftGF;
    MAIN.repIndex=repIndex;
    MAIN.lambdaMax=LAMBDA_MAX;
    MAIN.lambdaTol=LAMBDA_TOL;

    save(fullfile(outDir,'Fig4C_MAIN_FINAL_data.mat'), ...
        'MAIN','PRIM_CF','PRIM_GF','BOOT_CF','BOOT_GF', ...
        'MARGIN_CF','MARGIN_GF','LAMBDA_CF','LAMBDA_GF', ...
        'BOOT_LAM_CF','BOOT_LAM_GF', ...
        'overlapRefCF','overlapRefGF','overlapBestCF','overlapBestGF', ...
        '-v7.3');
end

fprintf('\n=================== FIG.4C MAIN FINAL SUMMARY ===================\n');
fprintf('CF Pcommit: T-only %.4f | Original %.4f | Flat %.4f | Random %.4f | Aligned %.4f\n', ...
    PRIM_CF.Pcommit);
fprintf('GF Pcommit: T-only %.4f | Original %.4f | Flat %.4f | Random %.4f | Aligned %.4f\n', ...
    PRIM_GF.Pcommit);
fprintf('CF r(lambda orig,flat)=%.4f | r_s(margin,lambda)=%.4f | P(lambda*<=1)=%.4f\n', ...
    LAMBDA_CF.corrOrigFlat,LAMBDA_CF.spearmanMarginLambda,LAMBDA_CF.fracAtOneOriginal);
fprintf('GF r(lambda orig,flat)=%.4f | r_s(margin,lambda)=%.4f | P(lambda*<=1)=%.4f\n', ...
    LAMBDA_GF.corrOrigFlat,LAMBDA_GF.spearmanMarginLambda,LAMBDA_GF.fracAtOneOriginal);
fprintf('Representative CF pair index = %d\n',repIndex);
fprintf('Outputs: %s\n',outDir);
fprintf('=================================================================\n');

%% =========================================================================
% Local functions
% =========================================================================
function path=resolveInputFile(scriptDir,name)
    candidates = {fullfile(pwd,name),fullfile(scriptDir,name), ...
        fullfile(scriptDir,'Fig4B_fate_routes_final',name)};
    path='';
    for k=1:numel(candidates)
        if isfile(candidates{k})
            path=candidates{k};
            return;
        end
    end
end

function E=positiveExcess(I)
    E=max(I-mean(I,2),0);
end

function [bestShift,score0,bestScore]=maxOverlapShift(A,B)
    [n,m]=size(A);
    if ~isequal(size(A),size(B)), error('A and B must have same size.'); end
    den=sqrt(sum(A.^2,2).*sum(B.^2,2));
    zeroDen=den<=eps;
    den=max(den,eps);
    score0=sum(A.*B,2)./den;
    score0(zeroDen)=0;
    bestScore=-Inf(n,1);
    bestShift=zeros(n,1);
    for s=0:m-1
        Bs=circshift(B,[0 s]);
        sc=sum(A.*Bs,2)./den;
        sc(zeroDen)=0;
        improve=sc>bestScore;
        bestScore(improve)=sc(improve);
        bestShift(improve)=s;
    end
end

function Y=rowCircularShift(X,shiftBins)
    [n,m]=size(X);
    s=mod(round(shiftBins(:)),m);
    if numel(s)==1, s=repmat(s,n,1); end
    if numel(s)~=n, error('shiftBins must be scalar or one value per row.'); end
    rows=repmat((1:n)',1,m);
    cols=mod(repmat(0:m-1,n,1)-s,m)+1;
    lin=sub2ind([n,m],rows,cols);
    Y=X(lin);
end

function X=initializeVariantResult(nPairs,nVariant)
    X.FPT_TS=inf(nPairs,nVariant);
    X.route=zeros(nPairs,nVariant,'uint8');
    X.Pcommit=zeros(1,nVariant);
    X.RMST=zeros(1,nVariant);
    X.jointFrac=zeros(1,nVariant);
end

function X=fillVariant(X,v,T,S,meso,dt,Tobs)
    R=evaluateCondition(T,S,meso,dt,Tobs,false);
    X.FPT_TS(:,v)=R.FPT_TS;
    X.route(:,v)=uint8(R.routeCode);
    X.Pcommit(v)=R.Pcommit;
    X.RMST(v)=R.RMST;
    X.jointFrac(v)=R.routeFrac(4);
end

function B=bootstrapVariants(X,idxT,idxS,nTD,nSD,Tobs,nBoot,alpha,seed)
    rng(seed,'twister');
    nVariant=size(X.FPT_TS,2);
    bP=nan(nBoot,nVariant);
    bRMST=nan(nBoot,nVariant);
    bJ=nan(nBoot,nVariant);

    for b=1:nBoot
        w=twoWayWeights(idxT,idxS,nTD,nSD);
        if sum(w)<=0, continue; end
        for v=1:nVariant
            f=X.FPT_TS(:,v);
            rc=X.route(:,v);
            bP(b,v)=weightedMean(f<=Tobs,w);
            bRMST(b,v)=weightedMean(min(f,Tobs),w);
            bJ(b,v)=weightedMean(rc==4,w);
        end
        if mod(b,max(1,round(nBoot/10)))==0
            fprintf('  bootstrap %d / %d\n',b,nBoot);
        end
    end

    B.P_CI=ci2D(bP,alpha);
    B.RMST_CI=ci2D(bRMST,alpha);
    B.Joint_CI=ci2D(bJ,alpha);
end

function M=temporalMarginOnly(T,meso,dt,Tobs)
    [FT,LT]=propagateLPOExact(meso.omegaT*T, ...
        meso.rhoL,meso.ML,meso.L0,dt,Tobs,true);
    maxLT=max(LT,[],2);
    M=struct();
    M.FPT_T=FT;
    M.maxLT=maxLT;
    M.margin=max(meso.ML-maxLT,0);
    M.Tinsufficient=FT>Tobs;
end

function X=criticalLambdaAnalysis(T,S,SFlat,MARGIN,meso,dt,Tobs,lambdaMax,tol,nIter)
    n=size(T,1); %#ok<NASGU>

    FS=propagateLPOExact(meso.omegaS*S,meso.rhoL,meso.ML,meso.L0,dt,Tobs,false);
    FSflat=propagateLPOExact(meso.omegaS*SFlat,meso.rhoL,meso.ML,meso.L0,dt,Tobs,false);

    Tinsufficient=MARGIN.Tinsufficient;
    [lamOrig,resOrig]=solveCriticalLambda(T,S,meso,dt,Tobs,lambdaMax,tol,nIter);
    [lamFlat,resFlat]=solveCriticalLambda(T,SFlat,meso,dt,Tobs,lambdaMax,tol,nIter);

    Rorig=evaluateCondition(T,S,meso,dt,Tobs,false);
    Rflat=evaluateCondition(T,SFlat,meso,dt,Tobs,false);

    predOrig=Tinsufficient & isfinite(lamOrig) & lamOrig<=1+2*tol;
    predFlat=Tinsufficient & isfinite(lamFlat) & lamFlat<=1+2*tol;
    actualOrig=Tinsufficient & Rorig.FPT_TS<=Tobs;
    actualFlat=Tinsufficient & Rflat.FPT_TS<=Tobs;

    finiteBoth=Tinsufficient & resOrig & resFlat & isfinite(lamOrig) & isfinite(lamFlat);
    finiteMargin=Tinsufficient & resOrig & isfinite(lamOrig) & isfinite(MARGIN.margin);

    X=struct();
    X.lambdaOriginal=lamOrig;
    X.lambdaFlat=lamFlat;
    X.resolvedOriginal=resOrig;
    X.resolvedFlat=resFlat;
    X.Tinsufficient=Tinsufficient;
    X.SinsufficientOriginal=FS>Tobs;
    X.SinsufficientFlat=FSflat>Tobs;
    X.routeOriginal=Rorig.routeCode;
    X.routeFlat=Rflat.routeCode;
    X.temporalMargin=MARGIN.margin;

    X.auditOriginal=classificationAccuracy(predOrig,actualOrig,Tinsufficient);
    X.auditFlat=classificationAccuracy(predFlat,actualFlat,Tinsufficient);

    X.medianOriginal=median(lamOrig(Tinsufficient & isfinite(lamOrig)),'omitnan');
    X.medianFlat=median(lamFlat(Tinsufficient & isfinite(lamFlat)),'omitnan');
    X.medianAbsDifference=median(abs(lamFlat(finiteBoth)-lamOrig(finiteBoth)),'omitnan');

    if sum(finiteBoth)>=3
        C=corrcoef(lamOrig(finiteBoth),lamFlat(finiteBoth));
        X.corrOrigFlat=C(1,2);
    else
        X.corrOrigFlat=NaN;
    end

    X.spearmanMarginLambda=spearmanLocal(MARGIN.margin(finiteMargin),lamOrig(finiteMargin));
    X.fracAtOneOriginal=mean(lamOrig(Tinsufficient)<=1+2*tol);
    X.fracAtOneFlat=mean(lamFlat(Tinsufficient)<=1+2*tol);
end

function [lambdaStar,resolved]=solveCriticalLambda(T,S,meso,dt,Tobs,lambdaMax,tol,nIter)
    n=size(T,1);
    lambdaStar=inf(n,1);
    resolved=false(n,1);

    F0=propagateLPOExact(meso.omegaT*T, ...
        meso.rhoL,meso.ML,meso.L0,dt,Tobs,false);
    cross0=F0<=Tobs;
    lambdaStar(cross0)=0;
    resolved(cross0)=true;

    active=find(~cross0);
    if isempty(active), return; end

    hi0=lambdaMax*ones(numel(active),1);
    crossHi=crossAtScaledSpatial(T(active,:),S(active,:),hi0,meso,dt,Tobs);

    solvable=active(crossHi);
    unresolved=active(~crossHi);
    resolved(solvable)=true;
    lambdaStar(unresolved)=Inf;

    if isempty(solvable), return; end

    lo=zeros(numel(solvable),1);
    hi=lambdaMax*ones(numel(solvable),1);
    Ts=T(solvable,:);
    Ss=S(solvable,:);

    for it=1:nIter
        mid=0.5*(lo+hi);
        crossMid=crossAtScaledSpatial(Ts,Ss,mid,meso,dt,Tobs);
        hi(crossMid)=mid(crossMid);
        lo(~crossMid)=mid(~crossMid);
        if max(hi-lo)<=tol, break; end
    end

    lambdaStar(solvable)=hi;
end

function cross=crossAtScaledSpatial(T,S,lambda,meso,dt,Tobs)
    lambda=lambda(:);
    if numel(lambda)==1
        lambda=repmat(lambda,size(T,1),1);
    end
    U=meso.omegaT*T + meso.omegaS*(S.*lambda);
    F=propagateLPOExact(U,meso.rhoL,meso.ML,meso.L0,dt,Tobs,false);
    cross=F<=Tobs;
end

function B=bootstrapLambdaSummary(X,idxT,idxS,nTD,nSD,nBoot,alpha,seed)
    rng(seed,'twister');
    bMedOrig=nan(nBoot,1);
    bMedFlat=nan(nBoot,1);
    bFracOrig=nan(nBoot,1);
    bFracFlat=nan(nBoot,1);

    validOrig=X.Tinsufficient & isfinite(X.lambdaOriginal);
    validFlat=X.Tinsufficient & isfinite(X.lambdaFlat);

    for b=1:nBoot
        w=twoWayWeights(idxT,idxS,nTD,nSD);
        bMedOrig(b)=weightedQuantile(X.lambdaOriginal(validOrig),w(validOrig),0.5);
        bMedFlat(b)=weightedQuantile(X.lambdaFlat(validFlat),w(validFlat),0.5);
        ii=X.Tinsufficient;
        bFracOrig(b)=weightedMean(X.lambdaOriginal(ii)<=1,w(ii));
        bFracFlat(b)=weightedMean(X.lambdaFlat(ii)<=1,w(ii));
    end

    B.medianOriginal_CI=quantileLocal(bMedOrig,[alpha/2 1-alpha/2]);
    B.medianFlat_CI=quantileLocal(bMedFlat,[alpha/2 1-alpha/2]);
    B.fracAtOneOriginal_CI=quantileLocal(bFracOrig,[alpha/2 1-alpha/2]);
    B.fracAtOneFlat_CI=quantileLocal(bFracFlat,[alpha/2 1-alpha/2]);
end

function R=evaluateCondition(T,S,meso,dt,Tobs,returnHist)
    if nargin<6, returnHist=false; end
    [FT,LT]=propagateLPOExact(meso.omegaT*T, ...
        meso.rhoL,meso.ML,meso.L0,dt,Tobs,returnHist);
    [FS,LS]=propagateLPOExact(meso.omegaS*S, ...
        meso.rhoL,meso.ML,meso.L0,dt,Tobs,returnHist);
    [FTS,LTS]=propagateLPOExact(meso.omegaT*T+meso.omegaS*S, ...
        meso.rhoL,meso.ML,meso.L0,dt,Tobs,returnHist);

    cT=FT<=Tobs; cS=FS<=Tobs; cTS=FTS<=Tobs;
    if any((cT|cS)&~cTS)
        error('Comparison-principle violation in T+S.');
    end

    rc=zeros(size(FTS),'uint8');
    rc(cT & ~cS)=1;
    rc(~cT & cS)=2;
    rc(cT & cS)=3;
    rc(~cT & ~cS & cTS)=4;
    rc(~cTS)=5;

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
    L=L0*ones(n,1);
    FPT=inf(n,1);
    if returnHistory, Lhist=zeros(n,nBins); else, Lhist=[]; end

    for b=1:nBins
        Ub=max(U(:,b),0);
        a=Ub+rho;
        Lnext=L;
        ok=a>eps;
        if any(ok)
            av=a(ok); uv=Ub(ok);
            Lss=uv./av;
            Lnext(ok)=Lss+(L(ok)-Lss).*exp(-av*dt);
        end
        Lnext=min(max(Lnext,0),1);

        cross=isinf(FPT)&(L<ML)&(Lnext>=ML);
        if any(cross)
            ii=find(cross);
            av=a(ii); uv=Ub(ii);
            Lss=uv./av;
            ratio=(ML-Lss)./(L(ii)-Lss);
            ratio=max(min(ratio,1),realmin);
            tc=-log(ratio)./av;
            tc=max(min(tc,dt),0);
            FPT(ii)=(b-1)*dt+tc;
        end

        L=Lnext;
        if returnHistory, Lhist(:,b)=L; end
    end
    FPT(FPT>Tobs)=inf;
end

function plotCriticalMap(ax,X,cMain,maxPoints,lambdaMax)
    valid=X.Tinsufficient & X.resolvedOriginal & ...
        isfinite(X.lambdaOriginal) & isfinite(X.temporalMargin);
    ii=find(valid);
    if numel(ii)>maxPoints
        ii=ii(randperm(numel(ii),maxPoints));
    end

    low=ii(X.lambdaOriginal(ii)<=1);
    high=ii(X.lambdaOriginal(ii)>1);

    scatter(ax,X.temporalMargin(high),X.lambdaOriginal(high),7, ...
        [0.72 0.72 0.72],'filled','MarkerFaceAlpha',0.23);
    scatter(ax,X.temporalMargin(low),X.lambdaOriginal(low),8, ...
        [0.93 0.66 0.20],'filled','MarkerFaceAlpha',0.48);

    yline(ax,1,':','Color',[0.15 0.15 0.15],'LineWidth',0.9);

    q=equalCountBins(X.temporalMargin,valid,4,'ascend');
    xm=nan(1,4); ym=nan(1,4);
    for k=1:4
        jj=q==k;
        xm(k)=median(X.temporalMargin(jj),'omitnan');
        ym(k)=median(X.lambdaOriginal(jj),'omitnan');
    end
    plot(ax,xm,ym,'o-','Color',cMain,'MarkerFaceColor','w', ...
        'LineWidth',1.25,'MarkerSize',4.3);

    xlim(ax,[0 max(0.50,max(X.temporalMargin(valid))*1.02)]);
    ymax=min(lambdaMax,max(2.25,quantileLocal(X.lambdaOriginal(valid),0.99)*1.05));
    ylim(ax,[0 ymax]);
end

function bins=equalCountBins(score,eligible,nBin,order)
    bins=zeros(size(score));
    idx=find(eligible & isfinite(score));
    if isempty(idx), return; end
    [~,ord]=sort(score(idx),order);
    idx=idx(ord);
    N=numel(idx);
    edges=round(linspace(0,N,nBin+1));
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
    x=double(x);
    ok=isfinite(x)&isfinite(w)&w>0;
    if ~any(ok), m=NaN; else, m=sum(w(ok).*x(ok))/sum(w(ok)); end
end

function q=weightedQuantile(x,w,p)
    x=x(:); w=w(:);
    ok=isfinite(x)&isfinite(w)&w>0;
    x=x(ok); w=w(ok);
    if isempty(x), q=NaN; return; end
    [x,ord]=sort(x); w=w(ord);
    cw=cumsum(w);
    k=find(cw>=p*cw(end),1,'first');
    q=x(k);
end

function CI=ci2D(X,alpha)
    nCol=size(X,2);
    CI=nan(nCol,2);
    for j=1:nCol
        x=X(:,j); x=x(isfinite(x));
        if ~isempty(x), CI(j,:)=quantileLocal(x,[alpha/2 1-alpha/2]); end
    end
end

function q=quantileLocal(x,p)
    x=x(isfinite(x)); x=sort(x(:));
    n=numel(x); p=p(:).'; q=nan(size(p));
    if n==0, return; end
    for j=1:numel(p)
        pp=min(max(p(j),0),1);
        pos=1+pp*(n-1);
        lo=floor(pos); hi=ceil(pos);
        if lo==hi
            q(j)=x(lo);
        else
            q(j)=x(lo)+(pos-lo)*(x(hi)-x(lo));
        end
    end
end

function r=spearmanLocal(x,y)
    x=x(:); y=y(:);
    ok=isfinite(x)&isfinite(y);
    x=x(ok); y=y(ok);
    if numel(x)<3, r=NaN; return; end
    rx=rankAverageTies(x); ry=rankAverageTies(y);
    C=corrcoef(rx,ry); r=C(1,2);
end

function r=rankAverageTies(x)
    [xs,ord]=sort(x);
    n=numel(x); rs=zeros(n,1);
    i=1;
    while i<=n
        j=i;
        while j<n && xs(j+1)==xs(i), j=j+1; end
        rs(i:j)=mean(i:j);
        i=j+1;
    end
    r=zeros(n,1); r(ord)=rs;
end

function a=classificationAccuracy(pred,actual,eligible)
    ii=find(eligible);
    if isempty(ii), a=NaN; else, a=mean(pred(ii)==actual(ii)); end
end

function T=buildLambdaTable(X,condition)
    n=numel(X.lambdaOriginal);
    T=table(repmat(string(condition),n,1),(1:n)', ...
        X.Tinsufficient,X.SinsufficientOriginal, ...
        X.lambdaOriginal,X.lambdaFlat,X.resolvedOriginal,X.resolvedFlat, ...
        X.temporalMargin,X.routeOriginal,X.routeFlat, ...
        'VariableNames',{'Condition','PairIndex','Tinsufficient', ...
        'SinsufficientOriginal','LambdaStarOriginal','LambdaStarFlat', ...
        'ResolvedOriginal','ResolvedFlat','TemporalMargin', ...
        'RouteOriginal','RouteFlat'});
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
