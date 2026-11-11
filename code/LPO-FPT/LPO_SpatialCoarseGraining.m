%% Fig4_PreC_SpatialCoarseGraining_Audit.m
% =========================================================================
% Pre-Fig.4C diagnostic:
% Does the weak Full-vs-Uniform effect at the mesoscopic level arise from
% the chosen spatial coarse-graining parameters?
%
% This diagnostic varies TWO interface choices while leaving the underlying
% SD-RDME physics unchanged:
%
%   (1) Local-window radius:
%         Rw = multiplier * Rw_baseline
%
%   (2) Whole-cell pooling exponent p:
%         I_S^(p) = [ mean_k( j_k^p ) ]^(1/p)
%
%       p=1   : spatial mean
%       p=2   : current RMS primary interface
%       p=4,8 : progressively hotspot-sensitive
%       p=Inf : maximum local-window exposure
%
% For EACH (Rw,p), an independent Uniform-RDME calibration gives B_S(Rw,p),
% which is then FIXED for both Uniform and Full. This prevents an operator
% from winning merely because it has a different numerical scale.
%
% omega_T, omega_S, M_L and tau_L are kept at the FINAL Fig.4B values.
% Thus this is an interface audit, not a post-hoc re-tuning of biology.
%
% Outputs:
%   - F/U ratios for pooled mean, pooled Q95 and integrated spatial dose
%   - mesoscopic Full-vs-Uniform delta RMST in CIR/GammaOU backgrounds
%   - Full-vs-Uniform change in Joint-required fraction
%
% Interpretation:
%   If smaller Rw / larger p strongly restores U-vs-F separation, the weak
%   main-text spatial effect is primarily caused by whole-cell averaging.
%   If separation remains weak across the grid, attenuation is intrinsic to
%   the current cross-scale mapping rather than one arbitrary RMS setting.
%
% NOTE: This diagnostic reruns SD-RDME. Start with RUN_MODE='pilot'.
% =========================================================================
clear; close all; clc;
clear functions;

RUN_MODE='pilot';   % 'smoke' | 'pilot' | 'final'

FINAL_MESO_NAME='Fig4B_fateRoutes_data.mat';
SD_DATA_NAME='Fig3BC_predictorIndependent_data.mat';

radiusMultipliers=[0.50 0.75 1.00 1.50 2.00];
pList=[1 2 4 8 Inf];

USE_PARALLEL=true;
AUTO_START_POOL=true;
N_WORKERS=8;

SAVE_FIG=true;
SAVE_CSV=true;
EXPORT_DPI=600;

switch lower(RUN_MODE)
    case 'smoke'
        nCal=8;  nSD=20;  nPairs=1000;
    case 'pilot'
        nCal=30; nSD=80;  nPairs=5000;
    case 'final'
        nCal=60; nSD=200; nPairs=12000;
    otherwise
        error('Unknown RUN_MODE.');
end

seedCal=20260901;
seedEval=20260902;
seedPairT=20260903;
seedPairS=20260904;

%% ------------------------------ Load files -------------------------------
scriptDir=fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir=pwd; end
mesoFile=resolveInputFile(scriptDir,FINAL_MESO_NAME);
sdFile=resolveInputFile(scriptDir,SD_DATA_NAME);
if isempty(mesoFile), error('Cannot find %s',FINAL_MESO_NAME); end
if isempty(sdFile), error('Cannot find %s',SD_DATA_NAME); end

M=load(mesoFile);
SD=load(sdFile);

IT_C=M.inputs.IT_C;
IT_G=M.inputs.IT_G;
meso=M.meta.meso;
dt=M.meta.dtBin_min;
Tobs=M.meta.Tobs_min;
tEdges=0:dt:Tobs;

% Locked geometry and fields from the current Fig.3 model.
xCoord=SD.x(:).'; yCoord=SD.y(:).';
mask=logical(SD.mask);
[Xg,Yg]=meshgrid(xCoord,yCoord);
K=nnz(mask);
Psd=SD.meta.parameters;

D0_Fe2=Psd.D0_Fe2;
if isfield(Psd,'D0_Fe3'), D0_Fe3=Psd.D0_Fe3; else, D0_Fe3=0.70*D0_Fe2; end
c2_0=Psd.c2_0;
FeTot=Psd.FeTot_per_voxel;
fracFe2=Psd.fracFe2_0;
alphaMatch=SD.meta.alphaMatch;
R0=SD.meta.windowRadius;

G=SD.G;
c1FullMap=SD.c1_sd;
D1FullMap=SD.D_Fe2_sd;
D2FullMap=D0_Fe3.*G; D2FullMap(~mask)=NaN;

c1F=c1FullMap(mask);
c2F=c2_0*ones(K,1);
D1F=D1FullMap(mask);
D2F=D2FullMap(mask);

D1Umap=NaN(size(D1FullMap));
D2Umap=NaN(size(D2FullMap));
D1Umap(mask)=mean(D1F);
D2Umap(mask)=mean(D2F);
c1U=alphaMatch*mean(c1F)*ones(K,1);
c2U=c2_0*ones(K,1);

idxMap=buildMaskIndexMap(mask);
dx=mean(diff(xCoord)); dy=mean(diff(yCoord));
[fromF1,toF1,qF1]=buildDirectedEdgesHarmonic(mask,idxMap,D1FullMap,dx,dy);
[fromF2,toF2,qF2]=buildDirectedEdgesHarmonic(mask,idxMap,D2FullMap,dx,dy);
[fromU1,toU1,qU1]=buildDirectedEdgesHarmonic(mask,idxMap,D1Umap,dx,dy);
[fromU2,toU2,qU2]=buildDirectedEdgesHarmonic(mask,idxMap,D2Umap,dx,dy);
outF=buildOutgoingLists(K,fromF1,toF1,qF1,fromF2,toF2,qF2);
outU=buildOutgoingLists(K,fromU1,toU1,qU1,fromU2,toU2,qU2);

N0=FeTot*ones(K,1);
x1_0=round(fracFe2*N0);
x2_0=N0-x1_0;

coordX=Xg(mask); coordY=Yg(mask);
nR=numel(radiusMultipliers);
wins=cell(nR,1);
for ir=1:nR
    wins{ir}=buildCircularWindowsSparse(coordX,coordY,R0*radiusMultipliers(ir));
end

usePar=configureParallel(USE_PARALLEL,AUTO_START_POOL,N_WORKERS);

outDir=fullfile(scriptDir,['Fig4_PreC_spatial_coarsegrain_' lower(RUN_MODE)]);
if ~exist(outDir,'dir'), mkdir(outDir); end

%% ------------------------ Simulate calibration/evaluation ----------------
fprintf('\nUniform calibration bank...\n');
CAL=simulateSpatialMultiBank(c1U,c2U,outU,wins,pList,x1_0,x2_0,tEdges,nCal,seedCal,usePar);

fprintf('\nUniform evaluation bank...\n');
Uraw=simulateSpatialMultiBank(c1U,c2U,outU,wins,pList,x1_0,x2_0,tEdges,nSD,seedEval,usePar);

fprintf('\nFull evaluation bank...\n');
Fraw=simulateSpatialMultiBank(c1F,c2F,outF,wins,pList,x1_0,x2_0,tEdges,nSD,seedEval,usePar);

% Shape: nRep x nBins x nR x nP
nP=numel(pList);
BS=zeros(nR,nP);
for ir=1:nR
    for ip=1:nP
        tmp=CAL(:,:,ir,ip);
        BS(ir,ip)=mean(tmp(:));
    end
end

%% ------------------------------ Pair TD/SD -------------------------------
idxT=balancedSourceIndices(size(IT_C,1),nPairs,seedPairT);
idxS=balancedSourceIndices(nSD,nPairs,seedPairS);
ITCp=IT_C(idxT,:); ITGp=IT_G(idxT,:);

meanFU=zeros(nR,nP);
q95FU=zeros(nR,nP);
doseFU=zeros(nR,nP);
pairedQ95=zeros(nR,nP);
dRMST_C=zeros(nR,nP);
dRMST_G=zeros(nR,nP);
dJoint_C=zeros(nR,nP);
dJoint_G=zeros(nR,nP);

for ir=1:nR
    for ip=1:nP
        IU=squeeze(Uraw(:,:,ir,ip))/BS(ir,ip);
        IF=squeeze(Fraw(:,:,ir,ip))/BS(ir,ip);

        meanFU(ir,ip)=mean(IF(:))/max(mean(IU(:)),eps);
        q95FU(ir,ip)=quantileLocal(IF(:),0.95)/max(quantileLocal(IU(:),0.95),eps);
        doseU=sum(IU,2)*dt;
        doseF=sum(IF,2)*dt;
        doseFU(ir,ip)=mean(doseF)/max(mean(doseU),eps);

        qU=rowQuantileLocal(IU,0.95);
        qF=rowQuantileLocal(IF,0.95);
        pairedQ95(ir,ip)=mean(qF>qU);

        IUp=IU(idxS,:);
        IFp=IF(idxS,:);
        R=evaluateFourBackgrounds(ITCp,ITGp,IUp,IFp,meso,dt,Tobs);

        % Full minus Uniform. Negative delta RMST means Full commits earlier.
        dRMST_C(ir,ip)=R(3).RMST60-R(1).RMST60;
        dRMST_G(ir,ip)=R(4).RMST60-R(2).RMST60;
        dJoint_C(ir,ip)=R(3).routeFrac(4)-R(1).routeFrac(4);
        dJoint_G(ir,ip)=R(4).routeFrac(4)-R(2).routeFrac(4);
    end
end

%% -------------------------------- Figures --------------------------------
set(groot,'defaultFigureColor','w');
set(groot,'defaultAxesFontName','Arial');

xlabels=cell(1,nP);
for ip=1:nP
    if isinf(pList(ip)), xlabels{ip}='Inf'; else, xlabels{ip}=num2str(pList(ip)); end
end
ylabels=arrayfun(@(x)sprintf('%.2g x R_0',x),radiusMultipliers,'UniformOutput',false);

fig=figure('Units','centimeters','Position',[1 1 18.3 12.0],'Color','w');
tl=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');

plotHeat(nexttile,meanFU,xlabels,ylabels,'Mean spatial input, F/U');
plotHeat(nexttile,q95FU,xlabels,ylabels,'Q_{0.95}(I_S), F/U');
plotHeat(nexttile,doseFU,xlabels,ylabels,'Integrated dose, F/U');
plotHeat(nexttile,dRMST_C,xlabels,ylabels,'\DeltaRMST Full-Uniform | CIR (min)');
plotHeat(nexttile,dRMST_G,xlabels,ylabels,'\DeltaRMST Full-Uniform | GammaOU (min)');
plotHeat(nexttile,0.5*(dJoint_C+dJoint_G),xlabels,ylabels,'Mean \Delta joint-required fraction');

if SAVE_FIG
    exportgraphics(fig,fullfile(outDir,'SpatialCoarseGraining_Audit.png'),'Resolution',EXPORT_DPI);
    exportgraphics(fig,fullfile(outDir,'SpatialCoarseGraining_Audit.pdf'),'ContentType','vector');
end

%% -------------------------------- Export ---------------------------------
rows=[];
for ir=1:nR
    for ip=1:nP
        rows=[rows; table(radiusMultipliers(ir),pList(ip),BS(ir,ip), ...
            meanFU(ir,ip),q95FU(ir,ip),doseFU(ir,ip),pairedQ95(ir,ip), ...
            dRMST_C(ir,ip),dRMST_G(ir,ip),dJoint_C(ir,ip),dJoint_G(ir,ip), ...
            'VariableNames',{'RadiusMultiplier','PoolingP','B_S', ...
            'MeanFU','Q95FU','DoseFU','PairedFractionQ95FgtU', ...
            'DeltaRMST_C','DeltaRMST_G','DeltaJoint_C','DeltaJoint_G'})]; %#ok<AGROW>
    end
end
if SAVE_CSV
    writetable(rows,fullfile(outDir,'SpatialCoarseGraining_Audit.csv'));
end

fprintf('\nBaseline cell (R=R0, p=2):\n');
[~,ir0]=min(abs(radiusMultipliers-1));
[~,ip0]=min(abs(pList-2));
fprintf('  Mean F/U       = %.4f\n',meanFU(ir0,ip0));
fprintf('  Q95 F/U        = %.4f\n',q95FU(ir0,ip0));
fprintf('  Dose F/U       = %.4f\n',doseFU(ir0,ip0));
fprintf('  DeltaRMST CIR  = %.4f min\n',dRMST_C(ir0,ip0));
fprintf('  DeltaRMST GOU  = %.4f min\n',dRMST_G(ir0,ip0));
fprintf('Outputs saved in %s\n',outDir);

%% ============================= Local functions ============================
function path=resolveInputFile(scriptDir,name)
    cands = {fullfile(pwd,name),fullfile(scriptDir,name), ...
        fullfile(scriptDir,'Fig4B_fate_routes_final',name), ...
        fullfile(fileparts(scriptDir),'SD-RDME',name)};
    path = '';
    for k=1:numel(cands), if isfile(cands{k}),path=cands{k};return;end,end
end

function usePar=configureParallel(requested,autoStart,nWorkers)
    usePar=false;
    if ~requested, return; end
    try
        paths={pwd,fileparts(mfilename('fullpath'))};
        if any(cellfun(@hasNonASCII,paths))
            warning('Non-ASCII project path: using serial loop for this diagnostic.');
            return;
        end
        if ~license('test','Distrib_Computing_Toolbox'), return; end
        p=gcp('nocreate');
        if isempty(p)&&autoStart, p=parpool('local',nWorkers); end
        usePar=~isempty(p);
    catch
        usePar=false;
    end
end

function tf=hasNonASCII(s)
    tf=any(double(char(s))>127);
end

function A=simulateSpatialMultiBank(c1,c2,out,wins,pList,x10,x20,tEdges,nRep,seed,usePar)
    nBins=numel(tEdges)-1; nR=numel(wins); nP=numel(pList);
    A=zeros(nRep,nBins,nR,nP);
    if usePar
        parfor r=1:nRep
            rng(seed+r,'twister');
            A(r,:,:,:)=simulateSpatialMultiOne(c1,c2,out,wins,pList,x10,x20,tEdges);
        end
    else
        for r=1:nRep
            rng(seed+r,'twister');
            A(r,:,:,:)=simulateSpatialMultiOne(c1,c2,out,wins,pList,x10,x20,tEdges);
            if mod(r,max(1,round(nRep/10)))==0
                fprintf('  finished %d / %d trajectories\n',r,nRep);
            end
        end
    end
end

function A=simulateSpatialMultiOne(c1,c2,out,wins,pList,x10,x20,tEdges)
    K=numel(c1); nBins=numel(tEdges)-1; dt=tEdges(2)-tEdges(1);
    nR=numel(wins); nP=numel(pList);
    x1=x10(:); x2=x20(:);

    counts=zeros(K,nR);
    A=zeros(1,nBins,nR,nP);

    rate=x1.*(c1+out.q1sum)+x2.*(c2+out.q2sum);
    tree=fenwickBuild(rate);

    t=tEdges(1); ib=1; Tend=tEdges(end);
    while t<Tend
        a0=fenwickTotal(tree);
        if ~(a0>0)&&isfinite(a0), break; end
        if ~(a0>0)||~isfinite(a0), break; end
        tau=-log(rand)/a0; tNext=t+tau;

        while ib<=nBins && tEdges(ib+1)<=min(tNext,Tend)+1e-14
            for ir=1:nR
                dens=counts(:,ir)./(wins{ir}.size*dt);
                for ip=1:nP
                    pp=pList(ip);
                    if isinf(pp), val=max(dens);
                    else, val=(mean(dens.^pp))^(1/pp); end
                    A(1,ib,ir,ip)=val;
                end
            end
            counts(:)=0; ib=ib+1;
        end
        if tNext>Tend, break; end

        k=fenwickFind(tree,rand*a0);
        aR1=c1(k)*x1(k); aR2=c2(k)*x2(k);
        aD1=out.q1sum(k)*x1(k); aD2=out.q2sum(k)*x2(k);
        aLocal=aR1+aR2+aD1+aD2;
        if aLocal<=0
            old=rate(k); rate(k)=0; tree=fenwickAdd(tree,k,-old); t=tNext; continue;
        end

        u=rand*aLocal; touched=k;
        if u<=aR1
            x1(k)=x1(k)-1; x2(k)=x2(k)+1;
            for ir=1:nR
                centers=wins{ir}.centersByVoxel{k};
                if ~isempty(centers), counts(centers,ir)=counts(centers,ir)+1; end
            end
        elseif u<=aR1+aR2
            x2(k)=x2(k)-1; x1(k)=x1(k)+1;
        elseif u<=aR1+aR2+aD1
            j=chooseDestination(out.to1{k},out.q1{k});
            x1(k)=x1(k)-1; x1(j)=x1(j)+1; touched=[k;j];
        else
            j=chooseDestination(out.to2{k},out.q2{k});
            x2(k)=x2(k)-1; x2(j)=x2(j)+1; touched=[k;j];
        end

        [rate,tree]=refreshVoxelRates(unique(touched),rate,tree,x1,x2,c1,c2,out.q1sum,out.q2sum);
        t=tNext;
    end

    while ib<=nBins
        for ir=1:nR
            dens=counts(:,ir)./(wins{ir}.size*dt);
            for ip=1:nP
                pp=pList(ip);
                if isinf(pp), val=max(dens);
                else, val=(mean(dens.^pp))^(1/pp); end
                A(1,ib,ir,ip)=val;
            end
        end
        counts(:)=0; ib=ib+1;
    end
end

function R=evaluateFourBackgrounds(IT_C,IT_G,IS_U,IS_F,meso,dt,Tobs)
    UTC=meso.omegaT*IT_C; UTG=meso.omegaT*IT_G;
    USU=meso.omegaS*IS_U; USF=meso.omegaS*IS_F;
    FTC=propagateLPOExact(UTC,meso.rhoL,meso.ML,meso.L0,dt,Tobs);
    FTG=propagateLPOExact(UTG,meso.rhoL,meso.ML,meso.L0,dt,Tobs);
    FSU=propagateLPOExact(USU,meso.rhoL,meso.ML,meso.L0,dt,Tobs);
    FSF=propagateLPOExact(USF,meso.rhoL,meso.ML,meso.L0,dt,Tobs);
    R(1)=assemble(FTC,FSU,propagateLPOExact(UTC+USU,meso.rhoL,meso.ML,meso.L0,dt,Tobs),Tobs);
    R(2)=assemble(FTG,FSU,propagateLPOExact(UTG+USU,meso.rhoL,meso.ML,meso.L0,dt,Tobs),Tobs);
    R(3)=assemble(FTC,FSF,propagateLPOExact(UTC+USF,meso.rhoL,meso.ML,meso.L0,dt,Tobs),Tobs);
    R(4)=assemble(FTG,FSF,propagateLPOExact(UTG+USF,meso.rhoL,meso.ML,meso.L0,dt,Tobs),Tobs);
end

function R=assemble(FT,FS,FTS,Tobs)
    cT=FT<=Tobs;cS=FS<=Tobs;cTS=FTS<=Tobs;
    rc=zeros(size(FTS));
    rc(cT&~cS)=1;rc(~cT&cS)=2;rc(cT&cS)=3;
    rc(~cT&~cS&cTS)=4;rc(~cTS)=5;
    R.FPT_T=FT;R.FPT_S=FS;R.FPT_TS=FTS;
    R.routeFrac=arrayfun(@(z)mean(rc==z),1:5);
    R.Pcross60=mean(cTS);R.RMST60=mean(min(FTS,Tobs));
end

function FPT=propagateLPOExact(U,rho,ML,L0,dt,Tobs)
    n=size(U,1);L=L0*ones(n,1);FPT=inf(n,1);
    for b=1:size(U,2)
        Ub=max(U(:,b),0);a=Ub+rho;Ln=L;ok=a>eps;
        if any(ok)
            av=a(ok);uv=Ub(ok);Lss=uv./av;
            Ln(ok)=Lss+(L(ok)-Lss).*exp(-av*dt);
        end
        Ln=min(max(Ln,0),1);
        cr=isinf(FPT)&L<ML&Ln>=ML;
        if any(cr)
            ii=find(cr);av=a(ii);uv=Ub(ii);Lss=uv./av;
            rr=(ML-Lss)./(L(ii)-Lss);rr=max(min(rr,1),realmin);
            tc=-log(rr)./av;tc=max(min(tc,dt),0);
            FPT(ii)=(b-1)*dt+tc;
        end
        L=Ln;
    end
    FPT(FPT>Tobs)=inf;
end

function idx=balancedSourceIndices(nBank,nPairs,seed)
    rng(seed,'twister');idx=zeros(nPairs,1);pos=1;
    while pos<=nPairs
        p=randperm(nBank);take=min(nBank,nPairs-pos+1);
        idx(pos:pos+take-1)=p(1:take);pos=pos+take;
    end
end

function idxMap=buildMaskIndexMap(mask)
    idxMap=zeros(size(mask));idxMap(mask)=1:nnz(mask);
end

function [from,to,q]=buildDirectedEdgesHarmonic(mask,idxMap,Dfield,dx,dy)
    [nr,nc]=size(mask);from=[];to=[];q=[];dirs=[0 1;0 -1;1 0;-1 0];
    for r=1:nr
        for c=1:nc
            if ~mask(r,c),continue;end
            k=idxMap(r,c);
            for d=1:4
                rr=r+dirs(d,1);cc=c+dirs(d,2);
                if rr<1||rr>nr||cc<1||cc>nc||~mask(rr,cc),continue;end
                l=idxMap(rr,cc);Dk=Dfield(r,c);Dl=Dfield(rr,cc);
                Dij=2*Dk*Dl/max(Dk+Dl,eps);
                if dirs(d,1)==0,h2=dx^2;else,h2=dy^2;end
                from(end+1,1)=k;to(end+1,1)=l;q(end+1,1)=Dij/h2;
            end
        end
    end
end

function out=buildOutgoingLists(K,from1,to1,q1,from2,to2,q2)
    out.to1=cell(K,1);out.q1=cell(K,1);out.to2=cell(K,1);out.q2=cell(K,1);
    for e=1:numel(q1),k=from1(e);out.to1{k}(end+1,1)=to1(e);out.q1{k}(end+1,1)=q1(e);end
    for e=1:numel(q2),k=from2(e);out.to2{k}(end+1,1)=to2(e);out.q2{k}(end+1,1)=q2(e);end
    out.q1sum=zeros(K,1);out.q2sum=zeros(K,1);
    for k=1:K,out.q1sum(k)=sum(out.q1{k});out.q2sum(k)=sum(out.q2{k});end
end

function win=buildCircularWindowsSparse(x,y,radius)
    K=numel(x);I=cell(K,1);J=cell(K,1);
    for k=1:K
        idx=find((x-x(k)).^2+(y-y(k)).^2<=radius^2+1e-14);
        I{k}=repmat(k,numel(idx),1);J{k}=idx(:);
    end
    I=vertcat(I{:});J=vertcat(J{:});W=sparse(I,J,1,K,K);
    win.size=full(sum(W,2));win.centersByVoxel=cell(K,1);
    for j=1:K,win.centersByVoxel{j}=find(W(:,j));end
end

function [rate,tree]=refreshVoxelRates(touched,rate,tree,x1,x2,c1,c2,q1s,q2s)
    for ii=1:numel(touched)
        k=touched(ii);new=max(x1(k)*(c1(k)+q1s(k))+x2(k)*(c2(k)+q2s(k)),0);
        tree=fenwickAdd(tree,k,new-rate(k));rate(k)=new;
    end
end

function j=chooseDestination(to,q)
    cs=cumsum(q);j=to(find(cs>=rand*cs(end),1,'first'));
end

function tree=fenwickBuild(v)
    tree=zeros(numel(v),1);for i=1:numel(v),tree=fenwickAdd(tree,i,v(i));end
end

function tree=fenwickAdd(tree,idx,delta)
    n=numel(tree);
    while idx<=n
        tree(idx)=tree(idx)+delta;lb=idx-bitand(idx,idx-1);idx=idx+lb;
    end
end

function s=fenwickTotal(tree)
    idx=numel(tree);s=0;
    while idx>0,s=s+tree(idx);lb=idx-bitand(idx,idx-1);idx=idx-lb;end
end

function idx=fenwickFind(tree,target)
    n=numel(tree);if target<=0,target=realmin;end;idx=0;bit=2^floor(log2(n));
    while bit>0
        nxt=idx+bit;
        if nxt<=n&&tree(nxt)<target,idx=nxt;target=target-tree(nxt);end
        bit=floor(bit/2);
    end
    idx=min(idx+1,n);
end

function q=rowQuantileLocal(X,p)
    X=sort(X,2);n=size(X,2);pos=1+p*(n-1);lo=floor(pos);hi=ceil(pos);
    if lo==hi,q=X(:,lo);else,q=X(:,lo)+(pos-lo)*(X(:,hi)-X(:,lo));end
end

function q=quantileLocal(x,p)
    x=sort(x(:));n=numel(x);p=p(:).';q=nan(size(p));
    for j=1:numel(p)
        pos=1+min(max(p(j),0),1)*(n-1);lo=floor(pos);hi=ceil(pos);
        if lo==hi,q(j)=x(lo);else,q(j)=x(lo)+(pos-lo)*(x(hi)-x(lo));end
    end
end

function plotHeat(ax,M,xlabs,ylabs,ttl)
    imagesc(ax,M);axis(ax,'tight');set(ax,'YDir','normal');
    set(ax,'XTick',1:numel(xlabs),'XTickLabel',xlabs, ...
        'YTick',1:numel(ylabs),'YTickLabel',ylabs);
    xlabel(ax,'Pooling exponent, p');ylabel(ax,'Local-window radius');
    title(ax,ttl,'FontWeight','normal');colorbar(ax);
    set(ax,'FontName','Arial','FontSize',7,'TickDir','out','Box','on');
end
