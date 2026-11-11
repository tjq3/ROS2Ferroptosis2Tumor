function SI_M5_Q1_StatisticalValidation_FINAL()
% Supplementary statistical / numerical support for Fig.5B Q1 inference.
% Plotting only; no new SSA.

clc; close all;

scriptDir=fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir=pwd; end

dataFile=fullfile(scriptDir,'Macro_Q1_TargetedRefinement_FINAL_OUT', ...
    'Macro_Q1_TargetedRefinement_FINAL.mat');
assert(isfile(dataFile),'Missing targeted-Q1 MAT file.');

S=load(dataFile,'Q1R');
Q=S.Q1R;

outDir=fullfile(scriptDir,'SI_M5_Q1_StatisticalValidation_FINAL_OUT');
if ~exist(outDir,'dir'), mkdir(outDir); end

outPNG=fullfile(outDir,'SI_M5_Q1_StatisticalValidation_FINAL.png');
outPDF=fullfile(outDir,'SI_M5_Q1_StatisticalValidation_FINAL.pdf');
outFIG=fullfile(outDir,'SI_M5_Q1_StatisticalValidation_FINAL.fig');
outMAT=fullfile(outDir,'SI_M5_Q1_StatisticalValidation_FINAL_data.mat');

chiF=double(Q.chiF(:)).';
chiMu=double(Q.chiMu(:)).';

Praw=double(Q.Pext);
Plo=double(Q.PextLo);
Phi=double(Q.PextHi);
Piso=double(Q.PextIsotonic);

pLevels=double(Q.boundary.pLevels(:)).';
B=double(Q.boundary.pExt);
BLo=double(Q.boundary.pExtLo);
BHi=double(Q.boundary.pExtHi);

gArrest=double(Q.boundary.growthArrest(:));
gArrestLo=double(Q.boundary.growthArrestLo(:));
gArrestHi=double(Q.boundary.growthArrestHi(:));

isBenchmarkMu=logical(Q.isOriginalMuBenchmark(:));

iMuRef=double(Q.referenceMu.index);
chiMuRef=double(Q.referenceMu.chiMu);

nRep=reshape(double(Q.tables.metrics.nRep),numel(chiF),numel(chiMu)).';
mcseNV=reshape(double(Q.tables.metrics.relMCSE_AUTC_NV),numel(chiF),numel(chiMu)).';
lateStart=double(Q.lateStartSurvivalFraction);

isoCorrection=Piso-Praw;

iP50=find(abs(pLevels-0.50)<1e-12,1);
iP90=find(abs(pLevels-0.90)<1e-12,1);

seqStops={'#eaf7e5','#d5eece','#bbe3be','#90d4bd', ...
    '#64c0cf','#3e9fca','#1a78b8','#085293'};
cmapSeq=interpHexPalette(seqStops,256);

divStops={'#2c5aa0','#8ca9d0','#f4f4f2','#eaa091','#b23a34'};
cmapDiv=interpHexPalette(divStops,257);

cGrowth=hex2rgb('#6D5A9B');
cP50=hex2rgb('#1A78B8');
cP90=hex2rgb('#C85A4A');

oldFont=get(groot,'defaultAxesFontName');
oldText=get(groot,'defaultTextFontName');
set(groot,'defaultAxesFontName','Arial');
set(groot,'defaultTextFontName','Arial');
set(groot,'defaultAxesFontSize',9);
set(groot,'defaultTextFontSize',9);
set(groot,'defaultAxesLineWidth',0.8);
set(groot,'defaultAxesTickDir','out');
cleanupObj=onCleanup(@()restoreDefaults(oldFont,oldText)); %#ok<NASGU>

fig=figure('Color','w','Position',[70 40 1500 850]);
tl=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');

% a
ax=nexttile(tl,1); hold(ax,'on');
errLo=Praw(iMuRef,:)-Plo(iMuRef,:);
errHi=Phi(iMuRef,:)-Praw(iMuRef,:);

errorbar(ax,chiF,Praw(iMuRef,:),errLo,errHi,'o', ...
    'Color',[0.35 0.35 0.35],'MarkerFaceColor','w', ...
    'CapSize',2.5,'LineWidth',0.8,'DisplayName','raw + Wilson 95% CI');

plot(ax,chiF,Piso(iMuRef,:),'-','Color',cP50,'LineWidth',1.8, ...
    'DisplayName','weighted isotonic fit');

yline(ax,0.10,':','Color',[0.65 0.65 0.65],'HandleVisibility','off');
yline(ax,0.50,':','Color',[0.65 0.65 0.65],'HandleVisibility','off');
yline(ax,0.90,':','Color',[0.65 0.65 0.65],'HandleVisibility','off');

xlabel(ax,'\chi_F = \lambda_F/r');
ylabel(ax,'P(V extinction by 168 h)');
title(ax,sprintf('a  Raw vs isotonic at \\chi_\\mu=%.3f',chiMuRef), ...
    'FontWeight','normal','HorizontalAlignment','left');
legend(ax,'Location','southeast','Box','off','FontSize',7.5);
ylim(ax,[-0.03 1.03]);
styleAxis(ax);

% b
ax=nexttile(tl,2);
imagesc(ax,chiF,chiMu,isoCorrection);
set(ax,'YDir','normal');
colormap(ax,cmapDiv);

lim=max(abs(isoCorrection(:)));
if lim==0, lim=0.01; end
caxis(ax,[-lim lim]);

xlabel(ax,'\chi_F = \lambda_F/r');
ylabel(ax,'\chi_\mu = \mu/r');
title(ax,'b  Raw-to-isotonic correction','FontWeight','normal', ...
    'HorizontalAlignment','left');
cb=colorbar(ax); cb.Label.String='P_{iso} - P_{raw}';
styleAxis(ax);

% c
ax=nexttile(tl,3); hold(ax,'on');
plotRibbon(ax,chiMu,gArrest,gArrestLo,gArrestHi,cGrowth,':','growth arrest');
plotRibbon(ax,chiMu,B(:,iP50),BLo(:,iP50),BHi(:,iP50),cP50,'-','P_{ext}=0.50');
plotRibbon(ax,chiMu,B(:,iP90),BLo(:,iP90),BHi(:,iP90),cP90,'--','P_{ext}=0.90');

scatter(ax,chiMu(isBenchmarkMu),B(isBenchmarkMu,iP50),28,cP50,'filled', ...
    'MarkerEdgeColor','w','LineWidth',0.4,'HandleVisibility','off');
scatter(ax,chiMu(~isBenchmarkMu),B(~isBenchmarkMu,iP50),28, ...
    'MarkerFaceColor','w','MarkerEdgeColor',cP50, ...
    'LineWidth',0.8,'HandleVisibility','off');

xlabel(ax,'\chi_\mu = \mu/r');
ylabel(ax,'Boundary \chi_F');
title(ax,'c  Bootstrap boundary precision','FontWeight','normal', ...
    'HorizontalAlignment','left');
legend(ax,'Location','northwest','Box','off','FontSize',7.5);
styleAxis(ax);

% d
ax=nexttile(tl,4);
imagesc(ax,chiF,chiMu,nRep);
set(ax,'YDir','normal');
colormap(ax,cmapSeq);

xlabel(ax,'\chi_F = \lambda_F/r');
ylabel(ax,'\chi_\mu = \mu/r');
title(ax,'d  Adaptive Monte-Carlo allocation','FontWeight','normal', ...
    'HorizontalAlignment','left');
cb=colorbar(ax); cb.Label.String='replicates per state';
styleAxis(ax);

% e
ax=nexttile(tl,5);
imagesc(ax,chiF,chiMu,100*mcseNV);
set(ax,'YDir','normal');
colormap(ax,cmapSeq);

xlabel(ax,'\chi_F = \lambda_F/r');
ylabel(ax,'\chi_\mu = \mu/r');
title(ax,'e  Precision of integrated viable burden','FontWeight','normal', ...
    'HorizontalAlignment','left');
cb=colorbar(ax); cb.Label.String='relative MCSE of AUTC_{NV} (%)';
styleAxis(ax);

% f
ax=nexttile(tl,6);
contourf(ax,chiF,chiMu,lateStart,0:0.1:1,'LineStyle','none');
colormap(ax,cmapSeq); hold(ax,'on');
plot(ax,B(:,iP90),chiMu,'--','Color',cP90,'LineWidth',1.2);

xlabel(ax,'\chi_F = \lambda_F/r');
ylabel(ax,'\chi_\mu = \mu/r');
title(ax,'f  Viable survival into the late window','FontWeight','normal', ...
    'HorizontalAlignment','left');
cb=colorbar(ax); cb.Label.String='P(N_V(120 h)>0)';
styleAxis(ax);

sgtitle(tl,'Statistical and numerical support for the refined Q1 boundaries', ...
    'FontWeight','normal','FontSize',11);

try
    axs=findall(fig,'Type','axes');
    for k=1:numel(axs), axtoolbar(axs(k),{}); end
catch
end

exportgraphics(fig,outPNG,'Resolution',600);
exportgraphics(fig,outPDF,'ContentType','vector');
savefig(fig,outFIG);

SAVE=struct();
SAVE.chiF=chiF; SAVE.chiMu=chiMu;
SAVE.Praw=Praw; SAVE.Piso=Piso;
SAVE.isotonicCorrection=isoCorrection;
SAVE.nRep=nRep; SAVE.mcseNV=mcseNV;
SAVE.lateStartSurvival=lateStart;
save(outMAT,'SAVE');

fprintf('\nSI-M5 statistical validation generated; no new simulation.\n');
fprintf('Saved:\n%s\n%s\n%s\n%s\n',outPNG,outPDF,outFIG,outMAT);

end

function styleAxis(ax)
ax.FontName='Arial'; ax.FontSize=9; ax.LineWidth=0.8;
ax.TickDir='out'; ax.TickLength=[0.015 0.015];
ax.Box='off'; ax.Layer='top'; grid(ax,'off');
end

function plotRibbon(ax,x,y,lo,hi,c,lineStyle,labelText)
x=x(:); y=y(:); lo=lo(:); hi=hi(:);
mask=isfinite(x)&isfinite(y)&isfinite(lo)&isfinite(hi);

if nnz(mask)>=2
    xx=x(mask); yy=y(mask); ll=lo(mask); hh=hi(mask);
    patch(ax,[xx;flipud(xx)],[ll;flipud(hh)],c, ...
        'FaceAlpha',0.12,'EdgeColor','none','HandleVisibility','off');
    plot(ax,xx,yy,'LineStyle',lineStyle,'Color',c, ...
        'LineWidth',1.6,'DisplayName',labelText);
end
end

function cmap=interpHexPalette(hexStops,n)
rgb=zeros(numel(hexStops),3);
for i=1:numel(hexStops), rgb(i,:)=hex2rgb(hexStops{i}); end
x=linspace(0,1,size(rgb,1)); xi=linspace(0,1,n);
cmap=zeros(n,3);
for j=1:3, cmap(:,j)=interp1(x,rgb(:,j),xi,'pchip'); end
cmap=max(0,min(1,cmap));
end

function rgb=hex2rgb(h)
h=char(h); if h(1)=='#', h=h(2:end); end
rgb=[hex2dec(h(1:2)) hex2dec(h(3:4)) hex2dec(h(5:6))]/255;
end

function restoreDefaults(oldFont,oldText)
set(groot,'defaultAxesFontName',oldFont);
set(groot,'defaultTextFontName',oldText);
end
