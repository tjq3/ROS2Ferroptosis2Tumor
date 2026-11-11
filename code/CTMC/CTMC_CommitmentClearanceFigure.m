function CTMC_CommitmentClearanceFigure()
% =========================================================================
% CTMC_CommitmentClearanceFigure
%
% Revised FINAL main-text Q1 plotting script.
%
% Layout:
%   B-i   Extinction landscape (full height)
%   B-ii  Representative growth phenotypes (upper ~2/3 of middle column)
%         + compact phenotype-translation strip (lower ~1/3)
%   B-iii Clearance modulation of transition boundaries (full height)
%
% The compact strip restores the scientifically useful content from the old
% phenotype-translation right panel, but in a non-redundant form:
%   - four mechanistic sectors
%   - five representative phenotype states
%   - three formal transition boundaries
% No state-specific g/g0 or Pext text is repeated in the strip; this avoids
% the inconsistent/overcrowded labeling problem. Those metrics are saved to
% CSV for the figure legend / SI table.
%
% Reads Macro_Q1_TargetedRefinement_FINAL.mat only.
% No new SSA simulation.
% =========================================================================

clc; close all;

%% 0. Locate FINAL Q1 file
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir = pwd; end

candidateFiles = { ...
    fullfile(scriptDir,'Macro_Q1_TargetedRefinement_FINAL.mat'), ...
    fullfile(scriptDir,'Macro_Q1_TargetedRefinement_FINAL_OUT', ...
        'Macro_Q1_TargetedRefinement_FINAL.mat')};

dataFile = '';
for k = 1:numel(candidateFiles)
    if isfile(candidateFiles{k})
        dataFile = candidateFiles{k};
        break;
    end
end
assert(~isempty(dataFile),'Cannot find Macro_Q1_TargetedRefinement_FINAL.mat.');

S = load(dataFile,'Q1R','RES');
assert(isfield(S,'Q1R') && isfield(S,'RES'), ...
    'FINAL MAT must contain Q1R and RES.');
Q = S.Q1R;
RES = S.RES;

%% 1. Read and verify FINAL quantities
chiF  = double(Q.chiF(:)).';
chiMu = double(Q.chiMu(:)).';
nChi = numel(chiF);
nMu  = numel(chiMu);

PextIso = orientStateMatrix(double(Q.PextIsotonic),nMu,nChi,'PextIsotonic');
PextRaw = orientStateMatrix(double(Q.Pext),nMu,nChi,'Pext');
gLate   = orientStateMatrix(double(Q.survivorAwareLateGrowth),nMu,nChi,'gLate');
gRatio  = orientStateMatrix(double(Q.survivorLateGrowthRatio),nMu,nChi,'gRatio');

pLevels = double(Q.boundary.pLevels(:)).';
B   = orientBoundaryMatrix(double(Q.boundary.pExt),nMu,numel(pLevels),'pExt');
BLo = orientBoundaryMatrix(double(Q.boundary.pExtLo),nMu,numel(pLevels),'pExtLo');
BHi = orientBoundaryMatrix(double(Q.boundary.pExtHi),nMu,numel(pLevels),'pExtHi');

gArrest   = double(Q.boundary.growthArrest(:));
gArrestLo = double(Q.boundary.growthArrestLo(:));
gArrestHi = double(Q.boundary.growthArrestHi(:));
isBenchmarkMu = logical(Q.isOriginalMuBenchmark(:));

iP50 = find(abs(pLevels-0.50)<1e-12,1);
iP90 = find(abs(pLevels-0.90)<1e-12,1);
assert(~isempty(iP50) && ~isempty(iP90), ...
    'Pext=0.50 and/or Pext=0.90 boundary is missing.');

iMuRef = round(double(Q.referenceMu.index));
chiMuRef = double(Q.referenceMu.chiMu);
assert(iMuRef>=1 && iMuRef<=nMu,'Invalid referenceMu.index.');

chiArrest = gArrest(iMuRef);
chiP50 = B(iMuRef,iP50);
chiP90 = B(iMuRef,iP90);

assert(isfinite(chiArrest) && isfinite(chiP50) && isfinite(chiP90), ...
    'Reference-mu boundaries must be finite.');
assert(chiArrest < chiP50 && chiP50 < chiP90, ...
    'Expected hierarchy arrest < P50 < P90 was not found.');

realChiF = double(Q.realMesoscopic.chiF);
realChiMin = min(realChiF(:));

%% 2. Select three main-text trajectory exemplars
% Growth-dominant
iGrowth = nearestIndex(chiF,0);

% Indolent-like: positive late growth, low extinction, below arrest boundary
targetInd = (2/3)*chiArrest;
eligibleInd = PextRaw(iMuRef,:) < 0.10 & ...
              isfinite(gLate(iMuRef,:)) & ...
              gLate(iMuRef,:) > 0 & ...
              chiF < chiArrest;

iInd = nearestEligibleIndex(chiF,targetInd,eligibleInd);
if isempty(iInd)
    iInd = nearestIndex(chiF,targetInd);
end

% Near-arrest / bounded-like
iBound = nearestIndex(chiF,chiArrest);

mainIdx = [iGrowth iInd iBound];
mainChiF = chiF(mainIdx);
mainPext = PextRaw(iMuRef,mainIdx);
mainGRatio = gRatio(iMuRef,mainIdx);

mainLabels = { ...
    'Growth-dominant (aggressive-like)', ...
    'Indolent-like persistence', ...
    'Near-arrest / bounded-like'};

assert(numel(unique(mainIdx))==3, ...
    'Three main phenotype representatives are not unique.');

%% 3. Select five states for the compact translation strip
% 1 growth-dominant
iAgg = iGrowth;

% 2 indolent-like
iInd5 = iInd;

% 3 near-arrest
iNear = iBound;

% 4 regression-like: midpoint between arrest and P50
targetReg = 0.5*(chiArrest + chiP50);
eligibleReg = chiF > chiArrest & chiF < chiP50;
iReg = nearestEligibleIndex(chiF,targetReg,eligibleReg);
if isempty(iReg)
    iReg = nearestIndex(chiF,targetReg);
end

% 5 extinction-dominated: nearest P90 boundary
iExt = nearestIndex(chiF,chiP90);

stripIdx = [iAgg iInd5 iNear iReg iExt];
stripChiF = chiF(stripIdx);
stripPext = PextRaw(iMuRef,stripIdx);
stripGRatio = gRatio(iMuRef,stripIdx);

assert(numel(unique(stripIdx))==5, ...
    'Five phenotype-translation representatives are not unique.');

stripNames = { ...
    'Growth-dominant', ...
    'Indolent-like', ...
    'Near-arrest', ...
    'Regression', ...
    'Extinction-dominated'};

%% 4. Read model-derived trajectories
assert(iscell(RES),'RES must be a cell array.');

if isequal(size(RES),[nMu nChi])
    RESm = RES;
elseif isequal(size(RES),[nChi nMu])
    RESm = RES.';
else
    error('Unexpected RES size [%d %d].',size(RES,1),size(RES,2));
end

t_h = double(RESm{iMuRef,mainIdx(1)}.t_h(:));
NV = nan(numel(t_h),3);

for k=1:3
    Rk = RESm{iMuRef,mainIdx(k)};
    assert(isfield(Rk,'t_h') && isfield(Rk,'meanNVcurve'), ...
        'Representative RES cell lacks t_h or meanNVcurve.');
    assert(isequal(double(Rk.t_h(:)),t_h), ...
        'Representative trajectory time grids differ.');
    NV(:,k) = double(Rk.meanNVcurve(:));
end

NVnorm = NV ./ NV(1,:);

%% 5. Output directory
outDir = fullfile(scriptDir,'Fig5B_Q1_Main_REVISED_WITH_TRANSLATION_FINAL_OUT');
if ~exist(outDir,'dir'), mkdir(outDir); end

outPNG = fullfile(outDir,'Fig5B_Q1_Main_REVISED_WITH_TRANSLATION_FINAL.png');
outPDF = fullfile(outDir,'Fig5B_Q1_Main_REVISED_WITH_TRANSLATION_FINAL.pdf');
outFIG = fullfile(outDir,'Fig5B_Q1_Main_REVISED_WITH_TRANSLATION_FINAL.fig');
outBoundCSV = fullfile(outDir,'Fig5B_Q1_Main_REVISED_Boundaries.csv');
outPhenCSV = fullfile(outDir,'Fig5B_Q1_Main_REVISED_Phenotypes.csv');
outMAT = fullfile(outDir,'Fig5B_Q1_Main_REVISED_WITH_TRANSLATION_FINAL_data.mat');

%% 6. Journal-style colors
seqStops = {'#eaf7e5','#d5eece','#bbe3be','#90d4bd', ...
    '#64c0cf','#3e9fca','#1a78b8','#085293'};
cmapSeq = interpHexPalette(seqStops,256);

cGrowth = hex2rgb('#6D5A9B');
cP50    = hex2rgb('#1A78B8');
cP90    = hex2rgb('#C85A4A');

% Main trajectory colors
cTraj = [ ...
    hex2rgb('#D24B40'); ...
    hex2rgb('#4E9A96'); ...
    hex2rgb('#8B8B8B')];

% Five phenotype state colors
cState = [ ...
    hex2rgb('#D24B40'); ...
    hex2rgb('#4E9A96'); ...
    hex2rgb('#8B8B8B'); ...
    hex2rgb('#7562A8'); ...
    hex2rgb('#9E3D32')];

% Compact sector fills
cSectorGrowth = hex2rgb('#EAF7E5');
cSectorReg    = hex2rgb('#DDEBE9');
cSectorStoch  = hex2rgb('#F4E8DD');
cSectorExt    = hex2rgb('#F2DDDA');

%% 7. Global graphics
oldFont = get(groot,'defaultAxesFontName');
oldText = get(groot,'defaultTextFontName');

set(groot,'defaultAxesFontName','Arial');
set(groot,'defaultTextFontName','Arial');
set(groot,'defaultAxesFontSize',9);
set(groot,'defaultTextFontSize',9);
set(groot,'defaultAxesLineWidth',0.8);
set(groot,'defaultAxesTickDir','out');

cleanupObj = onCleanup(@()restoreDefaults(oldFont,oldText)); %#ok<NASGU>

%% 8. Figure layout
% 3 rows x 3 columns:
%   left  : B-i spans all 3 rows
%   middle: B-ii trajectory spans top 2 rows
%           phenotype strip occupies bottom 1 row (< 1/2 height)
%   right : B-iii spans all 3 rows

fig = figure('Color','w','Position',[60 90 1640 610]);
tl = tiledlayout(fig,3,3,'TileSpacing','compact','Padding','compact');

%% B-i: extinction landscape
ax1 = nexttile(tl,1,[3 1]);

contourf(ax1,chiF,chiMu,PextIso,0:0.1:1,'LineStyle','none');
colormap(ax1,cmapSeq);
caxis(ax1,[0 1]);
hold(ax1,'on');

plot(ax1,gArrest,chiMu,':','LineWidth',1.7,'Color',cGrowth);
plot(ax1,B(:,iP50),chiMu,'-','LineWidth',1.8,'Color',cP50);
plot(ax1,B(:,iP90),chiMu,'--','LineWidth',1.6,'Color',cP90);

xlabel(ax1,'Relative ferroptotic commitment, \chi_F = \lambda_F/r');
ylabel(ax1,'Relative clearance, \chi_\mu = \mu/r');
title(ax1,'b-i  Ferroptotic commitment defines the extinction landscape', ...
    'FontWeight','normal','HorizontalAlignment','left');

cb = colorbar(ax1);
cb.Label.String = 'P(V extinction by 168 h)';

xlim(ax1,[chiF(1) chiF(end)]);
ylim(ax1,[chiMu(1) chiMu(end)]);

text(ax1,0.07,0.95,'growth-arrest boundary', ...
    'Units','normalized','Color',cGrowth,'FontSize',7.8);
text(ax1,0.07,0.885,'P_{ext}(168 h)=0.50', ...
    'Units','normalized','Color',cP50,'FontSize',7.8);
text(ax1,0.07,0.82,'P_{ext}(168 h)=0.90', ...
    'Units','normalized','Color',cP90,'FontSize',7.8);

text(ax1,0.98,0.055, ...
    sprintf('pro-ferroptotic mesoscopic states: \\chi_F \\geq %.2f (off-scale)',realChiMin), ...
    'Units','normalized','HorizontalAlignment','right','VerticalAlignment','bottom', ...
    'FontSize',7.2,'Color',[0.28 0.28 0.28]);

styleAxis(ax1);

%% B-ii upper: three representative normalized trajectories
ax2 = nexttile(tl,2,[2 1]);
hold(ax2,'on');

for k=1:3
    plot(ax2,t_h,NVnorm(:,k), ...
        'LineWidth',2.1,'Color',cTraj(k,:));
end

yline(ax2,1,':','Color',[0.72 0.72 0.72], ...
    'LineWidth',0.8,'HandleVisibility','off');

xlabel(ax2,'Time (h)');
ylabel(ax2,'Normalized viable burden, N_V(t)/N_V(0)');
title(ax2,'b-ii  Representative growth phenotypes under ferroptotic pressure', ...
    'FontWeight','normal','HorizontalAlignment','left');

xlim(ax2,[0 max(t_h)]);
ylim(ax2,[0 max(NVnorm,[],'all')*1.08]);

xLab = 0.70*max(t_h);
for k=1:3
    [~,ix] = min(abs(t_h-xLab));
    yk = NVnorm(ix,k);

    if k==1
        dy = +0.10;
    elseif k==2
        dy = +0.03;
    else
        dy = -0.07;
    end

    text(ax2,xLab,yk+dy, ...
        sprintf('%s\n\\chi_F=%.1f',mainLabels{k},mainChiF(k)), ...
        'Color',cTraj(k,:), ...
        'FontSize',7.7, ...
        'BackgroundColor','w','Margin',1);
end

text(ax2,0.03,0.04, ...
    sprintf('reference clearance: \\chi_\\mu = %.3f',chiMuRef), ...
    'Units','normalized','VerticalAlignment','bottom', ...
    'FontSize',7.3,'Color',[0.35 0.35 0.35]);

styleAxis(ax2);

%% B-ii lower: compact phenotype-translation strip
axStrip = nexttile(tl,8);
hold(axStrip,'on');

xMin = 0;
xMax = max([chiF(end),chiP90+0.25]);

% Four mechanistic sectors
patch(axStrip,[xMin chiArrest chiArrest xMin], ...
    [0 0 1 1],cSectorGrowth,'EdgeColor','none','HandleVisibility','off');
patch(axStrip,[chiArrest chiP50 chiP50 chiArrest], ...
    [0 0 1 1],cSectorReg,'EdgeColor','none','HandleVisibility','off');
patch(axStrip,[chiP50 chiP90 chiP90 chiP50], ...
    [0 0 1 1],cSectorStoch,'EdgeColor','none','HandleVisibility','off');
patch(axStrip,[chiP90 xMax xMax chiP90], ...
    [0 0 1 1],cSectorExt,'EdgeColor','none','HandleVisibility','off');

% Three formal boundaries
xline(axStrip,chiArrest,':','Color',cGrowth,'LineWidth',1.4, ...
    'HandleVisibility','off');
xline(axStrip,chiP50,'-','Color',cP50,'LineWidth',1.5, ...
    'HandleVisibility','off');
xline(axStrip,chiP90,'--','Color',cP90,'LineWidth',1.4, ...
    'HandleVisibility','off');

% Compact sector headings
text(axStrip,(xMin+chiArrest)/2,0.92,'growth-positive', ...
    'HorizontalAlignment','center','VerticalAlignment','top', ...
    'FontSize',6.8,'Color',[0.28 0.28 0.28]);

text(axStrip,(chiArrest+chiP50)/2,0.92,'regression', ...
    'HorizontalAlignment','center','VerticalAlignment','top', ...
    'FontSize',6.8,'Color',[0.28 0.28 0.28]);

text(axStrip,(chiP50+chiP90)/2,0.92,'stochastic extinction', ...
    'HorizontalAlignment','center','VerticalAlignment','top', ...
    'FontSize',6.5,'Color',[0.28 0.28 0.28]);

text(axStrip,(chiP90+xMax)/2,0.92,'extinction-dominated', ...
    'HorizontalAlignment','center','VerticalAlignment','top', ...
    'FontSize',6.5,'Color',[0.28 0.28 0.28]);

% Five representative states
yDot = 0.57;
scatter(axStrip,stripChiF,yDot*ones(size(stripChiF)),48,cState,'filled', ...
    'MarkerEdgeColor','w','LineWidth',0.6);

% Labels: only phenotype name + chi_F, consistently for ALL five states
% (g/g0 and Pext are intentionally moved to CSV / caption to avoid clutter.)
nameY = [0.35 0.28 0.35 0.28 0.35];
for k=1:5
    text(axStrip,stripChiF(k),nameY(k), ...
        sprintf('%s\n\\chi_F=%.1f',stripNames{k},stripChiF(k)), ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','top', ...
        'FontSize',6.4, ...
        'Color',cState(k,:), ...
        'FontWeight','bold');
end

% Full boundary definitions, with values explicitly tied to the reference chi_mu
text(axStrip,0.01,0.04, ...
    sprintf(['at \\chi_\\mu=%.3f:  ', ...
             '\\chi_F^{arrest}=%.2f (g_{late}=0);  ', ...
             '\\chi_F^{50}=%.2f [P_{ext}(168 h)=0.50];  ', ...
             '\\chi_F^{90}=%.2f [P_{ext}(168 h)=0.90]'], ...
             chiMuRef,chiArrest,chiP50,chiP90), ...
    'Units','normalized','VerticalAlignment','bottom', ...
    'FontSize',6.2,'Color',[0.25 0.25 0.25]);

xlabel(axStrip,'Relative ferroptotic commitment, \chi_F = \lambda_F/r');
ylabel(axStrip,'');
title(axStrip,'Phenotype translation along the Q1 control axis', ...
    'FontWeight','normal','FontSize',8);

xlim(axStrip,[xMin xMax]);
ylim(axStrip,[0 1]);
yticks(axStrip,[]);
styleAxis(axStrip);
axStrip.FontSize = 7.5;

%% B-iii: clearance modulation of formal boundaries
ax3 = nexttile(tl,3,[3 1]);
hold(ax3,'on');

plotRibbon(ax3,chiMu,gArrest,gArrestLo,gArrestHi,cGrowth,':', ...
    'Growth-arrest boundary');
plotRibbon(ax3,chiMu,B(:,iP50),BLo(:,iP50),BHi(:,iP50),cP50,'-', ...
    'P_{ext}(168 h)=0.50');
plotRibbon(ax3,chiMu,B(:,iP90),BLo(:,iP90),BHi(:,iP90),cP90,'--', ...
    'P_{ext}(168 h)=0.90');

scatter(ax3,chiMu(isBenchmarkMu),gArrest(isBenchmarkMu),26,cGrowth,'filled', ...
    'MarkerEdgeColor','w','LineWidth',0.4,'HandleVisibility','off');
scatter(ax3,chiMu(isBenchmarkMu),B(isBenchmarkMu,iP50),26,cP50,'filled', ...
    'MarkerEdgeColor','w','LineWidth',0.4,'HandleVisibility','off');
scatter(ax3,chiMu(isBenchmarkMu),B(isBenchmarkMu,iP90),26,cP90,'filled', ...
    'MarkerEdgeColor','w','LineWidth',0.4,'HandleVisibility','off');

xlabel(ax3,'Relative clearance, \chi_\mu = \mu/r');
ylabel(ax3,'Boundary in \chi_F');
title(ax3,'b-iii  Clearance modulates extinction thresholds', ...
    'FontWeight','normal','HorizontalAlignment','left');

legend(ax3,'Location','northwest','Box','off','FontSize',7.7);
xlim(ax3,[min(chiMu) max(chiMu)]);
ylim(ax3,[0 max(B(:,iP90))*1.13]);
styleAxis(ax3);

%% 9. Export
try
    axs=findall(fig,'Type','axes');
    for k=1:numel(axs)
        axtoolbar(axs(k),{});
    end
catch
end

exportgraphics(fig,outPNG,'Resolution',600);
exportgraphics(fig,outPDF,'ContentType','vector');
savefig(fig,outFIG);

%% 10. Save tables / data
Tbound = table( ...
    chiMu(:),gArrest,gArrestLo,gArrestHi, ...
    B(:,iP50),BLo(:,iP50),BHi(:,iP50), ...
    B(:,iP90),BLo(:,iP90),BHi(:,iP90),isBenchmarkMu, ...
    'VariableNames',{ ...
    'chiMu','chiF_growthArrest','chiF_growthArrest_lo','chiF_growthArrest_hi', ...
    'chiF_Pext50','chiF_Pext50_lo','chiF_Pext50_hi', ...
    'chiF_Pext90','chiF_Pext90_lo','chiF_Pext90_hi','isOriginalMuBenchmark'});
writetable(Tbound,outBoundCSV);

Phenotype = string(stripNames(:));
chiF_rep = stripChiF(:);
lateGrowth_vs_noF = stripGRatio(:);
PextV_168h = stripPext(:);

Tphen = table(Phenotype,chiF_rep,lateGrowth_vs_noF,PextV_168h, ...
    'VariableNames',{ ...
    'Phenotype','chiF_lambda_over_r','lateGrowth_vs_noF','PextV_168h'});
writetable(Tphen,outPhenCSV);

OUT = struct();
OUT.sourceFile = dataFile;
OUT.referenceChiMu = chiMuRef;

OUT.boundary.growthArrest = chiArrest;
OUT.boundary.Pext50 = chiP50;
OUT.boundary.Pext90 = chiP90;

OUT.mainPhenotype.labels = mainLabels;
OUT.mainPhenotype.chiF = mainChiF;
OUT.mainPhenotype.Pext = mainPext;
OUT.mainPhenotype.lateGrowthRatio = mainGRatio;
OUT.mainPhenotype.t_h = t_h;
OUT.mainPhenotype.NVnorm = NVnorm;

OUT.translation.labels = stripNames;
OUT.translation.chiF = stripChiF;
OUT.translation.Pext = stripPext;
OUT.translation.lateGrowthRatio = stripGRatio;

save(outMAT,'OUT');

%% 11. Console summary
fprintf('\n============================================================\n');
fprintf('REVISED Q1 MAIN FIGURE WITH COMPACT PHENOTYPE TRANSLATION\n');
fprintf('Reference chi_mu = %.6f\n',chiMuRef);
fprintf('Boundaries at reference chi_mu:\n');
fprintf('  chi_F^arrest = %.4f\n',chiArrest);
fprintf('  chi_F^50     = %.4f\n',chiP50);
fprintf('  chi_F^90     = %.4f\n',chiP90);

fprintf('\nMain B-ii trajectory exemplars:\n');
for k=1:3
    fprintf('  %-34s chi_F=%.2f | g/g0=%+.3f | Pext=%.3f\n', ...
        mainLabels{k},mainChiF(k),mainGRatio(k),mainPext(k));
end

fprintf('\nCompact translation strip states:\n');
for k=1:5
    fprintf('  %-22s chi_F=%.2f | g/g0=%+.3f | Pext=%.3f\n', ...
        stripNames{k},stripChiF(k),stripGRatio(k),stripPext(k));
end

fprintf('\nNo new SSA simulation was performed.\n');
fprintf('Saved:\n');
fprintf('  %s\n',outPNG);
fprintf('  %s\n',outPDF);
fprintf('  %s\n',outFIG);
fprintf('  %s\n',outBoundCSV);
fprintf('  %s\n',outPhenCSV);
fprintf('  %s\n',outMAT);
fprintf('============================================================\n');

end

%% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================
function M = orientStateMatrix(M,nMu,nChi,label)
if isequal(size(M),[nMu nChi])
    return;
elseif isequal(size(M),[nChi nMu])
    M = M.';
else
    error('%s has unexpected size [%d %d].',label,size(M,1),size(M,2));
end
end

function B = orientBoundaryMatrix(B,nMu,nLevel,label)
if isequal(size(B),[nMu nLevel])
    return;
elseif isequal(size(B),[nLevel nMu])
    B = B.';
else
    error('%s has unexpected size [%d %d].',label,size(B,1),size(B,2));
end
end

function idx = nearestIndex(x,target)
[~,idx] = min(abs(x-target));
end

function idx = nearestEligibleIndex(x,target,eligible)
cand = find(logical(eligible(:)).');
if isempty(cand)
    idx = [];
    return;
end
[~,k] = min(abs(x(cand)-target));
idx = cand(k);
end

function styleAxis(ax)
ax.FontName='Arial';
ax.FontSize=9;
ax.LineWidth=0.8;
ax.TickDir='out';
ax.TickLength=[0.015 0.015];
ax.Box='off';
ax.Layer='top';
grid(ax,'off');
end

function plotRibbon(ax,x,y,lo,hi,c,lineStyle,labelText)
x=x(:); y=y(:); lo=lo(:); hi=hi(:);
mask=isfinite(x)&isfinite(y)&isfinite(lo)&isfinite(hi);
if nnz(mask)>=2
    xx=x(mask); yy=y(mask); ll=lo(mask); hh=hi(mask);
    patch(ax,[xx;flipud(xx)],[ll;flipud(hh)],c, ...
        'FaceAlpha',0.12,'EdgeColor','none','HandleVisibility','off');
    plot(ax,xx,yy,'LineStyle',lineStyle,'Color',c, ...
        'LineWidth',1.7,'DisplayName',labelText);
end
end

function cmap=interpHexPalette(hexStops,n)
rgb=zeros(numel(hexStops),3);
for i=1:numel(hexStops)
    rgb(i,:)=hex2rgb(hexStops{i});
end
x=linspace(0,1,size(rgb,1));
xi=linspace(0,1,n);
cmap=zeros(n,3);
for j=1:3
    cmap(:,j)=interp1(x,rgb(:,j),xi,'pchip');
end
cmap=max(0,min(1,cmap));
end

function rgb=hex2rgb(h)
h=char(h);
if h(1)=='#', h=h(2:end); end
rgb=[hex2dec(h(1:2)) hex2dec(h(3:4)) hex2dec(h(5:6))]/255;
end

function restoreDefaults(oldFont,oldText)
set(groot,'defaultAxesFontName',oldFont);
set(groot,'defaultTextFontName',oldText);
end
