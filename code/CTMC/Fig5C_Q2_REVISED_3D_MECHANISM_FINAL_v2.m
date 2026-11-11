function Fig5C_Q2_REVISED_3D_MECHANISM_FINAL_v2()
% =========================================================================
% Fig5C_Q2_REVISED_3D_MECHANISM_FINAL_v2
%
% FIXED VERSION: uses the exact field names saved by the completed Q2 FINAL
% simulations. No fuzzy field-name searching and no new simulation.
%
% REQUIRED INPUTS
%   Macro_Q2_SpatialArrangement_FINAL_OUT/Macro_Q2_SpatialArrangement_FINAL.mat
%   Macro_Q2_NeighbourhoodAccess_Mechanism_FINAL_OUT/
%       Macro_Q2_NeighbourhoodAccess_Mechanism_FINAL.mat
%
% FIGURE LOGIC
%   C-i   3D Dispersed / Random / Clustered matched spatial counterfactuals
%   C-ii  Viable tumour burden dynamics + AUTC inset
%   C-iii Clearance-generated vacancy accessibility + low-lineage expansion
% =========================================================================

clc; close all;

%% 0. LOCATE FINAL FILES
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir = pwd; end

spatialFile = findFirst({ ...
    fullfile(scriptDir,'Macro_Q2_SpatialArrangement_FINAL_OUT', ...
        'Macro_Q2_SpatialArrangement_FINAL.mat'), ...
    fullfile(scriptDir,'Macro_Q2_SpatialArrangement_FINAL.mat')});

mechanismFile = findFirst({ ...
    fullfile(scriptDir,'Macro_Q2_NeighbourhoodAccess_Mechanism_FINAL_OUT', ...
        'Macro_Q2_NeighbourhoodAccess_Mechanism_FINAL.mat'), ...
    fullfile(scriptDir,'Macro_Q2_NeighbourhoodAccess_Mechanism_FINAL.mat')});

assert(~isempty(spatialFile),'Cannot find Macro_Q2_SpatialArrangement_FINAL.mat.');
assert(~isempty(mechanismFile),'Cannot find Macro_Q2_NeighbourhoodAccess_Mechanism_FINAL.mat.');

S = load(spatialFile,'Q2');
M = load(mechanismFile,'AUDIT');
assert(isfield(S,'Q2'),'Spatial FINAL MAT must contain Q2.');
assert(isfield(M,'AUDIT'),'Mechanism FINAL MAT must contain AUDIT.');
Q2 = S.Q2;
AUDIT = M.AUDIT;

assert(Q2.finalDecision.phenomenonPASS,'Q2 FINAL spatial phenomenon is not marked PASS.');
assert(strcmpi(AUDIT.mechanismDecision,'FINAL_MECHANISM_CONFIRMED'), ...
    'Q2 FINAL neighbourhood-access mechanism is not confirmed.');

%% 1. READ EXACT SAVED FIELDS
E  = Q2.ensemble;
P  = Q2.patternExample;
ME = AUDIT.ensemble;

% Exact fields written by Macro_Q2_SpatialArrangement_FINAL.m
assert(isfield(E,'t_h'),'Q2.ensemble.t_h is missing.');
assert(isfield(E,'NVcurve'),'Q2.ensemble.NVcurve is missing.');
assert(isfield(E,'AUTC_NV'),'Q2.ensemble.AUTC_NV is missing.');

% Exact fields written into Q2.patternExample
assert(isfield(P,'mapDispersed'),'Q2.patternExample.mapDispersed is missing.');
assert(isfield(P,'mapRandom'),'Q2.patternExample.mapRandom is missing.');
assert(isfield(P,'mapClustered'),'Q2.patternExample.mapClustered is missing.');
assert(isfield(P,'initialMask'),'Q2.patternExample.initialMask is missing.');
assert(isfield(P,'metrics'),'Q2.patternExample.metrics is missing.');

% Exact fields written by Macro_Q2_NeighbourhoodAccess_Mechanism_FINAL.m
assert(isfield(ME,'t_h'),'AUDIT.ensemble.t_h is missing.');
assert(isfield(ME,'highClearAccessibleFrac'),'AUDIT.ensemble.highClearAccessibleFrac is missing.');
assert(isfield(ME,'cumLowBirth'),'AUDIT.ensemble.cumLowBirth is missing.');
assert(isfield(ME,'nLowInitial'),'AUDIT.ensemble.nLowInitial is missing.');

arr4 = cellstr(Q2.arrangementNames);
arr3 = cellstr(AUDIT.arrangementNames);
assert(numel(arr4)==4 && strcmpi(arr4{1},'Homogeneous') && strcmpi(arr4{2},'Dispersed') && ...
    strcmpi(arr4{3},'Random') && strcmpi(arr4{4},'Clustered'));
assert(numel(arr3)==3 && strcmpi(arr3{1},'Dispersed') && strcmpi(arr3{2},'Random') && ...
    strcmpi(arr3{3},'Clustered'));

t = double(E.t_h(:));
NV = double(E.NVcurve);            % nTime x nRep x 4
AUTC = double(E.AUTC_NV);          % nRep x 4

tm = double(ME.t_h(:));
Hacc = double(ME.highClearAccessibleFrac); % nRep x 3
cumBirth = double(ME.cumLowBirth);         % nTime x nRep x 3
nLow0 = double(ME.nLowInitial);             % nRep x 3

assert(size(NV,1)==numel(t) && size(NV,3)==4, ...
    'Unexpected Q2.ensemble.NVcurve size; expected nTime x nRep x 4.');
assert(size(AUTC,2)==4,'Unexpected Q2.ensemble.AUTC_NV size; expected nRep x 4.');
assert(size(cumBirth,1)==numel(tm) && size(cumBirth,3)==3, ...
    'Unexpected AUDIT.ensemble.cumLowBirth size; expected nTime x nRep x 3.');
assert(isequal(size(Hacc),size(nLow0)), ...
    'highClearAccessibleFrac and nLowInitial must both be nRep x 3.');
assert(numel(t)==numel(tm) && max(abs(t-tm))<1e-12, ...
    'Spatial and mechanism FINAL time grids do not match.');

%% 2. STYLE
cHom  = hex2rgb('#3F3F3F');
cDisp = hex2rgb('#6AA6A1');
cRand = hex2rgb('#D3A536');
cClus = hex2rgb('#B4534B');
cLow  = hex2rgb('#76AAA5');
cHigh = hex2rgb('#C95D52');
cc4 = [cHom;cDisp;cRand;cClus];
cc3 = [cDisp;cRand;cClus];

%% 3. FIGURE
fig = figure('Color','w','Position',[40 40 1680 900]);
tl = tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');

% -------------------------------------------------------------------------
% C-i: REAL 3D matched spatial counterfactuals
% -------------------------------------------------------------------------
mapNames = {'mapDispersed','mapRandom','mapClustered'};
metricNames = {'dispersed','random','clustered'};
displayNames = {'Dispersed','Random','Clustered'};

for a=1:3
    ax = nexttile(tl,a);
    plotReal3DPattern(ax,double(P.(mapNames{a})),logical(P.initialMask), ...
        double(Q2.hazard.chiLow),double(Q2.hazard.chiHigh),cLow,cHigh);
    moranI = double(P.metrics.(metricNames{a}).moranI);
    title(ax,sprintf('%s\nMoran''s I = %+.2f',displayNames{a},moranI), ...
        'FontName','Arial','FontWeight','normal','FontSize',11);
    if a==1
        text(ax,-0.10,1.06,'C-i','Units','normalized','FontName','Arial', ...
            'FontSize',13,'FontWeight','bold');
    end
end

% -------------------------------------------------------------------------
% C-ii: macroscopic consequence
% -------------------------------------------------------------------------
ax4 = nexttile(tl,4); hold(ax4,'on');
for a=1:4
    plotRibbon95(ax4,t,squeeze(NV(:,:,a)),cc4(a,:),arr4{a});
end
xlabel(ax4,'Time (h)');
ylabel(ax4,'Viable tumour cells, $N_V$','Interpreter','latex');
title(ax4,'C-ii  Spatial organization reshapes viable-tumour expansion', ...
    'FontWeight','normal','FontSize',11);
legend(ax4,'Location','northwest','Box','off','FontSize',8);
styleAxis(ax4);

% AUTC inset, paired normalization to homogeneous
pos = ax4.Position;
axI = axes(fig,'Position',[pos(1)+0.58*pos(3),pos(2)+0.12*pos(4),0.34*pos(3),0.31*pos(4)]);
hold(axI,'on');
Arel = AUTC ./ AUTC(:,1);
for a=1:4
    [m,ci] = bootMeanCI(Arel(:,a),5000,8100+a);
    errorbar(axI,a,m,m-ci(1),ci(2)-m,'o','Color',cc4(a,:), ...
        'MarkerFaceColor',cc4(a,:),'MarkerEdgeColor','w', ...
        'LineWidth',1,'CapSize',4,'MarkerSize',5.5);
end
yline(axI,1,':','Color',[0.55 0.55 0.55]);
xlim(axI,[0.5 4.5]); xticks(axI,1:4); xticklabels(axI,{'H','D','R','C'});
ylabel(axI,'AUTC / H','FontSize',7.2);
title(axI,'Integrated viable burden','FontWeight','normal','FontSize',7.7);
styleAxisSmall(axI);

% -------------------------------------------------------------------------
% C-iii left: clearance-generated vacancy accessibility
% -------------------------------------------------------------------------
ax5 = nexttile(tl,5); hold(ax5,'on');
for a=1:3
    [m,ci] = bootMeanCI(Hacc(:,a),5000,8200+a);
    errorbar(ax5,a,m,m-ci(1),ci(2)-m,'o','Color',cc3(a,:), ...
        'MarkerFaceColor',cc3(a,:),'MarkerEdgeColor','w', ...
        'LineWidth',1.25,'CapSize',7,'MarkerSize',7.5);
end
xlim(ax5,[0.55 3.45]); xticks(ax5,1:3); xticklabels(ax5,arr3);
ylabel(ax5,{'Fraction of high-$\chi_F$ clearances','adjacent to low-$\chi_F$ viable cells'}, ...
    'Interpreter','latex');
title(ax5,'C-iii  Clearance-generated vacancy accessibility', ...
    'FontWeight','normal','FontSize',11);
styleAxis(ax5);

% -------------------------------------------------------------------------
% C-iii right: low-susceptibility lineage expansion
% -------------------------------------------------------------------------
ax6 = nexttile(tl,6); hold(ax6,'on');
for a=1:3
    Y = squeeze(cumBirth(:,:,a));       % nTime x nRep
    denom = nLow0(:,a).';               % 1 x nRep
    Ynorm = Y ./ denom;
    plotRibbon95(ax6,t,Ynorm,cc3(a,:),arr3{a});
end
xlabel(ax6,'Time (h)');
ylabel(ax6,'Cumulative low-lineage births / $N_{L,0}$','Interpreter','latex');
title(ax6,'Low-susceptibility survivor-lineage expansion', ...
    'FontWeight','normal','FontSize',11);
styleAxis(ax6);
text(ax6,0.04,0.94, ...
    'clustering  \rightarrow  less vacancy access  \rightarrow  less lineage expansion', ...
    'Units','normalized','FontName','Arial','FontSize',8.5, ...
    'Color',[0.32 0.32 0.32],'VerticalAlignment','top');

sgtitle(tl, ...
    'Q2: spatial organization of ferroptotic susceptibility controls tumour expansion through vacancy access', ...
    'FontName','Arial','FontWeight','normal','FontSize',14);

%% 4. EXPORT
outDir = fullfile(scriptDir,'Fig5C_Q2_REVISED_3D_MECHANISM_FINAL_v2_OUT');
if ~exist(outDir,'dir'), mkdir(outDir); end
base = fullfile(outDir,'Fig5C_Q2_REVISED_3D_MECHANISM_FINAL_v2');
removeAxesToolbars(fig);
exportgraphics(fig,[base '.png'],'Resolution',600);
exportgraphics(fig,[base '.pdf'],'ContentType','vector');
try, exportgraphics(fig,[base '.svg'],'ContentType','vector'); catch, end
savefig(fig,[base '.fig']);

fprintf('\n============================================================\n');
fprintf('FIG5C REVISED 3D/MECHANISM LAYOUT: SUCCESS\n');
fprintf('NVcurve size     : %s\n',mat2str(size(NV)));
fprintf('AUTC size        : %s\n',mat2str(size(AUTC)));
fprintf('cumLowBirth size : %s\n',mat2str(size(cumBirth)));
fprintf('Saved to         : %s\n',outDir);
fprintf('============================================================\n');
end

%% =========================================================================
function f=findFirst(cands)
f='';
for k=1:numel(cands)
    if isfile(cands{k}), f=cands{k}; return; end
end
end

function plotReal3DPattern(ax,map,mask,chiLow,chiHigh,cLow,cHigh)
% Uses the ACTUAL full 3D susceptibility map saved by Q2.patternExample.
assert(isequal(size(map),size(mask)), ...
    'patternExample map and initialMask have different sizes.');
idx=find(mask); vals=map(idx);
dLow=abs(vals-chiLow); dHigh=abs(vals-chiHigh);
isHigh=dHigh<dLow; isLow=~isHigh;
[x,y,z]=ind2sub(size(map),idx);
ctr=(size(map)+1)/2;
xc=x-ctr(1); yc=y-ctr(2); zc=z-ctr(3);
hold(ax,'on');
scatter3(ax,yc(isLow),xc(isLow),zc(isLow),34,'o','filled', ...
    'MarkerFaceColor',cLow,'MarkerEdgeColor','w','LineWidth',0.25, ...
    'MarkerFaceAlpha',0.70,'MarkerEdgeAlpha',0.25);
scatter3(ax,yc(isHigh),xc(isHigh),zc(isHigh),42,'o','filled', ...
    'MarkerFaceColor',cHigh,'MarkerEdgeColor','w','LineWidth',0.30, ...
    'MarkerFaceAlpha',0.90,'MarkerEdgeAlpha',0.35);
hold(ax,'off');
axis(ax,'equal'); axis(ax,'off'); view(ax,[-38 22]); camproj(ax,'perspective');
lim=max(abs([xc;yc;zc]))+0.8;
xlim(ax,[-lim lim]); ylim(ax,[-lim lim]); zlim(ax,[-lim lim]);
end

function plotRibbon95(ax,t,Y,c,name)
m=mean(Y,2,'omitnan'); n=sum(isfinite(Y),2);
se=std(Y,0,2,'omitnan')./sqrt(max(n,1));
lo=m-1.96*se; hi=m+1.96*se;
fill(ax,[t;flipud(t)],[lo;flipud(hi)],c,'FaceAlpha',0.10, ...
    'EdgeColor','none','HandleVisibility','off');
plot(ax,t,m,'Color',c,'LineWidth',1.7,'DisplayName',name);
end

function [m,ci]=bootMeanCI(x,nBoot,seed)
x=double(x(:)); x=x(isfinite(x)); m=mean(x);
rng(seed,'twister'); n=numel(x); b=zeros(nBoot,1);
for k=1:nBoot
    ii=randi(n,n,1); b(k)=mean(x(ii));
end
b=sort(b);
ci=[b(max(1,round(0.025*nBoot))) b(min(nBoot,round(0.975*nBoot)))];
end

function styleAxis(ax)
ax.FontName='Arial'; ax.FontSize=8.6; ax.LineWidth=0.85;
ax.TickDir='out'; ax.TickLength=[0.014 0.014]; ax.Box='off'; ax.Layer='top';
grid(ax,'off');
end

function styleAxisSmall(ax)
ax.FontName='Arial'; ax.FontSize=7.2; ax.LineWidth=0.70;
ax.TickDir='out'; ax.TickLength=[0.018 0.018]; ax.Box='off'; ax.Layer='top';
grid(ax,'off');
end

function removeAxesToolbars(fig)
try
    axs=findall(fig,'Type','axes');
    for k=1:numel(axs), axtoolbar(axs(k),{}); end
catch
end
end

function rgb=hex2rgb(h)
h=char(h); if h(1)=='#', h=h(2:end); end
rgb=[hex2dec(h(1:2)) hex2dec(h(3:4)) hex2dec(h(5:6))]/255;
end
