function SI_Q2_1_SpatialArrangement_Robustness_FINAL()
% Curated Q2 SI-1. NO NEW SIMULATION.
% a Moran manipulation; b AUTC; c roughness; d matched-burden cavity.

clc;close all;
root=fileparts(mfilename('fullpath'));if isempty(root),root=pwd;end
sf=findFile(root,{fullfile('Macro_Q2_SpatialArrangement_FINAL_OUT','Macro_Q2_SpatialArrangement_FINAL.mat'),'Macro_Q2_SpatialArrangement_FINAL.mat'});
mf=findFile(root,{fullfile('Macro_Q2_NeighbourhoodAccess_Mechanism_FINAL_OUT','Macro_Q2_NeighbourhoodAccess_Mechanism_FINAL.mat'),'Macro_Q2_NeighbourhoodAccess_Mechanism_FINAL.mat'});
S=load(sf,'Q2');M=load(mf,'AUDIT');Q=S.Q2;A=M.AUDIT;
assert(Q.finalDecision.phenomenonPASS);assert(strcmpi(A.mechanismDecision,'FINAL_MECHANISM_CONFIRMED'));

E=Q.ensemble;arr=cellstr(Q.arrangementNames);AUTC=double(E.AUTC_NV);Arel=AUTC./AUTC(:,1);
rough=double(E.rough120);cav=double(A.ensemble.matchedCavity);arr3=cellstr(A.arrangementNames);

c=[rgb('#3D3D3D');rgb('#5B9E99');rgb('#C99A28');rgb('#B85249')];
c3=c(2:4,:);

fig=figure('Color','w','Position',[80 60 1180 760]);tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

ax=nexttile(tl,1);hold(ax,'on');
for a=2:4,[m,ci]=mci(E.initialMoran(:,a),5000,6100+a);errorbar(ax,a-1,m,m-ci(1),ci(2)-m,'o','Color',c(a,:), ...
'MarkerFaceColor',c(a,:),'MarkerEdgeColor','w','LineWidth',1.2,'CapSize',6);end
yline(ax,0,':','Color',[.55 .55 .55]);xlim(ax,[.5 3.5]);xticks(ax,1:3);xticklabels(ax,arr(2:4));
ylabel(ax,'Initial Moran''s I of $\chi_F$ marks','Interpreter','latex');title(ax,'a  Spatial-manipulation audit');sty(ax);

ax=nexttile(tl,2);hold(ax,'on');
for a=1:4,[m,ci]=mci(Arel(:,a),5000,6200+a);errorbar(ax,a,m,m-ci(1),ci(2)-m,'o','Color',c(a,:), ...
'MarkerFaceColor',c(a,:),'MarkerEdgeColor','w','LineWidth',1.2,'CapSize',6);end
yline(ax,1,':','Color',[.55 .55 .55]);xlim(ax,[.5 4.5]);xticks(ax,1:4);xticklabels(ax,arr);
ylabel(ax,'$AUTC_{NV}$ / paired homogeneous','Interpreter','latex');title(ax,'b  Integrated viable burden');sty(ax);

ax=nexttile(tl,3);hold(ax,'on');
for a=1:4,[m,ci]=mci(rough(:,a),5000,6300+a);errorbar(ax,a,m,m-ci(1),ci(2)-m,'o','Color',c(a,:), ...
'MarkerFaceColor',c(a,:),'MarkerEdgeColor','w','LineWidth',1.2,'CapSize',6);end
xlim(ax,[.5 4.5]);xticks(ax,1:4);xticklabels(ax,arr);ylabel(ax,'Boundary roughness at 120 h (sites)');
title(ax,'c  Secondary boundary-morphology diagnostic');sty(ax);

ax=nexttile(tl,4);hold(ax,'on');
for a=1:3,[m,ci]=mci(cav(:,a),5000,6400+a);errorbar(ax,a,m,m-ci(1),ci(2)-m,'o','Color',c3(a,:), ...
'MarkerFaceColor',c3(a,:),'MarkerEdgeColor','w','LineWidth',1.2,'CapSize',6);end
xlim(ax,[.5 3.5]);xticks(ax,1:3);xticklabels(ax,arr3);
ylabel(ax,sprintf('Cavity fraction at first $N_V\\leq%d$',A.settings.targetNV),'Interpreter','latex');
title(ax,'d  Fragmentation at matched viable burden');sty(ax);

sgtitle(tl,'Q2 spatial-arrangement robustness and morphology controls','FontWeight','normal','FontSize',12);
out=fullfile(root,'SI_Q2_1_SpatialArrangement_Robustness_FINAL_OUT');if ~exist(out,'dir'),mkdir(out);end
base=fullfile(out,'SI_Q2_1_SpatialArrangement_Robustness_FINAL');clean(fig);
exportgraphics(fig,[base '.png'],'Resolution',600);exportgraphics(fig,[base '.pdf'],'ContentType','vector');
try,exportgraphics(fig,[base '.svg'],'ContentType','vector');catch,end;savefig(fig,[base '.fig']);
fprintf('Saved SI-Q2-1 to %s\n',out);
end
function f=findFile(r,c),f='';for i=1:numel(c),p=fullfile(r,c{i});if isfile(p),f=p;return;end,end,error('Required FINAL MAT not found.');end
function [m,ci]=mci(x,B,s),x=double(x(:));x=x(isfinite(x));m=mean(x);rng(s,'twister');n=numel(x);z=zeros(B,1);for k=1:B,ii=randi(n,n,1);z(k)=mean(x(ii));end;z=sort(z);ci=[z(round(.025*B)) z(round(.975*B))];end
function sty(ax),ax.FontName='Arial';ax.FontSize=8.5;ax.LineWidth=.8;ax.TickDir='out';ax.Box='off';ax.Layer='top';grid(ax,'off');end
function clean(fig),try,axs=findall(fig,'Type','axes');for k=1:numel(axs),axtoolbar(axs(k),{});end,catch,end,end
function z=rgb(h),h=char(h);if h(1)=='#',h=h(2:end);end;z=[hex2dec(h(1:2)) hex2dec(h(3:4)) hex2dec(h(5:6))]/255;end
