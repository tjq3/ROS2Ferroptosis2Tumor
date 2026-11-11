function CTMC_SpatialArrangementGenerality()
% Curated Q2 SI-3 condition/generalization figure. NO NEW SIMULATION.
% a Delta AUTC; b Delta Pext; c Delta vacancy access; d cross-condition alignment.
% Discrete cells only: no interpolated contour surface.

clc;close all;
root=fileparts(mfilename('fullpath'));if isempty(root),root=pwd;end
mf=findFile(root,{fullfile('Macro_Q2_TargetedConditionMap_FINAL_OUT','Macro_Q2_TargetedConditionMap_FINAL.mat'),'Macro_Q2_TargetedConditionMap_FINAL.mat'});
S=load(mf,'MAP');MAP=S.MAP;x=double(MAP.grid.meanChi(:)).';y=double(MAP.grid.deltaF(:)).';M=MAP.summary;
ax0=double(MAP.grid.anchorChiBar);ay0=double(MAP.grid.anchorDeltaF);
cmap=divmap(257);cseq=seqmap(257);

fig=figure('Color','w','Position',[80 60 1180 760]);tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

ax=nexttile(tl,1);heat(ax,x,y,double(M.dAUTC_pct).',cmap);hold(ax,'on');robust(ax,x,y,M.robustAUTC);
plot(ax,ax0,ay0,'p','MarkerSize',11,'MarkerFaceColor','k','MarkerEdgeColor','w');
xlabel(ax,'Mean commitment pressure, $\bar{\chi}_F$','Interpreter','latex');ylabel(ax,'Susceptibility contrast, $\delta_F$','Interpreter','latex');
title(ax,'a  Spatial effect on integrated viable burden');cb=colorbar(ax);cb.Label.String='$\Delta AUTC_{NV}^{C-R}$ (\%)';cb.Label.Interpreter='latex';sty(ax);

ax=nexttile(tl,2);heat(ax,x,y,double(M.dPext).',cmap);hold(ax,'on');plot(ax,ax0,ay0,'p','MarkerSize',11,'MarkerFaceColor','k','MarkerEdgeColor','w');
xlabel(ax,'Mean commitment pressure, $\bar{\chi}_F$','Interpreter','latex');ylabel(ax,'Susceptibility contrast, $\delta_F$','Interpreter','latex');
title(ax,'b  Spatial effect on extinction probability');cb=colorbar(ax);cb.Label.String='$\Delta P_{ext}^{C-R}(168\,h)$';cb.Label.Interpreter='latex';sty(ax);

ax=nexttile(tl,3);heat(ax,x,y,double(M.dAccess_pct).',cmap);hold(ax,'on');plot(ax,ax0,ay0,'p','MarkerSize',11,'MarkerFaceColor','k','MarkerEdgeColor','w');
xlabel(ax,'Mean commitment pressure, $\bar{\chi}_F$','Interpreter','latex');ylabel(ax,'Susceptibility contrast, $\delta_F$','Interpreter','latex');
title(ax,'c  Spatial effect on low-lineage vacancy access');cb=colorbar(ax);cb.Label.String='$\Delta access^{C-R}$ (\%)';cb.Label.Interpreter='latex';sty(ax);

ax=nexttile(tl,4);hold(ax,'on');
D=repmat(y,numel(x),1);mask=isfinite(M.dAccess_pct)&isfinite(M.dAUTC_pct)&D>0;
X=double(M.dAccess_pct(mask));Y=double(M.dAUTC_pct(mask));C=double(D(mask));
scatter(ax,X,Y,48,C,'filled','MarkerEdgeColor','w','LineWidth',.4);colormap(ax,cseq);
cb=colorbar(ax);cb.Label.String='$\delta_F$';cb.Label.Interpreter='latex';
r0=MAP.mechanismAlignment.r_access_to_AUTC;s=MAP.mechanismAlignment.slope_access_to_AUTC;b=mean(Y)-s*mean(X);
xx=linspace(min(X),max(X),120);plot(ax,xx,b+s*xx,'-','Color',[.15 .15 .15],'LineWidth',1.3);
xlabel(ax,'$\Delta$ low-lineage vacancy access$^{C-R}$ (\%)','Interpreter','latex');
ylabel(ax,'$\Delta AUTC_{NV}^{C-R}$ (\%)','Interpreter','latex');
title(ax,sprintf('d  Mechanism-outcome alignment across conditions  (r=%.2f)',r0));sty(ax);

sgtitle(tl,'Q2 generality across mean ferroptotic pressure and susceptibility contrast','FontWeight','normal','FontSize',12);
out=fullfile(root,'SI_Q2_3_ConditionGenerality_FINAL_OUT');if ~exist(out,'dir'),mkdir(out);end
base=fullfile(out,'SI_Q2_3_ConditionGenerality_FINAL');clean(fig);
exportgraphics(fig,[base '.png'],'Resolution',600);exportgraphics(fig,[base '.pdf'],'ContentType','vector');
try,exportgraphics(fig,[base '.svg'],'ContentType','vector');catch,end;savefig(fig,[base '.fig']);
fprintf('Saved SI-Q2-3 to %s\n',out);
end
function heat(ax,x,y,Z,cmap)
xe=edges(x);ye=edges(y);hold(ax,'on');v=Z(isfinite(Z));mx=max(abs(v),[],'omitnan');if isempty(mx)||~isfinite(mx)||mx<eps,mx=1;end
for j=1:numel(y),for i=1:numel(x),if isfinite(Z(j,i)),patch(ax,[xe(i) xe(i+1) xe(i+1) xe(i)],[ye(j) ye(j) ye(j+1) ye(j+1)],Z(j,i),'EdgeColor','none');end,end,end
xlim(ax,[xe(1) xe(end)]);ylim(ax,[ye(1) ye(end)]);caxis(ax,[-mx mx]);colormap(ax,cmap);
end
function e=edges(c),c=c(:).';m=(c(1:end-1)+c(2:end))/2;e=[c(1)-(m(1)-c(1)) m c(end)+(c(end)-m(end))];end
function robust(ax,x,y,R),for j=1:numel(y),for i=1:numel(x),if R(i,j)==-1,plot(ax,x(i),y(j),'o','MarkerSize',4,'MarkerFaceColor',[.1 .1 .1],'MarkerEdgeColor','w','LineWidth',.4);elseif R(i,j)==1,plot(ax,x(i),y(j),'^','MarkerSize',4.5,'MarkerFaceColor',[.1 .1 .1],'MarkerEdgeColor','w','LineWidth',.4);end,end,end,end
function f=findFile(r,c),f='';for i=1:numel(c),p=fullfile(r,c{i});if isfile(p),f=p;return;end,end,error('Required FINAL MAT not found.');end
function z=divmap(n),a=rgb('#2F7EA8');m=rgb('#F7F5EF');b=rgb('#C85C4A');n1=ceil(n/2);n2=n-n1+1;u=zeros(n1,3);v=zeros(n2,3);for k=1:3,u(:,k)=linspace(a(k),m(k),n1);v(:,k)=linspace(m(k),b(k),n2);end;z=[u;v(2:end,:)];end
function z=seqmap(n),a=rgb('#EAF4EF');m=rgb('#7FC8C6');b=rgb('#1F6FA8');n1=ceil(n/2);n2=n-n1+1;u=zeros(n1,3);v=zeros(n2,3);for k=1:3,u(:,k)=linspace(a(k),m(k),n1);v(:,k)=linspace(m(k),b(k),n2);end;z=[u;v(2:end,:)];end
function sty(ax),ax.FontName='Arial';ax.FontSize=8.5;ax.LineWidth=.8;ax.TickDir='out';ax.Box='off';ax.Layer='top';grid(ax,'off');end
function clean(fig),try,axs=findall(fig,'Type','axes');for k=1:numel(axs),axtoolbar(axs(k),{});end,catch,end,end
function z=rgb(h),h=char(h);if h(1)=='#',h=h(2:end);end;z=[hex2dec(h(1:2)) hex2dec(h(3:4)) hex2dec(h(5:6))]/255;end
