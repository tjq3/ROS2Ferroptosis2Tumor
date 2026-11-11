function SI_Q2_2_NeighbourhoodAccess_Mechanism_FINAL()
% Curated Q2 SI-2 mechanism figure. NO NEW SIMULATION.
% a low-lineage access; b low-lineage births; c high-clearance accessibility;
% d access-burden association.

clc;close all;
root=fileparts(mfilename('fullpath'));if isempty(root),root=pwd;end
mf=findFile(root,{fullfile('Macro_Q2_NeighbourhoodAccess_Mechanism_FINAL_OUT','Macro_Q2_NeighbourhoodAccess_Mechanism_FINAL.mat'),'Macro_Q2_NeighbourhoodAccess_Mechanism_FINAL.mat'});
M=load(mf,'AUDIT');A=M.AUDIT;assert(strcmpi(A.mechanismDecision,'FINAL_MECHANISM_CONFIRMED'));
E=A.ensemble;arr=cellstr(A.arrangementNames);t=double(E.t_h(:));
c=[rgb('#5B9E99');rgb('#C99A28');rgb('#B85249')];
access=double(E.integratedLowEmptyEdgeHours)./double(E.nLowInitial);AUTC=double(E.AUTC_NV);

fig=figure('Color','w','Position',[80 60 1180 760]);tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

ax=nexttile(tl,1);hold(ax,'on');
for a=1:3,ribbon(ax,t,double(E.lowEmptyPerLow(:,:,a)),c(a,:),arr{a});end
xlabel(ax,'Time (h)');ylabel(ax,'Empty neighbours per low-$\chi_F$ viable cell','Interpreter','latex');
title(ax,'a  Low-lineage access to proliferative space');legend(ax,'Location','best','Box','off','FontSize',7.5);sty(ax);

ax=nexttile(tl,2);hold(ax,'on');
for a=1:3,Y=double(E.cumLowBirth(:,:,a))./double(E.nLowInitial(:,a)).';ribbon(ax,t,Y,c(a,:),arr{a});end
xlabel(ax,'Time (h)');ylabel(ax,'Cumulative low-lineage births / $N_{L,0}$','Interpreter','latex');
title(ax,'b  Expansion of the low-susceptibility lineage');sty(ax);

ax=nexttile(tl,3);hold(ax,'on');
for a=1:3,[m,ci]=mci(E.highClearAccessibleFrac(:,a),5000,7100+a);errorbar(ax,a,m,m-ci(1),ci(2)-m,'o','Color',c(a,:), ...
'MarkerFaceColor',c(a,:),'MarkerEdgeColor','w','LineWidth',1.2,'CapSize',6);end
xlim(ax,[.5 3.5]);xticks(ax,1:3);xticklabels(ax,arr);
ylabel(ax,{'Fraction of high-$\chi_F$ clearances','adjacent to low-$\chi_F$ viable cells'},'Interpreter','latex');
title(ax,'c  Clearance-generated vacancy accessibility');sty(ax);

ax=nexttile(tl,4);hold(ax,'on');
for a=1:3,scatter(ax,access(:,a),AUTC(:,a),22,c(a,:),'filled','MarkerFaceAlpha',.55,'MarkerEdgeColor','none','DisplayName',arr{a});end
r0=A.association.accessAUTC.r;s=A.association.accessAUTC.slope;b=A.association.accessAUTC.intercept;
xx=linspace(min(access(:),[],'omitnan'),max(access(:),[],'omitnan'),120);plot(ax,xx,b+s*xx,'-','Color',[.15 .15 .15],'LineWidth',1.3,'HandleVisibility','off');
xlabel(ax,'Integrated low$\rightarrow$empty edge-hours / $N_{L,0}$','Interpreter','latex');ylabel(ax,'$AUTC_{NV}$','Interpreter','latex');
title(ax,sprintf('d  Vacancy access predicts viable burden  (r=%.2f)',r0));legend(ax,'Location','best','Box','off','FontSize',7.5);sty(ax);

sgtitle(tl,'Q2 neighbourhood-access mechanism','FontWeight','normal','FontSize',12);
out=fullfile(root,'SI_Q2_2_NeighbourhoodAccess_Mechanism_FINAL_OUT');if ~exist(out,'dir'),mkdir(out);end
base=fullfile(out,'SI_Q2_2_NeighbourhoodAccess_Mechanism_FINAL');clean(fig);
exportgraphics(fig,[base '.png'],'Resolution',600);exportgraphics(fig,[base '.pdf'],'ContentType','vector');
try,exportgraphics(fig,[base '.svg'],'ContentType','vector');catch,end;savefig(fig,[base '.fig']);
fprintf('Saved SI-Q2-2 to %s\n',out);
end
function ribbon(ax,t,Y,c,n),m=mean(Y,2,'omitnan');se=std(Y,0,2,'omitnan')./sqrt(sum(isfinite(Y),2));lo=m-1.96*se;hi=m+1.96*se;fill(ax,[t;flipud(t)],[lo;flipud(hi)],c,'FaceAlpha',.1,'EdgeColor','none','HandleVisibility','off');plot(ax,t,m,'Color',c,'LineWidth',1.7,'DisplayName',n);end
function f=findFile(r,c),f='';for i=1:numel(c),p=fullfile(r,c{i});if isfile(p),f=p;return;end,end,error('Required FINAL MAT not found.');end
function [m,ci]=mci(x,B,s),x=double(x(:));x=x(isfinite(x));m=mean(x);rng(s,'twister');n=numel(x);z=zeros(B,1);for k=1:B,ii=randi(n,n,1);z(k)=mean(x(ii));end;z=sort(z);ci=[z(round(.025*B)) z(round(.975*B))];end
function sty(ax),ax.FontName='Arial';ax.FontSize=8.5;ax.LineWidth=.8;ax.TickDir='out';ax.Box='off';ax.Layer='top';grid(ax,'off');end
function clean(fig),try,axs=findall(fig,'Type','axes');for k=1:numel(axs),axtoolbar(axs(k),{});end,catch,end,end
function z=rgb(h),h=char(h);if h(1)=='#',h=h(2:end);end;z=[hex2dec(h(1:2)) hex2dec(h(3:4)) hex2dec(h(5:6))]/255;end
