%% LPO_ThresholdAndWeightRobustness.m
% =========================================================================
% Mesoscopic SI robustness:
% Fate-boundary M_L x relative source weighting eta = omega_S / omega_T
%
% PURPOSE
% -------
% Test whether the main mesoscopic conclusions depend on the single
% operating point M_L = 0.5 and omega_S/omega_T = 1.
%
% DESIGN
% ------
% 1) Reuse FINAL Fig.4B locked microscopic trajectory banks.
% 2) Reuse FINAL Fig.4B locked TD-SD pairing: idxT_pair / idxS_pair.
% 3) Hold total source coupling fixed:
%       Omega = omega_T + omega_S
%    and vary only relative weighting eta = omega_S/omega_T:
%       omega_T = Omega/(1+eta)
%       omega_S = eta*Omega/(1+eta)
% 4) Keep rho_L, L0, dt and Tobs fixed.
% 5) Use exact within-bin LPO/FPT crossing.
%
% OUTPUT PANELS
% -------------
% A. DeltaRMST_T: >0 means CIR commits earlier than GammaOU.
% B. DeltaRMST_S: >0 means Full commits earlier than Uniform.
% C. Mean Joint-required fraction in Full background.
% D. Mean Resistant fraction in Full background.
% =========================================================================

clear; clc; close all;

%% ============================== USER SETTINGS =============================

scriptDir = fileparts(mfilename('fullpath'));
matFile = fullfile(scriptDir,'Fig4B_fate_routes_final','Fig4B_fateRoutes_data.mat');

outDir = fullfile(fileparts(matFile), ...
    'SI_Meso_ML_eta_robustness_FINAL');

if ~exist(outDir,'dir')
    mkdir(outDir);
end

Mvec   = 0.40:0.025:0.60;
etaVec = 2.^linspace(-2,2,9);

Mref   = 0.50;
etaRef = 1.00;

EXPORT_DPI = 600;

%% ================================ STYLE ==================================

set(groot,'defaultFigureColor','w');
set(groot,'defaultAxesFontName','Arial');
set(groot,'defaultTextFontName','Arial');
set(groot,'defaultAxesFontSize',8);
set(groot,'defaultTextFontSize',8);
set(groot,'defaultAxesLineWidth',0.75);
set(groot,'defaultAxesTickDir','out');
set(groot,'defaultAxesTickLength',[0.014 0.014]);
set(groot,'defaultAxesBox','off');
set(groot,'defaultAxesLayer','top');

%% ================================ LOAD ===================================

if ~isfile(matFile)
    error('Cannot find input MAT file:\n%s',matFile);
end

fprintf('Loading:\n%s\n\n',matFile);
D = load(matFile);

if ~isfield(D,'inputs')
    error('The MAT file does not contain top-level field D.inputs.');
end
if ~isfield(D,'meta')
    error('The MAT file does not contain top-level field D.meta.');
end

requiredInputs = {'IT_C','IT_G','IS_U','IS_F','idxT_pair','idxS_pair'};
for ii = 1:numel(requiredInputs)
    if ~isfield(D.inputs,requiredInputs{ii})
        error('Missing D.inputs.%s.',requiredInputs{ii});
    end
end

if ~isfield(D.meta,'meso')
    error('Missing D.meta.meso.');
end
if ~isfield(D.meta,'dtBin_min')
    error('Missing D.meta.dtBin_min.');
end
if ~isfield(D.meta,'Tobs_min')
    error('Missing D.meta.Tobs_min.');
end

%% -------------------------- RAW LOCKED BANKS -----------------------------

IT_C_bank = D.inputs.IT_C;
IT_G_bank = D.inputs.IT_G;
IS_U_bank = D.inputs.IS_U;
IS_F_bank = D.inputs.IS_F;

idxT_pair = D.inputs.idxT_pair(:);
idxS_pair = D.inputs.idxS_pair(:);

%% ----------------------------- META --------------------------------------

dt    = D.meta.dtBin_min;
Tobs  = D.meta.Tobs_min;
meso0 = D.meta.meso;

omegaT0 = getStructScalar(meso0,{'omegaT','omega_T'},0.0475);
omegaS0 = getStructScalar(meso0,{'omegaS','omega_S'},0.0475);
rho0    = getStructScalar(meso0,{'rhoL','rho_L'},0.1);
L0      = getStructScalar(meso0,{'L0','L_0'},0);

if isfield(meso0,'ML')
    ML0 = meso0.ML;
elseif isfield(meso0,'M_L')
    ML0 = meso0.M_L;
else
    ML0 = 0.5;
end

fprintf('Locked source-bank dimensions:\n');
fprintf('  TD CIR bank        = %d x %d\n',size(IT_C_bank,1),size(IT_C_bank,2));
fprintf('  TD GammaOU bank    = %d x %d\n',size(IT_G_bank,1),size(IT_G_bank,2));
fprintf('  SD Uniform bank    = %d x %d\n',size(IS_U_bank,1),size(IS_U_bank,2));
fprintf('  SD Full bank       = %d x %d\n',size(IS_F_bank,1),size(IS_F_bank,2));

fprintf('\nLocked pairing / mesoscopic configuration:\n');
fprintf('  nPairs             = %d\n',numel(idxT_pair));
fprintf('  dt                 = %.3f min\n',dt);
fprintf('  Tobs               = %.1f min\n',Tobs);
fprintf('  M_L (reference)    = %.4f\n',ML0);
fprintf('  rho_L              = %.5f min^-1\n',rho0);
fprintf('  omega_T(reference) = %.6f min^-1\n',omegaT0);
fprintf('  omega_S(reference) = %.6f min^-1\n',omegaS0);
fprintf('  L0                 = %.4f\n\n',L0);

%% -------------------------- INDEX AUDIT ----------------------------------

if numel(idxT_pair) ~= numel(idxS_pair)
    error('idxT_pair and idxS_pair have different lengths.');
end

if any(idxT_pair<1) || any(idxT_pair>size(IT_C_bank,1)) || ...
        any(idxT_pair>size(IT_G_bank,1))
    error('idxT_pair exceeds the TD-bank dimensions.');
end

if any(idxS_pair<1) || any(idxS_pair>size(IS_U_bank,1)) || ...
        any(idxS_pair>size(IS_F_bank,1))
    error('idxS_pair exceeds the SD-bank dimensions.');
end

%% ==================== RECONSTRUCT LOCKED CELL PAIRING ====================

IT_C = IT_C_bank(idxT_pair,:);
IT_G = IT_G_bank(idxT_pair,:);
IS_U = IS_U_bank(idxS_pair,:);
IS_F = IS_F_bank(idxS_pair,:);

nPairs = size(IT_C,1);

nBin = min([size(IT_C,2),size(IT_G,2),size(IS_U,2),size(IS_F,2)]);
IT_C = IT_C(:,1:nBin);
IT_G = IT_G(:,1:nBin);
IS_U = IS_U(:,1:nBin);
IS_F = IS_F(:,1:nBin);

if abs(nBin*dt - Tobs) > 1e-8
    warning(['nBin*dt = %.6f min, whereas saved Tobs = %.6f min. ' ...
             'The saved Tobs is retained.'],nBin*dt,Tobs);
end

fprintf('Locked Fig.4B cell pairing restored successfully:\n');
fprintf('  IT_C = %d x %d\n',size(IT_C,1),size(IT_C,2));
fprintf('  IT_G = %d x %d\n',size(IT_G,1),size(IT_G,2));
fprintf('  IS_U = %d x %d\n',size(IS_U,1),size(IS_U,2));
fprintf('  IS_F = %d x %d\n\n',size(IS_F,1),size(IS_F,2));

%% ======================= FIX TOTAL COUPLING BUDGET =======================

Omega = omegaT0 + omegaS0;

fprintf('Fixed total source-coupling budget:\n');
fprintf('  Omega = omega_T + omega_S = %.6f min^-1\n\n',Omega);

nM   = numel(Mvec);
nEta = numel(etaVec);

% condition order: CU, GU, CF, GF
Pcommit   = nan(nM,nEta,4);
RMST      = nan(nM,nEta,4);
routeFrac = nan(nM,nEta,4,5);

DeltaRMST_T = nan(nM,nEta);  % >0: CIR earlier than GammaOU
DeltaRMST_S = nan(nM,nEta);  % >0: Full earlier than Uniform
JointFull   = nan(nM,nEta);
ResFull     = nan(nM,nEta);

%% =============================== MAIN SCAN ===============================

fprintf('Starting M_L x eta robustness scan...\n\n');

for iM = 1:nM
    ML = Mvec(iM);

    for ie = 1:nEta
        eta = etaVec(ie);

        % Keep TOTAL source coupling fixed; vary only T:S allocation.
        omegaT = Omega/(1+eta);
        omegaS = eta*Omega/(1+eta);

        fprintf(['M_L = %.2f | eta = %.2f | ' ...
                 'omega_T = %.5f | omega_S = %.5f\n'], ...
                 ML,eta,omegaT,omegaS);

        RCU = evaluateRoutesExact(IT_C,IS_U,omegaT,omegaS,rho0,ML,L0,dt,Tobs);
        RGU = evaluateRoutesExact(IT_G,IS_U,omegaT,omegaS,rho0,ML,L0,dt,Tobs);
        RCF = evaluateRoutesExact(IT_C,IS_F,omegaT,omegaS,rho0,ML,L0,dt,Tobs);
        RGF = evaluateRoutesExact(IT_G,IS_F,omegaT,omegaS,rho0,ML,L0,dt,Tobs);

        RR = {RCU,RGU,RCF,RGF};

        for c = 1:4
            Pcommit(iM,ie,c) = RR{c}.Pcommit;
            RMST(iM,ie,c)    = RR{c}.RMST;
            routeFrac(iM,ie,c,:) = reshape(RR{c}.routeFrac,1,1,1,5);
        end

        % >0 means CIR commits earlier than GammaOU.
        DeltaRMST_T(iM,ie) = 0.5 * ...
            ((RGU.RMST-RCU.RMST) + (RGF.RMST-RCF.RMST));

        % >0 means Full commits earlier than Uniform.
        DeltaRMST_S(iM,ie) = 0.5 * ...
            ((RCU.RMST-RCF.RMST) + (RGU.RMST-RGF.RMST));

        JointFull(iM,ie) = 0.5 * (RCF.routeFrac(4) + RGF.routeFrac(4));
        ResFull(iM,ie)   = 0.5 * (RCF.routeFrac(5) + RGF.routeFrac(5));
    end
end

%% ====================== REFERENCE-POINT REPRODUCTION =====================

[~,iMref]   = min(abs(Mvec-Mref));
[~,iEtaref] = min(abs(etaVec-etaRef));

fprintf('\n================ REFERENCE-POINT AUDIT ================\n');
fprintf('Closest grid point: M_L=%.3f, eta=%.3f\n',Mvec(iMref),etaVec(iEtaref));

for c = 1:4
    fprintf('Condition %d | Pcommit=%.4f | RMST=%.4f min\n', ...
        c,Pcommit(iMref,iEtaref,c),RMST(iMref,iEtaref,c));
end

fprintf('DeltaRMST_T = %+0.4f min  (>0: CIR earlier)\n',DeltaRMST_T(iMref,iEtaref));
fprintf('DeltaRMST_S = %+0.4f min  (>0: Full earlier)\n',DeltaRMST_S(iMref,iEtaref));
fprintf('JointFull   = %.4f\n',JointFull(iMref,iEtaref));
fprintf('ResFull     = %.4f\n',ResFull(iMref,iEtaref));

% Optional exact reproduction audit against saved Fig.4B result.
if isfield(D,'R') && numel(D.R)>=4 && ...
        abs(Mvec(iMref)-ML0)<1e-12 && ...
        abs(etaVec(iEtaref)-omegaS0/omegaT0)<1e-12

    fprintf('\nSaved Fig.4B reproduction differences:\n');
    for c = 1:4
        dP = Pcommit(iMref,iEtaref,c) - D.R(c).Pcross60;
        dR = RMST(iMref,iEtaref,c)    - D.R(c).RMST60;
        rf = squeeze(routeFrac(iMref,iEtaref,c,:)).';
        dRoute = max(abs(rf-D.R(c).routeFrac));
        fprintf('  condition %d: dP=%+.3e | dRMST=%+.3e | max dRoute=%.3e\n', ...
            c,dP,dR,dRoute);
    end
end
fprintf('=======================================================\n');

%% =========================== GRID-LEVEL SUMMARY ==========================

fracTemporalCIREarlier = mean(DeltaRMST_T(:)>0);
fracSpatialFullEarlier = mean(DeltaRMST_S(:)>0);
fracJointPresent       = mean(JointFull(:)>0.05);

fprintf('\n================ ROBUSTNESS SUMMARY ====================\n');
fprintf('Grid fraction with CIR earlier than GammaOU = %.3f\n',fracTemporalCIREarlier);
fprintf('Grid fraction with Full earlier than Uniform = %.3f\n',fracSpatialFullEarlier);
fprintf('Grid fraction with Joint-required > 5%%       = %.3f\n',fracJointPresent);
fprintf('========================================================\n');

%% ================================ PLOT ===================================

fig = figure('Color','w','Units','centimeters','Position',[1 1 18.3 14.0]);
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

xIdx = 1:nEta;
yIdx = 1:nM;
xRef = iEtaref;
yRef = iMref;

% ---------------- A: temporal contrast -----------------------------------
ax = nexttile(tl,1);
imagesc(ax,xIdx,yIdx,DeltaRMST_T);
set(ax,'YDir','normal');
limA = max(abs(DeltaRMST_T(:))); if limA<1e-12, limA=1; end
caxis(ax,[-limA limA]);
colormap(ax,divergingMap(256));
cb=colorbar(ax); cb.Label.String='\DeltaRMST_T (min)'; cb.FontSize=7;
hold(ax,'on');
plot(ax,xRef,yRef,'p','MarkerSize',10,'MarkerFaceColor','w', ...
    'MarkerEdgeColor','k','LineWidth',1.2);
title(ax,'Temporal-path ordering','FontWeight','normal');
xlabel(ax,'Relative source weight, \eta=\omega_S/\omega_T');
ylabel(ax,'Ferroptotic boundary, M_L');
applyGridTicks(ax,etaVec,Mvec); styleAxes(ax);
text(ax,0.02,1.03,'>0: CIR commits earlier','Units','normalized', ...
    'FontSize',6.8,'Color',[0.25 0.25 0.25]);

% ---------------- B: spatial contrast ------------------------------------
ax = nexttile(tl,2);
imagesc(ax,xIdx,yIdx,DeltaRMST_S);
set(ax,'YDir','normal');
limB = max(abs(DeltaRMST_S(:))); if limB<1e-12, limB=1; end
caxis(ax,[-limB limB]);
colormap(ax,divergingMap(256));
cb=colorbar(ax); cb.Label.String='\DeltaRMST_S (min)'; cb.FontSize=7;
hold(ax,'on');
plot(ax,xRef,yRef,'p','MarkerSize',10,'MarkerFaceColor','w', ...
    'MarkerEdgeColor','k','LineWidth',1.2);
title(ax,'Spatial-structure effect','FontWeight','normal');
xlabel(ax,'Relative source weight, \eta=\omega_S/\omega_T');
ylabel(ax,'Ferroptotic boundary, M_L');
applyGridTicks(ax,etaVec,Mvec); styleAxes(ax);
text(ax,0.02,1.03,'>0: Full commits earlier','Units','normalized', ...
    'FontSize',6.8,'Color',[0.25 0.25 0.25]);

% ---------------- C: joint-required fraction -----------------------------
ax = nexttile(tl,3);
imagesc(ax,xIdx,yIdx,JointFull);
set(ax,'YDir','normal'); caxis(ax,[0 1]);
colormap(ax,sequentialMap(256));
cb=colorbar(ax); cb.Label.String='Joint-required fraction'; cb.FontSize=7;
hold(ax,'on');
plot(ax,xRef,yRef,'p','MarkerSize',10,'MarkerFaceColor','w', ...
    'MarkerEdgeColor','k','LineWidth',1.2);
title(ax,'Cooperative commitment | Full background','FontWeight','normal');
xlabel(ax,'Relative source weight, \eta=\omega_S/\omega_T');
ylabel(ax,'Ferroptotic boundary, M_L');
applyGridTicks(ax,etaVec,Mvec); styleAxes(ax);

% ---------------- D: resistant fraction ----------------------------------
ax = nexttile(tl,4);
imagesc(ax,xIdx,yIdx,ResFull);
set(ax,'YDir','normal'); caxis(ax,[0 1]);
colormap(ax,sequentialMap(256));
cb=colorbar(ax); cb.Label.String='Resistant fraction'; cb.FontSize=7;
hold(ax,'on');
plot(ax,xRef,yRef,'p','MarkerSize',10,'MarkerFaceColor','w', ...
    'MarkerEdgeColor','k','LineWidth',1.2);
title(ax,'Non-committed population | Full background','FontWeight','normal');
xlabel(ax,'Relative source weight, \eta=\omega_S/\omega_T');
ylabel(ax,'Ferroptotic boundary, M_L');
applyGridTicks(ax,etaVec,Mvec); styleAxes(ax);

% panel labels
annotation(fig,'textbox',[0.015 0.955 0.03 0.03],'String','A', ...
    'EdgeColor','none','FontName','Arial','FontSize',11,'FontWeight','bold');
annotation(fig,'textbox',[0.505 0.955 0.03 0.03],'String','B', ...
    'EdgeColor','none','FontName','Arial','FontSize',11,'FontWeight','bold');
annotation(fig,'textbox',[0.015 0.485 0.03 0.03],'String','C', ...
    'EdgeColor','none','FontName','Arial','FontSize',11,'FontWeight','bold');
annotation(fig,'textbox',[0.505 0.485 0.03 0.03],'String','D', ...
    'EdgeColor','none','FontName','Arial','FontSize',11,'FontWeight','bold');

%% =============================== SAVE ====================================

stem = fullfile(outDir,'SI_Meso_ML_eta_robustness_FINAL');
exportgraphics(fig,[stem '.png'],'Resolution',EXPORT_DPI);
exportgraphics(fig,[stem '.tif'],'Resolution',EXPORT_DPI);
exportgraphics(fig,[stem '.pdf'],'ContentType','vector');
savefig(fig,[stem '.fig']);

%% ------------------------------ TABLE ------------------------------------

rows = {};
for iM = 1:nM
    for ie = 1:nEta
        eta = etaVec(ie);
        omegaT = Omega/(1+eta);
        omegaS = eta*Omega/(1+eta);

        rows(end+1,:) = { ...
            Mvec(iM),eta,omegaT,omegaS, ...
            DeltaRMST_T(iM,ie),DeltaRMST_S(iM,ie), ...
            JointFull(iM,ie),ResFull(iM,ie), ...
            Pcommit(iM,ie,1),Pcommit(iM,ie,2), ...
            Pcommit(iM,ie,3),Pcommit(iM,ie,4), ...
            RMST(iM,ie,1),RMST(iM,ie,2), ...
            RMST(iM,ie,3),RMST(iM,ie,4)};
    end
end

summaryTable = cell2table(rows,'VariableNames',{ ...
    'M_L','eta','omegaT','omegaS', ...
    'DeltaRMST_T_GminusC','DeltaRMST_S_UminusF', ...
    'JointFull','ResistantFull', ...
    'Pcommit_CU','Pcommit_GU','Pcommit_CF','Pcommit_GF', ...
    'RMST_CU','RMST_GU','RMST_CF','RMST_GF'});

writetable(summaryTable,fullfile(outDir, ...
    'SI_Meso_ML_eta_robustness_FINAL_summary.csv'));

save(fullfile(outDir,'SI_Meso_ML_eta_robustness_FINAL_data.mat'), ...
    'Mvec','etaVec','Omega','omegaT0','omegaS0','rho0','L0','ML0', ...
    'Pcommit','RMST','routeFrac','DeltaRMST_T','DeltaRMST_S', ...
    'JointFull','ResFull','summaryTable', ...
    'fracTemporalCIREarlier','fracSpatialFullEarlier','fracJointPresent', ...
    '-v7.3');

fprintf('\nCompleted successfully.\n');
fprintf('Outputs saved in:\n%s\n',outDir);

%% =========================================================================
%% LOCAL FUNCTIONS
%% =========================================================================

function R = evaluateRoutesExact(IT,IS,omegaT,omegaS,rho,ML,L0,dt,Tobs)

    assert(isequal(size(IT),size(IS)), ...
        'IT and IS must have identical paired-cell dimensions.');

    UT  = omegaT*IT;
    US  = omegaS*IS;
    UTS = UT+US;

    FPT_T  = propagateLPOExact(UT, rho,ML,L0,dt,Tobs);
    FPT_S  = propagateLPOExact(US, rho,ML,L0,dt,Tobs);
    FPT_TS = propagateLPOExact(UTS,rho,ML,L0,dt,Tobs);

    crossT  = FPT_T <= Tobs;
    crossS  = FPT_S <= Tobs;
    crossTS = FPT_TS<= Tobs;

    if any((crossT | crossS) & ~crossTS)
        error('Pathwise monotonicity violation in route classification.');
    end

    routeCode = zeros(size(FPT_TS));
    routeCode(crossT & ~crossS) = 1;                 % T-sufficient
    routeCode(~crossT & crossS) = 2;                 % S-sufficient
    routeCode(crossT & crossS) = 3;                  % Redundant
    routeCode(~crossT & ~crossS & crossTS) = 4;      % Joint-required
    routeCode(~crossTS) = 5;                         % Resistant

    if any(routeCode==0)
        error('Unclassified fate route detected.');
    end

    R = struct();
    R.FPT_T  = FPT_T;
    R.FPT_S  = FPT_S;
    R.FPT_TS = FPT_TS;
    R.crossT  = crossT;
    R.crossS  = crossS;
    R.crossTS = crossTS;
    R.routeCode = routeCode;
    R.routeFrac = arrayfun(@(rr) mean(routeCode==rr),1:5);
    R.Pcommit = mean(crossTS);
    R.RMST = mean(min(FPT_TS,Tobs));
end

function FPT = propagateLPOExact(U,rho,ML,L0,dt,Tobs)

    [n,nBins] = size(U);
    L = L0*ones(n,1);
    FPT = inf(n,1);

    for b = 1:nBins
        Ub = max(U(:,b),0);
        a = Ub + rho;

        Lnext = L;
        active = a > eps;

        if any(active)
            aa = a(active);
            uu = Ub(active);
            Lss = uu./aa;
            E = exp(-aa*dt);
            Lnext(active) = Lss + (L(active)-Lss).*E;
        end

        Lnext = min(max(Lnext,0),1);

        crossing = isinf(FPT) & (L < ML) & (Lnext >= ML);

        if any(crossing)
            idx = find(crossing);
            aa = a(idx);
            uu = Ub(idx);
            Lss = uu./aa;
            ratio = (ML-Lss)./(L(idx)-Lss);
            ratio = max(min(ratio,1),realmin);
            tauCross = -log(ratio)./aa;
            tauCross = max(min(tauCross,dt),0);
            FPT(idx) = (b-1)*dt + tauCross;
        end

        L = Lnext;
    end

    FPT(FPT>Tobs) = inf;
end

function v = getStructScalar(S,names,defaultValue)
    v = defaultValue;
    for k = 1:numel(names)
        if isfield(S,names{k})
            x = S.(names{k});
            if isnumeric(x) && isscalar(x)
                v = x;
                return;
            end
        end
    end
end

function applyGridTicks(ax,etaVec,Mvec)
    set(ax, ...
        'XTick',1:numel(etaVec), ...
        'YTick',1:numel(Mvec), ...
        'XTickLabel',compose('%g',etaVec), ...
        'YTickLabel',compose('%.2f',Mvec));
    xlim(ax,[0.5 numel(etaVec)+0.5]);
    ylim(ax,[0.5 numel(Mvec)+0.5]);
end

function styleAxes(ax)
    set(ax, ...
        'FontName','Arial', ...
        'FontSize',8, ...
        'LineWidth',0.75, ...
        'TickDir','out', ...
        'TickLength',[0.014 0.014], ...
        'Box','off', ...
        'Layer','top');
end

function cmap = divergingMap(n)
    if nargin<1, n=256; end
    blue  = [0.20 0.43 0.69];
    white = [0.97 0.97 0.96];
    red   = [0.80 0.32 0.24];
    n1=floor(n/2); n2=n-n1;
    left=zeros(n1,3); right=zeros(n2,3);
    for k=1:3
        left(:,k)=linspace(blue(k),white(k),n1);
        right(:,k)=linspace(white(k),red(k),n2);
    end
    cmap=[left;right];
end

function cmap = sequentialMap(n)
    if nargin<1, n=256; end
    anchors = [ ...
        0.95 0.96 0.97; ...
        0.74 0.84 0.89; ...
        0.43 0.68 0.76; ...
        0.27 0.52 0.64; ...
        0.16 0.34 0.50];
    x=linspace(0,1,size(anchors,1));
    xi=linspace(0,1,n);
    cmap=zeros(n,3);
    for k=1:3
        cmap(:,k)=interp1(x,anchors(:,k),xi,'pchip');
    end
    cmap=min(max(cmap,0),1);
end
