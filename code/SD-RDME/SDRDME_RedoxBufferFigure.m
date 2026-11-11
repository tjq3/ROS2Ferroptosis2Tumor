%% SDRDME_RedoxBufferFigure.m
% Publication-style post-processing for the FINAL local redox-buffer scan.
%
% This script DOES NOT rerun the SD-RDME simulation. It reads the FINAL
% summary tables produced by SDRDME_RedoxBufferScan.m and creates
% a mechanism-led four-panel figure:
%   A, local Fe2+ source substitution;
%   B, rebuilding of the local Fe2+ reserve;
%   C, buffering of the transport-induced burden loss;
%   D, release of the high-threshold first-passage-time delay.
%
% Required input folder:
%   condition_redox_buffer_results/final
%
% Main output folder:
%   condition_redox_buffer_results/final/highlight

clear; close all; clc;

%% User-adjustable presentation settings
PRIMARY_GAMMA = 1;       % Finalized physical diffusion field
STRESS_GAMMA = 2;        % Stronger heterogeneity used as a stress test
HIGH_THRESHOLD = 1.4;    % Pre-specified high-damage threshold multiplier
EXPORT_DPI = 600;

scriptDir = fileparts(mfilename('fullpath'));
resultDir = fullfile(scriptDir,'condition_redox_buffer_results','final');
inputCSV = fullfile(resultDir,'Fig3H_redoxBuffer_results.csv');
inputThresholdCSV = fullfile(resultDir, ...
    'Fig3H_redoxBuffer_thresholdSensitivity.csv');
outDir = fullfile(resultDir,'highlight');
if ~exist(outDir,'dir'), mkdir(outDir); end

assert(exist(inputCSV,'file')==2, ...
    'FINAL summary CSV is missing: %s',inputCSV);
assert(exist(inputThresholdCSV,'file')==2, ...
    'FINAL threshold-sensitivity CSV is missing: %s',inputThresholdCSV);

T = readtable(inputCSV);
TT = readtable(inputThresholdCSV);

requiredMain = { ...
    'RecyclingCapacityScale','DiffusionHeterogeneityGamma', ...
    'RecyclingFraction_R2overR1','DiffusiveSupplyFraction_overR1', ...
    'Fe2FinalMean','LocalBurdenRatio','LocalBurdenRatio_CI_L', ...
    'LocalBurdenRatio_CI_U'};
requiredThreshold = { ...
    'RecyclingCapacityScale','DiffusionHeterogeneityGamma', ...
    'ThresholdMultiplier','PiTransport_RMSTDifference', ...
    'PiTransport_CI_L','PiTransport_CI_U'};
assertColumns(T,requiredMain,'FINAL main table');
assertColumns(TT,requiredThreshold,'FINAL threshold table');

rho = sort(unique(T.RecyclingCapacityScale));
gamma = sort(unique(T.DiffusionHeterogeneityGamma));
idxGamma0 = exactIndex(gamma,0,'gamma=0 reference');
idxGamma1 = exactIndex(gamma,PRIMARY_GAMMA,'primary gamma');
idxGamma2 = exactIndex(gamma,STRESS_GAMMA,'stress-test gamma');

assert(numel(rho)>=4, ...
    'At least four recycling-capacity levels are required.');
assert(all(rho>0),'Recycling-capacity values must be positive.');

recycle = matrixFromTable(T,'RecyclingFraction_R2overR1',rho,gamma);
diffusive = matrixFromTable(T,'DiffusiveSupplyFraction_overR1',rho,gamma);
fe2 = matrixFromTable(T,'Fe2FinalMean',rho,gamma);
burdenRatio = matrixFromTable(T,'LocalBurdenRatio',rho,gamma);
burdenL = matrixFromTable(T,'LocalBurdenRatio_CI_L',rho,gamma);
burdenU = matrixFromTable(T,'LocalBurdenRatio_CI_U',rho,gamma);

thresholdMask = abs(TT.ThresholdMultiplier-HIGH_THRESHOLD)<1e-10;
THigh = TT(thresholdMask,:);
assert(~isempty(THigh), ...
    'Threshold multiplier %.3g is absent from the FINAL table.', ...
    HIGH_THRESHOLD);
piHigh = matrixFromTable(THigh,'PiTransport_RMSTDifference',rho,gamma);
piHighL = matrixFromTable(THigh,'PiTransport_CI_L',rho,gamma);
piHighU = matrixFromTable(THigh,'PiTransport_CI_U',rho,gamma);

assert(all(isfinite(recycle(:))) && all(isfinite(diffusive(:))) && ...
       all(isfinite(fe2(:))) && all(isfinite(burdenRatio(:))), ...
    'The FINAL main table contains missing/non-finite grid cells.');
assert(all(isfinite(piHigh(:))) && all(isfinite(piHighL(:))) && ...
       all(isfinite(piHighU(:))), ...
    'The selected FINAL threshold table contains missing/non-finite cells.');
assert(all(burdenL(:)<=burdenRatio(:)+1e-12) && ...
       all(burdenU(:)>=burdenRatio(:)-1e-12), ...
    'Invalid burden confidence interval detected.');
assert(all(piHighL(:)<=piHigh(:)+1e-12) && ...
       all(piHighU(:)>=piHigh(:)-1e-12), ...
    'Invalid first-passage-time confidence interval detected.');

% Data-defined landmarks used for annotation rather than pilot hard-coding.
sourceSwitchIdx = find(recycle(:,idxGamma1)>=1,1,'first');
if isempty(sourceSwitchIdx)
    sourceSwitchRho = NaN;
else
    sourceSwitchRho = rho(sourceSwitchIdx);
end

% "Buffer onset" is the first FINAL recycling level at which the 95% CI of
% the physical-field FPT delay includes zero.
bufferOnsetIdx = find(piHighL(:,idxGamma1)<=0 & ...
                     piHighU(:,idxGamma1)>=0,1,'first');
if isempty(bufferOnsetIdx)
    bufferOnsetRho = NaN;
else
    bufferOnsetRho = rho(bufferOnsetIdx);
end

%% Save the exact values behind the highlight figure
summaryTable = table( ...
    rho, ...
    recycle(:,idxGamma1),diffusive(:,idxGamma1), ...
    recycle(:,idxGamma2),diffusive(:,idxGamma2), ...
    fe2(:,idxGamma0),fe2(:,idxGamma1),fe2(:,idxGamma2), ...
    burdenRatio(:,idxGamma1),burdenL(:,idxGamma1), ...
    burdenU(:,idxGamma1), ...
    burdenRatio(:,idxGamma2),burdenL(:,idxGamma2), ...
    burdenU(:,idxGamma2), ...
    piHigh(:,idxGamma1),piHighL(:,idxGamma1), ...
    piHighU(:,idxGamma1), ...
    piHigh(:,idxGamma2),piHighL(:,idxGamma2), ...
    piHighU(:,idxGamma2), ...
    'VariableNames',{ ...
    'RecyclingCapacityScale', ...
    'R2overR1_PrimaryGamma','DiffusionOverR1_PrimaryGamma', ...
    'R2overR1_StressGamma','DiffusionOverR1_StressGamma', ...
    'Fe2_Gamma0','Fe2_PrimaryGamma','Fe2_StressGamma', ...
    'BurdenRatio_PrimaryGamma','BurdenRatio_PrimaryGamma_CI_L', ...
    'BurdenRatio_PrimaryGamma_CI_U', ...
    'BurdenRatio_StressGamma','BurdenRatio_StressGamma_CI_L', ...
    'BurdenRatio_StressGamma_CI_U', ...
    'HighThresholdPi_PrimaryGamma','HighThresholdPi_PrimaryGamma_CI_L', ...
    'HighThresholdPi_PrimaryGamma_CI_U', ...
    'HighThresholdPi_StressGamma','HighThresholdPi_StressGamma_CI_L', ...
    'HighThresholdPi_StressGamma_CI_U'});
writetable(summaryTable,fullfile(outDir, ...
    'Fig3H_redoxBuffer_highlight_FINAL_summary.csv'));

%% Figure style
oldDefaults = captureGraphicsDefaults();
cleanupDefaults = onCleanup(@()restoreGraphicsDefaults(oldDefaults)); %#ok<NASGU>

set(groot,'defaultFigureColor','w');
set(groot,'defaultAxesFontName','Arial');
set(groot,'defaultTextFontName','Arial');
set(groot,'defaultAxesFontSize',6.9);
set(groot,'defaultTextFontSize',6.9);
set(groot,'defaultAxesLineWidth',0.75);
set(groot,'defaultAxesTickDir','out');

% Two-root palette with line style and marker redundancy for grayscale use.
purple = [0.40 0.20 0.58];
orange = [0.86 0.36 0.10];
blue = [0.12 0.42 0.68];
charcoal = [0.24 0.25 0.27];
lightGray = [0.86 0.87 0.88];
midGray = [0.48 0.49 0.51];

fig = figure('Units','centimeters','Position',[1 1 18.8 11.4], ...
    'Color','w','Renderer','painters');
tl = tiledlayout(fig,2,2,'TileSpacing','loose','Padding','loose');

%% A. Local Fe2+ source substitution in the finalized physical field
ax = nexttile(tl,1); hold(ax,'on');
hRecycle = semilogx(ax,rho,recycle(:,idxGamma1),'-o', ...
    'Color',purple,'MarkerFaceColor','w','MarkerEdgeColor',purple, ...
    'MarkerSize',4.2,'LineWidth',1.2);
hDiff = semilogx(ax,rho,diffusive(:,idxGamma1),'--s', ...
    'Color',blue,'MarkerFaceColor','w','MarkerEdgeColor',blue, ...
    'MarkerSize',4.1,'LineWidth',1.15);
yline(ax,1,':','Color',midGray,'LineWidth',0.8, ...
    'HandleVisibility','off');
if isfinite(sourceSwitchRho)
    xline(ax,sourceSwitchRho,':','Color',charcoal,'LineWidth',0.8, ...
        'HandleVisibility','off');
    yText = max([recycle(:,idxGamma1); 1])*0.91;
    text(ax,sourceSwitchRho*1.07,yText,'source switch', ...
        'Color',charcoal,'FontSize',6.0, ...
        'VerticalAlignment','top','HorizontalAlignment','left');
end
xlabel(ax,'Recycling capacity, c_2/c_2^0');
ylabel(ax,'Local Fe^{2+} source / R1');
title(ax,'Local Fe^{2+} source substitution (\gamma_D=1)', ...
    'FontSize',8.4,'FontWeight','normal');
legend(ax,[hRecycle hDiff], ...
    {'local regeneration, R2/R1','net diffusion / R1'}, ...
    'Location','east','Box','off','FontSize',6.2);
setLogTicks(ax,rho); panelLabel(ax,'A'); formatAxes(ax);

%% B. The source switch rebuilds the local Fe2+ reserve
ax = nexttile(tl,2); hold(ax,'on');
h0 = semilogx(ax,rho,fe2(:,idxGamma0),':^', ...
    'Color',midGray,'MarkerFaceColor','w','MarkerEdgeColor',midGray, ...
    'MarkerSize',4.0,'LineWidth',1.0);
h1 = semilogx(ax,rho,fe2(:,idxGamma1),'-o', ...
    'Color',purple,'MarkerFaceColor','w','MarkerEdgeColor',purple, ...
    'MarkerSize',4.2,'LineWidth',1.2);
h2 = semilogx(ax,rho,fe2(:,idxGamma2),'--s', ...
    'Color',orange,'MarkerFaceColor','w','MarkerEdgeColor',orange, ...
    'MarkerSize',4.1,'LineWidth',1.15);
xlabel(ax,'Recycling capacity, c_2/c_2^0');
ylabel(ax,'Final local Fe^{2+} (per voxel)');
title(ax,'Recycling rebuilds the hotspot Fe^{2+} reserve', ...
    'FontSize',8.4,'FontWeight','normal');
legend(ax,[h0 h1 h2], ...
    {'\gamma_D=0 reference','\gamma_D=1 physical field', ...
     '\gamma_D=2 stress test'}, ...
    'Location','southeast','Box','off','FontSize',6.0);
setLogTicks(ax,rho); panelLabel(ax,'B'); formatAxes(ax);

%% C. Burden penalty is buffered as recycling increases
ax = nexttile(tl,3); hold(ax,'on');
h1 = plotErrorSeries(ax,rho,burdenRatio(:,idxGamma1), ...
    burdenL(:,idxGamma1),burdenU(:,idxGamma1), ...
    purple,'-','o');
h2 = plotErrorSeries(ax,rho,burdenRatio(:,idxGamma2), ...
    burdenL(:,idxGamma2),burdenU(:,idxGamma2), ...
    orange,'--','s');
yline(ax,1,':','Color',midGray,'LineWidth',0.8, ...
    'HandleVisibility','off');
xlabel(ax,'Recycling capacity, c_2/c_2^0');
ylabel(ax,'Burden ratio, B_A(\gamma_D)/B_A(0)');
title(ax,'Recycling buffers transport-induced burden loss', ...
    'FontSize',8.4,'FontWeight','normal');
legend(ax,[h1 h2], ...
    {'\gamma_D=1 physical field','\gamma_D=2 stress test'}, ...
    'Location','southeast','Box','off','FontSize',6.1);
setLogTicks(ax,rho); panelLabel(ax,'C'); formatAxes(ax);

%% D. High-threshold FPT reveals the statistically defined buffer onset
ax = nexttile(tl,4); hold(ax,'on');
h1 = plotErrorSeries(ax,rho,piHigh(:,idxGamma1), ...
    piHighL(:,idxGamma1),piHighU(:,idxGamma1), ...
    purple,'-','o');
h2 = plotErrorSeries(ax,rho,piHigh(:,idxGamma2), ...
    piHighL(:,idxGamma2),piHighU(:,idxGamma2), ...
    orange,'--','s');
yline(ax,0,':','Color',midGray,'LineWidth',0.8, ...
    'HandleVisibility','off');
if isfinite(bufferOnsetRho)
    xline(ax,bufferOnsetRho,'-.','Color',charcoal,'LineWidth',0.9, ...
        'HandleVisibility','off');
    yl = paddedLimits([piHighL(:,idxGamma1);piHighU(:,idxGamma1); ...
                       piHighL(:,idxGamma2);piHighU(:,idxGamma2)],0.08);
    ylim(ax,yl);
    text(ax,bufferOnsetRho*1.07,yl(2)-0.35*diff(yl), ...
        sprintf('95%% CI includes 0 at %.2g',bufferOnsetRho), ...
        'Color',charcoal,'FontSize',6.0, ...
        'VerticalAlignment','top','HorizontalAlignment','left');
end
xlabel(ax,'Recycling capacity, c_2/c_2^0');
ylabel(ax,sprintf('FPT delay at %.1fx threshold (min)', ...
    HIGH_THRESHOLD));
title(ax,'Local recycling releases the high-threshold delay', ...
    'FontSize',8.4,'FontWeight','normal');
legend(ax,[h1 h2], ...
    {'\gamma_D=1 physical field','\gamma_D=2 stress test'}, ...
    'Location','northeast','Box','off','FontSize',6.1);
setLogTicks(ax,rho); panelLabel(ax,'D'); formatAxes(ax);

sgtitle(tl,'Local iron recycling buffers diffusion-limited Fenton damage', ...
    'FontName','Arial','FontSize',9.5,'FontWeight','bold');

%% Export
outBase = fullfile(outDir,'Fig3H_redoxBuffer_highlight_FINAL');
exportgraphics(fig,[outBase '.pdf'],'ContentType','vector');
exportgraphics(fig,[outBase '.png'],'Resolution',EXPORT_DPI);
exportgraphics(fig,[outBase '.tif'],'Resolution',EXPORT_DPI);
savefig(fig,[outBase '.fig']);

fprintf('\nFINAL redox-buffer highlight figure generated.\n');
fprintf('  Input:  %s\n',resultDir);
fprintf('  Output: %s.[pdf/png/tif/fig]\n',outBase);
fprintf('  High-damage threshold: %.2fx\n',HIGH_THRESHOLD);
if isfinite(sourceSwitchRho)
    fprintf('  Physical-field source switch (R2/R1 >= 1): rho = %.3g\n', ...
        sourceSwitchRho);
end
if isfinite(bufferOnsetRho)
    fprintf(['  Physical-field buffer onset (FPT-delay 95%% CI first ' ...
        'includes zero): rho = %.3g\n'],bufferOnsetRho);
else
    fprintf(['  Physical-field buffer onset was not reached within the ' ...
        'FINAL scan.\n']);
end

%% Local functions
function assertColumns(T,names,tableLabel)
    missing = setdiff(names,T.Properties.VariableNames);
    assert(isempty(missing),'%s is missing columns: %s', ...
        tableLabel,strjoin(missing,', '));
end

function idx = exactIndex(values,target,label)
    idx = find(abs(values-target)<1e-10,1,'first');
    assert(~isempty(idx),'Required %s (%.3g) is absent.',label,target);
end

function M = matrixFromTable(T,varName,rowValues,colValues)
    M = nan(numel(rowValues),numel(colValues));
    hit = zeros(size(M));
    for i = 1:height(T)
        ir = find(abs(rowValues-T.RecyclingCapacityScale(i))<1e-10, ...
            1,'first');
        ic = find(abs(colValues-T.DiffusionHeterogeneityGamma(i))<1e-10, ...
            1,'first');
        assert(~isempty(ir) && ~isempty(ic), ...
            'Unrecognized grid coordinate in column %s.',varName);
        M(ir,ic) = T.(varName)(i);
        hit(ir,ic) = hit(ir,ic)+1;
    end
    assert(all(hit(:)==1), ...
        'Column %s does not contain exactly one row per grid cell.',varName);
end

function h = plotErrorSeries(ax,x,y,lo,hi,color,lineStyle,marker)
    lower = max(0,y-lo);
    upper = max(0,hi-y);
    h = errorbar(ax,x,y,lower,upper, ...
        'Color',color,'LineStyle',lineStyle,'Marker',marker, ...
        'MarkerFaceColor','w','MarkerEdgeColor',color, ...
        'MarkerSize',4.2,'LineWidth',1.15,'CapSize',2.5);
    set(ax,'XScale','log');
end

function setLogTicks(ax,rho)
    preferred = [0.25 0.5 1 2 4 8];
    tickValues = [];
    for i = 1:numel(preferred)
        if any(abs(rho-preferred(i))<1e-10)
            tickValues(end+1) = preferred(i); %#ok<AGROW>
        end
    end
    if numel(tickValues)<4
        tickValues = rho;
    end
    set(ax,'XScale','log','XTick',tickValues, ...
        'XTickLabel',compose('%g',tickValues));
    xlim(ax,[min(rho)/1.16 max(rho)*1.16]);
end

function panelLabel(ax,label)
    text(ax,0.018,0.97,label,'Units','normalized', ...
        'FontWeight','bold','FontSize',11.5, ...
        'VerticalAlignment','top','HorizontalAlignment','left');
end

function formatAxes(ax)
    set(ax,'Box','on','LineWidth',0.75,'TickDir','out','Layer','top', ...
        'FontSize',6.9,'XColor',[0.08 0.08 0.08], ...
        'YColor',[0.08 0.08 0.08]);
    ax.XMinorTick = 'off';
    ax.YMinorTick = 'off';
end

function lim = paddedLimits(values,fraction)
    values = values(isfinite(values));
    lo = min(values); hi = max(values);
    span = max(hi-lo,eps);
    lim = [lo-fraction*span,hi+fraction*span];
end

function defaults = captureGraphicsDefaults()
    defaults.figureColor = get(groot,'defaultFigureColor');
    defaults.axesFontName = get(groot,'defaultAxesFontName');
    defaults.textFontName = get(groot,'defaultTextFontName');
    defaults.axesFontSize = get(groot,'defaultAxesFontSize');
    defaults.textFontSize = get(groot,'defaultTextFontSize');
    defaults.axesLineWidth = get(groot,'defaultAxesLineWidth');
    defaults.axesTickDir = get(groot,'defaultAxesTickDir');
end

function restoreGraphicsDefaults(defaults)
    set(groot,'defaultFigureColor',defaults.figureColor);
    set(groot,'defaultAxesFontName',defaults.axesFontName);
    set(groot,'defaultTextFontName',defaults.textFontName);
    set(groot,'defaultAxesFontSize',defaults.axesFontSize);
    set(groot,'defaultTextFontSize',defaults.textFontSize);
    set(groot,'defaultAxesLineWidth',defaults.axesLineWidth);
    set(groot,'defaultAxesTickDir',defaults.axesTickDir);
end
