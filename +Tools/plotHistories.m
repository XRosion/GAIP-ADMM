function plotHistories(results, problem, figureTitle, savePrefix, options)
%PLOTHISTORIES Plot objective and feasibility histories for algorithms.
%   PLOTHISTORIES(RESULTS, PROBLEM, FIGURETITLE, SAVEPREFIX, OPTS) creates
%   publication-style semilog iteration figures for a collection of algorithm
%   runs. For TV experiments the function also saves PSNR curves; when CPU
%   time is recorded it creates the corresponding time-based comparisons.
%
%   Input:
%       results     - cell array returned by run_algorithm_set. Each entry
%                     must contain name and history fields.
%       problem     - problem structure. The field PhiStar is used in the
%                     axis label and should come from the problem's
%                     reference solve.
%       figureTitle - title prefix for generated figures.
%       savePrefix  - path prefix such as fullfile(rootDir,'Results','LAD').
%                     Empty input saves under the project Results folder.
%       options     - plotting and experiment options. The function accepts
%                     options.tolerance. Optional visualization fields include
%                     plotSavePdf, plotLegendLocation, plotMaxLinePoints,
%                     plotPngDpi, plotFigureSizeInches, plotUseMarkers,
%                     plotShowTitle, plotObjectiveEnhanced,
%                     plotObjectiveSmoothingWindow, plotObjectiveRawMarkers,
%                     plotObjectiveRawMarkerMaxPoints, plotObjectiveMinPoint,
%                     plotObjectiveAxesLeftMargin, plotObjectiveSuffix, and
%                     plotFeasibilitySuffix. plotShowToleranceLine controls
%                     the existing feasibility tolerance reference line.
%
%   Output:
%       No MATLAB output. The function saves objective and feasibility PNG
%       files and optional PDF files.
if nargin < 1
    error('plotHistories:NotEnoughInputs', 'At least the results input is required.');
end
if nargin < 2
    problem = struct();
end
if nargin < 3 || isempty(figureTitle)
    figureTitle = 'Algorithm histories';
end
if nargin < 4 || isempty(savePrefix)
    savePrefix = defaultSavePrefix('histories');
end
if nargin < 5
    options = struct();
end

if isLegacyNameList(problem)
    results = normalizeResults(results, problem);
    problem = struct();
elseif iscell(results)
    results = normalizeResults(results, {});
else
    error('plotHistories:InvalidResults', 'results must be a cell array of algorithm outputs.');
end

plotConfiguration = plotConfig(options);
names = cellfun(@(r) r.name, results, 'UniformOutput', false);

% --- 修改点 1：将 '_{spp}' 替换为 LaTeX 数学模式 '$_{spp}$' ---
names = cellfun(@(n) strrep(n, '_{spp}', '$_{\rm spp}$'), names, 'UniformOutput', false);

legendNames = cellfun(@latexText, names, 'UniformOutput', false);
styles = plotStyles(numel(results), plotConfiguration);
if isempty(savePrefix)
    plotConfiguration.visible = 'on';
    plotConfiguration.closeFigures = false;
else
    plotConfiguration.visible = 'off';
    plotConfiguration.closeFigures = true;
end

objectiveCurves = cell(numel(results), 1);
objectiveSuffix = getOption(options, 'plotObjectiveSuffix', 'objective');
objectiveTimeSuffix = 'objective_time';
objectiveTitle = [figureTitle ' relative objective error'];
objectiveTimeTitle = [figureTitle ' relative objective error CPU time'];
objectiveXLabel = 'Iteration';
objectiveYLabel = objectiveLabel(problem);
for i = 1:numel(results)
    objectiveCurves{i} = prepareObjectiveCurve(results{i}, problem);
end
plotCurveFigure([], objectiveCurves, legendNames, styles, ...
    objectiveTitle, objectiveXLabel, ...
    objectiveYLabel, ...
    savePrefix, objectiveSuffix, plotConfiguration, false, 'objective');

feasCurves = cell(numel(results), 1);
for i = 1:numel(results)
    feasCurves{i} = prepareLogCurve(results{i}.history.feas);
end
plotCurveFigure([], feasCurves, legendNames, styles, ...
    [figureTitle ' feasibility'], 'Iteration', ...
    feasibilityLabel(problem), ...
    savePrefix, getOption(options, 'plotFeasibilitySuffix', 'feasibility_iter'), ...
    plotConfiguration, true, 'default');

timeCurves = prepareTimeCurves(results);
if ~isempty(timeCurves)
    plotCurveFigure(timeCurves, objectiveCurves, legendNames, styles, ...
        objectiveTimeTitle, 'CPU time (s)', ...
        objectiveYLabel, ...
        savePrefix, objectiveTimeSuffix, plotConfiguration, false, 'objective');
    plotCurveFigure(timeCurves, feasCurves, legendNames, styles, ...
        [figureTitle ' feasibility CPU time'], 'CPU time (s)', ...
        feasibilityLabel(problem), ...
        savePrefix, 'feasibility_time', plotConfiguration, true, 'default');
end

if allHaveHistoryField(results, 'psnr')
    psnrCurves = cell(numel(results), 1);
    for i = 1:numel(results)
        psnrCurves{i} = results{i}.history.psnr(:);
    end
    plotLinearCurveFigure([], psnrCurves, legendNames, styles, ...
        [figureTitle ' PSNR'], 'Iteration', 'PSNR (dB)', ...
        savePrefix, 'psnr_iter', plotConfiguration, 'default');
    if ~isempty(timeCurves)
        plotLinearCurveFigure(timeCurves, psnrCurves, legendNames, styles, ...
            [figureTitle ' PSNR CPU time'], 'CPU time (s)', 'PSNR (dB)', ...
            savePrefix, 'psnr_time', plotConfiguration, 'default');
    end
end

end

function tf = isLegacyNameList(x)
%ISLEGACYNAMELIST True for plotHistories(results, {'A','B',...}).
tf = iscell(x) && all(cellfun(@(v) ischar(v) || isStringScalar(v), x));
end

function tf = isStringScalar(x)
%ISSTRINGSCALAR Compatible scalar string check for older MATLAB releases.
tf = isa(x, 'string') && isscalar(x);
end

function results = normalizeResults(results, names)
%NORMALIZERESULTS Accept both .history and older .hist result structures.
if ~iscell(results)
    error('plotHistories:InvalidResults', 'results must be a cell array.');
end
if ~isempty(names) && numel(names) ~= numel(results)
    error('plotHistories:NameCountMismatch', ...
        'The number of names must match the number of result entries.');
end

for i = 1:numel(results)
    if ~isstruct(results{i})
        error('plotHistories:InvalidResultEntry', ...
            'Each result entry must be an algorithm output structure.');
    end

    if ~isfield(results{i}, 'history') && isfield(results{i}, 'hist')
        results{i}.history = results{i}.hist;
    end
    if ~isfield(results{i}, 'history')
        error('plotHistories:MissingHistory', ...
            'Result entry %d must contain a history or hist field.', i);
    end
    validateHistory(results{i}.history, i);

    if ~isempty(names)
        results{i}.name = char(names{i});
    elseif ~isfield(results{i}, 'name') || isempty(results{i}.name)
        results{i}.name = sprintf('Algorithm %d', i);
    end
end
end

function validateHistory(history, index)
%VALIDATEHISTORY Ensure required histories are present before plotting.
requiredFields = {'obj', 'feas'};
for j = 1:numel(requiredFields)
    fieldName = requiredFields{j};
    if ~isfield(history, fieldName) || isempty(history.(fieldName))
        error('plotHistories:MissingHistoryField', ...
            'Result entry %d must contain a nonempty history.%s field.', ...
            index, fieldName);
    end
end
end

function rel = prepareObjectiveCurve(result, problem)
%PREPARE_OBJECTIVE_CURVE Return optimality-gap curves.
if isfield(result.history, 'err') && ~isempty(result.history.err)
    rel = result.history.err(:);
    rel = prepareLogCurve(rel);
    return;
end

if isTVProblem(problem) && isfield(result.history, 'objPrimal') && ...
        ~isempty(result.history.objPrimal) && isfield(problem, 'PhiStar') && ...
        isfinite(problem.PhiStar)
    rel = abs(result.history.objPrimal(:) - problem.PhiStar) ./ ...
        max(abs(problem.PhiStar), eps);
    rel = prepareLogCurve(rel);
    return;
end

if isfield(problem, 'PhiStar') && isfinite(problem.PhiStar)
    rel = abs(result.history.obj - problem.PhiStar) ./ max(abs(problem.PhiStar), eps);
else
    rel = abs(result.history.obj);
end
rel = prepareLogCurve(rel);
end

function tf = isTVProblem(problem)
tf = isfield(problem, 'type') && strcmpi(problem.type, 'TV');
end

function label = objectiveLabel(problem)
if isTVProblem(problem)
    if isfield(problem, 'PhiStar') && isfinite(problem.PhiStar)
        label = '$\frac{|\Phi(x^{k+1},w^{k+1})-\Phi^*|}{|\Phi^*|}$';
    else
        label = '$|F(x^k)-F^*|$';
    end
    return;
end

if isfield(problem, 'PhiStar') && isfinite(problem.PhiStar)
    label = '$\frac{|\Phi(x^{k+1},w^{k+1})-\Phi^*|}{|\Phi^*|}$';
else
    label = '$|\Phi(x^{k+1},w^{k+1})|$';
end
end

function label = feasibilityLabel(problem)
%FEASIBILITYLABEL Return the split residual label for known problem types.
if isfield(problem, 'type')
    switch lower(problem.type)
        case 'lad'
            label = '$\|Ax^{k+1} +B w^{k+1} - b\|_2$';
        case {'lasso', 'elastic'}
            label = '$\|Ax^{k+1} + Bw^{k+1} - b\|_2$';
        case 'tv'
            label = '$\|Ax^{k+1} + B w^{k+1} - b\|_2$';
        otherwise
            label = '$\mathrm{feasibility\ residual}$';
    end
else
    label = '$\mathrm{feasibility\ residual}$';
end
end

function plotCurveFigure(xCurves, yCurves, legendNames, styles, plotTitle, ...
    xLabelText, yLabelText, savePrefix, suffix, plotConfiguration, showTol, curveKind)
%PLOT_CURVE_FIGURE Draw and save one semilog curve figure.
if nargin < 12 || isempty(curveKind)
    curveKind = 'default';
end
figureHandle = squareFigure(plotTitle, plotConfiguration);
plotAxes = axes('Parent', figureHandle, ...
    'Position', curveAxesPosition(plotConfiguration, curveKind));
hold(plotAxes, 'on');

allX = [];
allY = [];
useObjectiveTrace = strcmpi(char(curveKind), 'objective') && ...
    plotConfiguration.objectiveEnhanced;
for i = 1:numel(yCurves)
    y = yCurves{i}(:);
    if isempty(xCurves)
        x = (1:numel(y))';
    else
        x = xCurves{i}(:);
    end

    validateCurveLengths(x, y, legendNames{i}, plotTitle);
    y = sanitizeSemilogyData(y, plotConfiguration.logFloor, ...
        legendNames{i}, plotTitle);
    if useObjectiveTrace
        plotObjectiveTrace(plotAxes, x, y, styles, i, plotConfiguration);
    else
        [xPlot, yPlot] = downsampleCurve(x, y, plotConfiguration.maxLinePoints);
        lineHandle = semilogy(plotAxes, xPlot, yPlot, ...
            'Color', styles.colors(i, :), ...
            'LineStyle', styles.lineStyles{mod(i - 1, numel(styles.lineStyles)) + 1}, ...
            'LineWidth', styles.lineWidths(i), ...
            'Marker', styles.markers{mod(i - 1, numel(styles.markers)) + 1}, ...
            'MarkerSize', plotConfiguration.markerSize, ...
            'MarkerFaceColor', 'w', ...
            'MarkerEdgeColor', styles.colors(i, :));
        applyMarkerIndices(lineHandle, xPlot, plotConfiguration);
    end

    allX = [allX; x(isfinite(x))]; %#ok<AGROW>
    allY = [allY; y(isfinite(y) & y > 0)]; %#ok<AGROW>
end

formatCurveAxes(plotAxes, plotConfiguration);
xlabel(plotAxes, xLabelText, 'Interpreter', 'latex', 'FontSize', plotConfiguration.labelFontSize);
ylabel(plotAxes, yLabelText, 'Interpreter', 'latex', 'FontSize', plotConfiguration.labelFontSize);
if plotConfiguration.showTitle
    title(plotAxes, latexText(plotTitle), 'Interpreter', 'latex', ...
        'FontSize', plotConfiguration.titleFontSize, 'FontWeight', 'normal');
end

applyXLimits(plotAxes, allX);
if showTol && plotConfiguration.showToleranceLine
    allY = [allY; plotConfiguration.tolerance];
end
applyYLimits(plotAxes, allY, plotConfiguration);
if showTol && plotConfiguration.showToleranceLine
    drawToleranceLine(plotAxes, plotConfiguration.tolerance);
end

legendHandle = legend(plotAxes, legendNames, 'Interpreter', 'latex', ...
    'Location', plotConfiguration.legendLocation, 'Box', 'off');
legendHandle.FontSize = plotConfiguration.legendFontSize;
if isprop(legendHandle, 'Color')
    legendHandle.Color = 'none';
end
if plotConfiguration.legendNumColumns > 1 && isprop(legendHandle, 'NumColumns')
    legendHandle.NumColumns = plotConfiguration.legendNumColumns;
end

savePlot(figureHandle, savePrefix, suffix, plotConfiguration);
if plotConfiguration.closeFigures
    close(figureHandle);
end
end

function timeCurves = prepareTimeCurves(results)
%PREPARETIMECURVES Return per-algorithm CPU-time curves when available.
timeCurves = cell(numel(results), 1);
for i = 1:numel(results)
    history = results{i}.history;
    if ~isfield(history, 'time') || isempty(history.time)
        timeCurves = {};
        return;
    end
    timeCurve = history.time(:);
    if numel(timeCurve) ~= numel(history.obj)
        timeCurves = {};
        return;
    end
    timeCurves{i} = timeCurve;
end
end

function tf = allHaveHistoryField(results, fieldName)
%ALLHAVEHISTORYFIELD True when every result has a nonempty history field.
tf = true;
for i = 1:numel(results)
    history = results{i}.history;
    if ~isfield(history, fieldName) || isempty(history.(fieldName))
        tf = false;
        return;
    end
    if numel(history.(fieldName)) ~= numel(history.obj)
        tf = false;
        return;
    end
end
end

function plotLinearCurveFigure(xCurves, yCurves, legendNames, styles, plotTitle, ...
    xLabelText, yLabelText, savePrefix, suffix, plotConfiguration, curveKind)
%PLOTLINEARCURVEFIGURE Draw and save one linear-y curve figure.
if nargin < 11 || isempty(curveKind)
    curveKind = 'default';
end
figureHandle = squareFigure(plotTitle, plotConfiguration);
plotAxes = axes('Parent', figureHandle, ...
    'Position', curveAxesPosition(plotConfiguration, curveKind));
hold(plotAxes, 'on');

allX = [];
allY = [];
for i = 1:numel(yCurves)
    y = yCurves{i}(:);
    if isempty(xCurves)
        x = (1:numel(y))';
    else
        x = xCurves{i}(:);
    end

    validateCurveLengths(x, y, legendNames{i}, plotTitle);
    [xPlot, yPlot] = downsampleCurve(x, y, plotConfiguration.maxLinePoints);
    lineHandle = plot(plotAxes, xPlot, yPlot, ...
        'Color', styles.colors(i, :), ...
        'LineStyle', styles.lineStyles{mod(i - 1, numel(styles.lineStyles)) + 1}, ...
        'LineWidth', styles.lineWidths(i), ...
        'Marker', styles.markers{mod(i - 1, numel(styles.markers)) + 1}, ...
        'MarkerSize', plotConfiguration.markerSize, ...
        'MarkerFaceColor', 'w', ...
        'MarkerEdgeColor', styles.colors(i, :));
    applyMarkerIndices(lineHandle, xPlot, plotConfiguration);

    allX = [allX; x(isfinite(x))]; %#ok<AGROW>
    allY = [allY; y(isfinite(y))]; %#ok<AGROW>
end

formatLinearCurveAxes(plotAxes, plotConfiguration);
xlabel(plotAxes, xLabelText, 'Interpreter', 'latex', 'FontSize', plotConfiguration.labelFontSize);
ylabel(plotAxes, yLabelText, 'Interpreter', 'latex', 'FontSize', plotConfiguration.labelFontSize);
if plotConfiguration.showTitle
    title(plotAxes, latexText(plotTitle), 'Interpreter', 'latex', ...
        'FontSize', plotConfiguration.titleFontSize, 'FontWeight', 'normal');
end

applyXLimits(plotAxes, allX);
applyLinearYLimits(plotAxes, allY);

legendHandle = legend(plotAxes, legendNames, 'Interpreter', 'latex', ...
    'Location', plotConfiguration.legendLocation, 'Box', 'off');
legendHandle.FontSize = plotConfiguration.legendFontSize;
if isprop(legendHandle, 'Color')
    legendHandle.Color = 'none';
end
if plotConfiguration.legendNumColumns > 1 && isprop(legendHandle, 'NumColumns')
    legendHandle.NumColumns = plotConfiguration.legendNumColumns;
end

savePlot(figureHandle, savePrefix, suffix, plotConfiguration);
if plotConfiguration.closeFigures
    close(figureHandle);
end
end

function plotConfiguration = plotConfig(options)
%PLOT_CONFIG Build reusable plotting parameters.
if nargin < 1 || isempty(options)
    options = struct();
end
plotConfiguration.tolerance = getOption(options, 'tolerance', 1e-6);
plotConfiguration.maxLinePoints = getOption(options, 'plotMaxLinePoints', 2000);
plotConfiguration.maxYTicks = getOption(options, 'plotMaxYTicks', 7);
plotConfiguration.fontName = getOption(options, 'plotFontName', 'Times New Roman');
plotConfiguration.fontSize = getOption(options, 'plotFontSize', 12);
plotConfiguration.labelFontSize = getOption(options, 'plotLabelFontSize', 14);
plotConfiguration.titleFontSize = getOption(options, 'plotTitleFontSize', 14);
plotConfiguration.legendFontSize = getOption(options, 'plotLegendFontSize', 11);
plotConfiguration.showTitle = getLogicalOption(options, 'plotShowTitle', false);
plotConfiguration.pngDpi = getOption(options, 'plotPngDpi', 300);
plotConfiguration.savePdf = getLogicalOption(options, 'plotSavePdf', false);
plotConfiguration.showToleranceLine = getLogicalOption(options, ...
    'plotShowToleranceLine', true);
plotConfiguration.figureSizeInches = getFigureSize(options);
plotConfiguration.figurePixels = round(plotConfiguration.figureSizeInches * plotConfiguration.pngDpi);
plotConfiguration.lineWidth = getOption(options, 'plotLineWidth', 2.1);
plotConfiguration.highlightLineWidth = getOption(options, 'plotHighlightLineWidth', 2.7);
plotConfiguration.axisLineWidth = getOption(options, 'plotAxisLineWidth', 1.0);
plotConfiguration.markerSize = getOption(options, 'plotMarkerSize', 5.0);
plotConfiguration.useMarkers = getLogicalOption(options, 'plotUseMarkers', false);
plotConfiguration.markerCount = getOption(options, 'plotMarkerCount', 9);
plotConfiguration.objectiveEnhanced = getLogicalOption(options, 'plotObjectiveEnhanced', true);
plotConfiguration.objectiveSmoothingWindow = sanitizePositiveInteger( ...
    getOption(options, {'plotObjectiveSmoothingWindow', ...
    'plotRelativeErrorSmoothingWindow'}, 101), 101);
plotConfiguration.objectiveRawMarkers = getLogicalOption(options, ...
    'plotObjectiveRawMarkers', true);
plotConfiguration.objectiveRawMarkerMaxPoints = sanitizePositiveInteger( ...
    getOption(options, 'plotObjectiveRawMarkerMaxPoints', 220), 220);
plotConfiguration.objectiveRawMarkerStride = sanitizeOptionalPositiveInteger( ...
    getOption(options, 'plotObjectiveRawMarkerStride', []));
plotConfiguration.objectiveRawMarkerSize = getOption(options, ...
    'plotObjectiveRawMarkerSize', 11);
plotConfiguration.objectiveRawMarkerAlpha = getOption(options, ...
    'plotObjectiveRawMarkerAlpha', 0.24);
plotConfiguration.objectiveMinPoint = getLogicalOption(options, ...
    'plotObjectiveMinPoint', true);
plotConfiguration.objectiveMinMarkerSize = getOption(options, ...
    'plotObjectiveMinMarkerSize', 9.0);
plotConfiguration.objectiveMinMarkerLineWidth = getOption(options, ...
    'plotObjectiveMinMarkerLineWidth', 0.9);
plotConfiguration.objectiveAxesLeftMargin = getOption(options, ...
    'plotObjectiveAxesLeftMargin', 0.18);
plotConfiguration.logFloor = getOption(options, 'plotLogFloor', 1e-16);
if ~isnumeric(plotConfiguration.logFloor) || ~isscalar(plotConfiguration.logFloor) || ...
        ~isfinite(plotConfiguration.logFloor) || plotConfiguration.logFloor <= 0
    plotConfiguration.logFloor = 1e-16;
end
plotConfiguration.objectiveRawMarkerSize = sanitizePositiveScalar( ...
    plotConfiguration.objectiveRawMarkerSize, 11);
plotConfiguration.objectiveRawMarkerAlpha = sanitizeFraction( ...
    plotConfiguration.objectiveRawMarkerAlpha, 0.24);
plotConfiguration.objectiveMinMarkerSize = sanitizePositiveScalar( ...
    plotConfiguration.objectiveMinMarkerSize, 9.0);
plotConfiguration.objectiveMinMarkerLineWidth = sanitizePositiveScalar( ...
    plotConfiguration.objectiveMinMarkerLineWidth, 0.9);
plotConfiguration.objectiveAxesLeftMargin = sanitizeAxesMargin( ...
    plotConfiguration.objectiveAxesLeftMargin, 0.18);
plotConfiguration.legendLocation = getOption(options, 'plotLegendLocation', 'northeast');
plotConfiguration.legendNumColumns = getOption(options, 'plotLegendNumColumns', 1);
plotConfiguration.axesPosition = legendAwareAxesPosition(plotConfiguration.legendLocation, ...
    plotConfiguration.showTitle);
plotConfiguration.visible = 'off';
plotConfiguration.closeFigures = true;
end

function figureSize = getFigureSize(options)
%GETFIGURESIZE Return square, journal-friendly figure dimensions in inches.
figureSize = getOption(options, 'plotFigureSizeInches', []);
if isempty(figureSize)
    width = getOption(options, 'plotFigureWidthInches', []);
    height = getOption(options, 'plotFigureHeightInches', []);
    if isempty(width) && isempty(height)
        figureSize = [6.2 6.2];
    elseif isempty(height)
        figureSize = [width width];
    elseif isempty(width)
        figureSize = [height height];
    else
        figureSize = [width height];
    end
end
if ~isnumeric(figureSize) || numel(figureSize) ~= 2 || any(~isfinite(figureSize)) || ...
        any(figureSize <= 0)
    figureSize = [6.2 6.2];
end
figureSize = double(figureSize(:))';
end

function y = prepareLogCurve(y)
%PREPARE_LOG_CURVE Keep raw values so semilogy diagnostics expose bad data.
y = y(:);
end

function validateCurveLengths(x, y, curveName, plotTitle)
%VALIDATECURVELENGTHS Make curve/data mismatches explicit.
if numel(x) ~= numel(y)
    error('plotHistories:CurveLengthMismatch', ...
        'Curve "%s" in "%s" has %d x values and %d y values.', ...
        curveName, plotTitle, numel(x), numel(y));
end
end

function yPlot = sanitizeSemilogyData(y, logFloor, curveName, plotTitle)
%SANITIZESEMILOGYDATA Keep log plots drawable while exposing bad values.
yPlot = y(:);
positive = isfinite(yPlot) & yPlot > 0;
if ~any(positive)
    error('plotHistories:NoPositiveSemilogyData', ...
        'Curve "%s" in "%s" has no finite positive y values for a semilog-y plot.', ...
        curveName, plotTitle);
end

zeroOrNegative = isfinite(yPlot) & yPlot <= 0;
if any(zeroOrNegative)
    yPlot(zeroOrNegative) = logFloor;
end

nonFinite = ~isfinite(yPlot);
if any(nonFinite)
    yPlot(nonFinite) = NaN;
end
end

function styles = plotStyles(n, plotConfiguration)
%PLOT_STYLES Return a publication palette with Klein red first.
baseColors = [
    179  24  45
    225 107  36
    0 119 182
    0 150 136
    106  61 154
    232 171  45
    53  55  58
    86 180 233
    196  78 126
    ] / 255;
styles.colors = baseColors(mod(0:n - 1, size(baseColors, 1)) + 1, :);
styles.lineStyles = {'-', '-', '-', '-', '-', '--', '-.', '--', '-.'};
styles.markers = {'o', 's', '^', 'd', 'v', '>', '<'};
styles.lineWidths = repmat(plotConfiguration.lineWidth, n, 1);
if n >= 1
    styles.lineWidths(1) = plotConfiguration.highlightLineWidth;
end
end

function plotObjectiveTrace(plotAxes, x, y, styles, index, plotConfiguration)
%PLOTOBJECTIVETRACE Draw smoothed relative-error trend plus raw samples.
color = styles.colors(index, :);
marker = styles.markers{mod(index - 1, numel(styles.markers)) + 1};

if plotConfiguration.objectiveRawMarkers
    rawIdx = sampledPointIndices(numel(y), ...
        plotConfiguration.objectiveRawMarkerStride, ...
        plotConfiguration.objectiveRawMarkerMaxPoints);
    rawX = x(rawIdx);
    rawY = y(rawIdx);
    validRaw = isfinite(rawX) & isfinite(rawY) & rawY > 0;
    if any(validRaw)
        fadedColor = blendWithWhite(color, plotConfiguration.objectiveRawMarkerAlpha);
        scatter(plotAxes, rawX(validRaw), rawY(validRaw), ...
            plotConfiguration.objectiveRawMarkerSize, fadedColor, ...
            'filled', 'Marker', marker, ...
            'MarkerFaceColor', fadedColor, ...
            'MarkerEdgeColor', fadedColor, ...
            'HandleVisibility', 'off');
    end
end

smoothY = logMovingAverage(y, plotConfiguration.objectiveSmoothingWindow);
[xPlot, yPlot] = downsampleCurve(x, smoothY, plotConfiguration.maxLinePoints);
semilogy(plotAxes, xPlot, yPlot, ...
    'Color', color, ...
    'LineStyle', '-', ...
    'LineWidth', styles.lineWidths(index), ...
    'Marker', 'none');

if plotConfiguration.objectiveMinPoint
    plotMinimumPointMarker(plotAxes, x, y, color, plotConfiguration);
end
end

function ySmooth = logMovingAverage(y, windowSize)
%LOGMOVINGAVERAGE Smooth positive semilog data in log10 space.
ySmooth = y(:);
windowSize = min(sanitizePositiveInteger(windowSize, 1), numel(ySmooth));
if windowSize <= 1
    return;
end

logY = log10(ySmooth);
finiteLog = isfinite(logY);
logY(~finiteLog) = 0;
kernel = ones(windowSize, 1);
sumLog = conv(logY, kernel, 'same');
countLog = conv(double(finiteLog), kernel, 'same');
valid = countLog > 0;
ySmooth(:) = NaN;
ySmooth(valid) = 10 .^ (sumLog(valid) ./ countLog(valid));
end

function plotMinimumPointMarker(plotAxes, x, y, color, plotConfiguration)
%PLOTMINIMUMPOINTMARKER Mark the best observed relative error.
valid = isfinite(x) & isfinite(y) & y > 0;
if ~any(valid)
    return;
end
validIndex = find(valid);
[minY, localIndex] = min(y(valid));
minIndex = validIndex(localIndex);
plot(plotAxes, x(minIndex), minY, 'p', ...
    'Color', [0.08 0.08 0.08], ...
    'MarkerFaceColor', color, ...
    'MarkerEdgeColor', [0.08 0.08 0.08], ...
    'MarkerSize', plotConfiguration.objectiveMinMarkerSize, ...
    'LineWidth', plotConfiguration.objectiveMinMarkerLineWidth, ...
    'HandleVisibility', 'off');
end

function indexSet = sampledPointIndices(pointCount, stride, maxPoints)
%SAMPLEDPOINTINDICES Return sparse raw-data marker locations.
if pointCount <= 0
    indexSet = [];
    return;
end
if isempty(stride)
    stride = max(1, ceil(pointCount / max(maxPoints, 1)));
else
    stride = max(1, stride);
end
indexSet = 1:stride:pointCount;
if indexSet(end) ~= pointCount
    indexSet = [indexSet pointCount];
end
end

function fadedColor = blendWithWhite(color, alpha)
%BLENDWITHWHITE Use stable white blending instead of renderer transparency.
alpha = sanitizeFraction(alpha, 0.24);
fadedColor = 1 - alpha * (1 - color);
fadedColor = min(max(fadedColor, 0), 1);
end

function savePrefix = defaultSavePrefix(baseName)
%DEFAULTSAVEPREFIX Save figures to the project-level Results folder.
toolDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(toolDir);
savePrefix = fullfile(projectDir, 'Results', baseName);
end

function figureHandle = squareFigure(name, plotConfiguration)
%SQUARE_FIGURE Create a hidden square figure suitable for exported plots.
visibleState = plotConfiguration.visible;
positionPixels = [100 100 plotConfiguration.figurePixels];
paperPosition = [0 0 plotConfiguration.figureSizeInches];
figureHandle = figure('Color', 'w', 'Name', name, ...
    'Visible', visibleState, ...
    'Units', 'pixels', 'Position', positionPixels, ...
    'PaperUnits', 'inches', 'PaperPosition', paperPosition, ...
    'PaperSize', plotConfiguration.figureSizeInches, ...
    'InvertHardcopy', 'off');
end

function formatCurveAxes(plotAxes, plotConfiguration)
%FORMAT_CURVE_AXES Apply consistent publication-style axis formatting.
grid(plotAxes, 'on');
box(plotAxes, 'on');
plotAxes.YScale = 'log';
plotAxes.Color = 'w';
plotAxes.FontName = plotConfiguration.fontName;
plotAxes.FontSize = plotConfiguration.fontSize;
plotAxes.LineWidth = plotConfiguration.axisLineWidth;
plotAxes.TickLabelInterpreter = 'latex';
plotAxes.TickDir = 'in';
plotAxes.Layer = 'top';
plotAxes.GridAlpha = 0.18;
plotAxes.MinorGridAlpha = 0.08;
plotAxes.XMinorGrid = 'off';
plotAxes.YMinorGrid = 'off';
plotAxes.XColor = [0.10 0.10 0.10];
plotAxes.YColor = [0.10 0.10 0.10];
end

function formatLinearCurveAxes(plotAxes, plotConfiguration)
%FORMATLINEARCURVEAXES Apply the same style with a linear y-axis.
grid(plotAxes, 'on');
box(plotAxes, 'on');
plotAxes.YScale = 'linear';
plotAxes.Color = 'w';
plotAxes.FontName = plotConfiguration.fontName;
plotAxes.FontSize = plotConfiguration.fontSize;
plotAxes.LineWidth = plotConfiguration.axisLineWidth;
plotAxes.TickLabelInterpreter = 'latex';
plotAxes.TickDir = 'in';
plotAxes.Layer = 'top';
plotAxes.GridAlpha = 0.18;
plotAxes.MinorGridAlpha = 0.08;
plotAxes.XMinorGrid = 'off';
plotAxes.YMinorGrid = 'off';
plotAxes.XColor = [0.10 0.10 0.10];
plotAxes.YColor = [0.10 0.10 0.10];
end

function applyXLimits(plotAxes, allX)
%APPLY_X_LIMITS Set a stable x-axis range from valid curve data.
if isempty(allX)
    return;
end
xMin = min(allX);
xMax = max(allX);
if xMin == xMax
    xMax = xMin + 1;
end
padding = 0.015 * max(xMax - xMin, 1);
xlim(plotAxes, [xMin, xMax + padding]);
end

function applyYLimits(plotAxes, allY, plotConfiguration)
%APPLY_Y_LIMITS Set data-adaptive logarithmic y-axis limits and ticks.
validY = allY(isfinite(allY) & allY > 0);
if isempty(validY)
    return;
end
minExp = floor(log10(min(validY)));
maxExp = ceil(log10(max(validY)));
if minExp == maxExp
    minExp = minExp - 1;
    maxExp = maxExp + 1;
end
yMin = 10 ^ minExp;
yMax = 10 ^ maxExp;
ylim(plotAxes, [yMin, yMax]);
applyAdaptiveYTicks(plotAxes, minExp, maxExp, plotConfiguration.maxYTicks);
end

function applyLinearYLimits(plotAxes, allY)
%APPLYLINEARYLIMITS Add a small padding for linear diagnostic curves.
validY = allY(isfinite(allY));
if isempty(validY)
    return;
end
yMin = min(validY);
yMax = max(validY);
if yMin == yMax
    padding = max(abs(yMin) * 0.05, 1);
else
    padding = 0.06 * (yMax - yMin);
end
ylim(plotAxes, [yMin - padding, yMax + padding]);
end

function applyAdaptiveYTicks(plotAxes, minExp, maxExp, maxTicks)
%APPLY_ADAPTIVE_Y_TICKS Keep semilog tick precision aligned with data scale.
span = maxExp - minExp;
step = max(1, ceil(span / maxTicks));
tickExp = minExp:step:maxExp;
labels = arrayfun(@(e) sprintf('$10^{%d}$', e), tickExp, ...
    'UniformOutput', false);
set(plotAxes, 'YTick', 10 .^ tickExp, 'YTickLabel', labels);
end

function drawToleranceLine(plotAxes, tolerance)
%DRAW_TOLERANCE_LINE Draw a horizontal stopping-tolerance reference.
yl = ylim(plotAxes);
if tolerance < yl(1) || tolerance > yl(2)
    return;
end
xl = xlim(plotAxes);
plot(plotAxes, xl, [tolerance tolerance], ':', ...
    'Color', [0.25 0.25 0.25], 'LineWidth', 1.0, ...
    'HandleVisibility', 'off');
text(plotAxes, xl(1) + 0.03 * diff(xl), tolerance * 1.12, '$\varepsilon$', ...
    'Interpreter', 'latex', 'FontSize', 9, ...
    'Color', [0.25 0.25 0.25], ...
    'VerticalAlignment', 'bottom');
end

function applyMarkerIndices(lineHandle, xPlot, plotConfiguration)
%APPLYMARKERINDICES Optionally add sparse markers without cluttering curves.
if ~plotConfiguration.useMarkers
    lineHandle.Marker = 'none';
    return;
end
if isprop(lineHandle, 'MarkerIndices') && ~isempty(xPlot)
    markerCount = min(plotConfiguration.markerCount, numel(xPlot));
    lineHandle.MarkerIndices = unique(round(linspace(1, numel(xPlot), markerCount)));
elseif numel(xPlot) > max(plotConfiguration.markerCount, 1)
    lineHandle.Marker = 'none';
end
end

function [xPlot, yPlot] = downsampleCurve(x, y, maxPoints)
%DOWNSAMPLE_CURVE Keep long histories readable and export files compact.
if numel(y) <= maxPoints
    xPlot = x;
    yPlot = y;
    return;
end
indexSet = unique(round(linspace(1, numel(y), maxPoints)));
xPlot = x(indexSet);
yPlot = y(indexSet);
end

function value = sanitizePositiveInteger(value, defaultValue)
%SANITIZEPOSITIVEINTEGER Keep plotting options in a valid integer range.
if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) || value < 1
    value = defaultValue;
end
value = max(1, round(double(value)));
end

function value = sanitizeOptionalPositiveInteger(value)
%SANITIZEOPTIONALPOSITIVEINTEGER Return [] or a positive integer.
if isempty(value)
    return;
end
if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) || value < 1
    value = [];
else
    value = max(1, round(double(value)));
end
end

function value = sanitizePositiveScalar(value, defaultValue)
%SANITIZEPOSITIVESCALAR Keep scalar plotting sizes positive and finite.
if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) || value <= 0
    value = defaultValue;
end
value = double(value);
end

function value = sanitizeFraction(value, defaultValue)
%SANITIZEFRACTION Clamp white-blend strength to [0, 1].
if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value)
    value = defaultValue;
end
value = min(max(double(value), 0), 1);
end

function savePlot(figureHandle, savePrefix, suffix, plotConfiguration)
%SAVE_PLOT Save a figure as PNG and optionally PDF.
if isempty(savePrefix)
    return;
end
outDir = fileparts(savePrefix);
if ~isempty(outDir) && ~exist(outDir, 'dir')
    mkdir(outDir);
end
set(figureHandle, 'Renderer', 'painters');
pngFile = [savePrefix '_' suffix '.png'];
pdfFile = [savePrefix '_' suffix '.pdf'];
if canUseExportGraphics()
    exportgraphics(figureHandle, pngFile, ...
        'Resolution', plotConfiguration.pngDpi, 'BackgroundColor', 'white');
else
    print(figureHandle, pngFile, '-dpng', ['-r' num2str(plotConfiguration.pngDpi)]);
end
if plotConfiguration.savePdf
    if canUseExportGraphics()
        exportgraphics(figureHandle, pdfFile, 'ContentType', 'vector', ...
            'BackgroundColor', 'white');
    else
        print(figureHandle, pdfFile, '-dpdf', '-painters');
    end
end
end

function tf = canUseExportGraphics()
%CANUSEEXPORTGRAPHICS True for MATLAB releases that provide exportgraphics.
tf = exist('exportgraphics', 'file') == 2 || exist('exportgraphics', 'builtin') == 5;
end

function value = getOption(options, names, defaultValue)
%GETOPTION Read the first available field from a struct.
if ischar(names) || isStringScalar(names)
    names = {char(names)};
end
value = defaultValue;
for i = 1:numel(names)
    fieldName = names{i};
    if isfield(options, fieldName) && ~isempty(options.(fieldName))
        value = options.(fieldName);
        return;
    end
end
end

function value = getLogicalOption(options, name, defaultValue)
%GETLOGICALOPTION Accept logical, numeric, and on/off style option values.
rawValue = getOption(options, name, defaultValue);
if islogical(rawValue) && isscalar(rawValue)
    value = rawValue;
elseif isnumeric(rawValue) && isscalar(rawValue)
    value = rawValue ~= 0;
elseif ischar(rawValue) || isStringScalar(rawValue)
    rawValue = lower(strtrim(char(rawValue)));
    value = any(strcmp(rawValue, {'true', 'on', 'yes', '1'}));
else
    value = defaultValue;
end
end

function position = legendAwareAxesPosition(location, showTitle)
%LEGENDAWAREAXESPOSITION Leave room when the legend is outside the axes.
location = lower(char(location));
if showTitle
    topMargin = 0.12;
else
    topMargin = 0.07;
end
if hasText(location, 'outside')
    if hasText(location, 'east') || hasText(location, 'west')
        position = [0.13 0.13 0.66 1 - 0.13 - topMargin];
    elseif hasText(location, 'south')
        position = [0.13 0.23 0.82 1 - 0.23 - topMargin];
    else
        leftMargin = 0.20;      % 左侧留出 20% 空间，给长标签腾地
        bottomMargin = 0.13;
        width = 1 - leftMargin - 0.05;   % 右侧留 5% 避免图例/边缘裁剪
        height = 1 - bottomMargin - topMargin;
        position = [leftMargin, bottomMargin, width, height];
    end
else
    position = [0.13 0.13 0.82 1 - 0.13 - topMargin];
end
end

function position = curveAxesPosition(plotConfiguration, curveKind)
%CURVEAXESPOSITION Give the objective label enough room after rotation.
position = plotConfiguration.axesPosition;
if ~strcmpi(char(curveKind), 'objective')
    return;
end
rightEdge = position(1) + position(3);
position(1) = max(position(1), plotConfiguration.objectiveAxesLeftMargin);
position(3) = max(0.10, rightEdge - position(1));
end

function value = sanitizeAxesMargin(value, defaultValue)
%SANITIZEAXESMARGIN Keep normalized axes margins inside a usable range.
if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value)
    value = defaultValue;
end
value = min(max(double(value), 0.05), 0.45);
end

function tf = hasText(textValue, pattern)
%HASTEXT Compatibility wrapper for simple substring checks.
tf = ~isempty(strfind(textValue, pattern)); %#ok<STREMP>
end

% --- 修正后的 latexText：保护已有的数学模式 ---
function s = latexText(s)
%LATEX_TEXT Escape plain text for MATLAB LaTeX interpreter.
%   If the string already contains a '$' (i.e. it already has
%   inline math), it is left untouched; otherwise backslashes
%   and underscores are escaped.
if contains(s, '$')
    % 已经含有数学模式，直接返回
    return;
end
s = strrep(s, '\', '\textbackslash{}');
s = strrep(s, '_', '\_');
end
