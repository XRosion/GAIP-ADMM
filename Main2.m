% Main2: GAIP-ADMM one-parameter-at-a-time sweeps for LASSO.
% This script solves
%   min_x eta*||x||_1 + 0.5*||K*x-b||_2^2
% while varying theta_x, theta_y, and theta_lambda separately.

clc; clear;
import Algorithms.*
import Functions.*
import Problems.*
import Tools.*

rootDir = fileparts(mfilename('fullpath'));
resultsDir = fullfile(rootDir, 'Results');
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

fprintf('========== LASSO problem generation ==========\n');
observationDimension = 640;
variableDimension = 2000;
sparsityLevel = 100;
Problem2 = Problems.lassoProblem( ...
    observationDimension, variableDimension, sparsityLevel);
PhiStar = Problem2.PhiStar;
fprintf('CVX reference optimum Phi* = %.12g\n', PhiStar);
fprintf('LASSO data scale: ||K|| = %.4g, true-x objective = %.12g\n', ...
    Problem2.normK, ...
    Functions.objectiveValue(Problem2, Problem2.xtrue, Problem2.xtrue));

%% Experiment and linear-solve options
options.maxIterations = 5000;
options.tolerance = 1e-6;
options.keepOnlyVisualizationData = true;

% Match the Elastic-net publication-style visualization: square 300-DPI
% figures, an internal legend, and a smoothed relative-error trend with
% faint raw samples retained to show the unsmoothed behavior.
options.plotFigureSizeInches = [6.2, 6.2];
options.plotPngDpi = 300;
options.plotLegendLocation = 'northeast';
options.plotLegendNumColumns = 1;
options.plotShowTitle = false;
options.plotUseMarkers = false;
options.plotMaxLinePoints = 2000;
options.plotObjectiveEnhanced = true;
options.plotObjectiveSmoothingWindow = 101;
options.plotObjectiveRawMarkers = true;
options.plotObjectiveRawMarkerMaxPoints = 180;
options.plotObjectiveRawMarkerAlpha = 0.22;
options.plotObjectiveMinPoint = false;
options.plotSavePdf = false;
options.plotObjectiveSuffix = 'relative_error';
options.plotFeasibilitySuffix = 'constraint_violation';

% The smooth second-block CG solve starts from y^k in GAIPADMM.lasso.
% It stops at (2.10a), or when the absolute residual reaches 1e-10,
% and performs at most 300 inner iterations.
options.lassoCgMaxIter = 300;
options.lassoCgAbsTol = 1e-10;

thetaXValues = [-0.99, -0.8, -0.618, -0.3, 0, 0.3, 0.618, 0.8, 0.99];
thetaYValues = [-0.99, -0.8, -0.618, -0.3, 0, 0.3, 0.618, 0.8, 0.99];
% With the fixed alpha=1.618 and eps3=0, these values satisfy vartheta > 0.
thetaLambdaValues = [0.4, 0.44, 0.5, 0.7, 0.9];

% Baseline parameters from the LASSO experiment specification. Each sweep
% changes exactly one theta parameter while all other values remain fixed.
baselineParameters.beta = 1;
baselineParameters.tau1 = 0;
baselineParameters.tau2 = 0.01;
baselineParameters.eps1 = 0;
baselineParameters.eps2 = 0.99;
baselineParameters.eps3 = 0;
baselineParameters.theta_x = -0.99;
baselineParameters.theta_y = -0.99;
baselineParameters.theta_lambda = 0.44;
baselineParameters.alpha = 1.618;

% Vary one parameter at a time so all values of the same parameter can be
% compared in one pair of iteration-history figures.
numThetaXRuns = numel(thetaXValues);
numThetaYRuns = numel(thetaYValues);
numThetaLambdaRuns = numel(thetaLambdaValues);
thetaXRunIndices = 1:numThetaXRuns;
thetaYRunIndices = numThetaXRuns + (1:numThetaYRuns);
thetaLambdaRunIndices = numThetaXRuns + numThetaYRuns + ...
    (1:numThetaLambdaRuns);

thetaXList = [thetaXValues(:); ...
    repmat(baselineParameters.theta_x, numThetaYRuns + numThetaLambdaRuns, 1)];
thetaYList = [repmat(baselineParameters.theta_y, numThetaXRuns, 1); ...
    thetaYValues(:); ...
    repmat(baselineParameters.theta_y, numThetaLambdaRuns, 1)];
thetaLambdaList = [repmat(baselineParameters.theta_lambda, ...
    numThetaXRuns + numThetaYRuns, 1); thetaLambdaValues(:)];
numRuns = numel(thetaXList);
alphaList = repmat(baselineParameters.alpha, numRuns, 1);

varthetaList = 1 - alphaList .* (1 - thetaLambdaList) - ...
    baselineParameters.eps3;
if any(varthetaList <= 0)
    invalidIndex = find(varthetaList <= 0, 1);
    error('Main2:InvalidMultiplierExtrapolation', ...
        ['Run %d violates vartheta > 0: alpha=%g, theta_lambda=%g, ' ...
        'eps3=%g, vartheta=%g.'], ...
        invalidIndex, alphaList(invalidIndex), ...
        thetaLambdaList(invalidIndex), baselineParameters.eps3, ...
        varthetaList(invalidIndex));
end

sweepGroups = repmat(struct( ...
    'parameterName', '', 'values', [], 'runIndices', []), 3, 1);
sweepGroups(1).parameterName = 'theta_x';
sweepGroups(1).values = thetaXValues;
sweepGroups(1).runIndices = thetaXRunIndices;
sweepGroups(2).parameterName = 'theta_y';
sweepGroups(2).values = thetaYValues;
sweepGroups(2).runIndices = thetaYRunIndices;
sweepGroups(3).parameterName = 'theta_lambda';
sweepGroups(3).values = thetaLambdaValues;
sweepGroups(3).runIndices = thetaLambdaRunIndices;

runParameters = repmat(baselineParameters, numRuns, 1);
for i = 1:numRuns
    runParameters(i).theta_x = thetaXList(i);
    runParameters(i).theta_y = thetaYList(i);
    runParameters(i).theta_lambda = thetaLambdaList(i);
    runParameters(i).alpha = alphaList(i);
end

sweepDefinition = struct( ...
    'sweepMode', 'one-parameter-at-a-time', ...
    'thetaXValues', thetaXValues, ...
    'thetaYValues', thetaYValues, ...
    'thetaLambdaValues', thetaLambdaValues, ...
    'baselineParameters', baselineParameters, ...
    'maxOuterIterations', options.maxIterations, ...
    'cgMaxIterations', options.lassoCgMaxIter, ...
    'cgAbsoluteResidualTolerance', options.lassoCgAbsTol);

summaryTable = initializeSummaryTable( ...
    thetaXList, thetaYList, thetaLambdaList, alphaList);
results = cell(numRuns, 1);
elapsedTimes = nan(numRuns, 1);
finishedRuns = false(numRuns, 1);

checkpointFile = fullfile(resultsDir, 'lasso_parameter_sweep_checkpoint.mat');
summaryCsvFile = fullfile(resultsDir, 'lasso_parameter_sweep_summary.csv');
checkpointEvery = 5;

if exist(checkpointFile, 'file')
    checkpoint = load(checkpointFile);
    if ~isfield(checkpoint, 'sweepDefinition') || ...
            ~isequaln(checkpoint.sweepDefinition, sweepDefinition)
        error('Main2:CheckpointMismatch', ...
            ['The existing checkpoint was created with a different sweep ' ...
            'configuration: %s'], checkpointFile);
    end

    summaryTable = checkpoint.summaryTable;
    results = checkpoint.results;
    elapsedTimes = checkpoint.elapsedTimes;
    finishedRuns = checkpoint.finishedRuns;
    fprintf('Resumed %d/%d attempted runs from checkpoint.\n', ...
        nnz(finishedRuns), numRuns);
end

fprintf('\n========== GAIP-ADMM parameter sweep ==========\n');
fprintf(['Runs: %d = %d theta_x values + %d theta_y values + ' ...
    '%d theta_lambda values (one parameter at a time)\n'], ...
    numRuns, numel(thetaXValues), numel(thetaYValues), ...
    numel(thetaLambdaValues));
fprintf('Every run uses %d outer iterations and at most %d CG iterations.\n', ...
    options.maxIterations, options.lassoCgMaxIter);

attemptsSinceCheckpoint = 0;
for i = 1:numRuns
    if finishedRuns(i)
        continue;
    end

    parameters = runParameters(i);
    runName = parameterRunName(parameters);
    fprintf('\n>>> Run %d/%d: %s, alpha=%.8g\n', ...
        i, numRuns, runName, parameters.alpha);

    runTimer = tic;
    try
        out = Algorithms.GAIPADMM.lasso(Problem2, options, parameters);
        elapsedTimes(i) = toc(runTimer);

        out.name = runName;
        out.displayName = runName;
        out.parameters = parameters;

        history = out.history;
        if numel(history.obj) ~= options.maxIterations
            error('Main2:UnexpectedOuterIterationCount', ...
                'Expected %d outer iterations, but the solver returned %d.', ...
                options.maxIterations, numel(history.obj));
        end
        results{i} = out;
        relativeErrors = abs(history.obj - PhiStar) / (abs(PhiStar) + 1e-12);
        [minimumRelativeError, bestIteration] = min(relativeErrors);

        summaryTable.OuterIterations(i) = numel(history.obj);
        summaryTable.FinalObjective(i) = history.obj(end);
        summaryTable.RelativeObjectiveError(i) = relativeErrors(end);
        summaryTable.MinimumRelativeObjectiveError(i) = minimumRelativeError;
        summaryTable.BestObjectiveIteration(i) = bestIteration;
        summaryTable.FinalFeasibility(i) = history.feas(end);
        summaryTable.MeetsFeasibilityTolerance(i) = ...
            history.feas(end) <= options.tolerance;
        summaryTable.ElapsedSeconds(i) = elapsedTimes(i);
        summaryTable.Status{i} = 'completed';
        summaryTable.ErrorMessage{i} = '';

        fprintf(['    objective=%.6e, rel.error=%.3e, feasibility=%.3e, ' ...
            'time=%.2fs\n'], ...
            summaryTable.FinalObjective(i), ...
            summaryTable.RelativeObjectiveError(i), ...
            summaryTable.FinalFeasibility(i), elapsedTimes(i));
    catch runError
        elapsedTimes(i) = toc(runTimer);
        summaryTable.ElapsedSeconds(i) = elapsedTimes(i);
        summaryTable.Status{i} = 'failed';
        summaryTable.ErrorMessage{i} = getReport( ...
            runError, 'extended', 'hyperlinks', 'off');
        warning('Main2:RunFailed', ...
            'Run %d/%d failed: %s', i, numRuns, runError.message);
    end

    finishedRuns(i) = true;
    attemptsSinceCheckpoint = attemptsSinceCheckpoint + 1;
    if attemptsSinceCheckpoint >= checkpointEvery || all(finishedRuns)
        saveSweepCheckpoint(checkpointFile, sweepDefinition, summaryTable, ...
            results, elapsedTimes, finishedRuns);
        writetable(summaryTable, summaryCsvFile);
        attemptsSinceCheckpoint = 0;
        fprintf('    checkpoint saved (%d/%d attempted).\n', ...
            nnz(finishedRuns), numRuns);
    end
end

%% Save the comparison and one curve pair per swept parameter
completedMask = strcmp(summaryTable.Status, 'completed');
completedSummary = summaryTable(completedMask, :);
feasibleSummary = sortrows( ...
    completedSummary(completedSummary.MeetsFeasibilityTolerance, :), ...
    {'RelativeObjectiveError', 'FinalFeasibility'});
infeasibleSummary = sortrows( ...
    completedSummary(~completedSummary.MeetsFeasibilityTolerance, :), ...
    {'FinalFeasibility', 'RelativeObjectiveError'});
rankedSummary = [feasibleSummary; infeasibleSummary];
rankedCsvFile = fullfile(resultsDir, 'lasso_parameter_sweep_ranked.csv');
writetable(rankedSummary, rankedCsvFile);

finalMatFile = fullfile(resultsDir, 'lasso_parameter_sweep_results.mat');
save(finalMatFile, 'Problem2', 'sweepDefinition', 'summaryTable', ...
    'rankedSummary', 'results', 'elapsedTimes', '-v7.3');

saveParameterSweepCurves(results, completedMask, sweepGroups, Problem2, ...
    options, resultsDir);

fprintf('\n========== Sweep complete ==========\n');
fprintf('Completed runs: %d/%d; failed runs: %d.\n', ...
    nnz(completedMask), numRuns, nnz(strcmp(summaryTable.Status, 'failed')));
if ~isempty(rankedSummary)
    fprintf(['\nTop parameter combinations (feasibility tolerance first, ' ...
        'then final relative objective error):\n']);
    disp(rankedSummary(1:min(10, height(rankedSummary)), ...
        {'ThetaX', 'ThetaY', 'ThetaLambda', 'Alpha', ...
        'RelativeObjectiveError', 'FinalFeasibility', 'ElapsedSeconds'}));
end
fprintf('Results saved in:\n  %s\n', resultsDir);

%% Local helper functions
function summaryTable = initializeSummaryTable( ...
    thetaX, thetaY, thetaLambda, alpha)
numRuns = numel(thetaX);
Run = (1:numRuns)';
ThetaX = thetaX;
ThetaY = thetaY;
ThetaLambda = thetaLambda;
Alpha = alpha;
OuterIterations = zeros(numRuns, 1);
FinalObjective = nan(numRuns, 1);
RelativeObjectiveError = nan(numRuns, 1);
MinimumRelativeObjectiveError = nan(numRuns, 1);
BestObjectiveIteration = nan(numRuns, 1);
FinalFeasibility = nan(numRuns, 1);
MeetsFeasibilityTolerance = false(numRuns, 1);
ElapsedSeconds = nan(numRuns, 1);
Status = repmat({'pending'}, numRuns, 1);
ErrorMessage = repmat({''}, numRuns, 1);

summaryTable = table(Run, ThetaX, ThetaY, ThetaLambda, Alpha, ...
    OuterIterations, FinalObjective, RelativeObjectiveError, ...
    MinimumRelativeObjectiveError, BestObjectiveIteration, ...
    FinalFeasibility, MeetsFeasibilityTolerance, ElapsedSeconds, ...
    Status, ErrorMessage);
end

function name = parameterRunName(parameters)
name = sprintf('theta_x=%+.3f, theta_y=%+.3f, theta_lambda=%.1f', ...
    parameters.theta_x, parameters.theta_y, parameters.theta_lambda);
end

function saveSweepCheckpoint(fileName, sweepDefinition, summaryTable, ...
    results, elapsedTimes, finishedRuns)
save(fileName, 'sweepDefinition', 'summaryTable', 'results', ...
    'elapsedTimes', 'finishedRuns', '-v7.3');
end

function saveParameterSweepCurves(results, completedMask, sweepGroups, ...
    problem, options, outputDir)
if ~any(completedMask)
    return;
end

fprintf('\nGenerating six parameter-sweep figures...\n');
for groupIndex = 1:numel(sweepGroups)
    group = sweepGroups(groupIndex);
    groupResults = {};
    for valueIndex = 1:numel(group.runIndices)
        runIndex = group.runIndices(valueIndex);
        if ~completedMask(runIndex) || isempty(results{runIndex})
            warning('Main2:MissingParameterCurve', ...
                'Skipping %s=%g because run %d did not complete.', ...
                group.parameterName, group.values(valueIndex), runIndex);
            continue;
        end

        plotResult = results{runIndex};
        plotResult.name = parameterLegendName( ...
            group.parameterName, group.values(valueIndex));
        if isfield(plotResult.history, 'time')
            plotResult.history = rmfield(plotResult.history, 'time');
        end
        groupResults{end + 1, 1} = plotResult; %#ok<AGROW>
    end

    if isempty(groupResults)
        warning('Main2:EmptyParameterSweep', ...
            'No completed results are available for %s.', group.parameterName);
        continue;
    end

    figureTitle = sprintf('LASSO %s sweep', group.parameterName);
    savePrefix = fullfile(outputDir, group.parameterName);
    Tools.plotHistories(groupResults, problem, figureTitle, ...
        savePrefix, options);
    fprintf('    saved relative-error and constraint-violation figures for %s.\n', ...
        group.parameterName);
end
end

function name = parameterLegendName(parameterName, value)
switch parameterName
    case 'theta_x'
        name = sprintf('$\\theta_x = %.3g$', value);
    case 'theta_y'
        name = sprintf('$\\theta_w = %.3g$', value);
    case 'theta_lambda'
        name = sprintf('$\\theta_\\lambda = %.3g$', value);
    otherwise
        name = sprintf('%s = %.3g', parameterName, value);
end
end
