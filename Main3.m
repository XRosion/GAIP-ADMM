% Main3: Elastic-net numerical experiment main script.
% This script compares ADMM variants on
%   min_x eta*||x||_1 + gamma*||x||_2^2 + 0.5*||K*x-b||_2^2.
%
% The elastic-net implementations used here follow the split
%   min_{x,y} gamma*||x||_2^2 + 0.5*||K*x-b||_2^2 + eta*||y||_1,
%   s.t. x - y = 0.

clc; clear all;
import Algorithms.*
import Functions.*
import Problems.*
import Tools.*

rootDir = fileparts(mfilename('fullpath'));
resultsDir = fullfile(rootDir, 'Results');
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

fprintf('========== Elastic-net problem generation ==========\n');
observationDimension = 640;
variableDimension = 2000;
sparsityLevel = 100;
Problem3 = Problems.elasticProblem(observationDimension, variableDimension, sparsityLevel);
PhiStar = Problem3.PhiStar;
fprintf('CVX reference optimum Phi* = %.12g\n', PhiStar);
fprintf('Elastic-net data scale: ||K|| = %.4g, true-x objective = %.12g\n', ...
        Problem3.normK, Functions.objectiveValue(Problem3, Problem3.xtrue, Problem3.xtrue));
fprintf('Regularization parameters: eta = %.4g, gamma = %.4g\n', ...
        Problem3.eta, Problem3.gamma);

options.maxIterations = 3000;
options.tolerance = 1e-6;
options.keepOnlyVisualizationData = true;       % 仅保留绘图和汇总所需历史，降低长迭代实验的内存占用。
options.plotObjectiveSmoothingWindow = 101;     % 相对误差图的 log10 移动平均窗口，奇数窗口便于居中平滑。
options.plotObjectiveRawMarkerMaxPoints = 250;  % 原始相对误差散点最多显示数量，用稀疏散点体现真实波动。
options.plotObjectiveRawMarkerAlpha = 0.4;     % 原始散点淡化强度，数值越小越接近白色、视觉干扰越弱。
options.plotObjectiveMinPoint = false;           % 用特殊标记突出每个算法观测到的最小相对误差点。

% GAIPADMM and AGAIPADMM update the nonsmooth auxiliary block first and
% solve the second, smooth signal block by CG.
options.elasticCgMaxIter = 300;
options.elasticCgAbsTol = 1e-10;

algorithmParameters.alpha = 1.618;
algorithmParameters.theta_lambda = 0.44;
algorithmParameters.beta = 0.8;
algorithmParameters.theta_x = -0.618;
algorithmParameters.theta_y = -0.618;

% Block parameters follow Algorithm 1: first=x (nonsmooth auxiliary),
% second=y (smooth signal).
algorithmParameters.tau1 = 0;
algorithmParameters.tau2 = 0.01;
algorithmParameters.eps2 = 0.99;
algorithmParameters.eps3 = 0;

acceleratedParameters = algorithmParameters;
acceleratedParameters.theta_gamma = 1;
acceleratedParameters.beta = 0.2;

Algorithms = struct('Name', {}, 'Func', {}, 'Type', {}, 'Parameters', {}, 'Description', {});
algorithmIndex = 0;

algorithmIndex = algorithmIndex + 1;
Algorithms(algorithmIndex).Name = 'AGAIP-ADMM_{spp}';
Algorithms(algorithmIndex).Func = @Algorithms.AGAIPADMM.elastic;
Algorithms(algorithmIndex).Type = 1;
Algorithms(algorithmIndex).Parameters = acceleratedParameters;
Algorithms(algorithmIndex).Description = 'AGAIP-ADMM_spp with increasing penalty';


algorithmIndex = algorithmIndex + 1;
Algorithms(algorithmIndex).Name = 'GAIP-ADMM_{spp}';
Algorithms(algorithmIndex).Func = @Algorithms.GAIPADMM.elastic;
Algorithms(algorithmIndex).Type = 1;
Algorithms(algorithmIndex).Parameters = algorithmParameters;
Algorithms(algorithmIndex).Description = 'Fixed-penalty GAIP-ADMM_spp';


algorithmIndex = algorithmIndex + 1;
Algorithms(algorithmIndex).Name = 'IADMM';
Algorithms(algorithmIndex).Func = @Algorithms.IADMM.elastic;
Algorithms(algorithmIndex).Type = 0;
Algorithms(algorithmIndex).Description = 'ADMM with exact smooth elastic-net block';

algorithmIndex = algorithmIndex + 1;
Algorithms(algorithmIndex).Name = 'IGADMM';
Algorithms(algorithmIndex).Func = @Algorithms.IGADMM.elastic;
Algorithms(algorithmIndex).Type = 0;
Algorithms(algorithmIndex).Description = 'Generalized ADMM with exact elastic-net block';

fprintf('\n========== Start solving ==========\n');
numAlgorithms = numel(Algorithms);
results = cell(numAlgorithms, 1);
elapsedTimes = zeros(numAlgorithms, 1);

for i = 1:numAlgorithms
    alg = Algorithms(i);
    fprintf('\n>>> Algorithm %d/%d : %s (%s)\n', ...
        i, numAlgorithms, alg.Name, alg.Description);

    tic;
    if alg.Type == 0
        out = alg.Func(Problem3, options);
    else
        out = alg.Func(Problem3, options, alg.Parameters);
    end
    elapsedTimes(i) = toc;

    out.name = alg.Name;
    out.displayName = sprintf('%s (%s)', alg.Name, alg.Description);
    printSolverStatus(out, PhiStar, elapsedTimes(i), options.tolerance);
    results{i} = out;
end

fprintf('\n========== Performance summary ==========\n');
fprintf('%-12s | %8s | %12s | %12s | %10s | %8s\n', ...
        'Algorithm', 'Iter', 'ObjValue', 'RelError', 'Feasibility', 'Time(s)');
fprintf('%-12s-+-%8s-+-%12s-+-%12s-+-%10s-+-%8s\n', ...
        repmat('-',1,12), repmat('-',1,8), repmat('-',1,12), ...
        repmat('-',1,12), repmat('-',1,10), repmat('-',1,8));

for i = 1:numel(results)
    printTableRowWithTime(results{i}.name, results{i}, elapsedTimes(i), PhiStar);
end

plotPrefix = fullfile(resultsDir, 'ElasticNet');
Tools.plotHistories(results, Problem3, 'Elastic-net', plotPrefix, options);
fprintf('\nFigures saved to:\n');
fprintf('  %s\n', [plotPrefix '_objective.png']);
fprintf('  %s\n', [plotPrefix '_feasibility_iter.png']);

function printSolverStatus(out, PhiStar, elapsed, tolerance)
history = getHistory(out);
iterationCount = length(history.obj);
finalObjective = history.obj(end);
relativeError = relativeObjectiveError(finalObjective, PhiStar);
finalFeasibility = history.feas(end);
converged = (finalFeasibility <= tolerance) || ...
    (iterationCount >= 2 && abs(history.obj(end)-history.obj(end-1)) <= tolerance);

fprintf('  iterations : %d\n', iterationCount);
fprintf('  objective  : %.12g\n', finalObjective);
fprintf('  relative optimality error : %.2e\n', relativeError);
fprintf('  feasibility residual : %.2e\n', finalFeasibility);
printRelativeErrorStatus(history);
fprintf('  elapsed time : %.3f seconds\n', elapsed);
if converged
    fprintf('  status : converged by tolerance\n');
else
    fprintf('  status : reached max iterations\n');
end
end

function printRelativeErrorStatus(history)
if isfield(history, 'relErr2Ok') && ~isempty(history.relErr2Ok)
    satisfiedCount = nnz(history.relErr2Ok);
    totalCount = numel(history.relErr2Ok);
    ratios = history.relErr2Left ./ max(history.relErr2Right, eps);
    fprintf('  second-block relative criterion : %d/%d satisfied, max left/right = %.2e\n', ...
            satisfiedCount, totalCount, max(ratios));
end
if isfield(history, 'relErrYOk') && ~isempty(history.relErrYOk)
    satisfiedCount = nnz(history.relErrYOk);
    totalCount = numel(history.relErrYOk);
    if isfield(history, 'relErrYLeft') && isfield(history, 'relErrYRight')
        ratios = history.relErrYLeft ./ max(history.relErrYRight, eps);
        fprintf('  nonsmooth-block relative criterion : %d/%d satisfied, max left/right = %.2e\n', ...
                satisfiedCount, totalCount, max(ratios));
    else
        fprintf('  nonsmooth-block relative criterion : %d/%d satisfied\n', ...
                satisfiedCount, totalCount);
    end
end
if isfield(history, 'relErrXOk') && ~isempty(history.relErrXOk)
    satisfiedCount = nnz(history.relErrXOk);
    totalCount = numel(history.relErrXOk);
    if isfield(history, 'relErrXLeft') && isfield(history, 'relErrXRight')
        ratios = history.relErrXLeft ./ max(history.relErrXRight, eps);
        fprintf('  smooth-block relative criterion : %d/%d satisfied, max left/right = %.2e\n', ...
                satisfiedCount, totalCount, max(ratios));
    else
        fprintf('  smooth-block relative criterion : %d/%d satisfied\n', ...
                satisfiedCount, totalCount);
    end
end
end

function printTableRowWithTime(name, out, elapsed, PhiStar)
history = getHistory(out);
iterationCount = length(history.obj);
finalObjective = history.obj(end);
relativeError = relativeObjectiveError(finalObjective, PhiStar);
finalFeasibility = history.feas(end);
fprintf('%-12s | %8d | %12.5e | %10.2e | %10.2e | %8.2f\n', ...
        name, iterationCount, finalObjective, relativeError, finalFeasibility, elapsed);
end

function history = getHistory(out)
if isfield(out, 'history')
    history = out.history;
elseif isfield(out, 'hist')
    history = out.hist;
else
    error('Output structure does not contain history.');
end
end

function relativeError = relativeObjectiveError(objectiveValue, PhiStar)
relativeError = abs(objectiveValue - PhiStar) / (abs(PhiStar) + 1e-12);
end
