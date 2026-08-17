% Main4: TV image recovery numerical experiment main script.
% Compares ADMM variants on
%
%   min_x 0.5*mu*||x-c||_2^2 + TV(x),
%   using split: D*x - y = 0.
%
% This script follows the same lightweight framework as Main3.m:
% experiment options, algorithm registry, unified calling loop, summary,
% and visualization are kept in one file with local helper functions.

clc; clear;
import Algorithms.*
import Functions.*
import Problems.*
import Tools.*

rootDir = fileparts(mfilename('fullpath'));
resultsDir = fullfile(rootDir, 'TVResults');
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end


%% 1. Experiment options
imageNames = {'cameraman','man'};
% imageNames = {'cameraman'};

options.maxIterations = 500;
options.tolerance = 1e-6;
options.imageSize = [256, 256];
options.noiseStd = 0.05;
options.mu = 30;
options.seed = 20260807;

% Linear-solve controls for TV smooth image subproblems. GAIP/AGAIP update
% this subproblem second; comparison algorithms retain their native order.
options.tvCgMaxIter = 200;
options.tvCgAbsTol = 1e-6;

% Plot options.
options.plotLegendLocation = 'northeast';
options.plotUseMarkers = false;
options.plotShowTitle = false;
options.plotObjectiveEnhanced = false;
options.plotObjectiveRawMarkers = false;
options.plotObjectiveMinPoint = false;
options.plotObjectiveSmoothingWindow = 101;     % 目标函数相对误差图的 log10 移动平均窗口。

%% 2. Algorithm registry
% Proposed GAIPADMM parameters. Only GAIPADMM and AGAIPADMM use this
% third input. Comparison algorithms below are called as Func(problem, opts).
% 耦合参数
Paras.alpha = 1.618;
Paras.theta_lambda = 0.44;
Paras.beta = 80;
Paras.theta_x = -0.618;
Paras.theta_y = -0.618;

% Block parameters follow Algorithm 1: first=x (TV auxiliary),
% second=y (image).
Paras.tau1 = 0.01;
Paras.tau2 = 10;
Paras.eps2 = 0.99;
Paras.eps3 = 0;

AccParas = Paras;
AccParas.beta = 30;
AccParas.theta_gamma = 1;

Algorithms = struct('Name', {}, 'Func', {}, 'Type', {}, ...
    'Parameters', {}, 'Description', {}, 'Enabled', {});
algorithmIndex = 0;

algorithmIndex = algorithmIndex + 1;
Algorithms(algorithmIndex).Name = 'AGAIP-ADMM_{spp}';
Algorithms(algorithmIndex).Func = @Algorithms.AGAIPADMM.tv;
Algorithms(algorithmIndex).Type = 1;
Algorithms(algorithmIndex).Parameters = AccParas;
Algorithms(algorithmIndex).Description = 'AGAIP-ADMM_spp with exact TV first block and inexact image second block';
Algorithms(algorithmIndex).Enabled = true;

algorithmIndex = algorithmIndex + 1;
Algorithms(algorithmIndex).Name = 'GAIP-ADMM_{spp}';
Algorithms(algorithmIndex).Func = @Algorithms.GAIPADMM.tv;
Algorithms(algorithmIndex).Type = 1;
Algorithms(algorithmIndex).Parameters = Paras;
Algorithms(algorithmIndex).Description = 'GAIP-ADMM_spp with exact TV first block and inexact image second block';
Algorithms(algorithmIndex).Enabled = true;

algorithmIndex = algorithmIndex + 1;
Algorithms(algorithmIndex).Name = 'IADMM';
Algorithms(algorithmIndex).Func = @Algorithms.IADMM.tv;
Algorithms(algorithmIndex).Type = 0;
Algorithms(algorithmIndex).Description = 'Inexact ADMM with relative-error CG x-step';
Algorithms(algorithmIndex).Enabled = true;

algorithmIndex = algorithmIndex + 1;
Algorithms(algorithmIndex).Name = 'IGADMM';
Algorithms(algorithmIndex).Func = @Algorithms.IGADMM.tv;
Algorithms(algorithmIndex).Type = 0;
Algorithms(algorithmIndex).Description = 'Generalized ADMM with y proximal term and inexact TV x-step';
Algorithms(algorithmIndex).Enabled = true;

Algorithms = Algorithms([Algorithms.Enabled]);

%% 3. Run experiments image by image
allResults = struct();

for iImg = 1:numel(imageNames)
    imageName = imageNames{iImg};
    fprintf('\n========== TV problem generation (%d/%d) ==========\n', ...
        iImg, numel(imageNames));

    problem = buildTvProblem(imageName, options, iImg);
    noisyObj = Functions.objectiveValue(problem, problem.c, problem.D(problem.c));
    fprintf('Image: %s, size = %d x %d\n', problem.imageName, problem.m, problem.n);
    fprintf('Regularization: mu = %.4g, noiseStd = %.4g\n', ...
        problem.mu, problem.noiseStd);
    fprintf('Noisy objective = %.12g, noisy PSNR = %.3f dB\n', ...
        noisyObj, Functions.psnrValue(problem, problem.c));
    if isfield(problem, 'cvxInfo')
        fprintf('Reference solve: solver=%s, status=%s\n', ...
            problem.cvxInfo.solver, problem.cvxInfo.status);
    end

    fprintf('\n========== Start solving: %s ==========\n', problem.imageName);
    [results, elapsedTimes] = Tools.runRegisteredAlgorithms( ...
        Algorithms, problem, options, ...
        @(out, elapsed, alg) printSolverStatus(out, problem, elapsed, options.tolerance));

    printPerformanceSummary(results, elapsedTimes, problem);
    saveTvFigures(results, problem, resultsDir, options);

    key = Tools.fileSafeName(problem.imageName);
    allResults.(key).problem = problem;
    allResults.(key).results = results;
    allResults.(key).elapsedTimes = elapsedTimes;
end

%% 4. Save summary
experiment.imageNames = imageNames;
experiment.options = options;
save(fullfile(resultsDir, 'tv_experiment_summary.mat'), ...
    'allResults', 'experiment', 'Algorithms');

fprintf('\nAll TV image recovery experiments completed. Results saved in:\n  %s\n', ...
    resultsDir);

%% ======================= Local helper functions =======================
function problem = buildTvProblem(imageName, options, imageIndex)
probOptions = struct();
probOptions.imageSize = options.imageSize;
probOptions.noiseStd = options.noiseStd;
probOptions.mu = options.mu;
probOptions.seed = options.seed + imageIndex - 1;
problem = Problems.tvProblem(imageName, probOptions);
end

function printSolverStatus(out, problem, elapsed, tolerance)
history = out.history;
iterationCount = numel(history.obj);
finalObjective = history.objPrimal(end);
finalFeasibility = history.feas(end);
finalStopMetric = history.outerStopMetric(end);
converged = finalStopMetric <= tolerance;

fprintf('  iterations: %d\n', iterationCount);
fprintf('  objective (primal): %.12g\n', finalObjective);
fprintf('  relative optimality error: %.2e\n', history.err(end));
fprintf('  feasibility residual: %.2e\n', finalFeasibility);
fprintf('  outer relative stop metric: %.2e\n', finalStopMetric);
fprintf('  PSNR: %.3f dB\n', Functions.psnrValue(problem, out.x));
printRelativeErrorStatus(history);
fprintf('  elapsed time: %.3f seconds\n', elapsed);
if converged
    fprintf('  status: converged by tolerance\n');
else
    fprintf('  status: reached max iterations\n');
end
end

function printRelativeErrorStatus(history)
if isfield(history, 'relErr2Ok') && ~isempty(history.relErr2Ok)
    satisfiedCount = nnz(history.relErr2Ok);
    totalCount = numel(history.relErr2Ok);
    ratios = history.relErr2Left ./ max(history.relErr2Right, eps);
    fprintf('  second-block relative criterion: %d/%d satisfied, max left/right = %.2e\n', ...
        satisfiedCount, totalCount, max(ratios));
end

if isfield(history, 'relErrYOk') && ~isempty(history.relErrYOk)
    satisfiedCount = nnz(history.relErrYOk);
    totalCount = numel(history.relErrYOk);
    ratios = history.relErrYLeft ./ max(history.relErrYRight, eps);
    fprintf('  y-block relative criterion: %d/%d satisfied, max left/right = %.2e\n', ...
        satisfiedCount, totalCount, max(ratios));
end

if isfield(history, 'relErrXOk') && ~isempty(history.relErrXOk)
    satisfiedCount = nnz(history.relErrXOk);
    totalCount = numel(history.relErrXOk);
    ratios = history.relErrXLeft ./ max(history.relErrXRight, eps);
    fprintf('  x-block relative criterion: %d/%d satisfied, max left/right = %.2e\n', ...
        satisfiedCount, totalCount, max(ratios));
end

if isfield(history, 'absStop2')
    fprintf('  second-block linear residual tolerance: %d/%d satisfied\n', ...
        nnz(history.absStop2), numel(history.absStop2));
elseif isfield(history, 'xAbsStop')
    fprintf('  x-block linear residual tolerance: %d/%d satisfied\n', ...
        nnz(history.xAbsStop), numel(history.xAbsStop));
end

fprintf('  cumulative CG iterations: %d\n', history.innerTotal(end));

if isfield(history, 'linearResidual2')
    fprintf('  final normalized second-block linear residual: %.2e\n', ...
        history.linearResidual2(end));
elseif isfield(history, 'xLinearResidual')
    fprintf('  final normalized x linear residual: %.2e\n', ...
        history.xLinearResidual(end));
end
end

function printPerformanceSummary(results, elapsedTimes, problem)
fprintf('\n========== Performance summary: %s ==========\n', problem.imageName);
fprintf('%-20s | %8s | %12s | %10s | %10s | %9s | %8s\n', ...
    'Algorithm', 'Iter', 'ObjValue', 'Feas', 'RelObjErr', 'PSNR(dB)', 'Time(s)');
fprintf('%s\n', repmat('-', 1, 92));

for iAlg = 1:numel(results)
    out = results{iAlg};
    history = out.history;
    iterationCount = numel(history.obj);
    primalObjective = history.objPrimal(end);
    relErr = abs(history.objPrimal(end) - problem.PhiStar) / max(abs(problem.PhiStar), eps);
    fprintf('%-20s | %8d | %12.5e | %10.2e | %10.2e | %9.3f | %8.2f\n', ...
        out.name, iterationCount, primalObjective, history.feas(end), ...
        relErr, ...
        Functions.psnrValue(problem, out.x), elapsedTimes(iAlg));
end
end

function saveTvFigures(results, problem, resultsDir, options)
plotPrefix = fullfile(resultsDir, Tools.fileSafeName(problem.imageName));
Tools.plotHistories(results, problem, ['TV ' problem.imageName], plotPrefix, options);
Tools.plotTvRecovery(problem, results, plotPrefix);

fprintf('\nFigures saved to:\n');
if strcmpi(problem.type, 'TV')
    fprintf('  %s\n', [plotPrefix '_objective.png']);
    fprintf('  %s\n', [plotPrefix '_objective_time.png']);
else
    fprintf('  %s\n', [plotPrefix '_objective.png']);
end
fprintf('  %s\n', [plotPrefix '_feasibility_iter.png']);
fprintf('  %s\n', [plotPrefix '_psnr_iter.png']);
fprintf('  %s\n', [plotPrefix '_psnr_time.png']);
fprintf('  %s\n', [plotPrefix '_recovery_comparison.png']);
end
