function out = lad(problem, options)
% APADMC.lad  加速 PADMC 算法求解 LAD 问题（强凸变体）
%   输入：
%       problem - 包含字段 .x0, .y0, .lambda0, .K, .c, .mu1, .mu2, .normH2
%       options    - 包含字段 .maxIterations, .tolerance
%   输出：
%       out - 结构体，含 .x, .y, .lambda, .solver, .history

maxIterations = options.maxIterations;
tolerance     = options.tolerance;
startTime  = tic;

if ~isfield(problem, 'mu2') || problem.mu2 <= 0
    error('APADMC:NotStronglyConvex', ...
        ['APADMC.lad requires a strongly convex LAD variant with ', ...
         'problem.mu2 > 0. The standard LAD problem is not covered.']);
end

% 设置问题变体标记（原代码用于历史记录）
problem.variant = 'strong';

% 算法固定参数
alphaX = 2;
alphaW = 1;
rxBase = 1.01 * problem.normH2;
sw     = 1;

% 变量初始化
x      = problem.x0;
y      = problem.y0;
lambda = problem.lambda0;
xExtrapolated     = x;
yExtrapolated     = y;

history = struct('obj', [], 'feas', [], 'time', [], 'err', []);

for k = 1:maxIterations
    % 迭代依赖的惩罚参数 beta
    beta = 2 * problem.mu2 * k / (3 * problem.normH2);
    xProximalCoefficient = rxBase / beta;

    % 加速变量更新
    xExtrapolated = alphaX * x + (1 - alphaX) * xExtrapolated;
    yExtrapolated = alphaW * y + (1 - alphaW) * yExtrapolated;

    % x-子问题
    xCenter = xExtrapolated - (problem.K' * (lambda + beta * (problem.K * xExtrapolated - y))) / xProximalCoefficient;
    x = Functions.softThreshold((xProximalCoefficient / (xProximalCoefficient + 2 * problem.mu2)) * xCenter, ...
                       problem.mu1 / (xProximalCoefficient + 2 * problem.mu2));

    % y-子问题
    yCoefficient = beta + sw / beta;
    yCenter = (beta * problem.K * x + lambda + (sw / beta) * yExtrapolated) / yCoefficient;
    y = problem.c + Functions.softThreshold(yCenter - problem.c, 1 / yCoefficient);

    % 对偶变量更新
    lambda = lambda + beta * (problem.K * x - y);

    % 记录历史并检查停机
    history = Tools.appendHistory(history, problem, x, y, k, startTime);
    if shouldStop(history, k, tolerance)
        break;
    end
end

out = Tools.packOut('APADMC.lad', x, y, lambda, history);
end

function tf = shouldStop(history, k, tolerance)
tf = k >= 1 && isfield(history, 'feas') && history.feas(k) <= tolerance;
end
