function out = lad(problem, options)
% AGrpADMM.lad  加速黄金比例ADMM求解LAD问题（强凸变体）
%   输入：
%       problem - 包含字段 .x0, .y0, .lambda0, .K, .c, .mu1, .mu2, .normH2
%       options    - 包含字段 .maxIterations, .tolerance
%   输出：
%       out - 结构体，含 .x, .y, .lambda, .solver, .history

maxIterations = options.maxIterations;
tolerance     = options.tolerance;
constraintResidual     = getOption(options, 'grpPsi', 1.618);
startTime  = tic;

if ~isfield(problem, 'mu2') || problem.mu2 <= 0
    error('AGrpADMM:NotStronglyConvex', ...
        ['AGrpADMM.lad requires the strongly convex LAD model with ', ...
         'problem.mu2 > 0, as in GrpADMM paper equation (43).']);
end

% 标记为强凸LAD变体（原代码调用 make_lad_strong）
problem.variant = 'strong';

% 算法参数（原硬编码值）
rho = 0.001;
xProximalWeight  = problem.mu2;      % 强凸参数
tw  = 0.01;
sigma = rho;
tau   = constraintResidual * xProximalWeight / (sigma * problem.normH2);
t     = 1;

% 变量初始化
x      = problem.x0;
y      = problem.y0;
lambda = problem.lambda0;
u      = x;

history = struct('obj', [], 'feas', [], 'time', [], 'err', []);

for k = 1:maxIterations
    % 惯性步骤
    u = ((constraintResidual - 1) / constraintResidual) * x + (1 / constraintResidual) * u;

    % x-子问题（带强凸近端）
    xProximalCoefficient = xProximalWeight / tau;
    xCenter = u - (problem.K' * lambda) / xProximalCoefficient;
    x = Functions.softThreshold((xProximalCoefficient / (xProximalCoefficient + 2 * problem.mu2)) * xCenter, ...
                       problem.mu1 / (xProximalCoefficient + 2 * problem.mu2));

    % y-子问题（时变近端 + 软阈值）
    yCoefficient = sigma + tw / t;
    yCenter = (sigma * problem.K * x + lambda + (tw / t) * y) / yCoefficient;
    y = problem.c + Functions.softThreshold(yCenter - problem.c, 1 / yCoefficient);

    % 对偶变量更新
    lambda = lambda + sigma * (problem.K * x - y);

    % 记录历史
    history = Tools.appendHistory(history, problem, x, y, k, startTime);
    if shouldStop(history, k, tolerance)
        break;
    end

    % 更新加速参数
    tNew = 0.5 * (1 + sqrt(1 + 4 * t^2));
    t = tNew;
    sigma = rho * t;
    tau   = constraintResidual * xProximalWeight / (sigma * problem.normH2);
end

out = Tools.packOut('AGrpADMM.lad', x, y, lambda, history);
end

function tf = shouldStop(history, k, tolerance)
tf = k >= 1 && isfield(history, 'feas') && history.feas(k) <= tolerance;
end

function value = getOption(options, fieldName, defaultValue)
if isfield(options, fieldName) && ~isempty(options.(fieldName))
    value = options.(fieldName);
else
    value = defaultValue;
end
end
