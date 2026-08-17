function out = tv(problem, options)
% GrpADMM.tv  求解全变分正则化问题
%    min  0.5*μ·||x - c||² + TV(x)
%   TV(x) = Σ_i ||D_i x||₂（各向同性全变分）
%   约束 y = Dx

maxIterations = options.maxIterations;
constraintResidual     = getOption(options, 'grpPsi', 1.618);
safety  = getOption(options, 'grpTauSafety', 1);

% 增广拉格朗日惩罚参数（TV 常用较大值，如 10）
sigma = 10;
% 步长需要 problem.normD2 = ||D||²
tau = safety * constraintResidual / (sigma * problem.normD2);

startTime = tic;

% --- 变量初始化 ------------------------------------------------------
x      = problem.x0;
y      = problem.y0;
lambda = problem.lambda0;
c      = problem.c;
mu     = problem.mu;
u      = x;

history = struct('obj', [], 'feas', [], 'time', [], 'err', []);

% --- 主迭代循环 ------------------------------------------------------
for k = 1:maxIterations
    % 惯性步骤
    u = ((constraintResidual - 1) / constraintResidual) * x + (1 / constraintResidual) * u;

    % x-子问题：GrpADMM 的 (7b)，S=I，未引入增广二次项。
    rightHandSide = mu * c + (1 / tau) * u - problem.DT(lambda);
    x = rightHandSide / (mu + 1 / tau);

    % y-子问题：组软阈值
    shrinkInput = problem.D(x) + lambda / sigma;
    y = Functions.groupSoftThreshold(shrinkInput, 1 / sigma);

    % 对偶变量更新
    lambda = lambda + sigma * (problem.D(x) - y);

    % 记录历史
    history = Tools.appendHistory(history, problem, x, y, k, startTime);
end

out = Tools.packOut('GrpADMM.tv', x, y, lambda, history);
end

function value = getOption(options, fieldName, defaultValue)
if isfield(options, fieldName) && ~isempty(options.(fieldName))
    value = options.(fieldName);
else
    value = defaultValue;
end
end
