function out = lad(problem, options)
% GrpADMM.lad  求解 LAD 问题（最小绝对偏差）
%   输入：
%       problem - 包含字段 .x0, .y0, .lambda0, .K, .c, .mu, .normK2
%       options    - 包含字段 .maxIterations，可选字段 .grpPsi, .grpTauSafety
%   输出：
%       out - 结构体，含 .x, .y, .lambda, .solver, .history


maxIterations = options.maxIterations;
constraintResidual     = getOption(options, 'grpPsi', 1.618);
safety  = getOption(options, 'grpTauSafety', 1);
% 增广拉格朗日惩罚参数（LAD 常用 1）
sigma = 1;
% 步长 tau：满足收敛条件 σ·τ·||K||² < ψ
tau = safety * constraintResidual / (sigma * problem.normK2);
startTime = tic;
% --- 变量初始化 ------------------------------------------------------
x      = problem.x0;
y      = problem.y0;
lambda = problem.lambda0;
K      = problem.K;
c      = problem.c;
mu     = problem.mu;
u      = x;                         % 惯性变量 u0 = x0

% --- 历史记录 --------------------------------------------------------
history = struct('obj', [], 'feas', [], 'time', [], 'err', []);

% --- 主迭代循环 ------------------------------------------------------
for k = 1:maxIterations
    % 1) 惯性步骤 (7a)
    u = ((constraintResidual - 1) / constraintResidual) * x + (1 / constraintResidual) * u;

    % 2) x-子问题：软阈值
    x = Functions.softThreshold(u - tau * (K' * lambda), tau * mu);

    % 3) y-子问题：软阈值
    y = c + Functions.softThreshold(K * x + lambda / sigma - c, 1 / sigma);

    % 4) 对偶变量更新 (7d)
    lambda = lambda + sigma * (K * x - y);

    % 5) 记录历史
    history = Tools.appendHistory(history, problem, x, y, k, startTime);
end

% --- 打包输出 --------------------------------------------------------
out = Tools.packOut('GrpADMM.lad', x, y, lambda, history);
end

function value = getOption(options, fieldName, defaultValue)
if isfield(options, fieldName) && ~isempty(options.(fieldName))
    value = options.(fieldName);
else
    value = defaultValue;
end
end
