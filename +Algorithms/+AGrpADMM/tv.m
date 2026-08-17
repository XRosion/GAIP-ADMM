function out = tv(problem, options)
% AGrpADMM.tv  加速黄金比例ADMM求解全变分去噪问题
%   原始问题：min_x 0.5*μ‖x - c‖² + TV(x), TV(x)=Σ_i ‖D_i x‖₂
%   拆分形式：y = Dx

maxIterations = options.maxIterations;
tolerance     = options.tolerance;
constraintResidual     = getOption(options, 'grpPsi', 1.618);
startTime  = tic;

% 算法参数（原硬编码值）
rho = 1;
xStrongModulus = problem.preStrongModulusG;
xProximalWeight  = xStrongModulus;
tw  = 0.1;
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

    % x-子问题（求解线性系统，假设 DᵀD 可对角化）
    rightHandSide = problem.mu * problem.c + (xProximalWeight / tau) * u - problem.DT(lambda);
    x = rightHandSide / (problem.mu + xProximalWeight / tau);

    % y-子问题（组软阈值，时变近端）
    yCoefficient = sigma + tw / t;
    yCenter = (sigma * problem.D(x) + lambda + (tw / t) * y) / yCoefficient;
    y = Functions.groupSoftThreshold(yCenter, 1 / yCoefficient);

    % 对偶变量更新
    lambda = lambda + sigma * (problem.D(x) - y);

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

out = Tools.packOut('AGrpADMM.tv', x, y, lambda, history);
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
