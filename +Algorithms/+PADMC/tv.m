function out = tv(problem, options)
% PADMC.tv  求解全变分去噪问题的 PADMC 算法
%   原始问题：min_x 0.5*μ‖x - c‖² + TV(x), TV(x)=Σ_i ‖D_i x‖₂
%   等价拆分：min_{x,y} μ‖x - c‖² + Σ_i ‖w_i‖₂, s.t. Dx = y

maxIterations = options.maxIterations;
startTime = tic;

% PADMC 特定参数（TV 分支原硬编码值）
beta   = 10;                         % 增广拉格朗日惩罚参数
alphaX = 2;                          % x 的外推系数
alphaW = 1;                          % y 的外推系数
xProximalCoefficient     = 1.01 * beta * problem.normD2;   % x-子问题的近端步长相关参数
sw     = 0;                          % y-子问题的额外近端正则化系数

% 变量初始化
x      = problem.x0;
y      = problem.y0;
lambda = problem.lambda0;
xExtrapolated     = x;
yExtrapolated     = y;

history = struct('obj', [], 'feas', [], 'time', [], 'err', []);

for k = 1:maxIterations
    % 外推中心更新
    xExtrapolated = alphaX * x + (1 - alphaX) * xExtrapolated;
    yExtrapolated = alphaW * y + (1 - alphaW) * yExtrapolated;

    % x-子问题：求解线性系统（使用傅里叶域近似解法）
    gradAtXb = problem.DT(lambda + beta * (problem.D(xExtrapolated) - y));
    x = (problem.mu * problem.c + xProximalCoefficient * xExtrapolated - gradAtXb) / (problem.mu + xProximalCoefficient);

    % y-子问题：组软阈值（对每个梯度向量联合收缩）
    yCoefficient = beta + sw;
    yCenter = (beta * problem.D(x) + lambda + sw * yExtrapolated) / yCoefficient;
    y = Functions.groupSoftThreshold(yCenter, 1 / yCoefficient);

    % 对偶变量更新
    lambda = lambda + beta * (problem.D(x) - y);

    % 记录历史
    history = Tools.appendHistory(history, problem, x, y, k, startTime);
end

out = Tools.packOut('PADMC.tv', x, y, lambda, history);
end
