function out = tv(problem, options)
% APADMC.tv  加速 PADMC 算法求解全变分去噪问题
%   原始问题：min_x 0.5*μ‖x - c‖² + TV(x), TV(x)=Σ_i ‖D_i x‖₂
%   拆分形式：y = Dx

maxIterations = options.maxIterations;
tolerance     = options.tolerance;
startTime  = tic;

% 算法固定参数
alphaX = 1;
alphaW = 2;
sw     = 1;                % y 子问题的平滑参数
xProximalWeight     = problem.mu;       % x 子问题的近端权重
muStrong = problem.preStrongModulusG;

% 变量初始化
x      = problem.x0;
y      = problem.y0;
lambda = problem.lambda0;
xExtrapolated     = x;
yExtrapolated     = y;

history = struct('obj', [], 'feas', [], 'time', [], 'err', []);

for k = 1:maxIterations
    % 迭代依赖的惩罚参数 beta
    beta = muStrong * k / (3 * problem.normH2);

    % 加速变量更新（注意顺序：先 yExtrapolated，后 xExtrapolated）
    yExtrapolated = alphaW * y + (1 - alphaW) * yExtrapolated;
    xExtrapolated = alphaX * x + (1 - alphaX) * xExtrapolated;

    % y-子问题（组软阈值）
    yCoefficient = beta + sw / beta;
    yCenter = (beta * problem.D(x) + lambda + (sw / beta) * yExtrapolated) / yCoefficient;
    y = Functions.groupSoftThreshold(yCenter, 1 / yCoefficient);

    % x-子问题（FFT 求解线性系统）
    rightHandSide = problem.mu * problem.c - problem.DT(lambda) + beta * problem.DT(y) + (xProximalWeight / beta) * xExtrapolated;
    xCoefficient = problem.mu + xProximalWeight / beta;
    rightHandSideImage = reshape(rightHandSide, problem.m, problem.n);
    denominator = xCoefficient + beta * problem.eigsDtD;
    xImage = real(ifft2(fft2(rightHandSideImage) ./ denominator));
    x = xImage(:);

    % 对偶变量更新
    lambda = lambda + beta * (problem.D(x) - y);

    % 记录历史并检查停机
    history = Tools.appendHistory(history, problem, x, y, k, startTime);
    if shouldStop(history, k, tolerance)
        break;
    end
end

out = Tools.packOut('APADMC.tv', x, y, lambda, history);
end

function tf = shouldStop(history, k, tolerance)
tf = k >= 1 && isfield(history, 'feas') && history.feas(k) <= tolerance;
end
