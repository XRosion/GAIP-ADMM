function out = lasso(problem, options)
% IGADMM.lasso  用 x-y=0 分裂形式求解 LASSO 问题。
%
%   min_x  eta*||x||_1 + 0.5*||K*x-b||^2
%
% 采用分裂形式：
%
%   min_{x,y}  0.5*||K*x-b||^2 + eta*||y||_1,  s.t. x-y=0.
%
% x 子问题为光滑二次问题，由 CG 按 IGADMM 相对误差准则非精确求解；
% y 子问题为 l1 近端问题，由软阈值算子精确求解。

maxIterations = options.maxIterations;
beta = 1;
alpha = -0.3;
sigma1 = 0.99;
yProximalParameter = 0.01;
if ~(beta > 0)
    error('IGADMM:InvalidBeta', 'beta must be positive.');
end
if ~(alpha > -1 && alpha < 1)
    error('IGADMM:InvalidAlpha', 'alpha must be in (-1, 1).');
end
if ~(sigma1 >= 0 && sigma1 < 1)
    error('IGADMM:InvalidSigma1', 'sigma1 must be in [0, 1).');
end
if ~(yProximalParameter >= 0)
    error('IGADMM:InvalidYProximalParameter', ...
        'yProximalParameter must be nonnegative.');
end
maxInnerIterations = options.lassoCgMaxIter;
xAbsTol = options.lassoCgAbsTol;
keepOnlyVisualizationData = Tools.keepOnlyVisualizationData(options);
startTime = tic;

K = problem.K;
b = problem.b;
eta = problem.eta;
n = size(K, 2);

x = problem.x0;
y = problem.y0;
lambda = problem.lambda0;
auxiliaryX = x;
normalRightBase = K' * b;

if keepOnlyVisualizationData
    history = struct('obj', [], 'feas', []);
else
    history = struct('obj', [], 'feas', [], 'time', [], 'inner', [], 'err', [], ...
                  'relErrXLeft', [], 'relErrXRight', [], 'relErrXOk', [], ...
                  'relErrYLeft', [], 'relErrYRight', [], 'relErrYOk', [], ...
                  'innerResNorm', []);
end

for iteration = 1:maxIterations
    xRightHandSide = normalRightBase + beta * y + lambda;
    xCandidate = x;
    d = K' * (K * xCandidate) + beta * xCandidate - xRightHandSide;
    relLeft = 2 * abs((auxiliaryX - xCandidate)' * d) + d' * d;
    relRight = (1 - alpha) * sigma1 * beta * norm(xCandidate - y)^2;
    relOk = relLeft <= relRight;
    residual = -d;
    direction = residual;
    residualSquared = residual' * residual;
    innerIteration = 0;

    while innerIteration < maxInnerIterations && ...
            ~(relOk || sqrt(residualSquared) <= xAbsTol)
        matrixDirection = K' * (K * direction) + beta * direction;
        denominator = direction' * matrixDirection;
        if denominator <= eps
            break;
        end

        stepSize = residualSquared / denominator;
        xCandidate = xCandidate + stepSize * direction;
        residual = residual - stepSize * matrixDirection;
        d = -residual;
        newResidualSquared = residual' * residual;
        innerIteration = innerIteration + 1;

        relLeft = 2 * abs((auxiliaryX - xCandidate)' * d) + d' * d;
        relRight = (1 - alpha) * sigma1 * beta * norm(xCandidate - y)^2;
        relOk = relLeft <= relRight;
        if relOk || sqrt(newResidualSquared) <= xAbsTol
            residualSquared = newResidualSquared;
            break;
        end

        direction = residual + (newResidualSquared / residualSquared) * direction;
        residualSquared = newResidualSquared;
    end
    relOk = relOk || sqrt(residualSquared) <= xAbsTol;
    xNew = xCandidate;

    constraintResidual = xNew - y;
    lambdaHalf = lambda - alpha * beta * constraintResidual;

    yCoefficient = beta + yProximalParameter;
    yCenter = (beta * xNew - lambdaHalf + yProximalParameter * y) / yCoefficient;
    yNew = Functions.softThreshold(yCenter, eta / yCoefficient);

    newConstraintResidual = xNew - yNew;
    lambdaNew = lambdaHalf - beta * newConstraintResidual;
    auxiliaryX = auxiliaryX - d;

    history = Tools.appendHistory(history, problem, xNew, yNew, iteration, startTime, options);
    if ~keepOnlyVisualizationData
        history.inner(iteration, 1) = innerIteration;
        history.relErrXLeft(iteration, 1) = relLeft;
        history.relErrXRight(iteration, 1) = relRight;
        history.relErrXOk(iteration, 1) = relOk;
        history.relErrYLeft(iteration, 1) = 0;
        history.relErrYRight(iteration, 1) = 0;
        history.relErrYOk(iteration, 1) = true;
        history.innerResNorm(iteration, 1) = norm([d; zeros(size(yNew))]);
    end

    x = xNew;
    y = yNew;
    lambda = lambdaNew;
end

out = Tools.packOut('IGADMM.lasso', x, y, lambda, history, options);
end
