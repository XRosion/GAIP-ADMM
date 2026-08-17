function out = elastic(problem, options)
% APADMC.elastic  Accelerated PADMC for the elastic-net x-y split.
%
%   min  gamma*||x||^2 + 0.5*||K*x-b||^2 + eta*||y||_1,
%   s.t. x-y=0.

maxIterations = options.maxIterations;
tolerance = options.tolerance;
keepOnlyVisualizationData = Tools.keepOnlyVisualizationData(options);
startTime = tic;

alphaX = getOption(options, 'apadmcAlphaX', 2);
alphaY = getOption(options, 'apadmcAlphaY', 1);
xProximalWeight = getOption(options, 'apadmcXProximalWeight', 100);
yProximalWeight = getOption(options, 'apadmcYProximalWeight', 50);
smoothStrongModulus = problem.preStrongModulusG;
ridgeCurvature = problem.ridgeCurvature;

K = problem.K;
b = problem.b;
eta = problem.eta;
n = size(K, 2);

x = problem.x0;
y = compatibleVector(problem, 'y0', x, n);
lambda = compatibleVector(problem, 'lambda0', zeros(n, 1), n);
xExtrapolated = x;
yExtrapolated = y;
normalRightBase = K' * b;

if keepOnlyVisualizationData
    history = struct('obj', [], 'feas', []);
else
    history = struct('obj', [], 'feas', [], 'time', [], 'err', [], 'inner', []);
end

for k = 1:maxIterations
    beta = smoothStrongModulus * k / 3;

    xExtrapolated = alphaX * x + (1 - alphaX) * xExtrapolated;
    yExtrapolated = alphaY * y + (1 - alphaY) * yExtrapolated;

    xProximalCoefficient = xProximalWeight / beta;
    xCoefficient = beta + xProximalCoefficient + ridgeCurvature;
    xRightHandSide = normalRightBase + beta * y + lambda + ...
        xProximalCoefficient * xExtrapolated;
    woodburyMatrix = eye(size(K, 1)) + (K * K') / xCoefficient;
    [woodburyFactor, woodburyFlag] = chol(woodburyMatrix);
    projectedRightHandSide = K * xRightHandSide;
    if woodburyFlag == 0
        correction = woodburyFactor \ (woodburyFactor' \ projectedRightHandSide);
    else
        correction = woodburyMatrix \ projectedRightHandSide;
    end
    x = xRightHandSide / xCoefficient - K' * correction / (xCoefficient^2);

    yProximalCoefficient = yProximalWeight / beta;
    yCoefficient = beta + yProximalCoefficient;
    yCenter = (beta * x - lambda + yProximalCoefficient * yExtrapolated) / yCoefficient;
    y = Functions.softThreshold(yCenter, eta / yCoefficient);

    lambda = lambda - beta * (x - y);

    history = Tools.appendHistory(history, problem, x, y, k, startTime, options);
    if ~keepOnlyVisualizationData
        history.inner(k, 1) = 0;
    end
    if shouldStop(history, k, tolerance)
        break;
    end
end

out = Tools.packOut('APADMC.elastic', x, y, lambda, history, options);
end

function tf = shouldStop(history, k, tolerance)
tf = k >= 1 && isfield(history, 'feas') && numel(history.feas) >= k && ...
    history.feas(k) <= tolerance;
end

function value = compatibleVector(problem, fieldName, defaultValue, n)
if isfield(problem, fieldName) && numel(problem.(fieldName)) == n
    value = problem.(fieldName);
else
    value = defaultValue;
end
end

function value = getOption(options, fieldName, defaultValue)
if isfield(options, fieldName) && ~isempty(options.(fieldName))
    value = options.(fieldName);
else
    value = defaultValue;
end
end
