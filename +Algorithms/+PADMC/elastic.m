function out = elastic(problem, options)
% PADMC.elastic  PADMC for the elastic-net x-y split.
%
%   min  gamma*||x||^2 + 0.5*||K*x-b||^2 + eta*||y||_1,
%   s.t. x-y=0.
%
% The x-block is the smooth elastic-net quadratic; the y-block is solved by
% soft-thresholding.

maxIterations = options.maxIterations;
keepOnlyVisualizationData = Tools.keepOnlyVisualizationData(options);
startTime = tic;

beta = getOption(options, 'padmcBeta', 0.005);
alphaX = getOption(options, 'padmcAlphaX', 2);
alphaY = getOption(options, 'padmcAlphaY', 1);
xProximalWeight = getOption(options, 'padmcXProximalWeight', 100);
yProximalWeight = getOption(options, 'padmcYProximalWeight', 50);

K = problem.K;
b = problem.b;
eta = problem.eta;
ridgeCurvature = problem.ridgeCurvature;
n = size(K, 2);

x = problem.x0;
y = compatibleVector(problem, 'y0', x, n);
lambda = compatibleVector(problem, 'lambda0', zeros(n, 1), n);
xExtrapolated = x;
yExtrapolated = y;

xCoefficient = beta + xProximalWeight + ridgeCurvature;
woodburyMatrix = eye(size(K, 1)) + (K * K') / xCoefficient;
[woodburyFactor, woodburyFlag] = chol(woodburyMatrix);
useWoodburyCholesky = woodburyFlag == 0;
normalRightBase = K' * b;

if keepOnlyVisualizationData
    history = struct('obj', [], 'feas', []);
else
    history = struct('obj', [], 'feas', [], 'time', [], 'err', [], 'inner', []);
end

for k = 1:maxIterations
    xExtrapolated = alphaX * x + (1 - alphaX) * xExtrapolated;
    yExtrapolated = alphaY * y + (1 - alphaY) * yExtrapolated;

    xRightHandSide = normalRightBase + beta * y + lambda + ...
        xProximalWeight * xExtrapolated;
    projectedRightHandSide = K * xRightHandSide;
    if useWoodburyCholesky
        correction = woodburyFactor \ (woodburyFactor' \ projectedRightHandSide);
    else
        correction = woodburyMatrix \ projectedRightHandSide;
    end
    x = xRightHandSide / xCoefficient - K' * correction / (xCoefficient^2);

    yCoefficient = beta + yProximalWeight;
    yCenter = (beta * x - lambda + yProximalWeight * yExtrapolated) / yCoefficient;
    y = Functions.softThreshold(yCenter, eta / yCoefficient);

    lambda = lambda - beta * (x - y);

    history = Tools.appendHistory(history, problem, x, y, k, startTime, options);
    if ~keepOnlyVisualizationData
        history.inner(k, 1) = 0;
    end
end

out = Tools.packOut('PADMC.elastic', x, y, lambda, history, options);
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
