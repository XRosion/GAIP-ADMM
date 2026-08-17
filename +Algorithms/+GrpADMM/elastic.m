function out = elastic(problem, options)
% GrpADMM.elastic  Golden-ratio ADMM for the elastic-net x-y split.
%
%   min  gamma*||x||^2 + 0.5*||K*x-b||^2 + eta*||y||_1,
%   s.t. x-y=0.
%
% The project convention is that x is the smooth quadratic block and y is
% the nonsmooth l1 block. This matches objectiveValue/feasibilityNorm.

maxIterations = options.maxIterations;
psi = getOption(options, 'grpPsi', 1.618);
sigma = getOption(options, 'grpSigma', 1);
keepOnlyVisualizationData = Tools.keepOnlyVisualizationData(options);
startTime = tic;

K = problem.K;
b = problem.b;
eta = problem.eta;
ridgeCurvature = problem.ridgeCurvature;
n = size(K, 2);

x = problem.x0;
y = compatibleVector(problem, 'y0', x, n);
lambda = compatibleVector(problem, 'lambda0', zeros(n, 1), n);
yGolden = y;

xCoefficient = sigma + ridgeCurvature;
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
    yGolden = ((psi - 1) / psi) * y + (1 / psi) * yGolden;

    xRightHandSide = normalRightBase + sigma * yGolden + lambda;
    projectedRightHandSide = K * xRightHandSide;
    if useWoodburyCholesky
        correction = woodburyFactor \ (woodburyFactor' \ projectedRightHandSide);
    else
        correction = woodburyMatrix \ projectedRightHandSide;
    end
    x = xRightHandSide / xCoefficient - K' * correction / (xCoefficient^2);

    y = Functions.softThreshold(x - lambda / sigma, eta / sigma);
    lambda = lambda - sigma * (x - y);

    history = Tools.appendHistory(history, problem, x, y, k, startTime, options);
    if ~keepOnlyVisualizationData
        history.inner(k, 1) = 0;
    end
end

out = Tools.packOut('GrpADMM.elastic', x, y, lambda, history, options);
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
