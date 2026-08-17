function out = elastic(problem, options)
% AGrpADMM.elastic  Accelerated golden-ratio ADMM for elastic net.
%
%   min  gamma*||x||^2 + 0.5*||K*x-b||^2 + eta*||y||_1,
%   s.t. x-y=0.
%
% The accelerated variable is the strongly convex smooth x-block.

maxIterations = options.maxIterations;
tolerance = options.tolerance;
psi = getOption(options, 'grpPsi', 1.618);
keepOnlyVisualizationData = Tools.keepOnlyVisualizationData(options);
startTime = tic;

rho = getOption(options, 'agrpRho', 0.01);
xStrongModulus = problem.preStrongModulusG;
ridgeCurvature = problem.ridgeCurvature;
yProximalWeight = getOption(options, 'agrpYProximalWeight', 0.01);
sigma = rho;
tau = psi * xStrongModulus / (sigma * problem.constraintNorm2);
t = 1;

K = problem.K;
b = problem.b;
eta = problem.eta;
n = size(K, 2);

x = problem.x0;
y = compatibleVector(problem, 'y0', x, n);
lambda = compatibleVector(problem, 'lambda0', zeros(n, 1), n);
xGolden = x;
normalRightBase = K' * b;

if keepOnlyVisualizationData
    history = struct('obj', [], 'feas', []);
else
    history = struct('obj', [], 'feas', [], 'time', [], 'err', [], 'inner', []);
end

for k = 1:maxIterations
    xGolden = ((psi - 1) / psi) * x + (1 / psi) * xGolden;

    xProximalCoefficient = xStrongModulus / tau;
    xCoefficient = ridgeCurvature + xProximalCoefficient;
    xRightHandSide = normalRightBase - lambda + xProximalCoefficient * xGolden;
    woodburyMatrix = eye(size(K, 1)) + (K * K') / xCoefficient;
    [woodburyFactor, woodburyFlag] = chol(woodburyMatrix);
    projectedRightHandSide = K * xRightHandSide;
    if woodburyFlag == 0
        correction = woodburyFactor \ (woodburyFactor' \ projectedRightHandSide);
    else
        correction = woodburyMatrix \ projectedRightHandSide;
    end
    x = xRightHandSide / xCoefficient - K' * correction / (xCoefficient^2);

    yProximalCoefficient = yProximalWeight / t;
    yCoefficient = sigma + yProximalCoefficient;
    yCenter = (sigma * x + lambda + yProximalCoefficient * y) / yCoefficient;
    y = Functions.softThreshold(yCenter, eta / yCoefficient);

    lambda = lambda + sigma * (x - y);

    history = Tools.appendHistory(history, problem, x, y, k, startTime, options);
    if ~keepOnlyVisualizationData
        history.inner(k, 1) = 0;
    end
    if shouldStop(history, k, tolerance)
        break;
    end

    t = 0.5 * (1 + sqrt(1 + 4 * t^2));
    sigma = rho * t;
    tau = psi * xStrongModulus / (sigma * problem.constraintNorm2);
end

out = Tools.packOut('AGrpADMM.elastic', x, y, lambda, history, options);
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
