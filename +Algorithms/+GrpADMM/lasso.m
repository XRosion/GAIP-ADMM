function out = lasso(problem, options)
% GrpADMM.lasso  Golden-ratio style ADMM for the x-y split LASSO problem.
%
%   min  0.5*||K*x-b||^2 + eta*||y||_1,  s.t. x-y=0.
%
% The smooth x-block is solved by a cached Woodbury factorization and the
% nonsmooth y-block is solved by soft-thresholding.

maxIterations = options.maxIterations;
psi = getOption(options, 'grpPsi', 1.618);
sigma = getOption(options, 'grpSigma', 1);
keepOnlyVisualizationData = Tools.keepOnlyVisualizationData(options);
startTime = tic;

K = problem.K;
b = problem.b;
eta = problem.eta;
n = size(K, 2);

x = problem.x0;
y = compatibleVector(problem, 'y0', x, n);
lambda = compatibleVector(problem, 'lambda0', zeros(n, 1), n);
yGolden = y;

woodburyMatrix = eye(size(K, 1)) + (K * K') / sigma;
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
    x = xRightHandSide / sigma - K' * correction / (sigma^2);

    y = Functions.softThreshold(x - lambda / sigma, eta / sigma);
    lambda = lambda - sigma * (x - y);

    history = Tools.appendHistory(history, problem, x, y, k, startTime, options);
    if ~keepOnlyVisualizationData
        history.inner(k, 1) = 0;
    end
end

out = Tools.packOut('GrpADMM.lasso', x, y, lambda, history, options);
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
