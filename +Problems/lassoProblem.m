function problem = lassoProblem(m, n, s)
%LASSOPROBLEM Generate a LASSO sparse recovery instance.
%   PROBLEM = LASSOPROBLEM(M, N, S) creates data for
%
%       min_x eta*||x||_1 + 0.5*||K*x-b||_2^2.
%
%   The data generation follows the LASSO experiment in Zhou-Xu-Yang
%   (PADMC, Sec. 5.1): a correlated Gaussian sensing matrix, an s-sparse
%   ground-truth vector, and additive N(0,0.01) Gaussian noise.
%
%   The ADMM implementation uses the split form
%
%       min_{x,y} eta*||y||_1 + 0.5*||K*x - b||_2^2,  s.t. x - y = 0.
%
%   Input:
%       m - number of observations (number of rows in K).
%       n - number of variables (dimension of x and y).
%       s - sparsity level.
%
%   Output:
%       problem - structure with data, split initial values, and reference
%                 objective value computed by CVX when available.

rng(33);
nu = 0.5;
noiseVariance = 0.01;
regularizationParameter = 0.01;

gaussianInnovations = randn(m, n);
K = zeros(m, n);
K(:, 1) = gaussianInnovations(:, 1) / sqrt(1 - nu^2);
for j = 2:n
    K(:, j) = nu * K(:, j - 1) + gaussianInnovations(:, j);
end

xtrue = zeros(n, 1);
indexSet = randperm(n, s);
xtrue(indexSet) = randn(s, 1);
b = K * xtrue + sqrt(noiseVariance) * randn(m, 1);

problem.type = 'lasso';
problem.K = K;
problem.b = b;
problem.xtrue = xtrue;
problem.eta = regularizationParameter;
problem.x0 = zeros(n, 1);
problem.y0 = zeros(n, 1);      % y is also n-dimensional (split: x - y = 0)
problem.lambda0 = zeros(n, 1); % multiplier has same dimension as x and y
problem.normK = norm(K, 2);
problem.normK2 = problem.normK^2;
problem.constraintNorm = 1;
problem.constraintNorm2 = 1;
% Algorithm 1 order: l1 block first, least-squares block second.
problem.leastSquaresPreStrongModulus = ...
    Functions.leastSquaresPreStrongModulus(K);
problem.preStrongModulusF = 0;
problem.preStrongModulusG = problem.leastSquaresPreStrongModulus;
problem.m = m;
problem.n = n;
problem.s = s;

[problem.PhiStar, problem.xStar, problem.yStar, problem.cvxInfo] = ...
    solveLassoPhiStarCvx(K, b, problem.eta);
end

function [phiStar, xStar, yStar, cvxInfo] = solveLassoPhiStarCvx(K, b, regularizationParameter)
%SOLVELASSOPHISTARCVX Compute the reference LASSO optimum by CVX.
%   For the split min eta*||y||_1 + 0.5*||Kx-b||^2  s.t. x = y,
%   the optimal y equals the optimal x.

n = size(K, 2);

cvx_begin quiet
    variable x(n)
    minimize( regularizationParameter * norm(x, 1) + 0.5 * sum_square(K * x - b) )
cvx_end

if strcmp(cvx_status, 'Solved') || strcmp(cvx_status, 'Optimal')
    xStar = x;
    yStar = x;   % Under the constraint x = y, the optimal y is x
    phiStar = cvx_optval;
    cvxInfo = struct('solver', cvx_solver, ...
                     'status', cvx_status, ...
                     'precision', 'default');
else
    error('lassoProblem:CvxFailed', ...
        'CVX failed to solve the LASSO reference problem. Status: %s', cvx_status);
end
end
