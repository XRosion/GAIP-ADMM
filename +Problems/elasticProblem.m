function problem = elasticProblem(m, n, s)
%ELASTICPROBLEM Generate an elastic-net sparse recovery instance.
%   PROBLEM = ELASTICPROBLEM(M, N, S) creates data for
%
%       min_x eta*||x||_1 + gamma*||x||_2^2 + 0.5*||K*x-b||_2^2.
%
%   The data generation follows the LASSO experiment in Zhou-Xu-Yang
%   (PADMC, Sec. 5.1): a correlated Gaussian sensing matrix, an s-sparse
%   ground-truth vector, and additive N(0,0.01) Gaussian noise. The
%   elastic-net experiment uses the split form
%
%       min_{x,y} eta*||x||_1 + gamma*||y||_2^2
%                 + 0.5*||K*y-b||_2^2,  s.t. x - y = 0.
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
eta = 0.01;
gamma = 0.005;

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

problem.type = 'elastic';
problem.K = K;
problem.b = b;
problem.xtrue = xtrue;
problem.eta = eta;
problem.gamma = gamma;
problem.x0 = zeros(n, 1);
problem.y0 = zeros(n, 1);
problem.lambda0 = zeros(n, 1);
problem.normK = norm(K, 2);
problem.normK2 = problem.normK^2;
problem.constraintNorm = 1;
problem.constraintNorm2 = 1;
problem.normH = problem.constraintNorm;
problem.normH2 = problem.constraintNorm2;
% Algorithm 1 order: l1 block first, smooth elastic-net block second.
problem.ridgeCurvature = 2 * gamma;
problem.leastSquaresPreStrongModulus = ...
    Functions.leastSquaresPreStrongModulus(K);
problem.preStrongModulusF = 0;
problem.preStrongModulusG = problem.ridgeCurvature + ...
    problem.leastSquaresPreStrongModulus;
problem.m = m;
problem.n = n;
problem.s = s;
problem.title = 'Elastic-net regularization';

[problem.PhiStar, problem.xStar, problem.yStar, problem.cvxInfo] = ...
    solveElasticPhiStarCvx(K, b, eta, gamma);
end

function [phiStar, xStar, yStar, cvxInfo] = solveElasticPhiStarCvx(K, b, eta, gamma)
%SOLVEELASTICPHISTARCVX Compute the reference elastic-net optimum by CVX.
%   For the split min eta*||x||_1 + gamma*||y||^2 + 0.5*||K*y-b||^2
%   s.t. x = y, the optimal x and y coincide.

n = size(K, 2);

cvx_begin quiet
    variable x(n)
    minimize( eta * norm(x, 1) + gamma * sum_square(x) + 0.5 * sum_square(K * x - b) )
cvx_end

if strcmp(cvx_status, 'Solved') || strcmp(cvx_status, 'Optimal')
    xStar = x;
    yStar = x;
    phiStar = cvx_optval;
    cvxInfo = struct('solver', cvx_solver, ...
                     'status', cvx_status, ...
                     'precision', 'default');
else
    error('elasticProblem:CvxFailed', ...
        'CVX failed to solve the elastic-net reference problem. Status: %s', cvx_status);
end
end
