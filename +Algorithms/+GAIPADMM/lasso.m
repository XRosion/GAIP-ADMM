function out = lasso(problem, options, parameters)
%GAIPADMM.LASSO Solve LASSO with the nonsmooth block updated first.
%
% Public representation:
%
%   min_{x,y} 0.5*||K*x-b||^2 + eta*||y||_1,  s.t. x-y=0.
%
% Paper's Algorithm 1 order (restored to the public order on output):
%
%   first block  x := public y,  f(x) = eta*||x||_1,
%   second block y := public x,  g(y) = 0.5*||K*y-b||^2,
%   constraint  psi(x,y) = -x+y = 0.
%
% The x-step has a closed form. The y-step is solved inexactly by CG.

maxIter = options.maxIterations;

K = problem.K;
b = problem.b;
eta = problem.eta;

maxCgIter = options.lassoCgMaxIter;
cgAbsTol = options.lassoCgAbsTol;
compactHistory = Tools.keepOnlyVisualizationData(options);

alpha = parameters.alpha;
beta = parameters.beta;
theta_x = parameters.theta_x;
theta_y = parameters.theta_y;
theta_lambda = parameters.theta_lambda;
epsilon2 = parameters.eps2;
epsilon3 = parameters.eps3;
tau1 = parameters.tau1;
tau2 = parameters.tau2;

% Paper variables x/y are the public auxiliary/signal variables, respectively.
xk = problem.y0;
yk = problem.x0;
lambdak = problem.lambda0;

bar_xkm1 = xk;
bar_ykm1 = yk;
bar_lambdakm1 = lambdak;
w2k = yk;

Ktb = K' * b;

if compactHistory
    history = struct('obj', [], 'feas', []);
else
    history = struct('obj', [], 'feas', [], 'time', [], 'err', [], 'inner', [], ...
        'relErr2Left', [], 'relErr2Right', [], 'relErr2Ok', [], ...
        'innerResNorm', []);
end
startTime = tic;

% MATLAB iter=1 corresponds to k=0 in Algorithm 1.
for iter = 1:maxIter
    bar_xk = (1 - theta_x) * xk + theta_x * bar_xkm1;
    bar_yk = (1 - theta_y) * yk + theta_y * bar_ykm1;
    bar_lambdak = (1 - theta_lambda) * lambdak + ...
        theta_lambda * bar_lambdakm1;

    % Step 3 (paper x): exact l1 proximal update.
    coef_x = beta + tau1 / beta;
    center_x = (beta * yk - bar_lambdak + ...
        (tau1 / beta) * bar_xk) / coef_x;
    xkk = Functions.softThreshold(center_x, eta / coef_x);

    % Step 4 (paper y): inexact smooth quadratic update.
    coef_y = beta + tau2 / beta;
    rhs_y = Ktb + beta * xkk + bar_lambdak + ...
        (tau2 / beta) * bar_yk;

    ykk = yk;
    d2kk = K' * (K * ykk) + coef_y * ykk - rhs_y;

    lhs2 = 2 * abs((w2k - ykk)' * d2kk) + ...
        alpha * beta * (d2kk' * d2kk);
    rhs2 = epsilon2 * (1 + theta_y) * (tau2 / beta) * ...
        norm(ykk - bar_yk)^2 + epsilon3 * beta * norm(ykk - xkk)^2;
    relerr2_ok = lhs2 <= rhs2;

    r_cg = -d2kk;
    p_cg = r_cg;
    r_cg_norm2 = r_cg' * r_cg;
    cgIter = 0;

    while cgIter < maxCgIter && ~(relerr2_ok || sqrt(r_cg_norm2) <= cgAbsTol)
        Ap_cg = K' * (K * p_cg) + coef_y * p_cg;
        pAp = p_cg' * Ap_cg;
        if pAp <= eps
            break;
        end

        alpha_cg = r_cg_norm2 / pAp;
        ykk = ykk + alpha_cg * p_cg;
        r_cg = r_cg - alpha_cg * Ap_cg;
        d2kk = -r_cg;
        r_cg_norm2_new = r_cg' * r_cg;
        cgIter = cgIter + 1;

        lhs2 = 2 * abs((w2k - ykk)' * d2kk) + ...
            alpha * beta * (d2kk' * d2kk);
        rhs2 = epsilon2 * (1 + theta_y) * (tau2 / beta) * ...
            norm(ykk - bar_yk)^2 + epsilon3 * beta * norm(ykk - xkk)^2;
        relerr2_ok = lhs2 <= rhs2;

        if relerr2_ok || sqrt(r_cg_norm2_new) <= cgAbsTol
            r_cg_norm2 = r_cg_norm2_new;
            break;
        end

        p_cg = r_cg + (r_cg_norm2_new / r_cg_norm2) * p_cg;
        r_cg_norm2 = r_cg_norm2_new;
    end

    relerr2_ok = relerr2_ok || sqrt(r_cg_norm2) <= cgAbsTol;

    psikk = ykk - xkk;
    lambdakk = bar_lambdak - alpha * beta * psikk;
    w2kk = w2k - alpha * beta * d2kk;

    % Public order is signal y first, sparse auxiliary x second.
    history = Tools.appendHistory( ...
        history, problem, ykk, xkk, iter, startTime, options);
    if ~compactHistory
        history.inner(iter, 1) = cgIter;
        history.relErr2Left(iter, 1) = lhs2;
        history.relErr2Right(iter, 1) = rhs2;
        history.relErr2Ok(iter, 1) = relerr2_ok;
        history.innerResNorm(iter, 1) = norm(d2kk);
    end

    xk = xkk;
    yk = ykk;
    lambdak = lambdakk;
    bar_xkm1 = bar_xk;
    bar_ykm1 = bar_yk;
    bar_lambdakm1 = bar_lambdak;
    w2k = w2kk;
end

out = Tools.packOut('GAIPADMM.lasso', yk, xk, lambdak, history, options);
end
