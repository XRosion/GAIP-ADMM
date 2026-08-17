function out = elastic(problem, options, parameters)
%AGAIPADMM.ELASTIC AGAIP-ADMM_spp with the nonsmooth block first.
%
% Public representation:
%
%   min_{x,y} gamma*||x||^2 + 0.5*||K*x-b||^2 + eta*||y||_1,
%   s.t.      x-y=0.
%
% Paper's Algorithm 1 order (restored to the public order on output):
%
%   first block  x := public y,  f(x) = eta*||x||_1,
%   second block y := public x,  g(y) = gamma*||y||^2
%                                        + 0.5*||K*y-b||^2,
%   constraint  psi(x,y) = -x+y = 0.
%
% The x-step is exact and the strongly convex y-step is solved by CG.

maxIter = options.maxIterations;

K = problem.K;
b = problem.b;
eta = problem.eta;

two_gamma = problem.ridgeCurvature;
gamma_f = problem.preStrongModulusF;
gamma_g = problem.preStrongModulusG;
norm_B2 = problem.constraintNorm2;

maxCgIter = options.elasticCgMaxIter;
cgAbsTol = options.elasticCgAbsTol;
compactHistory = Tools.keepOnlyVisualizationData(options);

alpha = parameters.alpha;
betakm1 = parameters.beta;
theta_x = parameters.theta_x;
theta_y = parameters.theta_y;
theta_lambda = parameters.theta_lambda;
theta_gamma = parameters.theta_gamma;
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

gamma_hat = (1 - theta_gamma) * gamma_f + theta_gamma * gamma_g;
delta_beta = gamma_hat / (3 * (1 + norm_B2));

Ktb = K' * b;

if compactHistory
    history = struct('obj', [], 'feas', []);
else
    history = struct('obj', [], 'feas', [], 'time', [], 'err', [], 'inner', [], ...
        'beta', [], ...
        'relErr2Left', [], 'relErr2Right', [], 'relErr2Ok', [], ...
        'innerResNorm', []);
end
startTime = tic;

% MATLAB iter=1 corresponds to k=0 in Algorithm 1.
for iter = 1:maxIter
    betak = betakm1 + delta_beta;

    bar_xk = (1 - theta_x) * xk + theta_x * bar_xkm1;
    bar_yk = (1 - theta_y) * yk + theta_y * bar_ykm1;
    bar_lambdak = (1 - theta_lambda) * lambdak + ...
        theta_lambda * bar_lambdakm1;

    % Step 3 (paper x): exact l1 proximal update.
    coef_x = betak + tau1 / betak;
    center_x = (betak * yk - bar_lambdak + ...
        (tau1 / betak) * bar_xk) / coef_x;
    xkk = Functions.softThreshold(center_x, eta / coef_x);

    % Step 4 (paper y): inexact strongly convex quadratic update.
    coef_y = two_gamma + betak + tau2 / betak;
    rhs_y = Ktb + betak * xkk + bar_lambdak + ...
        (tau2 / betak) * bar_yk;

    ykk = yk;
    d2kk = K' * (K * ykk) + coef_y * ykk - rhs_y;

    lhs2 = 2 * abs((w2k - ykk)' * d2kk) + ...
        alpha * betak * (d2kk' * d2kk);
    rhs2 = epsilon2 * (1 + theta_y) * (tau2 / betak) * ...
        norm(ykk - bar_yk)^2 + epsilon3 * betak * norm(ykk - xkk)^2;
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
            alpha * betak * (d2kk' * d2kk);
        rhs2 = epsilon2 * (1 + theta_y) * (tau2 / betak) * ...
            norm(ykk - bar_yk)^2 + epsilon3 * betak * norm(ykk - xkk)^2;
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
    lambdakk = bar_lambdak - alpha * betak * psikk;
    w2kk = w2k - alpha * betak * d2kk;

    % Public order is signal y first, sparse auxiliary x second.
    history = Tools.appendHistory( ...
        history, problem, ykk, xkk, iter, startTime, options);
    if ~compactHistory
        history.inner(iter, 1) = cgIter;
        history.beta(iter, 1) = betak;
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
    betakm1 = betak;
end

out = Tools.packOut('AGAIPADMM.elastic', yk, xk, lambdak, history, options);
out.parameters = parameters;
out.penaltyUpdate = struct( ...
    'gammaF', gamma_f, ...
    'gammaG', gamma_g, ...
    'gammaHat', gamma_hat, ...
    'constraintNorm2', norm_B2, ...
    'betaIncrement', delta_beta);
end
