function out = tv(problem, options, parameters)
%AGAIPADMM.TV AGAIP-ADMM_spp with the TV block first.
%
% Public representation:
%
%   min_{x,y} 0.5*mu*||x-c||^2 + sum_i ||y_i||_2,
%   s.t.      D*x-y=0.
%
% Paper's Algorithm 1 order (restored to the public order on output):
%
%   first block  x := public y,  f(x) = sum_i ||x_i||_2,
%   second block y := public x,  g(y) = 0.5*mu*||y-c||^2,
%   constraint  psi(x,y) = -x+D*y = 0.
%
% The x-step has a closed form. The strongly convex y-step uses CG.

mu = problem.mu;
c = problem.c;
D = problem.D;
DT = problem.DT;

gamma_f = problem.preStrongModulusF;
gamma_g = problem.preStrongModulusG;
norm_B2 = problem.constraintNorm2;

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

maxIter = options.maxIterations;
outerTolerance = options.tolerance;
maxCgIter = options.tvCgMaxIter;
cgAbsTol = options.tvCgAbsTol;
compactHistory = Tools.keepOnlyVisualizationData(options);

% Paper variables x/y are the public TV auxiliary/image variables, respectively.
xk = problem.y0;
yk = problem.x0;
lambdak = problem.lambda0;

bar_xkm1 = xk;
bar_ykm1 = yk;
bar_lambdakm1 = lambdak;
w2k = yk;

gamma_hat = (1 - theta_gamma) * gamma_f + theta_gamma * gamma_g;
delta_beta = gamma_hat / (3 * (1 + norm_B2));

mu_c = mu * c;

if compactHistory
    history = struct('obj', [], 'feas', []);
else
    history = struct('obj', [], 'feas', [], 'time', [], 'err', [], 'inner', [], ...
        'innerTotal', [], 'beta', [], 'alpha', [], ...
        'theta_x', [], 'theta_y', [], 'theta_lambda', [], 'theta_gamma', [], ...
        'tau1', [], 'tau2', [], ...
        'relErr2Left', [], 'relErr2Right', [], 'relErr2Ok', [], ...
        'absStop2', [], 'linearResidual2', [], 'innerResNorm', [], ...
        'primalResidualRel', [], 'dualResidualRel', [], 'outerStopMetric', []);
end

startTime = tic;

% MATLAB iter=1 corresponds to k=0 in Algorithm 1.
for iter = 1:maxIter
    betak = betakm1 + delta_beta;

    bar_xk = (1 - theta_x) * xk + theta_x * bar_xkm1;
    bar_yk = (1 - theta_y) * yk + theta_y * bar_ykm1;
    bar_lambdak = (1 - theta_lambda) * lambdak + ...
        theta_lambda * bar_lambdakm1;

    % Step 3 (paper x): exact isotropic-TV proximal update.
    Dyk = D(yk);
    coef_x = betak + tau1 / betak;
    center_x = (betak * Dyk - bar_lambdak + ...
        (tau1 / betak) * bar_xk) / coef_x;
    xkk = Functions.groupSoftThreshold(center_x, 1 / coef_x);

    % Step 4 (paper y): inexact strongly convex image-space update.
    coef_y = mu + tau2 / betak;
    rhs_y = mu_c + DT(bar_lambdak + betak * xkk) + ...
        (tau2 / betak) * bar_yk;
    H_y = @(z) coef_y * z + betak * DT(D(z));

    ykk = yk;
    d2kk = H_y(ykk) - rhs_y;

    lhs2 = 2 * abs((w2k - ykk)' * d2kk) + ...
        alpha * betak * (d2kk' * d2kk);
    rhs2 = epsilon2 * (1 + theta_y) * (tau2 / betak) * ...
        norm(ykk - bar_yk)^2 + epsilon3 * betak * norm(D(ykk) - xkk)^2;
    relerr2_ok = lhs2 <= rhs2;

    r_cg = -d2kk;
    p_cg = r_cg;
    r_cg_norm2 = r_cg' * r_cg;
    cgIter = 0;

    while cgIter < maxCgIter && ~(relerr2_ok || sqrt(r_cg_norm2) <= cgAbsTol)
        Ap_cg = H_y(p_cg);
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
            norm(ykk - bar_yk)^2 + epsilon3 * betak * norm(D(ykk) - xkk)^2;
        relerr2_ok = lhs2 <= rhs2;

        if relerr2_ok || sqrt(r_cg_norm2_new) <= cgAbsTol
            r_cg_norm2 = r_cg_norm2_new;
            break;
        end

        p_cg = r_cg + (r_cg_norm2_new / r_cg_norm2) * p_cg;
        r_cg_norm2 = r_cg_norm2_new;
    end

    absstop2 = sqrt(r_cg_norm2) <= cgAbsTol;
    linear_residual2 = norm(d2kk) / max(norm(rhs_y), 1);

    Dykk = D(ykk);
    psikk = Dykk - xkk;
    lambdakk = bar_lambdak - alpha * betak * psikk;
    w2kk = w2k - alpha * betak * d2kk;

    % Public order is image y first, TV auxiliary x second.
    history = Tools.appendHistory( ...
        history, problem, ykk, xkk, iter, startTime, options);
    if ~compactHistory
        history.inner(iter, 1) = cgIter;
        history.innerTotal(iter, 1) = sum(history.inner(1:iter));
        history.beta(iter, 1) = betak;
        history.alpha(iter, 1) = alpha;
        history.theta_x(iter, 1) = theta_x;
        history.theta_y(iter, 1) = theta_y;
        history.theta_lambda(iter, 1) = theta_lambda;
        history.theta_gamma(iter, 1) = theta_gamma;
        history.tau1(iter, 1) = tau1;
        history.tau2(iter, 1) = tau2;
        history.relErr2Left(iter, 1) = lhs2;
        history.relErr2Right(iter, 1) = rhs2;
        history.relErr2Ok(iter, 1) = relerr2_ok;
        history.absStop2(iter, 1) = absstop2;
        history.linearResidual2(iter, 1) = linear_residual2;
        history.innerResNorm(iter, 1) = norm(d2kk);

        primal_res_rel = norm(psikk) / max([norm(Dykk), norm(xkk), 1]);
        dual_res = betak * norm(D(ykk - yk));
        dual_scale = max(1, norm(lambdak));
        dual_res_rel = dual_res / dual_scale;
        history.primalResidualRel(iter, 1) = primal_res_rel;
        history.dualResidualRel(iter, 1) = dual_res_rel;
        history.outerStopMetric(iter, 1) = max(primal_res_rel, dual_res_rel);
    end

    xk = xkk;
    yk = ykk;
    lambdak = lambdakk;
    bar_xkm1 = bar_xk;
    bar_ykm1 = bar_yk;
    bar_lambdakm1 = bar_lambdak;
    w2k = w2kk;
    betakm1 = betak;

    if ~compactHistory && history.outerStopMetric(iter) <= outerTolerance
        fields = fieldnames(history);
        for fieldIndex = 1:numel(fields)
            value = history.(fields{fieldIndex});
            if (isnumeric(value) || islogical(value)) && size(value, 1) >= iter
                history.(fields{fieldIndex}) = value(1:iter, :);
            end
        end
        break;
    end
end

out = Tools.packOut('AGAIPADMM.tv', yk, xk, lambdak, history, options);
out.parameters = parameters;
out.penaltyUpdate = struct( ...
    'gammaF', gamma_f, ...
    'gammaG', gamma_g, ...
    'gammaHat', gamma_hat, ...
    'constraintNorm2', norm_B2, ...
    'betaIncrement', delta_beta);
end
