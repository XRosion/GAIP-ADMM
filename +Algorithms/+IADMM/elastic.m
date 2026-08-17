function out = elastic(problem, options)
% IADMM.elastic  用 x-y=0 分裂形式求解弹性网络正则化问题。
%
%   min_z  eta*||z||_1 + gamma*||z||_2^2 + 0.5*||K*z-b||_2^2
%
% 采用分裂形式：
%
%   min_{x,y}  gamma*||x||_2^2 + 0.5*||K*x-b||_2^2 + eta*||y||_1,  s.t. x-y=0.
%
% x 子问题为光滑强凸二次问题，由 CG 按 IADMM 相对误差准则非精确求解；
% y 子问题为 l1 近端问题，由软阈值算子精确求解。
%
% 参数说明（通过 options 结构体传入）：
%   options.maxIterations       - 外迭代最大次数
%   options.elasticCgMaxIter    - CG 内迭代最大次数 (默认 100)
%   options.elasticCgAbsTol     - CG 绝对残差容差 (默认 1e-6)
%   keepOnlyVisualizationData   - 是否只保留绘图所需数据 (由 Tools 处理)

    maxIterations = options.maxIterations;
    beta = 1;
    alpha = 0.4;
    sigma = 0.99;
    gammaValue = 1 - alpha;
    tauUpperBound = (2 * gammaValue - 1 + sqrt(4 * gammaValue^2 + 1)) / (2 * gammaValue);
    tau = tauUpperBound - (tauUpperBound - 1) / 1000;
    maxInnerIterations = options.elasticCgMaxIter;
    xAbsTol = options.elasticCgAbsTol;
    
    keepOnlyVisualizationData = Tools.keepOnlyVisualizationData(options);
    startTime = tic;

    K = problem.K;
    b = problem.b;
    eta = problem.eta;
    ridgeCurvature = problem.ridgeCurvature;
    n = size(K, 2);

    x = problem.x0;
    y = problem.y0;
    lambda = problem.lambda0;
    auxiliaryX = x;

    normalRightBase = K' * b;

    % 初始化历史记录
    if keepOnlyVisualizationData
        history = struct('obj', [], 'feas', []);
    else
        history = struct('obj', [], 'feas', [], 'time', [], 'err', [], 'inner', [], ...
                         'relErrXLeft', [], 'relErrXRight', [], 'relErrXOk', [], ...
                         'relErrYLeft', [], 'relErrYRight', [], 'relErrYOk', [], ...
                         'innerResNorm', []);
    end

    for iteration = 1:maxIterations
        previousY = y;

        % ----- 1. 非精确求解 x 子问题 (CG 满足相对误差准则) -----
        xRightHandSide = normalRightBase + beta * y + lambda;
        xCandidate = x;
        grad_f = ridgeCurvature*xCandidate + K'*(K*xCandidate - b);
        d = grad_f - (lambda - beta*(xCandidate - previousY));
        relLeft = (2 / beta) * abs((auxiliaryX - xCandidate)' * d) + (d' * d);
        relRight = sigma * alpha * norm(xCandidate - previousY)^2;
        relOk = relLeft <= relRight;

        residual = -d;
        direction = residual;
        residualSquared = residual' * residual;
        innerIteration = 0;

        while innerIteration < maxInnerIterations && ...
                ~(relOk || sqrt(residualSquared) <= xAbsTol)
            Mv = K'*(K*direction) + ...
                (beta + ridgeCurvature) * direction;
            denominator = direction' * Mv;
            if denominator <= eps
                break;
            end

            stepSize = residualSquared / denominator;
            xCandidate = xCandidate + stepSize * direction;
            residual = residual - stepSize * Mv;
            d = -residual;
            newResidualSquared = residual' * residual;
            innerIteration = innerIteration + 1;

            relLeft = (2 / beta) * abs((auxiliaryX - xCandidate)' * d) + (d' * d);
            relRight = sigma * alpha * norm(xCandidate - previousY)^2;
            relOk = relLeft <= relRight;
            if relOk || sqrt(newResidualSquared) <= xAbsTol
                residualSquared = newResidualSquared;
                break;
            end

            direction = residual + (newResidualSquared / residualSquared) * direction;
            residualSquared = newResidualSquared;
        end
        relOk = relOk || sqrt(residualSquared) <= xAbsTol;
        x = xCandidate;

        % ----- 2. y 子问题精确求解 -----
        y = Functions.softThreshold(x - lambda / beta, eta / beta);

        % ----- 3. 乘子更新与辅助变量更新 -----
        constraintResidual = x - y;
        lambda = lambda - tau * beta * constraintResidual;
        auxiliaryX = auxiliaryX - beta * d;

        % ----- 记录历史 -----
        history = Tools.appendHistory(history, problem, x, y, iteration, startTime, options);
        if ~keepOnlyVisualizationData
            history.inner(iteration, 1) = innerIteration;
            history.relErrXLeft(iteration, 1) = relLeft;
            history.relErrXRight(iteration, 1) = relRight;
            history.relErrXOk(iteration, 1) = relOk;
            history.relErrYLeft(iteration, 1) = 0;
            history.relErrYRight(iteration, 1) = 0.99 * norm(y - previousY)^2;
            history.relErrYOk(iteration, 1) = true;
            history.innerResNorm(iteration, 1) = norm([d; zeros(size(y))]);
        end
    end

    out = Tools.packOut('IADMM.elastic', x, y, lambda, history, options);
end
