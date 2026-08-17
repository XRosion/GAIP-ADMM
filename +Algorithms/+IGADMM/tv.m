function out = tv(problem, options)
% IGADMM.TV  非精确广义 ADMM 求解各向同性 TV 去噪（精简版）
%
%   本函数实现 Wu 等 (2023) "Inexact generalized ADMM with relative error criteria"
%   中的 IGADMMa 算法 (Algorithm 1 + 准则 (11)(12))，用于求解 TV 去噪：
%        min_x  0.5*mu*||x-c||^2 + sum_i ||(Dx)_i||_2
%   采用分裂 Dx - y = 0，广义 ADMM 迭代包含：
%       1) 非精确求解 x（CG + 相对误差准则 (11)）
%       2) 松弛更新 lambdaHalf = lambda - alpha*beta*(Dx - y)
%       3) 精确求解 y（带近端项 T = yProx*I，闭式分组软阈值）
%       4) 更新 lambda 和辅助变量 v
%
%   固定算法参数（源自论文第 4.2 节图像复原实验）：
%       beta   = 1      (增广拉格朗日惩罚)
%       alpha  = -0.3   (松弛参数，∈(-1,1))
%       sigma1 = 0.99   (x-子问题相对误差因子)
%       yProx  = 0.01   (y-子问题近端正则系数，T = yProx*I)
%
%   输入：
%       problem : 结构体，必须包含 mu, c, D, DT, x0, y0, lambda0
%       options : 结构体，必须包含 maxIterations, tolerance, tvCgMaxIter, tvCgAbsTol
%                 可选 keepOnlyVisualizationData (默认 false)
%   输出：
%       out : 结构体，含 x, y, lambda, history, parameters, options

    % ---- 读取停止准则参数 ----
    maxIterations = options.maxIterations;
    outerTolerance = options.tolerance;
    tvCgMaxIter = options.tvCgMaxIter;
    tvCgAbsTol = options.tvCgAbsTol;
    if isfield(options, 'keepOnlyVisualizationData')
        keepOnlyVisualizationData = options.keepOnlyVisualizationData;
    else
        keepOnlyVisualizationData = false;
    end

    % ---- 固定算法参数 ----
    beta = 1;
    alpha = -0.3;
    sigma1 = 0.99;
    yProx = 0.01;

    % ---- 解包问题数据 ----
    mu = problem.mu;
    c = problem.c;
    D = problem.D;
    DT = problem.DT;

    % ---- 初始化变量 ----
    x = problem.x0;
    y = problem.y0;
    lambda = problem.lambda0;
    auxiliaryX = x;          % v1^k
    constantMuC = mu * c;

    % ---- 历史记录初始化 ----
    history = initializeHistory(keepOnlyVisualizationData);
    startTime = tic;

    % ---- 主循环 ----
    for iteration = 1:maxIterations
        previousY = y;
        previousLambda = lambda;

        % ------ 1. x 子问题（非精确 CG，满足相对误差准则 (11)） ------
        rightHandSide = constantMuC + DT(lambda + beta * y);
        linearOperator = @(v) mu * v + beta * DT(D(v));

        % CG 初始化
        x_cg = x;                     % warm-start
        d_cg = linearOperator(x_cg) - rightHandSide;   % Ax - b
        % 相对误差准则 (11) 右侧不包含 S 项，因为 S=0
        right_rhs = (1 - alpha) * sigma1 * beta * norm(D(x_cg) - y)^2;
        left_val = 2 * abs((auxiliaryX - x_cg)' * d_cg) + (d_cg' * d_cg);
        relOk = left_val <= right_rhs;
        residual_cg = -d_cg;
        direction_cg = residual_cg;
        resSq_cg = residual_cg' * residual_cg;
        innerIter = 0;

        % CG 迭代
        while innerIter < tvCgMaxIter && ~(relOk || sqrt(resSq_cg) <= tvCgAbsTol)
            Ad = linearOperator(direction_cg);
            denom = direction_cg' * Ad;
            if denom <= eps
                break;
            end
            step = resSq_cg / denom;
            x_cg = x_cg + step * direction_cg;
            residual_cg = residual_cg - step * Ad;
            d_cg = -residual_cg;
            newResSq = residual_cg' * residual_cg;
            innerIter = innerIter + 1;

            % 重新计算相对误差
            right_rhs = (1 - alpha) * sigma1 * beta * norm(D(x_cg) - y)^2;
            left_val = 2 * abs((auxiliaryX - x_cg)' * d_cg) + (d_cg' * d_cg);
            relOk = left_val <= right_rhs;

            if relOk || sqrt(newResSq) <= tvCgAbsTol
                resSq_cg = newResSq;
                break;
            end
            direction_cg = residual_cg + (newResSq / resSq_cg) * direction_cg;
            resSq_cg = newResSq;
        end

        xNew = x_cg;
        d = d_cg;
        relLeft = left_val;
        relRight = right_rhs;
        absStop = sqrt(resSq_cg) <= tvCgAbsTol;
        linearResidual = norm(d) / max(norm(rightHandSide), 1);

        DxNew = D(xNew);

        % ------ 2. 中间 lambda 更新（松弛步） ------
        lambdaHalf = lambda - alpha * beta * (DxNew - y);

        % ------ 3. y 子问题（带近端项 T = yProx*I，精确闭式解） ------
        % 近端广义 ADMM 中 y 更新：min_y  L_beta + 0.5*||y-y^k||_T^2
        % 对于各向同性 TV，其闭式解为分组软阈值
        yCoeff = beta + yProx;
        yCenter = (beta * DxNew - lambdaHalf + yProx * y) / yCoeff;
        yNew = Functions.groupSoftThreshold(yCenter, 1 / yCoeff);

        % ------ 4. lambda 最终更新 和 辅助变量更新 ------
        lambdaNew = lambdaHalf - beta * (DxNew - yNew);
        auxiliaryX = auxiliaryX - d;   % v1^{k+1} = v1^k - d1^{k+1}

        % ------ 记录历史 ------
        history = Tools.appendHistory(history, problem, xNew, yNew, iteration, startTime, options);
        if ~keepOnlyVisualizationData
            history.inner(iteration, 1) = innerIter;
            history.innerTotal(iteration, 1) = sum(history.inner(1:iteration));
            history.relErrXLeft(iteration, 1) = relLeft;
            history.relErrXRight(iteration, 1) = relRight;
            history.relErrXOk(iteration, 1) = relOk;
            history.xAbsStop(iteration, 1) = absStop;
            history.xLinearResidual(iteration, 1) = linearResidual;
            % y 子问题精确，故相对误差恒满足
            history.relErrYLeft(iteration, 1) = 0;
            history.relErrYRight(iteration, 1) = 0;
            history.relErrYOk(iteration, 1) = true;
            history.innerResNorm(iteration, 1) = norm([d; zeros(size(yNew))]);
            history.beta(iteration, 1) = beta;
            history.alpha(iteration, 1) = alpha;
            history.yProx(iteration, 1) = yProx;

            % 原始残差相对值
            primalRel = norm(DxNew - yNew) / max([norm(DxNew), norm(yNew), 1]);
            % 对偶残差相对值
            dualRes = beta * norm(DT(yNew - previousY));
            dualScale = max(1, norm(DT(previousLambda)));
            dualRel = dualRes / dualScale;
            history.primalResidualRel(iteration, 1) = primalRel;
            history.dualResidualRel(iteration, 1) = dualRel;
            history.outerStopMetric(iteration, 1) = max(primalRel, dualRel);
        end

        % ------ 外层停止判断 ------
        if ~keepOnlyVisualizationData && history.outerStopMetric(iteration) <= outerTolerance
            history = trimHistory(history, iteration);
            x = xNew; y = yNew; lambda = lambdaNew;
            break;
        end

        % 更新变量供下一轮
        x = xNew;
        y = yNew;
        lambda = lambdaNew;
    end

    % ---- 打包输出 ----
    out = Tools.packOut('IGADMM.tv', x, y, lambda, history, options);
    out.parameters = struct('beta', beta, 'alpha', alpha, 'sigma1', sigma1, 'yProx', yProx);
end

% ========= 局部辅助函数（仅历史记录相关） =========

function history = initializeHistory(keepOnlyVisualizationData)
    if keepOnlyVisualizationData
        history = struct('obj', [], 'feas', []);
    else
        history = struct('obj', [], 'feas', [], 'time', [], 'err', [], ...
            'inner', [], 'innerTotal', [], ...
            'relErrXLeft', [], 'relErrXRight', [], 'relErrXOk', [], ...
            'relErrYLeft', [], 'relErrYRight', [], 'relErrYOk', [], ...
            'xAbsStop', [], 'xLinearResidual', [], 'innerResNorm', [], ...
            'beta', [], 'alpha', [], 'yProx', [], ...
            'primalResidualRel', [], 'dualResidualRel', [], 'outerStopMetric', []);
    end
end

function history = trimHistory(history, iteration)
    fields = fieldnames(history);
    for i = 1:numel(fields)
        f = fields{i};
        val = history.(f);
        if (isnumeric(val) || islogical(val)) && size(val, 1) >= iteration
            history.(f) = val(1:iteration, :);
        end
    end
end