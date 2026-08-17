function out = tv(problem, options)
% IADMM.TV  基于非精确 ADMM 的各向同性 TV 图像去噪（精简内联版）
%
%   本函数实现了文献 "On inexact ADMMs with relative error criteria"
%   (Jiaxin Xie, 2018) 中提出的部分非精确 ADMM (Algorithm 1)，
%   用于求解各向同性全变分 (TV) 图像去噪问题：
%
%        min_x  0.5*mu*||x - c||_2^2  +  sum_i ||(Dx)_i||_2
%
%   其中 D 为离散梯度算子（各向同性 TV），c 为含噪观测图像。
%   采用分裂变量 Dx - y = 0，ADMM 迭代如下：
%       1) 非精确求解 x 子问题（使用预条件共轭梯度 CG），
%          满足相对误差准则 (论文式 11)；
%       2) 精确求解 y 子问题（分组软阈值）；
%       3) 以松弛步长 tau 更新对偶变量 lambda。
%
%   算法参数（beta, alpha, sigma, tau）固定为实验最优值，
%   以简化对比实验的调用接口。
%
%   参考文献：
%       Xie, J. (2018). On inexact ADMMs with relative error criteria.
%       Computational Optimization and Applications, 71(3), 743-765.
%
%   输入参数：
%       problem : 结构体，包含问题数据，必须有以下字段：
%           .mu        - 保真项权重（标量）
%           .c         - 观测图像，列向量 (N×1)
%           .D         - 梯度算子函数句柄  @(x) Dx
%           .DT        - 梯度伴随算子  @(y) D^T y
%           .x0        - 初始 x（列向量，N×1）
%           .y0        - 初始 y（列向量，2N×1，分组结构）
%           .lambda0   - 初始对偶变量（列向量，2N×1）
%
%       options : 结构体，用于控制停止准则，必须包含以下字段：
%           .maxIterations      - 最大外层迭代次数（正整数）
%           .tolerance          - 外层停止容差（标量，>0）
%           .tvCgMaxIter        - CG 最大迭代次数（正整数）
%           .tvCgAbsTol         - CG 绝对残差容差（标量，>0）
%           .keepOnlyVisualizationData (可选，逻辑值) -
%                 若为 true，只记录目标值和可行性，节省内存；
%                 默认 false。
%
%   输出参数：
%       out : 结构体，包含以下字段：
%           .name       - 算法名称 'IADMM.tv'
%           .x          - 最终恢复图像（列向量）
%           .y          - 最终辅助变量
%           .lambda     - 最终对偶变量
%           .history    - 迭代历史结构体（含目标值、残差、内层迭代数等）
%           .parameters - 所用算法参数结构体 (beta, alpha, sigma, tau)
%           .options    - 传入的 options 结构体
%
%   固定算法参数（依据论文实验选择）：
%       beta  = 30   : 增广拉格朗日惩罚参数
%       alpha = 0.4  : 松弛参数（控制误差容限与步长折中）
%       sigma = 0.99 : 相对误差准则中的收缩因子
%       tau   ~ 1.468: 对偶步长（由 alpha 自动计算，严格小于理论上限）
%
%   注意事项：
%       1) 所有向量均为列向量，尺寸为 N = m*n。
%       2) 函数依赖于外部工具箱：Functions.groupSoftThreshold、
%          Tools.appendHistory、Tools.packOut，请确保路径已添加。
%       3) 本实现针对各向同性 TV，梯度 D 返回 2N×1 的向量。
%       4) CG 采用 warm-start（以上次 x 为初始点），且使用相对误差
%          准则作为主要停止条件，绝对容差为辅助保护。

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
    beta = 30;
    alpha = 0.4;
    sigma = 0.99;
    % 计算 tau（基于 alpha）
    if alpha >= 1
        tau = 1;
    else
        gammaVal = 1 - alpha;
        upper = (2 * gammaVal - 1 + sqrt(4 * gammaVal^2 + 1)) / (2 * gammaVal);
        scaled = 1000 * upper;
        if abs(scaled - round(scaled)) <= 10 * eps(max(abs(scaled), 1))
            tau = (round(scaled) - 1) / 1000;
        else
            tau = floor(scaled) / 1000;
        end
    end

    % ---- 解包问题数据 ----
    mu = problem.mu;
    c = problem.c;
    D = problem.D;
    DT = problem.DT;

    % ---- 初始化变量 ----
    x = problem.x0;
    y = problem.y0;
    lambda = problem.lambda0;
    auxiliaryX = x;
    constantMuC = mu * c;

    % ---- 历史记录初始化 ----
    history = initializeHistory(keepOnlyVisualizationData);
    startTime = tic;

    % ---- 主循环 ----
    for iteration = 1:maxIterations
        previousY = y;
        previousLambda = lambda;

        % ------ x 子问题（非精确 CG，完全内联） ------
        rightHandSide = constantMuC + DT(lambda + beta * y);
        initialX = x;
        linearOperator = @(v) mu * v + beta * DT(D(v));

        % CG 初始化
        x_cg = initialX;
        d_cg = linearOperator(x_cg) - rightHandSide;   % Ax - b
        % 相对误差准则计算（内联）
        left_cg = (2 / beta) * abs((auxiliaryX - x_cg)' * d_cg) + (d_cg' * d_cg);
        right_cg = sigma * alpha * norm(D(x_cg) - previousY)^2;
        relOk = left_cg <= right_cg;
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
            alphaStep = resSq_cg / denom;
            x_cg = x_cg + alphaStep * direction_cg;
            residual_cg = residual_cg - alphaStep * Ad;
            d_cg = -residual_cg;
            newResSq = residual_cg' * residual_cg;
            innerIter = innerIter + 1;

            % 重新计算相对误差准则
            left_cg = (2 / beta) * abs((auxiliaryX - x_cg)' * d_cg) + (d_cg' * d_cg);
            right_cg = sigma * alpha * norm(D(x_cg) - previousY)^2;
            relOk = left_cg <= right_cg;

            if relOk || sqrt(newResSq) <= tvCgAbsTol
                resSq_cg = newResSq;
                break;
            end
            direction_cg = residual_cg + (newResSq / resSq_cg) * direction_cg;
            resSq_cg = newResSq;
        end

        % 最终 CG 结果
        x = x_cg;
        d = d_cg;
        relLeft = left_cg;
        relRight = right_cg;
        absStop = sqrt(resSq_cg) <= tvCgAbsTol;
        linearResidual = norm(d) / max(norm(rightHandSide), 1);

        Dx = D(x);

        % ------ y 子问题（精确分组软阈值） ------
        y = Functions.groupSoftThreshold(Dx - lambda / beta, 1 / beta);

        % ------ 约束残差与对偶更新 ------
        constraintResidual = Dx - y;
        lambda = lambda - tau * beta * constraintResidual;

        % ------ 辅助变量更新 ------
        auxiliaryX = auxiliaryX - beta * d;

        % ------ 记录历史 ------
        history = Tools.appendHistory(history, problem, x, y, iteration, startTime, options);
        if ~keepOnlyVisualizationData
            history.inner(iteration, 1) = innerIter;
            history.innerTotal(iteration, 1) = sum(history.inner(1:iteration));
            history.relErrXLeft(iteration, 1) = relLeft;
            history.relErrXRight(iteration, 1) = relRight;
            history.relErrXOk(iteration, 1) = relOk;
            history.xAbsStop(iteration, 1) = absStop;
            history.xLinearResidual(iteration, 1) = linearResidual;
            history.relErrYOk(iteration, 1) = true;
            history.relErrYLeft(iteration, 1) = 0;
            history.relErrYRight(iteration, 1) = 0;
            history.innerResNorm(iteration, 1) = norm([d; zeros(size(y))]);
            history.beta(iteration, 1) = beta;
            history.tau(iteration, 1) = tau;
            history.alpha(iteration, 1) = alpha;

            % 原始残差相对值（内联）
            primalRel = norm(Dx - y) / max([norm(Dx), norm(y), 1]);
            % 对偶残差相对值（内联）
            dualRes = beta * norm(DT(y - previousY));
            dualScale = max(1, norm(DT(previousLambda)));
            dualRel = dualRes / dualScale;
            history.primalResidualRel(iteration, 1) = primalRel;
            history.dualResidualRel(iteration, 1) = dualRel;
            history.outerStopMetric(iteration, 1) = max(primalRel, dualRel);
        end

        % ------ 外层停止判断 ------
        if ~keepOnlyVisualizationData && history.outerStopMetric(iteration) <= outerTolerance
            for f = fieldnames(history)'
                val = history.(f{1});
                if isnumeric(val) || islogical(val)
                    if size(val, 1) >= iteration
                        history.(f{1}) = val(1:iteration, :);
                    end
                end
            end
            break;
        end
    end

    % ---- 打包输出 ----
    out = Tools.packOut('IADMM.tv', x, y, lambda, history, options);
    out.parameters = struct('beta', beta, 'alpha', alpha, 'sigma', sigma, 'tau', tau);
end

% ========= 保留的辅助函数 =========

function history = initializeHistory(keepOnlyVisualizationData)
    if keepOnlyVisualizationData
        history = struct('obj', [], 'feas', []);
    else
        history = struct('obj', [], 'feas', [], 'time', [], 'err', [], ...
            'inner', [], 'innerTotal', [], ...
            'relErrXLeft', [], 'relErrXRight', [], 'relErrXOk', [], ...
            'relErrYLeft', [], 'relErrYRight', [], 'relErrYOk', [], ...
            'xAbsStop', [], 'xLinearResidual', [], 'innerResNorm', [], ...
            'beta', [], 'tau', [], 'alpha', [], ...
            'primalResidualRel', [], 'dualResidualRel', [], 'outerStopMetric', []);
    end
end