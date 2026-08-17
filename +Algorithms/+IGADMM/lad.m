function out = lad(problem, options)
% IGADMM.lad  广义 ADMM 格式求解 LAD 问题（闭式软阈值子问题）
%   问题：min μ‖x‖₁ + ‖Kx - c‖₁  <=>  min μ‖x‖₁ + ‖y-c‖₁  s.t.  Kx = y
%
%   输入：
%       problem - 包含 .K, .c, .mu, .x0, .y0, .lambda0, .normK2
%       options    - 包含 .maxIterations
%   输出：
%       out - 结构体 .x, .y, .lambda, .solver, .history

    %% ----- 算法参数（内部固定，参考论文推荐值）-----
    maxIterations = options.maxIterations; % 最大外迭代次数，可由 options 覆盖
    beta     = 1.0;           % 惩罚参数
    alpha    = -0.3;          % 广义 ADMM 参数，文献图像实验取 -0.3
    proximalParameter    = 0.01;          % T = proximalParameter * I  (T > 0)

    %% ----- 问题数据 -----
    K   = problem.K;
    c   = problem.c;
    mu  = problem.mu;
    x0  = problem.x0;
    y0  = problem.y0;        % 原辅助变量 y
    initialLambda= problem.lambda0;

    %% ----- 线性化近端参数 -----
    rhoX = 1.05 * beta * problem.normK2;   % 需大于 beta*||K||^2

    %% ----- 初始化 -----
    x   = x0;
    y   = y0;
    lambda = initialLambda;

    history = struct('obj', [], 'feas', [], 'time', [], 'inner', [], 'err', []);
    startTime = tic;

    %% ----- 主循环 -----
    for k = 1:maxIterations
        % ----- x^{k+1}：线性化 LAD 子问题的精确软阈值解 -----
        gradientX = beta * K' * (K * x - y - lambda / beta);
        xNew = Functions.softThreshold(x - gradientX / rhoX, mu / rhoX);
        KxNew = K * xNew;

        % ----- λ^{k+1/2} -----
        constraintResidual = KxNew - y;
        lambdaHalf = lambda - alpha * beta * constraintResidual;

        % ----- y^{k+1}：带 T=tI 近端项的平移软阈值 -----
        r = KxNew - lambdaHalf / beta;
        yCoefficient = beta + proximalParameter;
        yCenter = (beta * r + proximalParameter * y) / yCoefficient;
        yNew = c + Functions.softThreshold(yCenter - c, 1 / yCoefficient);

        % ----- 更新 λ^{k+1} -----
        newConstraintResidual = KxNew - yNew;
        lambdaNew = lambdaHalf - beta * newConstraintResidual;

        % ----- 记录历史，统一使用分裂目标 Phi(x,y) 和可行性残差 -----
        history = Tools.appendHistory(history, problem, xNew, yNew, k, startTime);
        history.inner(k, 1) = 1;

        if ~isfinite(history.obj(k)) || ~isfinite(history.feas(k))
            warning('IGADMM:NonFiniteHistory', ...
                'Non-finite value detected at iter %d; stopping early.', k);
            break;
        end

        % ----- 更新历史变量（为下一轮准备）-----
        x = xNew;           % x^{k+1}
        y = yNew;           % y^{k+1}
        lambda = lambdaNew;
    end

    out.x = x;
    out.y = y;
    out.lambda = lambda;
	    out.solver = 'IGADMMb';
	    out.history = history;
	end
