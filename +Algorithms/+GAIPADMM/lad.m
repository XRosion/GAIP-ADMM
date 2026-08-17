function out = lad(problem, options, parameters)
% GAIP_ADMM  求解 LAD 问题的广义加速近端 ADMM
%   输入：
%       problem - 包含 .K, .c, .mu, .x0, .y0, .lambda0, .normK2
%       options    - 必须包含以下字段：
%           .maxIterations
%       parameters   - 必须包含算法参数 .beta, .alpha, .theta_x, .theta_y,
%                 .theta_l, .eps3, .t_param, .rhoFactor。
%   输出：
%       out - 结构体 .x, .y, .lambda, .solver, .history

% ----- 直接从 parameters 读取算法参数（无默认值）-----
    maxIterations = options.maxIterations;
    beta    = parameters.beta;
    alpha   = parameters.alpha;
    theta_x = parameters.theta_x;
    theta_y = parameters.theta_y;
    theta_l = parameters.theta_l;
    epsilon3 = parameters.eps3;
    proximalParameter = parameters.t_param;

    rhoFactor = parameters.rhoFactor;
    
    % ----- 问题数据 -----
    K   = problem.K;
    c   = problem.c;
    mu  = problem.mu;
    x0  = problem.x0;
    initialY = problem.y0;
    initialLambda = problem.lambda0;
    normK2 = problem.normK2;
    
    % ----- 线性化步长（x子问题）-----
    maxEigenvalueS = 1;            % S=I 时最大特征值=1
    smoothLipschitz = beta * normK2 + (1/beta) * maxEigenvalueS;
    stepSize = 1 / (rhoFactor * smoothLipschitz);
    
    % ----- 初始化 -----
    x   = x0;
    y = initialY;
    lambda = initialLambda;
    x_bar_prev = x;
    y_bar_prev = y;
    lambda_bar_prev = lambda;
    
    % ----- 历史记录 -----
    history = struct('obj', [], 'feas', [], 'time', [], 'err', [], 'inner', []);
    startTime = tic;
    
    % ----- 主循环 -----
    for iteration = 1:maxIterations
        % 1. 广义加速点
        x_bar = (1 - theta_x) * x + theta_x * x_bar_prev;
        y_bar = (1 - theta_y) * y + theta_y * y_bar_prev;
        lambda_bar = (1 - theta_l) * lambda + theta_l * lambda_bar_prev;
        
        % 2. x-子问题：线性化 LAD 子问题的精确软阈值解
        gradientX = beta * K' * (K * x - y - lambda_bar / beta) ...
	               + (1 / beta) * (x - x_bar);
        x = Functions.softThreshold(x - stepSize * gradientX, mu * stepSize);
        Kx = K * x;
        
        % 3. y-子问题：带 T=(t/β)I 近端项的平移软阈值
        yCoefficient = beta + proximalParameter / beta;
        yCenter = (beta * Kx - lambda_bar + (proximalParameter / beta) * y_bar) / yCoefficient;
        y = c + Functions.softThreshold(yCenter - c, 1 / yCoefficient);
        
        % 4. 约束残差
        constraintResidual = Kx - y;
        
        % 5. 更新乘子
        lambda = lambda_bar - alpha * beta * constraintResidual;
        
        % 6. 记录历史，统一使用分裂目标 Phi(x,y) 和可行性残差
        history = Tools.appendHistory(history, problem, x, y, iteration, startTime);
        history.inner(iteration, 1) = 1;
        
        % 7. 更新加速历史
        x_bar_prev = x_bar;
        y_bar_prev = y_bar;
        lambda_bar_prev = lambda_bar;
        
    end
    
    % ----- 输出 -----
    out.x = x;
    out.y = y;
    out.lambda = lambda;
    out.solver = 'GAIP-ADMM';
    out.history = history;
end
