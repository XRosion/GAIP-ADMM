function out = lad(problem, options)
% PADMC.lad  使用 PADMC 文献格式求解 LAD 问题
%   LAD 拆分形式：
%       min mu*||x||_1 + ||y-c||_1,  s.t. K*x - y = 0.
%
%   对应 PADMC (4a)-(4d)：
%       A = K, B = -I, b = 0,
%       S_x = r*I - beta*K'*K, S_y = s_y*I.
%   输入：
%       problem - 包含字段 .x0, .y0, .lambda0, .K, .c, .mu, .normK2
%       options    - 包含字段 .maxIterations, .tolerance
%   输出：
%       out - 结构体，含 .x, .y, .lambda, .solver, .history

maxIterations = options.maxIterations;
startTime = tic;

% PADMC 文献数值实验参数
beta = 0.005;                         % 增广拉格朗日罚参数 beta
alpha_x = 2;                          % x 凸组合中心参数，文献实验取 2
alpha_y = 1;                          % y 凸组合中心参数
r = 1.001 * beta * problem.normK2;    % S_x = rI - beta*K'*K
s_y = 100;                            % S_y = s_y I

% 注：理论陈述中 alpha_x, alpha_y 属于 (0,2)。截图数值实验取 alpha_x=2，
% 这里按文献实验设置保留边界值。

% 变量初始化
x = problem.x0;
y = problem.y0;
lambda = problem.lambda0;
bar_x = x;
bar_y = y;
K = problem.K;
c = problem.c;
mu = problem.mu;

history = struct('obj', [], 'feas', [], 'time', [], 'err', []);

for k = 1:maxIterations
    % (4a) 凸组合近端中心
    bar_x = alpha_x * x + (1 - alpha_x) * bar_x;
    bar_y = alpha_y * y + (1 - alpha_y) * bar_y;

    % (4b) x 子问题：
    % S_x = rI - beta*K'*K 时，二次项合并为 r/2 ||x-v_x||^2。
    Kbarx = K * bar_x;
    v_x = bar_x - (K' * (lambda + beta * (Kbarx - y))) / r;
    x = Functions.softThreshold(v_x, mu / r);
    Kx = K * x;

    % (4c) y 子问题：S_y = s_y I
    a_y = beta + s_y;
    v_y = (beta * Kx + lambda + s_y * bar_y) / a_y;
    y = c + Functions.softThreshold(v_y - c, 1 / a_y);

    % (4d) 乘子更新，Ax+By-b = Kx-y
    lambda = lambda + beta * (Kx - y);

    % 记录历史
    history = Tools.appendHistory(history, problem, x, y, k, startTime);
end

out = Tools.packOut('PADMC.lad', x, y, lambda, history);
end
