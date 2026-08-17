function out = lad(problem, options)
% IADMM.lad  Linearized exact IADMM for LAD.
%   问题形式：min μ‖x‖₁ + ‖Kx - c‖₁，
%   等价于 min μ‖x‖₁ + ‖y - c‖₁ s.t. Kx - y = 0.
%
%   文献符号对应：
%       A = K, B = -I, b = 0, f(x)=μ‖x‖₁, g(y)=‖y-c‖₁.
%
%   本函数不是 Xie Algorithm 2 的非精确相对误差实现，而是 LAD 的
%   线性化精确闭式版本。由于 x/y 子问题都用软阈值闭式求解，不使用
%   文献中的非精确误差变量 w_1, w_2, d_1, d_2，也不检查 (30)、(32)
%   的相对误差准则。
%
%   输入：
%       problem - 包含 .x0, .y0, .lambda0, .K, .c, .mu, .normK2
%       options    - 包含 .maxIterations, .tolerance
%   输出：
%       out - 结构体 .x, .y, .lambda, .solver, .history
%
%   x 子问题使用线性化二次项的精确近端解：
%       min_x μ‖x‖₁ + <g_x, x-x^k> + (ρ_x/2)‖x-x^k‖²,
%       g_x = K'*(β*(K*x^k-y^k)-λ^k)
%   因而
%       x^{k+1} = shrink(x^k - g_x/ρ_x, μ/ρ_x).
%   y 子问题为平移 L1 近端映射：
%       min_y ‖y-c‖₁ + (β/2)‖y-(K*x^{k+1}-λ^k/β)‖²,
%       y^{k+1} = c + shrink(K*x^{k+1}-λ^k/β-c, 1/β).


maxIterations = options.maxIterations;
tau     = 1;
beta    = 1;
startTime  = tic;

% 变量初始化
x      = problem.x0;
y      = problem.y0;
lambda = problem.lambda0;

rhoX = 1.05 * beta * problem.normK2;   % 线性化近端参数，需大于 beta*||K||^2
history = struct('obj', [], 'feas', [], 'time', [], 'err', [], 'inner', []);
K = problem.K;  % 为了后续代码简洁，提前定义 K
mu = problem.mu;  % 提前定义 μ
c = problem.c;  % 提前定义 c
for k = 1:maxIterations
    % Step 1: x 子问题。线性化 L_beta(x,y^k;lambda^k) 后精确软阈值求解。
    gradientX = K' * (beta * (K * x - y) - lambda);
    x = Functions.softThreshold(x - gradientX / rhoX, mu / rhoX);
    
    % Step 2: y 子问题。精确求解 L_beta(x^{k+1},y;lambda^k)。
    y = c + Functions.softThreshold(K * x - lambda / beta - c, 1 / beta);
    
    % Step 3: 乘子更新。此处 A*x+B*y-b = K*x-y。
    lambda = lambda - tau * beta * (K * x - y);
    
    % ---------- 记录历史 ----------
    history = Tools.appendHistory(history, problem, x, y, k, startTime);
    history.inner(k, 1) = 1;
end

out = Tools.packOut('IADMM.lad', x, y, lambda, history);
end
