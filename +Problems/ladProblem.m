function problem = ladProblem(p, q, s)
%LADPROBLEM 生成一个一般的凸LAD稀疏恢复问题。
%   PROBLEM = LADPROBLEM(P, Q, S) 创建问题
%
%       min_x  mu*||x||_1 + ||K*x - c||_1
%
%   等价拆分形式为：
%
%       min_{x,y}  mu*||x||_1 + ||y - c||_1,  满足  K*x - y = 0.
%
%   通过 CVX 计算最优目标值并存储在 PROBLEM.PhiStar 中。
%
%   输入：
%       p - 测量数（K的行数）
%       q - 变量数（x的维度）
%       s - 真实解的稀疏度（非零元个数）
%
%   输出：
%       problem - 包含问题数据、最优值、解等信息的结构体

    if nargin ~= 3
        error('ladProblem:InvalidInput', 'ladProblem需要三个输入：p, q, s');
    end
    validateattributes(p, {'numeric'}, {'scalar', 'integer', 'positive'});
    validateattributes(q, {'numeric'}, {'scalar', 'integer', 'positive'});
    validateattributes(s, {'numeric'}, {'scalar', 'integer', 'positive', '<=', q});

    rng(11);                        % 固定随机种子，保证结果可重复
    % 数值复现实验采用 1/sqrt(p) 尺度归一化。
    % GrpADMM 论文文字写作 K_ij ~ N(0,1)，但归一化后的结果与文献曲线更一致；
    % 这里保留该方案以避免 ||K||^2 过大导致线性化步长过小。
    K = randn(p, q) / sqrt(p);      % 随机测量矩阵

    % 构造稀疏真实解 xtrue
    xtrue = zeros(q, 1);
    indexSet = randperm(q, s);
    xtrue(indexSet) = randn(s, 1);
    clean = K * xtrue;
    c = clean;

    % 加入离群点（outliers）
    numBad = max(1, round(0.05 * p));
    bad = randperm(p, numBad);
    sgn = 2 * (rand(numBad, 1) > 0.5) - 1;
    c(bad) = sgn * norm(clean, inf);

    % 填充输出结构体
    problem.type = 'LAD';
    problem.matrixScaling = 'randn/sqrt(p)';
    problem.K = K;
    problem.H = K;
    problem.c = c;
    problem.xtrue = xtrue;
    problem.mu = 0.05;
    problem.mu1 = problem.mu;
    problem.mu2 = 0;
    problem.p = p;
    problem.q = q;
    problem.s = s;
    problem.x0 = zeros(q, 1);
    problem.y0 = zeros(p, 1);
    problem.lambda0 = -c;
    problem.normK = norm(K, 2);
    problem.normK2 = problem.normK^2;
    problem.normH = problem.normK;
    problem.normH2 = problem.normK2;
    problem.title = 'LAD Sparse Recovery';

    % 调用 CVX 计算最优目标值和解
    [problem.PhiStar, problem.xStar, problem.yStar, problem.cvxInfo] = ...
        solveLadPhiStarCvx(K, c, problem.mu);
end

%% ------------------------------------------------------------------------
%  辅助函数：使用 CVX 求解 LAD 问题（直接形式）
% ------------------------------------------------------------------------
function [phiStar, xStar, wStar, cvxInfo] = solveLadPhiStarCvx(K, c, mu)
%SOLVE_LAD_PHI_STAR_CVX 通过 CVX 计算 LAD 问题的最优目标值。
%   输入：
%       K  - 测量矩阵 (p×q)
%       c  - 观测向量 (p×1)
%       mu - L1 正则化参数
%   输出：
%       phiStar   - 最优目标值
%       xStar     - 最优原始变量 x
%       wStar     - 最优辅助变量 y = K*xStar
%       cvxInfo   - 求解器信息结构体


    [p, q] = size(K);

    % 使用 CVX 求解原始问题：min mu*norm(x,1) + norm(K*x - c,1)
    % 注意：CVX 默认使用 SeDuMi / SDPT3，可自动选择。
    cvx_begin quiet                         % quiet 模式关闭详细输出
        variable x(q)
        minimize( mu * norm(x,1) + norm(K*x - c,1) )
    cvx_end

    % 处理求解结果
    if strcmp(cvx_status, 'Solved') || strcmp(cvx_status, 'Optimal')
        xStar = x;
        wStar = K * xStar;                  % 满足等式约束 y = K*x
        phiStar = cvx_optval;               % 最优目标值
        cvxInfo = struct('solver', cvx_solver, ...
                         'status', cvx_status, ...
                         'precision', 'default');
    else
        error('ladProblem:CvxFailed', ...
            'CVX 求解失败，状态为: %s', cvx_status);
    end
end
