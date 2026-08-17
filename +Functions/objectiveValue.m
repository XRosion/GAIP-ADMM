function value = objectiveValue(problem, x, y)
%OBJECTIVE_VALUE 计算分裂目标函数 Phi(x, y) 的值。
%   VAL = OBJECTIVE_VALUE(PROBLEM, X, W, VARIANT) 根据问题类型和变体，
%   返回当前分裂变量 (x, y) 对应的目标函数值。该值用于记录算法历史，
%   以及验证与 CPLEX 等精确求解器所得最优值的相对间隙。
%
%   输入参数：
%       problem - 结构体，必须包含字段 .type（'LAD'、'elastic' 或 'TV'），
%                 以及对应问题类型的正则化参数和数据。
%       x, y    - 分裂原始变量。y 是辅助变量，在 Phi(x, y) 表示中对应 y。
%
%   输出参数：
%       value   - 标量，目标函数值 Phi(x, y)。



    % 根据问题类型分别计算目标值
    switch problem.type

        case 'LAD'   % 最小绝对偏差（L1 正则化最小二乘）问题、
            % 标准 LAD: Phi(x,y)=mu*||x||_1 + ||y-c||_1。
            % 强凸 LAD 变体: Phi(x,y)=mu1*||x||_1 + mu2*||x||_2^2 + ||y-c||_1。
            if isfield(problem, 'mu2') && problem.mu2 > 0
	                value = problem.mu1 * norm(x, 1) + ...
	                      problem.mu2 * norm(x)^2 + ...
	                      norm(y - problem.c, 1);
	            else
	                value = problem.mu * norm(x, 1) + ...
	                      norm(y - problem.c, 1);
	            end
            

        case 'elastic'   % 弹性网问题（L1+L2 正则化 + 最小二乘拟合）
            % 目标函数：
            %   Φ(x, y) = γ * ||x||_2^2 + ½ * ||K x - b||_2^2
            %           + η * ||y||_1
            % 当前 elastic-net 采用与 LASSO 一致的 x-y=0 分裂：
            % 光滑二次块作用于 x，非光滑 L1 块作用于 y。
	            value = problem.gamma * sum(x.^2) + ...
	                  0.5 * norm(problem.K * x - problem.b)^2 + ...
	                  problem.eta * norm(y, 1);

        case 'lasso'     % LASSO 问题（x-y=0 分裂）
            % 当前 LASSO 采用分裂形式：
            %
            %   min_{x,y}  1/2 * ||K*x - b||_2^2 + eta * ||y||_1
            %   s.t.       x - y = 0.
            %
            % 因此目标函数为：
            %
            %   Phi(x,y) = 1/2 * ||K*x - b||_2^2 + eta * ||y||_1.
	            value = 0.5 * norm(problem.K * x - problem.b)^2 + ...
	                  problem.eta * norm(y, 1);

        case 'TV'        % 全变分（TV）正则化问题
            % 目标函数：
            %   Φ(x, y) = ½ * μ * ||x - c||_2^2 + Σ_g ||w_g||_2
            % 其中 y 是梯度场（通常 y = D x），第二项是各向同性 TV：
            %   group_norm_sum(y) = ∑_i sqrt( (D_x x)_i^2 + (D_y x)_i^2 )
	            value = 0.5 * problem.mu * norm(x - problem.c)^2 + ...  % 保真项（L2 距离）
	                  Functions.tvGroupNormSum(y);                      % 组 L2 范数之和（TV 项）
    end
end
