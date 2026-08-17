function history = appendHistory(history, problem, x, y, iteration, startTime, options)
%APPEND_HISTORY 将一次迭代的记录追加到算法历史结构中。
%   HISTORY = APPEND_HISTORY(HISTORY, PROBLEM, X, Y, ITERATION, STARTTIME)
%   计算当前迭代的目标函数值 Phi(x^k, y^k)、可行性残差、已耗CPU时间，
%   并在参考数据可用时（如最优值、真实图像）额外计算信号恢复诊断指标。
%
%   输入参数：
%       history    - 由算法维护的历史结构体。函数会向其中追加 obj,
%                 feas, time, err, psnr 等字段（只要所需参考数据存在）。
%       problem - 问题结构体。若包含 problem.PhiStar，则计算相对最优性误差 err。
%       x, y    - 当前迭代的原始分裂变量（通常为主变量和辅助变量）。
%       iteration - 正整数迭代索引。
%       startTime - 算法循环开始前由 tic 创建的计时器句柄。
%
%   输出参数：
%       history    - 更新后的历史结构体，包含本次迭代的各项记录。
%
%   注意事项：
%       - 函数不会检查输入维度的合法性，假定调用方已保证正确性。
%       - 若 problem 中缺少某些参考字段（PhiStar, xtrue, type），
%         则对应的诊断字段不会被创建或追加。
%       - 所有分母均使用 max(denominator, eps) 避免除零。

    if nargin < 7
        options = struct();
    end

    storeFinalOnly = Tools.keepOnlyFinalVisualizationData(options);
    if storeFinalOnly && isfield(options, 'maxIterations') && ...
            iteration < options.maxIterations
        return;
    end

    % 1. 计算当前目标函数值 Phi(x,y)
    %    Functions.objectiveValue 根据 problem 和 variant 返回具体的目标值
    phi = Functions.objectiveValue(problem, x, y);
    
    if storeFinalOnly
        storageIndex = 1;
    else
        storageIndex = iteration;
    end

    % 2. 记录目标函数值（列向量，行数为迭代次数）
    history.obj(storageIndex, 1) = phi;
    
    % 3. 计算可行性残差（例如约束违反的范数）
    %    feasibilityNorm 返回非负标量，衡量当前点满足约束的程度
    history.feas(storageIndex, 1) = Functions.feasibilityNorm(problem, x, y);

    % 4. 若提供了最优值 PhiStar，则计算相对最优性误差（归一化）
    if isfield(problem, 'PhiStar') && isfinite(problem.PhiStar)
        if strcmp(problem.type, 'TV')
            Dx = problem.D(x);
            primalObjective = 0.5 * problem.mu * norm(x - problem.c)^2 + ...
                Functions.tvGroupNormSum(Dx);
        else
            primalObjective = phi;
        end
        history.err(storageIndex, 1) = abs(primalObjective - problem.PhiStar) / ...
            max(abs(problem.PhiStar), eps);
    end

    if Tools.keepOnlyVisualizationData(options)
        return;
    end
    
    % 5. 记录从 startTime 开始已流逝的 CPU 时间（秒）
    history.time(iteration, 1) = toc(startTime);

    % 6. 针对 'TV'（全变分）类型问题，记录原问题目标和 PSNR
    if strcmp(problem.type, 'TV')
        if ~exist('Dx', 'var') || isempty(Dx)
            Dx = problem.D(x);
        end
        history.objSplit(iteration, 1) = phi;
        history.objPrimal(iteration, 1) = Functions.objectiveValue(problem, x, Dx);
        history.dataTerm(iteration, 1) = 0.5 * problem.mu * norm(x - problem.c)^2;
        history.tvDx(iteration, 1) = Functions.tvGroupNormSum(Dx);
        history.tvY(iteration, 1) = Functions.tvGroupNormSum(y);

        if isfield(problem, 'xtrue')
            history.psnr(iteration, 1) = Functions.psnrValue(problem, x);
        end
    end
    
end
