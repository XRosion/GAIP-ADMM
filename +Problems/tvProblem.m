function problem = tvProblem(imageName, options)
% TVPROBLEM 生成各向同性全变分图像去噪问题
%   PROBLEM = TVPROBLEM(IMAGENAME, OPTIONS) 构建如下形式的 TV 去噪问题：
%       min_x 0.5*mu*||x-c||_2^2 + sum_i ||(Dx)_i||_2
%
%   输入参数：
%       IMAGENAME : 字符串，图像名称（不含路径和扩展名，例如 'cameraman'）
%       OPTIONS   : 结构体，必须包含以下字段（无默认值，缺失将报错）：
%           imageSize  : [m, n] 目标图像尺寸（像素）
%           noiseStd   : 高斯噪声标准差
%           mu         : 保真项权重
%           seed       : 随机数种子
%
%   输出结构体 PROBLEM 包含算法所需的全部字段，并附加 CVX 参考解：
%       PhiStar     : 最优目标函数值（由 CVX 计算）
%       xStar       : 最优图像向量
%       yStar       : 对应梯度域变量 (D*xStar)
%       cvxInfo     : CVX 求解信息结构体

% 直接读取参数（无默认值，缺失字段将导致 MATLAB 报错）
m = options.imageSize(1);
n = options.imageSize(2);
noiseStd = options.noiseStd;
mu = options.mu;

% 设置随机种子
rng(options.seed);

% 图像路径：假定 Images 文件夹位于当前工作目录下
imgExts = {'.png', '.tif', '.jpg', '.bmp', '.pgm'};
imgPath = [];
for ext = imgExts
    fullPath = fullfile('Images', [imageName, ext{1}]);
    if exist(fullPath, 'file')
        imgPath = fullPath;
        break;
    end
end

% 读取或生成图像
if isempty(imgPath)
    warning('tvProblem:FileNotFound', ...
        '找不到图像文件 %s，使用合成图像替代。', imageName);
    raw = Tools.syntheticImage(m, n);
else
    raw = imread(imgPath);
    raw = double(raw) / 255;
    % 多通道取第一通道
    if size(raw, 3) > 1
        raw = raw(:, :, 1);
    end
    % 缩放到目标尺寸
    if size(raw,1) ~= m || size(raw,2) ~= n
        rr = round(linspace(1, size(raw,1), m));
        cc = round(linspace(1, size(raw,2), n));
        raw = raw(rr, cc);
    end
end

% 构建问题数据
xtrue = raw(:);
c = min(max(xtrue + noiseStd * randn(size(xtrue)), 0), 1);

% 输出结构体
problem.type = 'TV';
problem.m = m;
problem.n = n;
problem.N = m * n;
problem.xtrue = xtrue;
problem.c = c;
problem.pixelMax = 1;
problem.mu = mu;
problem.noiseStd = noiseStd;
problem.imageName = imageName;
problem.x0 = c;
problem.y0 = zeros(2 * m * n, 1);
problem.lambda0 = zeros(2 * m * n, 1);
problem.normH = sqrt(8);
problem.normH2 = 8;
problem.normD = problem.normH;
problem.normD2 = problem.normH2;
problem.constraintNorm = problem.normD;
problem.constraintNorm2 = problem.normD2;
% Algorithm 1 order: TV auxiliary block first, image block second.
problem.preStrongModulusF = 0;
problem.preStrongModulusG = mu;
problem.eigsDtD = Functions.tvLaplacianEigs(m, n);
problem.D = @(x) Functions.tvForward(x, m, n);
problem.DT = @(y) Functions.tvAdjoint(y, m, n);
problem.solveD = @(rhs, coeffI, coeffDtD) tvSolvePeriodic(rhs, coeffI, coeffDtD, problem);
problem.title = sprintf('TV denoising (%s)', imageName);

% 使用 CVX 计算参考最优解
[problem.PhiStar, problem.xStar, problem.yStar, problem.cvxInfo] = ...
    solveTvPhiStarCvx(m, n, mu, c);
end

% ========================================================================
% 子函数：周期边界下的快速求解 (I + c*D^T D) x = rhs
% ========================================================================
function x = tvSolvePeriodic(rightHandSide, coeffI, coeffDtD, problem)
    R = reshape(rightHandSide, problem.m, problem.n);
    den = coeffI + coeffDtD * problem.eigsDtD;
    x = real(ifft2(fft2(R) ./ den));
    x = x(:);
end

% ========================================================================
% 子函数：使用 CVX 计算 TV 去噪问题的精确最优解（参考值）
% ========================================================================
function [phiStar, xStar, yStar, cvxInfo] = solveTvPhiStarCvx(m, n, mu, c)
%SOLVETVPHISTARCVX 用 CVX 求解各向同性 TV 去噪，返回最优值及最优解。
%   输出：
%       phiStar : 最优目标函数值
%       xStar   : 最优图像向量（列）
%       yStar   : 最优梯度域变量，y = D*xStar，维度 2*m*n
%       cvxInfo : CVX 信息结构体

N = m * n;

cvx_begin quiet
    variable x(N)
    % 重塑为图像
    x_img = reshape(x, m, n);
    
    % 周期边界梯度
    % 水平梯度：x(i,j+1) - x(i,j)，列周期延拓
    Dx_h = x_img(:, [2:end, 1]) - x_img;
    % 垂直梯度：x(i+1,j) - x(i,j)，行周期延拓
    Dx_v = x_img([2:end, 1], :) - x_img;
    
    % 各向同性全变分项：逐像素 2-范数
    isotropic_norms = norms([Dx_h(:), Dx_v(:)], 2, 2);
    total_var = sum(isotropic_norms);
    
    % 目标：0.5*mu*||x - c||^2 + TV(x)
    minimize( 0.5 * mu * sum_square(x - c) + total_var )
cvx_end

if strcmp(cvx_status, 'Solved') || strcmp(cvx_status, 'Optimal')
    xStar = x;
    yStar = [Dx_h(:); Dx_v(:)];   % 与 problem.D 的输出顺序一致
    phiStar = cvx_optval;
    cvxInfo = struct('solver', cvx_solver, ...
                     'status', cvx_status, ...
                     'precision', 'default');
else
    error('tvProblem:CvxFailed', ...
        'CVX 无法求解 TV 参考问题。状态: %s', cvx_status);
end
end
