function H = walshHadamardMatrix(n)
%WALSH_HADAMARD_MATRIX 构建西尔维斯特型哈达玛矩阵。
%   H = WALSH_HADAMARD_MATRIX(N) 返回一个 N×N 的哈达玛矩阵，元素为 +1 或 -1。
%   N 必须是 2 的正整数次幂。
%
%   输入参数:
%       n - 矩阵阶数，限制为 2 的幂（如 1, 2, 4, 8, ...）。
%
%   输出参数:
%       H - N×N 的西尔维斯特构造的哈达玛矩阵。

% 检查 n 是否合法：必须大于等于 1 且是 2 的整数次幂。
% log2(n) 计算以 2 为底的对数，round(log2(n)) 四舍五入到最近整数，
% 若两者差的绝对值大于 eps（浮点精度），则 n 不是 2 的幂。
% 同时 n < 1 的情况也被拒绝。
if n < 1 || abs(log2(n) - round(log2(n))) > eps
    error('n must be a positive power of two.');
end

% 初始化 1×1 的哈达玛矩阵（定义 H_1 = 1）
H = 1;

% 通过递归克罗内克积构造更高阶矩阵。
% 当 H 的当前尺寸小于目标阶数 n 时，持续扩展：
%   已知 H_k 是 2^k 阶哈达玛矩阵，通过 [H_k  H_k; H_k  -H_k] 得到 2^{k+1} 阶矩阵。
% 该构造方法被称为西尔维斯特构造法，每次迭代将矩阵尺寸翻倍。
while size(H, 1) < n
    H = [H H; H -H]; %#ok<AGROW> 抑制 MATLAB 编辑器关于循环内数组增长的警告，这里是有意为之的预分配优化豁免。
end

% 循环结束后 H 即为 n×n 的哈达玛矩阵，其行（和列）相互正交，满足 H*H' = n*I。
end
