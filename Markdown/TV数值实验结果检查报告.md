# TV 数值实验结果检查报告

生成日期：2026-06-17  
状态说明：这是 2026-06-17 的历史检查记录。2026-06-18 后 TV 图像恢复程序已调整为“各算法内部非精确 CG 求解 x 子问题”，并移除了 `solveTvNormalEquation.m`、`solveTvQuadraticSubproblem.m` 以及额外的 IADMM 主程序入口；当前代码结构以 `Markdown/TV图像恢复程序整理说明.md` 为准。

检查对象：

- `Main4.m`
- `+Problems/tvProblem.m`
- `+Functions/tvForward.m`
- `+Functions/tvAdjoint.m`
- `+Functions/tvLaplacianEigs.m`
- `+Functions/solveTvNormalEquation.m`
- `+Functions/objectiveValue.m`
- `+Functions/feasibilityNorm.m`
- `+Algorithms/+AGAIPADMM/tv.m`
- `+Algorithms/+GAIPADMM/tv.m`
- `+Algorithms/+IADMM/tv.m`
- `+Algorithms/+IGADMM/tv.m`
- `TVResults/tv_experiment_summary.mat`
- `TVResults/cameraman_*.png`

## 1. 结论摘要

当前 TV 去噪实验的“怪”主要不是 TV 差分算子或 FFT 解线性系统写错造成的。静态检查和数值复核显示：

1. `tvForward`、`tvAdjoint`、`tvLaplacianEigs` 与 `solveTvNormalEquation` 在周期边界设定下是自洽的。用 Python 复核随机向量的伴随关系，得到相对误差约 `9.65e-15`。
2. 当前模型
   \[
   \min_x \frac{\mu}{2}\|x-c\|_2^2 + \mathrm{TV}(x)
   \]
   在 `mu = 11`、`noiseStd = 0.05` 时，本身会给出明显过平滑的解。一个普通精确 ADMM 基线也收敛到 `F(x) ≈ 2751.7`、`SNR ≈ 15.49 dB`，与当前结果非常接近。
3. 已保存的 `.mat` 结果和当前源码不完全一致。当前 `AGAIPADMM.tv` 和 `GAIPADMM.tv` 的 x 子问题使用 FFT 精确求解，代码中 `innerIter = 0`；但 `tv_experiment_summary.mat` 里 `AGAIP-ADMM_{ssp}` 的 `inner` 范围是 `3..79`，`GAIP-ADMM_{ssp}` 的 `inner` 范围是 `2..4`。说明现有图片/summary 很可能由旧版本源码生成，不能直接代表当前函数。
4. 目标函数图使用的是分裂目标 \(\Phi(x,y)\)，当 `D(x)-y` 还不小的时候，\(\Phi(x,y)\) 可能显著低估真实原问题目标 \(F(x)=\frac{\mu}{2}\|x-c\|^2+\mathrm{TV}(Dx)\)。因此 `GAIP/IADMM/IGADMM` 在 500 步未充分可行时，目标曲线和算法优劣判断有误导性。
5. 以当前 `tolerance = 1e-5` 判断，四个算法在保存结果里都没有真正达到容差。`AGAIP` 最接近，最终 `feas = 3.71e-05`；其余三个仍在 `7.96e-02` 到 `1.78e-01`。

因此，建议先不要把现象解释为“算法理论失效”。更合理的处理顺序是：重新用当前源码跑一次，改进指标记录，再对 `mu`、停止准则和比较公平性做系统扫描。

## 2. 实验结果复盘

实验参数来自 `Main4.m`：

| 参数 | 当前值 |
|---|---:|
| image | `cameraman` |
| imageSize | `256 x 256` |
| noiseStd | `0.05` |
| mu | `11` |
| maxIterations | `500` |
| tolerance | `1e-5` |
| seed | `33` |

从 `tv_experiment_summary.mat` 读取到的基础值：

| 对象 | 目标值/原问题目标 | SNR(dB) | 相对误差 | TV 项 | 数据项 |
|---|---:|---:|---:|---:|---:|
| noisy `c` | `7481.609` | `14.035` | `0.0924` | `7481.609` | `0` |
| true image | `4515.925` | 理论无穷大 | `0` | `3665.209` | `850.716` |

注意：真实干净图在该模型下的目标值 `4515.9` 高于算法最终约 `2751.7`。这意味着该正则化参数下，模型最优解本来就比真实图更倾向于低 TV、强平滑。

最终结果：

| 算法 | iter | time(s) | split obj \(\Phi(x,y)\) | 原问题目标 \(F(x)\) | feas | SNR(dB) | relerr |
|---|---:|---:|---:|---:|---:|---:|---:|
| AGAIP-ADMM\_{ssp} | 500 | `26.540` | `2751.687` | `2751.689` | `3.706e-05` | `15.488` | `0.07817` |
| GAIP-ADMM\_{ssp} | 500 | `4.067` | `2731.311` | `2761.813` | `1.713e-01` | `15.502` | `0.07804` |
| IADMM | 500 | `3.092` | `2743.724` | `2756.006` | `7.955e-02` | `15.498` | `0.07808` |
| IGADMM | 500 | `2.879` | `2730.391` | `2762.283` | `1.777e-01` | `15.503` | `0.07804` |

可行性阈值首次达到的迭代数：

| 算法 | `1e-1` | `1e-2` | `1e-3` | `1e-4` | `1e-5` |
|---|---:|---:|---:|---:|---:|
| AGAIP-ADMM\_{ssp} | 31 | 65 | 143 | 332 | 未达到 |
| GAIP-ADMM\_{ssp} | 未达到 | 未达到 | 未达到 | 未达到 | 未达到 |
| IADMM | 423 | 未达到 | 未达到 | 未达到 | 未达到 |
| IGADMM | 未达到 | 未达到 | 未达到 | 未达到 | 未达到 |

SNR 峰值也值得注意：

| 算法 | 最终 SNR | 最大 SNR | 最大 SNR 位置 |
|---|---:|---:|---:|
| AGAIP-ADMM\_{ssp} | `15.488` | `16.392` | iter 1 |
| GAIP-ADMM\_{ssp} | `15.502` | `16.489` | iter 3 |
| IADMM | `15.498` | `16.263` | iter 2 |
| IGADMM | `15.503` | `16.460` | iter 2 |

这说明 SNR 不是随着优化目标继续改善而改善；继续降低 TV 目标会让图像更平滑，反而牺牲恢复质量。

## 3. 函数正确性检查

### 3.1 TV 差分和伴随算子

`tvForward.m` 使用周期前向差分：

\[
D_v x(i,j)=x(i+1,j)-x(i,j),\quad
D_h x(i,j)=x(i,j+1)-x(i,j).
\]

`tvAdjoint.m` 使用：

\[
D_v^\top p = \mathrm{circshift}(p,1,1)-p,\quad
D_h^\top q = \mathrm{circshift}(q,1,2)-q.
\]

这与前向差分互为伴随。数值复核：

\[
\frac{|\langle Dx,y\rangle-\langle x,D^\top y\rangle|}{\max(1,|\langle Dx,y\rangle|,|\langle x,D^\top y\rangle|)}
\approx 9.65\times 10^{-15}.
\]

结论：没有发现符号或堆叠顺序错误。

### 3.2 \(D^\top D\) 特征值和 FFT 解法

`tvLaplacianEigs.m` 使用

\[
\lambda_{ij}=4-2\cos(2\pi i/m)-2\cos(2\pi j/n),
\]

这对应周期前向差分的 \(D^\top D\)。`solveTvNormalEquation.m` 解

\[
(aI+bD^\top D)x = rhs
\]

时使用 `fft2(rhs)./(a+b*eigsDtD)`。在 `a>0,b>=0` 时分母正，和当前模型一致。

结论：FFT 正规方程求解器形式正确。

### 3.3 TV 目标函数和 feasibility

`objectiveValue.m` 对 TV 记录：

\[
\Phi(x,y)=\frac{\mu}{2}\|x-c\|^2+\sum_i\sqrt{y_{1,i}^2+y_{2,i}^2}.
\]

`feasibilityNorm.m` 记录：

\[
\|D(x)-y\|_2.
\]

这两个函数对分裂模型是正确的。但在算法还不可行时，\(\Phi(x,y)\) 不是原问题目标 \(F(x)\)。例如最终：

- `GAIP`: \(\Phi=2731.31\)，但 \(F(x)=2761.81\)，相差 `30.50`。
- `IGADMM`: \(\Phi=2730.39\)，但 \(F(x)=2762.28\)，相差 `31.89`。

因此当前 objective 图不适合单独判断恢复质量或原问题收敛。

### 3.4 组软阈值

`groupSoftThreshold.m` 对 `[v1; v2]` 做逐像素二维块软阈值：

\[
\max(1-\tau/\sqrt{v_{1,i}^2+v_{2,i}^2},0)(v_{1,i},v_{2,i}).
\]

这是各向同性 TV 的正确 prox。

### 3.5 ADMM 乘子符号

代码统一采用如下符号约定：

\[
L=f(x)+g(y)-\langle\lambda,Dx-y\rangle+\frac{\beta}{2}\|Dx-y\|^2.
\]

在这个约定下：

- x 子问题右端出现 `DT(lambda + beta*y)`；
- y 子问题中心是 `Dx - lambda/beta`；
- 乘子更新是 `lambda = lambda - beta*(Dx-y)`。

这套符号内部一致。若论文公式使用 \(+\langle\lambda,Dx-y\rangle\)，则代码里的 `lambda` 等价于论文乘子的相反数。报告公式和代码对照时必须统一这个符号，否则容易误判。

## 4. 主要可疑点和风险

### 4.1 保存结果与当前源码不一致

当前源码：

- `+Algorithms/+AGAIPADMM/tv.m` 第 91 行使用 FFT 解 x 子问题；
- 第 101 行设置 `innerIter = 0`；
- `+Algorithms/+GAIPADMM/tv.m` 第 111 行使用 FFT 解 x 子问题；
- 第 120 行设置 `innerIter = 0`。

但 `.mat` 中：

| 算法 | `.mat` 中 inner 范围 | 当前源码预期 |
|---|---:|---:|
| AGAIP-ADMM\_{ssp} | `3..79` | 全部为 `0` |
| GAIP-ADMM\_{ssp} | `2..4` | 全部为 `0` |

这表明图片和 summary 很可能由旧版本函数生成。任何基于这些图的结论都需要先重新运行当前 `Main4.m`。

### 4.2 `mu = 11` 导致过平滑

用同一数据和普通精确 ADMM 做近似扫描：

| mu | 近似 SNR(dB) | relerr | 目标值 | TV 项 | 数据项 |
|---:|---:|---:|---:|---:|---:|
| 2 | `9.207` | `0.1611` | `1185.4` | `575.1` | `610.3` |
| 5 | `12.279` | `0.1131` | `1881.5` | `955.9` | `925.6` |
| 8 | `14.186` | `0.0908` | `2363.3` | `1231.0` | `1132.3` |
| 11 | `15.488` | `0.0782` | `2751.8` | `1437.4` | `1314.3` |
| 15 | `16.741` | `0.0677` | `3191.6` | `1665.3` | `1526.2` |
| 20 | `17.871` | `0.0594` | `3660.5` | `1928.5` | `1732.0` |
| 30 | `19.006` | `0.0521` | `4405.8` | `2505.6` | `1900.2` |
| 50 | `18.570` | `0.0548` | `5340.4` | `3652.8` | `1687.5` |
| 80 | `17.233` | `0.0639` | `6037.8` | `4767.6` | `1270.3` |
| 120 | `16.234` | `0.0717` | `6480.9` | `5555.8` | `925.2` |
| 200 | `15.369` | `0.0793` | `6863.2` | `6271.6` | `591.7` |

当前 `mu=11` 明显偏向低 TV 解。若实验目标是图像恢复质量，建议先在 `mu=20..50` 做精细扫描；`mu≈30` 在这次近似测试中 SNR 最高。

### 4.3 目标曲线不应只画 \(|\Phi(x,y)|\)

`plotHistories.m` 在没有 `PhiStar` 时画的是 `abs(history.obj)`，纵轴标签是 `|\Phi(x^{k+1},y^{k+1})|`。对 TV 来说：

1. 该值是分裂目标，不是原问题目标；
2. 它越小不一定图像越好；
3. 在不可行时，不同算法的 \(\Phi(x,y)\) 可比性较弱；
4. 目前图中目标从很小值上升到 2700 多，容易被误读成“目标变差”，其实是 `y` 从零场逐步接近 `D(x)` 后 TV 项被补回来。

建议同时记录并绘制：

- `history.obj_split = Phi(x,y)`;
- `history.obj_primal = Phi(x,Dx)`;
- `history.dataTerm = 0.5*mu*norm(x-c)^2`;
- `history.tvDx = TV(Dx)`;
- `history.tvY = TV(y)`;
- `history.feas = norm(Dx-y)`;
- `history.snr` 和 `history.err`。

### 4.4 只用 feasibility 停止不够

`Main4.m` 用 `finalFeasibility <= tolerance` 判断收敛。对分裂 TV 问题，仅 primal feasibility 不足以判断 ADMM 收敛，还应包含：

- primal residual：\(\|Dx-y\|\)；
- dual residual：例如 \(\beta D^\top(y^k-y^{k-1})\) 的范数；
- 目标变化或 KKT 残差；
- 对图像恢复实验，可额外记录 best SNR iteration。

当前结果中 `AGAIP` feasibility 很低，但 SNR 从早期峰值下降；这正说明优化收敛和恢复指标不等价。

### 4.5 `IADMM.tv` 和 `IGADMM.tv` 的 x 子问题非常粗

`.mat` 中 IADMM/IGADMM 的平均内迭代数只有 `1.028` 和 `0.762`，大多数步只做 0 到 2 次 CG。相对误差条件全部满足，但从公平比较角度看，当前设定等价于很粗的非精确 x 步，和 AGAIP/GAIP 的 FFT 精确步不是同一计算精度。

建议至少增加一个精确 ADMM 或 FFT-ADMM 基线，并在图中明确区分：

- 每步精确 x 解；
- 相对误差非精确 x 解；
- 固定 CG 次数/固定残差容差。

## 5. 建议修改

### 5.1 先重新生成结果

因为保存结果与当前源码不一致，建议先在 MATLAB 中重新运行：

```matlab
Main4
```

然后检查 `tv_experiment_summary.mat` 中：

```matlab
allResults.cameraman.results{1}.history.inner
allResults.cameraman.results{2}.history.inner
```

若当前源码生效，`AGAIP` 和 `GAIP` 应全部为 0。

### 5.2 增加原问题目标记录

建议在 `+Tools/appendHistory.m` 的 TV 分支增加：

```matlab
if strcmp(problem.type, 'TV')
    Dx = problem.D(x);
    history.obj_primal(iteration, 1) = Functions.objectiveValue(problem, x, Dx);
    history.obj_split(iteration, 1) = phi;
    history.dataTerm(iteration, 1) = 0.5 * problem.mu * norm(x - problem.c)^2;
    history.tvDx(iteration, 1) = tvGroupNormSum(Dx);
    history.tvY(iteration, 1) = tvGroupNormSum(y);
end
```

其中 `tvGroupNormSum` 可以复用 `objectiveValue.m` 中的 TV 组范数逻辑，最好抽成公开工具函数，避免重复实现。

### 5.3 增加 dual residual 和 KKT 残差

对 ADMM 类算法，建议每步记录：

```matlab
history.primalResidual(iteration,1) = norm(Dx - y);
history.dualResidual(iteration,1) = beta * norm(problem.DT(y - yPrevious));
```

若 `beta` 是变参数，则用当前步对应的 `betaCurrent`。对带近端项或广义松弛的算法，可另外记录算法自己的理论残差，但至少应有一个统一的可比较 residual。

### 5.4 扫描 `mu`

当前 `mu=11` 对 `noiseStd=0.05` 的 cameraman 过平滑。建议新建一个敏感性实验：

```matlab
muGrid = [15 20 25 30 35 40 50];
```

每个 `mu` 记录：

- final SNR；
- best SNR；
- best SNR iteration；
- final feasibility；
- final primal objective；
- restored image。

如果论文实验需要固定某个文献参数，也建议把“文献参数”和“恢复质量最优参数”分开报告。

### 5.5 绘图建议

建议替换/补充现有图：

1. `objective_primal_iter.png`：画 \(F(x^k)\)，不要只画 \(\Phi(x^k,y^k)\)。
2. `objective_split_iter.png`：保留分裂目标，但标题/标签说明是 split objective。
3. `residuals_iter.png`：同时画 primal/dual residual。
4. `snr_iter.png`：增加 best SNR 标记。
5. `mu_scan_snr.png`：展示参数选择对恢复质量的影响。

## 6. 对当前图片的解释

### feasibility 图

`AGAIP-ADMM_{ssp}` 可行性下降最快，最终到 `3.7e-05`，但仍未达到 `1e-5`。其他三个算法 500 步时仍在 `1e-1` 量级附近，因此不应按 `tolerance=1e-5` 宣称收敛。

### objective 图

目标曲线从很小值上升到约 2700，不代表算法在最大化目标。主要原因是初始 `y=0`，split TV 项很小；随着 `y` 接近 `D(x)`，TV 项被恢复，split objective 上升。应改看原问题目标 \(F(x)\) 和 residual。

### recovery comparison 图

恢复图偏平滑，与 `mu=11` 下 TV 模型最优解一致。近似 `mu` 扫描显示，把 `mu` 增大到 `20..50` 会显著提高 SNR，`mu≈30` 附近在本次检查中达到约 `19.0 dB`。

### SNR 图

SNR 在前几步达到峰值后下降并稳定，说明“优化目标继续下降/可行性继续改善”和“接近真实图像”不是同一个目标。论文若报告图像恢复质量，建议报告 best SNR 和 final SNR 两列。

## 7. 建议下一步实验清单

1. 用当前源码重新运行 `Main4.m`，确认 `.mat` 和源码一致。
2. 加入标准精确 ADMM/FFT-ADMM 基线，作为正确性锚点。
3. 在 history 中同时记录 split objective 和 primal objective。
4. 加入 dual residual。
5. 扫描 `mu`，先定位恢复质量合理区间。
6. 对每个算法统一停止到同一 residual 标准，而不是固定 500 步后直接比较。
7. 若要比较 CPU 时间，确保 MATLAB 环境、FFT/CG 精度、history 记录量一致，并多次运行取中位数。

## 8. 总体判断

未发现 TV 基础函数的确定性错误。当前数值效果奇怪，主要来自三件事叠加：

1. 保存结果很可能不是当前源码生成；
2. `mu=11` 对当前噪声水平和图像来说偏向过平滑；
3. objective 图画的是 split objective，在不可行阶段非常容易误导。

优先修复记录和绘图口径，再重新跑实验。若重新运行当前源码后仍出现 AGAIP/GAIP 的 `inner` 非零、或 FFT 解残差不接近机器精度，再继续检查 MATLAB 路径污染或旧函数遮蔽问题。
