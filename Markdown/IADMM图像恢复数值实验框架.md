# IADMM 文献图像恢复数值实验框架整理

依据文献：Jiaxin Xie, *On inexact ADMMs with relative error criteria*, Computational Optimization and Applications, 2018。  
本整理重点对应文献第 3 节 Algorithm 1 和第 4.1 节 TV-based image denoising。

## 1. 核心定位

文献的图像恢复实验不是把 x 子问题精确求解后当作 IADMM 特例，而是专门展示：

1. x 子问题是大规模线性系统，用 CG 非精确求解；
2. CG 内层不以固定高精度残差为主要标准，而是检查相对误差准则；
3. y 子问题可以用 shrinkage 精确求解；
4. 通过参数 `alpha` 牺牲一部分 x 子问题求解精度，换取更大的乘子步长 `tau`；
5. 最终目标是用更少 CPU 时间和更少累计内迭代数，达到与精确 ADMM 近似相同的 PSNR。

因此，如果按你的算法逻辑复现该文献，`IADMM.tv` 的 x 子问题应该保留 CG 非精确求解，而不是直接改成 FFT/Cholesky 精确解。FFT 或 Cholesky 可以作为精确 ADMM 对照组。

## 2. 文献模型

文献先考虑 analysis sparse recovery 形式：

$$
\min_x \|Dx\|_1 + \frac{1}{2\mu}\|Ax-c\|_2^2.
$$

引入辅助变量 \(y=Dx\)，得到分裂模型：

$$
\min_{x,y}\left\{
\|y\|_1 + \frac{1}{2\mu}\|Ax-c\|_2^2
\;:\;
Dx-y=0
\right\}.
$$

在 ADMM 统一形式中：

$$
f(x)=\frac{1}{2\mu}\|Ax-c\|_2^2,\quad
g(y)=\|y\|_1,\quad
\mathcal A=D,\quad
\mathcal B=-I,\quad
b=0.
$$

对 TV 图像去噪，文献使用：

$$
\min_x \|x\|_{\mathrm{TV}} + \frac{1}{2\mu}\|x-c\|_2^2.
$$

其中图像向量化为 \(x\in\mathbb R^{mn}\)，TV 写作：

$$
\|x\|_{\mathrm{TV}}
= \sum_{i,j}\|\nabla x_{i,j}\|_2,
\quad
\nabla x_{i,j}
=
\begin{bmatrix}
(D_1x)_{i,j}\\
(D_2x)_{i,j}
\end{bmatrix}.
$$

差分算子为：

$$
(D_1x)_{i,j}=x_{i+1,j}-x_{i,j},
\quad
(D_2x)_{i,j}=x_{i,j+1}-x_{i,j}.
$$

记 \(D=[D_1,D_2]^T\)，则 \(\|x\|_{\mathrm{TV}}=\|Dx\|\)，这里的 \(\|Dx\|\) 是逐像素二维梯度块范数求和。

## 3. 与当前代码模型的参数换算

你的当前 TV 代码模型通常写成：

$$
\min_x \frac{\mu_{\mathrm{code}}}{2}\|x-c\|_2^2+\mathrm{TV}(x).
$$

文献写成：

$$
\min_x \mathrm{TV}(x)+\frac{1}{2\mu_{\mathrm{paper}}}\|x-c\|_2^2.
$$

二者等价关系是：

$$
\mu_{\mathrm{code}}=\frac{1}{\mu_{\mathrm{paper}}}.
$$

例如：

| 文献参数 | 对应当前代码参数 |
|---:|---:|
| `mu_paper = 0.05` | `mu_code = 20` |
| `mu_paper = 0.125` | `mu_code = 8` |
| `mu_code = 11` | `mu_paper ≈ 0.0909` |

这个换算很重要。若直接把文献中的 `mu=0.05` 填到你的 `options.mu`，模型强度会完全反过来。

## 4. Algorithm 1: 单个子问题非精确的 IADMM

文献图像实验主要采用 Algorithm 1，即 y 子问题精确、x 子问题非精确。

给定参数：

$$
\beta>0,\quad
\sigma\in[0,1),\quad
\alpha\in(0,1),\quad
\gamma=1-\alpha.
$$

乘子步长满足：

$$
0<\tau<
\frac{2\gamma-1+\sqrt{4\gamma^2+1}}{2\gamma}.
$$

每轮外迭代：

### 4.1 x 子问题非精确求解

计算

$$
x^{k+1}\approx
\arg\min_x L_\beta(x,y^k;\lambda^k),
$$

但不要求精确最小化，只要求候选点 \(x^{k+1}\) 和误差向量 \(d^{k+1}\) 满足：

$$
\frac{2}{\beta}
\left|\langle w^k-x^{k+1},d^{k+1}\rangle\right|
+\|d^{k+1}\|_2^2
\le
\sigma\alpha\|Dx^{k+1}-y^k\|_2^2.
$$

其中

$$
d^{k+1}
\in
\partial f(x^{k+1})
-D^T\left[
\lambda^k-\beta(Dx^{k+1}-y^k)
\right].
$$

对 TV 去噪，\(f(x)\) 是二次函数，因此 \(d^{k+1}\) 就是 x 子问题一阶最优性残差。

按你的代码模型

$$
f(x)=\frac{\mu_{\mathrm{code}}}{2}\|x-c\|_2^2,
$$

可写为：

$$
d^{k+1}
=
\mu_{\mathrm{code}}(x^{k+1}-c)
-D^T\lambda^k
+\beta D^T(Dx^{k+1}-y^k).
$$

因此 x 子问题的线性系统是：

$$
(\mu_{\mathrm{code}}I+\beta D^TD)x^{k+1}
=
\mu_{\mathrm{code}}c+D^T(\lambda^k+\beta y^k).
$$

文献模型下等价写法是：

$$
(A^TA+\mu_{\mathrm{paper}}\beta D^TD)x^{k+1}
=
\mu_{\mathrm{paper}}D^T(\lambda^k+\beta y^k)+A^Tc.
$$

图像去噪时 \(A=I\)。

### 4.2 y 子问题精确求解

y 子问题是：

$$
y^{k+1}
=
\arg\min_y L_\beta(x^{k+1},y;\lambda^k).
$$

它有 shrinkage 闭式解。各向同性 TV 时应使用二维组 shrinkage：

$$
y_i^{k+1}
=
\max\left(1-\frac{1/\beta}{\|v_i\|_2},0\right)v_i,
\quad
v=Dx^{k+1}-\lambda^k/\beta.
$$

若是一维 \(\ell_1\) analysis 模型，则是逐分量 soft shrinkage。

### 4.3 乘子和辅助变量更新

文献采用的符号约定为：

$$
\lambda^{k+1}
=
\lambda^k-\tau\beta(Dx^{k+1}-y^{k+1}),
$$

$$
w^{k+1}=w^k-\beta d^{k+1}.
$$

注意：这与你现有代码中的 `lambda = lambda - tau*beta*residual` 方向一致。

## 5. CG 内层策略

文献对 x 子问题采用 CG，核心不是每次都解到高精度，而是：

1. 外层第 \(k+1\) 次求 x 时，用上一轮 \(x^k\) 作为 CG 初值；
2. 每做若干 CG 步后计算当前候选 \(x\) 的一阶残差 \(d\)；
3. 检查相对误差准则是否成立；
4. 一旦成立，就停止 CG，并接受该 x 候选；
5. 记录本次 CG 迭代数，累计为 `Inniter`。

文献也设置了精确 ADMM 对照：

- `AD-CG`：仍用 CG 解 x 线性系统，但解到高精度；
- `AD-Chol`：用稀疏 Cholesky 分解精确求解 x 线性系统；
- `iAD-CG.4`：使用 Algorithm 1，`alpha=0.4`，x 子问题按相对误差准则非精确求解；
- `iAD-CG1`：单位步长的非精确 ADMM，对应 `alpha=1`、`tau=1`。

精确 ADMM 中 CG 的高精度停止准则为：

$$
\frac{\|Cx-b\|_2}{\max\{\|b\|_2,1\}}<10^{-8},
$$

其中

$$
C=A^TA+\mu_{\mathrm{paper}}\beta D^TD,
\quad
b=\mu_{\mathrm{paper}}D^T(\lambda^k+\beta y^k)+A^Tc.
$$

但对 `iAD-CG.4`，真正体现算法逻辑的是相对误差准则，不是上面这个高精度准则。

## 6. `alpha` 和 `tau` 的实验逻辑

Algorithm 1 中：

$$
\gamma=1-\alpha,
\quad
\tau_{\max}
=
\frac{2\gamma-1+\sqrt{4\gamma^2+1}}{2\gamma}.
$$

文献实际取 \(\tau\) 为略小于 \(\tau_{\max}\) 的三位小数。若 \(1000\tau_{\max}\) 不是整数，则取：

$$
\tau=\frac{\lfloor1000\tau_{\max}\rfloor}{1000}.
$$

如果刚好是整数，则再减去 `0.001`。

典型例子：

| alpha | gamma | tau |
|---:|---:|---:|
| `0.4` | `0.6` | `1.468` |
| `1` | `0` | `1.000`，按文献作为单位步长旧式非精确 ADMM |

解释：

- `alpha` 太小：允许的 `tau` 上界更大，但误差准则右侧有 `sigma*alpha`，会变得更严格，CG 内迭代很多；
- `alpha=1`：误差容忍相对宽，但只能用 `tau=1`；
- `alpha=0.4`：文献实验中 CPU 时间最低，是折中点。

## 7. 图像实验流程

文献第 4.1 节的 TV 图像恢复实验流程如下：

1. 选取测试图像 `Pirate` 和 `Moon`。
2. 将像素缩放到 `[0,1]`。
3. 加入零均值白高斯噪声，噪声水平记为 `delta`。
4. 初始化：

   $$
   x^0=c,\quad y^0=0,\quad \lambda^0=0,\quad w^0=0.
   $$

5. 外层迭代中：
   - 用 CG 非精确求解 x 线性系统；
   - CG 初值用上一轮外迭代的 \(x^k\)；
   - 用相对误差准则决定何时停止 CG；
   - 用 shrinkage 精确更新 y；
   - 用较大步长更新乘子；
   - 更新辅助变量 w。
6. 手动调节 `mu_paper` 以获得较高 PSNR。
7. 对 `beta` 试若干正值，选取表现较好的值。
8. 记录外层迭代数、累计内层 CG 次数、CPU 时间和 PSNR。

## 8. 外层停止准则

文献使用 primal residual 和 dual residual 的归一化最大值：

$$
\max\left\{
\frac{\|Dx^k-y^k\|_2}
{\max\{\|Dx^k\|_2,\|y^k\|_2,1\}},
\frac{\|\beta D^T(y^{k+1}-y^k)\|_2}
{\max\{1,\|D^T\lambda^k\|_2\}}
\right\}
<\eta.
$$

实践中更自然的实现方式是在得到 \(x^{k+1},y^{k+1}\) 后记录：

$$
r_{\mathrm{pri}}^{k+1}
=
\frac{\|Dx^{k+1}-y^{k+1}\|_2}
{\max\{\|Dx^{k+1}\|_2,\|y^{k+1}\|_2,1\}},
$$

$$
r_{\mathrm{dual}}^{k+1}
=
\frac{\|\beta D^T(y^{k+1}-y^k)\|_2}
{\max\{1,\|D^T\lambda^k\|_2\}}.
$$

然后用

$$
\max\{r_{\mathrm{pri}}^{k+1},r_{\mathrm{dual}}^{k+1}\}<\eta
$$

作为停止条件。

这比只看 \(\|Dx-y\|\) 更符合文献实验。

## 9. 图像质量指标

文献使用 PSNR：

$$
\mathrm{PSNR}
=
10\log_{10}\left(\frac{x_{\max}^2}{\mathrm{MSE}}\right),
$$

其中

$$
\mathrm{MSE}
=
\frac{1}{mn}\sum_{i=1}^m\sum_{j=1}^n
(\hat x_{i,j}-x_{i,j})^2.
$$

若图像已缩放到 `[0,1]`，则 \(x_{\max}=1\)。

你的代码里当前使用的是基于图像均值中心化信号能量的 SNR。若要严格复现文献，应额外记录 PSNR，而不是只记录 SNR。

## 10. 文献具体实验设置和结果

### 10.1 alpha 扫描实验

目的：展示 `alpha` 控制误差容忍和步长上界的折中。

设置：

| 项目 | Pirate | Moon |
|---|---:|---:|
| noise level `delta` | `0.01` | `0.01` |
| `mu_paper` | `0.05` | `0.05` |
| 对应 `mu_code` | `20` | `20` |
| beta | `20` | `60` |
| alpha grid | `0.005, 0.1, 0.2, ..., 1` | 同左 |

观察：

- CPU 时间与累计内层迭代数高度相关；
- `alpha=0.005` 太小，误差准则过严，内层迭代多；
- `alpha=0.4` 在 Pirate 和 Moon 上 CPU 时间最低；
- `alpha=0.1..0.9` 多数比 `alpha=1` 更快，说明较大乘子步长确实有收益。

### 10.2 表 1 对比实验

| Image | Algorithm | Step length | Iter | Inniter | Time(s) | PSNR |
|---|---|---:|---:|---:|---:|---:|
| Pirate | AD-Chol | 1.618 | 68 | - | 49.47 | 27.40 |
| Pirate | AD-CG | 1.618 | 69 | 967 | 9.52 | 27.40 |
| Pirate | iAD-CG.4 | 1.468 | 75 | 357 | 8.10 | 27.40 |
| Pirate | iAD-CG1 | 1.000 | 110 | 403 | 10.18 | 27.40 |
| Pirate | iAD-CG.4 | 1.000 | 110 | 466 | 10.91 | 27.40 |
| Pirate | AD-CG | 1.000 | 110 | 1413 | 14.42 | 27.40 |
| Moon | AD-Chol | 1.618 | 116 | - | 43.78 | 28.58 |
| Moon | AD-CG | 1.618 | 116 | 2764 | 19.75 | 28.58 |
| Moon | iAD-CG.4 | 1.468 | 120 | 1069 | 12.79 | 28.58 |
| Moon | iAD-CG1 | 1.000 | 175 | 1231 | 17.30 | 28.58 |
| Moon | iAD-CG.4 | 1.000 | 175 | 1419 | 19.10 | 28.58 |
| Moon | AD-CG | 1.000 | 175 | 3836 | 24.00 | 28.58 |

结论：

- `iAD-CG.4` 的 PSNR 与精确 ADMM 基本相同；
- `iAD-CG.4` 的累计内层迭代数显著少于 `AD-CG`；
- `iAD-CG.4` 通常比单位步长的 `iAD-CG1` 更快；
- 精确 ADMM 的 Cholesky 版本在这些测试中明显慢于 CG 版本。

### 10.3 restored image 展示实验

文献还展示了 `Pirate` 在较高噪声下的恢复图：

| 项目 | 值 |
|---|---:|
| noise level `delta` | `0.03` |
| `mu_paper` | `0.125` |
| 对应 `mu_code` | `8` |
| beta | `30` |
| noisy PSNR | `15.62 dB` |
| restored PSNR | `26.22 dB` |

图中报告时间：

| Algorithm | PSNR | Time(s) |
|---|---:|---:|
| AD-Chol | 26.22 | 68.20 |
| AD-CG | 26.22 | 24.13 |
| iAD-CG.4 | 26.22 | 16.60 |
| iAD-CG1 | 26.22 | 20.16 |

这组图的目的不是证明图像视觉上差异很大，而是说明同等 PSNR 下，`iAD-CG.4` 更省时间。

## 11. 对你的实现最重要的复现要点

### 11.1 IADMM.tv 应保留非精确 x 步

若目标是复现 Xie 文献，`IADMM.tv` 应采用：

```matlab
xCandidate = xOld;          % warm start
for inner = 1:maxInner
    % CG step for (mu_code*I + beta*D'*D)x = mu_code*c + D'*(lambda + beta*yOld)
    % update xCandidate
    d = mu_code*(xCandidate - c) ...
        - D' * lambda ...
        + beta * D' * (D*xCandidate - yOld);
    left = (2/beta) * abs((w - xCandidate)' * d) + norm(d)^2;
    right = sigma * alpha * norm(D*xCandidate - yOld)^2;
    if left <= right
        break;
    end
end
```

接受当前 `xCandidate` 后再更新 y、lambda、w。

### 11.2 精确解只能作为对照组

可以保留一个 `AD-CG` 或 `AD-FFT` 作为精确 ADMM baseline，但不应把它混作 IADMM 主算法。建议命名区分：

| 算法名 | x 子问题 |
|---|---|
| `AD-FFT` 或 `AD-Chol` | 精确或高精度 |
| `AD-CG` | 高精度 CG |
| `iAD-CG.4` | 相对误差准则非精确 CG，`alpha=0.4,tau=1.468` |
| `iAD-CG1` | 单位步长非精确 CG，`tau=1` |

### 11.3 history 应记录这些字段

为了对齐文献和方便调试，建议至少记录：

| 字段 | 含义 |
|---|---|
| `inner` | 每轮 CG 内迭代数 |
| `innerTotal` | 累计 CG 内迭代数 |
| `relErrLeft` | 相对误差准则左端 |
| `relErrRight` | 相对误差准则右端 |
| `relErrOk` | 是否由相对误差准则停止 |
| `xLinearResidual` | x 子问题最终归一化线性残差 |
| `primalResidualRel` | 外层归一化 primal residual |
| `dualResidualRel` | 外层归一化 dual residual |
| `psnr` | 文献 PSNR |
| `snr` | 你当前使用的 SNR，可保留 |
| `objPrimal` | \(f(x)+TV(Dx)\) |
| `objSplit` | \(f(x)+TV(y)\) |

### 11.4 参数建议

按文献逻辑复现时，可先设置：

| 参数 | 建议 |
|---|---|
| `sigma` | `0.99` |
| `alpha` | `0.4` |
| `tau` | `1.468` |
| `beta` | 先扫若干正值 |
| `mu_code` | 用 `1/mu_paper` 换算，文献 `0.05` 对应 `20` |
| `x0` | noisy image `c` |
| `y0, lambda0, w0` | zero |

若当前图像是 `cameraman`，不能直接照搬 Pirate/Moon 的最优 `beta`，但可以照搬实验思路：固定噪声水平后，对 `mu_code` 和 `beta` 做网格扫描，以 PSNR 和 CPU 时间共同选取。

## 12. 推荐复现实验结构

建议把 TV 图像实验拆成四组：

1. **正确性基线**：`AD-FFT` 或 `AD-CG-highacc`，确认目标解和 PSNR 水平。
2. **文献主算法**：`iAD-CG.4`，使用相对误差准则和 `tau=1.468`。
3. **单位步长对照**：`iAD-CG1`，使用非精确 CG 但 `tau=1`。
4. **步长隔离对照**：`iAD-CG.4-tau1`，保留 `alpha=0.4` 的误差准则但强制 `tau=1`。

这样可以分别回答：

- 非精确 x 步是否减少内层成本？
- 较大乘子步长是否减少外层迭代？
- `alpha=0.4` 的优势来自误差准则，还是来自 `tau=1.468`？
- 同等 PSNR 下 CPU 时间是否更低？

## 13. 与当前 TV 实验的关系

你当前实验里若把 GAIP/AGAIP 的 x 子问题改成 FFT 精确解，它们就不再体现“非精确 x 子问题”的思想。若你的新算法逻辑要求 x 子问题非精确，则应像 IADMM 文献一样：

1. 保留内层迭代；
2. 用相对误差准则，而不是单纯固定 CG 容差；
3. 记录内层是否真的满足准则；
4. 以累计内层迭代数和 CPU 时间作为主要效率指标；
5. 用精确 ADMM 只作 baseline，不作你的算法主体。

从文献角度看，实验的核心叙事不是“谁的最终 PSNR 更高”，而是：

> 在达到几乎相同恢复质量的前提下，非精确 ADMM 通过相对误差准则和较大乘子步长显著减少计算时间。

这也是你后续整理数值实验时最值得保留的主线。
