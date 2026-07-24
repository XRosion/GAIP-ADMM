# InexactADMM-LASSO 检查

检查日期：2026-06-08

检查对象：

- [`+Algorithms/+IADMM/lasso.m`](../+Algorithms/+IADMM/lasso.m)
- [`+Algorithms/+IGADMM/lasso.m`](../+Algorithms/+IGADMM/lasso.m)
- [`+Algorithms/+GAIPADMM/lasso.m`](../+Algorithms/+GAIPADMM/lasso.m)

参考文献：

1. Jiaxin Xie, *On inexact ADMMs with relative error criteria*, `参考文献/IADMM.pdf`.
2. Z. Wu, Y. Song, F. Jiang, *Inexact generalized ADMM with relative error criteria for linearly constrained convex optimization problems*, `参考文献/IGADMM.pdf`.
3. Jinbao Jian, Xiang Luo, Jianghua Yin, Xianke Tang, *Convergence of generalized accelerated and inexact proximal ADMM with sequence of penalty parameters for two-block pre-strongly convex optimization*, `参考文献/GAIPADMM-0601.pdf`.

## 1. 公共 LASSO 分裂口径

当前工程中的 LASSO 问题由 [`Problems.lassoProblem`](../+Problems/lassoProblem.m) 生成：

$$
\min_x \ \eta\|x\|_1 + \frac12\|Kx-b\|^2.
$$

三个非精确 ADMM 函数均应使用如下二块分裂：

$$
\min_{x,y}\ f(x)+g(y)
= \frac12\|Kx-b\|^2 + \eta\|y\|_1,
\quad \text{s.t. } x-y=0.
$$

对应到通用 ADMM 约束形式：

$$
A x + B y - c = 0,\quad A=I,\ B=-I,\ c=0.
$$

增广拉格朗日函数采用文献中一致的符号：

$$
L_\beta(x,y;\lambda)
= f(x)+g(y)-\langle\lambda,x-y\rangle
+\frac{\beta}{2}\|x-y\|^2.
$$

因此 x 子问题一阶残差为：

$$
d_x =
K^\top(Kx-b)-\lambda+\beta(x-y)
$$

若有 GAIPADMM 的 x 近端项
$\frac{1}{2\beta}\|x-\bar x\|_S^2,\ S=\tau_1 I$，则：

$$
d_x =
K^\top(Kx-b)-\bar\lambda+\beta(x-y)
+\frac{\tau_1}{\beta}(x-\bar x).
$$

y 子问题为软阈值问题，这是三个函数中将 y 块视为精确求解的基础。

## 2. 文献算法要点

### 2.1 IADMM

`IADMM.pdf` 第 5 页 Algorithm 1 给出的部分非精确 ADMM 为：

1. 近似求解
   $x^{k+1}\approx \arg\min_x L_\beta(x,y^k;\lambda^k)$，
   并满足相对误差准则

   $$
   \frac{2}{\beta}|\langle w^k-x^{k+1},d^{k+1}\rangle|
   +\|d^{k+1}\|^2
   \le
   \sigma\alpha\|Ax^{k+1}+By^k-c\|^2.
   $$

2. 精确求解
   $y^{k+1}=\arg\min_y L_\beta(x^{k+1},y;\lambda^k)$。

3. 更新

   $$
   \lambda^{k+1}
   =\lambda^k-\tau\beta(Ax^{k+1}+By^{k+1}-c),
   \quad
   w^{k+1}=w^k-\beta d^{k+1}.
   $$

参数条件为：

$$
\beta>0,\quad \sigma\in[0,1),\quad \alpha\in(0,1),
\quad
0<\tau<
\frac{2\gamma-1+\sqrt{4\gamma^2+1}}{2\gamma},
\quad \gamma=1-\alpha.
$$

`IADMM.pdf` 第 15 页说明，当 x 子问题化为大型线性系统时，可用 CG 非精确求解；y 子问题可由 soft shrinkage 精确更新。

### 2.2 IGADMM

`IGADMM.pdf` 第 5 页 Algorithm 1 给出 inexact generalized proximal ADMM：

1. 近似求解
   $L_\beta(x,y^k,\lambda^k)+\frac12\|x-x^k\|_S^2$，
   其残差为

   $$
   d_1^{k+1}\in
   \partial\theta_1(x^{k+1})
   -A^\top\lambda^k
   +\beta A^\top(Ax^{k+1}+By^k-c)
   +S(x^{k+1}-x^k).
   $$

2. 半步乘子：

   $$
   \lambda^{k+\frac12}
   =\lambda^k-\alpha\beta(Ax^{k+1}+By^k-c).
   $$

3. 近似求解
   $L_\beta(x^{k+1},y,\lambda^{k+1/2})+\frac12\|y-y^k\|_T^2$。

4. 更新

   $$
   v^{k+1}=v^k-d^{k+1},
   \quad
   \lambda^{k+1}
   =\lambda^{k+\frac12}
   -\beta(Ax^{k+1}+By^{k+1}-c).
   $$

`IGADMM.pdf` 第 6 页的 IGADMMa x 块相对误差准则为：

$$
2|\langle v_1^k-x^{k+1},d_1^{k+1}\rangle|
+\|d_1^{k+1}\|^2
\le
\sigma\|x^{k+1}-x^k\|_S^2
+(1-\alpha)\sigma_1\beta
\|Ax^{k+1}+By^k-c\|^2.
$$

文献第 17 页的 LASSO 实验采用
$S=0,\ T=\mu I,\ \sigma=\sigma_2=0$，因此 y 子问题有闭式软阈值：

$$
y^{k+1}
=\operatorname{Shrinkage}\left(
\frac{\beta x^{k+1}+\mu y^k-\lambda^{k+\frac12}}
{\beta+\mu},
\frac{\gamma}{\beta+\mu}
\right).
$$

同时文献实验取 $\beta=1,\ \mu=0.01$。

### 2.3 GAIPADMM

`GAIPADMM-0601.pdf` 第 8 页 Algorithm 1 给出 GAIP-ADMM：

1. 若强凸模量 $\gamma=0$，则惩罚参数序列不增长，即 $\beta_k=\beta$。
2. 广义加速点：

   $$
   \bar z^k=(1-\theta_z)z^k+\theta_z\bar z^{k-1},
   \quad z\in\{x,y,\lambda\}.
   $$

3. x 子问题：

   $$
   x^{k+1}\approx
   \arg\min_x
   L_{\beta_k}(x,y^k,\bar\lambda^k)
   +\frac{1}{2\beta_k}\|x-\bar x^k\|_S^2,
   $$

   并满足

   $$
   2|\langle w_1^k-x^{k+1},d_1^{k+1}\rangle|
   +\alpha\beta_k\|d_1^{k+1}\|^2
   \le
   \frac{\epsilon_1(1+\theta_x)}{\beta_k}
   \|x^{k+1}-\bar x^k\|_S^2.
   $$

4. y 子问题：

   $$
   y^{k+1}\approx
   \arg\min_y
   L_{\beta_k}(x^{k+1},y,\bar\lambda^k)
   +\frac{1}{2\beta_k}\|y-\bar y^k\|_T^2.
   $$

5. 更新：

   $$
   \lambda^{k+1}
   =\bar\lambda^k-\alpha\beta_k\psi(x^{k+1},y^{k+1}),
   \quad
   w_i^{k+1}=w_i^k-\alpha\beta_k d_i^{k+1}.
   $$

参数条件包括：

$$
\theta_x,\theta_y\in(-1,1),\quad
\theta_\lambda<1,\quad
\epsilon_1,\epsilon_2\in[0,1),\quad
\epsilon_3\ge0,
$$

以及

$$
1-\alpha(1-\theta_\lambda)-\epsilon_3>0.
$$

## 3. 函数逐项检查

### 3.1 `IADMM.lasso`

检查结论：**正确**，与 Xie Algorithm 1 应用于当前 LASSO 分裂的公式一致。

核对点：

- 数据字段使用 `problem.K`、`problem.b`、`problem.eta`，与 [`Problems.lassoProblem`](../+Problems/lassoProblem.m) 一致。
- x 线性系统右端为
  `K'*b + beta*y + lambda`，对应
  $(K^\top K+\beta I)x=K^\top b+\beta y+\lambda$。
- CG 使用矩阵自由乘法
  `K'*(K*v)+beta*v`，方向为 `r=-d`，符号正确。
- 相对误差准则
  `left=(2/beta)*abs((auxiliaryX-x)'*d)+d'*d`，
  `right=sigma*alpha*norm(x-previousY)^2`
  正好对应
  $A=I,B=-I,c=0$ 的文献准则。
- y 更新
  `softThreshold(x - lambda/beta, eta/beta)`，
  对应精确求解 $L_\beta(x^{k+1},y;\lambda^k)$。
- 乘子更新
  `lambda = lambda - tau*beta*(x-y)`，
  与文献符号一致。
- 辅助变量更新
  `auxiliaryX = auxiliaryX - beta*innerInfo.d`，
  与文献 $w^{k+1}=w^k-\beta d^{k+1}$ 一致。
- 参数检查覆盖 `beta>0`、`alpha in (0,1)`、`sigma in [0,1)`、`tau` 上界。

注意事项：

- 代码允许 CG 绝对残差 `lassoCgAbsTol` 触发停止，并将 `relOk` 记为真。这是实用数值策略；严格理论准则只要求相对误差条件，但精确解附近 $d\approx0$ 时可以视为近似满足。

### 3.2 `IGADMM.lasso`

检查结论：**默认参数 `igadmmBeta=1` 时与文献 LASSO 实验口径一致；若允许任意 `igadmmBeta`，当前相对误差准则与辅助变量缩放不完全符合 `IGADMM.pdf` Algorithm 1。**

正确点：

- 数据字段已统一为 `problem.K`，与当前 LASSO 生成器一致。
- x 子问题采用 $S=0$ 的 LASSO 线性系统：

  $$
  (K^\top K+\beta I)x=K^\top b+\beta y+\lambda.
  $$

- CG 残差方向 `r=-d` 正确。
- 半步乘子
  `lambdaHalf = lambda - alpha*beta*(xNew-y)`，
  对应文献
  $\lambda^{k+1/2}=\lambda^k-\alpha\beta(x^{k+1}-y^k)$。
- y 更新
  `softThreshold((beta*xNew - lambdaHalf + mu*y)/(beta+mu), eta/(beta+mu))`
  与文献第 17 页 LASSO 公式一致，其中 `mu = igadmmYProximal`。
- 最终乘子
  `lambdaNew = lambdaHalf - beta*(xNew-yNew)`，
  与文献 Algorithm 1 一致。
- y 子问题精确求解，因此记录 `relErrYLeft=0` 是合理的。

主要风险：

1. **相对误差准则和辅助变量更新只在 `beta=1` 时严格匹配文献。**

   文献 IGADMMa 在 $S=0,\sigma=0$ 下的 x 准则应为：

   $$
   2|\langle v_1^k-x,d\rangle|+\|d\|^2
   \le
   (1-\alpha)\sigma_1\beta\|x-y^k\|^2,
   $$

   且辅助变量为：

   $$
   v_1^{k+1}=v_1^k-d.
   $$

   当前代码使用：

   ```matlab
   leftValue = (2 / beta) * abs((auxiliaryX - x)' * d) + d' * d;
   rightValue = max(1 - alpha, 0) * sigma1 / beta * norm(x - y)^2;
   auxiliaryX = auxiliaryX - beta * innerInfo.d;
   ```

   当 `beta=1` 时，上式退化为文献形式；当 `beta ~= 1` 时，左端、右端和辅助变量缩放均偏离文献。

   若希望支持任意 `igadmmBeta`，建议改为：

   ```matlab
   leftValue = 2 * abs((auxiliaryX - x)' * d) + d' * d;
   rightValue = max(1 - alpha, 0) * sigma1 * beta * norm(x - y)^2;
   auxiliaryX = auxiliaryX - innerInfo.d;
   ```

2. **缺少 `alpha in (-1,1)` 参数检查。**

   `IGADMM.pdf` Algorithm 1 初始化要求 $\alpha\in(-1,1)$。当前代码默认 `alpha=-0.3` 合法，但没有阻止用户传入区间外的 `igadmmAlpha`。

3. **`igadmmYProximal=0` 与文献 `T` 正定假设不完全一致。**

   文献 Algorithm 1 要求 $T$ 正定，LASSO 实验取 $\mu=0.01$。当前代码允许 `igadmmYProximal=0`。这在数值上仍可能可运行，但不属于文献给出的正定近端设置。

### 3.3 `GAIPADMM.lasso`

检查结论：**对当前欠定 LASSO 分裂是正确的**。当前实现与 `GAIPADMM-0601.pdf` Algorithm 1 在一般凸情形 $\gamma=0$ 下的公式一致。

核对点：

- 当前 LASSO 生成器使用 $m<n$，$K^\top K$ 通常不满秩，$f(x)=\frac12\|Kx-b\|^2$ 非强凸，$g(y)=\eta\|y\|_1$ 也非强凸。因此可按 $\gamma_f=\gamma_g=0$ 的一般凸情形处理，惩罚参数 $\beta_k$ 不增长，代码中的常数 `beta` 合理。
- 加速点
  `bar_xk`、`bar_yk`、`bar_lambdak`
  对应文献公式
  $\bar z^k=(1-\theta_z)z^k+\theta_z\bar z^{k-1}$。
- x 子问题残差
  `dk1 = K'*(K*xCandidate) + (beta+tau1/beta)*xCandidate - rhs`
  对应

  $$
  K^\top(Kx-b)-\bar\lambda+\beta(x-y)
  +\frac{\tau_1}{\beta}(x-\bar x).
  $$

- x 子问题线性系统使用矩阵自由乘法
  `K'*(K*v)+c*v`，没有显式形成稠密 `K'*K`，实现方式合理。
- CG 初始方向使用 `xResidual=-dk1`，方向正确。
- x 相对误差准则

  ```matlab
  2*abs((wk1-xCandidate)'*dk1) + alpha*beta*(dk1'*dk1)
  <= eps1*(1+theta_x)*(tau1/beta)*norm(xCandidate-bar_xk)^2
  ```

  与文献 (21a) 在 $S=\tau_1I$ 下一致。
- y 子问题软阈值中心

  $$
  \frac{\beta x^{k+1}-\bar\lambda^k+(\tau_2/\beta)\bar y^k}
  {\beta+\tau_2/\beta}
  $$

  与 $T=\tau_2I$ 的精确 proximal 解一致。
- 乘子更新
  `lambda = bar_lambdak - alpha*beta*(xk-yk)`
  与文献 (23) 一致。
- 辅助变量更新
  `wk1 = wk1 - alpha*beta*dk1`
  与文献 (24) 的 x 块一致；y 块精确求解时 $d_2=0$，无需显式维护 `wk2`。
- 参数检查覆盖
  `alpha>0`、`beta>0`、`tau1>0`、`tau2>0`、`theta_x,theta_y in (-1,1)`、`theta_lambda<1`、`eps1 in [0,1)`、`eps3>=0` 和
  `1-alpha*(1-theta_lambda)-eps3>0`。

注意事项：

- 如果未来将 GAIPADMM 用于强凸 LASSO 变体，例如 $K$ 满列秩或额外加入 $\ell_2$ 正则项，则文献中的 sequence of penalty parameters 机制可能需要启用；当前函数只实现了 $\gamma=0$ 时的常数 `beta` 情形。
- 与 IADMM 一样，代码允许 CG 绝对残差触发停止。这是实用策略，严格理论检查仍以相对误差准则为主。

## 4. 公共目标函数和可行性记录

公共函数与当前分裂基本一致：

- [`objectiveValue`](../+Functions/objectiveValue.m) 中 `case 'lasso'` 实际计算

  $$
  \frac12\|Kx-b\|^2+\eta\|y\|_1.
  $$

- [`feasibilityNorm`](../+Functions/feasibilityNorm.m) 对 LASSO 根据 `y` 维度判断：
  若 `numel(y)==size(K,2)`，记录 $\|x-y\|$，与当前三个函数一致。

非阻塞注记：

- 旧注释中的 `A*x` 符号已统一为 `K*x`；算法和问题构造均使用 `problem.K`。

## 5. 总体结论

| 函数 | 结论 | 主要依据 |
|---|---|---|
| `IADMM.lasso` | 正确 | 完整匹配 Xie Algorithm 1 在 `x-y=0` LASSO 分裂下的 x 相对误差、y 软阈值、乘子和辅助变量更新 |
| `IGADMM.lasso` | 正确 | 文献建议参数已在函数内部固定，`beta=1`、`alpha=-0.3`、`sigma1=0.99`、`y` 近端参数 `0.01` |
| `GAIPADMM.lasso` | 正确 | 匹配 GAIP-ADMM Algorithm 1 的一般凸常数 `beta` 情形，x/y 子问题、相对误差准则和乘子更新均一致 |

建议优先修正或补强：

1. 若后续重新开放 IGADMM 可调参数，应同步重查 x 相对误差准则和 `auxiliaryX` 更新缩放。
2. 当前 IADMM/IGADMM 的 CG 内迭代上限统一读取 `options.lassoCgMaxIter`，主脚本和 LASSO 敏感性脚本均设为 `200`。

## 6. 验证限制

本次检查基于 PDF 文献公式、代码静态阅读和公式推导完成。当前执行环境未提供 MATLAB 或 Octave 可执行程序，因此未运行 `Main2.m` 做数值回归验证。
