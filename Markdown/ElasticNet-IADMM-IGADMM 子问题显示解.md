# Elastic-net IADMM/IGADMM/GAIPADMM/AGAIPADMM 子问题显示解

## 1. 问题构造检查

截图中的 elastic-net 原问题写为

$$
\min_z\ \eta\|z\|_1+\gamma\|z\|_2^2+\frac12\|Kz-b\|_2^2,\quad \eta,\gamma>0.
$$

当前代码采用与 LASSO 一致的分裂方式：

$$
\min_{x,y}\ F(x,y)
:=\gamma\|x\|_2^2+\frac12\|Kx-b\|_2^2+\eta\|y\|_1,
\quad \text{s.t. }x-y=0.
$$

这个构造是正确的：可行约束给出 `x=y=z`，代回即得到原始 elastic-net 目标。该分裂把光滑强凸二次项放在 `x` 块，把非光滑 `l1` 项放在 `y` 块；因此 `x` 块是二次系统，`y` 块是软阈值近端。

检查结论：

| 项目 | 检查结果 |
|---|---|
| 数据矩阵符号 | 数据矩阵统一为 `problem.K`，不再保留 `problem.A` 兼容别名。 |
| 数据生成 | `Problems.elasticProblem` 沿用 LASSO 的相关高斯矩阵和噪声框架生成 `K,b`。 |
| 目标函数记录 | `Functions.objectiveValue` 的 elastic 分支已改为 `gamma*||x||^2 + 0.5*||K*x-b||^2 + eta*||y||_1`。 |
| 约束矩阵范数 | elastic-net 拆分约束是 `x-y=0`，约束算子范数平方为 `1`；`constraintNorm2=1`，`normH2=1`。 |

## 2. 公共记号

负号拉格朗日约定为

$$
L_\beta(x,y;\lambda)
=\gamma\|x\|^2+\frac12\|Kx-b\|^2+\eta\|y\|_1
-\langle\lambda,x-y\rangle+\frac{\beta}{2}\|x-y\|^2.
$$

软阈值算子记为

$$
\operatorname{shrink}(v,\kappa)
=\operatorname{sign}(v)\odot\max(|v|-\kappa,0).
$$

对任意 `a>0` 和右端 `r`，二次系统

$$
(K^TK+aI)x=r
$$

的显示解为

$$
x=(K^TK+aI)^{-1}r.
$$

代码在 `m<=n` 时使用等价 Woodbury 形式

$$
x=\frac1a\left[r-K^T(KK^T+aI)^{-1}Kr\right].
$$

## 3. IADMM.elastic

对应文件：`+Algorithms/+IADMM/elastic.m`。

### 3.1 x 子问题

固定 `y^k,lambda^k`：

$$
x^{k+1}
=\arg\min_x\ \gamma\|x\|^2+\frac12\|Kx-b\|^2
-\langle\lambda^k,x-y^k\rangle
+\frac{\beta}{2}\|x-y^k\|^2.
$$

一阶条件为

$$
(K^TK+(2\gamma+\beta)I)x^{k+1}
=K^Tb+\beta y^k+\lambda^k.
$$

### 3.2 y 子问题

固定 `x^{k+1},lambda^k`：

$$
y^{k+1}
=\arg\min_y\ \eta\|y\|_1
+\langle\lambda^k,y\rangle
+\frac{\beta}{2}\|x^{k+1}-y\|^2.
$$

因此

$$
y^{k+1}
=\operatorname{shrink}\left(x^{k+1}-\lambda^k/\beta,\eta/\beta\right).
$$

### 3.3 乘子更新

$$
\lambda^{k+1}
=\lambda^k-\tau\beta(x^{k+1}-y^{k+1}).
$$

## 4. IGADMM.elastic

对应文件：`+Algorithms/+IGADMM/elastic.m`。文献建议参数在代码中直接固定，设 `t=0.01`。

### 4.1 x 子问题

$$
(K^TK+(2\gamma+\beta)I)x^{k+1}
=K^Tb+\beta y^k+\lambda^k.
$$

### 4.2 半步乘子

$$
\lambda^{k+\frac12}
=\lambda^k-\alpha\beta(x^{k+1}-y^k).
$$

### 4.3 y 子问题

$$
y^{k+1}
=\operatorname{shrink}\left(
\frac{\beta x^{k+1}-\lambda^{k+\frac12}+t y^k}{\beta+t},
\frac{\eta}{\beta+t}
\right).
$$

### 4.4 乘子更新

$$
\lambda^{k+1}
=\lambda^{k+\frac12}-\beta(x^{k+1}-y^{k+1}).
$$

## 5. GAIPADMM.elastic

对应文件：`+Algorithms/+GAIPADMM/elastic.m`。

参数沿用 `GAIPADMM.lasso`：

| 参数 | 含义 |
|---|---|
| `alpha,beta` | 乘子步长与增广拉格朗日罚参数 |
| `theta_x,theta_y,theta_lambda` | 广义加速参数 |
| `tau1,tau2` | x/y 块近端参数 |
| `eps1,eps3` | x 块相对误差和理论参数条件 |

广义加速点：

$$
\bar x^k=(1-\theta_x)x^k+\theta_x\bar x^{k-1},\quad
\bar y^k=(1-\theta_y)y^k+\theta_y\bar y^{k-1},
$$

$$
\bar\lambda^k=(1-\theta_\lambda)\lambda^k+\theta_\lambda\bar\lambda^{k-1}.
$$

### 5.1 x 子问题

带近端项的 x 块一阶条件为

$$
\left(K^TK+\left(2\gamma+\beta+\tau_1/\beta\right)I\right)x^{k+1}
=K^Tb+\beta y^k+\bar\lambda^k+(\tau_1/\beta)\bar x^k.
$$

代码用 CG 求解该 SPD 系统，并记录相对误差条件。

### 5.2 y 子问题

$$
y^{k+1}
=\operatorname{shrink}\left(
\frac{\beta x^{k+1}-\bar\lambda^k+(\tau_2/\beta)\bar y^k}
{\beta+\tau_2/\beta},
\frac{\eta}{\beta+\tau_2/\beta}
\right).
$$

### 5.3 乘子更新

$$
\lambda^{k+1}
=\bar\lambda^k-\alpha\beta(x^{k+1}-y^{k+1}).
$$

## 6. AGAIPADMM.elastic

对应文件：`+Algorithms/+AGAIPADMM/elastic.m`。

elastic-net 的 `x` 块

$$
\gamma\|x\|^2+\frac12\|Kx-b\|^2
$$

至少具有强凸模量 `2*gamma`。GAIP-ADMMspp 要求强凸模量最大的块作为第二个子问题，因此代码在算法顺序上采用 `(y,x)`：先更新非光滑 `y` 块，再更新强凸 `x` 块；输出仍保持 `out.x` 为二次系统变量、`out.y` 为软阈值变量。

参数：

| 参数 | 含义 |
|---|---|
| `alpha,beta` | 乘子步长与初始罚参数 |
| `theta_x,theta_y,theta_lambda` | 广义加速参数 |
| `theta_gamma` | 罚参数序列中强凸模量的凸组合权重 |
| `tau1,tau2` | x/y 块近端参数 |
| `eps1,eps2,eps3` | y/x 块相对误差与可行性残差误差参数 |

罚参数序列：

$$
\beta_{-1}=\beta,\quad
\beta_k=\beta_{k-1}+\frac{\theta_\gamma(2\gamma)}{3(1+\|I\|^2)}
=\beta_{k-1}+\frac{\theta_\gamma\gamma}{3}.
$$

广义加速点：

$$
\bar x^k=(1-\theta_x)x^k+\theta_x\bar x^{k-1},\quad
\bar y^k=(1-\theta_y)y^k+\theta_y\bar y^{k-1},
$$

$$
\bar\lambda^k=(1-\theta_\lambda)\lambda^k+\theta_\lambda\bar\lambda^{k-1}.
$$

### 6.1 y 子问题

第一块为 `y`，固定 `x^k,\bar\lambda^k` 后

$$
y^{k+1}
=\operatorname{shrink}\left(
\frac{\beta_kx^k-\bar\lambda^k+(\tau_2/\beta_k)\bar y^k}
{\beta_k+\tau_2/\beta_k},
\frac{\eta}{\beta_k+\tau_2/\beta_k}
\right).
$$

该子问题精确求解，因此 `d_y=0`。

### 6.2 x 子问题

第二块为强凸 `x`，固定 `y^{k+1},\bar\lambda^k` 后一阶条件为

$$
\left(K^TK+\left(2\gamma+\beta_k+\tau_1/\beta_k\right)I\right)x^{k+1}
=K^Tb+\beta_ky^{k+1}+\bar\lambda^k+(\tau_1/\beta_k)\bar x^k.
$$

代码用 CG 求解该 SPD 系统，并按第二块相对误差条件记录：

$$
2|\langle w_x^k-x^{k+1},d_x^{k+1}\rangle|
+\alpha\beta_k\|d_x^{k+1}\|^2
\le
\frac{\epsilon_2(1+\theta_x)}{\beta_k}
\|x^{k+1}-\bar x^k\|_{\tau_1 I}^2
+\epsilon_3\beta_k\|x^{k+1}-y^{k+1}\|^2.
$$

### 6.3 乘子更新

$$
\lambda^{k+1}
=\bar\lambda^k-\alpha\beta_k(x^{k+1}-y^{k+1}).
$$
