# TV 图像恢复程序整理说明

整理日期：2026-06-18

本次整理按照你的原始算法框架处理：每个算法的子问题都保留在算法函数内部求解，不再抽象成公共二次子问题求解器。原因是 `IADMM`、`IGADMM`、`GAIPADMM`、`AGAIPADMM` 的 x 子问题虽然都表现为 TV 二次系统，但右端项、近端项、辅助变量和相对误差准则并不相同。

## 1. 主程序入口

图像恢复数值实验主程序只保留：

```text
Main4.m
```

已删除额外入口：

```text
Main4_IADMM_TV.m
Main4_IADMM_TV_AlphaScan.m
```

`Main4.m` 中只保留非精确 x-CG 控制参数：

```matlab
options.tvCgMaxIter = 100;
options.tvCgAbsTol = 1e-6;
```

不再提供 `xStepMode`、`highAccuracyCg`、`fftExact` 等精确/高精度对照入口。

## 2. 子问题求解框架

四个 TV 算法均采用同一原则：

1. 在算法主循环内部构造本算法自己的 x 子问题；
2. 在算法文件内部用局部函数 `solveXByInexactCg` 求解；
3. x 子问题只做非精确 CG；
4. CG 停机由本算法自己的 `relativeCriterion` 决定；
5. y 子问题继续用 `Functions.groupSoftThreshold` 精确求解。

已删除公共求解器：

```text
+Functions/solveTvQuadraticSubproblem.m
+Functions/solveTvNormalEquation.m
+Tools/plotIadmmTvDiagnostics.m
```

## 3. 各算法的 x 子问题

### 3.1 `Algorithms.IADMM.tv`

x 子问题：

$$
(\mu I+\beta D^TD)x=\mu c+D^T(\lambda+\beta y).
$$

内部非精确准则：

$$
\frac{2}{\beta}|\langle w-x,d\rangle|+\|d\|^2
\le
\sigma\alpha\|Dx-y^k\|^2.
$$

辅助变量更新：

```matlab
auxiliaryX = auxiliaryX - beta*d;
```

### 3.2 `Algorithms.IGADMM.tv`

x 子问题仍为：

$$
(\mu I+\beta D^TD)x=\mu c+D^T(\lambda+\beta y).
$$

但相对误差准则改为 IGADMM 自己的形式：

$$
2|\langle w-x,d\rangle|+\|d\|^2
\le
(1-\alpha)\sigma_1\beta\|Dx-y\|^2.
$$

并保留中间乘子与 y 近端项：

```matlab
lambdaHalf = lambda - alpha*beta*(DxNew-y);
yCenter = (beta*DxNew - lambdaHalf + yProx*y)/(beta + yProx);
yNew = Functions.groupSoftThreshold(yCenter, 1/(beta + yProx));
lambdaNew = lambdaHalf - beta*(DxNew-yNew);
```

### 3.3 `Algorithms.GAIPADMM.tv`

先构造广义加速点：

```matlab
bar_xk = (1-theta_x)*xk + theta_x*bar_xk_prev;
bar_yk = (1-theta_y)*yk + theta_y*bar_yk_prev;
bar_lambdak = (1-theta_lambda)*lambdak + theta_lambda*bar_lambdak_prev;
```

x 子问题：

$$
\left((\mu+\tau_1/\beta)I+\beta D^TD\right)x
=
\mu c+D^T(\bar\lambda+\beta y)
+(\tau_1/\beta)\bar x.
$$

相对误差准则：

$$
2|\langle w-x,d\rangle|+\alpha\beta\|d\|^2
\le
\varepsilon_1(1+\theta_x)\frac{\tau_1}{\beta}\|x-\bar x\|^2.
$$

### 3.4 `Algorithms.AGAIPADMM.tv`

结构与 `GAIPADMM.tv` 相同，但每轮使用递增罚参数：

```matlab
betaCurrent = betaPrevious + betaIncrement;
```

x 子问题为：

$$
\left((\mu+\tau_1/\beta_k)I+\beta_kD^TD\right)x
=
\mu c+D^T(\bar\lambda+\beta_k y)
+(\tau_1/\beta_k)\bar x.
$$

相对误差准则使用当前 `betaCurrent`。

## 4. 诊断字段

四个 TV 算法保留以下 x 子问题诊断：

| 字段 | 含义 |
|---|---|
| `history.inner` | 每轮 x-CG 内迭代数 |
| `history.innerTotal` | 累计 x-CG 内迭代数 |
| `history.relErrXLeft` | x 相对误差准则左端 |
| `history.relErrXRight` | x 相对误差准则右端 |
| `history.relErrXOk` | x 相对误差准则是否满足 |
| `history.xAbsStop` | 是否由绝对线性残差阈值停机 |
| `history.xLinearResidual` | x 子问题最终归一化线性残差 |

注意：`xAbsStop` 只是数值保护，不是理论相对误差准则。若要分析非精确求解是否符合算法逻辑，应优先看 `relErrXOk` 和 `relErrXLeft ./ relErrXRight`。

## 5. 保留与删除的内容

保留：

```text
Main4.m
+Algorithms/+IADMM/tv.m
+Algorithms/+IGADMM/tv.m
+Algorithms/+GAIPADMM/tv.m
+Algorithms/+AGAIPADMM/tv.m
+Functions/tvGroupNormSum.m
```

删除：

```text
Main4_IADMM_TV.m
Main4_IADMM_TV_AlphaScan.m
+Functions/solveTvQuadraticSubproblem.m
+Functions/solveTvNormalEquation.m
+Tools/plotIadmmTvDiagnostics.m
```

`+Functions/tvLaplacianEigs.m` 仍保留，因为其他算法文件如 `APADMC.tv` 还在使用频域 TV 线性系统信息。

## 6. 后续检查建议

运行 `Main4.m` 后重点查看：

1. 每个算法的 `x-block relative criterion` 是否大多数迭代满足；
2. `max left/right` 是否长期大于 1；
3. `xAbsStop` 是否过多触发，若过多说明 `options.tvCgAbsTol` 可能干扰了理论准则；
4. `objPrimal`、`SNR`、`feasibility` 是否趋势一致。

当前环境未检测到 MATLAB 或 Octave，因此本次只完成静态结构整理，未执行数值运行。
