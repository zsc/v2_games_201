# 第五章：欧拉视角（2）：线性系统求解器

在欧拉视角的物理仿真中，我们经常需要求解大规模稀疏线性系统，特别是在压力投影步骤中出现的泊松方程。本章将深入探讨这些线性系统的特性以及高效求解它们的各种方法。从基础的迭代法到先进的多重网格方法，我们将系统地学习如何在保证精度的同时实现高性能计算。

## 5.1 稀疏矩阵与零空间(Nullspaces)

### 5.1.1 稀疏矩阵存储格式

在物理仿真中，离散化后的线性系统通常具有稀疏性——矩阵中大部分元素为零。例如，二维5点Laplace算子每行最多只有5个非零元素。高效存储这些稀疏矩阵对于节省内存和加速计算至关重要。

**压缩稀疏行格式(CSR)**是最常用的存储格式：
- `values[]`: 存储所有非零元素
- `col_indices[]`: 每个非零元素的列索引
- `row_ptrs[]`: 每行的起始位置在values数组中的索引

例如，对于矩阵：
$$A = \begin{bmatrix} 4 & -1 & 0 \\ -1 & 4 & -1 \\ 0 & -1 & 4 \end{bmatrix}$$

CSR表示为：
- `values = [4, -1, -1, 4, -1, -1, 4]`
- `col_indices = [0, 1, 0, 1, 2, 1, 2]`
- `row_ptrs = [0, 2, 5, 7]`

**坐标格式(COO)**适合构建阶段：
- `(row[], col[], value[])` 三元组表示每个非零元素

### 5.1.2 兼容性条件

对于奇异系统 $Ax = b$，解存在的必要条件是右端项 $b$ 必须垂直于矩阵 $A$ 的零空间。

考虑纯Neumann边界条件下的泊松方程：
$$\nabla^2 p = \nabla \cdot u^*$$

离散化后得到的线性系统具有一维零空间，其基向量为 $\mathbf{1} = [1, 1, ..., 1]^T$。这是因为压力场可以相差一个常数而不影响压力梯度。

兼容性条件要求：
$$\mathbf{1}^T b = \sum_i b_i = \sum_i (\nabla \cdot u^*)_i = 0$$

这在物理上对应于不可压缩条件——流入等于流出。

### 5.1.3 零空间投影

当系统具有非平凡零空间时，解不唯一。我们需要将解投影到零空间的正交补空间中。

对于压力泊松方程，常用的处理方法包括：

1. **固定一点压力**：设置 $p_0 = 0$，移除一个自由度
2. **投影法**：求解后减去平均值 $p = p - \frac{1}{n}\sum_i p_i$
3. **增广系统**：添加约束 $\sum_i p_i = 0$

投影操作可以表示为：
$$P = I - \frac{\mathbf{1}\mathbf{1}^T}{n}$$

其中 $P$ 是投影矩阵，将向量投影到零空间的正交补。

### 5.1.4 奇异系统的处理

处理奇异系统的实用策略：

1. **正则化**：添加小量 $\epsilon I$ 使矩阵非奇异
   $$(A + \epsilon I)x = b$$
   
2. **最小二乘解**：求解 $\min ||Ax - b||_2$，使用广义逆
   $$x = A^+ b$$
   
3. **兼容性修正**：强制右端项满足兼容性条件
   $$b' = b - \frac{\mathbf{1}^T b}{n}\mathbf{1}$$

4. **迭代求解器调整**：
   - 在CG中使用投影预条件
   - 在多重网格中特殊处理粗网格零空间

## 5.2 Krylov子空间求解器

### 5.2.1 Krylov子空间理论

Krylov子空间方法是求解大规模稀疏线性系统的主力军。对于系统 $Ax = b$，Krylov子空间定义为：

$$\mathcal{K}_m(A, r_0) = \text{span}\{r_0, Ar_0, A^2r_0, ..., A^{m-1}r_0\}$$

其中 $r_0 = b - Ax_0$ 是初始残差。

核心思想是在逐渐扩大的Krylov子空间中寻找最优近似解：
$$x_m \in x_0 + \mathcal{K}_m(A, r_0)$$

这类方法的优势：
- 只需要矩阵-向量乘积
- 内存需求低（通常是 $O(n)$）
- 可以利用矩阵的特殊结构

### 5.2.2 共轭梯度法(CG)

对于对称正定(SPD)矩阵，共轭梯度法是最优选择。算法通过构造一组 $A$-共轭的搜索方向 $\{p_k\}$：

$$p_i^T A p_j = 0, \quad i \neq j$$

**CG算法核心步骤**：
```
初始化: r_0 = b - Ax_0, p_0 = r_0
for k = 0, 1, 2, ... do
    α_k = (r_k^T r_k) / (p_k^T A p_k)     # 步长
    x_{k+1} = x_k + α_k p_k              # 更新解
    r_{k+1} = r_k - α_k A p_k            # 更新残差
    β_k = (r_{k+1}^T r_{k+1}) / (r_k^T r_k)  # 共轭系数
    p_{k+1} = r_{k+1} + β_k p_k          # 新搜索方向
end
```

理论上，CG在 $n$ 步内收敛到精确解（不考虑舍入误差）。实际收敛速度取决于条件数：

$$||e_k||_A \leq 2\left(\frac{\sqrt{\kappa} - 1}{\sqrt{\kappa} + 1}\right)^k ||e_0||_A$$

其中 $\kappa = \lambda_{\max}/\lambda_{\min}$ 是条件数。

### 5.2.3 BiCGSTAB方法

对于非对称矩阵，BiCGSTAB（双共轭梯度稳定化）是常用选择。它结合了BiCG的思想和稳定化技术：

```
初始化: r_0 = b - Ax_0, 选择 r̃_0 (通常 r̃_0 = r_0)
for k = 0, 1, 2, ... do
    ρ_k = r̃_0^T r_k
    β = (ρ_k/ρ_{k-1}) × (α/ω_{k-1})
    p_k = r_k + β(p_{k-1} - ω_{k-1}v_{k-1})
    v_k = Ap_k
    α = ρ_k / (r̃_0^T v_k)
    s = r_k - αv_k
    t = As
    ω_k = (t^T s) / (t^T t)
    x_{k+1} = x_k + αp_k + ω_k s
    r_{k+1} = s - ω_k t
end
```

BiCGSTAB的优点是避免了BiCG中的不规则收敛行为，缺点是每步需要两次矩阵-向量乘积。

### 5.2.4 收敛性分析

Krylov方法的收敛性主要受以下因素影响：

1. **谱分布**：特征值聚集程度比条件数更重要
2. **右端项**：在特征向量基下的分解影响收敛
3. **舍入误差**：可能导致正交性损失

**重启策略**：
对于GMRES等方法，存储需求随迭代次数增长。GMRES(m)每 $m$ 步重启：
- 优点：限制内存使用
- 缺点：可能减慢收敛甚至停滞

**实用收敛准则**：
- 相对残差：$||r_k||/||b|| < \epsilon$
- 相对改变：$||x_{k+1} - x_k||/||x_k|| < \epsilon$
- 组合准则：同时考虑残差和解的变化

## 5.3 预条件(Preconditioning)

### 5.3.1 预条件的作用

预条件是加速Krylov方法收敛的关键技术。基本思想是找到一个容易求逆的矩阵 $M \approx A$，将原系统转换为条件数更好的等价系统。

**左预条件**：
$$M^{-1}Ax = M^{-1}b$$

**右预条件**：
$$AM^{-1}y = b, \quad x = M^{-1}y$$

**分裂预条件**（对称情况）：
$$L^{-1}AL^{-T}\hat{x} = L^{-1}b, \quad x = L^{-T}\hat{x}$$

其中 $M = LL^T$。

理想的预条件器应该满足：
1. $M^{-1}A$ 的条件数远小于 $A$ 的条件数
2. $M^{-1}v$ 容易计算
3. $M$ 的构造和存储开销合理

预条件的效果可以通过谱分析理解。如果 $A$ 的特征值分布在 $[\lambda_{\min}, \lambda_{\max}]$，好的预条件器会使 $M^{-1}A$ 的特征值聚集在1附近。

### 5.3.2 Jacobi预条件

最简单的预条件是对角预条件（Jacobi预条件）：
$$M = \text{diag}(A) = \text{diag}(a_{11}, a_{22}, ..., a_{nn})$$

**优点**：
- 构造简单：$O(n)$ 时间和空间
- 完全并行：$M^{-1}v$ 的计算无数据依赖
- 适合GPU实现

**缺点**：
- 改善有限：通常只能将条件数减少常数倍
- 对病态问题效果差

**加权Jacobi**：
$$M = \omega \text{diag}(A)$$

其中 $\omega \in (0, 1]$ 是松弛因子。对于某些问题，$\omega < 1$ 可以改善稳定性。

### 5.3.3 不完全Cholesky分解

对于SPD矩阵，完全Cholesky分解 $A = LL^T$ 提供了完美预条件，但分解过程中的填充(fill-in)使其对大规模稀疏矩阵不实用。

不完全Cholesky分解(IC)通过限制填充来保持稀疏性：

**IC(0)**：零填充，只在 $A$ 的非零位置计算 $L$
```
for i = 1 to n do
    L_{ii} = sqrt(A_{ii} - sum_{k<i, L_{ik}≠0} L_{ik}^2)
    for j > i where A_{ij} ≠ 0 do
        L_{ji} = (A_{ji} - sum_{k<i, L_{jk}≠0, L_{ik}≠0} L_{jk}L_{ik}) / L_{ii}
    end
end
```

**IC(p)**：允许 $p$ 级填充
- 级别定义：$\text{level}(i,j) = \min_{路径} \sum \text{level}(边)$
- 只计算 $\text{level}(i,j) \leq p$ 的元素

**阈值IC(τ)**：基于数值大小的填充策略
- 如果 $|L_{ij}| > \tau \cdot ||L_{i,:}||$，则保留该元素

### 5.3.4 修正不完全Cholesky(MIC)

标准IC可能不稳定，特别是对于接近奇异的矩阵。修正IC通过保持某些数学性质来改善稳定性。

**MIC(0)**保持行和：
$$\sum_j L_{ij}L_{jk} = \sum_j A_{jk}$$

实现方式是将丢弃的填充值累加到对角元：
```
for i = 1 to n do
    dropped_sum = 0
    for j < i do
        if (i,j) not in pattern then
            dropped_sum += computed_value^2
        end
    end
    L_{ii} = sqrt(A_{ii} + dropped_sum - sum_{k<i} L_{ik}^2)
end
```

**AINV预条件**：近似逆预条件
直接近似 $A^{-1}$ 而不是分解 $A$：
- 计算稀疏矩阵 $Z \approx L^{-1}$
- $M^{-1} = Z^TZ \approx A^{-1}$

## 5.4 多重网格方法

### 5.4.1 误差的频率分析

多重网格方法的核心洞察是：简单迭代法（如Jacobi、Gauss-Seidel）对误差的不同频率成分有不同的效果。

考虑一维Poisson方程的误差传播。对于网格间距 $h$，误差可以分解为不同频率的Fourier模式：
$$e^{(k)} = \sum_{j=1}^{n-1} \alpha_j \sin(j\pi x/L)$$

**光滑性质**：
- 高频误差（$j > n/2$）：被Jacobi/GS快速衰减
- 低频误差（$j \leq n/2$）：衰减缓慢

误差衰减因子：
$$\mu_j = 1 - \frac{4\sin^2(j\pi h/2)}{4/h^2} = \cos^2(j\pi h/2)$$

- 高频模式（$j = n/2$）：$\mu_{n/2} = 0$（一步消除）
- 低频模式（$j = 1$）：$\mu_1 \approx 1 - (\pi h)^2/2$（衰减极慢）

**多重网格思想**：
在粗网格上，原来的低频误差变成高频误差，可以被有效消除。

### 5.4.2 限制与延拓算子

网格间的信息传递通过限制(restriction)和延拓(prolongation)算子实现。

**限制算子** $I_{2h}^h: \Omega_h \to \Omega_{2h}$（细到粗）：

1. **注入(Injection)**：
   $$u_{2h,i} = u_{h,2i}$$
   
2. **全权重(Full weighting)**（一维）：
   $$u_{2h,i} = \frac{1}{4}u_{h,2i-1} + \frac{1}{2}u_{h,2i} + \frac{1}{4}u_{h,2i+1}$$
   
3. **全权重(Full weighting)**（二维）：
   $$u_{2h,i,j} = \frac{1}{16}\begin{bmatrix}1 & 2 & 1\\2 & 4 & 2\\1 & 2 & 1\end{bmatrix} \cdot u_h$$

**延拓算子** $I_h^{2h}: \Omega_{2h} \to \Omega_h$（粗到细）：

1. **线性插值**（一维）：
   $$u_{h,2i} = u_{2h,i}$$
   $$u_{h,2i+1} = \frac{1}{2}(u_{2h,i} + u_{2h,i+1})$$
   
2. **双线性插值**（二维）：
   使用双线性基函数，粗网格节点直接复制，其他点插值

**Galerkin条件**：
为保证变分性质，通常选择：
$$I_{2h}^h = c(I_h^{2h})^T$$

### 5.4.3 V-cycle与W-cycle

多重网格通过在不同层级间递归来消除所有频率的误差。

**V-cycle算法**：
```
function V_cycle(A_h, u_h, f_h, level)
    if level == coarsest then
        u_h = A_h^{-1} f_h  // 直接求解
    else
        // 前光滑
        u_h = Smooth(A_h, u_h, f_h, ν₁)
        
        // 计算残差并限制
        r_h = f_h - A_h u_h
        r_{2h} = I_{2h}^h r_h
        
        // 粗网格修正
        e_{2h} = 0
        e_{2h} = V_cycle(A_{2h}, e_{2h}, r_{2h}, level+1)
        
        // 延拓并修正
        u_h = u_h + I_h^{2h} e_{2h}
        
        // 后光滑
        u_h = Smooth(A_h, u_h, f_h, ν₂)
    end
    return u_h
end
```

**W-cycle**：在每层递归调用两次
```
e_{2h} = W_cycle(A_{2h}, 0, r_{2h}, level+1)
e_{2h} = W_cycle(A_{2h}, e_{2h}, r_{2h}, level+1)
```

**F-cycle**：介于V和W之间的策略

**Full Multigrid (FMG)**：
使用粗网格解作为细网格初值：
1. 在最粗网格求解
2. 延拓到下一层作为初值
3. 执行V-cycle
4. 重复直到最细网格

### 5.4.4 粗网格修正

粗网格修正的数学原理基于误差方程：
$$A_h e_h = r_h$$

其中 $e_h = u_h^* - u_h$ 是误差，$r_h = f_h - A_h u_h$ 是残差。

**两网格算法分析**：
1. 光滑后的误差主要是低频成分
2. 限制到粗网格：$r_{2h} = I_{2h}^h r_h$
3. 求解粗网格误差方程：$A_{2h} e_{2h} = r_{2h}$
4. 延拓修正：$u_h^{new} = u_h + I_h^{2h} e_{2h}$

**收敛性分析**：
两网格收敛因子：
$$\rho_{TG} = ||(I - I_h^{2h} A_{2h}^{-1} I_{2h}^h A_h) S^{\nu}||$$

其中 $S$ 是光滑迭代矩阵。典型值：$\rho_{TG} \approx 0.1-0.2$。

## 5.5 几何多重网格

### 5.5.1 网格层次构建

几何多重网格需要构建一系列逐渐变粗的网格。对于结构化网格，这个过程相对简单。

**均匀粗化策略**：
- 一维：每隔一个点取一个（间距翻倍）
- 二维：每个方向都翻倍，4个细网格单元对应1个粗网格单元
- 三维：8个细网格单元对应1个粗网格单元

**网格层次示例**（二维）：
```
Level 0 (finest):   64×64   (h = 1/64)
Level 1:            32×32   (h = 1/32)
Level 2:            16×16   (h = 1/16)
Level 3:             8×8    (h = 1/8)
Level 4 (coarsest):  4×4    (h = 1/4)
```

**边界处理**：
- Dirichlet边界：粗网格继承细网格边界条件
- Neumann边界：需要特殊处理以保持通量守恒

**非均匀网格**：
对于自适应网格，粗化策略更复杂：
- 聚集(agglomeration)：将相邻细网格单元合并
- 确保粗网格的连通性和质量

### 5.5.2 Galerkin粗化

Galerkin粗化通过变分原理自动生成粗网格算子：
$$A_{2h} = I_{2h}^h A_h I_h^{2h}$$

**优点**：
- 保持对称性：如果 $A_h$ 对称且 $I_{2h}^h = (I_h^{2h})^T$，则 $A_{2h}$ 对称
- 保持正定性：如果 $A_h$ 正定，则 $A_{2h}$ 正定
- 自动处理复杂边界条件

**计算优化**：
直接矩阵三乘积计算开销大，可以优化：
1. 对于标准离散化，$A_{2h}$ 的模板可以预先推导
2. 利用稀疏性，只计算非零元素

**五点Laplace算子的Galerkin粗化**：
细网格算子：
$$A_h = \frac{1}{h^2}\begin{bmatrix}
0 & -1 & 0\\
-1 & 4 & -1\\
0 & -1 & 0
\end{bmatrix}$$

粗网格算子（全权重限制+双线性延拓）：
$$A_{2h} = \frac{1}{(2h)^2}\begin{bmatrix}
0 & -1 & 0\\
-1 & 4 & -1\\
0 & -1 & 0
\end{bmatrix}$$

注意：粗网格算子保持相同的模板形式！

### 5.5.3 光滑器选择

不同的光滑器对多重网格性能有重要影响。

**点Jacobi**：
$$u_i^{new} = \frac{1}{a_{ii}}(f_i - \sum_{j \neq i} a_{ij}u_j^{old})$$

- 完全并行，适合GPU
- 需要欠松弛 $\omega \approx 2/3$ 获得最佳光滑性
- 收敛较慢，通常需要更多迭代

**Gauss-Seidel**：
$$u_i^{new} = \frac{1}{a_{ii}}(f_i - \sum_{j < i} a_{ij}u_j^{new} - \sum_{j > i} a_{ij}u_j^{old})$$

- 串行依赖，并行性差
- 光滑效果好，收敛快
- 前向和后向GS的组合（对称GS）保持对称性

**Red-Black Gauss-Seidel**：
将网格点按棋盘模式分为红黑两组：
1. 更新所有红点（并行）
2. 更新所有黑点（并行）

- 保持GS的良好光滑性
- 实现并行化
- 特别适合规则网格

**线/面松弛**：
- **线松弛**：同时求解一条线上的所有未知数
- **面松弛**：同时求解一个面上的所有未知数
- 适用于各向异性问题（如 $\epsilon u_{xx} + u_{yy} = f$ 当 $\epsilon \ll 1$）

### 5.5.4 并行多重网格

多重网格的并行化面临独特挑战，特别是在粗网格层级。

**并行化策略**：

1. **网格分区**：
   - 细网格：良好的负载均衡
   - 粗网格：通信开销增加，计算/通信比下降
   
2. **粗网格并行性退化**：
   ```
   Level 0: 1024×1024 / 64 processors = 16384 points/proc
   Level 1: 512×512 / 64 processors = 4096 points/proc
   ...
   Level 6: 16×16 / 64 processors = 4 points/proc (!)
   ```

3. **解决方案**：
   - **聚集(Agglomeration)**：在粗网格上使用更少的处理器
   - **冗余计算**：每个处理器都计算粗网格问题
   - **混合方法**：粗网格切换到Krylov方法

**通信优化**：
- **重叠计算与通信**：在等待边界数据时计算内部点
- **消息合并**：将多个小消息合并为大消息
- **持久通信**：对于固定通信模式，使用MPI持久通信

**GPU实现考虑**：
- 细网格：高并行度，GPU效率高
- 粗网格：并行度低，可能需要CPU处理
- 使用统一内存简化数据传输

## 5.6 代数多重网格(AMG)

### 5.6.1 强连接与粗网格选择

代数多重网格不依赖几何信息，而是通过分析矩阵系数来构建网格层次。

**强连接定义**：
节点 $i$ 强连接到节点 $j$ 如果：
$$|a_{ij}| \geq \theta \max_{k \neq i} |a_{ik}|$$

典型选择 $\theta = 0.25$ 或 $0.5$。

**强连接的意义**：
- 大系数表示变量间的强耦合
- 光滑迭代难以消除强连接变量间的误差
- 需要在粗网格上同时处理强连接的变量

**C/F分裂算法**（Coarse/Fine splitting）：
```
1. 计算每个点的"影响度"：λᵢ = |{j: j强连接到i}|
2. 初始化：所有点标记为未定(U)
3. while 存在U点：
   a. 选择λ最大的U点i，标记为C点
   b. 将所有强连接到i的U点标记为F点
   c. 更新受影响点的λ值
4. 后处理：确保每个F点至少强连接到一个C点
```

**两遍算法**(Second pass)：
检查每个F点，如果它没有强连接到任何C点，则：
- 将其改为C点，或
- 将某个强连接的F点改为C点

### 5.6.2 插值算子构造

AMG的核心是构造好的插值(延拓)算子。基本原则是准确插值光滑误差。

**经典插值**：
对于F点 $i$，其值由强连接的C点插值：
$$e_i = \sum_{j \in C_i} w_{ij} e_j$$

其中 $C_i$ 是强连接到 $i$ 的C点集合。

**插值权重计算**：
基于光滑误差假设 $Ae \approx 0$：
$$a_{ii}e_i + \sum_{j \in N_i} a_{ij}e_j \approx 0$$

将邻居分为C点和F点：
$$e_i = -\frac{1}{a_{ii}} \left( \sum_{j \in C_i} a_{ij}e_j + \sum_{k \in F_i} a_{ik}e_k \right)$$

对F点的贡献进行近似（如分配到C点），得到插值权重。

**标准插值公式**：
$$w_{ij} = -\frac{a_{ij} + \sum_{k \in F_i} \frac{a_{ik}a_{kj}}{\sum_{m \in C_i} a_{km}}}{a_{ii} + \sum_{k \in F_i} \frac{a_{ik}a_{ki}}{a_{kk}}}$$

### 5.6.3 经典AMG与平滑聚合

**经典AMG特点**：
- 基于强连接的C/F分裂
- 直接插值构造
- 对各向异性问题效果好
- 构造复杂度较高

**平滑聚合(SA)AMG**：
不同于C/F分裂，SA使用聚合策略：

1. **聚合阶段**：
   - 将强连接的节点聚合成组
   - 每组形成一个粗网格节点
   - 初始插值：分片常数

2. **平滑阶段**：
   - 对初始插值应用光滑迭代
   - $P = S \tilde{P}$，其中 $S$ 是光滑算子
   - 改善插值质量

**聚合算法**：
```
1. 初始化：所有点未聚合
2. for each 未聚合点 i：
   if i 有足够的未聚合强连接邻居：
      创建新聚合，包含 i 和其强连接邻居
   else：
      将 i 加入邻近的聚合
3. 处理剩余点（可能需要放松强连接准则）
```

### 5.6.4 自适应AMG

自适应AMG通过测试向量自动改进插值算子。

**基本思想**：
使用实际问题的近零特征向量（代表光滑误差）来构造插值。

**自适应设置算法**：
1. 初始设置：使用标准AMG构造
2. 测试阶段：
   - 执行多重网格迭代
   - 收集收敛慢的误差分量
3. 改进插值：
   - 将测试向量加入插值的范围
   - 重新构造网格层次

**Bootstrap AMG**：
```
1. 从常数向量开始：v⁽⁰⁾ = 1
2. for k = 1 to K：
   a. 用当前向量构造AMG层次
   b. 应用AMG求解，记录收敛慢的分量
   c. 更新测试向量集：V⁽ᵏ⁾ = [v⁽⁰⁾, v⁽¹⁾, ..., v⁽ᵏ⁾]
3. 使用最终的向量集构造AMG
```

**优点**：
- 自动适应问题特性
- 对困难问题效果显著
- 可以处理多个近零特征向量

**缺点**：
- 设置阶段开销大
- 需要存储测试向量
- 对问题变化敏感

## 5.7 无矩阵(Matrix-free)方法

### 5.7.1 矩阵-向量乘积的隐式计算

在现代计算架构中，内存带宽往往是性能瓶颈。无矩阵方法通过直接计算矩阵-向量乘积而不存储矩阵来优化性能。

**动机**：
- 存储稀疏矩阵需要大量内存
- 矩阵元素的读取受内存带宽限制
- 重新计算往往比读取更快

**Laplace算子的无矩阵实现**：
```python
def laplace_matvec(u, h):
    # 不存储矩阵，直接计算 Au
    Au = zeros_like(u)
    inv_h2 = 1.0 / (h * h)
    
    for i in range(1, n-1):
        for j in range(1, n-1):
            Au[i,j] = inv_h2 * (
                4*u[i,j] - u[i-1,j] - u[i+1,j] 
                         - u[i,j-1] - u[i,j+1]
            )
    return Au
```

**性能分析**：
- 传统方法：读取5个矩阵元素 + 5个向量元素 = 10次内存访问
- 无矩阵方法：读取5个向量元素 = 5次内存访问
- 算术运算增加但通常"免费"（计算隐藏在内存访问延迟中）

### 5.7.2 内存带宽优化

现代处理器的关键特性：
- 计算能力 >> 内存带宽
- FLOP/字节比持续增加
- 缓存层次结构的重要性

**Roofline模型分析**：
算术强度(AI) = FLOPs / 内存字节数

对于Laplace算子：
- FLOPs: 5次乘法 + 4次加法 = 9 FLOPs
- 内存: 5次读取 × 8字节 = 40字节
- AI = 9/40 ≈ 0.225 FLOP/字节

这远低于现代CPU的平衡点（通常 > 10 FLOP/字节），说明是内存带宽受限。

**优化策略**：

1. **时间阻塞(Temporal blocking)**：
   ```python
   # 在缓存中多次使用数据
   for t_block in range(0, T, block_size):
       for t in range(t_block, min(t_block+block_size, T)):
           update_interior()
           exchange_boundaries()
   ```

2. **空间阻塞(Spatial blocking)**：
   ```python
   # 分块处理以提高缓存利用率
   for ii in range(0, n, block):
       for jj in range(0, n, block):
           for i in range(ii, min(ii+block, n)):
               for j in range(jj, min(jj+block, n)):
                   compute_stencil(i, j)
   ```

3. **数据布局优化**：
   - 使用SoA而非AoS提高向量化效率
   - 内存对齐提高SIMD效率

### 5.7.3 算子融合技术

算子融合通过合并多个操作减少内存访问次数。

**未融合的CG迭代**：
```python
# 5个独立的循环，5次内存遍历
q = A @ p           # 循环1
alpha = rTr / (p @ q)  # 循环2
x = x + alpha * p    # 循环3
r = r - alpha * q    # 循环4
rTr_new = r @ r      # 循环5
```

**融合后的实现**：
```python
def fused_cg_iteration(x, r, p, A_func):
    # 融合多个操作在一个循环中
    q = zeros_like(p)
    pq = 0.0
    rTr_new = 0.0
    
    for i in range(n):
        # 计算 q = Ap
        q[i] = A_func(p, i)
        # 累加 p·q
        pq += p[i] * q[i]
    
    alpha = rTr / pq
    
    for i in range(n):
        # 更新 x, r 并计算新的 r·r
        x[i] += alpha * p[i]
        r[i] -= alpha * q[i]
        rTr_new += r[i] * r[i]
    
    return x, r, q, rTr_new
```

**融合收益**：
- 减少内存遍历次数
- 提高时间局部性
- 更好的编译器优化机会

### 5.7.4 GPU实现策略

GPU上的无矩阵方法需要特殊考虑。

**线程映射策略**：
```cuda
__global__ void laplace_3d_kernel(
    const float* u, float* Au, 
    int nx, int ny, int nz, float inv_h2
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;
    
    if (i > 0 && i < nx-1 && 
        j > 0 && j < ny-1 && 
        k > 0 && k < nz-1) {
        
        int idx = i + nx*(j + ny*k);
        Au[idx] = inv_h2 * (
            6*u[idx] 
            - u[idx-1] - u[idx+1]
            - u[idx-nx] - u[idx+nx]
            - u[idx-nx*ny] - u[idx+nx*ny]
        );
    }
}
```

**共享内存优化**：
```cuda
__global__ void optimized_laplace_kernel(...) {
    __shared__ float tile[TILE_DIM+2][TILE_DIM+2];
    
    // 协作加载包含边界的数据块到共享内存
    load_tile_with_halo(tile, u, ...);
    __syncthreads();
    
    // 从共享内存计算
    if (threadIdx.x > 0 && threadIdx.x < TILE_DIM-1 &&
        threadIdx.y > 0 && threadIdx.y < TILE_DIM-1) {
        
        float result = inv_h2 * (
            4*tile[ty][tx]
            - tile[ty-1][tx] - tile[ty+1][tx]
            - tile[ty][tx-1] - tile[ty][tx+1]
        );
        
        // 写回全局内存
        Au[global_idx] = result;
    }
}
```

**性能考虑**：
1. **合并访问**：确保相邻线程访问相邻内存
2. **占用率**：平衡寄存器使用和块大小
3. **避免分支发散**：最小化条件语句
4. **使用纹理内存**：对于有空间局部性的访问模式

## 5.8 泊松方程的快速解法

### 5.8.1 格林函数与基本解

泊松方程 $\nabla^2 \phi = f$ 可以通过格林函数方法求解。

**基本解（三维）**：
$$G(x, y) = -\frac{1}{4\pi||x-y||}$$

**通解**：
$$\phi(x) = \int_\Omega G(x, y) f(y) dy + \int_{\partial\Omega} \left[ G(x, y)\frac{\partial\phi}{\partial n} - \phi(y)\frac{\partial G}{\partial n} \right] ds$$

对于自由空间（无边界）问题：
$$\phi(x) = -\frac{1}{4\pi} \int_\Omega \frac{f(y)}{||x-y||} dy$$

**离散形式**：
$$\phi_i = \sum_j G_{ij} f_j \Delta V_j$$

其中 $G_{ij} = -1/(4\pi||x_i - x_j||)$。

### 5.8.2 快速多极子方法(FMM)

直接计算所有粒子对相互作用需要 $O(N^2)$ 操作。FMM通过多极展开将复杂度降至 $O(N)$。

**核心思想**：
1. 远场使用多极展开近似
2. 近场直接计算
3. 层次结构加速计算

**多极展开**：
对于源点 $y$ 在原点附近，场点 $x$ 远离原点：
$$\frac{1}{||x-y||} \approx \sum_{l=0}^p \sum_{m=-l}^l \frac{Y_l^m(\hat{y})}{r^{l+1}} r_y^l Y_l^m(\hat{x})$$

其中 $Y_l^m$ 是球谐函数。

**FMM算法步骤**：
1. **构建树结构**：八叉树(3D)或四叉树(2D)
2. **上行遍历(Upward pass)**：
   - P2M: 粒子到多极
   - M2M: 多极到多极（子到父）
3. **下行遍历(Downward pass)**：
   - M2L: 多极到局部（远场作用）
   - L2L: 局部到局部（父到子）
4. **最终求值**：
   - L2P: 局部到粒子
   - P2P: 粒子到粒子（近场直接）

### 5.8.3 PPPM方法

Particle-Particle Particle-Mesh (PPPM) 方法结合了直接求和和网格方法。

**基本思想**：
$$\phi = \phi_{short} + \phi_{long}$$

- 短程部分：直接粒子-粒子相互作用
- 长程部分：通过FFT在网格上求解

**算法流程**：
1. **粒子到网格(P2G)**：
   $$\rho_{\text{grid}}(x) = \sum_i q_i W(x - x_i)$$
   
2. **网格上求解泊松方程**（FFT）：
   $$\hat{\phi}_k = \frac{\hat{\rho}_k}{k^2}$$
   
3. **网格到粒子(G2P)**：
   $$\phi_i^{\text{long}} = \sum_{\text{grid}} \phi_{\text{grid}} W(x_{\text{grid}} - x_i)$$
   
4. **短程修正**：
   $$\phi_i = \phi_i^{\text{long}} + \sum_{j \in \text{near}} \frac{q_j \text{erfc}(\alpha r_{ij})}{r_{ij}}$$

**参数选择**：
- $\alpha$：短程/长程分离参数
- 网格分辨率：平衡精度和效率
- 插值阶数：通常使用3次B样条

### 5.8.4 FFT方法(周期边界)

对于周期边界条件，泊松方程可以通过FFT高效求解。

**算法**：
1. **前向FFT**：
   $$\hat{f}_k = \text{FFT}(f)$$
   
2. **谱空间求解**：
   $$\hat{\phi}_k = \frac{\hat{f}_k}{-k^2}$$
   
   注意：$k = 0$ 模式需要特殊处理（设为0或平均值）
   
3. **逆向FFT**：
   $$\phi = \text{IFFT}(\hat{\phi})$$

**离散波数**：
对于 $N \times N$ 网格：
$$k_x = \frac{2\pi}{L} n_x, \quad n_x = 0, 1, ..., N-1$$

实际计算中使用：
$$k^2 = \frac{4}{h^2}\left[\sin^2\left(\frac{\pi n_x}{N}\right) + \sin^2\left(\frac{\pi n_y}{N}\right)\right]$$

**非周期边界的处理**：
1. **正弦变换**：适用于Dirichlet边界
2. **余弦变换**：适用于Neumann边界
3. **Chebyshev方法**：非均匀网格

**实现优化**：
```python
def poisson_fft_solve(f, h):
    # 前向FFT
    f_hat = np.fft.fftn(f)
    
    # 构造波数
    kx = np.fft.fftfreq(n, d=h) * 2 * np.pi
    ky = np.fft.fftfreq(n, d=h) * 2 * np.pi
    kx, ky = np.meshgrid(kx, ky, indexing='ij')
    
    # 离散Laplace算子的特征值
    k_squared = (2 - 2*np.cos(kx*h))/h**2 + \
                (2 - 2*np.cos(ky*h))/h**2
    
    # 避免除零
    k_squared[0, 0] = 1.0
    phi_hat = f_hat / (-k_squared)
    phi_hat[0, 0] = 0.0  # 设置平均值为0
    
    # 逆向FFT
    return np.fft.ifftn(phi_hat).real
```

## 本章小结

本章深入探讨了欧拉视角物理仿真中的核心计算问题——大规模稀疏线性系统的求解。我们学习了从基础迭代法到现代高性能算法的完整谱系：

**关键概念**：
1. **稀疏矩阵特性**：CSR/COO存储格式，零空间处理，兼容性条件
2. **Krylov子空间方法**：CG用于对称正定系统，BiCGSTAB用于非对称系统
3. **预条件技术**：Jacobi、IC、MIC等预条件器显著加速收敛
4. **多重网格方法**：利用误差的频率特性，实现 $O(n)$ 复杂度
5. **代数多重网格**：不依赖几何信息，自动构建网格层次
6. **无矩阵方法**：优化内存带宽使用，适应现代计算架构
7. **特殊快速算法**：FFT、FMM、PPPM等针对泊松方程的优化方法

**核心公式**：
- Krylov子空间：$\mathcal{K}_m(A, r_0) = \text{span}\{r_0, Ar_0, ..., A^{m-1}r_0\}$
- CG收敛率：$||e_k||_A \leq 2\left(\frac{\sqrt{\kappa} - 1}{\sqrt{\kappa} + 1}\right)^k ||e_0||_A$
- Galerkin粗化：$A_{2h} = I_{2h}^h A_h I_h^{2h}$
- 泊松方程基本解：$G(x, y) = -\frac{1}{4\pi||x-y||}$（3D）

## 练习题

### 基础题

**练习 5.1**：稀疏矩阵存储
给定5×5三对角矩阵：
$$A = \begin{bmatrix}
2 & -1 & 0 & 0 & 0\\
-1 & 2 & -1 & 0 & 0\\
0 & -1 & 2 & -1 & 0\\
0 & 0 & -1 & 2 & -1\\
0 & 0 & 0 & -1 & 2
\end{bmatrix}$$

写出其CSR格式表示。

<details>
<summary>提示</summary>
CSR需要三个数组：values（非零值）、col_indices（列索引）、row_ptrs（行指针）。
</details>

<details>
<summary>答案</summary>

CSR表示：
- values = [2, -1, -1, 2, -1, -1, 2, -1, -1, 2, -1, -1, 2]
- col_indices = [0, 1, 0, 1, 2, 1, 2, 3, 2, 3, 4, 3, 4]
- row_ptrs = [0, 2, 5, 8, 11, 13]

验证：第i行的非零元素在values[row_ptrs[i]:row_ptrs[i+1]]中。
</details>

**练习 5.2**：CG迭代
对于系统 $Ax = b$，其中：
$$A = \begin{bmatrix} 4 & 1 \\ 1 & 3 \end{bmatrix}, \quad b = \begin{bmatrix} 1 \\ 2 \end{bmatrix}$$

手工执行一步CG迭代，初始猜测 $x_0 = [0, 0]^T$。

<details>
<summary>提示</summary>
CG第一步：$r_0 = b - Ax_0$，$p_0 = r_0$，计算 $\alpha_0 = \frac{r_0^T r_0}{p_0^T A p_0}$。
</details>

<details>
<summary>答案</summary>

1. 初始残差：$r_0 = b - Ax_0 = [1, 2]^T$
2. 初始搜索方向：$p_0 = r_0 = [1, 2]^T$
3. 计算 $Ap_0 = [4·1 + 1·2, 1·1 + 3·2]^T = [6, 7]^T$
4. 步长：$\alpha_0 = \frac{r_0^T r_0}{p_0^T Ap_0} = \frac{1^2 + 2^2}{1·6 + 2·7} = \frac{5}{20} = 0.25$
5. 更新解：$x_1 = x_0 + \alpha_0 p_0 = [0, 0]^T + 0.25[1, 2]^T = [0.25, 0.5]^T$
6. 更新残差：$r_1 = r_0 - \alpha_0 Ap_0 = [1, 2]^T - 0.25[6, 7]^T = [-0.5, 0.25]^T$
</details>

**练习 5.3**：多重网格限制算子
对于一维网格上的向量 $u_h = [1, 4, 2, 5, 3]$（5个点），使用全权重限制算子计算粗网格向量 $u_{2h}$。

<details>
<summary>提示</summary>
一维全权重限制：$u_{2h,i} = \frac{1}{4}u_{h,2i-1} + \frac{1}{2}u_{h,2i} + \frac{1}{4}u_{h,2i+1}$。
</details>

<details>
<summary>答案</summary>

粗网格有3个点（索引0, 1, 2对应细网格的0, 2, 4）：
- $u_{2h,0} = u_{h,0} = 1$（边界点直接复制）
- $u_{2h,1} = \frac{1}{4}u_{h,1} + \frac{1}{2}u_{h,2} + \frac{1}{4}u_{h,3} = \frac{1}{4}·4 + \frac{1}{2}·2 + \frac{1}{4}·5 = 3.25$
- $u_{2h,2} = u_{h,4} = 3$（边界点直接复制）

因此 $u_{2h} = [1, 3.25, 3]$。
</details>

**练习 5.4**：AMG强连接
给定矩阵行：$a_i = [0, -0.1, 0, -0.8, 2.0, -0.5, 0, -0.2, -0.4]$，使用阈值 $\theta = 0.25$，确定节点 $i$ 强连接到哪些节点。

<details>
<summary>提示</summary>
节点 $j$ 是强连接如果 $|a_{ij}| \geq \theta \max_{k \neq i} |a_{ik}|$。
</details>

<details>
<summary>答案</summary>

1. 找出最大非对角元素：$\max_{k \neq i} |a_{ik}| = 0.8$（对应节点3）
2. 阈值：$\theta \max = 0.25 × 0.8 = 0.2$
3. 强连接判断：
   - 节点1：$|a_{i1}| = 0.1 < 0.2$ ✗
   - 节点3：$|a_{i3}| = 0.8 \geq 0.2$ ✓
   - 节点5：$|a_{i5}| = 0.5 \geq 0.2$ ✓
   - 节点7：$|a_{i7}| = 0.2 \geq 0.2$ ✓
   - 节点8：$|a_{i8}| = 0.4 \geq 0.2$ ✓

节点 $i$ 强连接到节点 {3, 5, 7, 8}。
</details>

### 挑战题

**练习 5.5**：条件数与CG收敛
证明：对于条件数为 $\kappa$ 的SPD矩阵，CG方法在 $k = O(\sqrt{\kappa} \log(1/\epsilon))$ 步内将相对误差降至 $\epsilon$。

<details>
<summary>提示</summary>
使用CG收敛估计 $||e_k||_A \leq 2\left(\frac{\sqrt{\kappa} - 1}{\sqrt{\kappa} + 1}\right)^k ||e_0||_A$。
</details>

<details>
<summary>答案</summary>

从收敛估计出发：
$$\frac{||e_k||_A}{||e_0||_A} \leq 2\left(\frac{\sqrt{\kappa} - 1}{\sqrt{\kappa} + 1}\right)^k$$

要使相对误差 $\leq \epsilon$，需要：
$$2\left(\frac{\sqrt{\kappa} - 1}{\sqrt{\kappa} + 1}\right)^k \leq \epsilon$$

取对数：
$$k \log\left(\frac{\sqrt{\kappa} - 1}{\sqrt{\kappa} + 1}\right) \leq \log(\epsilon/2)$$

当 $\kappa \gg 1$ 时：
$$\frac{\sqrt{\kappa} - 1}{\sqrt{\kappa} + 1} = \frac{1 - 1/\sqrt{\kappa}}{1 + 1/\sqrt{\kappa}} \approx 1 - \frac{2}{\sqrt{\kappa}}$$

因此：
$$\log\left(\frac{\sqrt{\kappa} - 1}{\sqrt{\kappa} + 1}\right) \approx \log\left(1 - \frac{2}{\sqrt{\kappa}}\right) \approx -\frac{2}{\sqrt{\kappa}}$$

代入得：
$$k \cdot \left(-\frac{2}{\sqrt{\kappa}}\right) \leq \log(\epsilon/2)$$
$$k \geq \frac{\sqrt{\kappa}}{2} \log(2/\epsilon)$$

因此 $k = O(\sqrt{\kappa} \log(1/\epsilon))$。
</details>

**练习 5.6**：多重网格复杂度分析
对于 $n × n$ 二维网格，证明V-cycle多重网格的计算复杂度是 $O(n^2)$。假设每层使用固定次数的Jacobi光滑。

<details>
<summary>提示</summary>
计算所有层级的工作量总和，利用几何级数求和。
</details>

<details>
<summary>答案</summary>

设最细网格有 $N = n^2$ 个点。

各层网格点数：
- Level 0: $n^2$
- Level 1: $(n/2)^2 = n^2/4$
- Level 2: $(n/4)^2 = n^2/16$
- ...
- Level L: $O(1)$

每层工作量正比于该层点数。设每个点的工作量为 $c$（光滑迭代）。

V-cycle总工作量：
$$W = c \sum_{l=0}^L \frac{n^2}{4^l} = cn^2 \sum_{l=0}^L \frac{1}{4^l}$$

这是几何级数，和为：
$$\sum_{l=0}^{\infty} \frac{1}{4^l} = \frac{1}{1-1/4} = \frac{4}{3}$$

因此：
$$W = cn^2 \cdot \frac{4}{3} = O(n^2)$$

这证明了V-cycle的线性复杂度（相对于未知数个数）。
</details>

**练习 5.7**：无矩阵方法的内存带宽分析
对于三维Laplace算子（7点模板），比较传统稀疏矩阵方法和无矩阵方法的内存带宽需求。假设使用双精度浮点数。

<details>
<summary>提示</summary>
计算每个矩阵-向量乘积需要的内存读写量。
</details>

<details>
<summary>答案</summary>

**传统稀疏矩阵方法**（CSR格式）：
- 读取：7个矩阵值 + 7个列索引 + 2个行指针 + 7个向量元素
- 写入：1个结果元素
- 总计：7×8 + 7×4 + 2×4 + 7×8 + 1×8 = 156字节/点

**无矩阵方法**：
- 读取：7个向量元素
- 写入：1个结果元素
- 总计：7×8 + 1×8 = 64字节/点

内存带宽减少：$\frac{156-64}{156} \approx 59\%$

**算术强度比较**：
- 传统：AI = 13 FLOPs / 156 bytes ≈ 0.083 FLOP/byte
- 无矩阵：AI = 13 FLOPs / 64 bytes ≈ 0.203 FLOP/byte

无矩阵方法的算术强度提高约2.4倍，更适合现代处理器。
</details>

**练习 5.8**：FFT求解泊松方程的误差分析
使用FFT方法求解周期边界的二维泊松方程时，为什么需要特殊处理 $k = 0$ 模式？这如何影响解的唯一性？

<details>
<summary>提示</summary>
考虑泊松方程的兼容性条件和周期边界下的零空间。
</details>

<details>
<summary>答案</summary>

**数学分析**：

泊松方程：$\nabla^2 u = f$

Fourier变换后：$-k^2 \hat{u}_k = \hat{f}_k$

当 $k = 0$ 时：
- 左边：$-0^2 \cdot \hat{u}_0 = 0$
- 右边：$\hat{f}_0 = \int_\Omega f dx$

**兼容性条件**：
有解的必要条件是 $\hat{f}_0 = 0$，即：
$$\int_\Omega f dx = 0$$

这对应于物理上的守恒条件（如不可压缩流体的质量守恒）。

**零空间**：
$k = 0$ 模式对应常数函数，是Laplace算子在周期边界下的零空间。

**处理方法**：
1. 检查兼容性：验证 $\sum_{ij} f_{ij} \approx 0$
2. 设置 $\hat{u}_0 = 0$：固定解的平均值为0
3. 或者：添加约束 $\int u dx = c$

**物理意义**：
在周期边界下，泊松方程只能确定势函数的梯度，不能确定绝对值。这反映了规范不变性——物理量（如速度、电场）只依赖于势的梯度。
</details>

## 常见陷阱与错误

1. **零空间处理不当**：
   - 错误：直接求解奇异系统导致发散
   - 正确：投影到零空间正交补或固定一点

2. **预条件器选择**：
   - 错误：对所有问题都用Jacobi预条件
   - 正确：根据问题特性选择合适预条件器

3. **多重网格粗网格**：
   - 错误：粗网格太粗导致不能表示低频误差
   - 正确：保证粗网格仍能解析问题的主要特征

4. **Krylov方法数值稳定性**：
   - 错误：忽视正交性损失导致停滞
   - 正确：使用重正交化或更稳定的变种

5. **无矩阵方法的边界处理**：
   - 错误：硬编码边界条件导致不灵活
   - 正确：设计通用的边界处理框架

## 最佳实践检查清单

### 求解器选择
- [ ] 矩阵是否对称正定？→ 使用CG
- [ ] 矩阵是否非对称？→ 使用BiCGSTAB或GMRES
- [ ] 问题规模是否很大？→ 考虑多重网格
- [ ] 是否有好的预条件器？→ 使用预条件Krylov方法
- [ ] 内存是否受限？→ 使用无矩阵方法

### 性能优化
- [ ] 是否分析了矩阵的稀疏模式？
- [ ] 是否测试了不同的预条件器？
- [ ] 是否优化了矩阵-向量乘积？
- [ ] 是否考虑了并行化？
- [ ] 是否利用了问题的特殊结构？

### 数值稳定性
- [ ] 是否处理了零空间？
- [ ] 是否设置了合理的收敛准则？
- [ ] 是否监控了迭代历史？
- [ ] 是否有备用求解策略？
- [ ] 是否验证了解的正确性？