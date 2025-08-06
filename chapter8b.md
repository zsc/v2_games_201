# 第八章：多重网格方法

多重网格方法是求解大规模线性系统最高效的算法之一，特别适用于椭圆型偏微分方程离散化产生的系统。本章将深入探讨多重网格的理论基础、实现细节和在物理仿真中的应用。通过学习本章，读者将掌握设计和实现高性能多重网格求解器的核心技术。

**学习目标**：
- 理解误差的频率分析和光滑性质
- 掌握几何多重网格和代数多重网格的实现
- 学会选择合适的光滑器和传输算子
- 能够在GPU上实现并行多重网格算法
- 将多重网格应用于流体压力投影和弹性力学问题

---

## 8.1 多重网格方法基础

### 8.1.1 迭代方法的局限性

考虑线性系统 $Au = f$，传统迭代方法如Jacobi和Gauss-Seidel的收敛速度依赖于误差的频率成分。定义误差 $e = u - u^*$，其中 $u^*$ 是精确解。

对于一维泊松方程 $-u_{xx} = f$，在均匀网格上离散化得到：
$$\frac{-u_{i-1} + 2u_i - u_{i+1}}{h^2} = f_i$$

误差可以展开为傅里叶级数：
$$e(x) = \sum_{k=1}^{n-1} \alpha_k \sin(k\pi x/L)$$

其中高频成分（大的 $k$）对应快速振荡的误差，低频成分对应缓慢变化的误差。

### 8.1.2 光滑性质

Jacobi迭代的误差传播算子为：
$$G_{Jacobi} = I - \omega D^{-1}A$$

对于模式 $k$，其衰减因子为：
$$\mu_k = 1 - \omega \frac{\lambda_k}{\lambda_{max}}$$

其中 $\lambda_k = 4\sin^2(k\pi h/2)/h^2$ 是矩阵 $A$ 的第 $k$ 个特征值。

**关键观察**：高频误差（$k$ 接近 $n/2$）快速衰减，但低频误差（小的 $k$）衰减缓慢。这就是光滑性质——简单迭代方法能快速消除高频误差，但对低频误差效果差。

### 8.1.3 粗网格修正原理

多重网格的核心思想：在粗网格上，原本的低频误差变成了相对高频的误差，因此可以被高效消除。

设细网格步长为 $h$，粗网格步长为 $2h$。在细网格上频率为 $k$ 的模式，在粗网格上对应频率 $2k$，相对频率提高了一倍。

### 8.1.4 两网格算法

基本的两网格算法流程：
1. 在细网格上进行 $\nu_1$ 次前光滑（如Gauss-Seidel）
2. 计算残差 $r^h = f^h - A^h u^h$
3. 限制残差到粗网格：$r^{2h} = I_{2h}^h r^h$
4. 在粗网格上求解：$A^{2h} e^{2h} = r^{2h}$
5. 延拓误差到细网格：$e^h = I_h^{2h} e^{2h}$
6. 修正解：$u^h \leftarrow u^h + e^h$
7. 在细网格上进行 $\nu_2$ 次后光滑

### 8.1.5 收敛性分析

两网格方法的误差传播算子：
$$M_{TG} = S^{\nu_2} (I - I_h^{2h} (A^{2h})^{-1} I_{2h}^h A^h) S^{\nu_1}$$

其中 $S$ 是光滑器的误差传播算子。

**收敛条件**：$\rho(M_{TG}) < 1$，其中 $\rho$ 表示谱半径。

对于适当选择的光滑器和传输算子，可以证明：
$$\|M_{TG}\| \leq C \cdot h^{\alpha}$$

其中 $\alpha > 0$，表明网格越细，收敛越快。

### 8.1.6 多重网格的动机

两网格方法需要在粗网格上精确求解，当问题规模大时仍然昂贵。解决方案：递归应用两网格思想，形成多重网格。

在最粗网格上，问题规模足够小，可以用直接法求解（如LU分解）。

---

## 8.2 几何多重网格

### 8.2.1 网格层次构建

对于 $d$ 维问题，构建网格序列：
$$\Omega^{h_0} \supset \Omega^{h_1} \supset ... \supset \Omega^{h_L}$$

其中 $h_{l+1} = 2h_l$，每个方向上网格点数减半。

**网格粗化策略**：
- 标准粗化：每个方向均匀粗化因子为2
- 半粗化：只在某些方向粗化（用于各向异性问题）
- 自适应粗化：根据问题特性选择粗化区域

### 8.2.2 限制算子

限制算子 $I_{2h}^h: \Omega^h \rightarrow \Omega^{2h}$ 将细网格函数映射到粗网格。

**1D情况**：
- 注入(Injection)：$u_{2i}^{2h} = u_{2i}^h$
- 全权重(Full weighting)：$u_i^{2h} = \frac{1}{4}(u_{2i-1}^h + 2u_{2i}^h + u_{2i+1}^h)$

**2D情况**（9点模板）：
$$u_{i,j}^{2h} = \frac{1}{16} \begin{bmatrix} 1 & 2 & 1 \\ 2 & 4 & 2 \\ 1 & 2 & 1 \end{bmatrix} \ast u^h$$

### 8.2.3 延拓算子

延拓算子 $I_h^{2h}: \Omega^{2h} \rightarrow \Omega^h$ 将粗网格函数插值到细网格。

**线性插值（1D）**：
- 粗网格点直接复制：$u_{2i}^h = u_i^{2h}$
- 中间点线性插值：$u_{2i+1}^h = \frac{1}{2}(u_i^{2h} + u_{i+1}^{2h})$

**双线性插值（2D）**：
对于细网格点 $(2i+p, 2j+q)$，其中 $p,q \in \{0,1\}$：
$$u_{2i+p,2j+q}^h = \sum_{k,l \in \{0,1\}} w_{p,k} w_{q,l} u_{i+k,j+l}^{2h}$$

其中权重 $w_{0,0} = w_{1,1} = 1$，$w_{0,1} = w_{1,0} = 0.5$。

### 8.2.4 Galerkin粗化

粗网格算子通过Galerkin条件构造：
$$A^{2h} = I_{2h}^h A^h I_h^{2h}$$

**性质**：
- 如果 $A^h$ 对称正定，则 $A^{2h}$ 也对称正定
- 保持变分性质
- 自动满足能量最小化原理

**计算优化**：对于标准离散化，可以直接推导粗网格模板，避免矩阵乘法。

### 8.2.5 边界条件处理

**Dirichlet边界**：
- 细网格边界值直接限制到粗网格
- 延拓时保持边界值不变

**Neumann边界**：
- 使用修正的限制/延拓模板
- 保证通量守恒

**混合边界**：分区域处理，确保一致性。

### 8.2.6 变系数问题

对于变系数泊松方程 $-\nabla \cdot (a(x)\nabla u) = f$：

**算术平均**：
$$a_{i+1/2}^{2h} = \frac{1}{2}(a_{2i+1/2}^h + a_{2i+3/2}^h)$$

**调和平均**（更适合强间断）：
$$\frac{1}{a_{i+1/2}^{2h}} = \frac{1}{2}\left(\frac{1}{a_{2i+1/2}^h} + \frac{1}{a_{2i+3/2}^h}\right)$$

---

## 8.3 代数多重网格(AMG)

### 8.3.1 AMG的动机

几何多重网格需要：
- 结构化网格
- 几何信息
- 规则的粗化策略

AMG只需要矩阵 $A$，适用于：
- 非结构网格
- 复杂几何
- 各向异性问题

### 8.3.2 强连接与影响

定义强连接：对于给定阈值 $\theta \in (0,1)$，如果
$$|a_{ij}| \geq \theta \max_{k \neq i} |a_{ik}|$$
则称 $i$ 强依赖于 $j$，记作 $j \in S_i$。

**强连接的意义**：
- 表示变量间的强耦合
- 指导粗网格点选择
- 决定插值权重

### 8.3.3 粗网格点选择

**Classical AMG的C/F分裂**：
1. 计算每个点的影响度量：$\lambda_i = |S_i^T|$（被多少点强依赖）
2. 选择 $\lambda$ 最大的点作为C点
3. 将其强连接的邻居标记为F点
4. 更新剩余点的 $\lambda$ 值
5. 重复直到所有点被分类

**性质要求**：
- F点应被足够的C点强连接
- C点集应该是极大独立集
- C点分布相对均匀

### 8.3.4 插值算子构造

**经典插值**：对F点 $i$，其插值公式：
$$u_i = \sum_{j \in C_i} w_{ij} u_j$$

其中 $C_i$ 是强连接到 $i$ 的C点集合。

**权重计算**：
$$w_{ij} = -\frac{a_{ij} + \sum_{k \in F_i^s} a_{ik} \bar{a}_{kj}}{a_{ii} + \sum_{k \in F_i^w} a_{ik}}$$

其中：
- $F_i^s$：强连接的F点
- $F_i^w$：弱连接的F点  
- $\bar{a}_{kj}$：从F点到C点的平均连接

### 8.3.5 平滑聚合AMG

**聚合阶段**：
1. 构建聚合 $\{A_k\}$，每个聚合包含强连接的点
2. 初始延拓算子：$(P^0)_{ij} = 1$ 如果 $i \in A_j$，否则为0

**平滑阶段**：
$$P = (I - \omega D^{-1}A) P^0$$

其中 $\omega = 2/3 \cdot 1/\rho(D^{-1}A)$。

**优点**：
- 构造简单
- 并行性好
- 对各向异性问题鲁棒

### 8.3.6 自适应AMG

**Bootstrap AMG**：
1. 用当前AMG求解测试问题
2. 分析收敛慢的误差分量
3. 将这些分量加入插值算子
4. 重新构造AMG层次

**收敛性监测**：
$$\rho_{est} = \|e^{(k)}\| / \|e^{(k-1)}\|$$

如果 $\rho_{est}$ 过大，触发自适应。

---

## 8.4 光滑器与限制/延拓算子

### 8.4.1 光滑器的选择准则

光滑器的作用是消除高频误差分量。理想的光滑器应该：
- 计算成本低
- 并行性好
- 对高频误差有强阻尼作用
- 不放大低频误差

**光滑因子**：定义为高频误差的最大放大率
$$\mu = \max_{k > n/2} |\mu_k|$$

其中 $\mu_k$ 是模式 $k$ 的衰减因子。

### 8.4.2 点光滑器

**Jacobi迭代**：
$$u_i^{new} = u_i^{old} + \omega \frac{r_i}{a_{ii}}$$

其中 $r_i = f_i - \sum_j a_{ij}u_j^{old}$。

- 最优松弛因子：$\omega = 2/3$（对于泊松方程）
- 光滑因子：$\mu \approx 0.6$
- 并行性：完美并行

**Gauss-Seidel迭代**：
$$u_i^{new} = \frac{1}{a_{ii}}\left(f_i - \sum_{j<i} a_{ij}u_j^{new} - \sum_{j>i} a_{ij}u_j^{old}\right)$$

- 光滑因子：$\mu \approx 0.5$
- 并行性：需要着色或分块策略

**红黑Gauss-Seidel**：
1. 更新红点（棋盘着色）
2. 更新黑点

优点：保持Gauss-Seidel的收敛性，实现并行。

### 8.4.3 块光滑器

**块Jacobi**：将未知量分组，每组内部用直接法求解
$$U_I^{new} = U_I^{old} + \omega A_{II}^{-1} R_I$$

其中 $I$ 表示第 $I$ 个块。

**线光滑器**：特别适合各向异性问题
- x-line：沿x方向的所有点组成一个块
- y-line：沿y方向的所有点组成一个块
- 交替方向：先x-line再y-line

**ILU光滑器**：不完全LU分解
$$A \approx LU$$
保持稀疏模式，作为预条件子。

### 8.4.4 多项式光滑器

**Chebyshev多项式**：
$$u^{(k+1)} = u^{(k)} + \alpha_k r^{(k)} + \beta_k (u^{(k)} - u^{(k-1)})$$

参数选择基于谱半径估计：
$$\alpha_k = \frac{4}{3\rho}, \quad \beta_k = \left(\frac{\rho - 2}{\rho + 2}\right)^2$$

优点：
- 不需要矩阵元素，只需矩阵向量乘积
- 可以精确控制阻尼特性
- 适合无矩阵方法

### 8.4.5 限制/延拓算子的优化

**能量最小化延拓**：
$$P = \arg\min_{\tilde{P}} \|(I - \tilde{P}A_c^{-1}\tilde{P}^T A)\|_A$$

其中 $\|\cdot\|_A$ 是能量范数。

**保持常数延拓**：确保 $P \mathbf{1} = \mathbf{1}$，保证常数函数精确传输。

**高阶延拓**：
- 三次插值：使用更多邻居点
- 保持多项式：精确传输低阶多项式

### 8.4.6 算子依赖的传输

**基于矩阵的插值**：
$$P_{ij} = \begin{cases}
1 & \text{if } i \text{ is coarse point } j \\
-\frac{\sum_{k \in N_i^C} a_{ik} P_{kj}}{a_{ii}} & \text{if } i \text{ is fine point}
\end{cases}$$

其中 $N_i^C$ 是 $i$ 的粗网格邻居。

**优点**：
- 自动适应系数变化
- 处理各向异性
- 保持矩阵性质（如M矩阵）

---

## 8.5 V循环、W循环与完全多重网格

### 8.5.1 多重网格循环策略

**递归定义**：
```
MGM(A_l, u_l, f_l, γ):
  if l == L (最粗层):
    直接求解 A_L u_L = f_L
  else:
    前光滑 ν₁ 次
    r_l = f_l - A_l u_l
    r_{l+1} = I_{l+1}^l r_l
    e_{l+1} = 0
    for i = 1 to γ:
      MGM(A_{l+1}, e_{l+1}, r_{l+1}, γ)
    e_l = I_l^{l+1} e_{l+1}
    u_l = u_l + e_l
    后光滑 ν₂ 次
```

其中 $\gamma$ 决定循环类型。

### 8.5.2 V循环（γ=1）

**计算复杂度**：设细网格有 $N$ 个未知量
- 工作量：$W_V = O(N)$（最优）
- 存储：$S_V = O(N)$

**收敛因子估计**：
$$\rho_V \approx 1 - O(h^2)$$

对于模型问题，典型值 $\rho_V \approx 0.1$。

**V循环的特点**：
- 每层只访问一次
- 适合作为预条件子
- 对初值敏感

### 8.5.3 W循环（γ=2）

**动机**：在粗网格上做更多工作，更彻底地消除误差。

**计算复杂度**：
- 2D：$W_W = O(N)$
- 3D：$W_W = O(N\log N)$

**收敛因子**：
$$\rho_W \approx \rho_V^2$$

更鲁棒但计算量更大。

### 8.5.4 F循环

F循环是V循环和W循环的折中：
```
第一次下降：V循环方式
之后：W循环方式
```

平衡了效率和鲁棒性。

### 8.5.5 完全多重网格（FMG）

**思想**：使用粗网格解作为细网格的初始猜测。

**算法**：
```
FMG(f):
  if 最粗层:
    直接求解
  else:
    f_{coarse} = restrict(f)
    u_{coarse} = FMG(f_{coarse})
    u = interpolate(u_{coarse})
    u = V-cycle(u, f)  // 一次或多次
  return u
```

**误差估计**：
$$\|u - u^*\|_A \leq C h^p$$

其中 $p$ 是离散化阶数。

### 8.5.6 循环策略选择

**V循环适用于**：
- 规则问题
- 作为Krylov方法的预条件子
- 内存受限情况

**W循环适用于**：
- 困难问题（各向异性、不连续系数）
- 需要高鲁棒性
- 粗网格代价低

**FMG适用于**：
- 需要高精度解
- 一次性求解（非迭代环境）
- 有好的初值估计

---

## 8.6 并行多重网格实现

### 8.6.1 并行化挑战

**细网格**：
- 大量并行度
- 通信/计算比低
- 负载均衡容易

**粗网格**：
- 并行度降低
- 通信开销相对增加
- 可能成为瓶颈

### 8.6.2 数据分布策略

**标准分布**：每层独立分区
```python
for level in range(L):
    partition[level] = domain_decomposition(grid[level])
```

**聚合策略**：粗网格聚合到部分处理器
```python
if grid_size[level] < threshold:
    use_subset_processors(level)
```

**优点**：减少粗网格通信，提高缓存利用率。

### 8.6.3 并行光滑器

**Jacobi**：自然并行
```python
@ti.kernel
def jacobi_smooth(u: ti.field, f: ti.field, omega: float):
    for i, j in u:
        if not_boundary(i, j):
            u_new[i,j] = u[i,j] + omega * residual(i,j) / a_diag
```

**多色Gauss-Seidel**：
```python
@ti.kernel  
def red_black_gs(u: ti.field, f: ti.field, color: int):
    for i, j in u:
        if (i + j) % 2 == color and not_boundary(i, j):
            u[i,j] = (f[i,j] - off_diagonal_sum(i,j)) / a_diag
```

**并行线光滑**：流水线方式或循环缩减。

### 8.6.4 GPU实现优化

**内存合并访问**：
```python
# 结构数组(SoA)布局
@ti.kernel
def restrict_2d(r_fine: ti.field, r_coarse: ti.field):
    for i, j in r_coarse:
        r_coarse[i,j] = 0.25 * r_fine[2*i, 2*j] + \
                        0.125 * (r_fine[2*i+1, 2*j] + r_fine[2*i-1, 2*j] + 
                                r_fine[2*i, 2*j+1] + r_fine[2*i, 2*j-1]) + \
                        0.0625 * (r_fine[2*i+1, 2*j+1] + r_fine[2*i-1, 2*j+1] +
                                 r_fine[2*i+1, 2*j-1] + r_fine[2*i-1, 2*j-1])
```

**共享内存利用**：
```python
@ti.kernel
def block_smooth(u: ti.field, f: ti.field):
    ti.block_local(u_local, BLOCK_SIZE+2, BLOCK_SIZE+2)
    # 加载到共享内存，包括halo
    # 在共享内存中迭代
    # 写回全局内存
```

### 8.6.5 通信优化

**非阻塞通信**：
```python
# 重叠计算和通信
def parallel_smooth():
    # 启动边界通信
    start_halo_exchange()
    
    # 计算内部点
    smooth_interior()
    
    # 等待通信完成
    wait_halo_exchange()
    
    # 计算边界点
    smooth_boundary()
```

**通信避免**：
- 冗余计算换通信
- 使用更宽的halo
- 粗网格复制

### 8.6.6 负载均衡

**动态负载均衡**：
```python
def adaptive_distribution(level):
    work_estimate = estimate_work(level)
    if work_estimate < threshold:
        # 聚合到更少的处理器
        new_procs = min(nprocs, work_estimate // min_work)
        redistribute(level, new_procs)
```

**空间填充曲线**：使用Hilbert或Z-order曲线保持局部性。

---

## 8.7 物理仿真中的应用

### 8.7.1 不可压缩流体的压力投影

**泊松方程**：在Chorin投影法中，压力满足：
$$\nabla^2 p = \frac{\rho}{\Delta t} \nabla \cdot u^*$$

**离散化**：在MAC网格上
$$\frac{p_{i+1,j} - 2p_{i,j} + p_{i-1,j}}{\Delta x^2} + \frac{p_{i,j+1} - 2p_{i,j} + p_{i,j-1}}{\Delta y^2} = \frac{\rho}{\Delta t} \left(\frac{u^*_{i+1/2,j} - u^*_{i-1/2,j}}{\Delta x} + \frac{v^*_{i,j+1/2} - v^*_{i,j-1/2}}{\Delta y}\right)$$

**多重网格设置**：
- 光滑器：红黑Gauss-Seidel（保持对称性）
- 限制/延拓：全权重限制，线性延拓
- 边界条件：Neumann边界（$\partial p/\partial n = 0$）

**奇异性处理**：
- 纯Neumann边界导致奇异系统
- 解决方案：固定一点压力值或投影到零均值空间

### 8.7.2 弹性力学问题

**线弹性方程**：
$$-\nabla \cdot \sigma = f$$
$$\sigma = 2\mu\epsilon + \lambda(\nabla \cdot u)I$$
$$\epsilon = \frac{1}{2}(\nabla u + \nabla u^T)$$

**块结构**：2D问题每个节点有2个自由度$(u,v)$
$$\begin{bmatrix} A_{uu} & A_{uv} \\ A_{vu} & A_{vv} \end{bmatrix} \begin{bmatrix} u \\ v \end{bmatrix} = \begin{bmatrix} f_x \\ f_y \end{bmatrix}$$

**多重网格考虑**：
- 使用块光滑器保持耦合
- Vanka光滑器：同时更新一个单元的所有自由度
- 各向异性问题需要线光滑或ILU

### 8.7.3 热传导与扩散

**隐式时间离散**：
$$\frac{u^{n+1} - u^n}{\Delta t} = \alpha \nabla^2 u^{n+1} + f$$

整理得：
$$(I - \alpha \Delta t \nabla^2) u^{n+1} = u^n + \Delta t f$$

**多重网格效率**：
- 大时间步（$\Delta t$ 大）使系统更接近泊松方程
- 小时间步系统接近恒等算子，简单迭代即可
- 典型策略：V(2,2)循环作为时间步进器

### 8.7.4 相场方法

**Cahn-Hilliard方程**：
$$\frac{\partial \phi}{\partial t} = \nabla \cdot (M \nabla \mu)$$
$$\mu = f'(\phi) - \epsilon^2 \nabla^2 \phi$$

**算子分裂**：
1. 求解化学势：$(I - \epsilon^2 \nabla^2)\mu = f'(\phi)$
2. 更新相场：$\phi^{n+1} = \phi^n + \Delta t \nabla \cdot (M \nabla \mu)$

两步都可用多重网格加速。

### 8.7.5 复杂边界处理

**浸入边界法**：
- 使用笛卡尔网格，边界切割单元
- 修正离散化模板近边界处
- 多重网格需要特殊处理切割单元

**自适应网格**：
- 局部加密需要的区域
- 多重网格在非均匀网格上的修正
- 限制/延拓算子的特殊处理

### 8.7.6 多物理场耦合

**流固耦合**：
```
while not converged:
    # 固体子问题（多重网格）
    solve_elasticity(u_solid, f_interface)
    
    # 流体子问题（多重网格）  
    solve_fluid(u_fluid, p, bc_from_solid)
    
    # 界面力更新
    update_interface_forces()
```

**单片式求解**：将耦合系统整体用多重网格
- 需要合适的块光滑器
- 物理量尺度差异需要预条件

---

## 8.8 高级主题与优化

### 8.8.1 自适应多重网格

**误差估计**：
$$\eta_K = h_K \|r_K\|_{L^2(K)}$$

其中 $K$ 是单元，$r_K$ 是局部残差。

**自适应策略**：
1. 计算误差指示子
2. 标记需要加密的区域
3. 局部加密网格
4. 更新多重网格层次

**挑战**：
- 非嵌套网格的传输算子
- 悬挂节点处理
- 负载均衡

### 8.8.2 无矩阵多重网格

**矩阵向量乘积**：
```python
@ti.kernel
def apply_stencil(u: ti.field, Au: ti.field):
    for i, j in u:
        Au[i,j] = stencil_center * u[i,j] + \
                  stencil_x * (u[i+1,j] + u[i-1,j]) + \
                  stencil_y * (u[i,j+1] + u[i,j-1])
```

**优势**：
- 减少内存需求
- 提高缓存效率
- 适合GPU实现

**几何信息**：需要保存
- 网格层次结构
- 系数场（如果是变系数）
- 边界标记

### 8.8.3 多重网格预条件

**与Krylov方法结合**：
```python
def pcg_with_multigrid(A, b, x0):
    r = b - A @ x0
    z = multigrid_v_cycle(r)  # 预条件
    p = z
    
    while norm(r) > tol:
        Ap = A @ p
        alpha = dot(r, z) / dot(p, Ap)
        x = x + alpha * p
        r_new = r - alpha * Ap
        z_new = multigrid_v_cycle(r_new)
        beta = dot(r_new, z_new) / dot(r, z)
        p = z_new + beta * p
        r = r_new
        z = z_new
```

**优点**：
- 结合多重网格的最优复杂度
- Krylov方法的鲁棒性
- 可处理非对称系统（BiCGSTAB+AMG）

### 8.8.4 非线性多重网格（FAS）

**Full Approximation Scheme**：
处理非线性问题 $N(u) = f$

```
FAS(level, u, f):
  if coarsest level:
    solve N(u) = f exactly
  else:
    # 前光滑
    smooth(N, u, f)
    
    # 限制
    u_coarse = restrict(u)
    r = f - N(u)
    f_coarse = restrict(r) + N_coarse(u_coarse)
    
    # 粗网格求解
    FAS(level+1, u_coarse, f_coarse)
    
    # 延拓修正
    u = u + prolongate(u_coarse - restrict(u))
    
    # 后光滑
    smooth(N, u, f)
```

### 8.8.5 时空多重网格

**并行时间积分**：
- 时间方向也进行多重网格
- Parareal算法的推广
- 适合长时间积分

**实现考虑**：
- 时间粗化策略
- 时空耦合的光滑器
- 并行效率分析

### 8.8.6 性能优化技巧

**缓存优化**：
```python
# 循环分块
@ti.kernel
def blocked_smooth(u: ti.field, f: ti.field):
    for bi, bj in ti.ndrange(N//B, N//B):
        for i, j in ti.ndrange(B, B):
            ii, jj = bi*B + i, bj*B + j
            if 0 < ii < N-1 and 0 < jj < N-1:
                u[ii,jj] = compute_update(ii, jj)
```

**向量化**：
- 使用SIMD指令
- 数据对齐
- 避免分支

**GPU特定优化**：
- Warp级别原语
- Texture内存用于插值
- 常量内存存储模板系数

---

## 本章小结

多重网格方法是求解大规模线性系统的最优算法之一，其核心思想是利用不同尺度网格的互补性：
- 细网格上的光滑迭代快速消除高频误差
- 粗网格将低频误差转换为相对高频并高效消除
- 递归应用形成O(N)复杂度的求解器

**关键概念回顾**：
1. **误差的频率分析**：理解光滑性质是多重网格成功的基础
2. **网格传输算子**：限制和延拓需要保持问题的物理性质
3. **循环策略**：V、W、FMG各有适用场景
4. **并行化**：粗网格是并行瓶颈，需要特殊处理
5. **代数多重网格**：只需矩阵信息，适用范围更广

**算法复杂度总结**：
- 计算复杂度：O(N)（2D和3D的V循环）
- 存储复杂度：O(N)
- 并行可扩展性：优秀（需要处理粗网格）

---

## 练习题

### 基础题

**练习8.1**：一维泊松方程的双网格分析
考虑一维泊松方程 $-u_{xx} = f$ 在 $[0,1]$ 上，使用中心差分离散。设细网格有 $n=7$ 个内部节点。
1. 写出细网格上的系数矩阵 $A^h$
2. 使用全权重限制和线性插值，构造粗网格（$n=3$）上的系数矩阵 $A^{2h}$
3. 计算并比较两个矩阵的特征值
4. 验证高频模式在细网格上的衰减率

<details>
<summary>提示</summary>

- 细网格矩阵是三对角矩阵，对角元素为 $2/h^2$，非对角元素为 $-1/h^2$
- 使用Galerkin条件 $A^{2h} = I_{2h}^h A^h I_h^{2h}$ 构造粗网格矩阵
- 三对角矩阵的特征值有解析表达式
</details>

<details>
<summary>答案</summary>

1. 细网格矩阵（$h=1/8$）：
$$A^h = \frac{64}{1} \begin{bmatrix}
2 & -1 & & & \\
-1 & 2 & -1 & & \\
& -1 & 2 & -1 & \\
& & \ddots & \ddots & \ddots \\
& & & -1 & 2
\end{bmatrix}_{7 \times 7}$$

2. 限制算子和延拓算子构造粗网格矩阵，得到：
$$A^{2h} = \frac{16}{1} \begin{bmatrix}
2 & -1 & 0 \\
-1 & 2 & -1 \\
0 & -1 & 2
\end{bmatrix}$$

3. 特征值：
- 细网格：$\lambda_k^h = 4\sin^2(k\pi/16)/h^2$，$k=1,...,7$
- 粗网格：$\lambda_k^{2h} = 4\sin^2(k\pi/8)/(2h)^2$，$k=1,2,3$

4. Jacobi迭代对模式$k=4,5,6,7$（高频）的衰减因子小于0.5
</details>

**练习8.2**：红黑Gauss-Seidel的并行实现
实现二维泊松方程的红黑Gauss-Seidel光滑器，要求：
1. 正确处理边界条件
2. 实现并行更新
3. 测试不同问题规模的加速比

<details>
<summary>提示</summary>

- 棋盘着色：$(i+j)\%2$ 决定红黑
- 红点和黑点可以分别并行更新
- 注意边界点的特殊处理
</details>

**练习8.3**：V循环收敛性分析
对于二维泊松方程，实现V循环并：
1. 测量不同光滑次数（$\nu_1, \nu_2$）的收敛因子
2. 绘制残差下降曲线
3. 比较与理论预测的差异

<details>
<summary>提示</summary>

- 收敛因子 = $\|r^{(k+1)}\| / \|r^{(k)}\|$
- 使用随机初值避免特殊情况
- 理论预测：$\rho \approx 0.1$ 对于 $\nu_1=\nu_2=1$
</details>

### 挑战题

**练习8.4**：各向异性问题的线光滑器
考虑各向异性泊松方程：
$$-\epsilon u_{xx} - u_{yy} = f$$
其中 $\epsilon = 0.001$。

1. 解释为什么标准点光滑器效果差
2. 实现x-line和y-line光滑器
3. 设计交替方向策略
4. 比较不同光滑器的效率

<details>
<summary>提示</summary>

- 各向异性导致x方向和y方向的耦合强度不同
- y-line光滑器对这个问题更有效
- 可以结合：先y-line再x-line
</details>

**练习8.5**：AMG的强连接分析
给定稀疏矩阵，实现AMG的C/F分裂算法：
1. 计算强连接矩阵（阈值$\theta=0.25$）
2. 实现经典C/F分裂
3. 构造插值算子
4. 验证粗网格的质量

<details>
<summary>提示</summary>

- 强连接：$|a_{ij}| \geq \theta \max_{k \neq i}|a_{ik}|$
- C点选择：最大化强影响
- 插值权重需要归一化
</details>

**练习8.6**：多重网格用于特征值问题
使用多重网格加速求解最小特征值：
$$Au = \lambda u$$

1. 推导Rayleigh商迭代的多重网格加速
2. 实现算法
3. 与幂法比较收敛速度

<details>
<summary>提示</summary>

- Rayleigh商：$\lambda = u^TAu / u^Tu$
- 每步需要求解 $(A-\sigma I)v = u$
- 多重网格作为内层求解器
</details>

**练习8.7**：非线性多重网格（FAS）
实现FAS求解非线性方程：
$$-\nabla^2 u + u^3 = f$$

1. 推导FAS的限制和延拓
2. 实现非线性光滑器
3. 与Newton-多重网格比较
4. 分析计算成本

<details>
<summary>提示</summary>

- 非线性残差：$r = f - (-\nabla^2 u + u^3)$
- 粗网格方程包含细网格解的信息
- 光滑器可用阻尼Newton
</details>

**练习8.8**：时空多重网格
对热方程实现时空多重网格：
$$u_t - \nabla^2 u = f$$

1. 设计时空网格的粗化策略
2. 推导时空耦合的光滑器
3. 分析并行效率
4. 与时间步进方法比较

<details>
<summary>提示</summary>

- 时间也可以看作一个维度
- 粗化可以在时间、空间或两者
- 注意因果性约束
</details>

---

## 常见陷阱与调试技巧

### 陷阱1：不当的边界条件处理
**问题**：边界条件在不同网格层不一致
**解决**：
- Dirichlet边界：细网格值直接复制到粗网格
- Neumann边界：使用修正的限制算子
- 混合边界：分别处理不同类型

### 陷阱2：错误的残差计算
**问题**：忘记更新残差或使用旧值
**调试**：
```python
# 正确的残差计算
r = f - A @ u  # 不要重用旧的r！
```

### 陷阱3：光滑不足
**症状**：V循环不收敛或收敛慢
**诊断**：
- 检查光滑后的残差频谱
- 增加光滑次数
- 尝试不同光滑器

### 陷阱4：粗网格过小
**问题**：直接求解器在"粗"网格上仍然昂贵
**解决**：
- 设置合理的粗网格大小（如 $4 \times 4$）
- 使用迭代法替代直接法
- 考虑聚合到单处理器

### 陷阱5：并行效率下降
**症状**：增加处理器数反而变慢
**分析**：
- Profile通信时间
- 检查负载均衡
- 优化粗网格策略

---

## 最佳实践检查清单

### 算法设计
- [ ] 根据问题特性选择几何或代数多重网格
- [ ] 光滑器与问题匹配（各向异性→线光滑）
- [ ] 合理的网格层数（通常5-10层）
- [ ] 适当的粗网格规模（避免过小或过大）

### 实现优化
- [ ] 数据结构支持高效的网格传输
- [ ] 避免不必要的内存分配
- [ ] 利用问题的对称性和稀疏性
- [ ] 预计算不变的系数

### 并行化
- [ ] 选择合适的数据分布策略
- [ ] 处理粗网格的并行瓶颈
- [ ] 最小化通信（宽halo、聚合等）
- [ ] 负载均衡（特别是自适应网格）

### 鲁棒性
- [ ] 处理奇异和近奇异系统
- [ ] 边界条件的一致处理
- [ ] 数值稳定性（避免除零等）
- [ ] 收敛性监测和自适应策略

### 调试验证
- [ ] 单元测试各组件（光滑器、传输等）
- [ ] 验证Galerkin条件
- [ ] 检查网格间的守恒性质
- [ ] 与已知解或其他方法对比
