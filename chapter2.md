# 第二章：拉格朗日视角（1）

本章将介绍拉格朗日视角下的物理仿真方法，从最基础的弹簧质点系统开始，逐步深入到现代流体仿真技术。我们将学习如何追踪材料粒子的运动轨迹，理解显式与隐式时间积分的权衡，掌握SPH和PBF等无网格方法，并实现高效的邻居搜索算法。通过本章学习，读者将建立起粒子系统仿真的完整知识体系，为后续的有限元和混合方法打下坚实基础。

## 2.1 弹簧质点系统

弹簧质点系统是物理仿真的基石，通过离散的质点和连接它们的弹簧来近似连续体的行为。这种方法直观简单，却能模拟出丰富的物理现象。

### 2.1.1 胡克定律与牛顿第二定律

弹簧的弹性力遵循胡克定律。对于连接质点$i$和$j$的弹簧，弹性力为：

$$\mathbf{f}_{ij} = -k(\|\mathbf{x}_i - \mathbf{x}_j\| - l_{ij})\frac{\mathbf{x}_i - \mathbf{x}_j}{\|\mathbf{x}_i - \mathbf{x}_j\|}$$

其中：
- $k$是弹簧刚度系数
- $l_{ij}$是弹簧的静止长度
- $\mathbf{x}_i, \mathbf{x}_j$是质点的位置向量

注意这个公式中的方向：当弹簧被拉伸时（$\|\mathbf{x}_i - \mathbf{x}_j\| > l_{ij}$），力的方向与$(\mathbf{x}_i - \mathbf{x}_j)$相反，将质点拉近；当弹簧被压缩时，力将质点推开。

每个质点的运动遵循牛顿第二定律：

$$m_i\frac{d\mathbf{v}_i}{dt} = \sum_j \mathbf{f}_{ij} + \mathbf{f}_i^{ext}$$

其中$\mathbf{f}_i^{ext}$包含重力、风力等外力。在实际计算中，我们通常会加入阻尼力来稳定系统：

$$\mathbf{f}_{ij}^{damping} = -c(\mathbf{v}_i - \mathbf{v}_j) \cdot \frac{\mathbf{x}_i - \mathbf{x}_j}{\|\mathbf{x}_i - \mathbf{x}_j\|} \cdot \frac{\mathbf{x}_i - \mathbf{x}_j}{\|\mathbf{x}_i - \mathbf{x}_j\|}$$

这是投影到弹簧方向的相对速度阻尼，$c$是阻尼系数。

### 2.1.2 弹簧刚度与阻尼系数

选择合适的弹簧参数对仿真稳定性至关重要。刚度系数$k$决定了系统的响应速度和数值稳定性条件。对于单个弹簧-质点系统，固有频率为：

$$\omega = \sqrt{\frac{k}{m}}$$

时间步长必须满足：

$$\Delta t < \frac{2}{\omega} = 2\sqrt{\frac{m}{k}}$$

阻尼系数的选择通常基于阻尼比$\zeta$：

$$c = 2\zeta\sqrt{km}$$

其中：
- $\zeta < 1$：欠阻尼，系统会振荡
- $\zeta = 1$：临界阻尼，最快达到平衡无振荡
- $\zeta > 1$：过阻尼，缓慢达到平衡

在布料仿真中，通常选择$\zeta \approx 0.1$保持一定的动态效果；在需要快速收敛的场合（如求解静态平衡），可以使用$\zeta \geq 1$。

### 2.1.3 系统的能量守恒

理想的弹簧质点系统应该保持总能量守恒。系统的总能量包括：

动能：
$$T = \sum_i \frac{1}{2}m_i\|\mathbf{v}_i\|^2$$

弹性势能：
$$U = \sum_{(i,j)} \frac{1}{2}k(\|\mathbf{x}_i - \mathbf{x}_j\| - l_{ij})^2$$

重力势能：
$$V = \sum_i m_i g h_i$$

理想情况下，$E = T + U + V$应该保持恒定。然而，数值积分会引入能量漂移：
- 显式积分通常会增加能量（数值不稳定）
- 隐式积分通常会耗散能量（数值阻尼）

辛积分器（如辛欧拉法）能更好地保持长时间仿真的能量守恒性。

### 2.1.4 多质点系统的矩阵表示

对于包含$n$个质点的系统，我们可以将所有方程组装成矩阵形式。定义：
- 位置向量：$\mathbf{x} = [\mathbf{x}_1^T, \mathbf{x}_2^T, ..., \mathbf{x}_n^T]^T \in \mathbb{R}^{3n}$
- 速度向量：$\mathbf{v} = [\mathbf{v}_1^T, \mathbf{v}_2^T, ..., \mathbf{v}_n^T]^T \in \mathbb{R}^{3n}$

系统的运动方程可以写成：

$$\mathbf{M}\frac{d\mathbf{v}}{dt} = -\mathbf{K}(\mathbf{x} - \mathbf{x}_0) - \mathbf{C}\mathbf{v} + \mathbf{f}^{ext}$$

其中：
- $\mathbf{M} = \text{diag}(m_1\mathbf{I}_3, m_2\mathbf{I}_3, ..., m_n\mathbf{I}_3)$是质量矩阵
- $\mathbf{K}$是刚度矩阵（注意：对于非线性弹簧，这是线性化的结果）
- $\mathbf{C}$是阻尼矩阵
- $\mathbf{x}_0$是静止位置

这种矩阵表示为后续的隐式积分方法奠定了基础。在Taichi中，我们通常使用`ti.Matrix`来高效地处理这些矩阵运算。

## 2.2 布料模拟

布料是弹簧质点系统的经典应用。通过在规则网格上布置质点，并用不同类型的弹簧连接，我们可以模拟出逼真的布料行为。

### 2.2.1 结构弹簧、剪切弹簧与弯曲弹簧

布料模型通常包含三种类型的弹簧：

**结构弹簧（Structural Springs）**：连接水平和垂直相邻的质点，维持布料的基本形状。对于网格位置$(i,j)$的质点，结构弹簧连接到$(i+1,j)$和$(i,j+1)$。

**剪切弹簧（Shearing Springs）**：连接对角相邻的质点，抵抗剪切变形。连接$(i,j)$到$(i+1,j+1)$和$(i+1,j-1)$。

**弯曲弹簧（Bending Springs）**：连接间隔一个质点的邻居，提供弯曲刚度。连接$(i,j)$到$(i+2,j)$和$(i,j+2)$。

不同弹簧的刚度系数反映了材料特性：
- 结构弹簧刚度$k_s$：决定拉伸强度
- 剪切弹簧刚度$k_{sh}$：通常$k_{sh} \approx 0.5k_s$
- 弯曲弹簧刚度$k_b$：通常$k_b \ll k_s$，使布料易弯曲

静止长度的设置：
- 结构弹簧：网格间距$l_0$
- 剪切弹簧：$\sqrt{2}l_0$
- 弯曲弹簧：$2l_0$

### 2.2.2 外力

布料仿真中常见的外力包括：

**重力**：
$$\mathbf{f}_i^{gravity} = m_i\mathbf{g}$$

其中$\mathbf{g} = [0, -9.8, 0]^T$ m/s²。

**风力**：风力的计算需要考虑布料的法向量和相对速度：
$$\mathbf{f}_{wind} = c_{wind}((\mathbf{v}_{wind} - \mathbf{v}) \cdot \mathbf{n})\mathbf{n}$$

其中：
- $\mathbf{v}_{wind}$是风速
- $\mathbf{v}$是布料表面速度
- $\mathbf{n}$是表面法向量
- $c_{wind}$是风力系数

对于三角形$(i,j,k)$，法向量计算为：
$$\mathbf{n} = \frac{(\mathbf{x}_j - \mathbf{x}_i) \times (\mathbf{x}_k - \mathbf{x}_i)}{\|(\mathbf{x}_j - \mathbf{x}_i) \times (\mathbf{x}_k - \mathbf{x}_i)\|}$$

**碰撞力**：当质点穿透障碍物时，施加排斥力：
$$\mathbf{f}_{collision} = k_{collision} \cdot d \cdot \mathbf{n}_{surface}$$

其中$d$是穿透深度，$\mathbf{n}_{surface}$是障碍物表面法向量。

### 2.2.3 约束处理

布料仿真经常需要处理各种约束：

**固定约束**：某些质点被固定在空间中（如晾衣绳上的衣服）。实现方法：
- 每个时间步后将固定点的位置和速度重置
- 或在求解时将对应自由度从系统中移除

**滑动约束**：质点被约束在曲线或曲面上滑动。处理方法：
1. 计算无约束的新位置$\mathbf{x}_{new}$
2. 将位置投影到约束流形上：$\mathbf{x}_{projected} = \text{Project}(\mathbf{x}_{new})$
3. 更新速度以保持一致性：$\mathbf{v}_{new} = (\mathbf{x}_{projected} - \mathbf{x}_{old})/\Delta t$

**距离约束**：保持两点间距离恒定（不可延展布料）。使用位置修正：
$$\Delta\mathbf{x}_i = -\frac{w_i}{w_i + w_j}\frac{\|\mathbf{x}_i - \mathbf{x}_j\| - l_{ij}}{\|\mathbf{x}_i - \mathbf{x}_j\|}(\mathbf{x}_i - \mathbf{x}_j)$$

其中$w_i = 1/m_i$是质点的逆质量。

### 2.2.4 自碰撞检测与响应

布料的自碰撞是仿真中的主要挑战之一。

**空间哈希加速**：将空间划分为网格，每个单元大小约为典型碰撞距离的2倍。对每个质点：
1. 计算其所在的网格单元：`cell = floor(position / cell_size)`
2. 检查相邻的27个单元中的所有质点
3. 只对距离小于阈值的质点对进行精确碰撞检测

**连续碰撞检测（CCD）**：对于高速运动，需要检测运动轨迹间的碰撞。给定两个三角形在$t$和$t+\Delta t$时刻的位置，求解：
$$\text{min } t \in [0,1] \text{ such that triangles intersect}$$

这通常需要求解一个三次方程。

**冲量响应**：当检测到碰撞时，施加冲量分离质点：
$$\mathbf{J} = \frac{2m_i m_j}{m_i + m_j}v_{rel} \cdot \mathbf{n}$$

其中$v_{rel}$是相对速度，$\mathbf{n}$是碰撞法向量。

为了稳定性，通常会加入一个小的分离距离$\epsilon$（如2mm），在质点距离小于$\epsilon$时就开始施加排斥力。

## 2.3 显式与隐式时间积分器

时间积分是将连续的运动方程离散化的关键步骤。不同的积分方法在精度、稳定性和计算效率上有不同的权衡。

### 2.3.1 前向欧拉法

前向欧拉法是最简单的显式积分方法：

$$\mathbf{v}_{t+\Delta t} = \mathbf{v}_t + \Delta t \frac{\mathbf{f}_t}{m}$$
$$\mathbf{x}_{t+\Delta t} = \mathbf{x}_t + \Delta t \mathbf{v}_t$$

注意速度更新使用的是时刻$t$的力，位置更新使用的是时刻$t$的速度。

**稳定性条件**：对于弹簧系统，前向欧拉法的稳定性要求：
$$\Delta t < \frac{2}{\omega_{max}} = 2\sqrt{\frac{m}{k_{max}}}$$

其中$\omega_{max}$是系统的最大固有频率。这意味着：
- 刚度越大，需要的时间步长越小
- 质量越小，需要的时间步长越小

**能量行为**：前向欧拉法会人工增加系统能量，导致"爆炸"不稳定性。考虑简谐振子$\ddot{x} + \omega^2 x = 0$，前向欧拉的放大因子为：
$$|1 + i\omega\Delta t| = \sqrt{1 + (\omega\Delta t)^2} > 1$$

### 2.3.2 辛欧拉法

辛欧拉法（Symplectic Euler）只需要交换更新顺序：

$$\mathbf{v}_{t+\Delta t} = \mathbf{v}_t + \Delta t \frac{\mathbf{f}_t}{m}$$
$$\mathbf{x}_{t+\Delta t} = \mathbf{x}_t + \Delta t \mathbf{v}_{t+\Delta t}$$

注意位置更新使用的是新的速度$\mathbf{v}_{t+\Delta t}$。

**辛性质**：辛欧拉法保持相空间体积（Liouville定理），这对长时间仿真很重要。对于哈密顿系统，辛积分器能更好地保持能量有界。

**稳定域**：辛欧拉法的稳定条件与前向欧拉相同，但能量行为更好——它不会单调增加或减少能量，而是在真实值附近振荡。

### 2.3.3 后向欧拉法

后向欧拉法是最基本的隐式方法：

$$\mathbf{v}_{t+\Delta t} = \mathbf{v}_t + \Delta t \frac{\mathbf{f}(\mathbf{x}_{t+\Delta t})}{m}$$
$$\mathbf{x}_{t+\Delta t} = \mathbf{x}_t + \Delta t \mathbf{v}_{t+\Delta t}$$

关键区别是力$\mathbf{f}$在新位置$\mathbf{x}_{t+\Delta t}$处计算，这导致需要求解非线性方程。

**无条件稳定性**：后向欧拉法对任意大的时间步长都是稳定的（A-稳定）。这允许使用大时间步长，特别适合刚性系统。

**数值阻尼**：后向欧拉会引入人工阻尼，导致能量耗散。振幅衰减因子约为：
$$\text{damping} \approx 1 - \frac{\omega^2\Delta t^2}{2}$$

**求解策略**：通常使用Newton-Raphson迭代或固定点迭代求解非线性系统。对于弹簧系统，一次Newton迭代通常就足够了。

### 2.3.4 中点法与RK4

**中点法（Midpoint Method）**：二阶精度的隐式方法：
$$\mathbf{k}_1 = \mathbf{v}_t + \frac{\Delta t}{2}\frac{\mathbf{f}(\mathbf{x}_t + \frac{\Delta t}{2}\mathbf{k}_1)}{m}$$
$$\mathbf{x}_{t+\Delta t} = \mathbf{x}_t + \Delta t \mathbf{k}_1$$

**通用$\beta$方法**：
$$\mathbf{x}_{t+\Delta t} = \mathbf{x}_t + \Delta t \mathbf{v}_t + \Delta t^2[(1-\beta)\mathbf{a}_t + \beta\mathbf{a}_{t+\Delta t}]$$

其中：
- $\beta = 0$：显式欧拉
- $\beta = 0.5$：中点法（二阶精度）
- $\beta = 1$：后向欧拉

**Runge-Kutta 4（RK4）**：四阶精度的显式方法：
$$\mathbf{k}_1 = \mathbf{f}(\mathbf{x}_t, t)$$
$$\mathbf{k}_2 = \mathbf{f}(\mathbf{x}_t + \frac{\Delta t}{2}\mathbf{k}_1, t + \frac{\Delta t}{2})$$
$$\mathbf{k}_3 = \mathbf{f}(\mathbf{x}_t + \frac{\Delta t}{2}\mathbf{k}_2, t + \frac{\Delta t}{2})$$
$$\mathbf{k}_4 = \mathbf{f}(\mathbf{x}_t + \Delta t\mathbf{k}_3, t + \Delta t)$$
$$\mathbf{x}_{t+\Delta t} = \mathbf{x}_t + \frac{\Delta t}{6}(\mathbf{k}_1 + 2\mathbf{k}_2 + 2\mathbf{k}_3 + \mathbf{k}_4)$$

RK4提供了精度和计算成本的良好平衡，每步需要4次力计算。

## 2.4 隐式积分的线性化与求解

隐式方法的核心挑战是求解非线性系统。本节介绍实用的线性化和迭代求解技术。

### 2.4.1 牛顿法线性化

对于隐式欧拉法，我们需要求解：
$$\mathbf{v}_{t+\Delta t} = \mathbf{v}_t + \Delta t \mathbf{M}^{-1}\mathbf{f}(\mathbf{x}_t + \Delta t \mathbf{v}_{t+\Delta t})$$

使用泰勒展开线性化力：
$$\mathbf{f}(\mathbf{x} + \Delta \mathbf{x}) \approx \mathbf{f}(\mathbf{x}) + \frac{\partial \mathbf{f}}{\partial \mathbf{x}}\Delta \mathbf{x}$$

定义雅可比矩阵$\mathbf{J} = \frac{\partial \mathbf{f}}{\partial \mathbf{x}}$，代入得到线性系统：
$$(\mathbf{I} - \Delta t^2 \mathbf{M}^{-1}\mathbf{J})\Delta \mathbf{v} = \Delta t \mathbf{M}^{-1}\mathbf{f}(\mathbf{x}_t)$$

对于弹簧力，雅可比矩阵的元素为：
$$\mathbf{J}_{ij} = -k\left[\frac{(\mathbf{x}_i - \mathbf{x}_j)(\mathbf{x}_i - \mathbf{x}_j)^T}{\|\mathbf{x}_i - \mathbf{x}_j\|^2} + \left(1 - \frac{l_{ij}}{\|\mathbf{x}_i - \mathbf{x}_j\|}\right)\mathbf{I}\right]$$

### 2.4.2 雅可比矩阵的构建

系统的完整雅可比矩阵是稀疏的，只在有弹簧连接的质点间有非零元素。

**矩阵结构**：对于$n$个质点的3D系统，雅可比矩阵大小为$3n \times 3n$。每个弹簧贡献4个$3 \times 3$的块：
- $\mathbf{J}_{ii}$：质点$i$对自身的导数
- $\mathbf{J}_{ij}$：质点$i$对质点$j$的导数
- $\mathbf{J}_{ji}$：质点$j$对质点$i$的导数
- $\mathbf{J}_{jj}$：质点$j$对自身的导数

**对称性**：由牛顿第三定律，$\mathbf{J}_{ij} = -\mathbf{J}_{ji}$，因此整体雅可比矩阵是对称的。

**正定性**：加入阻尼项和隐式积分的系数后，系统矩阵$\mathbf{A} = \mathbf{I} - \Delta t^2 \mathbf{M}^{-1}\mathbf{J}$通常是对称正定的，这保证了唯一解的存在。

### 2.4.3 Jacobi与Gauss-Seidel迭代

对于大规模系统，直接求解线性系统可能过于昂贵。迭代方法提供了实用的替代方案。

**Jacobi迭代**：
$$\mathbf{x}_i^{(k+1)} = \frac{1}{\mathbf{A}_{ii}}\left(\mathbf{b}_i - \sum_{j \neq i} \mathbf{A}_{ij}\mathbf{x}_j^{(k)}\right)$$

Jacobi方法的优点：
- 完全并行，适合GPU实现
- 实现简单，内存访问规则
- 对于对角占优矩阵收敛

**Gauss-Seidel迭代**：
$$\mathbf{x}_i^{(k+1)} = \frac{1}{\mathbf{A}_{ii}}\left(\mathbf{b}_i - \sum_{j < i} \mathbf{A}_{ij}\mathbf{x}_j^{(k+1)} - \sum_{j > i} \mathbf{A}_{ij}\mathbf{x}_j^{(k)}\right)$$

Gauss-Seidel使用最新的值，收敛速度约为Jacobi的两倍。

**松弛因子**：超松弛（SOR）可以加速收敛：
$$\mathbf{x}_i^{(k+1)} = (1-\omega)\mathbf{x}_i^{(k)} + \omega\mathbf{x}_i^{GS}$$

其中$\omega \in (0, 2)$是松弛因子。最优$\omega$依赖于矩阵谱半径，实践中通常选择$\omega \in [1.2, 1.8]$。

**收敛准则**：迭代直到残差足够小：
$$\|\mathbf{A}\mathbf{x}^{(k)} - \mathbf{b}\| < \epsilon\|\mathbf{b}\|$$

通常$\epsilon = 10^{-3}$到$10^{-5}$就足够了。

### 2.4.4 共轭梯度法简介

共轭梯度（CG）法是求解对称正定系统的首选方法。

**基本算法**：
1. 初始化：$\mathbf{r}_0 = \mathbf{b} - \mathbf{A}\mathbf{x}_0$，$\mathbf{p}_0 = \mathbf{r}_0$
2. 迭代：
   - $\alpha_k = \frac{\mathbf{r}_k^T\mathbf{r}_k}{\mathbf{p}_k^T\mathbf{A}\mathbf{p}_k}$
   - $\mathbf{x}_{k+1} = \mathbf{x}_k + \alpha_k\mathbf{p}_k$
   - $\mathbf{r}_{k+1} = \mathbf{r}_k - \alpha_k\mathbf{A}\mathbf{p}_k$
   - $\beta_k = \frac{\mathbf{r}_{k+1}^T\mathbf{r}_{k+1}}{\mathbf{r}_k^T\mathbf{r}_k}$
   - $\mathbf{p}_{k+1} = \mathbf{r}_{k+1} + \beta_k\mathbf{p}_k$

**收敛性**：CG理论上在$n$步内收敛（$n$是矩阵维度），但实际收敛速度取决于条件数$\kappa(\mathbf{A})$：
$$\|\mathbf{e}_k\|_{\mathbf{A}} \leq 2\left(\frac{\sqrt{\kappa}-1}{\sqrt{\kappa}+1}\right)^k\|\mathbf{e}_0\|_{\mathbf{A}}$$

**预条件**：使用预条件矩阵$\mathbf{M} \approx \mathbf{A}$可以显著加速收敛。对于弹簧系统，对角预条件（Jacobi预条件）通常就很有效：
$$\mathbf{M} = \text{diag}(\mathbf{A})$$

在Taichi中实现CG时，注意利用稀疏矩阵结构和向量化操作来优化性能。

## 2.5 光滑粒子流体动力学(SPH)

SPH是一种无网格的拉格朗日方法，通过粒子和核函数来离散化连续场。它特别适合处理大变形和自由表面流动。

### 2.5.1 核函数与粒子近似

SPH的核心思想是用离散粒子的加权和来近似连续场。对于任意物理量$A(\mathbf{x})$：

$$A(\mathbf{x}) = \int A(\mathbf{x}')W(\mathbf{x} - \mathbf{x}', h)d\mathbf{x}' \approx \sum_j A_j \frac{m_j}{\rho_j}W(\mathbf{x} - \mathbf{x}_j, h)$$

其中：
- $W(\mathbf{r}, h)$是核函数，$h$是光滑长度
- $m_j$是粒子$j$的质量
- $\rho_j$是粒子$j$处的密度

**核函数性质**：
1. 归一化：$\int W(\mathbf{r}, h)d\mathbf{r} = 1$
2. 紧支性：$W(\mathbf{r}, h) = 0$ for $\|\mathbf{r}\| > kh$（通常$k=2$或$3$）
3. 对称性：$W(\mathbf{r}, h) = W(-\mathbf{r}, h)$
4. 光滑性：至少$C^1$连续

**Cubic Spline核函数**（最常用）：
$$W(r, h) = \frac{\sigma_d}{h^d} \begin{cases}
1 - \frac{3}{2}q^2 + \frac{3}{4}q^3 & 0 \leq q \leq 1 \\
\frac{1}{4}(2-q)^3 & 1 < q \leq 2 \\
0 & q > 2
\end{cases}$$

其中$q = r/h$，$\sigma_d$是归一化常数（2D: $10/7\pi$，3D: $1/\pi$）。

**梯度计算**：SPH中梯度的近似特别重要：
$$\nabla A(\mathbf{x}) \approx \sum_j A_j \frac{m_j}{\rho_j}\nabla W(\mathbf{x} - \mathbf{x}_j, h)$$

注意梯度作用在核函数上，而不是物理量上。

### 2.5.2 密度计算与状态方程

**密度计算**是SPH的第一步：
$$\rho_i = \sum_j m_j W(\mathbf{x}_i - \mathbf{x}_j, h)$$

这保证了质量守恒：$\sum_i m_i = \sum_i \rho_i V_i$。

**状态方程**将密度与压力联系起来。对于弱可压缩流体，常用Tait方程：
$$p = B\left[\left(\frac{\rho}{\rho_0}\right)^\gamma - 1\right]$$

其中：
- $B = \frac{\rho_0 c_s^2}{\gamma}$，$c_s$是声速
- $\gamma = 7$（水）
- $\rho_0$是参考密度

选择合适的$B$值很重要：
- 太小：流体过于可压缩
- 太大：时间步长限制严格（$\Delta t \propto 1/c_s$）

实践中，通常选择$c_s = 10 v_{max}$，其中$v_{max}$是预期的最大速度。

### 2.5.3 压力梯度与粘性力

**压力梯度**的计算需要特别注意，直接使用梯度公式会导致动量不守恒。对称形式：
$$\nabla p_i = \rho_i \sum_j m_j \left(\frac{p_i}{\rho_i^2} + \frac{p_j}{\rho_j^2}\right)\nabla W_{ij}$$

这保证了作用力与反作用力相等。

**粘性力**模拟流体的内摩擦。物理粘性：
$$\mathbf{f}_i^{viscosity} = \mu \sum_j \frac{m_j}{\rho_j}\frac{\mathbf{v}_j - \mathbf{v}_i}{r_{ij}^2 + 0.01h^2}\mathbf{r}_{ij} \cdot \nabla W_{ij}$$

其中$\mu$是动力粘度，$0.01h^2$项防止除零。

**人工粘性**用于稳定数值解：
$$\Pi_{ij} = \begin{cases}
-\alpha \bar{c}_{ij} \frac{\mathbf{v}_{ij} \cdot \mathbf{r}_{ij}}{\bar{\rho}_{ij}(r_{ij}^2 + 0.01h^2)} & \mathbf{v}_{ij} \cdot \mathbf{r}_{ij} < 0 \\
0 & \text{otherwise}
\end{cases}$$

其中$\alpha \approx 0.01-0.1$，$\bar{c}_{ij}$和$\bar{\rho}_{ij}$是平均值。

### 2.5.4 表面张力模型

表面张力在小尺度流动中很重要。**CSF（Continuum Surface Force）方法**将表面张力转换为体积力。

首先定义颜色场：
$$c_i = \sum_j \frac{m_j}{\rho_j}W_{ij}$$

表面法向量：
$$\mathbf{n}_i = \nabla c_i = \sum_j \frac{m_j}{\rho_j}\nabla W_{ij}$$

曲率：
$$\kappa_i = -\nabla \cdot \hat{\mathbf{n}}_i = -\sum_j \frac{m_j}{\rho_j}\frac{\mathbf{n}_j - \mathbf{n}_i}{|\mathbf{n}_j - \mathbf{n}_i|} \cdot \nabla W_{ij}$$

表面张力：
$$\mathbf{f}_i^{surface} = \sigma \kappa_i \mathbf{n}_i$$

其中$\sigma$是表面张力系数。只在$|\mathbf{n}_i| > \epsilon$的粒子上施加此力，以识别表面粒子。

**算法流程总结**：
1. 邻居搜索
2. 计算密度：$\rho_i = \sum_j m_j W_{ij}$
3. 计算压力：$p_i = B[(\rho_i/\rho_0)^7 - 1]$
4. 计算加速度：
   - 压力：$-\sum_j m_j (p_i/\rho_i^2 + p_j/\rho_j^2)\nabla W_{ij}$
   - 粘性：粘性力公式
   - 重力：$\mathbf{g}$
   - 表面张力（如需要）
5. 时间积分更新位置和速度

## 2.6 基于位置的流体(PBF)

PBF将不可压缩性作为约束条件，通过迭代投影来满足。这种方法稳定性好，特别适合实时应用。

### 2.6.1 位置约束与拉格朗日乘子

PBF的核心是密度约束：
$$C_i(\mathbf{x}_1, ..., \mathbf{x}_n) = \frac{\rho_i}{\rho_0} - 1 = 0$$

使用拉格朗日乘子法，位置修正为：
$$\Delta \mathbf{x}_i = -\sum_j \lambda_j \nabla_{\mathbf{x}_i} C_j$$

其中$\lambda_j$是拉格朗日乘子。对于密度约束：
$$\nabla_{\mathbf{x}_i} C_j = \frac{1}{\rho_0}\nabla_{\mathbf{x}_i}\rho_j = \frac{1}{\rho_0}\sum_k m_k \nabla_{\mathbf{x}_i} W_{jk}$$

当$i = j$时：
$$\nabla_{\mathbf{x}_i} C_i = \frac{1}{\rho_0}\sum_{k \neq i} m_k \nabla W_{ik}$$

当$i \neq j$时：
$$\nabla_{\mathbf{x}_i} C_j = -\frac{m_i}{\rho_0}\nabla W_{ij}$$

### 2.6.2 不可压缩性约束

通过Newton迭代求解$\lambda_i$：
$$\lambda_i = -\frac{C_i(\mathbf{x})}{\sum_k |\nabla_{\mathbf{x}_k} C_i|^2 + \epsilon}$$

其中$\epsilon$是松弛参数，防止除零。

完整的分母展开：
$$\sum_k |\nabla_{\mathbf{x}_k} C_i|^2 = \frac{1}{\rho_0^2}\left[\left|\sum_{k \neq i} m_k \nabla W_{ik}\right|^2 + \sum_{k \neq i} m_k^2 |\nabla W_{ik}|^2\right]$$

位置修正：
$$\Delta \mathbf{x}_i = -\frac{1}{\rho_0}\sum_j \lambda_j m_j \nabla W_{ij}$$

注意这里使用了对称化的梯度，保证动量守恒。

### 2.6.3 人工压力与涡量补偿

**人工压力**防止粒子聚集。在计算$\lambda$时加入修正项：
$$s_{corr} = -k\left(\frac{W(\mathbf{x}_i - \mathbf{x}_j, h)}{W(\Delta \mathbf{q}, h)}\right)^n$$

其中$k = 0.1$，$n = 4$，$\Delta \mathbf{q} = 0.1h$。修正后的$\lambda$：
$$\lambda_i^{corr} = \lambda_i + s_{corr}$$

**涡量补偿**恢复数值耗散损失的旋转运动：
1. 计算涡量：$\boldsymbol{\omega}_i = \nabla \times \mathbf{v}_i = \sum_j \frac{m_j}{\rho_j}\mathbf{v}_{ij} \times \nabla W_{ij}$
2. 计算涡量力：$\mathbf{f}_i^{vorticity} = \epsilon(\mathbf{N}_i \times \boldsymbol{\omega}_i)$
3. 其中$\mathbf{N}_i = \nabla|\boldsymbol{\omega}|_i/|\nabla|\boldsymbol{\omega}|_i|$，$\epsilon \approx 0.01$

### 2.6.4 PCISPH与DFSPH方法

**PCISPH（Predictive-Corrective Incompressible SPH）**：
1. 预测速度和位置
2. 计算预测密度误差
3. 计算压力：$p_i = \delta \rho_i^{err} / (\Delta t^2 \beta)$
4. 其中$\beta = 2m^2(\sum_j \nabla W_{ij})^2/\rho_0^2$
5. 迭代直到密度误差小于阈值

**DFSPH（Divergence-Free SPH）**：
除了密度约束，还强制速度场散度为零：
$$\nabla \cdot \mathbf{v} = 0$$

这需要两个投影步骤：
1. 密度不变投影（确保$D\rho/Dt = 0$）
2. 散度自由投影（确保$\nabla \cdot \mathbf{v} = 0$）

DFSPH能更好地保持不可压缩性，减少密度振荡。

**PBF算法总结**：
```
1. 预测位置：x* = x + Δt·v + Δt²·f_ext
2. 邻居搜索
3. while 密度误差 > 阈值：
   a. 计算密度
   b. 计算λ_i
   c. 计算位置修正Δx
   d. 更新位置：x* += Δx
4. 更新速度：v = (x* - x)/Δt
5. 应用涡量补偿
6. 更新位置：x = x*
```

## 2.7 体素化：从三角网格生成粒子

将三角网格转换为粒子表示是初始化粒子系统的关键步骤。

### 2.7.1 点在多边形内测试

**射线法（Ray Casting）**：从测试点发射射线，计算与多边形边界的交点数：
- 奇数次相交：点在内部
- 偶数次相交：点在外部

实现细节：
1. 选择射线方向（通常是+X方向）
2. 处理特殊情况：
   - 射线穿过顶点：只在顶点Y坐标大于射线Y坐标时计数
   - 射线与边平行：忽略该边
   - 数值精度：使用$\epsilon$容差

**绕数（Winding Number）方法**：计算多边形相对于测试点的绕数：
$$w = \frac{1}{2\pi}\sum_{i=0}^{n-1} \theta_i$$

其中$\theta_i$是从点到边$(v_i, v_{i+1})$的有向角。$|w| \geq 1$表示点在内部。

**角度和方法**：对于凸多边形，计算点到所有顶点的角度和：
$$\sum_{i=0}^{n-1} \angle(v_i, p, v_{i+1}) = \begin{cases}
2\pi & \text{点在内部} \\
0 & \text{点在外部}
\end{cases}$$

### 2.7.2 射线-三角形相交

**Möller-Trumbore算法**是最高效的射线-三角形相交测试。给定射线$\mathbf{r}(t) = \mathbf{o} + t\mathbf{d}$和三角形顶点$\mathbf{v}_0, \mathbf{v}_1, \mathbf{v}_2$：

参数化三角形内的点：
$$\mathbf{p}(u,v) = (1-u-v)\mathbf{v}_0 + u\mathbf{v}_1 + v\mathbf{v}_2$$

求解相交：
$$\mathbf{o} + t\mathbf{d} = (1-u-v)\mathbf{v}_0 + u\mathbf{v}_1 + v\mathbf{v}_2$$

整理得到线性系统：
$$\begin{bmatrix}
-\mathbf{d} & \mathbf{v}_1-\mathbf{v}_0 & \mathbf{v}_2-\mathbf{v}_0
\end{bmatrix}
\begin{bmatrix}
t \\ u \\ v
\end{bmatrix} = \mathbf{o} - \mathbf{v}_0$$

使用Cramer法则求解：
```
E1 = v1 - v0
E2 = v2 - v0
P = d × E2
det = E1 · P
if |det| < ε: return no_intersection

T = o - v0
u = (T · P) / det
if u < 0 or u > 1: return no_intersection

Q = T × E1
v = (d · Q) / det
if v < 0 or u + v > 1: return no_intersection

t = (E2 · Q) / det
if t > 0: return intersection at t
```

### 2.7.3 有符号距离场(SDF)

SDF提供了更丰富的几何信息：
$$\phi(\mathbf{x}) = \begin{cases}
-d(\mathbf{x}, \partial\Omega) & \mathbf{x} \in \Omega \\
+d(\mathbf{x}, \partial\Omega) & \mathbf{x} \notin \Omega
\end{cases}$$

**快速行进法（FMM）构建SDF**：
1. 初始化：
   - 边界点：$\phi = 0$，标记为已知
   - 其他点：$\phi = \infty$，标记为远场
2. 维护一个优先队列（最小堆），包含所有试探点
3. 重复直到队列为空：
   - 取出最小距离的试探点，标记为已知
   - 更新其邻居的距离值
   - 将新的试探点加入队列

距离更新使用Eikonal方程的数值解：
$$|\nabla \phi| = 1$$

对于规则网格，使用迎风差分：
$$\max(D^{-x}\phi, -D^{+x}\phi, 0)^2 + \max(D^{-y}\phi, -D^{+y}\phi, 0)^2 = 1$$

### 2.7.4 自适应采样策略

均匀采样可能导致粒子分布不理想。自适应采样根据局部特征调整密度。

**基于曲率的采样**：
1. 计算表面曲率：$\kappa = \nabla \cdot (\nabla\phi/|\nabla\phi|)$
2. 采样密度：$\rho_{sample} = \rho_{min} + (\rho_{max} - \rho_{min})|\kappa|/\kappa_{max}$
3. 在高曲率区域放置更多粒子

**蓝噪声分布**产生视觉上更均匀的分布：
- 任意两个粒子的最小距离被最大化
- 避免聚集和空洞

**Poisson Disk采样**：
1. 初始化：随机选择第一个样本
2. 对每个现有样本：
   - 在环形区域$[r, 2r]$内生成$k$个候选点（通常$k=30$）
   - 检查每个候选点是否与现有点距离$> r$
   - 接受第一个满足条件的候选点
3. 使用空间数据结构（如网格）加速距离查询

**算法实现要点**：
```
function voxelize(mesh, particle_radius):
    # 1. 计算包围盒
    bbox = compute_bbox(mesh)
    
    # 2. 构建SDF（可选）
    sdf = build_sdf(mesh, bbox, grid_size=particle_radius/2)
    
    # 3. 生成候选点
    if use_adaptive_sampling:
        candidates = adaptive_sample(sdf, particle_radius)
    else:
        candidates = uniform_grid(bbox, particle_radius*2)
    
    # 4. 筛选内部点
    particles = []
    for p in candidates:
        if point_in_mesh(p, mesh):  # 使用射线法或SDF
            particles.append(p)
    
    # 5. 后处理（可选）
    particles = blue_noise_relaxation(particles)
    
    return particles
```

## 2.8 快速邻居搜索

高效的邻居搜索是粒子方法性能的关键。对于$n$个粒子，朴素的$O(n^2)$搜索是不可接受的。

### 2.8.1 均匀网格法

最简单有效的加速结构。将空间划分为规则网格，每个单元存储其内的粒子。

**网格大小选择**：设置为支持半径$h$，这样每个粒子只需检查$3^d$个相邻单元（$d$是维度）。

**实现步骤**：
1. 计算粒子所在单元：`cell = floor(position / h)`
2. 将粒子ID插入对应单元的列表
3. 查询时，遍历相邻单元的所有粒子

**数据结构选择**：
- 密集网格：直接分配所有单元，适合粒子分布均匀的情况
- 哈希表：只存储非空单元，适合稀疏分布

**并行构建**：
```
# 并行计数
for each particle i in parallel:
    cell = floor(x[i] / h)
    atomic_add(count[cell], 1)

# 前缀和计算偏移
offset[0] = 0
for i = 1 to num_cells:
    offset[i] = offset[i-1] + count[i-1]

# 并行插入
for each particle i in parallel:
    cell = floor(x[i] / h)
    idx = atomic_add(offset[cell], 1)
    particle_list[idx] = i
```

### 2.8.2 紧凑哈希(Compact Hashing)

对于稀疏分布，哈希表更节省内存。使用空间哈希函数将3D坐标映射到1D：

**Z-order哈希**（Morton码）：
```python
def morton_3d(x, y, z):
    x = spread_bits(x)  # 0b0000 -> 0b0000000
    y = spread_bits(y)  # 0b0000 -> 0b0000000
    z = spread_bits(z)  # 0b0000 -> 0b0000000
    return x | (y << 1) | (z << 2)

def spread_bits(v):
    v = (v | (v << 16)) & 0x030000FF
    v = (v | (v << 8))  & 0x0300F00F
    v = (v | (v << 4))  & 0x030C30C3
    v = (v | (v << 2))  & 0x09249249
    return v
```

**哈希冲突处理**：
- 开放寻址：线性探测或二次探测
- 链表法：每个槽存储链表
- Cuckoo hashing：保证最坏情况O(1)查询

### 2.8.3 Z-order与空间填充曲线

空间填充曲线将多维空间映射到一维，同时保持空间局部性。

**Z-order曲线优势**：
1. 保持空间局部性：相近的点在曲线上也相近（大部分情况）
2. 缓存友好：顺序访问时空间跳跃小
3. 易于计算：位操作即可

**其他空间填充曲线**：
- Hilbert曲线：局部性更好但计算复杂
- Peano曲线：类似Hilbert
- Gray码：减少位翻转

**分层Z-order**：
```
level_0: [    0    ] 整个空间
level_1: [0][1][2][3] 四个象限
level_2: 16个子区域
...
```

用于自适应结构，如八叉树。

### 2.8.4 并行邻居搜索算法

GPU上的高效实现需要特别考虑。

**原子操作构建**：
```cuda
__global__ void build_grid(float3* positions, int* grid, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    
    int3 cell = make_int3(positions[idx] / cell_size);
    int hash = hash_function(cell);
    
    // 原子操作插入
    int old = atomicExch(&grid[hash], idx);
    if (old != -1) {
        // 处理冲突：链表或其他方法
    }
}
```

**避免竞态条件**：
1. 两遍法：第一遍计数，第二遍填充
2. 原子操作：使用atomicAdd等
3. 排序法：先排序粒子，然后顺序构建

**查询优化**：
- 使用共享内存缓存邻居单元
- 合并内存访问
- 使用纹理内存（如果适用）

**性能考虑**：
- 负载均衡：每个线程处理相似数量的邻居
- 内存合并：连续线程访问连续内存
- 占用率：平衡寄存器使用和线程块大小

## 2.9 刚体模拟简介

刚体是不可变形的理想化物体，其内部任意两点的距离保持恒定。

### 2.9.1 刚体运动学

刚体的状态由质心位置和姿态完全确定。

**状态表示**：
- 位置：$\mathbf{x} \in \mathbb{R}^3$
- 姿态：四元数$\mathbf{q} = (w, x, y, z)$，满足$|\mathbf{q}| = 1$
- 线速度：$\mathbf{v} \in \mathbb{R}^3$
- 角速度：$\boldsymbol{\omega} \in \mathbb{R}^3$

**四元数优势**：
- 无奇异性（vs欧拉角）
- 紧凑表示（4个数vs9个旋转矩阵）
- 插值简单（SLERP）

**四元数运算**：
- 旋转向量：$\mathbf{v}' = \mathbf{q} \mathbf{v} \mathbf{q}^*$
- 四元数乘法：$\mathbf{q}_1 \mathbf{q}_2 = (w_1w_2 - \mathbf{v}_1 \cdot \mathbf{v}_2, w_1\mathbf{v}_2 + w_2\mathbf{v}_1 + \mathbf{v}_1 \times \mathbf{v}_2)$
- 角速度更新：$\dot{\mathbf{q}} = \frac{1}{2}\boldsymbol{\omega} \mathbf{q}$

**速度扭曲（Twist）表示**：
$$\boldsymbol{\xi} = \begin{bmatrix} \mathbf{v} \\ \boldsymbol{\omega} \end{bmatrix} \in \mathbb{R}^6$$

统一表示线速度和角速度，便于约束求解。

### 2.9.2 刚体动力学

**牛顿-欧拉方程**：
$$m\dot{\mathbf{v}} = \mathbf{F}$$
$$\mathbf{I}\dot{\boldsymbol{\omega}} + \boldsymbol{\omega} \times \mathbf{I}\boldsymbol{\omega} = \boldsymbol{\tau}$$

其中：
- $m$是质量
- $\mathbf{I}$是惯性张量（世界坐标系）
- $\mathbf{F}$是合外力
- $\boldsymbol{\tau}$是合外力矩

**惯性张量计算**：
对于离散质点系统：
$$\mathbf{I} = \sum_i m_i [(\mathbf{r}_i^T\mathbf{r}_i)\mathbf{I}_3 - \mathbf{r}_i\mathbf{r}_i^T]$$

对于连续体：
$$I_{ij} = \int_V \rho(x_k x_k \delta_{ij} - x_i x_j)dV$$

**世界坐标系惯性张量**：
$$\mathbf{I}_{world} = \mathbf{R}\mathbf{I}_{body}\mathbf{R}^T$$

其中$\mathbf{R}$是旋转矩阵，$\mathbf{I}_{body}$是物体坐标系下的惯性张量（常量）。

### 2.9.3 碰撞检测

**宽相位（Broad Phase）**：快速排除不可能碰撞的物体对。
- AABB（轴对齐包围盒）树
- Sweep and Prune
- 空间哈希

**窄相位（Narrow Phase）**：精确检测碰撞。

**GJK算法**（Gilbert-Johnson-Keerthi）：计算凸物体间的最近距离。
- 基于Minkowski差的原理
- 迭代搜索支撑点
- 可扩展为EPA（Expanding Polytope Algorithm）计算穿透深度

**SAT（分离轴定理）**：如果两个凸物体不相交，则存在一个轴使得两物体在该轴上的投影不重叠。
- 对于凸多面体，只需检查面法线和边叉积方向
- 对于盒子，只需检查15个轴（3+3+9）

### 2.9.4 碰撞响应

**冲量计算**：
碰撞点的相对速度：
$$\mathbf{v}_{rel} = (\mathbf{v}_A + \boldsymbol{\omega}_A \times \mathbf{r}_A) - (\mathbf{v}_B + \boldsymbol{\omega}_B \times \mathbf{r}_B)$$

法向冲量：
$$j_n = \frac{-(1+e)(\mathbf{v}_{rel} \cdot \mathbf{n})}{\frac{1}{m_A} + \frac{1}{m_B} + \mathbf{n} \cdot [(\mathbf{I}_A^{-1}(\mathbf{r}_A \times \mathbf{n})) \times \mathbf{r}_A + (\mathbf{I}_B^{-1}(\mathbf{r}_B \times \mathbf{n})) \times \mathbf{r}_B]}$$

其中$e$是恢复系数（0=完全非弹性，1=完全弹性）。

**摩擦处理**：
库仑摩擦模型：$|\mathbf{f}_t| \leq \mu|\mathbf{f}_n|$

切向冲量：
$$j_t = \min(\mu j_n, j_{t,stick})$$

其中$j_{t,stick}$是使相对切向速度为零所需的冲量。

**LCP（线性互补问题）求解**：
对于多接触点，需要同时求解所有约束：
$$\mathbf{A}\boldsymbol{\lambda} + \mathbf{b} \geq 0, \boldsymbol{\lambda} \geq 0, \boldsymbol{\lambda}^T(\mathbf{A}\boldsymbol{\lambda} + \mathbf{b}) = 0$$

使用投影Gauss-Seidel（PGS）迭代求解。

## 本章小结

本章介绍了拉格朗日视角下的粒子仿真方法，涵盖了从基础的弹簧质点系统到先进的流体仿真技术：

**核心概念**：
- **拉格朗日视角**：追踪材料粒子的运动轨迹，$\frac{D}{Dt} = \frac{\partial}{\partial t} + \mathbf{v} \cdot \nabla$
- **时间积分**：显式方法（前向欧拉、RK4）条件稳定但高效；隐式方法（后向欧拉）无条件稳定但需求解非线性系统
- **SPH核近似**：$A(\mathbf{x}) \approx \sum_j A_j \frac{m_j}{\rho_j}W(\mathbf{x} - \mathbf{x}_j, h)$
- **位置基方法**：将物理约束转化为几何约束，通过迭代投影满足

**关键公式**：
1. 弹簧力：$\mathbf{f}_{ij} = -k(\|\mathbf{x}_i - \mathbf{x}_j\| - l_{ij})\frac{\mathbf{x}_i - \mathbf{x}_j}{\|\mathbf{x}_i - \mathbf{x}_j\|}$
2. SPH密度：$\rho_i = \sum_j m_j W_{ij}$
3. SPH压力梯度：$\nabla p_i = \rho_i \sum_j m_j (\frac{p_i}{\rho_i^2} + \frac{p_j}{\rho_j^2})\nabla W_{ij}$
4. PBF密度约束：$C_i = \frac{\rho_i}{\rho_0} - 1 = 0$
5. 刚体角动量：$\mathbf{I}\dot{\boldsymbol{\omega}} + \boldsymbol{\omega} \times \mathbf{I}\boldsymbol{\omega} = \boldsymbol{\tau}$

**数值稳定性要点**：
- 显式积分时间步长限制：$\Delta t < 2\sqrt{m/k}$
- SPH中选择合适的声速：$c_s \approx 10 v_{max}$
- PBF中使用人工压力防止粒子聚集
- 邻居搜索使用空间数据结构将复杂度从$O(n^2)$降至$O(n)$

## 练习题

### 基础题

**习题2.1** 弹簧系统稳定性分析
考虑一个由两个质点组成的弹簧系统，质量分别为$m_1 = 1$kg和$m_2 = 2$kg，弹簧刚度$k = 100$N/m。
a) 使用前向欧拉法，计算保证数值稳定的最大时间步长。
b) 如果使用辛欧拉法，系统的总能量会如何变化？
c) 推导该系统的固有频率。

<details>
<summary>提示</summary>
考虑系统的最高频率模式，这决定了稳定性条件。对于两质点系统，可以先转换到质心坐标系。
</details>

<details>
<summary>答案</summary>

a) 转换到相对坐标：$\mu = \frac{m_1 m_2}{m_1 + m_2} = \frac{2}{3}$kg
   固有频率：$\omega = \sqrt{\frac{k}{\mu}} = \sqrt{\frac{100}{2/3}} = \sqrt{150} \approx 12.25$ rad/s
   最大时间步长：$\Delta t_{max} = \frac{2}{\omega} = \frac{2}{12.25} \approx 0.163$s

b) 辛欧拉法保持相空间体积，总能量在真实值附近振荡，振荡幅度约为$O(\Delta t^2)$。

c) 系统有两个模式：
   - 质心模式（频率为0）
   - 相对运动模式：$\omega = \sqrt{\frac{k(m_1 + m_2)}{m_1 m_2}} = \sqrt{150}$ rad/s
</details>

**习题2.2** SPH核函数性质
给定Cubic Spline核函数，证明：
a) 核函数满足归一化条件
b) 计算2D情况下的归一化常数$\sigma_2$
c) 推导核函数的梯度表达式

<details>
<summary>提示</summary>
使用极坐标进行积分，注意分段函数的处理。
</details>

<details>
<summary>答案</summary>

a) 2D极坐标积分：
   $\int_0^{2\pi} \int_0^{\infty} W(r,h) r dr d\theta = 2\pi \sigma_2 \int_0^{2h} W(r,h) r dr = 1$

b) 分段积分：
   $\int_0^h (1 - \frac{3r^2}{2h^2} + \frac{3r^3}{4h^3}) r dr + \int_h^{2h} \frac{1}{4}(\frac{2h-r}{h})^3 r dr = \frac{7h^2}{10\pi\sigma_2}$
   
   因此$\sigma_2 = \frac{10}{7\pi}$

c) 梯度：
   $\nabla W = \frac{\partial W}{\partial r} \frac{\mathbf{r}}{r}$，其中
   $\frac{\partial W}{\partial r} = \frac{\sigma_2}{h^2} \begin{cases}
   -\frac{3r}{h^2} + \frac{9r^2}{4h^3} & 0 \leq r \leq h \\
   -\frac{3(2h-r)^2}{4h^3} & h < r \leq 2h \\
   0 & r > 2h
   \end{cases}$
</details>

**习题2.3** PBF约束投影
在PBF算法中，考虑3个粒子的一维情况，位置分别为$x_1 = 0$, $x_2 = 0.8h$, $x_3 = 1.6h$，每个粒子质量为$m$。
a) 计算粒子2的密度（使用简化的线性核$W(r) = 1 - r/h$当$r < h$）
b) 计算密度约束的梯度
c) 求解拉格朗日乘子$\lambda_2$

<details>
<summary>提示</summary>
在一维情况下，梯度简化为导数。注意核函数的支持域。
</details>

<details>
<summary>答案</summary>

a) $\rho_2 = m[W(0) + W(0.8h) + W(0.8h)] = m[1 + 0.2 + 0.2] = 1.4m$

b) 梯度：
   - $\nabla_{x_1} C_2 = -\frac{m}{\rho_0 h} = -\frac{m}{\rho_0 h}$
   - $\nabla_{x_2} C_2 = \frac{m}{\rho_0 h} + \frac{m}{\rho_0 h} = \frac{2m}{\rho_0 h}$
   - $\nabla_{x_3} C_2 = -\frac{m}{\rho_0 h}$

c) $\lambda_2 = -\frac{C_2}{\sum |\nabla C_2|^2} = -\frac{1.4m/\rho_0 - 1}{6m^2/(\rho_0 h)^2}$
</details>

### 挑战题

**习题2.4** 隐式积分器设计
设计一个介于显式和隐式之间的积分器，使其在保持稳定性的同时减少能量耗散。
a) 提出一个混合方案，明确说明哪些项使用显式/隐式处理
b) 分析该方案的稳定性条件
c) 讨论如何自适应地调整混合比例

<details>
<summary>提示</summary>
考虑将力分解为线性部分（弹性力）和非线性部分（碰撞力等），对不同部分使用不同的处理方式。
</details>

<details>
<summary>答案</summary>

a) 混合方案：
   - 线性弹性力使用隐式：$\mathbf{f}_{elastic}^{n+1}$
   - 非线性力使用显式：$\mathbf{f}_{nonlinear}^n$
   - 更新方程：$\mathbf{v}^{n+1} = \mathbf{v}^n + \Delta t[\alpha \mathbf{f}_{elastic}^{n+1} + (1-\alpha)\mathbf{f}_{elastic}^n + \mathbf{f}_{nonlinear}^n]/m$

b) 稳定性分析：
   - 当$\alpha \geq 0.5$时，对线性系统无条件稳定
   - 非线性项仍需满足CFL条件
   - 有效时间步长：$\Delta t_{eff} = \min(\Delta t_{implicit}, \Delta t_{CFL})$

c) 自适应策略：
   - 监测能量变化率：$\alpha = \min(1, \max(0.5, 1 - |dE/dt|/E_{threshold}))$
   - 基于局部刚度：高刚度区域增大$\alpha$
   - 使用误差估计：比较显式和隐式预测的差异
</details>

**习题2.5** SPH表面张力优化
标准CSF方法在计算表面张力时存在数值噪声。设计一个改进的表面张力模型。
a) 分析CSF方法产生噪声的原因
b) 提出至少两种降噪策略
c) 设计一个自适应表面张力系数

<details>
<summary>提示</summary>
考虑法向量计算的数值误差，以及曲率计算的二阶导数性质。
</details>

<details>
<summary>答案</summary>

a) 噪声来源：
   - 颜色场梯度在内部粒子处接近零，导致法向量不稳定
   - 曲率计算涉及二阶导数，放大数值误差
   - 核函数在边界处的截断误差

b) 降噪策略：
   1. 法向量平滑：$\mathbf{n}_i^{smooth} = \sum_j \frac{m_j}{\rho_j} \mathbf{n}_j W_{ij} / \sum_j \frac{m_j}{\rho_j} W_{ij}$
   2. 使用高阶核函数计算曲率
   3. 只在$|\mathbf{n}| > \epsilon_{surface}$的粒子上施加表面张力
   4. 时间平均：$\kappa_i^{filtered} = \alpha \kappa_i^{new} + (1-\alpha)\kappa_i^{old}$

c) 自适应系数：
   $\sigma_{adaptive} = \sigma_0 \cdot f(|\mathbf{n}|) \cdot g(|\nabla \rho|/\rho)$
   其中$f(x) = \tanh(x/\epsilon)$过滤内部粒子，$g(x)$根据密度梯度调整强度
</details>

**习题2.6** 高效邻居搜索数据结构
设计一个针对非均匀粒子分布的自适应邻居搜索结构。
a) 提出数据结构设计
b) 分析时间和空间复杂度
c) 讨论并行化策略

<details>
<summary>提示</summary>
考虑结合多种数据结构的优点，如八叉树的自适应性和网格的简单性。
</details>

<details>
<summary>答案</summary>

a) 混合数据结构：
   - 顶层：松散八叉树，叶节点大小自适应
   - 叶节点：当粒子数>阈值时，使用局部均匀网格
   - 稀疏区域：直接存储粒子列表
   - 密集区域：Z-order哈希表

b) 复杂度分析：
   - 构建：$O(n \log n)$平均情况，$O(n^2)$最坏情况（所有粒子聚集）
   - 查询：$O(k)$平均情况，$k$是邻居数
   - 空间：$O(n + m)$，$m$是活跃节点数

c) 并行化：
   - 构建阶段：自顶向下并行分裂节点
   - 使用Morton码并行排序粒子
   - 查询阶段：每个粒子独立查询，无需同步
   - 动态更新：使用双缓冲避免读写冲突
</details>

**习题2.7** 刚体-流体耦合
设计一个刚体漂浮在SPH流体上的耦合算法。
a) 推导流体对刚体的力和力矩
b) 处理刚体边界条件
c) 保证动量守恒

<details>
<summary>提示</summary>
考虑使用虚拟粒子或者边界积分方法。注意作用力和反作用力。
</details>

<details>
<summary>答案</summary>

a) 力和力矩计算：
   - 压力：$\mathbf{F}_p = -\sum_i \rho_i V_i p_i \nabla W_{ib}$
   - 粘性：$\mathbf{F}_v = \mu \sum_i \frac{m_i}{\rho_i} (\mathbf{v}_b - \mathbf{v}_i) \nabla^2 W_{ib}$
   - 力矩：$\boldsymbol{\tau} = \sum_i (\mathbf{x}_i - \mathbf{x}_{cm}) \times \mathbf{f}_i$

b) 边界条件：
   - 虚拟粒子法：在刚体表面生成虚拟粒子
   - 虚拟粒子速度：$\mathbf{v}_{virtual} = \mathbf{v}_{rigid} + \boldsymbol{\omega} \times \mathbf{r}$
   - 压力镜像：$p_{virtual} = p_{fluid} + \rho g \Delta h$

c) 动量守恒：
   - 使用对称的力计算确保作用力等于反作用力
   - 同时更新流体和刚体：$m_{fluid}\Delta\mathbf{v}_{fluid} + m_{rigid}\Delta\mathbf{v}_{rigid} = 0$
   - 时间积分使用相同的方案
</details>

## 常见陷阱与错误

### 数值不稳定
1. **时间步长过大**：显式积分爆炸，表现为粒子飞散
   - 解决：使用CFL条件自动调整时间步长
   
2. **刚度矩阵病态**：隐式求解不收敛
   - 解决：添加正则化项或使用预条件

3. **粒子聚集**：SPH/PBF中粒子过度聚集
   - 解决：使用人工压力或XSPH粘性

### 性能陷阱
1. **邻居搜索低效**：使用$O(n^2)$的朴素搜索
   - 解决：实现空间数据结构

2. **缓存未命中**：随机内存访问模式
   - 解决：使用空间排序提高局部性

3. **过度同步**：GPU上频繁的全局同步
   - 解决：设计异步算法，减少同步点

### 物理失真
1. **能量不守恒**：长时间仿真后能量漂移
   - 解决：使用辛积分器或能量修正

2. **体积损失**：流体体积逐渐减少
   - 解决：使用DFSPH或体积修正

3. **穿透问题**：高速碰撞时物体穿透
   - 解决：使用CCD或减小时间步长

## 最佳实践检查清单

### 算法选择
- [ ] 根据刚度选择显式/隐式积分器
- [ ] 流体仿真考虑SPH vs PBF的权衡
- [ ] 大变形用拉格朗日，复杂边界用欧拉
- [ ] 实时应用优先考虑稳定性

### 参数调优
- [ ] 时间步长满足稳定性条件
- [ ] 粒子间距与核函数支持域匹配（通常$h = 2\Delta x$）
- [ ] 人工粘性/压力参数经过测试
- [ ] 迭代求解器的收敛阈值合理

### 数据结构
- [ ] 使用空间数据结构加速邻居搜索
- [ ] 考虑内存访问模式优化缓存
- [ ] 粒子数据使用SoA布局
- [ ] 预分配内存避免动态分配

### 并行优化
- [ ] 识别并行化机会（粒子更新、邻居搜索）
- [ ] 最小化原子操作和同步
- [ ] 负载均衡（每个线程处理相似工作量）
- [ ] 利用共享内存减少全局内存访问

### 调试验证
- [ ] 实现能量/动量监测
- [ ] 边界条件正确处理
- [ ] 单元测试核心算法
- [ ] 可视化中间结果
