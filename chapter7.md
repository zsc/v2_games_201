# 第七章：混合欧拉-拉格朗日视角（1）

在前面的章节中，我们分别深入探讨了拉格朗日视角和欧拉视角的物理仿真方法。拉格朗日方法擅长追踪材料的运动轨迹，保持材料属性，但在处理大变形时容易产生网格扭曲；欧拉方法在固定网格上求解，便于处理边界条件和不可压缩性约束，但对流项容易引入数值耗散。本章将介绍结合两种视角优势的混合方法，重点讲解粒子元胞法(PIC)、流体隐粒子(FLIP)、仿射粒子元胞法(APIC)等经典算法。

## 7.1 混合方法概述

混合欧拉-拉格朗日方法的核心思想是：使用拉格朗日粒子携带和传输材料信息，使用欧拉网格进行力的计算和不可压缩性投影。这种设计充分利用了两种视角的优势，避免了各自的缺陷。

### 7.1.1 拉格朗日与欧拉的优缺点

**拉格朗日视角的优势：**
- 自然处理对流：粒子随流体运动，无需显式求解对流项
- 保持材料属性：每个粒子携带密度、温度等属性，无数值耗散
- 易于处理自由表面：粒子的存在即定义了流体区域
- 历史依赖性：可以追踪材料的应变历史（如塑性变形）

**拉格朗日视角的劣势：**
- 大变形问题：网格可能严重扭曲，需要重新网格化
- 不可压缩性投影：在粒子上直接投影困难，需要构建连接关系
- 邻居搜索开销：SPH等方法需要频繁的邻居搜索

**欧拉视角的优势：**
- 固定网格：无需处理网格扭曲，数值格式成熟
- 高效投影：在规则网格上求解泊松方程效率高
- 边界处理简单：固定网格便于施加各种边界条件

**欧拉视角的劣势：**
- 数值耗散：对流项离散化引入人工粘性
- 界面追踪困难：需要额外的Level Set或VOF方法
- 分辨率限制：细节受网格分辨率限制

### 7.1.2 物理量守恒问题

在设计混合方法时，需要特别关注各种物理量的守恒性：

**动量守恒：**
粒子-网格传输过程中必须保证总动量守恒。使用一致的插值核函数，确保：
$$\sum_p m_p \mathbf{v}_p = \sum_i m_i \mathbf{v}_i$$

**角动量守恒：**
简单的速度插值可能不保持角动量。APIC方法通过引入仿射速度场解决了这个问题。

**体积守恒（不可压缩性）：**
通过网格上的压力投影强制满足 $\nabla \cdot \mathbf{u} = 0$。

**能量守恒：**
完全的能量守恒很难实现，通常追求低耗散而非严格守恒。

### 7.1.3 数值耗散与噪声

混合方法需要在数值耗散和噪声之间找到平衡：

- **PIC方法**：过度耗散但稳定，适合需要稳定性的应用
- **FLIP方法**：保持能量但产生噪声，适合湍流等高能量流动
- **PIC/FLIP混合**：通过线性组合在两者间权衡

耗散的来源主要是信息在粒子和网格间传输时的损失。考虑2D情况：
- 粒子自由度：$2N_p$（N_p个粒子，每个2个速度分量）
- 网格自由度：$2N_g$（N_g个网格点）

通常 $N_p \gg N_g$，因此P2G过程必然损失信息。

### 7.1.4 混合策略设计

成功的混合方法需要精心设计以下组件：

1. **插值核函数**：决定粒子和网格间的信息传输
   - 线性：简单高效，但可能产生噪声
   - 二次B样条：更平滑，计算成本适中
   - 三次B样条：最平滑，但计算开销大

2. **传输策略**：
   - 直接传输（PIC）：传输速度值
   - 增量传输（FLIP）：传输速度变化量
   - 仿射传输（APIC）：传输速度及其梯度

3. **时间积分方案**：
   - 显式：简单但受CFL条件限制
   - 隐式：稳定但需要求解线性系统

4. **压力求解**：
   - 在网格上构建和求解压力泊松方程
   - 使用多重网格或共轭梯度法加速

## 7.2 粒子元胞法(PIC)

粒子元胞法(Particle-In-Cell)最早由Harlow在1964年提出，用于求解流体动力学方程。其核心思想是使用拉格朗日粒子携带物理量，在欧拉网格上进行动量方程求解。

### 7.2.1 P2G与G2P传输

PIC的核心是粒子到网格(Particle-to-Grid, P2G)和网格到粒子(Grid-to-Particle, G2P)的传输过程。

**P2G传输（Particle to Grid）：**
将粒子上的物理量传输到网格节点。对于速度场：

$$m_i = \sum_p w_{ip} m_p$$
$$m_i \mathbf{v}_i = \sum_p w_{ip} m_p \mathbf{v}_p$$

其中 $w_{ip}$ 是插值权重，由核函数决定：
$$w_{ip} = N\left(\frac{\mathbf{x}_p - \mathbf{x}_i}{h}\right)$$

**G2P传输（Grid to Particle）：**
将网格上更新后的速度传回粒子：

$$\mathbf{v}_p^{new} = \sum_i w_{ip} \mathbf{v}_i^{new}$$

注意这里直接覆盖粒子速度，这是PIC产生耗散的主要原因。

### 7.2.2 核函数选择

核函数的选择对算法性能和数值性质有重要影响：

**线性核函数（tent function）：**
$$N_1(x) = \begin{cases}
1 - |x| & |x| < 1 \\
0 & \text{otherwise}
\end{cases}$$

- 计算效率高，每个粒子只影响 $2^d$ 个网格点
- 可能产生粒子聚集和速度场不连续

**二次B样条：**
$$N_2(x) = \begin{cases}
\frac{3}{4} - x^2 & |x| < \frac{1}{2} \\
\frac{1}{2}(|x| - \frac{3}{2})^2 & \frac{1}{2} \leq |x| < \frac{3}{2} \\
0 & \text{otherwise}
\end{cases}$$

- 一阶连续，速度场更平滑
- 每个粒子影响 $3^d$ 个网格点

**三次B样条：**
$$N_3(x) = \begin{cases}
\frac{1}{6}(2 - |x|)^3 & 1 < |x| < 2 \\
\frac{2}{3} - x^2 + \frac{1}{2}|x|^3 & |x| \leq 1 \\
0 & \text{otherwise}
\end{cases}$$

- 二阶连续，最平滑
- 计算成本最高，每个粒子影响 $4^d$ 个网格点

### 7.2.3 信息损失分析

PIC方法的主要问题是信息损失导致的数值耗散。考虑一维情况：

设有 $N_p$ 个粒子和 $N_g$ 个网格点。粒子携带 $N_p$ 个速度值，而网格只能存储 $N_g$ 个速度值。当 $N_p > N_g$ 时（通常情况），P2G过程必然损失信息。

**频谱分析：**
PIC相当于一个低通滤波器，每次P2G-G2P循环都会衰减高频成分。设传输算子为 $T$，则：
$$\hat{v}_{k}^{n+1} = T(k) \hat{v}_{k}^n$$

其中 $T(k) < 1$ 对于高波数 $k$，导致高频信息的指数衰减。

**能量耗散率：**
对于二次B样条核函数，单次传输的能量保留率约为：
$$\frac{E_{after}}{E_{before}} \approx 1 - C \left(\frac{\Delta x}{L}\right)^2$$

其中 $L$ 是流动特征长度，$C$ 是依赖于核函数的常数。

### 7.2.4 动量守恒性

PIC方法的一个重要优点是精确的动量守恒。证明如下：

**P2G阶段：**
$$\sum_i m_i \mathbf{v}_i = \sum_i \sum_p w_{ip} m_p \mathbf{v}_p = \sum_p m_p \mathbf{v}_p \sum_i w_{ip}$$

由于核函数的归一化性质 $\sum_i w_{ip} = 1$，因此：
$$\sum_i m_i \mathbf{v}_i = \sum_p m_p \mathbf{v}_p$$

**G2P阶段：**
类似地，如果使用相同的核函数和权重，G2P过程也保持动量守恒。

这种动量守恒性对于长时间仿真的稳定性至关重要。

## 7.3 流体隐粒子(FLIP)

为了解决PIC的过度耗散问题，Brackbill和Ruppel在1986年提出了流体隐粒子方法(Fluid Implicit Particle)。FLIP的关键创新是传输速度的增量而非绝对值。

### 7.3.1 增量式速度更新

FLIP的核心思想是保留粒子的速度历史，只从网格传输速度的变化量：

**标准FLIP更新：**
$$\mathbf{v}_p^{n+1} = \mathbf{v}_p^n + \sum_i w_{ip} (\mathbf{v}_i^{n+1} - \mathbf{v}_i^n)$$

这里 $\mathbf{v}_i^n$ 是P2G后网格上的速度，$\mathbf{v}_i^{n+1}$ 是经过力的作用和压力投影后的速度。

**物理解释：**
- 粒子保持自己的速度"记忆"
- 网格只提供速度的修正量
- 避免了高频信息的损失

### 7.3.2 噪声vs耗散权衡

FLIP解决了耗散问题，但引入了新的挑战：

**噪声来源：**
1. 粒子速度不再通过网格"平均化"
2. 舍入误差和离散化误差会累积
3. 粒子可能发展出网格无法解析的小尺度运动

**噪声特征：**
- 高频、小尺度的速度扰动
- 在剪切流中特别明显
- 可能导致粒子分布不均匀

**能量守恒分析：**
理想情况下，FLIP保持动能：
$$E_{kinetic} = \frac{1}{2} \sum_p m_p |\mathbf{v}_p|^2$$

但实际上，压力投影和力的离散化仍会引入一些耗散。

### 7.3.3 PIC/FLIP混合

为了在稳定性和精度之间取得平衡，通常使用PIC和FLIP的线性组合：

$$\mathbf{v}_p^{n+1} = (1-\alpha) \mathbf{v}_p^{PIC} + \alpha \mathbf{v}_p^{FLIP}$$

其中：
- $\mathbf{v}_p^{PIC} = \sum_i w_{ip} \mathbf{v}_i^{n+1}$ （PIC更新）
- $\mathbf{v}_p^{FLIP} = \mathbf{v}_p^n + \sum_i w_{ip} (\mathbf{v}_i^{n+1} - \mathbf{v}_i^n)$ （FLIP更新）

**混合比例选择：**
- $\alpha = 0$：纯PIC，最稳定但耗散大
- $\alpha = 1$：纯FLIP，最少耗散但可能不稳定
- $\alpha = 0.95-0.99$：常用值，保持大部分能量同时提供一定稳定性

### 7.3.4 自适应混合比例

更先进的方法是根据局部流动特征动态调整混合比例：

**基于涡度的调整：**
$$\alpha = \alpha_{base} + (1 - \alpha_{base}) \cdot \min\left(1, \frac{|\omega|}{|\omega|_{max}}\right)$$

其中 $\omega = \nabla \times \mathbf{v}$ 是涡度。高涡度区域使用更多FLIP以保持旋转运动。

**基于应变率的调整：**
在高剪切区域减少FLIP比例以抑制噪声：
$$\alpha = \alpha_{base} \cdot \exp\left(-\beta \frac{|\mathbf{S}|}{|\mathbf{S}|_{avg}}\right)$$

其中 $\mathbf{S} = \frac{1}{2}(\nabla \mathbf{v} + \nabla \mathbf{v}^T)$ 是应变率张量。

## 7.4 仿射粒子元胞法(APIC)

仿射粒子元胞法(Affine Particle-In-Cell)由Jiang等人在2015年提出，是PIC/FLIP的重要改进。APIC的核心创新是让每个粒子携带局部的仿射速度场，而不仅仅是速度值。

### 7.4.1 仿射速度场

APIC中，每个粒子 $p$ 的速度场定义为：

$$\mathbf{v}(\mathbf{x}) = \mathbf{v}_p + \mathbf{C}_p (\mathbf{x} - \mathbf{x}_p)$$

其中 $\mathbf{C}_p$ 是 $d \times d$ 的矩阵，表示速度梯度。

**物理意义：**
- $\mathbf{v}_p$：粒子中心的速度
- $\mathbf{C}_p$：速度的空间变化率
- 包含了旋转（反对称部分）和变形（对称部分）

**初始化：**
通常初始化 $\mathbf{C}_p = \mathbf{0}$，在仿真过程中自然演化。

### 7.4.2 角动量守恒

APIC的一个关键优势是严格的角动量守恒。考虑单个粒子的角动量：

$$\mathbf{L}_p = \mathbf{x}_p \times m_p \mathbf{v}_p + \mathbf{I}_p \boldsymbol{\omega}_p$$

其中第二项来自粒子的"内部"角动量，$\boldsymbol{\omega}_p$ 与 $\mathbf{C}_p$ 的反对称部分相关。

**守恒性证明要点：**
1. P2G传输保持总角动量
2. 网格上的力（如压力梯度）不产生净力矩
3. G2P传输将角动量正确分配回粒子

数学上可以证明，APIC的传输过程满足：
$$\sum_p \mathbf{x}_p \times m_p \mathbf{v}_p = \text{constant}$$

### 7.4.3 C矩阵的物理意义

$\mathbf{C}_p$ 矩阵可以分解为对称和反对称部分：

$$\mathbf{C}_p = \mathbf{S}_p + \mathbf{W}_p$$

其中：
- $\mathbf{S}_p = \frac{1}{2}(\mathbf{C}_p + \mathbf{C}_p^T)$：应变率（拉伸/压缩）
- $\mathbf{W}_p = \frac{1}{2}(\mathbf{C}_p - \mathbf{C}_p^T)$：旋转率

**与流体力学的联系：**
在连续介质中，速度梯度 $\nabla \mathbf{v}$ 有相同的分解。APIC粒子携带的 $\mathbf{C}_p$ 是局部速度梯度的离散近似。

**C矩阵的更新：**
G2P过程中，$\mathbf{C}_p$ 通过以下公式更新：

$$\mathbf{C}_p = \frac{4}{h^2} \sum_i w_{ip} \mathbf{v}_i (\mathbf{x}_i - \mathbf{x}_p)^T$$

这个公式确保了角动量守恒和二阶精度。

### 7.4.4 APIC的稳定性

APIC相比PIC/FLIP具有更好的稳定性：

**能量分析：**
APIC的能量耗散介于PIC和FLIP之间：
- 比PIC保持更多能量（因为保留了速度梯度信息）
- 比FLIP更稳定（因为通过网格传输提供了隐式滤波）

**数值实验表明：**
1. APIC不需要显式的PIC/FLIP混合
2. 在各种流动条件下都表现稳定
3. 特别适合有旋转的流动（如涡旋）

**与FLIP的关系：**
可以证明，当 $\mathbf{C}_p = \mathbf{0}$ 且使用特定的传输公式时，APIC退化为FLIP。这说明APIC是FLIP的推广。

## 7.5 PolyPIC方法

PolyPIC（Polynomial PIC）是PIC方法的高阶推广，由Fu等人在2017年提出。它使用多项式而非仿射函数来表示粒子周围的速度场，实现了粒子和网格间的无损信息传输。

### 7.5.1 高阶多项式表示

PolyPIC中，每个粒子携带的速度场表示为：

$$\mathbf{v}(\mathbf{x}) = \sum_{|\alpha| \leq k} \mathbf{a}_{\alpha} \frac{(\mathbf{x} - \mathbf{x}_p)^{\alpha}}{\alpha!}$$

其中 $\alpha$ 是多重指标，$k$ 是多项式阶数。

**2D二阶展开示例：**
$$v_x(x,y) = a_{00} + a_{10}(x-x_p) + a_{01}(y-y_p) + a_{20}\frac{(x-x_p)^2}{2} + a_{11}(x-x_p)(y-y_p) + a_{02}\frac{(y-y_p)^2}{2}$$

每个速度分量需要6个系数，2D速度场共需12个系数。

### 7.5.2 Taylor展开传输

PolyPIC的关键是设计P2G和G2P传输，使得信息传输是可逆的。

**自由度匹配：**
- 2D二阶PolyPIC：每个粒子18个自由度（2个速度分量×9个多项式系数）
- 3×3网格模板：18个自由度（9个网格点×2个速度分量）

当自由度匹配时，理论上可以实现无损传输。

**P2G传输：**
$$m_i \mathbf{v}_i = \sum_p \int_{\Omega_p} \rho(\mathbf{x}) \mathbf{v}_p(\mathbf{x}) N_i(\mathbf{x}) d\mathbf{x}$$

其中积分在粒子的影响域 $\Omega_p$ 上进行。

**G2P传输：**
通过最小二乘拟合恢复多项式系数：
$$\min_{\mathbf{a}} \sum_i w_{ip} \left| \mathbf{v}_i - \sum_{|\alpha| \leq k} \mathbf{a}_{\alpha} \frac{(\mathbf{x}_i - \mathbf{x}_p)^{\alpha}}{\alpha!} \right|^2$$

### 7.5.3 精度与效率分析

**精度优势：**
1. 高阶精度：$k$ 阶PolyPIC具有 $(k+1)$ 阶空间精度
2. 保持更多物理特征：可以准确表示涡旋、剪切等复杂流动
3. 减少数值耗散：高阶方法天然具有更低的数值耗散

**计算成本：**
- 存储开销：$O(k^d)$ 个系数per粒子
- P2G成本：需要计算高阶矩 $\int x^{\alpha} N_i(x) dx$
- G2P成本：求解最小二乘系统，可以预计算伪逆

**效率对比（2D情况）：**
| 方法 | 存储/粒子 | P2G FLOPs | G2P FLOPs |
|------|-----------|-----------|-----------|
| PIC | 2 | ~20 | ~20 |
| APIC | 6 | ~60 | ~80 |
| PolyPIC-2 | 18 | ~200 | ~400 |

### 7.5.4 实现细节

**数值积分：**
P2G中的积分通常使用高斯求积：
$$\int_{\Omega_p} f(\mathbf{x}) d\mathbf{x} \approx \sum_{q} w_q f(\mathbf{x}_q)$$

对于二阶多项式，需要至少3×3个积分点。

**条件数问题：**
G2P的最小二乘系统可能病态，特别是当网格点接近共线时。解决方法：
1. 使用SVD求解，截断小奇异值
2. 添加正则化项
3. 使用正交多项式基

**边界处理：**
粒子靠近边界时，影响域可能只覆盖部分网格点，导致欠定系统。处理策略：
1. 降低多项式阶数
2. 使用虚拟网格点
3. 特殊的边界核函数

## 7.6 粒子-网格传输细节

高质量的粒子-网格传输是混合方法成功的关键。本节深入讨论各种实现细节和优化技术。

### 7.6.1 B样条核函数

B样条是最常用的插值核函数，具有良好的数学性质。

**递归定义：**
$$N_n(x) = \frac{x}{n-1} N_{n-1}(x) + \frac{n-x}{n-1} N_{n-1}(x-1)$$

初始条件：$N_1(x) = \mathbb{1}_{[0,1)}(x)$

**重要性质：**
1. 紧支撑：$N_n(x) = 0$ for $x \notin [0, n]$
2. 非负性：$N_n(x) \geq 0$
3. 归一化：$\sum_{i} N_n(x-i) = 1$
4. 光滑性：$N_n \in C^{n-2}$

**导数计算：**
$$N_n'(x) = N_{n-1}(x) - N_{n-1}(x-1)$$

这在计算速度梯度（如APIC的C矩阵）时很有用。

### 7.6.2 二次vs三次插值

选择插值阶数需要在精度和效率间权衡：

**二次B样条（推荐用于大多数应用）：**
- 支撑域：3×3×3（3D）
- C¹连续
- 计算成本适中
- 足够的光滑性

**三次B样条：**
- 支撑域：4×4×4（3D）
- C²连续
- 计算成本高（64vs27个网格点）
- 用于需要高光滑性的情况

**数值比较：**
对于典型的流体仿真，二次和三次的视觉差异通常很小，但三次的计算成本显著更高。

### 7.6.3 核函数的紧支性

紧支性（compact support）是核函数的关键特性，影响算法效率。

**紧支域大小：**
- 线性：$2^d$ 个网格点
- 二次：$3^d$ 个网格点
- 三次：$4^d$ 个网格点

**优化存储访问：**
```
// 2D二次B样条的高效实现
for (int i = base_i; i < base_i + 3; i++) {
    float wi = N2((xp - i*dx) / dx);
    for (int j = base_j; j < base_j + 3; j++) {
        float wj = N2((yp - j*dy) / dy);
        float w = wi * wj;
        // 执行P2G或G2P操作
    }
}
```

**边界处理：**
当粒子靠近域边界时，部分支撑域可能在边界外：
1. 夹clamp索引到有效范围
2. 使用周期边界条件
3. 扩展虚拟网格点

### 7.6.4 并行传输算法

P2G过程涉及多个粒子写入同一网格点，需要处理竞态条件。

**原子操作方法：**
```
// CUDA风格伪代码
atomicAdd(&grid_mass[i], wp * particle_mass);
atomicAdd(&grid_momentum[i], wp * particle_mass * particle_velocity);
```

优点：实现简单，正确性保证
缺点：原子操作开销，可能成为瓶颈

**分块并行（Tiling）：**
1. 将粒子分配到不重叠的块
2. 每个块独立处理其粒子
3. 合并各块的贡献

**排序方法：**
1. 按空间位置排序粒子
2. 相邻粒子可能影响相同网格点
3. 使用共享内存减少全局内存访问

**GPU优化策略：**
- 使用texture内存加速网格数据读取
- Warp级别的协作，减少原子操作
- 混合精度：位置用float，累加用double

## 本章小结

本章介绍了混合欧拉-拉格朗日方法的基本概念和经典算法。关键要点包括：

1. **混合方法的动机**：结合拉格朗日的材料追踪能力和欧拉的高效投影算法
2. **PIC方法**：简单稳定但存在数值耗散，适合需要稳定性的应用
3. **FLIP方法**：通过增量更新减少耗散，但可能产生噪声
4. **APIC方法**：引入仿射速度场，实现角动量守恒，是现代MPM的基础
5. **PolyPIC方法**：使用高阶多项式实现更高精度，但计算成本显著增加
6. **核函数选择**：二次B样条通常是精度和效率的最佳平衡

这些混合方法为下一章的物质点法（MPM）奠定了理论基础。MPM可以看作是将这些概念扩展到固体力学的自然推广。

## 练习题

### 基础题

**练习7.1**：推导二维情况下线性插值核函数的显式表达式，并验证其满足归一化条件。

<details>
<summary>提示</summary>

考虑粒子位于 $(x_p, y_p)$，网格间距为 $h$。对于网格点 $(i,j)$，权重为：
$$w_{ij} = N\left(\frac{x_p - ih}{h}\right) N\left(\frac{y_p - jh}{h}\right)$$
</details>

<details>
<summary>答案</summary>

线性核函数为：
$$N(x) = \max(0, 1 - |x|)$$

二维权重：
$$w_{ij} = \max\left(0, 1 - \left|\frac{x_p - ih}{h}\right|\right) \max\left(0, 1 - \left|\frac{y_p - jh}{h}\right|\right)$$

归一化验证：粒子最多影响4个网格点，设粒子在网格 $(i_0, j_0)$ 和 $(i_0+1, j_0+1)$ 之间，局部坐标 $(\alpha, \beta) \in [0,1]^2$，则：
$$\sum_{i,j} w_{ij} = (1-\alpha)(1-\beta) + \alpha(1-\beta) + (1-\alpha)\beta + \alpha\beta = 1$$
</details>

**练习7.2**：分析PIC方法在一维周期域上的数值耗散。设粒子均匀分布，使用线性插值，计算单次P2G-G2P循环后的能量损失。

<details>
<summary>提示</summary>

考虑正弦速度场 $v(x) = A\sin(kx)$，分析不同波数 $k$ 的衰减率。
</details>

<details>
<summary>答案</summary>

对于波数 $k$ 的正弦模式，传输算子的特征值为：
$$\lambda(k) = \frac{\sin^2(kh/2)}{(kh/2)^2}$$

能量保留率：
$$\frac{E_{after}}{E_{before}} = \lambda^2(k)$$

对于 $kh = \pi$（Nyquist频率），$\lambda = 4/\pi^2 \approx 0.405$，能量保留率仅为 16.4%。
</details>

**练习7.3**：实现FLIP的速度更新公式，并分析在什么情况下FLIP和PIC会给出相同的结果。

<details>
<summary>提示</summary>

考虑特殊的初始条件或流动状态。
</details>

<details>
<summary>答案</summary>

FLIP更新：$\mathbf{v}_p^{n+1} = \mathbf{v}_p^n + \sum_i w_{ip}(\mathbf{v}_i^{n+1} - \mathbf{v}_i^n)$

PIC更新：$\mathbf{v}_p^{n+1} = \sum_i w_{ip}\mathbf{v}_i^{n+1}$

两者相同当且仅当：
$$\mathbf{v}_p^n = \sum_i w_{ip}\mathbf{v}_i^n$$

即粒子速度已经与其插值到网格的速度一致。这发生在：
1. 稳态流动
2. 粒子速度场已经充分光滑
3. 第一个时间步（如果粒子初始化为网格插值）
</details>

### 挑战题

**练习7.4**：证明APIC方法保持角动量守恒。考虑2D情况，写出P2G和G2P过程的详细公式，并证明总角动量不变。

<details>
<summary>提示</summary>

角动量 $L = \sum_p \mathbf{x}_p \times m_p \mathbf{v}_p$。需要证明P2G和G2P过程都保持这个量。
</details>

<details>
<summary>答案</summary>

P2G传输：
$$m_i = \sum_p w_{ip} m_p$$
$$m_i \mathbf{v}_i = \sum_p w_{ip} m_p[\mathbf{v}_p + \mathbf{C}_p(\mathbf{x}_i - \mathbf{x}_p)]$$

关键观察：$\sum_i w_{ip}(\mathbf{x}_i - \mathbf{x}_p) = \mathbf{0}$（一阶矩为零）

因此网格上的总动量：
$$\sum_i m_i \mathbf{v}_i = \sum_p m_p \mathbf{v}_p$$

网格上的角动量：
$$\sum_i \mathbf{x}_i \times m_i \mathbf{v}_i = \sum_p \mathbf{x}_p \times m_p \mathbf{v}_p + \sum_p m_p \sum_i w_{ip}(\mathbf{x}_i - \mathbf{x}_p) \times \mathbf{C}_p(\mathbf{x}_i - \mathbf{x}_p)$$

第二项可以证明为零（使用C矩阵的特殊结构）。G2P过程类似，从而总角动量守恒。
</details>

**练习7.5**：设计一个自适应的PIC/FLIP混合策略，根据局部流动特征（如涡度、散度、应变率）动态调整混合比例。给出具体的公式和物理justification。

<details>
<summary>提示</summary>

考虑不同流动特征对数值稳定性和精度的要求。
</details>

<details>
<summary>答案</summary>

自适应混合策略：
$$\alpha = \alpha_{base} \cdot f_{vorticity} \cdot f_{strain} \cdot f_{divergence}$$

其中：

1. 涡度因子（保持旋转）：
$$f_{vorticity} = \min\left(1, \frac{|\omega|/|\omega|_{ref} + \epsilon}{1 + \epsilon}\right)$$

2. 应变因子（抑制剪切噪声）：
$$f_{strain} = \exp\left(-\beta \frac{||\mathbf{S}||}{||\mathbf{S}||_{avg}}\right)$$

3. 散度因子（压缩区域稳定性）：
$$f_{divergence} = \begin{cases}
1 & \nabla \cdot \mathbf{v} \geq 0 \\
1 - \gamma|\nabla \cdot \mathbf{v}|/|\mathbf{v}|_{max} & \nabla \cdot \mathbf{v} < 0
\end{cases}$$

物理动机：
- 高涡度区需要FLIP保持能量
- 高应变率区需要PIC抑制噪声
- 压缩区域需要更多PIC保证稳定性
</details>

**练习7.6**：推导PolyPIC的最优G2P重构公式。给定网格速度 $\mathbf{v}_i$ 和权重 $w_{ip}$，如何计算粒子的多项式系数以最小化重构误差？

<details>
<summary>提示</summary>

这是一个加权最小二乘问题。构建法方程并求解。
</details>

<details>
<summary>答案</summary>

最小化目标函数：
$$J = \sum_i w_{ip} \left|\mathbf{v}_i - \sum_{|\alpha| \leq k} \mathbf{a}_{\alpha} \frac{(\mathbf{x}_i - \mathbf{x}_p)^{\alpha}}{\alpha!}\right|^2$$

令 $\mathbf{\Phi}$ 为基函数矩阵，$\mathbf{\Phi}_{i,\alpha} = \frac{(\mathbf{x}_i - \mathbf{x}_p)^{\alpha}}{\alpha!}$

法方程：
$$(\mathbf{\Phi}^T \mathbf{W} \mathbf{\Phi}) \mathbf{a} = \mathbf{\Phi}^T \mathbf{W} \mathbf{v}$$

其中 $\mathbf{W} = \text{diag}(w_{ip})$。

解：
$$\mathbf{a} = (\mathbf{\Phi}^T \mathbf{W} \mathbf{\Phi})^{-1} \mathbf{\Phi}^T \mathbf{W} \mathbf{v}$$

实践中，预计算 $(\mathbf{\Phi}^T \mathbf{W} \mathbf{\Phi})^{-1} \mathbf{\Phi}^T \mathbf{W}$ 对于规则网格。
</details>

**练习7.7**：分析三维APIC的存储和计算复杂度。与PIC和FLIP相比，额外的开销是什么？这些开销带来了什么好处？

<details>
<summary>提示</summary>

考虑C矩阵的存储、P2G/G2P的额外计算等。
</details>

<details>
<summary>答案</summary>

存储复杂度：
- PIC/FLIP: 3个浮点数（速度）+ 3个浮点数（位置）= 6个浮点数/粒子
- APIC: 额外9个浮点数（3×3 C矩阵）= 15个浮点数/粒子
- 存储增加: 2.5倍

计算复杂度（使用二次B样条）：
- PIC P2G: 27次权重计算 + 27×3次加法 ≈ 108 FLOPs
- APIC P2G: 额外27×9次乘加（C矩阵贡献）≈ 351 FLOPs
- 计算增加: 约3.25倍

好处：
1. 角动量守恒（无需额外修正）
2. 更好的稳定性（无需PIC/FLIP混合）
3. 更准确的涡旋保持
4. 为MPM提供自然的应力计算框架

开销是值得的，特别是对于需要准确旋转运动的应用。
</details>

**练习7.8**：设计一个实验来定量比较PIC、FLIP、APIC和PolyPIC的性能。应该测试哪些场景？如何量化各方法的优缺点？

<details>
<summary>提示</summary>

考虑不同的流动特征和评价指标。
</details>

<details>
<summary>答案</summary>

测试场景：

1. **Taylor-Green涡旋**：测试能量守恒和涡旋保持
   - 初始条件：$u = \sin(x)\cos(y)$, $v = -\cos(x)\sin(y)$
   - 度量：动能随时间的衰减、涡度场的演化

2. **旋转圆盘**：测试角动量守恒
   - 初始条件：刚体旋转速度场
   - 度量：总角动量误差、速度profile变形

3. **Kelvin-Helmholtz不稳定性**：测试对复杂流动的解析能力
   - 初始条件：剪切层with小扰动
   - 度量：涡旋卷起的时间、小尺度结构

4. **静水平衡**：测试数值噪声
   - 初始条件：静止流体
   - 度量：速度场的RMS误差

评价指标：
- 定量：能量/动量/角动量误差、L2/L∞范数误差
- 定性：视觉质量、稳定性
- 性能：运行时间、内存使用

预期结果：
- PIC：最稳定但耗散最大
- FLIP：能量守恒最好但可能不稳定
- APIC：平衡性能，角动量守恒
- PolyPIC：最高精度但计算成本最高
</details>

## 常见陷阱与错误（Gotchas）

1. **核函数归一化错误**
   - 错误：假设权重自动归一化
   - 正确：某些情况（如边界附近）需要显式归一化
   - 调试：检查 $\sum_i w_{ip}$ 是否为1

2. **P2G原子操作竞态**
   - 错误：非原子累加导致结果不确定
   - 正确：使用原子操作或其他并行策略
   - 调试：单线程vs多线程结果对比

3. **FLIP噪声累积**
   - 错误：纯FLIP在长时间仿真后变得不稳定
   - 正确：使用少量PIC混合或定期重新初始化
   - 调试：监控速度场的高频分量

4. **APIC的C矩阵对称性**
   - 错误：假设C矩阵是对称的
   - 正确：C包含反对称部分（旋转）
   - 调试：分别检查对称和反对称部分

5. **边界粒子处理**
   - 错误：边界附近粒子的权重和不为1
   - 正确：特殊处理或使用ghost cells
   - 调试：可视化边界附近的粒子行为

## 最佳实践检查清单

### 算法选择
- [ ] 稳定性优先→选择PIC或高比例PIC混合
- [ ] 精度优先→选择FLIP或PolyPIC
- [ ] 平衡方案→选择APIC或PIC/FLIP混合
- [ ] 角动量重要→必须选择APIC或PolyPIC

### 实现细节
- [ ] 核函数选择匹配应用需求（通常二次B样条）
- [ ] 正确处理边界条件（ghost cells或特殊权重）
- [ ] 并行化策略避免竞态条件
- [ ] 定期验证守恒量（质量、动量、角动量）

### 性能优化
- [ ] 内存布局优化（AoS vs SoA）
- [ ] 减少原子操作（通过排序或分块）
- [ ] 利用空间局部性（粒子排序）
- [ ] 考虑混合精度计算

### 数值稳定性
- [ ] 监控能量和守恒量
- [ ] 自适应时间步长（CFL条件）
- [ ] 处理退化情况（如粒子聚集）
- [ ] 定期重新分布粒子（如需要）

### 调试技巧
- [ ] 简单案例验证（如平移、旋转）
- [ ] 可视化中间结果（网格速度、粒子分布）
- [ ] 对比不同方法的结果
- [ ] 使用解析解验证（如Taylor-Green涡旋）
