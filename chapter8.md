# 第八章：混合欧拉-拉格朗日视角（2）：物质点法

物质点法（Material Point Method, MPM）是一种强大的数值方法，结合了拉格朗日粒子和欧拉网格的优势。本章将深入探讨MPM的理论基础、经典算法、现代改进以及在各种材料模拟中的应用。我们将从基础的MPM算法开始，逐步深入到MLS-MPM、本构模型、断裂模拟等高级主题，并介绍Taichi中实现MPM的高级特性。

## 8.1 物质点法(MPM)基础

### 8.1.1 MPM的历史与发展

物质点法最初由Sulsky和Schreyer在1996年提出，作为有限元方法(FEM)的扩展来处理大变形问题。传统FEM在处理极大变形时会遇到网格扭曲问题，而MPM通过使用拉格朗日粒子（物质点）携带材料信息，欧拉背景网格进行动量方程求解，巧妙地避免了这个问题。

MPM的发展历程可以分为几个重要阶段：

**早期发展（1994-2000）**：
- 1994年：Sulsky, Chen和Schreyer提出MPM的原始形式，用于固体力学中的冲击和穿透问题
- 1995年：引入GIMP（Generalized Interpolation Material Point）方法，改善了数值稳定性
- 1999年：Bardenhagen等人改进了MPM的动量守恒性质

**理论完善（2000-2010）**：
- 2004年：Steffen等人提出了CPDI（Convected Particle Domain Interpolation），减少了格子噪声
- 2008年：Sadeghirad等人发展了CPDI2，进一步提高了大变形下的精度
- 2010年：引入了双网格MPM，分离了动量和应力的计算

**图形学应用（2013至今）**：
2013年，Stomakhin等人将MPM引入计算机图形学领域，首次实现了雪的真实感模拟，在迪士尼动画《冰雪奇缘》中得到应用。这项工作展示了MPM在处理相变、断裂等复杂物理现象上的独特优势。

关键突破包括：
- 2014年：Jiang等人提出了APIC方法，实现了角动量守恒
- 2016年：Klár等人用MPM模拟沙子，引入了Drucker-Prager塑性模型
- 2016年：Daviet和Bertails-Descoubes提出了implicit MPM
- 2017年：Gao等人发展了adaptive MPM，支持自适应网格细化
- 2018年：Hu等人提出MLS-MPM，将实现简化到88行代码
- 2019年：引入了基于神经网络的本构模型
- 2020年：发展了GPU优化的MPM，实现了实时模拟

MPM的核心思想是将连续介质离散为一系列携带质量、动量、应力等物理量的粒子，而背景网格仅用于计算内力和更新动量。这种双重表示方式使得MPM既保留了拉格朗日方法追踪材料历史的能力，又具备了欧拉方法处理碰撞和自碰撞的便利性。

**应用领域扩展**：
- 地质工程：滑坡、土壤液化、基础沉降
- 生物力学：软组织变形、手术模拟
- 制造业：3D打印、粉末冶金、增材制造
- 影视特效：雪崩、沙尘暴、泥石流、爆炸效果
- 游戏引擎：可破坏环境、流体-固体交互

### 8.1.2 与FEM的关系

MPM可以视为无网格Galerkin方法的一种，特别是属于无单元Galerkin（Element-Free Galerkin, EFG）方法家族。从数学角度看，MPM和FEM都基于连续介质力学的弱形式：

$$\int_\Omega \rho \mathbf{a} \cdot \mathbf{w} \, dV = -\int_\Omega \boldsymbol{\sigma} : \nabla \mathbf{w} \, dV + \int_\Omega \rho \mathbf{b} \cdot \mathbf{w} \, dV + \int_{\partial\Omega_t} \mathbf{t} \cdot \mathbf{w} \, dA$$

其中$\mathbf{w}$是测试函数，$\mathbf{a}$是加速度，$\boldsymbol{\sigma}$是柯西应力张量，$\mathbf{b}$是体力，$\mathbf{t}$是表面力。

**积分方式的本质区别**：

MPM与FEM的主要区别在于积分方式和材料点的处理：

| 特性 | FEM | MPM |
|------|-----|-----|
| 积分点 | 高斯点，位置固定在单元内 | 物质点，可自由移动 |
| 网格作用 | 承载所有信息 | 仅用于动量更新 |
| 拓扑变化 | 需要重新网格化 | 自然处理 |
| 历史变量 | 存储在高斯点 | 存储在粒子上 |
| 大变形 | 网格扭曲问题 | 无网格扭曲 |
| 接触处理 | 需要显式接触算法 | 自动处理 |

**数学等价性**：

在小变形情况下，MPM可以完全等价于FEM。考虑一个单元内有$n_g$个高斯点的FEM和有$n_p$个粒子的MPM：

FEM的积分：
$$\int_{\Omega_e} f(\mathbf{x}) \, dV \approx \sum_{g=1}^{n_g} w_g J_g f(\mathbf{x}_g)$$

MPM的积分：
$$\int_{\Omega_e} f(\mathbf{x}) \, dV \approx \sum_{p=1}^{n_p} V_p f(\mathbf{x}_p)$$

当粒子初始位置与高斯点重合，且$V_p = w_g J_g$时，两者完全等价。

**MPM作为FEM的推广**：

MPM可以理解为使用特殊积分规则的FEM：
1. **动态积分点**：积分点（粒子）随材料移动，自动追踪材料历史
2. **单点积分**：每个粒子使用单点积分，$V_p$是积分权重
3. **无单元结构**：不需要维护单元连接关系，粒子通过背景网格交互

**理论基础的继承**：

MPM继承了FEM的许多理论结果：
- **收敛性**：在网格加密和粒子加密时，MPM收敛到连续解
- **误差估计**：$||u - u_h|| \leq Ch^k$，其中$k$取决于形函数阶数
- **稳定性条件**：CFL条件类似，$\Delta t \leq C\frac{h}{c}$，$c = \sqrt{E/\rho}$
- **守恒性质**：质量、动量自动守恒，APIC/MLS-MPM还保证角动量守恒

**优势与局限**：

MPM相对于FEM的优势：
- 处理超大变形无需重网格化
- 自动处理拓扑变化（断裂、合并）
- 材料界面追踪精确
- 历史相关本构模型实现简单

MPM的局限性：
- 计算成本通常高于FEM（需要更多粒子）
- 边界条件施加不如FEM精确
- 对于小变形问题，FEM更高效
- 可能出现粒子聚集或空洞

这种关系意味着MPM继承了FEM的坚实理论基础，同时又具有处理大变形的灵活性，使其成为极端变形问题的理想选择。

### 8.1.3 弱形式推导

从动量守恒方程出发：

$$\rho \frac{D\mathbf{v}}{Dt} = \nabla \cdot \boldsymbol{\sigma} + \rho \mathbf{b}$$

乘以测试函数$\mathbf{w}$并在域$\Omega$上积分：

$$\int_\Omega \rho \frac{D\mathbf{v}}{Dt} \cdot \mathbf{w} \, dV = \int_\Omega (\nabla \cdot \boldsymbol{\sigma}) \cdot \mathbf{w} \, dV + \int_\Omega \rho \mathbf{b} \cdot \mathbf{w} \, dV$$

对应力项使用分部积分：

$$\int_\Omega (\nabla \cdot \boldsymbol{\sigma}) \cdot \mathbf{w} \, dV = -\int_\Omega \boldsymbol{\sigma} : \nabla \mathbf{w} \, dV + \int_{\partial\Omega} (\boldsymbol{\sigma} \cdot \mathbf{n}) \cdot \mathbf{w} \, dA$$

在MPM中，我们使用粒子近似积分：

$$\int_\Omega f(\mathbf{x}) \, dV \approx \sum_p V_p f(\mathbf{x}_p)$$

其中$V_p$是粒子$p$的体积，$\mathbf{x}_p$是粒子位置。

### 8.1.4 空间离散化

MPM的空间离散化采用混合方法，结合了拉格朗日粒子和欧拉网格的优势：

**1. 粒子表示（拉格朗日部分）**

材料被离散为$N_p$个粒子，每个粒子$p$携带完整的材料状态：

| 变量 | 符号 | 物理意义 | 维度 |
|------|------|----------|------|
| 位置 | $\mathbf{x}_p$ | 当前构型中的位置 | $\mathbb{R}^d$ |
| 速度 | $\mathbf{v}_p$ | 材料点速度 | $\mathbb{R}^d$ |
| 质量 | $m_p$ | 粒子质量（守恒） | $\mathbb{R}$ |
| 体积 | $V_p$ | 当前体积 | $\mathbb{R}$ |
| 初始体积 | $V_p^0$ | 参考构型体积 | $\mathbb{R}$ |
| 变形梯度 | $\mathbf{F}_p$ | 变形梯度张量 | $\mathbb{R}^{d\times d}$ |
| 应力 | $\boldsymbol{\sigma}_p$ | 柯西应力张量 | $\mathbb{R}^{d\times d}$ |
| 仿射速度 | $\mathbf{C}_p$ | APIC/MLS-MPM速度梯度 | $\mathbb{R}^{d\times d}$ |
| 塑性变形 | $\mathbf{F}_p^p$ | 塑性变形部分 | $\mathbb{R}^{d\times d}$ |
| 内部变量 | $\boldsymbol{\alpha}_p$ | 硬化参数等 | 问题相关 |

粒子密度的选择原则：
- 2D：每个网格单元4-9个粒子
- 3D：每个网格单元8-27个粒子
- 自适应：基于变形程度动态调整

**2. 背景网格（欧拉部分）**

使用规则的背景网格进行动量更新，网格节点$i$存储：

| 变量 | 符号 | 作用 | 生命周期 |
|------|------|------|----------|
| 位置 | $\mathbf{x}_i$ | 节点坐标（固定） | 永久 |
| 速度 | $\mathbf{v}_i$ | 节点速度 | 单个时间步 |
| 质量 | $m_i$ | 节点质量 | 单个时间步 |
| 动量 | $\mathbf{p}_i$ | $= m_i \mathbf{v}_i$ | 单个时间步 |
| 力 | $\mathbf{f}_i$ | 节点受力 | 单个时间步 |

网格类型选择：
- **均匀笛卡尔网格**：实现简单，GPU友好
- **MAC网格**：交错网格，用于流体
- **自适应网格**：AMR或八叉树结构

**3. 形函数（插值核）**

形函数$N_i(\mathbf{x})$定义了粒子和网格之间的映射关系。常用B样条基函数：

**线性B样条（tent函数）**：
$$N^1(x) = \begin{cases}
1 - |x| & |x| \leq 1 \\
0 & |x| > 1
\end{cases}$$

导数：
$$\frac{dN^1}{dx} = \begin{cases}
-\text{sign}(x) & |x| < 1 \\
0 & |x| \geq 1
\end{cases}$$

**二次B样条**：
$$N^2(x) = \begin{cases}
\frac{3}{4} - x^2 & |x| \leq \frac{1}{2} \\
\frac{1}{2}(\frac{3}{2} - |x|)^2 & \frac{1}{2} < |x| \leq \frac{3}{2} \\
0 & |x| > \frac{3}{2}
\end{cases}$$

导数：
$$\frac{dN^2}{dx} = \begin{cases}
-2x & |x| \leq \frac{1}{2} \\
-(\frac{3}{2} - |x|)\text{sign}(x) & \frac{1}{2} < |x| \leq \frac{3}{2} \\
0 & |x| > \frac{3}{2}
\end{cases}$$

**三次B样条**：
$$N^3(x) = \begin{cases}
\frac{1}{2}|x|^3 - x^2 + \frac{2}{3} & |x| \leq 1 \\
-\frac{1}{6}|x|^3 + x^2 - 2|x| + \frac{4}{3} & 1 < |x| \leq 2 \\
0 & |x| > 2
\end{cases}$$

**多维形函数**：

使用张量积构造：
$$N_i(\mathbf{x}_p) = \prod_{d=1}^{\text{dim}} N^{1D}\left(\frac{x_p^d - x_i^d}{\Delta x}\right)$$

**4. 离散化误差分析**

空间离散化引入的误差主要来源于：

1. **积分误差**：$O(h^{k+1})$，其中$k$是形函数阶数
2. **插值误差**：$O(h^{k})$，影响P2G/G2P传输
3. **格子噪声**：粒子穿越网格单元边界时的误差

**5. 粒子-网格映射关系**

定义权重函数：
$$w_{ip} = N_i(\mathbf{x}_p)$$

梯度权重：
$$\nabla w_{ip} = \nabla N_i(\mathbf{x}_p)$$

这些权重满足：
- 分割统一性：$\sum_i w_{ip} = 1$
- 紧支性：只有有限个$w_{ip} \neq 0$
- 光滑性：$C^{k-1}$连续，其中$k$是B样条阶数

空间离散化的关键是粒子-网格传输（P2G）和网格-粒子传输（G2P）操作，这些传输保证了动量守恒和数值稳定性。

## 8.2 经典MPM算法

### 8.2.1 应力更新(USL vs USF)

MPM算法的一个关键选择是应力更新的时机，主要有两种策略：

**Update Stress Last (USL)**：
1. 粒子到网格传输（P2G）
2. 网格动量更新
3. 网格到粒子传输（G2P）
4. 在粒子上更新应力

USL算法流程：
```
for p in particles:
    F_p^{n+1} = (I + dt * ∇v^n) * F_p^n
    σ_p^{n+1} = constitutive_model(F_p^{n+1})
```

**Update Stress First (USF)**：
1. 在粒子上更新应力
2. 粒子到网格传输（P2G）
3. 网格动量更新
4. 网格到粒子传输（G2P）

USF算法流程：
```
for p in particles:
    σ_p^n = constitutive_model(F_p^n)
    # 然后进行P2G传输
```

USF通常具有更好的能量守恒性质，但USL在某些情况下数值稳定性更好。选择取决于具体的应用场景和材料模型。

### 8.2.2 形函数选择

形函数的选择对MPM的精度和稳定性有重要影响，需要在计算效率、数值精度和稳定性之间权衡：

**形函数对比**：

| 特性 | 线性B样条 | 二次B样条 | 三次B样条 |
|------|-----------|-----------|------------|
| 支持域 | $[-1, 1]$ | $[-1.5, 1.5]$ | $[-2, 2]$ |
| 影响节点(2D) | 4 | 9 | 16 |
| 影响节点(3D) | 8 | 27 | 64 |
| 连续性 | $C^0$ | $C^1$ | $C^2$ |
| 计算复杂度 | 低 | 中 | 高 |
| 格子噪声 | 严重 | 轻微 | 几乎无 |
| 内存占用 | 小 | 中 | 大 |

**线性B样条（tent函数）**：
- 优点：计算效率高，实现简单，内存占用小
- 缺点：严重的格子噪声（grid crossing error），梯度不连续
- 应用：快速原型开发，实时应用，GPU实现
- 数学表达：$N^1(x) = \max(0, 1 - |x|)$

**二次B样条**：
- 优点：$C^1$连续，格子噪声较小，精度-效率平衡好
- 缺点：计算量是线性的3倍左右
- 应用：大多数生产级MPM实现的默认选择
- 特殊性质：满足再生条件，可精确再生线性场

**三次B样条**：
- 优点：$C^2$连续，几乎无格子噪声，高精度
- 缺点：计算成本高，内存占用大，实现复杂
- 应用：高精度科学计算，小规模精确模拟
- 特殊性质：可精确再生二次多项式场

**格子噪声（Grid Crossing Error）分析**：

当粒子穿越网格单元边界时，权重函数的不连续变化导致的数值噪声：
$$\text{Error} \propto \frac{\partial^{k+1} N}{\partial x^{k+1}}$$

其中$k$是B样条的阶数。高阶B样条具有更高阶的连续性，因此格子噪声更小。

**形函数计算优化**：

多维形函数使用张量积：
$$N_i(\mathbf{x}_p) = \prod_{d=1}^{\text{dim}} N^{1D}\left(\frac{x_p^d - x_i^d}{\Delta x}\right)$$

梯度计算（利用链式法则）：
$$\nabla N_i(\mathbf{x}_p) = \begin{bmatrix}
\frac{\partial N^{1D}_x}{\partial x} N^{1D}_y N^{1D}_z \\
N^{1D}_x \frac{\partial N^{1D}_y}{\partial y} N^{1D}_z \\
N^{1D}_x N^{1D}_y \frac{\partial N^{1D}_z}{\partial z}
\end{bmatrix}$$

**实现技巧**：
1. 预计算权重表：对于固定网格，预计算并存储权重
2. SIMD优化：利用向量指令并行计算多个权重
3. 稀疏性利用：只计算非零权重（利用紧支性）
4. 缓存优化：按粒子空间位置排序，提高缓存命中率

### 8.2.3 积分点与背景网格

MPM中的数值积分使用单点积分规则，每个粒子作为一个积分点，这是MPM区别于传统FEM的关键特征：

**积分近似**：
$$\int_\Omega f(\mathbf{x}) \, dV \approx \sum_p V_p^0 f(\mathbf{x}_p)$$

其中$V_p^0$是粒子的初始（参考）体积，满足：
$$\sum_p V_p^0 = V_{\text{total}}$$

**粒子作为积分点的特性**：
1. **单点积分**：每个粒子使用单点求积规则，可能引入零能模式
2. **权重更新**：$V_p = J_p V_p^0$，其中$J_p = \det(\mathbf{F}_p)$
3. **守恒性**：质量守恒自动满足，$m_p = \rho_0 V_p^0$保持不变
4. **精度分析**：积分精度取决于粒子密度和分布均匀性

**背景网格类型详解**：

**1. 均匀笛卡尔网格**：
- 结构：规则的立方体单元，间距$\Delta x$
- 索引：直接映射$(i,j,k) \rightarrow i + j \cdot N_x + k \cdot N_x \cdot N_y$
- 优点：实现简单，缓存友好，GPU高效
- 缺点：内存浪费（空区域也分配）
- 适用：中小规模、密集型模拟

**2. SPGrid（Sparse Paged Grid）**：
- 结构：虚拟内存页（通常512³或1024³）
- 原理：利用OS的虚拟内存管理，未使用页不占物理内存
- 优点：自动内存管理，支持超大域
- 实现：
```python
# 概念性实现
page_mask = 0xFFFFF000  # 4KB页
offset_mask = 0x00000FFF
page_table = {}  # 稀疏页表
```
- 适用：大规模稀疏模拟（如烟雾）

**3. OpenVDB风格层次网格**：
- 结构：B+树，典型配置5-4-3（根-内部-叶子）
- 分辨率：支持$2^{30}$级别的虚拟分辨率
- 优点：极度稀疏时内存效率最高，支持自适应
- 缺点：随机访问开销大，实现复杂
- 适用：电影级特效，极大规模模拟

**4. 自适应网格（AMR）**：
- 结构：八叉树（3D）或四叉树（2D）
- 细化准则：基于粒子密度、变形梯度或误差估计
- 优点：计算资源集中在关键区域
- 挑战：P2G/G2P跨级别传输复杂

**粒子分布策略**：

**初始分布模式**：
1. **规则分布**：
   - 2D：$2\times2$或$3\times3$每单元
   - 3D：$2\times2\times2$或$3\times3\times3$每单元
   - 优点：均匀，易实现
   - 缺点：可能产生各向异性

2. **随机扰动**：
   ```python
   # Jittered sampling
   x_p = x_regular + random.uniform(-0.5, 0.5) * dx * jitter_factor
   ```
   - jitter_factor通常取0.1-0.3
   - 减少规则分布的伪影

3. **Poisson盘采样**：
   - 保证最小间距$r_{\min}$
   - 蓝噪声特性，更均匀
   - 实现复杂但质量最高

**粒子密度准则**：

| 维度 | 最小密度 | 推荐密度 | 高质量密度 |
|------|----------|----------|------------|
| 1D | 2 ppc | 3 ppc | 4 ppc |
| 2D | 4 ppc | 9 ppc | 16 ppc |
| 3D | 8 ppc | 27 ppc | 64 ppc |

（ppc = particles per cell）

**积分精度与粒子数关系**：
$$\text{Error} = O\left(\frac{1}{\sqrt{N_p}}\right) + O(h^k)$$

第一项是统计误差，第二项是离散化误差。

### 8.2.4 边界条件处理

MPM中的边界条件在网格级别施加，这是欧拉部分处理的优势。边界条件的正确施加对模拟的物理真实性至关重要：

**边界条件类型**：

**1. 粘性边界（Sticky/No-slip）**：
物理意义：完全粘附，模拟粗糙表面或胶合界面
$$\mathbf{v}_i = \mathbf{v}_{\text{wall}} \quad \text{if } \mathbf{x}_i \in \partial\Omega_{\text{sticky}}$$

通常$\mathbf{v}_{\text{wall}} = 0$（静止壁面），但也可以是运动边界。

**2. 滑动边界（Slip/Free-slip）**：
物理意义：无摩擦滑动，只约束法向分量
$$\begin{cases}
\mathbf{v}_i \cdot \mathbf{n} = \mathbf{v}_{\text{wall}} \cdot \mathbf{n} \\
\mathbf{v}_i^{\text{new}} = \mathbf{v}_i - (\mathbf{v}_i \cdot \mathbf{n} - \mathbf{v}_{\text{wall}} \cdot \mathbf{n})\mathbf{n}
\end{cases}$$

**3. 分离边界（Separable/One-way）**：
物理意义：单向约束，允许分离但阻止穿透
$$\mathbf{v}_i^{\text{new}} = \begin{cases}
\mathbf{v}_i - \min(0, \mathbf{v}_i \cdot \mathbf{n})\mathbf{n} & \text{if approaching} \\
\mathbf{v}_i & \text{if separating}
\end{cases}$$

**4. 摩擦边界（Frictional）**：
结合库仑摩擦模型：
$$\begin{cases}
\mathbf{v}_n = -e \cdot \mathbf{v}_i \cdot \mathbf{n} \cdot \mathbf{n} & \text{(法向，e是恢复系数)} \\
\mathbf{v}_t = \max(0, 1 - \mu \frac{|\mathbf{v}_n|}{|\mathbf{v}_t|}) \mathbf{v}_t & \text{(切向，μ是摩擦系数)}
\end{cases}$$

**实现策略**：

```python
@ti.kernel
def apply_boundary_conditions():
    for i, j, k in grid_v:
        # 检查边界
        if i < boundary_width or i >= res_x - boundary_width:
            apply_bc_x(i, j, k)
        if j < boundary_width or j >= res_y - boundary_width:
            apply_bc_y(i, j, k)
        if k < boundary_width or k >= res_z - boundary_width:
            apply_bc_z(i, j, k)

@ti.func
def apply_bc_x(i, j, k):
    if boundary_type == STICKY:
        grid_v[i, j, k] = vec3(0)
    elif boundary_type == SLIP:
        normal = vec3(1, 0, 0) if i < boundary_width else vec3(-1, 0, 0)
        vn = grid_v[i, j, k].dot(normal)
        grid_v[i, j, k] -= vn * normal
    elif boundary_type == SEPARATE:
        normal = vec3(1, 0, 0) if i < boundary_width else vec3(-1, 0, 0)
        vn = grid_v[i, j, k].dot(normal)
        if vn * normal.x < 0:  # 朝向边界
            grid_v[i, j, k] -= vn * normal
    elif boundary_type == FRICTION:
        apply_friction_bc(i, j, k)
```

**复杂边界几何**：

**1. 隐式表面（Level Set）**：
使用有符号距离场$\phi(\mathbf{x})$表示边界：
```python
@ti.func
def apply_sdf_boundary(x, v):
    phi = sample_sdf(x)
    if phi < 0:  # 在物体内部
        normal = compute_sdf_gradient(x).normalized()
        vn = v.dot(normal)
        if boundary_type == STICKY:
            v = vec3(0)
        elif boundary_type == SLIP:
            v -= vn * normal
    return v
```

**2. 解析边界**：
球体、平面、盒子等简单几何：
```python
@ti.func
def sphere_boundary(x, v, center, radius):
    dist = (x - center).norm()
    if dist < radius:
        normal = (x - center).normalized()
        # 应用边界条件
        v = apply_bc(v, normal)
    return v
```

**3. 三角网格边界**：
复杂几何使用三角网格表示，需要加速结构（BVH、空间哈希）。

**边界条件的时机**：

正确的施加顺序：
1. P2G传输
2. 网格动量更新（重力、内力）
3. **施加边界条件** ← 关键时机
4. G2P传输

**常见问题与解决**：

1. **粒子穿透**：
   - 原因：时间步过大或边界太薄
   - 解决：CFL条件、多层边界、连续碰撞检测

2. **粘附伪影**：
   - 原因：数值粘性
   - 解决：高阶形函数、FLIP混合

3. **边界层分离**：
   - 原因：边界处粒子稀疏
   - 解决：边界附近增加粒子密度

4. **动量不守恒**：
   - 原因：边界力未正确计算
   - 解决：记录边界冲量，用于后处理分析

## 8.3 移动最小二乘MPM(MLS-MPM)

### 8.3.1 MLS插值理论

MLS-MPM是Hu等人在2018年提出的MPM简化版本，核心思想是使用移动最小二乘（Moving Least Squares）插值替代传统的B样条插值。

MLS插值的目标是找到一个多项式$p(\mathbf{x})$，最小化加权误差：

$$\min_p \sum_i w_i(\mathbf{x}) [p(\mathbf{x}_i) - f_i]^2$$

其中$w_i(\mathbf{x})$是权重函数。

对于一阶MLS（线性），我们寻找：
$$p(\mathbf{x}) = \mathbf{a}^T \mathbf{x} + b$$

解得系数后，MLS形函数为：
$$\phi_i(\mathbf{x}) = w_i(\mathbf{x}) \mathbf{q}^T(\mathbf{x}) \mathbf{M}^{-1}(\mathbf{x}) \mathbf{q}(\mathbf{x}_i)$$

其中$\mathbf{q}(\mathbf{x}) = [1, x, y, z]^T$是基函数向量。

### 8.3.2 88行实现解析

MLS-MPM的一个重要贡献是极简的实现。以下是核心算法的伪代码解析：

```python
# 初始化
for p in particles:
    x_p = initial_position[p]
    v_p = initial_velocity[p]
    F_p = I  # 变形梯度初始为单位矩阵
    C_p = 0  # 仿射速度场矩阵

# 主循环
for step in time_steps:
    # 清空网格
    grid_v = 0
    grid_m = 0
    
    # P2G: 粒子到网格
    for p in particles:
        base = (x_p / dx - 0.5).int()
        fx = x_p / dx - base
        
        # 二次B样条权重
        w = [0.5 * (1.5 - fx)^2,
             0.75 - (fx - 1)^2,
             0.5 * (fx - 0.5)^2]
        
        # 计算应力（本构模型）
        stress = constitutive_model(F_p)
        affine = stress + mass * C_p
        
        # 传输到网格
        for offset in 3x3x3:
            i = base + offset
            weight = w[offset.x] * w[offset.y] * w[offset.z]
            grid_v[i] += weight * (mass * v_p + affine * (x_i - x_p))
            grid_m[i] += weight * mass
    
    # 网格更新
    for i in grid:
        if grid_m[i] > 0:
            grid_v[i] /= grid_m[i]  # 动量转速度
            grid_v[i] += dt * gravity  # 重力
            # 施加边界条件
            apply_boundary_conditions(grid_v[i])
    
    # G2P: 网格到粒子
    for p in particles:
        v_p = 0
        C_p = 0
        
        for offset in 3x3x3:
            i = base + offset
            weight = w[offset.x] * w[offset.y] * w[offset.z]
            v_p += weight * grid_v[i]
            C_p += 4 * weight * grid_v[i] * (x_i - x_p)^T / dx^2
        
        # 更新位置和变形梯度
        x_p += dt * v_p
        F_p = (I + dt * C_p) * F_p
```

关键简化：
1. 直接使用$\mathbf{C}_p$矩阵近似速度梯度$\nabla \mathbf{v}$
2. 变形梯度更新：$\mathbf{F}^{n+1} = (\mathbf{I} + \Delta t \mathbf{C}^n) \mathbf{F}^n$
3. 体积更新：$J_{p} *= 1 + \Delta t \cdot \text{trace}(\mathbf{C})$

### 8.3.3 性能优势分析

MLS-MPM相比传统MPM的性能优势：

1. **计算复杂度降低**：
   - 传统MPM需要计算$\nabla N_i$：约30 FLOPs/粒子
   - MLS-MPM直接使用$\mathbf{C}_p$：约15 FLOPs/粒子
   - FLOPs减少约50%

2. **内存访问优化**：
   - 减少了梯度计算的内存访问
   - 更好的缓存局部性
   - 适合GPU实现

3. **数值稳定性**：
   - 自然满足角动量守恒
   - 减少了数值耗散
   - 更稳定的大时间步长

4. **实现简洁性**：
   - 代码行数大幅减少
   - 易于理解和调试
   - 便于扩展到不同材料模型

### 8.3.4 与APIC的关系

MLS-MPM本质上是APIC（Affine Particle-in-Cell）方法在弹性固体上的应用：

**APIC回顾**：
- 粒子携带仿射速度场：$\mathbf{v}(\mathbf{x}) = \mathbf{v}_p + \mathbf{C}_p(\mathbf{x} - \mathbf{x}_p)$
- 保证角动量守恒

**MLS-MPM的创新**：
1. 将APIC的动量传输用于MPM
2. 在P2G阶段融合了弹性力计算
3. 统一了流体和固体的处理框架

数学关系：
$$\mathbf{C}_p^{APIC} = \mathbf{C}_p^{MLS-MPM} \quad \text{(速度梯度矩阵相同)}$$

$$\text{MLS-MPM} = \text{APIC} + \text{弹性力} + \text{本构模型}$$

这种统一使得MLS-MPM可以无缝处理流固耦合问题。

## 8.4 MPM中的本构模型

### 8.4.1 弹性固体(Neo-Hookean, Corotated)

MPM可以轻松处理各种超弹性材料模型。最常用的是Neo-Hookean和Corotated模型。

**Neo-Hookean模型**：
应变能密度函数：
$$\psi(\mathbf{F}) = \frac{\mu}{2}(\text{tr}(\mathbf{F}^T\mathbf{F}) - d) - \mu\ln(J) + \frac{\lambda}{2}\ln^2(J)$$

其中$J = \det(\mathbf{F})$，$d$是空间维度，$\mu$和$\lambda$是Lamé参数：
$$\mu = \frac{E}{2(1+\nu)}, \quad \lambda = \frac{E\nu}{(1+\nu)(1-2\nu)}$$

第一Piola-Kirchhoff应力：
$$\mathbf{P}(\mathbf{F}) = \mu(\mathbf{F} - \mathbf{F}^{-T}) + \lambda\ln(J)\mathbf{F}^{-T}$$

**Corotated模型**：
使用极分解$\mathbf{F} = \mathbf{R}\mathbf{S}$，其中$\mathbf{R}$是旋转矩阵，$\mathbf{S}$是对称矩阵。

应变能密度函数：
$$\psi(\mathbf{F}) = \mu\sum_i(\sigma_i - 1)^2 + \frac{\lambda}{2}(J - 1)^2$$

其中$\sigma_i$是$\mathbf{F}$的奇异值。

第一Piola-Kirchhoff应力：
$$\mathbf{P}(\mathbf{F}) = 2\mu(\mathbf{F} - \mathbf{R}) + \lambda(J - 1)J\mathbf{F}^{-T}$$

Corotated模型更适合大旋转小应变的情况，而Neo-Hookean更适合大变形。

### 8.4.2 弱可压缩流体

对于流体模拟，MPM可以使用简化的状态方程模型。

**状态方程**：
$$p = K(1 - J)$$

其中$K$是体积模量，$J = \det(\mathbf{F})$是体积比。

**流体的特殊处理**：
1. 只维护$J$而不是完整的$\mathbf{F}$矩阵
2. 避免数值不稳定性
3. 压力贡献：$\mathbf{P} = -pJ\mathbf{F}^{-T}$

**粘性项**：
可以添加粘性应力：
$$\boldsymbol{\tau}_{\text{viscous}} = \mu(\nabla\mathbf{v} + \nabla\mathbf{v}^T)$$

在P2G阶段直接使用$\mathbf{C}_p$计算粘性力。

### 8.4.3 弹塑性材料

弹塑性材料在超过屈服极限后会产生永久变形。MPM中常用乘法分解：

$$\mathbf{F} = \mathbf{F}^e \mathbf{F}^p$$

其中$\mathbf{F}^e$是弹性变形，$\mathbf{F}^p$是塑性变形。

**屈服准则**：
判断材料是否进入塑性状态。常用的有：

1. **von Mises准则**（金属）：
   $$f = \sqrt{\frac{3}{2}\mathbf{s}:\mathbf{s}} - \sigma_Y \leq 0$$
   其中$\mathbf{s}$是偏应力张量，$\sigma_Y$是屈服应力。

2. **Drucker-Prager准则**（砂土）：
   $$f = \sqrt{J_2} + \alpha I_1 - k \leq 0$$
   其中$I_1$是第一应力不变量，$J_2$是第二偏应力不变量。

**返回映射算法**：
当应力超过屈服面时，需要将应力投影回屈服面上：

```python
# SVD分解
U, Sigma, V = svd(F_elastic)

# von Mises返回映射
for i in range(dim):
    if Sigma[i] > yield_surface:
        Sigma[i] = yield_surface
    elif Sigma[i] < -yield_surface:
        Sigma[i] = -yield_surface

# 重构弹性变形梯度
F_elastic = U @ diag(Sigma) @ V.T
```

### 8.4.4 屈服准则与塑性流动

塑性流动的方向由流动法则决定：

**关联流动法则**：
$$\dot{\boldsymbol{\varepsilon}}^p = \dot{\lambda} \frac{\partial f}{\partial \boldsymbol{\sigma}}$$

其中$\dot{\lambda}$是塑性乘子，$f$是屈服函数。

**具体模型实例**：

1. **雪的模型**（Stomakhin 2013）：
   - 使用可变的临界压缩和拉伸比
   - 硬化参数随塑性变形增加
   - $\theta_c = \theta_{c0}(1 + \xi\max(0, \varepsilon_p))$

2. **沙子模型**（Klár 2016）：
   - Drucker-Prager屈服准则
   - 体积保持：$\det(\mathbf{F}^p) = 1$
   - 摩擦角$\phi$和内聚力$c$

3. **泥土模型**（Cam-Clay）：
   - 椭圆形屈服面
   - 硬化/软化行为
   - 临界状态线

**SVD夹持技巧**：
MPM中常用SVD分解处理塑性：

$$\mathbf{F} = \mathbf{U}\boldsymbol{\Sigma}\mathbf{V}^T$$

通过夹持奇异值$\boldsymbol{\Sigma}$到屈服面内，可以简单地实现各种屈服准则：

```python
# Box屈服准则（用于雪）
for i in range(dim):
    Sigma[i] = clamp(Sigma[i], 1-theta_c, 1+theta_s)

# Drucker-Prager（用于沙子）
if trace(Sigma) > 0:  # 压缩
    # 投影到屈服面
    project_to_yield_surface(Sigma)
```

## 8.5 MPM中的拉格朗日力

### 8.5.1 弹簧与阻尼器

MPM可以结合拉格朗日力来模拟额外的约束和相互作用。最常见的是弹簧-阻尼系统。

**弹簧力模型**：
两个粒子之间的弹簧力：
$$\mathbf{f}_{spring} = -k_s(||\mathbf{x}_i - \mathbf{x}_j|| - L_0)\frac{\mathbf{x}_i - \mathbf{x}_j}{||\mathbf{x}_i - \mathbf{x}_j||}$$

其中$k_s$是弹簧刚度，$L_0$是原始长度。

**阻尼力**：
$$\mathbf{f}_{damp} = -k_d(\mathbf{v}_i - \mathbf{v}_j) \cdot \frac{\mathbf{x}_i - \mathbf{x}_j}{||\mathbf{x}_i - \mathbf{x}_j||} \cdot \frac{\mathbf{x}_i - \mathbf{x}_j}{||\mathbf{x}_i - \mathbf{x}_j||}$$

**实现方式**：
1. 在P2G之前计算拉格朗日力
2. 将力转换为动量变化
3. 在P2G阶段传输到网格

```python
# 计算弹簧力
for spring in springs:
    i, j = spring.endpoints
    direction = (x[j] - x[i]).normalized()
    stretch = (x[j] - x[i]).norm() - spring.rest_length
    
    # 弹簧力
    f_spring = spring.stiffness * stretch * direction
    
    # 阻尼力
    v_rel = v[j] - v[i]
    f_damp = spring.damping * dot(v_rel, direction) * direction
    
    # 应用力
    forces[i] += f_spring + f_damp
    forces[j] -= f_spring + f_damp
```

### 8.5.2 薄壳与薄膜

MPM可以处理嵌入在高维空间中的低维结构。

**薄膜模型**：
使用面内应变能：
$$\psi_{membrane} = \frac{h}{2}\int_S \mathbf{E} : \mathbb{C} : \mathbf{E} \, dS$$

其中$h$是厚度，$\mathbf{E}$是Green应变张量，$\mathbb{C}$是弹性张量。

**弯曲能**：
对于薄壳，需要添加弯曲项：
$$\psi_{bending} = \frac{h^3}{24}\int_S \kappa^2 \, dS$$

其中$\kappa$是曲率。

**离散化方法**：
1. 使用三角网格表示薄壳
2. 每个三角形计算应变能
3. 将力贡献传输到MPM粒子

**耦合策略**：
- 单向耦合：薄壳影响流体但不受流体影响
- 双向耦合：通过网格交换动量

### 8.5.3 刚体约束

MPM可以与刚体动力学耦合，实现刚体-变形体相互作用。

**刚体表示**：
- 质心位置：$\mathbf{x}_c$
- 姿态四元数：$\mathbf{q}$
- 线速度：$\mathbf{v}_c$
- 角速度：$\boldsymbol{\omega}$

**约束施加**：
对于附着在刚体上的粒子：
$$\mathbf{v}_p = \mathbf{v}_c + \boldsymbol{\omega} \times (\mathbf{x}_p - \mathbf{x}_c)$$

**碰撞处理**：
1. 检测粒子与刚体的碰撞
2. 计算冲量响应
3. 更新粒子和刚体的速度

```python
# 刚体约束
for p in constrained_particles:
    # 计算刚体上对应点的速度
    r = x[p] - rigid_body.center
    v_rigid = rigid_body.velocity + cross(rigid_body.omega, r)
    
    # 约束粒子速度
    v[p] = v_rigid
```

### 8.5.4 接触力模型

MPM天然处理自碰撞，但有时需要额外的接触力模型。

**罚函数方法**：
当检测到穿透时，施加罚力：
$$\mathbf{f}_{penalty} = k_{penalty} \cdot d_{penetration} \cdot \mathbf{n}$$

**摩擦力**：
使用库仑摩擦模型：
$$||\mathbf{f}_t|| \leq \mu ||\mathbf{f}_n||$$

其中$\mu$是摩擦系数。

**粘附力**：
模拟材料之间的粘附：
$$\mathbf{f}_{adhesion} = -k_{adhesion} \cdot A_{contact} \cdot \mathbf{n}$$

## 8.6 数值断裂

### 8.6.1 连续损伤力学(CDM)

连续损伤力学通过引入损伤变量$d \in [0,1]$来描述材料的退化。

**损伤演化方程**：
$$\dot{d} = \frac{\langle f(\varepsilon) - \kappa \rangle}{\eta}$$

其中$f$是损伤加载函数，$\kappa$是损伤阈值，$\eta$是粘性参数。

**有效应力**：
$$\boldsymbol{\sigma}_{effective} = (1-d)\boldsymbol{\sigma}_{undamaged}$$

**实现步骤**：
1. 计算应变或应力指标
2. 更新损伤变量
3. 修正应力响应
4. 当$d \approx 1$时，粒子失效

```python
# 损伤更新
for p in particles:
    # 计算等效应变
    epsilon_eq = compute_equivalent_strain(F[p])
    
    # 更新损伤变量
    if epsilon_eq > damage_threshold:
        damage[p] = min(1.0, damage[p] + dt * damage_rate * (epsilon_eq - damage_threshold))
    
    # 修正应力
    stress[p] *= (1 - damage[p])
    
    # 检查失效
    if damage[p] > 0.99:
        particle_active[p] = False
```

### 8.6.2 相场断裂模型

相场方法使用连续场$\phi \in [0,1]$表示裂纹，避免了显式追踪裂纹表面。

**能量泛函**：
$$\Psi = \int_\Omega g(\phi)\psi^e(\boldsymbol{\varepsilon}) \, dV + G_c\int_\Omega \left(\frac{\phi^2}{2l} + \frac{l}{2}|\nabla\phi|^2\right) \, dV$$

其中$g(\phi) = (1-\phi)^2$是退化函数，$G_c$是断裂能，$l$是长度尺度。

**演化方程**：
$$\frac{\partial \phi}{\partial t} = -M\frac{\delta \Psi}{\delta \phi}$$

**MPM实现**：
1. 在粒子上存储相场值
2. 通过网格求解相场方程
3. 根据相场值修正材料响应

### 8.6.3 CPIC方法

Compatible PIC (CPIC)方法通过维护粒子之间的连接关系来处理断裂。

**连接矩阵**：
维护粒子对之间的连接强度：
$$C_{ij} \in [0,1]$$

**断裂准则**：
当应变超过临界值时，减少连接强度：
$$C_{ij} = \max(0, C_{ij} - \Delta t \cdot R(\varepsilon_{ij}))$$

**力的计算**：
只在连接的粒子之间传递力：
$$\mathbf{f}_i = \sum_j C_{ij} \mathbf{f}_{ij}$$

### 8.6.4 断裂准则

不同的断裂准则适用于不同的材料：

**最大主应力准则**：
$$\sigma_1 > \sigma_{critical}$$

适用于脆性材料。

**最大主应变准则**：
$$\varepsilon_1 > \varepsilon_{critical}$$

适用于延性材料。

**能量释放率准则**：
$$G > G_c$$

基于Griffith断裂理论。

**J积分准则**：
用于评估裂纹尖端的应力强度因子。

**实现示例**：
```python
# 基于主应力的断裂
for p in particles:
    # SVD分解得到主应力
    U, S, V = svd(stress[p])
    max_principal_stress = S.max()
    
    # 检查断裂准则
    if max_principal_stress > fracture_threshold:
        # 标记粒子断裂
        fractured[p] = True
        
        # 释放应力
        stress[p] *= stress_release_factor
        
        # 可选：生成裂纹粒子
        if generate_crack_particles:
            create_crack_particles(x[p], U[:, 0])  # 沿主应力方向
```

## 8.7 隐式MPM

### 8.7.1 隐式时间积分

当模拟刚性材料（高杨氏模量）或需要大时间步长时，显式MPM会遇到稳定性问题。隐式MPM通过求解非线性系统来保证无条件稳定性。

**隐式欧拉格式**：
$$\mathbf{v}^{n+1} = \mathbf{v}^n + \Delta t \mathbf{M}^{-1}[\mathbf{f}_{ext} - \mathbf{f}_{int}(\mathbf{x}^{n+1})]$$
$$\mathbf{x}^{n+1} = \mathbf{x}^n + \Delta t \mathbf{v}^{n+1}$$

这导致非线性系统：
$$\mathbf{g}(\mathbf{v}^{n+1}) = \mathbf{M}\mathbf{v}^{n+1} - \mathbf{M}\mathbf{v}^n - \Delta t[\mathbf{f}_{ext} - \mathbf{f}_{int}(\mathbf{x}^n + \Delta t \mathbf{v}^{n+1})] = 0$$

**Newton-Raphson求解**：
线性化系统：
$$\mathbf{J}\Delta\mathbf{v} = -\mathbf{g}(\mathbf{v}^k)$$

其中雅可比矩阵：
$$\mathbf{J} = \frac{\partial \mathbf{g}}{\partial \mathbf{v}} = \mathbf{M} + \Delta t^2 \frac{\partial \mathbf{f}_{int}}{\partial \mathbf{x}}$$

### 8.7.2 切线刚度矩阵

计算切线刚度矩阵是隐式MPM的关键。对于超弹性材料：

**应力导数**：
$$\frac{\partial \mathbf{P}}{\partial \mathbf{F}} = \frac{\partial^2 \psi}{\partial \mathbf{F} \partial \mathbf{F}} = \mathbb{C}$$

这是四阶张量，称为切线模量。

**Neo-Hookean的切线模量**：
$$\mathbb{C}_{ijkl} = \mu\delta_{ik}\delta_{jl} + (\mu - \lambda\ln J)F^{-1}_{ji}F^{-1}_{lk} + \lambda F^{-1}_{ji}F^{-1}_{lk}$$

**网格级别的刚度矩阵**：
$$\mathbf{K}_{IJ} = \sum_p V_p \nabla N_I(\mathbf{x}_p) : \mathbb{C}_p : \nabla N_J(\mathbf{x}_p)$$

### 8.7.3 Newton-Raphson迭代

隐式MPM的Newton迭代过程：

```python
def implicit_mpm_step():
    # 初始猜测
    v_new = v_old + dt * f_ext / mass
    
    for newton_iter in range(max_iterations):
        # 计算残差
        x_trial = x_old + dt * v_new
        f_int = compute_internal_forces(x_trial)
        residual = mass * (v_new - v_old) - dt * (f_ext - f_int)
        
        # 检查收敛
        if norm(residual) < tolerance:
            break
        
        # 计算切线刚度矩阵
        K = compute_tangent_stiffness(x_trial)
        J = mass + dt^2 * K
        
        # 求解线性系统
        delta_v = solve(J, -residual)
        
        # 线搜索（可选）
        alpha = line_search(v_new, delta_v)
        v_new += alpha * delta_v
    
    return v_new
```

**线搜索策略**：
为确保收敛，使用Armijo线搜索：
$$||\mathbf{g}(\mathbf{v} + \alpha\Delta\mathbf{v})|| < (1 - c\alpha)||\mathbf{g}(\mathbf{v})||$$

其中$c \in (0, 0.5)$是常数。

### 8.7.4 大变形处理

大变形下的隐式MPM需要特殊处理：

**几何非线性**：
变形梯度的增量更新：
$$\mathbf{F}^{n+1} = (\mathbf{I} + \Delta t \nabla \mathbf{v}^{n+1})\mathbf{F}^n$$

这是关于$\mathbf{v}^{n+1}$的非线性函数。

**共旋框架**：
使用共旋坐标系减少非线性：
1. 提取旋转：$\mathbf{F} = \mathbf{R}\mathbf{S}$
2. 在局部坐标系求解
3. 转换回全局坐标系

**增量形式**：
使用增量变形梯度：
$$\delta\mathbf{F} = \Delta t \nabla\delta\mathbf{v} \cdot \mathbf{F}^n$$

**自适应时间步长**：
根据Newton迭代的收敛性调整时间步长：
```python
if newton_iterations > max_iter * 0.8:
    dt *= 0.5  # 减小时间步长
elif newton_iterations < max_iter * 0.2:
    dt *= 1.2  # 增大时间步长
```

## 8.8 高级Taichi特性（2）

### 8.8.1 稀疏数据结构设计

MPM的计算域通常是稀疏的，使用稀疏数据结构可以大幅提升性能。

**SPGrid (Setaluri et al. 2014)**：
利用虚拟内存的分页机制：
```python
@ti.data_oriented
class SPGrid:
    def __init__(self, resolution):
        self.resolution = resolution
        # 使用位掩码管理稀疏块
        self.mask = ti.field(dtype=ti.i32, shape=resolution//block_size)
        self.data = ti.field(dtype=ti.f32)
        
        # 动态分配
        ti.root.pointer(ti.ijk, resolution//block_size).dense(ti.ijk, block_size).place(self.data)
```

**优势**：
- 内存自动管理
- 缓存友好的访问模式
- 支持动态拓扑变化

**OpenVDB风格的层次结构**：
```python
# Taichi中的层次稀疏网格
grid = ti.field(dtype=ti.f32)
ti.root.pointer(ti.ijk, 64).pointer(ti.ijk, 8).dense(ti.ijk, 8).place(grid)
```

三层结构：
1. 顶层：64³的指针网格
2. 中层：8³的指针块
3. 底层：8³的密集数据

### 8.8.2 动态内存分配

Taichi支持动态内存管理，适合粒子数量变化的场景。

**动态粒子数组**：
```python
@ti.data_oriented
class ParticleSystem:
    def __init__(self, max_particles):
        self.max_particles = max_particles
        self.n_particles = ti.field(ti.i32, shape=())
        
        # 动态数组
        self.x = ti.Vector.field(3, dtype=ti.f32)
        self.v = ti.Vector.field(3, dtype=ti.f32)
        self.F = ti.Matrix.field(3, 3, dtype=ti.f32)
        
        # 使用动态SNode
        self.particle = ti.root.dynamic(ti.i, self.max_particles)
        self.particle.place(self.x, self.v, self.F)
    
    @ti.kernel
    def add_particle(self, pos: ti.types.vector(3, ti.f32)):
        i = ti.atomic_add(self.n_particles[None], 1)
        self.x[i] = pos
        self.v[i] = ti.Vector([0.0, 0.0, 0.0])
        self.F[i] = ti.Matrix.identity(ti.f32, 3)
```

**内存池管理**：
```python
@ti.data_oriented
class MemoryPool:
    def __init__(self, chunk_size, max_chunks):
        self.chunk_size = chunk_size
        self.free_list = ti.field(ti.i32, shape=max_chunks)
        self.n_free = ti.field(ti.i32, shape=())
        
    @ti.func
    def allocate(self) -> ti.i32:
        n = ti.atomic_sub(self.n_free[None], 1)
        if n > 0:
            return self.free_list[n-1]
        return -1  # 分配失败
    
    @ti.func
    def deallocate(self, chunk_id: ti.i32):
        n = ti.atomic_add(self.n_free[None], 1)
        self.free_list[n] = chunk_id
```

### 8.8.3 层次化网格

层次化网格可以在不同尺度上高效处理物理现象。

**自适应网格细化(AMR)**：
```python
@ti.data_oriented
class AdaptiveGrid:
    def __init__(self):
        # 多层级网格
        self.level0 = ti.field(dtype=ti.f32, shape=(64, 64, 64))
        self.level1 = ti.field(dtype=ti.f32, shape=(128, 128, 128))
        self.level2 = ti.field(dtype=ti.f32, shape=(256, 256, 256))
        
        # 细化标记
        self.refine_flag = ti.field(dtype=ti.i32, shape=(64, 64, 64))
    
    @ti.kernel
    def adaptive_refine(self):
        for i, j, k in self.level0:
            # 检查细化准则（如梯度）
            if self.needs_refinement(i, j, k):
                self.refine_flag[i, j, k] = 1
                # 插值到细网格
                self.interpolate_to_fine(i, j, k)
    
    @ti.func
    def needs_refinement(self, i, j, k) -> ti.i32:
        # 基于梯度的细化准则
        grad = ti.abs(self.level0[i+1,j,k] - self.level0[i-1,j,k])
        return grad > self.refine_threshold
```

**多重网格加速**：
用于隐式求解器：
```python
@ti.kernel
def multigrid_vcycle(level: ti.i32):
    # 前光滑
    for _ in range(n_smooth):
        smooth(level)
    
    if level > 0:
        # 限制到粗网格
        restrict(level, level-1)
        
        # 递归求解
        multigrid_vcycle(level-1)
        
        # 延拓到细网格
        prolongate(level-1, level)
    
    # 后光滑
    for _ in range(n_smooth):
        smooth(level)
```

### 8.8.4 自定义数据布局

优化内存访问模式对性能至关重要。

**AoS vs SoA布局**：
```python
# Array of Structures (AoS)
@ti.dataclass
class Particle:
    x: ti.types.vector(3, ti.f32)
    v: ti.types.vector(3, ti.f32)
    m: ti.f32

particles_aos = Particle.field(shape=n_particles)

# Structure of Arrays (SoA)
x_soa = ti.Vector.field(3, dtype=ti.f32, shape=n_particles)
v_soa = ti.Vector.field(3, dtype=ti.f32, shape=n_particles)
m_soa = ti.field(dtype=ti.f32, shape=n_particles)
```

**性能考虑**：
- AoS：空间局部性好，适合访问单个粒子的所有属性
- SoA：向量化友好，适合批量处理同一属性

**自定义布局示例**：
```python
# 块状布局，结合AoS和SoA的优势
block_size = 32
x = ti.Vector.field(3, dtype=ti.f32)
v = ti.Vector.field(3, dtype=ti.f32)

# 二级结构：块内SoA，块间数组
ti.root.dense(ti.i, n_particles//block_size).dense(ti.j, block_size).place(x, v)
```

**内存对齐**：
```python
# 确保缓存行对齐
@ti.kernel
def aligned_access():
    # Taichi自动处理对齐
    for i in range(n_particles):
        # 连续访问，利用缓存
        x[i] += v[i] * dt
```

## 本章小结

物质点法(MPM)作为混合欧拉-拉格朗日方法的代表，成功结合了两种视角的优势。本章深入探讨了MPM的理论基础、算法实现和工程应用：

**核心概念**：
1. **MPM基础**：粒子携带材料信息，网格用于动量更新，避免网格扭曲和数值耗散
2. **MLS-MPM**：通过移动最小二乘简化实现，仅需88行代码，性能提升50%
3. **本构模型**：支持弹性、塑性、流体等多种材料，通过SVD实现屈服准则
4. **数值断裂**：连续损伤力学、相场方法、CPIC等断裂模拟技术
5. **隐式积分**：处理刚性材料和大时间步长，Newton-Raphson迭代求解

**关键公式**：
- 弱形式：$\int_\Omega \rho \mathbf{a} \cdot \mathbf{w} \, dV = -\int_\Omega \boldsymbol{\sigma} : \nabla \mathbf{w} \, dV + \int_\Omega \rho \mathbf{b} \cdot \mathbf{w} \, dV$
- 变形梯度更新：$\mathbf{F}^{n+1} = (\mathbf{I} + \Delta t \mathbf{C}^n) \mathbf{F}^n$
- Neo-Hookean应力：$\mathbf{P} = \mu(\mathbf{F} - \mathbf{F}^{-T}) + \lambda\ln(J)\mathbf{F}^{-T}$
- 塑性返回映射：通过SVD夹持奇异值到屈服面

**实现要点**：
- 形函数选择：二次B样条提供精度-效率平衡
- P2G/G2P传输：APIC保证角动量守恒
- 稀疏数据结构：SPGrid或OpenVDB处理大规模稀疏域
- 性能优化：SoA布局、向量化、多层级网格

## 练习题

### 基础题

**练习8.1**：推导二维情况下二次B样条基函数的显式表达式，并计算其梯度。

*提示*：从一维B样条$N(x)$出发，使用张量积构造二维基函数。

<details>
<summary>答案</summary>

一维二次B样条：
$$N(x) = \begin{cases}
\frac{3}{4} - x^2 & |x| \leq \frac{1}{2} \\
\frac{1}{2}(\frac{3}{2} - |x|)^2 & \frac{1}{2} < |x| \leq \frac{3}{2} \\
0 & |x| > \frac{3}{2}
\end{cases}$$

二维基函数：
$$N_{ij}(x, y) = N(x - x_i)N(y - y_j)$$

梯度：
$$\nabla N_{ij} = \begin{bmatrix}
N'(x - x_i)N(y - y_j) \\
N(x - x_i)N'(y - y_j)
\end{bmatrix}$$

其中$N'(x) = -2x$当$|x| \leq \frac{1}{2}$，$N'(x) = -(\frac{3}{2} - |x|)\text{sign}(x)$当$\frac{1}{2} < |x| \leq \frac{3}{2}$。
</details>

**练习8.2**：证明MLS-MPM中的$\mathbf{C}_p$矩阵确实近似了速度梯度$\nabla\mathbf{v}$。

*提示*：从APIC的推导出发，考虑仿射速度场$\mathbf{v}(\mathbf{x}) = \mathbf{v}_p + \mathbf{C}_p(\mathbf{x} - \mathbf{x}_p)$。

<details>
<summary>答案</summary>

仿射速度场的梯度：
$$\nabla\mathbf{v} = \nabla[\mathbf{v}_p + \mathbf{C}_p(\mathbf{x} - \mathbf{x}_p)] = \mathbf{C}_p$$

在G2P传输中：
$$\mathbf{C}_p = \sum_i \mathbf{v}_i \otimes \nabla N_i(\mathbf{x}_p)$$

这正是速度场在粒子位置的梯度的加权平均，因此$\mathbf{C}_p \approx \nabla\mathbf{v}|_{\mathbf{x}_p}$。
</details>

**练习8.3**：计算Neo-Hookean模型在单轴拉伸下的应力-应变关系。设拉伸比为$\lambda$，其他方向自由变形。

*提示*：利用不可压缩条件$\det(\mathbf{F}) = 1$。

<details>
<summary>答案</summary>

单轴拉伸的变形梯度：
$$\mathbf{F} = \text{diag}(\lambda, \lambda^{-1/2}, \lambda^{-1/2})$$

保证$\det(\mathbf{F}) = \lambda \cdot \lambda^{-1/2} \cdot \lambda^{-1/2} = 1$。

Neo-Hookean应力：
$$P_{11} = \mu(\lambda - \lambda^{-1})$$

工程应力：
$$\sigma = \frac{P_{11}}{\lambda} = \mu(1 - \lambda^{-2})$$
</details>

### 挑战题

**练习8.4**：设计一个自适应粒子细化算法，当检测到大变形时自动增加粒子密度。给出细化准则和粒子分裂策略。

*提示*：考虑基于变形梯度的行列式或最大主伸长比的准则。

<details>
<summary>答案</summary>

细化准则：
1. 体积变化：$|\det(\mathbf{F}) - 1| > \theta_V$
2. 最大主伸长：$\lambda_{\max} > \theta_{\lambda}$
3. 塑性应变：$\varepsilon_p > \theta_p$

分裂策略：
```python
def split_particle(p):
    # 获取主方向
    U, S, V = svd(F[p])
    
    # 沿最大伸长方向分裂
    direction = U[:, 0]
    offset = 0.25 * dx * direction
    
    # 创建子粒子
    for i in [-1, 1]:
        new_p = create_particle()
        x[new_p] = x[p] + i * offset
        v[new_p] = v[p]
        F[new_p] = F[p]
        m[new_p] = m[p] / 2
        V[new_p] = V[p] / 2
    
    # 删除原粒子
    delete_particle(p)
```
</details>

**练习8.5**：推导并实现各向异性材料的MPM本构模型，考虑纤维方向的影响。

*提示*：使用结构张量$\mathbf{M} = \mathbf{a}_0 \otimes \mathbf{a}_0$表示纤维方向。

<details>
<summary>答案</summary>

各向异性超弹性模型：
$$\psi = \psi_{iso}(\mathbf{F}) + \psi_{aniso}(\mathbf{F}, \mathbf{a}_0)$$

其中各向异性部分：
$$\psi_{aniso} = \frac{k_1}{2k_2}[\exp(k_2(I_4 - 1)^2) - 1]$$

$I_4 = \mathbf{a}_0 \cdot (\mathbf{C} \mathbf{a}_0)$是纤维伸长的平方。

应力：
$$\mathbf{P} = \mathbf{P}_{iso} + 2k_1(I_4 - 1)\exp(k_2(I_4 - 1)^2)\mathbf{F}\mathbf{a}_0 \otimes \mathbf{a}_0$$
</details>

**练习8.6**：分析MPM中的能量守恒性质。证明在无外力、无阻尼的情况下，APIC传输保证总角动量守恒。

*提示*：计算系统总角动量$\mathbf{L} = \sum_p m_p \mathbf{x}_p \times \mathbf{v}_p$的时间导数。

<details>
<summary>答案</summary>

系统总角动量：
$$\mathbf{L} = \sum_p m_p \mathbf{x}_p \times \mathbf{v}_p$$

P2G后网格角动量：
$$\mathbf{L}_g = \sum_i m_i \mathbf{x}_i \times \mathbf{v}_i$$

由于APIC传输包含了$\mathbf{C}_p$项：
$$\mathbf{v}_i = \frac{1}{m_i}\sum_p w_{ip}m_p[\mathbf{v}_p + \mathbf{C}_p(\mathbf{x}_i - \mathbf{x}_p)]$$

可以证明：
$$\mathbf{L}_g = \mathbf{L}_p$$

G2P后，由于使用相同的权重和$\mathbf{C}_p$更新，角动量继续守恒。
</details>

**练习8.7**：设计一个MPM-FEM耦合算法，在小变形区域使用FEM，大变形区域使用MPM。

*提示*：考虑过渡区域的处理和动量传递。

<details>
<summary>答案</summary>

耦合策略：
1. **区域划分**：基于变形梯度的行列式判断
2. **过渡区**：同时存在粒子和网格节点
3. **动量交换**：
   - FEM→MPM：在FEM边界生成虚拟粒子
   - MPM→FEM：粒子贡献作为FEM的边界力

算法框架：
```python
def coupled_step():
    # FEM区域
    K_fem = assemble_stiffness_matrix()
    f_fem = assemble_force_vector()
    
    # MPM贡献到FEM边界
    f_mpm_to_fem = compute_mpm_boundary_force()
    f_fem += f_mpm_to_fem
    
    # 求解FEM
    u_fem = solve(K_fem, f_fem)
    
    # FEM速度作为MPM边界条件
    v_fem_boundary = differentiate(u_fem) / dt
    
    # MPM步骤（使用FEM边界速度）
    mpm_step_with_boundary(v_fem_boundary)
```
</details>

**练习8.8**：实现一个基于机器学习的本构模型，使用神经网络替代解析的应力-应变关系。

*提示*：网络输入为$\mathbf{F}$或其不变量，输出为$\mathbf{P}$。

<details>
<summary>答案</summary>

神经网络本构模型：
```python
class NeuralConstitutive(nn.Module):
    def __init__(self):
        super().__init__()
        # 输入：F的不变量 (I1, I2, I3)
        self.net = nn.Sequential(
            nn.Linear(3, 64),
            nn.ReLU(),
            nn.Linear(64, 64),
            nn.ReLU(),
            nn.Linear(64, 9)  # 输出P的9个分量
        )
    
    def forward(self, F):
        # 计算不变量
        I1 = torch.trace(F.T @ F)
        I2 = 0.5 * (I1**2 - torch.trace((F.T @ F)**2))
        I3 = torch.det(F)**2
        
        invariants = torch.stack([I1, I2, I3])
        P_vec = self.net(invariants)
        P = P_vec.reshape(3, 3)
        
        # 确保物理一致性
        P = self.enforce_symmetry(P)
        return P
```

训练数据来自：
1. 解析模型生成
2. 实验测量
3. 高精度仿真
</details>

## 常见陷阱与错误 (Gotchas)

1. **粒子穿越网格边界**：粒子离开计算域时需要正确处理，否则会导致段错误
   - 解决：边界检查和粒子删除机制

2. **变形梯度退化**：$\det(\mathbf{F}) \leq 0$导致应力计算失败
   - 解决：使用可逆有限元技术或SVD夹持

3. **时间步长过大**：CFL条件$\Delta t < \frac{\Delta x}{c}$，其中$c = \sqrt{E/\rho}$
   - 解决：自适应时间步长或隐式积分

4. **能量不守恒**：数值耗散导致能量损失
   - 解决：使用APIC/MLS-MPM，辛积分器

5. **内存爆炸**：粒子数量动态增长导致内存不足
   - 解决：粒子数量上限，自适应粒子管理

6. **并行竞态**：P2G阶段的原子操作性能瓶颈
   - 解决：粒子排序，分块并行

7. **边界条件不一致**：网格边界条件与粒子运动冲突
   - 解决：统一的边界处理策略

8. **数值断裂不稳定**：断裂后应力集中导致数值爆炸
   - 解决：应力释放，损伤正则化

## 最佳实践检查清单

### 算法选择
- [ ] 小变形：考虑使用FEM而非MPM
- [ ] 流固耦合：MLS-MPM统一框架
- [ ] 断裂模拟：CPIC或相场方法
- [ ] 大时间步长：隐式MPM

### 参数设置
- [ ] 粒子密度：每个单元2-4个粒子（2D），4-8个（3D）
- [ ] 形函数阶数：二次B样条（默认），三次（高精度）
- [ ] 时间步长：满足CFL条件，考虑材料刚度
- [ ] 网格分辨率：捕捉最小特征尺寸

### 性能优化
- [ ] 数据布局：SoA用于批量计算，AoS用于随机访问
- [ ] 稀疏结构：SPGrid或OpenVDB处理大规模稀疏域
- [ ] 并行策略：粒子排序减少原子操作冲突
- [ ] 内存管理：预分配，避免动态分配

### 数值稳定性
- [ ] 变形梯度监控：检测并处理退化情况
- [ ] 能量监控：跟踪系统总能量变化
- [ ] 残差检查：隐式求解器的收敛性
- [ ] 边界处理：确保边界条件的一致性

### 验证与调试
- [ ] 单元测试：本构模型、传输算子
- [ ] 收敛性测试：网格细化研究
- [ ] 守恒性检查：质量、动量、角动量
- [ ] 基准测试：与解析解或实验对比
