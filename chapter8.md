# 第八章：混合欧拉-拉格朗日视角（2）：物质点法

物质点法（Material Point Method, MPM）是一种强大的数值方法，结合了拉格朗日粒子和欧拉网格的优势。本章将深入探讨MPM的理论基础、经典算法、现代改进以及在各种材料模拟中的应用。我们将从基础的MPM算法开始，逐步深入到MLS-MPM、本构模型、断裂模拟等高级主题，并介绍Taichi中实现MPM的高级特性。

## 8.1 物质点法(MPM)基础

### 8.1.1 MPM的历史与发展

物质点法最初由Sulsky和Schreyer在1996年提出，作为有限元方法(FEM)的扩展来处理大变形问题。传统FEM在处理极大变形时会遇到网格扭曲问题，而MPM通过使用拉格朗日粒子（物质点）携带材料信息，欧拉背景网格进行动量方程求解，巧妙地避免了这个问题。

2013年，Stomakhin等人将MPM引入计算机图形学领域，首次实现了雪的真实感模拟。这项工作展示了MPM在处理相变、断裂等复杂物理现象上的独特优势。随后，MPM在图形学领域得到快速发展，被广泛应用于各种材料的模拟，包括沙子、泥土、泡沫、布料等。

MPM的核心思想是将连续介质离散为一系列携带质量、动量、应力等物理量的粒子，而背景网格仅用于计算内力和更新动量。这种双重表示方式使得MPM既保留了拉格朗日方法追踪材料历史的能力，又具备了欧拉方法处理碰撞和自碰撞的便利性。

### 8.1.2 与FEM的关系

MPM可以视为无网格Galerkin方法的一种，特别是属于无单元Galerkin（Element-Free Galerkin, EFG）方法家族。从数学角度看，MPM和FEM都基于连续介质力学的弱形式：

$$\int_\Omega \rho \mathbf{a} \cdot \mathbf{w} \, dV = -\int_\Omega \boldsymbol{\sigma} : \nabla \mathbf{w} \, dV + \int_\Omega \rho \mathbf{b} \cdot \mathbf{w} \, dV + \int_{\partial\Omega_t} \mathbf{t} \cdot \mathbf{w} \, dA$$

其中$\mathbf{w}$是测试函数，$\mathbf{a}$是加速度，$\boldsymbol{\sigma}$是柯西应力张量，$\mathbf{b}$是体力，$\mathbf{t}$是表面力。

MPM与FEM的主要区别在于积分方式：
- FEM使用高斯积分点，这些点的位置在单元内固定
- MPM使用物质点作为积分点，这些点可以自由移动
- MPM中每个粒子可以视为一个积分点，其积分权重对应于粒子的体积

这种关系意味着MPM继承了FEM的许多理论基础，包括收敛性分析、误差估计等，同时又具有处理大变形的灵活性。

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

MPM的空间离散化采用混合方法：

1. **粒子表示**：材料被离散为$N_p$个粒子，每个粒子携带：
   - 位置 $\mathbf{x}_p$
   - 速度 $\mathbf{v}_p$
   - 质量 $m_p$
   - 体积 $V_p$
   - 变形梯度 $\mathbf{F}_p$
   - 应力 $\boldsymbol{\sigma}_p$

2. **背景网格**：使用规则的背景网格（通常是笛卡尔网格），网格节点$i$具有：
   - 位置 $\mathbf{x}_i$
   - 速度 $\mathbf{v}_i$
   - 质量 $m_i$
   - 动量 $\mathbf{p}_i = m_i \mathbf{v}_i$

3. **形函数**：使用形函数$N_i(\mathbf{x})$在粒子和网格之间插值，常用的是B样条基函数：
   - 线性（1次）：$N(x) = 1 - |x|$, $|x| \leq 1$
   - 二次：$N(x) = \begin{cases} \frac{3}{4} - x^2 & |x| \leq \frac{1}{2} \\ \frac{1}{2}(|x| - \frac{3}{2})^2 & \frac{1}{2} < |x| \leq \frac{3}{2} \end{cases}$
   - 三次：使用三次B样条，支持域更大但更光滑

空间离散化的关键是粒子-网格传输（P2G）和网格-粒子传输（G2P）操作，这些将在后续章节详细讨论。

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

形函数的选择对MPM的精度和稳定性有重要影响：

**线性B样条**：
- 计算效率高，每个粒子影响$2^d$个节点（$d$是维度）
- 可能产生"格子噪声"（grid crossing error）
- 适合快速原型开发

**二次B样条**：
- 更光滑，减少格子噪声
- 每个粒子影响$3^d$个节点
- 计算成本适中，是常用选择

**三次B样条**：
- 最光滑，几乎消除格子噪声
- 每个粒子影响$4^d$个节点
- 计算成本较高，用于高精度要求场景

形函数及其梯度的计算：
$$N_i(\mathbf{x}_p) = \prod_{d=1}^{dim} N^{1D}\left(\frac{x_p^d - x_i^d}{\Delta x}\right)$$

$$\nabla N_i(\mathbf{x}_p) = \nabla \left[ \prod_{d} N^{1D}_d \right] = \sum_{d} \frac{\partial N^{1D}_d}{\partial x^d} \prod_{k \neq d} N^{1D}_k$$

### 8.2.3 积分点与背景网格

MPM中的数值积分使用单点积分规则，每个粒子作为一个积分点：

$$\int_\Omega f(\mathbf{x}) \, dV \approx \sum_p V_p^0 f(\mathbf{x}_p)$$

其中$V_p^0$是粒子的初始体积。

**背景网格类型**：
1. **均匀笛卡尔网格**：最常用，实现简单，适合GPU并行
2. **SPGrid**：稀疏分页网格，利用虚拟内存管理稀疏数据
3. **OpenVDB**：层次化B+树结构，适合极度稀疏的大规模模拟

**粒子分布策略**：
- 每个网格单元通常放置$2^d$到$4^d$个粒子
- 粒子初始位置可以是规则分布或随机扰动
- 粒子密度影响数值精度和计算成本

### 8.2.4 边界条件处理

MPM中的边界条件在网格级别施加，主要有三种类型：

**粘性边界（Sticky）**：
粒子在边界上速度完全为零
$$\mathbf{v}_i = 0 \quad \text{if } \mathbf{x}_i \in \partial\Omega_{sticky}$$

**滑动边界（Slip）**：
只约束法向速度，允许切向滑动
$$\mathbf{v}_i \cdot \mathbf{n} = 0, \quad \mathbf{v}_i = \mathbf{v}_i - (\mathbf{v}_i \cdot \mathbf{n})\mathbf{n}$$

**分离边界（Separate）**：
单向约束，只在穿透时起作用
$$\text{if } \mathbf{v}_i \cdot \mathbf{n} < 0: \quad \mathbf{v}_i = \mathbf{v}_i - (\mathbf{v}_i \cdot \mathbf{n})\mathbf{n}$$

实现时，通常在网格动量更新后立即施加边界条件：
```python
# 网格动量更新后
if is_boundary(i, j):
    if boundary_type == STICKY:
        v[i, j] = 0
    elif boundary_type == SLIP:
        normal = get_boundary_normal(i, j)
        v[i, j] -= dot(v[i, j], normal) * normal
```

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
