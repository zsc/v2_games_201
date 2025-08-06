# 第四章：欧拉视角（1）

本章我们将从拉格朗日视角转向欧拉视角，探讨基于网格的流体仿真方法。欧拉方法在处理大变形流动、自由表面追踪和不可压缩性约束方面具有独特优势。我们将深入理解Navier-Stokes方程的数值求解，掌握压力投影法的核心思想，并学习如何处理各种边界条件。

## 4.1 欧拉vs拉格朗日描述

### 4.1.1 物质导数与对流项

在流体力学中，我们需要追踪流体质点的物理量变化。欧拉视角下，我们在固定的空间位置观察流体，这导致了物质导数的出现：

$$\frac{D}{Dt} = \frac{\partial}{\partial t} + \mathbf{u} \cdot \nabla$$

其中第一项 $\frac{\partial}{\partial t}$ 是局部时间导数，表示固定位置处物理量的时间变化率；第二项 $\mathbf{u} \cdot \nabla$ 是对流项，表示由于流体运动导致的物理量变化。

**物质导数的物理意义**：考虑一个流体质点沿轨迹 $\mathbf{x}(t)$ 运动，其携带的物理量 $\phi$ 的变化率为：
$$\frac{d\phi}{dt} = \frac{\partial \phi}{\partial t} + \frac{\partial \phi}{\partial x_i}\frac{dx_i}{dt} = \frac{\partial \phi}{\partial t} + u_i\frac{\partial \phi}{\partial x_i}$$

这里使用了爱因斯坦求和约定。物质导数连接了拉格朗日描述（跟随质点）和欧拉描述（固定空间点）。

例如，对于速度场的物质导数：
$$\frac{D\mathbf{u}}{Dt} = \frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u} \cdot \nabla)\mathbf{u}$$

展开对流项的分量形式（以2D为例）：
$$(\mathbf{u} \cdot \nabla)\mathbf{u} = \begin{pmatrix} u\frac{\partial u}{\partial x} + v\frac{\partial u}{\partial y} \\ u\frac{\partial v}{\partial x} + v\frac{\partial v}{\partial y} \end{pmatrix}$$

这个对流项 $(\mathbf{u} \cdot \nabla)\mathbf{u}$ 是非线性的，也是Navier-Stokes方程数值求解的主要挑战之一。它导致了湍流等复杂现象的产生。

**对流项的另一种理解**：使用向量恒等式，对流项可以改写为：
$$(\mathbf{u} \cdot \nabla)\mathbf{u} = \nabla(\frac{|\mathbf{u}|^2}{2}) - \mathbf{u} \times (\nabla \times \mathbf{u})$$

第一项是动能梯度，第二项包含涡量 $\boldsymbol{\omega} = \nabla \times \mathbf{u}$，体现了旋转效应对动量输送的影响。

### 4.1.2 不可压Navier-Stokes方程

不可压缩流体的运动由Navier-Stokes方程描述：

$$\rho \frac{D\mathbf{u}}{Dt} = -\nabla p + \mu \nabla^2 \mathbf{u} + \rho \mathbf{g}$$

配合不可压缩条件：
$$\nabla \cdot \mathbf{u} = 0$$

其中：
- $\rho$ 是流体密度（对于不可压流体为常数）
- $p$ 是压力（实际上是压力除以密度，有压力势的量纲）
- $\mu$ 是动力粘度
- $\nu = \mu/\rho$ 是运动粘度
- $\mathbf{g}$ 是重力加速度

展开物质导数后，动量方程变为：
$$\rho \left(\frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u} \cdot \nabla)\mathbf{u}\right) = -\nabla p + \mu \nabla^2 \mathbf{u} + \rho \mathbf{g}$$

**无量纲化与Reynolds数**：通过特征尺度 $L$、特征速度 $U$ 和特征时间 $L/U$ 无量纲化，可得：
$$\frac{\partial \mathbf{u}^*}{\partial t^*} + (\mathbf{u}^* \cdot \nabla^*)\mathbf{u}^* = -\nabla^* p^* + \frac{1}{Re} \nabla^{*2} \mathbf{u}^* + \frac{1}{Fr^2}\hat{\mathbf{g}}$$

其中 Reynolds数 $Re = \frac{UL}{\nu}$ 表征惯性力与粘性力之比，Froude数 $Fr = \frac{U}{\sqrt{gL}}$ 表征惯性力与重力之比。

**压力的作用**：在不可压缩流中，压力不是状态变量，而是一个拉格朗日乘子，用于强制满足不可压缩约束。压力瞬时调整以维持 $\nabla \cdot \mathbf{u} = 0$。

**涡量形式**：取动量方程的旋度，可得涡量输送方程：
$$\frac{D\boldsymbol{\omega}}{Dt} = (\boldsymbol{\omega} \cdot \nabla)\mathbf{u} + \nu \nabla^2 \boldsymbol{\omega}$$

在2D情况下，涡量拉伸项 $(\boldsymbol{\omega} \cdot \nabla)\mathbf{u}$ 消失，涡量仅通过对流和扩散演化。

### 4.1.3 算子分裂方法

直接求解耦合的速度-压力系统非常困难。Chorin和Temam提出的算子分裂方法将Navier-Stokes方程分解为几个子步骤：

1. **对流步**：求解 $\frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u} \cdot \nabla)\mathbf{u} = 0$
2. **外力步**：添加重力和粘性力 $\frac{\partial \mathbf{u}}{\partial t} = \nu \nabla^2 \mathbf{u} + \mathbf{g}$
3. **投影步**：通过压力投影使速度场满足不可压缩条件

**数学基础**：算子分裂基于Trotter-Lie公式。对于方程 $\frac{\partial u}{\partial t} = (A + B)u$，其解可以近似为：
$$u(t + \Delta t) \approx e^{\Delta t B} e^{\Delta t A} u(t) + O(\Delta t^2)$$

这是一阶分裂。二阶Strang分裂为：
$$u(t + \Delta t) \approx e^{\Delta t B/2} e^{\Delta t A} e^{\Delta t B/2} u(t) + O(\Delta t^3)$$

**具体算法步骤**：

1. **预测步**（忽略压力）：
   $$\mathbf{u}^* = \mathbf{u}^n + \Delta t \left[ -(\mathbf{u}^n \cdot \nabla)\mathbf{u}^n + \nu \nabla^2 \mathbf{u}^n + \mathbf{g} \right]$$

2. **压力求解**：
   $$\nabla^2 p^{n+1} = \frac{\rho}{\Delta t} \nabla \cdot \mathbf{u}^*$$

3. **速度修正**：
   $$\mathbf{u}^{n+1} = \mathbf{u}^* - \frac{\Delta t}{\rho} \nabla p^{n+1}$$

这种分裂使得每个子问题都可以高效求解。在Taichi中，典型的实现框架如下：

```python
@ti.kernel
def step():
    advect_velocity()      # 对流步
    apply_external_forces() # 外力步  
    project_velocity()      # 投影步
    advect_other_quantities() # 输送其他物理量
```

**分裂误差分析**：算子分裂引入的误差主要来自于忽略了各项之间的耦合。例如，压力梯度实际上会影响对流，但在分裂方法中这种影响被延迟到投影步。这种误差通常是 $O(\Delta t)$，但可以通过高阶分裂方案减小。

### 4.1.4 稳定性与CFL条件

显式时间积分的稳定性受CFL（Courant-Friedrichs-Lewy）条件限制：

$$\Delta t \leq C \frac{\Delta x}{|\mathbf{u}|_{\max}}$$

其中 $C$ 是CFL数，通常取0.5-1.0。这个条件确保在一个时间步内，流体粒子移动的距离不超过一个网格单元。

**CFL条件的推导**：考虑一维线性对流方程 $\frac{\partial u}{\partial t} + c\frac{\partial u}{\partial x} = 0$，使用前向时间、中心空间差分：
$$\frac{u_i^{n+1} - u_i^n}{\Delta t} + c\frac{u_{i+1}^n - u_{i-1}^n}{2\Delta x} = 0$$

通过von Neumann稳定性分析，代入 $u_i^n = \hat{u}^n e^{ik_xi\Delta x}$，可得放大因子：
$$G = 1 - i\nu \sin(k\Delta x)$$

其中 $\nu = c\Delta t/\Delta x$ 是CFL数。稳定性要求 $|G| \leq 1$，但对于中心差分这永远无法满足。使用迎风格式可得稳定条件 $\nu \leq 1$。

对于粘性项，如果使用显式积分，还需要满足扩散稳定性条件：
$$\Delta t \leq \frac{(\Delta x)^2}{2d\nu}$$

其中 $d$ 是空间维度。由于这个条件对小网格非常严格（$\Delta t \propto (\Delta x)^2$），粘性项通常使用隐式或半隐式方法处理。

**多物理场的综合稳定性条件**：

1. **对流CFL**：$\Delta t_{adv} = C_{adv} \frac{\Delta x}{|\mathbf{u}|_{\max}}$

2. **粘性稳定性**：$\Delta t_{visc} = C_{visc} \frac{(\Delta x)^2}{\nu}$

3. **表面张力**（如果存在）：$\Delta t_{\sigma} = C_{\sigma} \sqrt{\frac{\rho (\Delta x)^3}{\sigma}}$

4. **重力波**（自由表面）：$\Delta t_{grav} = C_{grav} \sqrt{\frac{\Delta x}{g}}$

最终时间步长：$\Delta t = \min(\Delta t_{adv}, \Delta t_{visc}, \Delta t_{\sigma}, \Delta t_{grav})$

其中 $C_{adv} \approx 0.5$, $C_{visc} \approx 0.25$, $C_{\sigma} \approx 0.5$, $C_{grav} \approx 0.5$ 是安全系数。

## 4.2 网格类型与离散化

### 4.2.1 同位网格vs交错网格(Staggered Grid)

**同位网格（Collocated Grid）**：
- 所有物理量（速度分量、压力等）存储在网格单元中心
- 实现简单，但可能产生棋盘格压力振荡
- 需要特殊处理来避免数值不稳定

**棋盘格问题的根源**：在同位网格上，压力梯度使用中心差分：
$$\left(\frac{\partial p}{\partial x}\right)_i = \frac{p_{i+1} - p_{i-1}}{2\Delta x}$$

这个离散算子无法感知棋盘格模式 $p_i = (-1)^i$，因为 $p_{i+1} - p_{i-1} = 0$。这导致压力泊松方程的零空间包含非物理的高频模式。

**交错网格（MAC Grid）**：
- 压力存储在单元中心
- 速度分量存储在对应的单元面中心
- $u$ 存储在 $(i+1/2, j)$ 位置
- $v$ 存储在 $(i, j+1/2)$ 位置

MAC网格的布局：
```
      v(i,j+1/2)
         |
u(i-1/2,j)---p(i,j)---u(i+1/2,j)
         |
      v(i,j-1/2)
```

**MAC网格的数学优势**：
1. **紧凑的梯度-散度对**：散度和梯度算子互为负转置：
   $$\langle \nabla \cdot \mathbf{u}, p \rangle = -\langle \mathbf{u}, \nabla p \rangle$$
   这保证了离散系统的能量守恒性。

2. **自然的通量计算**：速度直接定义在单元界面上，便于计算质量通量。

3. **最优的inf-sup条件**：MAC离散满足离散inf-sup（LBB）条件，保证了压力的唯一性。

**实现细节**：在Taichi中，MAC网格的索引约定：
```python
# 压力和标量场：单元中心
pressure = ti.field(float, shape=(nx, ny))
# u速度：垂直面中心
u = ti.field(float, shape=(nx+1, ny))
# v速度：水平面中心  
v = ti.field(float, shape=(nx, ny+1))
```

### 4.2.2 MAC网格的优势

1. **自然满足散度定理**：在单元边界上直接计算通量
2. **避免压力振荡**：压力和速度的耦合更紧密
3. **边界条件处理简单**：速度分量直接位于边界上
4. **守恒性好**：离散后仍保持质量、动量守恒

离散散度算子在MAC网格上特别简洁：
$$(\nabla \cdot \mathbf{u})_{i,j} = \frac{u_{i+1/2,j} - u_{i-1/2,j}}{\Delta x} + \frac{v_{i,j+1/2} - v_{i,j-1/2}}{\Delta y}$$

### 4.2.3 双线性插值

由于速度分量存储在不同位置，经常需要插值获取任意位置的速度。对于双线性插值：

$$f(x,y) = f_{00}(1-\alpha)(1-\beta) + f_{10}\alpha(1-\beta) + f_{01}(1-\alpha)\beta + f_{11}\alpha\beta$$

其中 $\alpha = (x - x_0)/\Delta x$，$\beta = (y - y_0)/\Delta y$ 是局部坐标。

在Taichi中实现双线性插值：
```python
@ti.func
def bilinear_interp(field: ti.template(), x: float, y: float) -> float:
    i, j = int(x), int(y)
    fx, fy = x - i, y - j
    return (field[i,j] * (1-fx) * (1-fy) +
            field[i+1,j] * fx * (1-fy) +
            field[i,j+1] * (1-fx) * fy +
            field[i+1,j+1] * fx * fy)
```

### 4.2.4 散度与梯度算子离散

**散度算子**（用于不可压缩约束）：
$$(\nabla \cdot \mathbf{u})_{i,j} = \frac{u_{i+1/2,j} - u_{i-1/2,j}}{\Delta x} + \frac{v_{i,j+1/2} - v_{i,j-1/2}}{\Delta y}$$

**梯度算子**（用于压力梯度）：
$$(\nabla p)_x|_{i+1/2,j} = \frac{p_{i+1,j} - p_{i,j}}{\Delta x}$$
$$(\nabla p)_y|_{i,j+1/2} = \frac{p_{i,j+1} - p_{i,j}}{\Delta y}$$

**Laplace算子**（用于粘性项和压力泊松方程）：
$$(\nabla^2 p)_{i,j} = \frac{p_{i+1,j} - 2p_{i,j} + p_{i-1,j}}{(\Delta x)^2} + \frac{p_{i,j+1} - 2p_{i,j} + p_{i,j-1}}{(\Delta y)^2}$$

## 4.3 半拉格朗日输送(Semi-Lagrangian Advection)

### 4.3.1 反向粒子追踪

半拉格朗日方法结合了拉格朗日和欧拉方法的优点。基本思想是：对于网格点 $\mathbf{x}$，追踪到达该点的流体粒子在上一时刻的位置 $\mathbf{x}_{prev}$：

$$\mathbf{x}_{prev} = \mathbf{x} - \Delta t \cdot \mathbf{u}(\mathbf{x})$$

然后通过插值获取该位置的物理量作为新时刻的值：
$$q^{n+1}(\mathbf{x}) = q^n(\mathbf{x}_{prev})$$

这种方法无条件稳定，允许使用大时间步长，但会引入数值耗散。

### 4.3.2 速度插值策略

简单的一阶精度追踪使用当前位置的速度，但可以通过Runge-Kutta方法提高精度：

**RK2（中点法）**：
```python
@ti.func
def backtrace_rk2(x: ti.math.vec2, dt: float) -> ti.math.vec2:
    v1 = sample_velocity(x)
    x_mid = x - 0.5 * dt * v1
    v2 = sample_velocity(x_mid)
    return x - dt * v2
```

**RK3**：
```python
@ti.func  
def backtrace_rk3(x: ti.math.vec2, dt: float) -> ti.math.vec2:
    v1 = sample_velocity(x)
    x1 = x - 0.5 * dt * v1
    v2 = sample_velocity(x1)
    x2 = x - 0.75 * dt * v2
    v3 = sample_velocity(x2)
    return x - dt * (2.0/9.0 * v1 + 1.0/3.0 * v2 + 4.0/9.0 * v3)
```

### 4.3.3 数值耗散问题

半拉格朗日方法的主要缺点是数值耗散，表现为：
- 涡旋快速衰减
- 细节特征模糊
- 能量不守恒

数值耗散的根源是插值过程中的平滑效应。每次插值相当于应用了一个低通滤波器，多次输送后高频信息逐渐丢失。

### 4.3.4 单调性保持

为了避免产生非物理的新极值，可以使用限制器：

```python
@ti.func
def clamp_extrema(q_new: float, x: ti.math.vec2) -> float:
    # 获取周围网格点的最大最小值
    q_min = ti.math.inf
    q_max = -ti.math.inf
    for di in range(-1, 2):
        for dj in range(-1, 2):
            q_neighbor = q_field[int(x[0])+di, int(x[1])+dj]
            q_min = min(q_min, q_neighbor)
            q_max = max(q_max, q_neighbor)
    return clamp(q_new, q_min, q_max)
```

## 4.4 高阶输送格式

### 4.4.1 MacCormack方法

MacCormack方法通过前向和后向输送估计误差并补偿：

1. 前向输送：$\phi^* = \text{SemiLagrangian}(\phi^n, \Delta t)$
2. 后向输送：$\phi^{**} = \text{SemiLagrangian}(\phi^*, -\Delta t)$  
3. 误差估计：$e = \frac{1}{2}(\phi^n - \phi^{**})$
4. 修正结果：$\phi^{n+1} = \phi^* + e$

实现时需要加入限制器防止过冲：
```python
@ti.kernel
def maccormack_advect(phi: ti.template(), phi_new: ti.template(), dt: float):
    # 前向输送
    for i, j in phi:
        x = ti.Vector([i, j]) + 0.5
        x_prev = backtrace(x, dt)
        phi_star[i,j] = bilinear_sample(phi, x_prev)
    
    # 后向输送  
    for i, j in phi:
        x = ti.Vector([i, j]) + 0.5
        x_next = backtrace(x, -dt)
        phi_nn[i,j] = bilinear_sample(phi_star, x_next)
    
    # 误差补偿
    for i, j in phi:
        phi_new[i,j] = phi_star[i,j] + 0.5 * (phi[i,j] - phi_nn[i,j])
        phi_new[i,j] = clamp_extrema(phi_new[i,j], i, j)
```

### 4.4.2 BFECC方法

BFECC (Back and Forth Error Compensation and Correction) 是另一种减少数值耗散的方法：

1. 前向输送：$\phi^{n+1/2} = A(\phi^n)$
2. 后向输送：$\phi^{n*} = A^{-1}(\phi^{n+1/2})$
3. 误差估计：$e = \frac{1}{2}(\phi^n - \phi^{n*})$
4. 修正输入：$\tilde{\phi}^n = \phi^n + e$
5. 最终输送：$\phi^{n+1} = A(\tilde{\phi}^n)$

其中 $A$ 表示输送算子，$A^{-1}$ 表示反向输送。

### 4.4.3 Runge-Kutta积分

对于速度场的高阶积分，常用RK4方法：

```python
@ti.func
def rk4_backtrace(x: ti.math.vec2, dt: float) -> ti.math.vec2:
    k1 = sample_velocity(x)
    k2 = sample_velocity(x - 0.5 * dt * k1)
    k3 = sample_velocity(x - 0.5 * dt * k2)
    k4 = sample_velocity(x - dt * k3)
    return x - dt * (k1 + 2*k2 + 2*k3 + k4) / 6.0
```

### 4.4.4 限制器与夹持(Clamping)

为了保持物理量的界限和避免振荡，需要使用各种限制器：

**MinMod限制器**：
$$\text{minmod}(a, b) = \begin{cases}
a & \text{if } |a| < |b| \text{ and } ab > 0 \\
b & \text{if } |b| < |a| \text{ and } ab > 0 \\
0 & \text{otherwise}
\end{cases}$$

**Superbee限制器**：
$$\text{superbee}(a, b) = \text{maxmod}(\text{minmod}(2a, b), \text{minmod}(a, 2b))$$

这些限制器确保高阶方法不会产生新的极值，保持解的物理性。

## 4.5 Chorin式压力投影

### 4.5.1 Helmholtz-Hodge分解

Helmholtz-Hodge定理指出，任意向量场可以唯一分解为无散场和无旋场之和：

$$\mathbf{u} = \mathbf{u}_{div-free} + \nabla \phi$$

其中 $\nabla \cdot \mathbf{u}_{div-free} = 0$。这是压力投影法的理论基础。

### 4.5.2 压力泊松方程推导

从动量方程出发：
$$\frac{\mathbf{u}^{n+1} - \mathbf{u}^*}{\Delta t} = -\frac{1}{\rho}\nabla p$$

其中 $\mathbf{u}^*$ 是不满足不可压缩条件的中间速度。

对两边取散度，并应用不可压缩条件 $\nabla \cdot \mathbf{u}^{n+1} = 0$：
$$\nabla \cdot \mathbf{u}^{n+1} - \nabla \cdot \mathbf{u}^* = -\frac{\Delta t}{\rho} \nabla^2 p$$

因此得到压力泊松方程：
$$\nabla^2 p = \frac{\rho}{\Delta t} \nabla \cdot \mathbf{u}^*$$

### 4.5.3 离散化与边界条件

在MAC网格上，压力泊松方程的离散形式为：

$$\frac{p_{i+1,j} - 2p_{i,j} + p_{i-1,j}}{(\Delta x)^2} + \frac{p_{i,j+1} - 2p_{i,j} + p_{i,j-1}}{(\Delta y)^2} = \frac{\rho}{\Delta t} \cdot divergence_{i,j}$$

**边界条件**：
- **固体边界**：法向压力梯度为零（Neumann条件）$\frac{\partial p}{\partial n} = 0$
- **自由表面**：压力为大气压（Dirichlet条件）$p = 0$
- **周期边界**：压力和梯度都周期延拓

离散系统可以写成矩阵形式 $\mathbf{A}\mathbf{p} = \mathbf{b}$，其中 $\mathbf{A}$ 是离散Laplace算子。

### 4.5.4 速度修正步

求解压力后，通过压力梯度修正速度场：

$$\mathbf{u}^{n+1} = \mathbf{u}^* - \frac{\Delta t}{\rho} \nabla p$$

在MAC网格上的具体形式：
$$u_{i+1/2,j}^{n+1} = u_{i+1/2,j}^* - \frac{\Delta t}{\rho} \frac{p_{i+1,j} - p_{i,j}}{\Delta x}$$
$$v_{i,j+1/2}^{n+1} = v_{i,j+1/2}^* - \frac{\Delta t}{\rho} \frac{p_{i,j+1} - p_{i,j}}{\Delta y}$$

投影后的速度场自动满足不可压缩条件。

## 4.6 固体边界条件

### 4.6.1 无滑移条件

对于粘性流体，在固体壁面上需要满足无滑移条件：
- 法向速度：$\mathbf{u} \cdot \mathbf{n} = 0$
- 切向速度：$\mathbf{u} \cdot \mathbf{t} = 0$

在MAC网格上，直接设置边界上的速度分量：
```python
@ti.kernel
def enforce_boundary():
    # 左右边界
    for j in range(ny):
        u[0, j] = 0.0      # 左边界
        u[nx, j] = 0.0     # 右边界
    # 上下边界  
    for i in range(nx):
        v[i, 0] = 0.0      # 下边界
        v[i, ny] = 0.0     # 上边界
```

### 4.6.2 自由滑移条件

自由滑移条件只约束法向速度，允许切向滑动：
- 法向速度：$\mathbf{u} \cdot \mathbf{n} = 0$
- 切向应力：$\tau \cdot \mathbf{t} = 0$

实现时，设置法向速度为零，但不修改切向速度。

### 4.6.3 浸入边界法

对于复杂几何形状，浸入边界法（Immersed Boundary Method）将边界力作为体积力添加到动量方程中：

$$\rho \frac{D\mathbf{u}}{Dt} = -\nabla p + \mu \nabla^2 \mathbf{u} + \mathbf{f}_{IB}$$

其中 $\mathbf{f}_{IB}$ 是边界力，通过拉格朗日乘子或罚函数方法计算。

### 4.6.4 切割单元法

切割单元法（Cut-cell Method）精确处理与网格不对齐的边界：

1. 计算每个单元被固体占据的体积分数
2. 修改离散算子以考虑部分单元
3. 在切割单元上应用特殊的插值和梯度计算

这种方法可以达到二阶精度，但实现相对复杂。

## 4.7 自由表面边界条件

### 4.7.1 运动学条件

自由表面必须满足运动学条件：表面上的粒子始终保持在表面上。用等势面函数 $\phi$ 表示界面，其满足：

$$\frac{D\phi}{Dt} = \frac{\partial \phi}{\partial t} + \mathbf{u} \cdot \nabla \phi = 0$$

这是一个Hamilton-Jacobi型方程，描述了界面的输送。

### 4.7.2 动力学条件  

在自由表面上，应力必须连续。忽略空气的影响，边界条件为：

$$p = p_{atm} - \sigma \kappa$$

其中：
- $p_{atm}$ 是大气压（通常设为0）
- $\sigma$ 是表面张力系数
- $\kappa = \nabla \cdot \mathbf{n}$ 是界面曲率

### 4.7.3 表面张力处理

**连续表面力（CSF）模型**将表面张力转换为体积力：

$$\mathbf{F}_{st} = \sigma \kappa \delta(\phi) \nabla \phi$$

其中 $\delta(\phi)$ 是Dirac delta函数的光滑近似。在实际计算中：

$$\kappa = \nabla \cdot \left(\frac{\nabla \phi}{|\nabla \phi|}\right)$$

为了数值稳定性，通常使用光滑的delta函数：
$$\delta_{\epsilon}(\phi) = \begin{cases}
\frac{1}{2\epsilon}\left(1 + \cos\left(\frac{\pi\phi}{\epsilon}\right)\right) & |\phi| < \epsilon \\
0 & \text{otherwise}
\end{cases}$$

### 4.7.4 Ghost Fluid方法

Ghost Fluid方法在界面两侧使用虚拟值来处理不连续：

1. 在界面附近定义ghost cells
2. 根据界面条件外推物理量到ghost cells
3. 使用标准离散格式，自动处理界面跳跃条件

例如，对于压力的ghost值：
```python
@ti.func
def extrapolate_pressure(i: int, j: int) -> float:
    if phi[i,j] * phi[i+1,j] < 0:  # 界面穿过单元
        # 线性外推压力值
        theta = phi[i,j] / (phi[i,j] - phi[i+1,j])
        return (1-theta) * p[i,j] + theta * p_air
    return p[i,j]
```

## 本章小结

本章介绍了欧拉视角下的流体仿真核心概念：

1. **物质导数**：$\frac{D}{Dt} = \frac{\partial}{\partial t} + \mathbf{u} \cdot \nabla$ 连接了欧拉和拉格朗日描述
2. **算子分裂**：将复杂的Navier-Stokes方程分解为对流、外力和投影三步
3. **MAC网格**：交错网格自然满足离散散度定理，避免压力振荡
4. **半拉格朗日输送**：无条件稳定但有数值耗散，需要高阶修正
5. **压力投影**：通过求解泊松方程 $\nabla^2 p = \frac{\rho}{\Delta t}\nabla \cdot \mathbf{u}^*$ 实现不可压缩性
6. **边界条件**：固体边界的无滑移/自由滑移，自由表面的运动学/动力学条件

掌握这些概念是实现稳定高效的流体求解器的基础。

## 练习题

### 基础题

**练习4.1** 推导二维情况下的物质导数展开式。对于标量场 $\phi(x,y,t)$，证明：
$$\frac{D\phi}{Dt} = \frac{\partial \phi}{\partial t} + u\frac{\partial \phi}{\partial x} + v\frac{\partial \phi}{\partial y}$$

<details>
<summary>提示</summary>
考虑流体质点的轨迹 $\mathbf{x}(t)$，使用链式法则。
</details>

<details>
<summary>答案</summary>
沿着流体质点的轨迹 $\mathbf{x}(t) = (x(t), y(t))$，有 $\frac{dx}{dt} = u$，$\frac{dy}{dt} = v$。

应用链式法则：
$$\frac{D\phi}{Dt} = \frac{d}{dt}\phi(x(t), y(t), t) = \frac{\partial \phi}{\partial t} + \frac{\partial \phi}{\partial x}\frac{dx}{dt} + \frac{\partial \phi}{\partial y}\frac{dy}{dt}$$

代入速度分量：
$$\frac{D\phi}{Dt} = \frac{\partial \phi}{\partial t} + u\frac{\partial \phi}{\partial x} + v\frac{\partial \phi}{\partial y} = \frac{\partial \phi}{\partial t} + \mathbf{u} \cdot \nabla\phi$$
</details>

**练习4.2** 在MAC网格上，给定速度场 $u_{i+1/2,j}$ 和 $v_{i,j+1/2}$，写出计算单元 $(i,j)$ 处散度的公式。如果 $\Delta x = \Delta y = h$，散度为零意味着什么？

<details>
<summary>提示</summary>
考虑流入和流出单元的通量平衡。
</details>

<details>
<summary>答案</summary>
散度公式：
$$(\nabla \cdot \mathbf{u})_{i,j} = \frac{u_{i+1/2,j} - u_{i-1/2,j}}{\Delta x} + \frac{v_{i,j+1/2} - v_{i,j-1/2}}{\Delta y}$$

当 $\Delta x = \Delta y = h$ 时：
$$(\nabla \cdot \mathbf{u})_{i,j} = \frac{1}{h}[(u_{i+1/2,j} - u_{i-1/2,j}) + (v_{i,j+1/2} - v_{i,j-1/2})]$$

散度为零意味着流入单元的质量通量等于流出的质量通量，满足质量守恒。
</details>

**练习4.3** 证明压力投影步骤后的速度场满足不可压缩条件。给定 $\mathbf{u}^{n+1} = \mathbf{u}^* - \frac{\Delta t}{\rho}\nabla p$，其中 $p$ 满足 $\nabla^2 p = \frac{\rho}{\Delta t}\nabla \cdot \mathbf{u}^*$。

<details>
<summary>提示</summary>
对速度更新公式两边取散度。
</details>

<details>
<summary>答案</summary>
对速度更新公式取散度：
$$\nabla \cdot \mathbf{u}^{n+1} = \nabla \cdot \mathbf{u}^* - \frac{\Delta t}{\rho}\nabla \cdot (\nabla p)$$

注意到 $\nabla \cdot (\nabla p) = \nabla^2 p$，代入压力泊松方程：
$$\nabla \cdot \mathbf{u}^{n+1} = \nabla \cdot \mathbf{u}^* - \frac{\Delta t}{\rho} \cdot \frac{\rho}{\Delta t}\nabla \cdot \mathbf{u}^*$$

$$\nabla \cdot \mathbf{u}^{n+1} = \nabla \cdot \mathbf{u}^* - \nabla \cdot \mathbf{u}^* = 0$$

因此投影后的速度场自动满足不可压缩条件。
</details>

### 挑战题

**练习4.4** 分析半拉格朗日方法的数值耗散。考虑一维线性对流方程 $\frac{\partial u}{\partial t} + c\frac{\partial u}{\partial x} = 0$，使用线性插值，证明数值解引入了人工扩散项。

<details>
<summary>提示</summary>
将插值误差用Taylor级数展开，分析截断误差的主导项。
</details>

<details>
<summary>答案</summary>
设网格间距为 $h$，时间步长为 $\Delta t$，CFL数 $\nu = c\Delta t/h$。

半拉格朗日更新：$u_i^{n+1} = (1-\alpha)u_{i-1}^n + \alpha u_i^n$，其中 $\alpha = 1 - \nu$ 是插值系数。

使用Taylor展开：
$$u_{i-1}^n = u_i^n - h\frac{\partial u}{\partial x} + \frac{h^2}{2}\frac{\partial^2 u}{\partial x^2} + O(h^3)$$

代入更新公式：
$$u_i^{n+1} = u_i^n - \nu h\frac{\partial u}{\partial x} + \frac{\nu h^2}{2}\frac{\partial^2 u}{\partial x^2} + O(h^3)$$

与精确解 $u_i^{n+1} = u_i^n - c\Delta t\frac{\partial u}{\partial x}$ 比较，得到修正方程：
$$\frac{\partial u}{\partial t} + c\frac{\partial u}{\partial x} = \frac{ch}{2}\nu(1-\nu)\frac{\partial^2 u}{\partial x^2}$$

右端项是数值扩散，扩散系数 $D_{num} = \frac{ch}{2}\nu(1-\nu)$ 总是正的，导致解的耗散。
</details>

**练习4.5** 设计一个测试案例验证你的流体求解器的收敛阶。考虑Taylor-Green涡旋：
$$u(x,y,t) = -\cos(x)\sin(y)e^{-2\nu t}$$
$$v(x,y,t) = \sin(x)\cos(y)e^{-2\nu t}$$
$$p(x,y,t) = -\frac{1}{4}[\cos(2x) + \cos(2y)]e^{-4\nu t}$$

<details>
<summary>提示</summary>
这是Navier-Stokes方程的精确解，可以计算不同网格分辨率下的L2误差。
</details>

<details>
<summary>答案</summary>
1. 验证此解满足Navier-Stokes方程（代入可验证）
2. 在 $[0, 2\pi]^2$ 域上使用周期边界条件
3. 计算L2误差：$e_h = \sqrt{\frac{1}{N}\sum_{i,j}(u_{computed} - u_{exact})^2}$
4. 对不同网格尺寸 $h = 2\pi/N$，计算误差
5. 收敛阶：$p = \log(e_h/e_{h/2})/\log(2)$
6. 预期：空间二阶精度 $e_h = O(h^2)$
</details>

**练习4.6** 实现一个自适应时间步长策略，基于CFL条件和粘性稳定性条件自动选择最大允许时间步长。

<details>  
<summary>提示</summary>
需要同时考虑对流CFL和扩散稳定性条件。
</details>

<details>
<summary>答案</summary>
```python
@ti.kernel
def compute_max_dt() -> float:
    # 对流CFL条件
    max_vel = 0.0
    for i, j in u:
        max_vel = max(max_vel, abs(u[i,j]))
    for i, j in v:
        max_vel = max(max_vel, abs(v[i,j]))
    
    dt_cfl = CFL * dx / (max_vel + 1e-8)
    
    # 粘性稳定性条件（如果使用显式粘性）
    dt_visc = 0.25 * dx * dx / (nu + 1e-8)
    
    # 表面张力稳定性条件（如果有）
    dt_tension = ti.sqrt(rho * dx**3 / (sigma + 1e-8))
    
    return min(dt_cfl, dt_visc, dt_tension)
```
</details>

**练习4.7** 推导并实现涡量-流函数形式的不可压缩流动。这种形式自动满足不可压缩条件，避免了压力投影。

<details>
<summary>提示</summary>  
涡量 $\omega = \nabla \times \mathbf{u}$，流函数 $\psi$ 满足 $u = \frac{\partial \psi}{\partial y}$，$v = -\frac{\partial \psi}{\partial x}$。
</details>

<details>
<summary>答案</summary>
二维涡量输送方程：
$$\frac{\partial \omega}{\partial t} + \mathbf{u} \cdot \nabla \omega = \nu \nabla^2 \omega$$

流函数泊松方程：
$$\nabla^2 \psi = -\omega$$

边界条件：
- 无滑移：$\psi = const$，$\frac{\partial \psi}{\partial n} = 0$
- 自由滑移：$\psi = const$，$\omega = 0$

算法步骤：
1. 输送涡量（使用半拉格朗日或高阶方法）
2. 求解流函数泊松方程
3. 从流函数计算速度场
4. 重复

这种方法的优点是自动满足不可压缩性，缺点是边界条件处理较复杂，且难以推广到三维。
</details>

## 常见陷阱与错误

1. **压力求解器不收敛**
   - 检查边界条件是否正确设置
   - 对于全Neumann边界，需要固定一个参考压力点
   - 确保离散算子的对称性

2. **速度场发散**
   - CFL条件可能太宽松
   - 检查边界条件实现
   - 压力投影可能不充分（迭代次数不够）

3. **MAC网格索引混淆**
   - 记住速度分量存储在不同位置
   - 插值时要使用正确的网格点
   - 绘图时需要将速度插值到单元中心

4. **数值耗散过大**
   - 使用高阶输送方法（MacCormack、BFECC）
   - 减小时间步长
   - 考虑使用涡量约束

5. **自由表面不稳定**
   - 表面张力的显式处理可能导致不稳定
   - 使用隐式或半隐式表面张力
   - 限制界面附近的时间步长

## 最佳实践检查清单

- [ ] 使用MAC网格避免压力振荡
- [ ] 实现自适应时间步长确保稳定性
- [ ] 对对流项使用至少二阶精度方法
- [ ] 压力求解器使用预条件共轭梯度法
- [ ] 实现质量守恒检查
- [ ] 边界条件处理要保持离散算子的对称性
- [ ] 使用高阶插值减少数值耗散
- [ ] 定期验证不可压缩条件的满足程度
- [ ] 实现能量或涡量监控以检测数值耗散
- [ ] 对复杂边界使用浸入边界法或切割单元法