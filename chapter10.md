# 第十章：可微编程与机器学习

可微编程正在革新物理仿真领域。传统的物理引擎是"前向"的——给定初始条件和参数，计算系统的演化。而可微编程使我们能够"反向"求解：给定期望的结果，自动推导出所需的参数或控制策略。这种端到端的优化能力为物理仿真开辟了全新的应用场景，从机器人控制到材料设计，从流体优化到动画生成。

本章将深入探讨可微物理仿真的核心技术，包括自动微分的实现原理、优化算法的选择策略、逆问题的求解框架，以及与机器学习的结合方法。我们将学习如何利用Taichi的自动微分功能构建可微仿真系统，并通过丰富的案例展示这些技术的实际应用。掌握这些内容后，读者将能够设计和实现智能化的物理系统，让仿真不仅能预测未来，更能优化决策。

## 10.1 可微仿真基础

可微仿真的核心思想是将整个物理系统视为一个可微函数，从输入参数到最终结果的每一步计算都保留梯度信息。这使得我们能够通过反向传播计算任意参数对结果的影响，进而优化系统行为。

### 10.1.1 端到端优化

考虑一个简单的弹道问题：给定目标位置，求解最优发射角度。传统方法需要多次试射或解析求解，而可微仿真可以直接优化：

$$\mathcal{L}(\theta) = \|x_{final}(\theta) - x_{target}\|^2$$

其中$x_{final}(\theta)$是发射角为$\theta$时的最终位置。通过自动微分计算$\frac{\partial \mathcal{L}}{\partial \theta}$，我们可以使用梯度下降找到最优角度。

更复杂的例子包括：
- **软体机器人控制**：优化驱动器序列使机器人到达目标位置
- **流体形状优化**：设计最小阻力的物体形状
- **材料参数识别**：从实验数据反推材料属性

端到端优化的关键是保持整个计算链的可微性。这要求我们：
1. 使用可微的数值格式（如可微的碰撞处理）
2. 避免不可微操作（如硬阈值、离散选择）
3. 处理数值稳定性问题（避免梯度爆炸/消失）

### 10.1.2 参数敏感度分析

在优化之前，理解参数如何影响系统行为至关重要。敏感度分析计算：

$$S_{ij} = \frac{\partial y_i}{\partial \theta_j}$$

其中$y_i$是第$i$个输出，$\theta_j$是第$j$个参数。敏感度矩阵$S$告诉我们：
- 哪些参数对结果影响最大
- 参数间是否存在耦合效应
- 优化问题是否良态（well-posed）

对于时间相关问题，我们需要考虑敏感度的时间演化：

$$\frac{d}{dt}\frac{\partial x}{\partial \theta} = \frac{\partial f}{\partial x}\frac{\partial x}{\partial \theta} + \frac{\partial f}{\partial \theta}$$

这是一个关于敏感度的ODE，可以与原系统一起求解。

### 10.1.3 伴随方法

当参数维度远大于目标维度时（如优化整个速度场），直接计算所有梯度的成本过高。伴随方法（Adjoint Method）提供了高效的解决方案。

考虑优化问题：
$$\min_\theta \mathcal{L}(x(T), \theta) \quad \text{s.t.} \quad \dot{x} = f(x, \theta, t)$$

引入拉格朗日乘子$\lambda(t)$（称为伴随变量），构造拉格朗日函数：

$$\mathcal{L}_{aug} = \mathcal{L}(x(T), \theta) + \int_0^T \lambda^T[\dot{x} - f(x, \theta, t)]dt$$

通过变分法，我们得到伴随方程：

$$\dot{\lambda} = -\left(\frac{\partial f}{\partial x}\right)^T \lambda, \quad \lambda(T) = \frac{\partial \mathcal{L}}{\partial x(T)}$$

梯度计算公式为：

$$\frac{d\mathcal{L}}{d\theta} = \int_0^T \lambda^T \frac{\partial f}{\partial \theta} dt$$

伴随方法的优势：
- 计算成本与参数数量无关
- 只需一次前向求解和一次反向求解
- 内存需求可通过Checkpointing控制

### 10.1.4 Checkpointing技术

可微仿真的主要挑战是内存需求——反向传播需要存储所有中间状态。对于长时间仿真，这可能需要TB级内存。Checkpointing通过时间-内存权衡解决这个问题。

**基本策略**：
1. 前向传播时只保存部分时刻的状态（checkpoints）
2. 反向传播时从最近的checkpoint重新计算所需状态

**最优Checkpoint分布**：
对于$n$个时间步和$c$个checkpoints，最优分布遵循对数规律：

$$t_i = T\left(\frac{i}{c}\right)^{p}, \quad p = \frac{\log c}{\log n}$$

**递归二分策略**：
```
checkpoint_binary(start, end, budget):
    if budget == 0:
        return []
    mid = (start + end) // 2
    left_budget = budget // 2
    right_budget = budget - left_budget - 1
    return checkpoint_binary(start, mid, left_budget) + [mid] + 
           checkpoint_binary(mid, end, right_budget)
```

这种策略保证重计算次数为$O(\log n)$。

**实践考虑**：
- Checkpoint粒度：太细增加开销，太粗增加重计算
- 异步Checkpoint：利用GPU计算时异步保存到CPU/磁盘
- 自适应策略：根据内存压力动态调整checkpoint密度

通过合理使用Checkpointing，我们可以在有限内存下处理任意长度的可微仿真。

## 10.2 自动微分实现

自动微分（Automatic Differentiation, AD）是可微编程的核心技术。与数值微分（有限差分）和符号微分不同，AD能够高效且精确地计算导数，没有截断误差，计算复杂度仅为原函数的常数倍。

### 10.2.1 计算图构建

自动微分的基础是将计算过程表示为有向无环图（DAG），其中节点代表变量，边代表操作。

**基本元素**：
- **变量节点**：输入参数、中间结果、输出值
- **操作节点**：基本运算（+、-、×、÷）、函数调用（sin、exp）
- **边**：数据依赖关系

**双数（Dual Number）表示**：
对于前向模式AD，每个变量同时存储值和导数：
$$\tilde{x} = (x, \dot{x}) = x + \dot{x}\epsilon$$

其中$\epsilon^2 = 0$。运算规则：
- $\tilde{x} + \tilde{y} = (x + y, \dot{x} + \dot{y})$
- $\tilde{x} \times \tilde{y} = (xy, x\dot{y} + y\dot{x})$
- $f(\tilde{x}) = (f(x), f'(x)\dot{x})$

**计算图示例**：
考虑函数$f(x, y) = x^2 + \sin(xy)$：
```
x ─┬─ [×] ─ x² ─┐
   │              ├─ [+] ─ f
   └─ [×] ─ xy ─ [sin] ─┘
     │
y ───┘
```

**动态图vs静态图**：
- **静态图**：编译时构建，优化机会多，如TensorFlow 1.x
- **动态图**：运行时构建，灵活性高，如PyTorch、Taichi

Taichi使用静态编译与动态执行的混合模式，在kernel编译时构建计算图，运行时高效执行。

### 10.2.2 前向vs反向传播

**前向模式（Forward Mode）**：
计算方向与原函数相同，适合输入维度小于输出维度的情况。

对于$y = f(x)$，前向模式计算：
$$\dot{y} = \frac{\partial f}{\partial x}\dot{x}$$

算法流程：
1. 初始化：$\dot{x}_i = e_i$（第$i$个单位向量）
2. 前向计算：对每个操作同时计算值和导数
3. 输出：得到第$i$列的Jacobian矩阵

**反向模式（Reverse Mode）**：
计算方向与原函数相反，适合输出维度小于输入维度的情况（最常见）。

引入伴随变量$\bar{v} = \frac{\partial \mathcal{L}}{\partial v}$，反向模式计算：
$$\bar{x} = \frac{\partial f}{\partial x}^T \bar{y}$$

算法流程：
1. 前向传播：计算所有中间值
2. 初始化：$\bar{y} = 1$
3. 反向传播：从输出到输入计算伴随变量
4. 累积梯度：$\bar{x}_i += \frac{\partial v}{\partial x_i}\bar{v}$

**基本操作的反向规则**：
- 加法：$z = x + y \Rightarrow \bar{x} += \bar{z}, \bar{y} += \bar{z}$
- 乘法：$z = xy \Rightarrow \bar{x} += y\bar{z}, \bar{y} += x\bar{z}$
- 除法：$z = x/y \Rightarrow \bar{x} += \bar{z}/y, \bar{y} -= x\bar{z}/y^2$
- 幂函数：$z = x^n \Rightarrow \bar{x} += nx^{n-1}\bar{z}$

### 10.2.3 高阶导数

计算Hessian矩阵$H_{ij} = \frac{\partial^2 f}{\partial x_i \partial x_j}$对于牛顿法等二阶优化算法至关重要。

**方法一：嵌套微分**
对梯度再次应用AD：
```python
# 一阶导数
def grad_f(x):
    return autograd(f, x)

# 二阶导数
def hessian_f(x):
    return autograd(grad_f, x)
```

**方法二：前向对反向（Forward-over-Reverse）**
1. 反向模式计算梯度$g(x) = \nabla f(x)$
2. 前向模式计算$Hv = \nabla g(x) \cdot v$

这种方法计算Hessian向量积的成本仅为$O(1)$次函数评估。

**方法三：边对边（Edge Pushing）**
直接在计算图上传播二阶导数信息：
$$\frac{\partial^2 f}{\partial x_i \partial x_j} = \sum_{k} \frac{\partial f}{\partial v_k} \frac{\partial^2 v_k}{\partial x_i \partial x_j} + \sum_{k,l} \frac{\partial^2 f}{\partial v_k \partial v_l} \frac{\partial v_k}{\partial x_i} \frac{\partial v_l}{\partial x_j}$$

**稀疏Hessian的利用**：
物理问题的Hessian通常稀疏（如FEM刚度矩阵）。可以通过图着色（Graph Coloring）减少计算：
- 将变量分组，同组变量的Hessian元素不重叠
- 每组只需一次前向传播
- 总计算量从$O(n)$降至$O(\chi)$，其中$\chi$是着色数

### 10.2.4 混合模式AD

实际应用中，纯前向或纯反向模式都不是最优的。混合模式根据问题结构选择最佳策略。

**矩阵链式法则**：
对于复合函数$h = f \circ g$，Jacobian矩阵满足：
$$J_h = J_f J_g$$

矩阵乘法的顺序影响计算量。设$J_f \in \mathbb{R}^{m \times k}$，$J_g \in \mathbb{R}^{k \times n}$：
- 左乘$(J_f J_g)$：$O(mkn)$
- 右乘$J_f (J_g)$：$O(mkn)$
- 但如果先计算$v^T J_f$：$O(mk + kn)$

**最优Jacobian累积**：
Bauer定理：寻找最优计算顺序是NP困难的。实用启发式：
1. **贪心选择**：每步选择计算量最小的操作
2. **动态规划**：对小规模问题求解最优顺序
3. **预设模式**：根据问题类型选择前向/反向/混合

**稀疏性利用**：
物理仿真中的Jacobian通常高度稀疏：
- **压缩存储**：只存储非零元素
- **稀疏AD**：只计算非零元素的导数
- **符号分析**：预先确定稀疏模式

**实现技巧**：
1. **表达式模板**：延迟计算，优化整个表达式
2. **算子融合**：将多个操作合并，减少内存访问
3. **内存池**：预分配内存，避免动态分配开销
4. **并行化**：独立的导数计算可以并行执行

通过合理使用这些技术，现代AD系统可以达到接近手写导数代码的性能，同时保持自动化的便利性。

## 10.3 优化算法

有了梯度信息后，选择合适的优化算法至关重要。物理仿真的优化问题通常具有特殊结构：高维、非凸、存在约束、计算昂贵。本节介绍各类优化算法及其在物理问题中的应用。

### 10.3.1 梯度下降变种

**基础梯度下降**：
$$\theta_{k+1} = \theta_k - \alpha \nabla f(\theta_k)$$

学习率$\alpha$的选择至关重要。过大导致振荡，过小收敛缓慢。

**动量法（Momentum）**：
引入速度变量累积历史梯度信息：
$$v_{k+1} = \beta v_k + (1-\beta) \nabla f(\theta_k)$$
$$\theta_{k+1} = \theta_k - \alpha v_{k+1}$$

物理解释：将优化过程视为阻尼系统中的质点运动，$\beta$是阻尼系数。

**Nesterov加速梯度（NAG）**：
先根据动量"展望"未来位置，再计算梯度：
$$v_{k+1} = \beta v_k + (1-\beta) \nabla f(\theta_k - \alpha \beta v_k)$$
$$\theta_{k+1} = \theta_k - \alpha v_{k+1}$$

收敛速度从$O(1/k)$提升到$O(1/k^2)$。

**Adam优化器**：
自适应学习率，对每个参数维护一阶矩和二阶矩估计：
$$m_{k+1} = \beta_1 m_k + (1-\beta_1) \nabla f(\theta_k)$$
$$v_{k+1} = \beta_2 v_k + (1-\beta_2) [\nabla f(\theta_k)]^2$$
$$\hat{m} = m_{k+1}/(1-\beta_1^{k+1}), \quad \hat{v} = v_{k+1}/(1-\beta_2^{k+1})$$
$$\theta_{k+1} = \theta_k - \alpha \hat{m}/(\sqrt{\hat{v}} + \epsilon)$$

默认参数：$\beta_1=0.9$，$\beta_2=0.999$，$\epsilon=10^{-8}$。

**RMSprop**：
仅使用二阶矩，适合非平稳目标：
$$v_{k+1} = \beta v_k + (1-\beta) [\nabla f(\theta_k)]^2$$
$$\theta_{k+1} = \theta_k - \alpha \nabla f(\theta_k)/\sqrt{v_{k+1} + \epsilon}$$

**AdaGrad**：
累积所有历史梯度平方：
$$G_{k+1} = G_k + [\nabla f(\theta_k)]^2$$
$$\theta_{k+1} = \theta_k - \alpha \nabla f(\theta_k)/\sqrt{G_{k+1} + \epsilon}$$

适合稀疏梯度，但学习率单调递减可能过早停止。

### 10.3.2 拟牛顿方法

牛顿法使用二阶信息加速收敛：
$$\theta_{k+1} = \theta_k - H^{-1} \nabla f(\theta_k)$$

其中$H$是Hessian矩阵。直接计算$H$成本高昂，拟牛顿方法通过历史信息近似。

**BFGS算法**：
维护Hessian逆的近似$B_k \approx H^{-1}$：
$$s_k = \theta_{k+1} - \theta_k, \quad y_k = \nabla f(\theta_{k+1}) - \nabla f(\theta_k)$$
$$B_{k+1} = (I - \rho_k s_k y_k^T) B_k (I - \rho_k y_k s_k^T) + \rho_k s_k s_k^T$$
其中$\rho_k = 1/(y_k^T s_k)$。

**L-BFGS（限制内存BFGS）**：
只存储最近$m$步的$(s_k, y_k)$对，通过两步递归计算搜索方向：
```
# 第一步：从右向左递归
q = ∇f(θ)
for i = k-1, ..., k-m:
    α[i] = ρ[i] * s[i]ᵀ * q
    q = q - α[i] * y[i]

# 初始化：使用简单的对角近似
r = H₀ * q

# 第二步：从左向右递归  
for i = k-m, ..., k-1:
    β = ρ[i] * y[i]ᵀ * r
    r = r + s[i] * (α[i] - β)

搜索方向 = -r
```

内存需求：$O(mn)$，其中$m$通常取3-20。

**线搜索策略**：
拟牛顿方法需要线搜索确定步长：
1. **Armijo条件**：$f(\theta + \alpha p) \leq f(\theta) + c_1 \alpha \nabla f(\theta)^T p$
2. **Wolfe条件**：加上曲率条件$\nabla f(\theta + \alpha p)^T p \geq c_2 \nabla f(\theta)^T p$
3. **强Wolfe条件**：$|\nabla f(\theta + \alpha p)^T p| \leq c_2 |\nabla f(\theta)^T p|$

典型参数：$c_1 = 10^{-4}$，$c_2 = 0.9$。

### 10.3.3 约束优化

物理问题常带约束，如不可压缩性、接触约束、关节限制等。

**问题形式**：
$$\min_\theta f(\theta) \quad \text{s.t.} \quad g_i(\theta) \leq 0, \quad h_j(\theta) = 0$$

**罚函数法**：
将约束转化为罚项：
$$f_{penalty}(\theta) = f(\theta) + \mu \sum_i \max(0, g_i(\theta))^2 + \mu \sum_j h_j(\theta)^2$$

随着$\mu \to \infty$，解收敛到原问题。但$\mu$过大导致病态。

**增广拉格朗日法（ALM）**：
结合拉格朗日乘子和罚函数：
$$\mathcal{L}_\mu(\theta, \lambda, \nu) = f(\theta) + \sum_i \lambda_i g_i(\theta) + \frac{\mu}{2}\sum_i \max(0, g_i(\theta))^2 + \sum_j \nu_j h_j(\theta) + \frac{\mu}{2}\sum_j h_j(\theta)^2$$

交替更新：
1. 固定$(\lambda, \nu)$，优化$\theta$
2. 更新乘子：$\lambda_i \leftarrow \max(0, \lambda_i + \mu g_i(\theta))$，$\nu_j \leftarrow \nu_j + \mu h_j(\theta)$

**内点法（Interior Point Method）**：
使用障碍函数处理不等式约束：
$$f_{barrier}(\theta) = f(\theta) - \tau \sum_i \log(-g_i(\theta))$$

随着$\tau \to 0$，解收敛。需要保持严格可行性$g_i(\theta) < 0$。

**投影梯度法**：
每步后投影到可行域：
$$\theta_{k+1} = \Pi_\mathcal{C}[\theta_k - \alpha \nabla f(\theta_k)]$$

对于盒约束$a \leq \theta \leq b$，投影简单：
$$[\Pi_\mathcal{C}(\theta)]_i = \max(a_i, \min(b_i, \theta_i))$$

### 10.3.4 随机优化

当梯度计算昂贵或存在噪声时，随机优化方法更实用。

**随机梯度下降（SGD）**：
使用梯度的无偏估计：
$$\theta_{k+1} = \theta_k - \alpha_k \tilde{\nabla} f(\theta_k)$$

其中$\mathbb{E}[\tilde{\nabla} f] = \nabla f$。学习率需满足Robbins-Monro条件：
$$\sum_k \alpha_k = \infty, \quad \sum_k \alpha_k^2 < \infty$$

**方差缩减技术**：
- **SVRG**：周期性计算完整梯度作为锚点
- **SAGA**：存储每个样本的历史梯度
- **SAG**：使用历史梯度的平均

**协方差矩阵自适应（CMA-ES）**：
维护高斯分布$\mathcal{N}(\mu, \sigma^2 C)$采样候选解：
1. 采样：$x_i \sim \mathcal{N}(\mu, \sigma^2 C)$
2. 评估并排序
3. 更新均值：$\mu \leftarrow \sum_{i=1}^{\lambda} w_i x_{i:λ}$
4. 更新协方差：考虑进化路径和秩-1、秩-$\mu$更新

适合非凸、多峰、噪声问题，但采样效率低。

**粒子群优化（PSO）**：
模拟鸟群觅食行为：
$$v_i \leftarrow \omega v_i + c_1 r_1 (p_i - x_i) + c_2 r_2 (g - x_i)$$
$$x_i \leftarrow x_i + v_i$$

其中$p_i$是粒子历史最优，$g$是全局最优。

**模拟退火（SA）**：
以概率接受劣解跳出局部最优：
$$P(accept) = \begin{cases}
1 & \text{if } f(x_{new}) < f(x_{current}) \\
\exp(-(f(x_{new}) - f(x_{current}))/T) & \text{otherwise}
\end{cases}$$

温度$T$逐渐降低，最终收敛到全局最优（理论上）。

**实践建议**：
1. **问题规模小（<100维）**：拟牛顿方法（L-BFGS）
2. **问题规模大但光滑**：Adam或其变种
3. **存在约束**：增广拉格朗日或投影方法
4. **非凸多峰**：CMA-ES或多起点局部搜索
5. **梯度噪声大**：方差缩减SGD或无梯度方法