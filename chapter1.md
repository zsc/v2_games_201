# 第一章：导论

本章将为读者建立物理仿真的基础知识框架，并深入介绍Taichi编程语言的核心特性。我们将从物理引擎的基本概念出发，逐步深入到高性能并行计算的实现细节。通过本章学习，读者将掌握使用Taichi进行物理仿真开发所需的全部基础知识。

## 1.1 基于物理的动画简介

### 1.1.1 物理引擎的定义与应用领域

物理引擎是一类专门用于模拟物理系统行为的计算机软件。它通过数值方法求解物理方程，为虚拟世界中的物体提供逼真的运动和交互效果。现代物理引擎通常包含三大核心组件：

1. **刚体动力学**：模拟不可变形物体的运动，处理碰撞检测与响应
2. **软体动力学**：模拟可变形物体如布料、橡胶、肌肉等的形变行为  
3. **流体动力学**：模拟液体和气体的流动，包括水波、烟雾等效果

物理引擎的应用领域极其广泛：
- **游戏产业**：从简单的《愤怒的小鸟》到复杂的《围攻》(Besiege)，物理引擎让游戏世界更加真实
- **影视特效**：迪士尼、皮克斯等顶级动画工作室依赖物理仿真创造震撼视觉效果
- **工程仿真**：汽车碰撞测试、建筑结构分析、流体机械设计等
- **虚拟现实**：提供真实的触觉反馈和物理交互体验
- **机器人学**：训练和验证机器人控制算法的虚拟环境

### 1.1.2 计算机图形学中的物理仿真

在计算机图形学的发展历程中，物理仿真扮演着越来越重要的角色。早期的动画完全依赖艺术家手工关键帧，而现代制作流程已经深度集成了物理仿真技术：

**传统动画 vs 物理动画**：
- 传统动画：艺术家控制每一帧，灵活但耗时
- 物理动画：定义初始条件和物理参数，自动生成逼真运动

**里程碑作品**：
- 1997年《泰坦尼克号》：大规模刚体破碎和流体仿真
- 2013年《冰雪奇缘》：基于物质点法(MPM)的雪景仿真
- 2016年《海洋奇缘》：高度逼真的水体和头发仿真

**技术演进**：
```
1980s: 简单粒子系统
1990s: 刚体动力学成熟
2000s: 流体仿真突破（SPH、FLIP）
2010s: 统一仿真框架（MPM）
2020s: 机器学习增强、实时光线追踪物理
```

### 1.1.3 工程CAE与游戏物理的差异

虽然都是物理仿真，但工程计算机辅助工程(CAE)与游戏物理在设计理念上存在本质差异：

| 特性 | 工程CAE | 游戏物理 |
|------|---------|----------|
| **精度要求** | 误差 < 1%，需要验证 | 视觉合理即可 |
| **实时性** | 可以离线计算数小时 | 必须 60 FPS |
| **稳定性** | 绝对稳定，不允许崩溃 | 偶尔穿模可接受 |
| **复杂度** | 百万级自由度 | 千级自由度 |
| **物理真实性** | 严格遵循物理定律 | 可以作弊提升体验 |

**工程CAE的典型流程**：
1. 几何建模（CAD）
2. 网格生成（六面体/四面体）
3. 边界条件设置
4. 求解器计算（可能需要数小时）
5. 后处理可视化

**游戏物理的设计原则**：
1. 性能优先：宁可牺牲精度也要保证帧率
2. 鲁棒性：处理各种极端情况不崩溃
3. 可控性：艺术家/设计师可以调整参数
4. 确定性：相同输入产生相同结果（重要于网络同步）

### 1.1.4 实时性与精确性的权衡

在物理仿真中，实时性和精确性往往不可兼得。理解这种权衡对于选择合适的算法至关重要：

**时间步长的选择**：
- 显式积分器的稳定性条件：$\Delta t \leq c\sqrt{\frac{m}{k}}$
- 隐式积分器：无条件稳定但每步计算量大
- 自适应时间步长：根据场景动态调整

**模型简化技术**：
1. **几何简化**：
   - 碰撞用简化凸包代替复杂网格
   - LOD (Level of Detail) 根据距离切换精度
   
2. **物理简化**：
   - 刚体近似：忽略小变形
   - 质点近似：远距离物体视为质点
   - 维度降低：3D布料用2D壳体模拟

3. **数值简化**：
   - 低阶积分格式
   - 稀疏更新（只更新活跃区域）
   - 预计算查找表

**性能优化层次**：
```
算法级: O(n²) → O(n log n) → O(n)
并行级: CPU多线程 → GPU大规模并行
近似级: 精确解 → 数值解 → 启发式解
```

## 1.2 Taichi编程语言简介

### 1.2.1 安装与环境配置

Taichi是一个专为高性能数值计算设计的编程语言，特别适合物理仿真。其最大特点是将Python的易用性与C++的性能完美结合。

**安装方式**：
```bash
# 推荐使用pip安装
pip install taichi

# 对于开发版本
pip install taichi-nightly

# 验证安装
python -c "import taichi as ti; ti.init(); print(ti.__version__)"
```

**支持平台**：
- 操作系统：Windows、Linux、macOS
- Python版本：3.6-3.10
- 后端支持：
  - CPU：x64、ARM64
  - GPU：CUDA（NVIDIA）、Vulkan（跨平台）、Metal（Apple）
  - WebGPU：浏览器端运行

**初始化配置**：
```python
import taichi as ti

# 基本初始化
ti.init(arch=ti.cuda)  # 选择CUDA后端

# 高级配置
ti.init(
    arch=ti.gpu,           # 自动选择可用GPU
    device_memory_GB=4,    # 限制GPU内存使用
    debug=True,            # 开启调试模式
    print_ir=True,         # 打印中间表示
    random_seed=42         # 固定随机种子
)
```

### 1.2.2 基本数据类型与张量操作

Taichi提供了丰富的数据类型和张量操作，专门针对科学计算优化：

**基本数据类型**：
```python
# 标量类型
ti.i32, ti.i64    # 32/64位整数
ti.f32, ti.f64    # 32/64位浮点数
ti.u8, ti.u16     # 无符号整数

# 类型推断
x = ti.field(dtype=ti.f32, shape=(100, 100))
y = ti.field(dtype=ti.i32, shape=100)
```

**张量（Tensor）定义**：
```python
# 标量场
temperature = ti.field(dtype=ti.f32, shape=(256, 256))

# 向量场
velocity = ti.Vector.field(2, dtype=ti.f32, shape=(100, 100))

# 矩阵场
stress = ti.Matrix.field(3, 3, dtype=ti.f32, shape=(50, 50, 50))

# 动态形状（高级用法）
particles = ti.Vector.field(3, dtype=ti.f32)
ti.root.dynamic(ti.i, 1024).place(particles)
```

**张量操作**：
```python
@ti.kernel
def tensor_operations():
    # 向量运算
    a = ti.Vector([1.0, 2.0, 3.0])
    b = ti.Vector([4.0, 5.0, 6.0])
    
    dot_product = a.dot(b)              # 点积: 32.0
    cross_product = a.cross(b)          # 叉积: [-3, 6, -3]
    norm = a.norm()                     # 范数: sqrt(14)
    normalized = a.normalized()         # 单位向量
    
    # 矩阵运算
    M = ti.Matrix([[1, 2], [3, 4]])
    M_inv = M.inverse()                 # 逆矩阵
    M_T = M.transpose()                 # 转置
    det = M.determinant()               # 行列式: -2
    
    # 特殊矩阵
    I = ti.Matrix.identity(ti.f32, 3)   # 3x3单位矩阵
    Z = ti.Matrix.zero(ti.f32, 2, 3)   # 2x3零矩阵
```

### 1.2.3 核函数(kernel)与函数(func)

Taichi的计算核心是kernel和func，它们定义了可以在GPU上高效执行的计算：

**Kernel函数**：
```python
@ti.kernel
def add_one(field: ti.template()):
    # Kernel是Taichi程序的入口点
    for i, j in field:
        field[i, j] += 1.0

# Kernel可以有返回值
@ti.kernel
def sum_field(field: ti.template()) -> ti.f32:
    total = 0.0
    for i, j in field:
        total += field[i, j]
    return total
```

**Func函数**：
```python
@ti.func
def complex_calculation(x: ti.f32, y: ti.f32) -> ti.f32:
    # Func只能被Kernel或其他Func调用
    # 会被强制内联，无函数调用开销
    return ti.sqrt(x * x + y * y)

@ti.kernel
def use_func(field: ti.template()):
    for i, j in field:
        field[i, j] = complex_calculation(i * 0.1, j * 0.1)
```

**类型提示与模板**：
```python
# 严格类型提示
@ti.kernel
def typed_kernel(x: ti.types.ndarray(dtype=ti.f32, ndim=2)):
    for i, j in ti.ndrange(x.shape[0], x.shape[1]):
        x[i, j] *= 2.0

# 模板参数
@ti.kernel
def generic_kernel(field: ti.template(), scalar: ti.f32):
    for I in ti.grouped(field):
        field[I] *= scalar
```

**编译时优化**：
- 常量折叠：编译时计算常量表达式
- 死代码消除：移除不可达代码
- 循环展开：小循环自动展开
- 向量化：自动SIMD指令生成

### 1.2.4 并行for循环与原子操作

Taichi的一个核心特性是自动并行化，开发者只需要写串行代码，编译器会自动并行执行：

**并行For循环**：
```python
@ti.kernel
def parallel_computation(n: ti.i32):
    # 这个循环会自动并行执行
    for i in range(n):
        # 每个迭代独立，无数据竞争
        result[i] = ti.sin(i * 0.1) + ti.cos(i * 0.2)
    
    # 多维并行循环
    for i, j in ti.ndrange(100, 200):
        matrix[i, j] = i + j
```

**原子操作**：
```python
@ti.kernel 
def atomic_operations():
    # 处理并发写入冲突
    for i in range(1000):
        # 原子加法，保证线程安全
        ti.atomic_add(total[0], 1.0)
        
        # 原子最大值
        ti.atomic_max(max_value[0], values[i])
        
        # 原子交换
        old = ti.atomic_sub(counter[0], 1)
    
    # 复合原子操作
    for i in range(n):
        # 自动识别并转换为原子操作
        histogram[values[i]] += 1
```

**竞态条件的避免**：
```python
# 错误示例 - 有竞态条件
@ti.kernel
def race_condition():
    for i in range(100):
        # 多个线程可能同时读写
        shared_sum[0] = shared_sum[0] + i  # 危险！

# 正确示例 - 使用原子操作
@ti.kernel 
def no_race_condition():
    for i in range(100):
        ti.atomic_add(shared_sum[0], i)  # 安全
```

**性能考虑**：
- 原子操作比普通操作慢10-100倍
- 尽量减少原子操作的使用
- 考虑使用归约(reduction)模式
- 合理的数据布局可以避免冲突

## 1.3 Taichi的自动并行化

### 1.3.1 并行计算基础概念

理解Taichi的自动并行化机制，首先需要掌握并行计算的基本概念：

**并行计算模型**：
1. **数据并行(Data Parallelism)**：
   - 相同操作应用于不同数据
   - SIMD (Single Instruction Multiple Data)
   - 适合规则的数组运算

2. **任务并行(Task Parallelism)**：
   - 不同任务同时执行
   - MIMD (Multiple Instruction Multiple Data)
   - 适合异构计算任务

**GPU执行模型(SIMT)**：
```
GPU执行层次：
Grid
├── Block (0,0)
│   ├── Warp 0 (32 threads)
│   ├── Warp 1 (32 threads)
│   └── ...
├── Block (0,1)
└── ...
```

**并行粒度**：
- 细粒度：每个数据元素一个线程
- 中粒度：数据块级别并行
- 粗粒度：任务级别并行

Taichi采用细粒度数据并行，自动将循环映射到GPU线程。

### 1.3.2 Range-for与Struct-for循环

Taichi提供两种主要的并行循环结构：

**Range-for循环**：
```python
@ti.kernel
def range_for_examples():
    # 一维循环
    for i in range(100):
        arr[i] = i * i
    
    # 多维循环 - 使用ti.ndrange
    for i, j in ti.ndrange(10, 20):
        matrix[i, j] = i * j
    
    # 带步长的循环
    for i in range(0, 100, 2):
        even_array[i // 2] = i
    
    # 嵌套循环优化
    for i, j, k in ti.ndrange(32, 32, 32):
        # 编译器会优化循环顺序以提高缓存局部性
        volume[i, j, k] = i + j + k
```

**Struct-for循环**：
```python
# 稀疏数据结构定义
block = ti.root.pointer(ti.ij, 32)
pixel = block.dense(ti.ij, 16)
pixel.place(sparse_field)

@ti.kernel
def struct_for_example():
    # 只遍历激活的元素
    for i, j in sparse_field:
        sparse_field[i, j] *= 2.0
    
    # 分层遍历
    for i, j in block:
        print(f"Active block at ({i}, {j})")
```

**ti.grouped优化**：
```python
@ti.kernel
def grouped_access():
    # 将多维索引打包为单个变量
    for I in ti.grouped(tensor_3d):
        # I是一个向量，包含所有维度的索引
        tensor_3d[I] = I.dot(I)  # 自动展开为 i*i + j*j + k*k
```

**并行化规则**：
1. 最外层循环自动并行
2. 内层循环保持串行
3. 循环迭代必须独立
4. 不支持break/continue

### 1.3.3 线程安全与竞态条件

在并行编程中，正确处理共享数据访问至关重要：

**竞态条件示例**：
```python
# 危险代码 - 竞态条件
@ti.kernel
def race_example():
    for i in range(1000000):
        # 多个线程同时读-修改-写
        counter[0] = counter[0] + 1  # 结果不确定！

# 安全代码 - 原子操作
@ti.kernel
def safe_example():
    for i in range(1000000):
        ti.atomic_add(counter[0], 1)  # 保证正确结果
```

**Taichi的安全保证**：
1. **自动原子化**：
   ```python
   # Taichi自动识别并转换
   a[i] += b[j]  # 自动转换为atomic_add
   a[i] = max(a[i], b[j])  # 自动转换为atomic_max
   ```

2. **支持的原子操作**：
   - `atomic_add`, `atomic_sub`
   - `atomic_min`, `atomic_max`
   - `atomic_and`, `atomic_or`, `atomic_xor`
   - `atomic_exchange` (原子交换)

3. **线程局部存储**：
   ```python
   @ti.kernel
   def thread_local_example():
       # 每个线程有自己的局部变量
       for i in range(n):
           local_sum = 0.0
           for j in range(m):
               local_sum += matrix[i, j]
           row_sums[i] = local_sum
   ```

**避免竞态的策略**：
1. 数据分区：每个线程处理独立数据
2. 归约模式：使用树形归约减少冲突
3. 双缓冲：读写分离到不同缓冲区
4. 原子操作：必要时使用但注意性能

### 1.3.4 性能分析与调优基础

理解和优化Taichi程序性能的关键工具和技术：

**内置性能分析器**：
```python
# 启用性能分析
ti.init(arch=ti.cuda, kernel_profiler=True)

@ti.kernel
def computation():
    for i in range(1000000):
        a[i] = ti.sin(b[i]) + ti.cos(c[i])

# 执行并打印性能报告
computation()
ti.print_kernel_profile_info()
```

**性能指标**：
1. **kernel时间**：GPU执行时间
2. **内存带宽**：数据传输速率
3. **占用率**：活跃线程比例
4. **缓存命中率**：数据局部性

**常见性能问题**：
```python
# 问题1：内存访问模式差
@ti.kernel
def bad_access_pattern():
    for i, j in ti.ndrange(1024, 1024):
        # 列优先访问，缓存不友好
        a[j, i] = b[j, i] + 1

# 优化：行优先访问
@ti.kernel 
def good_access_pattern():
    for i, j in ti.ndrange(1024, 1024):
        # 行优先访问，缓存友好
        a[i, j] = b[i, j] + 1

# 问题2：过多原子操作
@ti.kernel
def too_many_atomics():
    for i in range(n):
        ti.atomic_add(sum[data[i] % 10], 1)  # 高冲突

# 优化：线程局部累加
@ti.kernel
def optimized_histogram():
    # 使用线程局部数组减少冲突
    for i in range(n):
        local_hist = [0] * 10
        # 先本地累加
        for j in range(i, min(i + 1000, n)):
            local_hist[data[j] % 10] += 1
        # 再原子更新全局
        for k in range(10):
            ti.atomic_add(histogram[k], local_hist[k])
```

**优化技巧清单**：
1. 合并内存访问
2. 减少分支分歧
3. 使用共享内存
4. 循环展开
5. 避免存储bank冲突

## 1.4 Taichi程序的调试技巧

### 1.4.1 调试模式与边界检查

调试是开发高质量代码的关键环节，Taichi提供了丰富的调试工具：

**启用调试模式**：
```python
# 完整调试模式
ti.init(debug=True, arch=ti.cpu)

# 调试选项详解
ti.init(
    debug=True,           # 启用所有调试检查
    print_ir=True,        # 打印中间表示
    print_kernel_llvm_ir=True,  # 打印LLVM IR
    check_out_of_bound=True,    # 边界检查
)
```

**边界检查器**：
```python
@ti.kernel
def boundary_check_example():
    arr = ti.field(ti.f32, shape=100)
    
    # 调试模式下会检测越界
    for i in range(101):  # 错误：越界访问
        arr[i] = i  # 运行时错误：index 100 out of bound [0, 100)
    
    # 安全的访问方式
    for i in range(100):
        if 0 <= i < 100:  # 显式边界检查
            arr[i] = i
```

**断言使用**：
```python
@ti.kernel
def debug_with_assert():
    for i in range(n):
        # 运行时断言
        assert 0 <= values[i] <= 1.0, f"Value {values[i]} out of range"
        
        # 计算前置条件检查
        denominator = compute_denominator(i)
        assert denominator != 0, "Division by zero"
        result[i] = numerator[i] / denominator
```

### 1.4.2 打印调试与断言

Taichi支持在kernel中打印，这是调试的重要工具：

**打印功能**：
```python
@ti.kernel
def print_debugging():
    # 基本打印
    print("Starting computation")
    
    # 打印变量
    for i in range(5):
        print(f"i = {i}, value = {arr[i]}")
    
    # 条件打印
    for i in range(1000):
        if i % 100 == 0:
            print(f"Progress: {i}/1000")
    
    # 打印向量和矩阵
    vec = ti.Vector([1.0, 2.0, 3.0])
    mat = ti.Matrix([[1, 2], [3, 4]])
    print("Vector:", vec, "Matrix:", mat)
```

**打印限制与最佳实践**：
```python
@ti.kernel
def print_best_practices():
    # 限制打印数量避免输出爆炸
    print_count = 0
    for i in range(1000000):
        if arr[i] < 0 and print_count < 10:
            print(f"Negative value at {i}: {arr[i]}")
            print_count += 1
    
    # 使用统计而非逐个打印
    error_count = 0
    for i in range(n):
        if is_error(i):
            error_count += 1
    print(f"Total errors: {error_count}")
```

### 1.4.3 常见错误与解决方案

**错误类型1：越界访问**
```python
# 错误示例
@ti.kernel
def out_of_bound_error():
    arr = ti.field(ti.f32, shape=(10, 10))
    for i, j in ti.ndrange(11, 11):  # 错误！
        arr[i, j] = 0

# 解决方案
@ti.kernel
def safe_access():
    arr = ti.field(ti.f32, shape=(10, 10))
    for i, j in arr:  # 使用field自身作为范围
        arr[i, j] = 0
```

**错误类型2：类型不匹配**
```python
# 错误示例
@ti.kernel
def type_mismatch():
    int_field = ti.field(ti.i32, shape=100)
    for i in range(100):
        int_field[i] = 3.14  # 错误：浮点数赋值给整数

# 解决方案
@ti.kernel
def correct_types():
    int_field = ti.field(ti.i32, shape=100)
    for i in range(100):
        int_field[i] = int(3.14)  # 显式类型转换
```

**错误类型3：未初始化变量**
```python
# 错误示例
@ti.kernel
def uninitialized_error():
    sum = ti.field(ti.f32, shape=1)
    for i in range(100):
        sum[0] += arr[i]  # 错误：sum[0]未初始化

# 解决方案
@ti.kernel
def proper_initialization():
    sum = ti.field(ti.f32, shape=1)
    sum[0] = 0.0  # 显式初始化
    for i in range(100):
        sum[0] += arr[i]
```

**错误类型4：并行冲突**
```python
# 错误示例
@ti.kernel
def parallel_conflict():
    for i in range(n):
        # 多个线程写入相同位置
        shared_array[i % 10] = i  # 结果不确定

# 解决方案
@ti.kernel
def conflict_resolution():
    for i in range(n):
        # 使用原子操作
        ti.atomic_max(shared_array[i % 10], i)
```

### 1.4.4 性能瓶颈识别

识别和解决性能瓶颈是优化的关键：

**性能分析工具**：
```python
import time

# 手动计时
@ti.kernel
def timed_kernel():
    # kernel代码
    pass

start = time.time()
for _ in range(100):
    timed_kernel()
ti.sync()  # 重要：等待GPU完成
end = time.time()
print(f"Average time: {(end - start) / 100 * 1000:.2f} ms")

# 使用Taichi profiler
ti.init(kernel_profiler=True)
# 运行kernels
ti.print_kernel_profile_info("time")  # 按时间排序
```

**常见性能瓶颈**：

1. **内存带宽限制**：
```python
# 问题：随机内存访问
@ti.kernel
def random_access():
    for i in range(n):
        sum += arr[indices[i]]  # 缓存不友好

# 优化：顺序访问
@ti.kernel
def sequential_access():
    for i in range(n):
        sum += arr[i]  # 缓存友好
```

2. **计算瓶颈**：
```python
# 问题：复杂数学运算
@ti.kernel
def expensive_math():
    for i in range(n):
        result[i] = ti.sin(ti.cos(ti.tan(x[i])))

# 优化：查表或近似
@ti.kernel
def optimized_math():
    for i in range(n):
        # 使用泰勒展开或查找表
        result[i] = fast_approximation(x[i])
```

3. **核函数粒度**：
```python
# 问题：kernel太小，启动开销大
for i in range(1000):
    small_kernel(i)  # 1000次kernel启动

# 优化：批处理
large_kernel(0, 1000)  # 1次kernel启动
```

## 1.5 面向数据的编程(DOP)

### 1.5.1 AoS vs SoA布局

数据布局对性能有决定性影响，理解AoS(Array of Structures)和SoA(Structure of Arrays)的差异至关重要：

**AoS布局示例**：
```python
# Array of Structures - 数据按对象组织
@ti.data_oriented
class ParticleSystemAoS:
    def __init__(self, n):
        self.particles = ti.Struct.field({
            'x': ti.f32,
            'y': ti.f32,
            'vx': ti.f32,
            'vy': ti.f32,
            'mass': ti.f32
        }, shape=n)
    
    @ti.kernel
    def update_positions(self, dt: ti.f32):
        for i in self.particles:
            # 访问同一粒子的不同属性 - 缓存友好
            self.particles[i].x += self.particles[i].vx * dt
            self.particles[i].y += self.particles[i].vy * dt
```

**SoA布局示例**：
```python
# Structure of Arrays - 数据按属性组织
@ti.data_oriented
class ParticleSystemSoA:
    def __init__(self, n):
        self.x = ti.field(ti.f32, shape=n)
        self.y = ti.field(ti.f32, shape=n)
        self.vx = ti.field(ti.f32, shape=n)
        self.vy = ti.field(ti.f32, shape=n)
        self.mass = ti.field(ti.f32, shape=n)
    
    @ti.kernel
    def update_positions(self, dt: ti.f32):
        for i in range(self.x.shape[0]):
            # 访问不同数组的相同索引 - SIMD友好
            self.x[i] += self.vx[i] * dt
            self.y[i] += self.vy[i] * dt
```

**性能对比**：

| 操作类型 | AoS | SoA |
|---------|-----|-----|
| 单对象全属性访问 | 优秀 | 较差 |
| 批量单属性处理 | 较差 | 优秀 |
| SIMD向量化 | 困难 | 简单 |
| 缓存利用率 | 取决于访问模式 | 通常更好 |

**混合策略**：
```python
# AoSoA - Array of Structure of Arrays
# 结合两者优势
@ti.data_oriented
class HybridLayout:
    def __init__(self, n, block_size=32):
        self.blocks = ti.Struct.field({
            'x': ti.types.vector(block_size, ti.f32),
            'y': ti.types.vector(block_size, ti.f32),
            'vx': ti.types.vector(block_size, ti.f32),
            'vy': ti.types.vector(block_size, ti.f32)
        }, shape=n // block_size)
```

### 1.5.2 缓存友好的数据组织

现代处理器的性能瓶颈往往是内存访问而非计算，优化数据布局至关重要：

**缓存层级**：
```
寄存器: ~1 cycle
L1缓存: ~4 cycles (32KB)
L2缓存: ~12 cycles (256KB)
L3缓存: ~40 cycles (8MB)
主内存: ~200 cycles
```

**空间局部性优化**：
```python
# 差：跨步访问
@ti.kernel
def poor_spatial_locality():
    for i in range(n):
        for j in range(m):
            # 每次访问跳过m个元素
            sum += matrix[j * n + i]

# 好：连续访问
@ti.kernel  
def good_spatial_locality():
    for i in range(n):
        for j in range(m):
            # 连续内存访问
            sum += matrix[i * m + j]
```

**时间局部性优化**：
```python
# 差：反复加载相同数据
@ti.kernel
def poor_temporal_locality():
    for k in range(p):
        for i in range(n):
            for j in range(m):
                C[i, j] += A[i, k] * B[k, j]

# 好：分块矩阵乘法
@ti.kernel
def blocked_matmul():
    tile_size = 32
    for bi in range(0, n, tile_size):
        for bj in range(0, m, tile_size):
            for bk in range(0, p, tile_size):
                # 处理tile_size x tile_size的块
                for i in range(bi, min(bi + tile_size, n)):
                    for j in range(bj, min(bj + tile_size, m)):
                        for k in range(bk, min(bk + tile_size, p)):
                            C[i, j] += A[i, k] * B[k, j]
```

**数据对齐**：
```python
# 确保数据对齐到缓存行(64字节)
@ti.data_oriented
class AlignedData:
    def __init__(self):
        # Taichi自动处理对齐
        self.data = ti.field(ti.f32, shape=(1024, 1024))
        
        # 手动填充避免伪共享
        self.counters = ti.field(ti.i32, shape=16)
        self.padding = ti.field(ti.i32, shape=16)  # 填充到缓存行
```

### 1.5.3 内存访问模式优化

不同的内存访问模式对性能影响巨大：

**合并访问(Coalesced Access)**：
```python
# GPU上的合并访问
@ti.kernel
def coalesced_access():
    # 好：相邻线程访问相邻内存
    for i in range(n):
        arr[i] = i * 2.0
    
    # 差：跨步访问
    for i in range(n):
        arr[i * stride] = i * 2.0
```

**Bank冲突避免**：
```python
# 共享内存bank冲突
@ti.kernel
def bank_conflict_example():
    shared_mem = ti.field(ti.f32, shape=32)
    
    # 冲突：所有线程访问同一bank
    for i in range(32):
        val = shared_mem[i % 4]  # 8路冲突
    
    # 无冲突：添加偏移
    for i in range(32):
        val = shared_mem[(i + i // 4) % 32]
```

**预取策略**：
```python
@ti.kernel
def manual_prefetch():
    # 软件预取提示
    for i in range(n - prefetch_distance):
        # 预取未来要用的数据
        ti.prefetch(arr[i + prefetch_distance])
        # 处理当前数据
        result[i] = complex_computation(arr[i])
```

### 1.5.4 向量化编程技巧

利用SIMD指令提升计算密集型代码性能：

**自动向量化**：
```python
@ti.kernel
def auto_vectorization():
    # Taichi自动向量化简单循环
    for i in range(n):
        c[i] = a[i] + b[i]  # 自动使用SIMD
    
    # 复杂操作可能阻止向量化
    for i in range(n):
        if a[i] > 0:  # 分支阻止向量化
            c[i] = ti.sqrt(a[i])
```

**显式向量操作**：
```python
@ti.kernel
def explicit_vectorization():
    # 使用向量类型
    vec_size = 4
    n_vec = n // vec_size
    
    # 向量化主循环
    for i in range(n_vec):
        vec_a = ti.Vector([a[i*4+j] for j in range(4)])
        vec_b = ti.Vector([b[i*4+j] for j in range(4)])
        vec_c = vec_a + vec_b
        for j in range(4):
            c[i*4+j] = vec_c[j]
    
    # 处理剩余元素
    for i in range(n_vec * vec_size, n):
        c[i] = a[i] + b[i]
```

**向量化友好的算法**：
```python
# 归约操作的向量化
@ti.kernel
def vectorized_reduction() -> ti.f32:
    # 使用多个累加器减少依赖
    sum0 = 0.0
    sum1 = 0.0
    sum2 = 0.0
    sum3 = 0.0
    
    # 4路展开
    for i in range(n // 4):
        sum0 += arr[i * 4 + 0]
        sum1 += arr[i * 4 + 1]
        sum2 += arr[i * 4 + 2]
        sum3 += arr[i * 4 + 3]
    
    # 合并结果
    total = sum0 + sum1 + sum2 + sum3
    
    # 处理剩余
    for i in range((n // 4) * 4, n):
        total += arr[i]
    
    return total
```

## 1.6 面向对象的编程(OOP)

### 1.6.1 Taichi中的类与装饰器

Taichi支持面向对象编程，通过`@ti.data_oriented`装饰器实现：

**基本类定义**：
```python
@ti.data_oriented
class RigidBody:
    def __init__(self, mass, position, velocity):
        # 成员变量可以是Taichi fields
        self.mass = ti.field(ti.f32, shape=())
        self.position = ti.Vector.field(3, dtype=ti.f32, shape=())
        self.velocity = ti.Vector.field(3, dtype=ti.f32, shape=())
        self.force = ti.Vector.field(3, dtype=ti.f32, shape=())
        
        # 初始化
        self.mass[None] = mass
        self.position[None] = position
        self.velocity[None] = velocity
    
    @ti.kernel
    def update_physics(self, dt: ti.f32):
        # kernel可以作为成员函数
        acceleration = self.force[None] / self.mass[None]
        self.velocity[None] += acceleration * dt
        self.position[None] += self.velocity[None] * dt
    
    @ti.func
    def compute_kinetic_energy(self) -> ti.f32:
        # func也可以作为成员函数
        return 0.5 * self.mass[None] * self.velocity[None].norm_sqr()
```

**继承与多态**：
```python
@ti.data_oriented
class Particle:
    def __init__(self, n):
        self.x = ti.Vector.field(3, dtype=ti.f32, shape=n)
        self.v = ti.Vector.field(3, dtype=ti.f32, shape=n)
        self.m = ti.field(ti.f32, shape=n)
    
    @ti.kernel
    def step(self, dt: ti.f32):
        for i in self.x:
            self.x[i] += self.v[i] * dt

@ti.data_oriented
class ChargedParticle(Particle):
    def __init__(self, n):
        super().__init__(n)
        self.charge = ti.field(ti.f32, shape=n)
        self.E_field = ti.Vector.field(3, dtype=ti.f32, shape=n)
    
    @ti.kernel
    def step(self, dt: ti.f32):
        # 重写父类方法
        for i in self.x:
            # 电场力
            force = self.charge[i] * self.E_field[i]
            self.v[i] += force / self.m[i] * dt
            self.x[i] += self.v[i] * dt
```

### 1.6.2 ODOP设计模式

ODOP (Objective Data-Oriented Programming) 结合OOP的抽象能力和DOP的性能优势：

**ODOP原则**：
1. 数据和行为封装在类中
2. 使用SoA布局优化性能
3. 批量操作而非单对象操作
4. 最小化虚函数调用

**ODOP示例**：
```python
@ti.data_oriented
class ParticleSystem:
    def __init__(self, max_particles):
        # SoA布局
        self.position = ti.Vector.field(3, ti.f32, shape=max_particles)
        self.velocity = ti.Vector.field(3, ti.f32, shape=max_particles)
        self.life = ti.field(ti.f32, shape=max_particles)
        self.active = ti.field(ti.i32, shape=max_particles)
        self.count = ti.field(ti.i32, shape=())
        
    @ti.kernel
    def emit_particles(self, num: ti.i32, pos: ti.types.vector(3, ti.f32)):
        old_count = ti.atomic_add(self.count[None], num)
        for i in range(num):
            if old_count + i < self.position.shape[0]:
                idx = old_count + i
                self.position[idx] = pos + ti.Vector([
                    (ti.random() - 0.5) * 0.1,
                    (ti.random() - 0.5) * 0.1,
                    (ti.random() - 0.5) * 0.1
                ])
                self.velocity[idx] = ti.Vector([0.0, 1.0, 0.0])
                self.life[idx] = 1.0
                self.active[idx] = 1
    
    @ti.kernel
    def update(self, dt: ti.f32):
        # 批量更新所有粒子
        for i in range(self.count[None]):
            if self.active[i]:
                # 重力
                self.velocity[i].y -= 9.81 * dt
                # 位置更新
                self.position[i] += self.velocity[i] * dt
                # 生命值衰减
                self.life[i] -= dt * 0.5
                # 死亡检测
                if self.life[i] <= 0:
                    self.active[i] = 0
```

### 1.6.3 太阳系仿真案例

完整的面向对象太阳系仿真示例：

```python
import taichi as ti
import math

ti.init(arch=ti.gpu)

@ti.data_oriented
class CelestialBody:
    def __init__(self, N):
        self.N = N
        self.mass = ti.field(ti.f32, shape=N)
        self.position = ti.Vector.field(3, dtype=ti.f32, shape=N)
        self.velocity = ti.Vector.field(3, dtype=ti.f32, shape=N)
        self.force = ti.Vector.field(3, dtype=ti.f32, shape=N)
        self.radius = ti.field(ti.f32, shape=N)
        self.color = ti.Vector.field(3, dtype=ti.f32, shape=N)
        
    @ti.kernel
    def compute_force(self):
        G = 6.67430e-11  # 引力常数
        
        # 清零力
        for i in range(self.N):
            self.force[i] = ti.Vector([0.0, 0.0, 0.0])
        
        # 计算万有引力
        for i in range(self.N):
            for j in range(self.N):
                if i != j:
                    r = self.position[j] - self.position[i]
                    r_mag = r.norm()
                    if r_mag > 1e-6:  # 避免除零
                        F_mag = G * self.mass[i] * self.mass[j] / (r_mag * r_mag)
                        self.force[i] += F_mag * r.normalized()
    
    @ti.kernel
    def update_verlet(self, dt: ti.f32):
        # Velocity Verlet积分
        for i in range(self.N):
            # 更新位置
            self.position[i] += self.velocity[i] * dt + \
                               0.5 * self.force[i] / self.mass[i] * dt * dt
        
        # 重新计算力
        self.compute_force()
        
        # 更新速度
        for i in range(self.N):
            self.velocity[i] += self.force[i] / self.mass[i] * dt
    
    def initialize_solar_system(self):
        # 太阳
        self.mass[0] = 1.989e30
        self.position[0] = ti.Vector([0.0, 0.0, 0.0])
        self.velocity[0] = ti.Vector([0.0, 0.0, 0.0])
        self.radius[0] = 6.96e8
        self.color[0] = ti.Vector([1.0, 1.0, 0.0])
        
        # 地球
        self.mass[1] = 5.972e24
        self.position[1] = ti.Vector([1.496e11, 0.0, 0.0])
        self.velocity[1] = ti.Vector([0.0, 2.978e4, 0.0])
        self.radius[1] = 6.371e6
        self.color[1] = ti.Vector([0.0, 0.0, 1.0])
        
        # 更多行星...
```

### 1.6.4 代码复用与模块化

良好的OOP设计促进代码复用：

**组件化设计**：
```python
@ti.data_oriented
class Transform:
    """位置和方向组件"""
    def __init__(self, n):
        self.position = ti.Vector.field(3, ti.f32, shape=n)
        self.rotation = ti.Vector.field(4, ti.f32, shape=n)  # 四元数
        self.scale = ti.Vector.field(3, ti.f32, shape=n)

@ti.data_oriented  
class Mesh:
    """网格组件"""
    def __init__(self, vertices, indices):
        self.vertices = ti.Vector.field(3, ti.f32, shape=len(vertices))
        self.indices = ti.field(ti.i32, shape=len(indices))
        
@ti.data_oriented
class GameObject:
    """组合多个组件"""
    def __init__(self, n):
        self.transform = Transform(n)
        self.physics = RigidBody(n)
        self.renderer = MeshRenderer(n)
    
    @ti.kernel
    def update(self, dt: ti.f32):
        # 协调各组件更新
        self.physics.step(dt)
        self.sync_transform()
```

**接口设计**：
```python
class Integrator:
    """积分器接口"""
    @ti.kernel
    def step(self, state: ti.template(), dt: ti.f32):
        raise NotImplementedError

@ti.data_oriented
class EulerIntegrator(Integrator):
    @ti.kernel  
    def step(self, state: ti.template(), dt: ti.f32):
        for i in state.position:
            state.position[i] += state.velocity[i] * dt

@ti.data_oriented
class VerletIntegrator(Integrator):
    def __init__(self, n):
        self.prev_position = ti.Vector.field(3, ti.f32, shape=n)
        
    @ti.kernel
    def step(self, state: ti.template(), dt: ti.f32):
        for i in state.position:
            temp = state.position[i]
            state.position[i] = 2 * state.position[i] - \
                               self.prev_position[i] + \
                               state.force[i] / state.mass[i] * dt * dt
            self.prev_position[i] = temp
```

## 1.7 元编程(MP)

### 1.7.1 模板核函数

Taichi的元编程允许编写高度通用和优化的代码：

**基本模板**：
```python
@ti.kernel
def copy_field(src: ti.template(), dst: ti.template()):
    """通用的field复制函数"""
    for I in ti.grouped(src):
        dst[I] = src[I]

# 可用于任何维度和类型的field
field_2d = ti.field(ti.f32, shape=(100, 100))
field_3d = ti.field(ti.f32, shape=(50, 50, 50))
copy_field(field_2d, another_2d_field)
copy_field(field_3d, another_3d_field)
```

**参数化kernel**：
```python
@ti.kernel
def parametric_kernel(field: ti.template(), 
                     alpha: ti.template(),
                     beta: ti.template()):
    """编译时参数化"""
    for I in ti.grouped(field):
        # alpha和beta在编译时确定
        field[I] = alpha * field[I] + beta

# 为不同参数生成特化版本
parametric_kernel(field, 2.0, 3.0)  # 生成一个版本
parametric_kernel(field, 1.0, 0.0)  # 生成另一个版本
```

### 1.7.2 维度无关编程

编写可以在2D/3D甚至更高维度工作的代码：

```python
@ti.kernel
def laplacian(f: ti.template(), result: ti.template()):
    """计算任意维度的拉普拉斯算子"""
    for I in ti.grouped(f):
        result[I] = -2.0 * ti.static(len(I)) * f[I]
        
        # 静态循环，编译时展开
        for d in ti.static(range(len(I))):
            offset = ti.Vector.unit(len(I), d)
            if I[d] > 0:
                result[I] += f[I - offset]
            else:
                result[I] += f[I]  # 边界条件
                
            if I[d] < f.shape[d] - 1:
                result[I] += f[I + offset]
            else:
                result[I] += f[I]  # 边界条件

# 同一函数用于2D和3D
field_2d = ti.field(ti.f32, shape=(100, 100))
laplacian_2d = ti.field(ti.f32, shape=(100, 100))
laplacian(field_2d, laplacian_2d)

field_3d = ti.field(ti.f32, shape=(50, 50, 50))
laplacian_3d = ti.field(ti.f32, shape=(50, 50, 50))
laplacian(field_3d, laplacian_3d)
```

**维度通用的物理仿真**：
```python
@ti.data_oriented
class GenericParticleSystem:
    def __init__(self, n, dim):
        self.dim = dim
        self.position = ti.Vector.field(dim, ti.f32, shape=n)
        self.velocity = ti.Vector.field(dim, ti.f32, shape=n)
        self.force = ti.Vector.field(dim, ti.f32, shape=n)
        self.mass = ti.field(ti.f32, shape=n)
    
    @ti.kernel
    def compute_pairwise_force(self):
        for i in range(self.position.shape[0]):
            self.force[i] = ti.Vector.zero(ti.f32, self.dim)
            
        for i, j in ti.ndrange(self.position.shape[0], 
                               self.position.shape[0]):
            if i < j:
                r = self.position[j] - self.position[i]
                r_norm = r.norm()
                if r_norm > 1e-6:
                    # Lennard-Jones势
                    f_mag = 24 * (2 / r_norm**14 - 1 / r_norm**8)
                    f = f_mag * r.normalized()
                    self.force[i] += f
                    self.force[j] -= f
```

### 1.7.3 编译时分支与循环展开

利用`ti.static`进行编译时优化：

**静态分支**：
```python
@ti.kernel
def static_branching(field: ti.template(), mode: ti.template()):
    for I in ti.grouped(field):
        # 编译时决定分支
        if ti.static(mode == 'add'):
            field[I] += 1.0
        elif ti.static(mode == 'mul'):
            field[I] *= 2.0
        elif ti.static(mode == 'sin'):
            field[I] = ti.sin(field[I])

# 生成三个不同的kernel版本
static_branching(field, 'add')
static_branching(field, 'mul')
static_branching(field, 'sin')
```

**静态循环展开**：
```python
@ti.kernel
def matrix_vector_multiply(mat: ti.template(), 
                          vec: ti.template(),
                          result: ti.template(),
                          n: ti.template()):
    """编译时已知大小的矩阵向量乘法"""
    for i in range(result.shape[0]):
        result[i] = 0.0
        # 完全展开的循环
        for j in ti.static(range(n)):
            result[i] += mat[i, j] * vec[j]

# 为不同大小生成优化版本
mat4 = ti.Matrix.field(4, 4, ti.f32, shape=())
vec4 = ti.Vector.field(4, ti.f32, shape=())
res4 = ti.Vector.field(4, ti.f32, shape=())
matrix_vector_multiply(mat4, vec4, res4, 4)
```

**递归模板**：
```python
@ti.func
def sum_recursive(arr: ti.template(), l: ti.i32, r: ti.i32) -> ti.f32:
    """编译时递归模板"""
    if ti.static(r - l <= 4):
        # 基础情况：直接求和
        s = 0.0
        for i in range(l, r):
            s += arr[i]
        return s
    else:
        # 递归情况：分治
        mid = (l + r) // 2
        return sum_recursive(arr, l, mid) + sum_recursive(arr, mid, r)
```

### 1.7.4 静态类型推导

Taichi的类型系统支持编译时类型推导：

**自动类型推导**：
```python
@ti.kernel
def type_deduction_example():
    # 自动推导类型
    a = 42  # ti.i32
    b = 3.14  # ti.f32
    c = a + b  # ti.f32 (自动提升)
    
    # 向量类型推导
    v1 = ti.Vector([1, 2, 3])  # ti.types.vector(3, ti.i32)
    v2 = ti.Vector([1.0, 2.0, 3.0])  # ti.types.vector(3, ti.f32)
    v3 = v1 + v2  # ti.types.vector(3, ti.f32)
```

**泛型函数**：
```python
@ti.func
def generic_norm(vec: ti.template()) -> ti.template():
    """自动推导返回类型"""
    norm_sqr = 0
    for i in ti.static(range(vec.n)):
        norm_sqr += vec[i] * vec[i]
    return ti.sqrt(norm_sqr)

# 可用于不同精度
vec_f32 = ti.Vector([1.0, 2.0, 3.0])
norm_f32 = generic_norm(vec_f32)  # 返回f32

vec_f64 = ti.Vector([1.0, 2.0, 3.0], dt=ti.f64)  
norm_f64 = generic_norm(vec_f64)  # 返回f64
```

**类型约束**：
```python
@ti.kernel
def type_constrained(field: ti.types.ndarray(dtype=ti.f32, ndim=2)):
    """明确的类型约束"""
    for i, j in ti.ndrange(field.shape[0], field.shape[1]):
        field[i, j] = ti.sin(field[i, j])

# 只接受2D f32数组
import numpy as np
arr = np.random.rand(100, 100).astype(np.float32)
type_constrained(arr)  # OK

# arr_3d = np.random.rand(10, 10, 10).astype(np.float32)
# type_constrained(arr_3d)  # 编译错误：维度不匹配
```

## 本章小结

本章全面介绍了基于物理的动画和Taichi编程语言的核心概念：

**关键概念回顾**：
1. **物理引擎基础**：理解了物理仿真在游戏、电影和工程中的应用，以及实时性与精确性的权衡
2. **Taichi语言特性**：掌握了kernel/func、张量操作、自动并行化等核心功能
3. **编程范式**：
   - DOP：数据布局优化（AoS vs SoA）、缓存友好设计、向量化
   - OOP：类封装、ODOP模式、组件化设计
   - MP：模板编程、维度通用代码、编译时优化

**重要公式和算法**：
- 稳定性条件：$\Delta t \leq c\sqrt{\frac{m}{k}}$（显式积分器）
- 缓存访问时间：L1(~4 cycles) < L2(~12 cycles) < L3(~40 cycles) < 主存(~200 cycles)
- 并行效率：$E = \frac{S}{p} = \frac{T_1}{p \cdot T_p}$，其中S是加速比，p是处理器数

**性能优化要点**：
1. 数据局部性：行优先访问、分块算法、预取
2. 并行化：避免竞态条件、减少原子操作、负载均衡
3. 向量化：SIMD友好的数据布局、循环展开、避免分支

## 练习题

### 基础题（熟悉材料）

**练习1.1**：编写一个Taichi kernel，计算两个向量场的点积。要求支持任意维度。

<details>
<summary>提示</summary>

使用`ti.template()`和`ti.grouped()`实现维度无关的代码。

</details>

<details>
<summary>参考答案</summary>

```python
@ti.kernel
def vector_field_dot(a: ti.template(), b: ti.template()) -> ti.f32:
    result = 0.0
    for I in ti.grouped(a):
        result += a[I].dot(b[I])
    return result
```

关键点：使用`ti.grouped`遍历任意维度的field，利用向量的`dot`方法计算点积。

</details>

**练习1.2**：比较AoS和SoA布局在粒子系统更新中的性能差异。实现两个版本并测量执行时间。

<details>
<summary>提示</summary>

分别实现只更新位置和更新所有属性两种情况，观察哪种布局更优。

</details>

<details>
<summary>参考答案</summary>

AoS在访问单个粒子的所有属性时更优（如碰撞检测），SoA在批量更新单一属性时更优（如位置更新）。具体性能差异取决于：
- 缓存行大小（通常64字节）
- SIMD向量宽度
- 访问模式的规律性

测量结果通常显示SoA在纯位置更新时快2-4倍。

</details>

**练习1.3**：实现一个线程安全的直方图统计kernel，统计0-255范围内整数的出现频率。

<details>
<summary>提示</summary>

考虑使用原子操作，但要注意优化策略减少冲突。

</details>

<details>
<summary>参考答案</summary>

```python
@ti.kernel
def histogram(data: ti.template(), hist: ti.template()):
    # 清零
    for i in range(256):
        hist[i] = 0
    
    # 统计
    for i in data:
        ti.atomic_add(hist[data[i]], 1)
```

优化版本可以使用线程局部数组减少原子操作频率。

</details>

**练习1.4**：使用Taichi的调试功能找出以下代码的错误：

```python
@ti.kernel
def buggy_kernel(arr: ti.template()):
    for i in range(arr.shape[0]):
        if i > 0:
            arr[i] = arr[i-1] + arr[i+1]  # 有bug
```

<details>
<summary>提示</summary>

考虑边界条件和数组访问范围。

</details>

<details>
<summary>参考答案</summary>

当`i = arr.shape[0] - 1`时，`arr[i+1]`会越界。修正：
```python
@ti.kernel
def fixed_kernel(arr: ti.template()):
    for i in range(arr.shape[0]):
        if i > 0 and i < arr.shape[0] - 1:
            arr[i] = arr[i-1] + arr[i+1]
```

</details>

### 挑战题（深入思考）

**练习1.5**：设计并实现一个高效的2D空间哈希数据结构，支持动态插入和范围查询。要求：
- 支持百万级粒子
- 查询半径内的所有邻居
- 处理非均匀分布

<details>
<summary>提示</summary>

考虑使用多级网格或自适应网格大小，处理稀疏和密集区域。

</details>

<details>
<summary>参考答案</summary>

关键设计要点：
1. 使用质数表大小减少哈希冲突
2. 链表处理冲突，但要考虑内存局部性
3. 动态调整网格大小基于局部密度
4. 使用Morton编码保持空间局部性

性能优化：
- 预分配内存池避免动态分配
- 使用原子操作并行构建
- 定期重建以保持效率

</details>

**练习1.6**：实现一个通用的归约（reduction）框架，支持任意二元操作和数据类型。要求在GPU上达到接近峰值带宽。

<details>
<summary>提示</summary>

考虑warp-level primitives、shared memory使用、多级归约策略。

</details>

<details>
<summary>参考答案</summary>

高效归约的关键技术：
1. **Warp内归约**：利用warp shuffle指令避免共享内存
2. **多级策略**：warp级→block级→grid级
3. **负载均衡**：每个线程处理多个元素
4. **内存访问**：合并访问模式，避免bank冲突

理论分析：
- 带宽利用率 = 实际带宽 / 峰值带宽
- 目标：达到80%以上的带宽利用率
- 瓶颈：最后几步的并行度下降

</details>

**练习1.7**：使用元编程技术实现一个编译时计算的快速傅里叶变换（FFT）。支持2的幂次大小。

<details>
<summary>提示</summary>

使用`ti.static`展开蝶形操作，预计算旋转因子。

</details>

<details>
<summary>参考答案</summary>

元编程FFT的要点：
1. 编译时确定FFT大小和迭代次数
2. 静态展开所有循环层级
3. 预计算twiddle factors为常量
4. 使用位反转的查表实现

性能考虑：
- 寄存器压力vs循环开销权衡
- 共享内存bank冲突避免
- 适合小尺寸FFT（≤1024点）

</details>

**练习1.8**：设计一个物理引擎的整体架构，整合本章所学的所有概念。包括：
- 模块化设计（碰撞检测、约束求解、积分器）
- 数据结构选择（粒子、刚体、约束）
- 并行化策略
- 可扩展性考虑

<details>
<summary>提示</summary>

参考现代物理引擎如Bullet、PhysX的设计，但要考虑Taichi的特点。

</details>

<details>
<summary>参考答案</summary>

架构设计要点：

1. **核心模块**：
   - BroadPhase：空间划分，快速剔除
   - NarrowPhase：精确碰撞检测
   - ConstraintSolver：约束求解（PGS/SI）
   - Integrator：时间积分（插件式）

2. **数据组织**：
   - ECS模式：实体-组件-系统
   - SoA布局：性能优化
   - 双缓冲：避免竞态

3. **并行策略**：
   - 空间分区并行
   - 图着色避免冲突
   - 任务图调度

4. **扩展性**：
   - 插件式力场
   - 自定义约束
   - 脚本绑定

</details>

## 常见陷阱与错误

### 1. 并行编程陷阱

**陷阱1**：假设kernel内的执行顺序
```python
# 错误：依赖执行顺序
@ti.kernel
def wrong():
    for i in range(n):
        if i > 0:
            a[i] = a[i-1] + 1  # 并行执行，顺序不确定！
```

**解决**：使用串行循环或改变算法

**陷阱2**：过度使用原子操作
```python
# 低效：大量原子操作
@ti.kernel
def inefficient():
    for i in range(1000000):
        ti.atomic_add(sum[0], data[i])  # 严重的串行化
```

**解决**：使用归约或线程局部累加

### 2. 内存访问陷阱

**陷阱3**：列主序访问二维数组
```python
# 缓存不友好
for j in range(m):
    for i in range(n):
        process(a[i, j])  # 跨步访问
```

**解决**：调整循环顺序或转置数据

**陷阱4**：伪共享
```python
# 多线程写入相邻内存
counters = ti.field(ti.i32, shape=num_threads)
# 不同线程的counter可能在同一缓存行
```

**解决**：添加填充或使用线程局部变量

### 3. 类型系统陷阱

**陷阱5**：整数除法精度损失
```python
# 意外的整数除法
result = 3 / 2  # 在kernel中可能是1而非1.5
```

**解决**：显式使用浮点数或转换类型

### 4. 调试陷阱

**陷阱6**：GPU调试时忘记同步
```python
kernel()
print(field[0])  # 可能读到旧值！
```

**解决**：添加`ti.sync()`确保GPU完成

## 最佳实践检查清单

### 设计阶段
- [ ] 选择合适的数据布局（AoS/SoA/AoSoA）
- [ ] 识别并行化机会和数据依赖
- [ ] 考虑缓存友好的算法设计
- [ ] 规划内存使用和数据流

### 实现阶段
- [ ] 使用`ti.template()`编写通用代码
- [ ] 最小化原子操作使用
- [ ] 确保内存访问模式合理
- [ ] 适当使用`ti.static`优化

### 优化阶段
- [ ] 使用性能分析工具定位瓶颈
- [ ] 考虑向量化和循环展开
- [ ] 平衡计算与内存访问
- [ ] 验证并行正确性

### 调试阶段
- [ ] 启用调试模式进行边界检查
- [ ] 使用小数据集验证正确性
- [ ] 添加断言检查不变量
- [ ] 记录性能基准用于回归测试

### 代码质量
- [ ] 遵循ODOP设计模式
- [ ] 模块化和可测试性
- [ ] 文档化性能关键决策
- [ ] 考虑跨平台兼容性
