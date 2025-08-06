# 第九章：高性能计算

本章深入探讨物理引擎的性能优化技术，从底层硬件架构到高层算法优化。我们将学习如何充分利用现代处理器的并行计算能力，设计高效的数据结构，以及实现可扩展的并行算法。通过本章的学习，读者将掌握构建高性能物理仿真系统所需的核心技术。

## 9.1 处理器微架构

### 9.1.1 CPU流水线与超标量

现代CPU采用流水线技术将指令执行分解为多个阶段，典型的五级流水线包括：

1. **取指(Fetch)**：从内存读取指令
2. **译码(Decode)**：解析指令操作码和操作数
3. **执行(Execute)**：在ALU中执行运算
4. **访存(Memory)**：读写内存（如需要）
5. **回写(Write Back)**：将结果写回寄存器

超标量处理器能够在一个时钟周期内发射多条指令，现代CPU通常具有4-6个发射端口。指令级并行(ILP)的关键在于：

- **数据依赖性**：后续指令依赖前面指令的结果
- **控制依赖性**：分支指令影响程序流
- **结构冲突**：多条指令竞争相同的执行单元

为了提高ILP，编译器和程序员需要：
- 展开循环减少分支
- 重排指令减少数据依赖
- 使用SIMD指令增加吞吐量

### 9.1.2 SIMD指令集

单指令多数据(SIMD)是提升数值计算性能的关键技术。主要的SIMD扩展包括：

**SSE/AVX系列**：
- SSE：128位寄存器，4个float或2个double
- AVX：256位寄存器，8个float或4个double  
- AVX-512：512位寄存器，16个float或8个double

**常用SIMD操作**：
```
// 加法：c = a + b
__m256 c = _mm256_add_ps(a, b);

// 乘加：d = a * b + c
__m256 d = _mm256_fmadd_ps(a, b, c);

// 广播：将标量复制到向量所有位置
__m256 vec = _mm256_set1_ps(scalar);
```

**向量化考虑**：
- 数据对齐：使用`_mm_malloc`确保32字节对齐
- 避免分支：使用blend指令代替if-else
- 减少gather/scatter：连续内存访问效率最高

### 9.1.3 分支预测与投机执行

现代CPU使用分支预测器猜测条件跳转的方向，典型预测准确率达95%以上。分支预测失败会导致流水线清空，损失10-20个时钟周期。

**优化策略**：
1. **减少分支**：使用条件移动指令或位运算
2. **分支提示**：使用`__builtin_expect`提供预测信息
3. **循环展开**：减少循环边界检查
4. **模板化**：将运行时分支转为编译时决策

**投机执行风险**：
- Spectre/Meltdown等安全漏洞
- 需要考虑侧信道攻击防护

### 9.1.4 乱序执行

乱序执行允许CPU重排指令顺序以隐藏延迟：

**Tomasulo算法核心组件**：
- **保留站**：缓存等待执行的指令
- **重排缓冲区(ROB)**：维护程序顺序
- **寄存器重命名**：消除伪依赖

**性能影响因素**：
- 指令窗口大小（典型200-300条）
- 执行端口数量和类型
- 内存消歧(memory disambiguation)能力

**编程建议**：
- 增加指令间距离，给CPU更多重排空间
- 使用restrict关键字帮助编译器优化
- 避免长依赖链

## 9.2 内存层级

### 9.2.1 缓存结构与大小

现代CPU的缓存层次结构：

| 级别 | 大小 | 延迟 | 带宽 |
|------|------|------|------|
| L1-D | 32KB | 4-5 cycles | ~100 GB/s |
| L1-I | 32KB | 4-5 cycles | ~100 GB/s |
| L2 | 256KB-1MB | 12-15 cycles | ~50 GB/s |
| L3 | 8MB-32MB | 30-40 cycles | ~30 GB/s |
| DRAM | 16GB+ | 100-300 cycles | ~20 GB/s |

**缓存组织**：
- **组相联**：L1通常8路，L2/L3可达16-20路
- **缓存行**：64字节，是缓存传输的基本单位
- **包含性**：L3通常包含L1/L2的内容

**缓存优化原则**：
1. **空间局部性**：连续访问相邻数据
2. **时间局部性**：重复访问相同数据
3. **工作集大小**：保持热数据在缓存内

### 9.2.2 缓存行与伪共享

伪共享(False Sharing)是多核编程的常见性能陷阱：

```cpp
struct BadLayout {
    int thread1_data;  // 线程1写
    int thread2_data;  // 线程2写，但在同一缓存行！
};

struct GoodLayout {
    alignas(64) int thread1_data;  // 独占缓存行
    alignas(64) int thread2_data;  // 独占缓存行
};
```

**检测方法**：
- 使用perf监控缓存失效率
- Intel VTune的False Sharing检测器
- 手动padding或使用`alignas`

**缓存行优化技巧**：
- 将只读数据分离到独立缓存行
- 使用本地变量累积，最后写回
- 考虑数据布局的缓存友好性

### 9.2.3 TLB与页表

转换查找缓冲器(TLB)缓存虚拟地址到物理地址的映射：

**TLB层次**：
- L1 DTLB：64-128项，4KB页
- L1 ITLB：64-128项，指令页
- L2 TLB：512-1536项，统一
- 大页支持：2MB/1GB页减少TLB压力

**TLB优化**：
1. **使用大页**：`madvise(MADV_HUGEPAGE)`
2. **减少工作集**：避免稀疏访问模式
3. **页面着色**：控制物理页分配减少冲突

**页表遍历开销**：
- 4级页表需要4次内存访问
- 使用大页减少页表级数
- 考虑NUMA下的页面放置

### 9.2.4 NUMA架构

非统一内存访问(NUMA)是多插槽系统的标准架构：

**NUMA特征**：
- 每个CPU有本地内存控制器
- 访问远程内存延迟高50-100%
- QPI/UPI互连带宽有限

**NUMA感知编程**：
```cpp
// 绑定线程到NUMA节点
numa_run_on_node(node_id);

// 在指定节点分配内存
void* ptr = numa_alloc_onnode(size, node_id);

// 内存访问模式优化
#pragma omp parallel for schedule(static)
for (int i = 0; i < n; i++) {
    // 确保线程访问本地数据
    process_local_data(i);
}
```

**NUMA优化策略**：
- First-touch策略：首次访问的线程所在节点分配页面
- 显式内存绑定：使用`numactl`或API
- 避免远程访问：数据分区对齐NUMA拓扑

## 9.3 性能分析与优化

### 9.3.1 Roofline模型

Roofline模型将程序性能上限可视化为计算能力和内存带宽的函数：

**性能上限计算**：
$$\text{Performance} = \min(\text{Peak FLOPS}, \text{Arithmetic Intensity} \times \text{Bandwidth})$$

其中算术强度(Arithmetic Intensity) = FLOPs / Bytes

**典型算术强度**：
- DAXPY (y = ax + y): 2 FLOPs / 24 Bytes = 0.083
- 矩阵乘法: O(n³) FLOPs / O(n²) Bytes = O(n)
- 模板计算: 取决于模板大小和重用

**优化方向判断**：
- AI < 机器平衡点：内存带宽受限，需要数据重用
- AI > 机器平衡点：计算受限，需要向量化/并行化

### 9.3.2 带宽vs计算瓶颈

识别性能瓶颈的关键指标：

**内存带宽受限特征**：
- CPU利用率低（等待内存）
- L3缓存未命中率高
- 内存控制器接近饱和

**计算受限特征**：
- CPU利用率接近100%
- 指令吞吐量接近峰值
- 缓存命中率高

**平衡优化策略**：
1. **提高算术强度**：循环阻塞、数据重用
2. **减少内存流量**：压缩数据结构、精度权衡
3. **隐藏内存延迟**：预取、软件流水线

### 9.3.3 性能计数器

硬件性能计数器提供详细的执行信息：

**关键性能事件**：
```bash
# 使用perf统计缓存失效
perf stat -e cache-misses,cache-references ./program

# 监控分支预测
perf stat -e branch-misses,branch-instructions ./program

# 内存带宽测量
perf stat -e uncore_imc/data_reads/,uncore_imc/data_writes/ ./program
```

**常用计数器**：
- Instructions Per Cycle (IPC)
- Cache hit/miss rates
- TLB miss rates
- Branch misprediction rate
- Memory bandwidth utilization

**分析工具链**：
- Linux perf：轻量级，内置内核
- Intel VTune：详细的微架构分析
- AMD uProf：AMD平台优化
- NVIDIA Nsight：GPU性能分析

### 9.3.4 Profile工具使用

**采样分析(Sampling Profiler)**：
```bash
# 记录性能数据
perf record -g ./simulation

# 生成火焰图
perf script | flamegraph.pl > flame.svg

# 查看热点函数
perf report --stdio
```

**追踪分析(Tracing Profiler)**：
- 精确的函数调用时序
- 更高的开销
- 适合分析复杂交互

**优化工作流**：
1. **识别热点**：找到占用时间最多的函数
2. **分析瓶颈**：确定是计算、内存还是同步
3. **针对优化**：应用相应的优化技术
4. **验证效果**：对比优化前后性能

## 9.4 并行编程模型

### 9.4.1 共享内存并行(OpenMP)

OpenMP是共享内存并行编程的事实标准，通过编译器指令实现并行化：

**基本并行结构**：
```cpp
#pragma omp parallel for
for (int i = 0; i < n; i++) {
    process(data[i]);  // 自动分配给多个线程
}
```

**高级OpenMP特性**：

1. **调度策略**：
   - `static`：固定大小块分配，开销最小
   - `dynamic`：动态分配，适合负载不均
   - `guided`：递减块大小，平衡负载和开销
   - `runtime`：运行时通过环境变量决定

2. **归约操作**：
```cpp
double sum = 0.0;
#pragma omp parallel for reduction(+:sum)
for (int i = 0; i < n; i++) {
    sum += array[i];
}
```

3. **同步原语**：
   - `#pragma omp critical`：临界区
   - `#pragma omp atomic`：原子操作
   - `#pragma omp barrier`：栅栏同步

**OpenMP优化技巧**：
- 避免false sharing：使用`firstprivate`
- 减少同步开销：批量处理减少critical区域
- NUMA感知：使用`proc_bind`和`places`子句

### 9.4.2 分布式内存(MPI)

消息传递接口(MPI)是分布式内存并行的标准：

**基本通信模式**：
```cpp
// 点对点通信
MPI_Send(buffer, count, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
MPI_Recv(buffer, count, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);

// 集合通信
MPI_Bcast(buffer, count, MPI_DOUBLE, root, MPI_COMM_WORLD);
MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
```

**非阻塞通信**：
```cpp
MPI_Request request;
MPI_Isend(buffer, count, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD, &request);
// 计算与通信重叠
do_computation();
MPI_Wait(&request, &status);
```

**MPI优化策略**：
1. **消息聚合**：减少小消息的开销
2. **通信-计算重叠**：使用非阻塞通信
3. **拓扑感知**：匹配通信模式与网络拓扑
4. **负载均衡**：动态任务分配

### 9.4.3 任务并行vs数据并行

**数据并行**：
- 相同操作应用于不同数据
- SIMD、GPU天然适合
- 易于负载均衡
- 例：矩阵运算、图像处理

**任务并行**：
- 不同任务同时执行
- 依赖关系复杂
- 负载均衡困难
- 例：流水线、DAG调度

**混合并行模式**：
```cpp
#pragma omp parallel sections
{
    #pragma omp section
    { compute_forces(); }     // 任务1
    
    #pragma omp section  
    { update_positions(); }   // 任务2
}

#pragma omp parallel for     // 数据并行
for (int i = 0; i < n; i++) {
    apply_constraints(particles[i]);
}
```

### 9.4.4 负载均衡策略

**静态负载均衡**：
- 编译时或启动时分配
- 开销小但灵活性差
- 适合规则计算

**动态负载均衡**：
1. **工作池模式**：
```cpp
while (true) {
    Task* task = get_next_task();  // 原子获取
    if (!task) break;
    process_task(task);
}
```

2. **工作窃取(Work Stealing)**：
- 每个线程维护本地队列
- 空闲线程从忙碌线程窃取任务
- Intel TBB、Cilk Plus实现

3. **层次分解**：
- 递归划分问题
- 适合树形结构
- 自然负载均衡

**负载均衡度量**：
$$\text{Efficiency} = \frac{\text{Average Load}}{\text{Maximum Load}}$$

目标是使效率接近1.0

## 9.5 GPU编程

### 9.5.1 GPU架构概述

GPU采用SIMT(Single Instruction Multiple Thread)架构，与CPU的关键区别：

**架构特点**：
- 数千个简单核心 vs 几十个复杂核心
- 深度多线程隐藏内存延迟
- 高内存带宽（~1TB/s vs ~100GB/s）
- 有限的分支能力

**GPU内存层次**：
| 内存类型 | 大小 | 延迟 | 带宽 | 作用域 |
|---------|------|------|------|--------|
| 寄存器 | 64KB/SM | 1 cycle | 极高 | 线程私有 |
| 共享内存 | 48-96KB/SM | ~30 cycles | ~2TB/s | Block共享 |
| L1缓存 | 24-48KB/SM | ~30 cycles | ~1TB/s | 自动管理 |
| L2缓存 | 4-6MB | ~200 cycles | ~2TB/s | 全局共享 |
| 全局内存 | 8-80GB | ~400 cycles | ~1TB/s | 全局可见 |

### 9.5.2 线程组织

CUDA/OpenCL的三级线程层次：

**Grid → Block → Thread**：
```cuda
// Kernel启动配置
dim3 blockDim(256);  // 每个block 256线程
dim3 gridDim((n + 255) / 256);  // 足够的blocks
kernel<<<gridDim, blockDim>>>(data, n);
```

**Warp执行模型**：
- 32个线程组成一个warp
- Warp内线程锁步执行
- 分支导致线程发散(divergence)
- 合并内存访问提高带宽利用

**线程同步**：
```cuda
__syncthreads();  // Block内同步
__syncwarp();     // Warp内同步（隐式）
// Grid级同步需要kernel边界
```

### 9.5.3 内存合并访问

内存合并是GPU性能的关键：

**合并条件**：
1. 连续的线程访问连续的地址
2. 起始地址对齐到128字节
3. 访问大小为4、8或16字节

**合并访问模式**：
```cuda
// Good: 合并访问
data[threadIdx.x + blockIdx.x * blockDim.x] = value;

// Bad: 跨步访问
data[threadIdx.x * stride] = value;

// Bad: 随机访问  
data[indices[threadIdx.x]] = value;
```

**优化技巧**：
- 使用共享内存转置数据
- Structure of Arrays (SoA)布局
- 纹理内存处理非规则访问

### 9.5.4 占用率优化

占用率(Occupancy)表示活跃warp数与最大warp数的比例：

**影响因素**：
1. **寄存器使用**：每线程寄存器数×线程数 ≤ 总寄存器数
2. **共享内存**：每block共享内存 ≤ 总共享内存
3. **Block大小**：必须是warp大小的倍数

**优化策略**：
```cuda
// 使用launch bounds限制寄存器
__global__ __launch_bounds__(256, 4)
void kernel() {
    // 编译器将优化到指定的占用率
}
```

**动态并行**：
- 从kernel内启动新kernel
- 适合不规则并行
- 注意启动开销

**占用率计算**：
```python
def calc_occupancy(regs_per_thread, shared_per_block, threads_per_block):
    max_warps_per_sm = 64  # 架构相关
    max_regs_per_sm = 65536
    max_shared_per_sm = 49152
    
    warps_limited_by_regs = max_regs_per_sm // (regs_per_thread * 32)
    blocks_limited_by_shared = max_shared_per_sm // shared_per_block
    warps_limited_by_shared = blocks_limited_by_shared * (threads_per_block // 32)
    
    active_warps = min(warps_limited_by_regs, warps_limited_by_shared)
    return active_warps / max_warps_per_sm
```

## 9.6 稀疏数据结构

### 9.6.1 SPGrid与OpenVDB

**SPGrid特点**：
- 利用虚拟内存的稀疏性
- 页面级别的存在性跟踪
- 适合均匀稀疏的数据

**实现原理**：
```cpp
template<int dim>
class SPGrid {
    static constexpr int page_size = 4096;
    static constexpr int block_size = (dim == 2) ? 4 : 8;
    
    uint64_t* page_table;  // 页表跟踪
    T* data;               // 虚拟内存池
    
    size_t offset(const Vec<dim>& idx) {
        // Morton编码保持空间局部性
        return morton_encode(idx);
    }
};
```

**OpenVDB特点**：
- 层次化B+树结构
- 支持极度稀疏数据
- 工业标准（电影特效）

**树结构**：
```
Root (dynamic hash map)
  └── Internal Node (32³ or 16³)
       └── Internal Node (16³ or 8³)  
            └── Leaf Node (8³)
```

### 9.6.2 分层稀疏网格

多分辨率表示的优势：

**自适应细化**：
```cpp
struct HierarchicalGrid {
    struct Node {
        union {
            Node* children[8];  // 内部节点
            T data[8*8*8];      // 叶子节点
        };
        uint8_t mask;           // 子节点存在性
        
        bool is_leaf() const { return mask & 0x80; }
    };
    
    void refine(Node* node, const Vec3& pos) {
        if (need_refinement(node, pos)) {
            split_node(node);
        }
    }
};
```

**内存池管理**：
- 预分配大块内存
- 自定义分配器减少碎片
- 延迟释放提高重用率

### 9.6.3 位掩码与指针结构

**位掩码优化**：
```cpp
struct SparseMask {
    uint64_t mask[BLOCKS/64];  // 每bit表示一个block
    
    bool is_active(int idx) {
        return mask[idx/64] & (1ULL << (idx%64));
    }
    
    void activate(int idx) {
        mask[idx/64] |= (1ULL << (idx%64));
    }
    
    int popcount() {  // 活跃block数
        int count = 0;
        for (auto m : mask) {
            count += __builtin_popcountll(m);
        }
        return count;
    }
};
```

**指针压缩**：
- 使用相对偏移代替绝对地址
- 32位索引处理大部分情况
- 指针打包提高缓存效率

### 9.6.4 动态拓扑更新

**增量更新策略**：
1. **延迟分配**：首次写入时分配
2. **批量更新**：累积修改减少开销
3. **垃圾回收**：定期清理空块

**并发更新**：
```cpp
class ConcurrentSparseGrid {
    struct Block {
        std::atomic<uint32_t> version;
        T data[BLOCK_SIZE];
    };
    
    void update(const Vec3& pos, T value) {
        Block* block = get_or_create_block(pos);
        uint32_t ver = block->version.load();
        
        // 乐观并发控制
        block->data[local_idx(pos)] = value;
        
        if (!block->version.compare_exchange_weak(ver, ver+1)) {
            // 重试或使用其他策略
        }
    }
};
```

## 9.7 算法优化技术

### 9.7.1 循环优化

循环是数值计算的核心，优化循环对性能至关重要：

**循环展开(Loop Unrolling)**：
```cpp
// 原始循环
for (int i = 0; i < n; i++) {
    sum += a[i] * b[i];
}

// 4倍展开
for (int i = 0; i < n-3; i += 4) {
    sum += a[i] * b[i] + a[i+1] * b[i+1] + 
           a[i+2] * b[i+2] + a[i+3] * b[i+3];
}
// 处理剩余元素
for (int i = n - n%4; i < n; i++) {
    sum += a[i] * b[i];
}
```

**循环阻塞(Loop Blocking/Tiling)**：
```cpp
// 矩阵乘法的缓存优化版本
const int TILE = 64;  // 选择适合L1缓存的块大小

for (int ii = 0; ii < n; ii += TILE) {
    for (int jj = 0; jj < n; jj += TILE) {
        for (int kk = 0; kk < n; kk += TILE) {
            // 处理一个TILE×TILE的块
            for (int i = ii; i < min(ii+TILE, n); i++) {
                for (int j = jj; j < min(jj+TILE, n); j++) {
                    double sum = C[i][j];
                    for (int k = kk; k < min(kk+TILE, n); k++) {
                        sum += A[i][k] * B[k][j];
                    }
                    C[i][j] = sum;
                }
            }
        }
    }
}
```

**循环融合(Loop Fusion)**：
```cpp
// 分离的循环
for (int i = 0; i < n; i++) {
    a[i] = b[i] + c[i];
}
for (int i = 0; i < n; i++) {
    d[i] = a[i] * e[i];
}

// 融合后：减少内存流量
for (int i = 0; i < n; i++) {
    a[i] = b[i] + c[i];
    d[i] = a[i] * e[i];  // a[i]还在缓存中
}
```

**循环交换(Loop Interchange)**：
```cpp
// 原始：列优先访问（缓存不友好）
for (int j = 0; j < n; j++) {
    for (int i = 0; i < n; i++) {
        process(A[i][j]);
    }
}

// 交换后：行优先访问（缓存友好）
for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
        process(A[i][j]);
    }
}
```

### 9.7.2 数据重排与预取

**AoS到SoA转换**：
```cpp
// Array of Structures (AoS) - 缓存不友好
struct Particle_AoS {
    float x, y, z;
    float vx, vy, vz;
};
Particle_AoS particles[N];

// Structure of Arrays (SoA) - SIMD友好
struct Particles_SoA {
    float x[N], y[N], z[N];
    float vx[N], vy[N], vz[N];
};
```

**软件预取**：
```cpp
void process_with_prefetch(float* data, int n) {
    const int PREFETCH_DISTANCE = 8;
    
    for (int i = 0; i < n; i++) {
        // 预取未来的数据
        if (i + PREFETCH_DISTANCE < n) {
            __builtin_prefetch(&data[i + PREFETCH_DISTANCE], 0, 3);
        }
        
        // 处理当前数据
        expensive_computation(data[i]);
    }
}
```

**数据打包优化**：
```cpp
// 压缩存储减少内存带宽
struct CompressedParticle {
    // 位置使用相对坐标（16位整数）
    int16_t dx, dy, dz;
    // 速度使用半精度浮点
    half vx, vy, vz;
};

// 解压时转换回全精度
void decompress(CompressedParticle& cp, Particle& p, Vec3 base) {
    p.x = base.x + cp.dx * SCALE;
    p.y = base.y + cp.dy * SCALE;
    p.z = base.z + cp.dz * SCALE;
    p.vx = float(cp.vx);
    p.vy = float(cp.vy);
    p.vz = float(cp.vz);
}
```

### 9.7.3 向量化技巧

**手动向量化with Intrinsics**：
```cpp
// 标量版本
void add_scalar(float* a, float* b, float* c, int n) {
    for (int i = 0; i < n; i++) {
        c[i] = a[i] + b[i];
    }
}

// AVX向量化版本
void add_avx(float* a, float* b, float* c, int n) {
    int i = 0;
    // 主循环：8个元素一组
    for (; i < n - 7; i += 8) {
        __m256 va = _mm256_load_ps(&a[i]);
        __m256 vb = _mm256_load_ps(&b[i]);
        __m256 vc = _mm256_add_ps(va, vb);
        _mm256_store_ps(&c[i], vc);
    }
    // 处理剩余元素
    for (; i < n; i++) {
        c[i] = a[i] + b[i];
    }
}
```

**数据对齐**：
```cpp
// 确保32字节对齐（AVX需要）
float* aligned_alloc_float(size_t n) {
    void* ptr = nullptr;
    if (posix_memalign(&ptr, 32, n * sizeof(float)) != 0) {
        return nullptr;
    }
    return static_cast<float*>(ptr);
}

// 使用aligned属性
struct alignas(32) AlignedData {
    float values[8];
};
```

**避免Gather/Scatter**：
```cpp
// Bad: gather操作性能差
__m256 gather_bad(float* data, int* indices) {
    return _mm256_i32gather_ps(data, _mm256_load_si256((__m256i*)indices), 4);
}

// Good: 重组数据避免gather
void transpose_block_8x8(float* src, float* dst, int stride) {
    __m256 row0 = _mm256_load_ps(src + 0*stride);
    __m256 row1 = _mm256_load_ps(src + 1*stride);
    // ... 加载所有行
    
    // 使用shuffle和blend进行转置
    // ... 转置操作
    
    _mm256_store_ps(dst + 0*8, col0);
    _mm256_store_ps(dst + 1*8, col1);
    // ... 存储所有列
}
```

### 9.7.4 通信隐藏

**计算与通信重叠**：
```cpp
// MPI非阻塞通信示例
void overlap_compute_comm(float* local_data, float* halo_data, int n) {
    MPI_Request requests[4];
    
    // 启动非阻塞接收
    MPI_Irecv(halo_left, size, MPI_FLOAT, left_rank, 0, 
              MPI_COMM_WORLD, &requests[0]);
    MPI_Irecv(halo_right, size, MPI_FLOAT, right_rank, 1, 
              MPI_COMM_WORLD, &requests[1]);
    
    // 启动非阻塞发送
    MPI_Isend(send_left, size, MPI_FLOAT, left_rank, 1, 
              MPI_COMM_WORLD, &requests[2]);
    MPI_Isend(send_right, size, MPI_FLOAT, right_rank, 0, 
              MPI_COMM_WORLD, &requests[3]);
    
    // 计算内部区域（不依赖halo）
    compute_interior(local_data, n);
    
    // 等待通信完成
    MPI_Waitall(4, requests, MPI_STATUSES_IGNORE);
    
    // 计算边界区域（使用halo数据）
    compute_boundary(local_data, halo_left, halo_right, n);
}
```

**双缓冲Pipeline**：
```cpp
template<typename T>
class DoubleBuffer {
    T* buffers[2];
    int current = 0;
    
public:
    void swap() { current = 1 - current; }
    T* read_buffer() { return buffers[current]; }
    T* write_buffer() { return buffers[1 - current]; }
};

void pipeline_processing() {
    DoubleBuffer<Data> buffer;
    
    // 启动第一次读取
    async_read(buffer.write_buffer());
    
    for (int i = 0; i < iterations; i++) {
        buffer.swap();
        
        // 并行：处理当前数据 + 读取下一批
        auto future = std::async(std::launch::async, [&]() {
            async_read(buffer.write_buffer());
        });
        
        process(buffer.read_buffer());
        
        future.wait();
    }
}
```

**GPU异步传输**：
```cuda
void multi_stream_processing(float* h_data, float* d_data, int n, int chunks) {
    const int chunk_size = n / chunks;
    cudaStream_t streams[chunks];
    
    // 创建流
    for (int i = 0; i < chunks; i++) {
        cudaStreamCreate(&streams[i]);
    }
    
    // 异步处理
    for (int i = 0; i < chunks; i++) {
        int offset = i * chunk_size;
        
        // 异步拷贝到GPU
        cudaMemcpyAsync(d_data + offset, h_data + offset, 
                       chunk_size * sizeof(float), 
                       cudaMemcpyHostToDevice, streams[i]);
        
        // 在流中启动kernel
        process_kernel<<<blocks, threads, 0, streams[i]>>>
                      (d_data + offset, chunk_size);
        
        // 异步拷贝回CPU
        cudaMemcpyAsync(h_data + offset, d_data + offset, 
                       chunk_size * sizeof(float), 
                       cudaMemcpyDeviceToHost, streams[i]);
    }
    
    // 同步所有流
    for (int i = 0; i < chunks; i++) {
        cudaStreamSynchronize(streams[i]);
        cudaStreamDestroy(streams[i]);
    }
}
```

## 9.8 多GPU并行

### 9.8.1 域分解策略

**1D分解**：
```cpp
struct Domain1D {
    int global_size;
    int num_gpus;
    int gpu_id;
    
    int local_size() const {
        return (global_size + num_gpus - 1) / num_gpus;
    }
    
    int local_start() const {
        return gpu_id * local_size();
    }
    
    int local_end() const {
        return min(local_start() + local_size(), global_size);
    }
};
```

**2D/3D分解**：
```cpp
struct Domain3D {
    int3 global_size;
    int3 gpu_grid;    // GPU的3D排列
    int3 gpu_coords;  // 当前GPU的坐标
    
    int3 local_size() const {
        return make_int3(
            (global_size.x + gpu_grid.x - 1) / gpu_grid.x,
            (global_size.y + gpu_grid.y - 1) / gpu_grid.y,
            (global_size.z + gpu_grid.z - 1) / gpu_grid.z
        );
    }
    
    // 计算表面积/体积比
    float surface_volume_ratio() const {
        int3 ls = local_size();
        int surface = 2 * (ls.x*ls.y + ls.y*ls.z + ls.z*ls.x);
        int volume = ls.x * ls.y * ls.z;
        return float(surface) / volume;
    }
};
```

**最优分解选择**：
- 最小化通信量：减小表面积/体积比
- 考虑问题的各向异性
- 匹配硬件拓扑（如NVLink连接）

### 9.8.2 Halo区域交换

**Ghost Cell实现**：
```cpp
class HaloExchange {
    int3 local_size;
    int3 halo_width;
    int3 total_size;  // local + 2*halo
    
    float* data;
    
    // 6个方向的邻居GPU
    int neighbors[6];  // -x, +x, -y, +y, -z, +z
    
    void exchange_halos() {
        // X方向交换
        exchange_x_direction();
        
        // Y方向交换（包含X方向的halo）
        exchange_y_direction();
        
        // Z方向交换（包含X,Y方向的halo）
        exchange_z_direction();
    }
    
    void exchange_x_direction() {
        size_t yz_slice_size = halo_width.x * local_size.y * local_size.z;
        
        // 准备发送缓冲区
        pack_x_minus(send_buffer_xm);
        pack_x_plus(send_buffer_xp);
        
        // 非阻塞P2P传输
        cudaMemcpyPeerAsync(recv_xm, gpu_xm, send_xp, gpu_id, 
                           yz_slice_size * sizeof(float), stream_xm);
        cudaMemcpyPeerAsync(recv_xp, gpu_xp, send_xm, gpu_id, 
                           yz_slice_size * sizeof(float), stream_xp);
        
        // 同步
        cudaStreamSynchronize(stream_xm);
        cudaStreamSynchronize(stream_xp);
        
        // 解包到halo区域
        unpack_x_minus(recv_buffer_xm);
        unpack_x_plus(recv_buffer_xp);
    }
};
```

**GPUDirect P2P优化**：
```cpp
void enable_peer_access(int num_gpus) {
    for (int i = 0; i < num_gpus; i++) {
        cudaSetDevice(i);
        for (int j = 0; j < num_gpus; j++) {
            if (i != j) {
                int can_access;
                cudaDeviceCanAccessPeer(&can_access, i, j);
                if (can_access) {
                    cudaDeviceEnablePeerAccess(j, 0);
                }
            }
        }
    }
}
```

**重叠计算和通信**：
```cpp
void timestep_with_overlap() {
    // 1. 启动halo交换（非阻塞）
    start_halo_exchange();
    
    // 2. 计算内部区域（不需要halo）
    dim3 interior_blocks(
        (local_size.x - 2*halo_width.x + 255) / 256,
        local_size.y - 2*halo_width.y,
        local_size.z - 2*halo_width.z
    );
    compute_interior<<<interior_blocks, 256, 0, compute_stream>>>(
        data, local_size, halo_width
    );
    
    // 3. 等待halo交换完成
    finish_halo_exchange();
    
    // 4. 计算边界区域
    compute_boundary_x<<<...>>>();
    compute_boundary_y<<<...>>>();
    compute_boundary_z<<<...>>>();
    
    // 5. 同步所有计算
    cudaDeviceSynchronize();
}
```

### 9.8.3 异步传输与计算重叠

**CUDA流管理**：
```cpp
class MultiGPUManager {
    struct GPUContext {
        int device_id;
        cudaStream_t compute_stream;
        cudaStream_t transfer_stream[2];  // 双缓冲传输
        
        float* d_data;      // 设备内存
        float* h_buffer[2]; // 钉扎内存缓冲区
    };
    
    std::vector<GPUContext> contexts;
    
    void process_multi_gpu() {
        // 为每个GPU创建上下文
        for (int i = 0; i < num_gpus; i++) {
            cudaSetDevice(i);
            contexts[i].compute_stream = create_stream();
            contexts[i].transfer_stream[0] = create_stream();
            contexts[i].transfer_stream[1] = create_stream();
            
            // 分配钉扎内存加速传输
            cudaHostAlloc(&contexts[i].h_buffer[0], buffer_size, 
                         cudaHostAllocPortable);
            cudaHostAlloc(&contexts[i].h_buffer[1], buffer_size, 
                         cudaHostAllocPortable);
        }
        
        // Pipeline处理
        for (int iter = 0; iter < num_iterations; iter++) {
            for (int gpu = 0; gpu < num_gpus; gpu++) {
                cudaSetDevice(gpu);
                auto& ctx = contexts[gpu];
                
                int buf_idx = iter % 2;
                
                // 异步传输下一批数据
                if (iter + 1 < num_iterations) {
                    prepare_next_data(ctx.h_buffer[1-buf_idx], iter+1);
                    cudaMemcpyAsync(ctx.d_data, ctx.h_buffer[1-buf_idx],
                                   data_size, cudaMemcpyHostToDevice,
                                   ctx.transfer_stream[1-buf_idx]);
                }
                
                // 在当前数据上计算
                compute_kernel<<<grid, block, 0, ctx.compute_stream>>>
                              (ctx.d_data, params);
                
                // 等待上一次传输完成
                cudaStreamSynchronize(ctx.transfer_stream[buf_idx]);
            }
        }
    }
};
```

**事件同步**：
```cpp
class EventSync {
    cudaEvent_t events[MAX_GPUS][MAX_EVENTS];
    
    void cross_gpu_dependency() {
        // GPU 0完成计算
        cudaSetDevice(0);
        kernel_gpu0<<<...>>>();
        cudaEventRecord(events[0][0]);
        
        // GPU 1等待GPU 0
        cudaSetDevice(1);
        cudaStreamWaitEvent(stream_gpu1, events[0][0]);
        kernel_gpu1<<<..., stream_gpu1>>>();
        
        // 可以形成复杂的依赖DAG
    }
};
```

### 9.8.4 负载动态迁移

**工作窃取实现**：
```cpp
class WorkStealingScheduler {
    struct WorkQueue {
        std::atomic<int> head;
        std::atomic<int> tail;
        Task* tasks[QUEUE_SIZE];
        
        bool try_push(Task* task) {
            int h = head.load(std::memory_order_acquire);
            int t = tail.load(std::memory_order_relaxed);
            
            if (t - h >= QUEUE_SIZE) return false;
            
            tasks[t % QUEUE_SIZE] = task;
            tail.store(t + 1, std::memory_order_release);
            return true;
        }
        
        Task* try_pop() {
            int h = head.load(std::memory_order_relaxed);
            int t = tail.load(std::memory_order_acquire);
            
            if (h >= t) return nullptr;
            
            Task* task = tasks[h % QUEUE_SIZE];
            head.store(h + 1, std::memory_order_release);
            return task;
        }
        
        Task* try_steal() {
            int h = head.load(std::memory_order_acquire);
            int t = tail.load(std::memory_order_acquire);
            
            if (h >= t) return nullptr;
            
            Task* task = tasks[h % QUEUE_SIZE];
            if (head.compare_exchange_strong(h, h+1)) {
                return task;
            }
            return nullptr;
        }
    };
    
    WorkQueue queues[MAX_GPUS];
    
    void gpu_worker(int gpu_id) {
        cudaSetDevice(gpu_id);
        WorkQueue& my_queue = queues[gpu_id];
        
        while (!done) {
            Task* task = my_queue.try_pop();
            
            if (!task) {
                // 尝试从其他GPU窃取
                for (int victim = 0; victim < num_gpus; victim++) {
                    if (victim != gpu_id) {
                        task = queues[victim].try_steal();
                        if (task) break;
                    }
                }
            }
            
            if (task) {
                execute_task(task, gpu_id);
            } else {
                std::this_thread::yield();
            }
        }
    }
};
```

**动态负载监控**：
```cpp
class LoadMonitor {
    struct GPULoad {
        std::atomic<float> utilization;
        std::atomic<int> task_count;
        std::atomic<int64_t> total_time;
    };
    
    GPULoad loads[MAX_GPUS];
    
    void update_load(int gpu_id, float util, int64_t task_time) {
        loads[gpu_id].utilization.store(util);
        loads[gpu_id].task_count.fetch_add(1);
        loads[gpu_id].total_time.fetch_add(task_time);
    }
    
    int find_least_loaded_gpu() {
        float min_load = 1.0f;
        int best_gpu = 0;
        
        for (int i = 0; i < num_gpus; i++) {
            float load = loads[i].utilization.load();
            if (load < min_load) {
                min_load = load;
                best_gpu = i;
            }
        }
        
        return best_gpu;
    }
    
    void rebalance_work() {
        // 计算平均负载
        float avg_load = 0.0f;
        for (int i = 0; i < num_gpus; i++) {
            avg_load += loads[i].utilization.load();
        }
        avg_load /= num_gpus;
        
        // 识别过载和空闲GPU
        std::vector<int> overloaded, underloaded;
        for (int i = 0; i < num_gpus; i++) {
            float load = loads[i].utilization.load();
            if (load > avg_load * 1.2f) {
                overloaded.push_back(i);
            } else if (load < avg_load * 0.8f) {
                underloaded.push_back(i);
            }
        }
        
        // 迁移任务
        for (int src : overloaded) {
            for (int dst : underloaded) {
                migrate_tasks(src, dst);
            }
        }
    }
};
```

**NCCL集合通信**：
```cpp
void nccl_allreduce_example() {
    ncclComm_t comms[num_gpus];
    
    // 初始化NCCL
    ncclCommInitAll(comms, num_gpus, device_ids);
    
    // 每个GPU执行allreduce
    #pragma omp parallel num_threads(num_gpus)
    {
        int gpu = omp_get_thread_num();
        cudaSetDevice(gpu);
        
        // 局部计算
        compute_local<<<...>>>(d_data[gpu]);
        
        // 全局归约
        ncclAllReduce(d_data[gpu], d_data[gpu], data_size,
                     ncclFloat, ncclSum, comms[gpu], 
                     streams[gpu]);
        
        // 等待完成
        cudaStreamSynchronize(streams[gpu]);
    }
    
    // 清理
    for (int i = 0; i < num_gpus; i++) {
        ncclCommDestroy(comms[i]);
    }
}
```

## 本章小结

本章深入探讨了高性能计算在物理引擎中的应用，从底层硬件架构到高层算法优化：

1. **处理器微架构**：理解CPU流水线、SIMD指令集、分支预测和乱序执行对于编写高效代码至关重要。

2. **内存层级**：缓存优化、NUMA感知编程和TLB优化能够显著提升内存密集型应用的性能。

3. **性能分析**：Roofline模型帮助识别性能瓶颈，性能计数器和分析工具指导优化方向。

4. **并行编程模型**：OpenMP简化共享内存并行，MPI支持分布式计算，合理的负载均衡策略确保高效率。

5. **GPU编程**：理解GPU架构特点，优化内存访问模式和占用率，充分利用GPU的并行计算能力。

6. **稀疏数据结构**：SPGrid和OpenVDB等技术处理大规模稀疏数据，动态拓扑更新支持自适应仿真。

7. **算法优化**：循环优化、向量化、通信隐藏等技术在各个层面提升性能。

8. **多GPU并行**：域分解、halo交换、异步传输和动态负载均衡实现可扩展的多GPU计算。

## 练习题

### 基础题

1. **缓存行分析**
   给定以下代码，分析其缓存行为并优化：
   ```cpp
   struct Particle {
       float x, y, z;
       float vx, vy, vz;
       float mass;
       int type;
   };
   
   void update_positions(Particle* p, int n, float dt) {
       for (int i = 0; i < n; i++) {
           p[i].x += p[i].vx * dt;
           p[i].y += p[i].vy * dt;
           p[i].z += p[i].vz * dt;
       }
   }
   ```
   
   <details>
   <summary>提示</summary>
   考虑AoS vs SoA布局，以及缓存行的利用率。
   </details>
   
   <details>
   <summary>答案</summary>
   
   原始代码问题：
   - Particle结构体大小32字节，一个缓存行(64字节)只能容纳2个粒子
   - 更新位置时只使用了6/8的数据，缓存利用率75%
   
   优化方案1 - SoA布局：
   ```cpp
   struct ParticlesSoA {
       float* x;  float* y;  float* z;
       float* vx; float* vy; float* vz;
       float* mass;
       int* type;
   };
   
   void update_positions_soa(ParticlesSoA& p, int n, float dt) {
       for (int i = 0; i < n; i++) {
           p.x[i] += p.vx[i] * dt;
           p.y[i] += p.vy[i] * dt;
           p.z[i] += p.vz[i] * dt;
       }
   }
   ```
   
   优化方案2 - 数据分离：
   ```cpp
   struct Position { float x, y, z; };
   struct Velocity { float vx, vy, vz; };
   
   void update_positions_split(Position* pos, Velocity* vel, int n, float dt) {
       for (int i = 0; i < n; i++) {
           pos[i].x += vel[i].vx * dt;
           pos[i].y += vel[i].vy * dt;
           pos[i].z += vel[i].vz * dt;
       }
   }
   ```
   </details>

2. **Roofline模型应用**
   某矩阵乘法核心循环执行2n³次浮点运算，读取2n²个浮点数。机器峰值性能1 TFLOPS，内存带宽100 GB/s。问：
   - a) 计算算术强度
   - b) n至少为多少时才能达到计算受限？
   - c) 如何通过分块提高性能？
   
   <details>
   <summary>提示</summary>
   算术强度 = FLOPs / Bytes，机器平衡点 = Peak FLOPS / Bandwidth
   </details>
   
   <details>
   <summary>答案</summary>
   
   a) 算术强度计算：
   - FLOPs = 2n³
   - Bytes = 2n² × 4 = 8n²（假设float）
   - AI = 2n³ / 8n² = n/4
   
   b) 机器平衡点：
   - Balance = 1000 GFLOPS / 100 GB/s = 10 FLOP/Byte
   - 需要 n/4 ≥ 10，即 n ≥ 40
   
   c) 分块优化：
   ```cpp
   // 原始：AI = n/4
   for(i) for(j) for(k)
       C[i][j] += A[i][k] * B[k][j];
   
   // 分块后：AI = TILE/4（每个块内）
   for(ii) for(jj) for(kk)  // 块级循环
       for(i) for(j) for(k)  // 块内循环
           C[ii+i][jj+j] += A[ii+i][kk+k] * B[kk+k][jj+j];
   ```
   选择TILE=64可以充分利用L1缓存，提高数据重用。
   </details>

3. **GPU内存合并访问**
   分析以下GPU kernel的内存访问模式：
   ```cuda
   __global__ void transpose(float* out, float* in, int n) {
       int x = blockIdx.x * blockDim.x + threadIdx.x;
       int y = blockIdx.y * blockDim.y + threadIdx.y;
       
       if (x < n && y < n) {
           out[y * n + x] = in[x * n + y];
       }
   }
   ```
   
   <details>
   <summary>提示</summary>
   考虑warp内线程的访问模式，读和写的合并情况。
   </details>
   
   <details>
   <summary>答案</summary>
   
   问题分析：
   - 读取：相邻线程读取跨步为n的元素，非合并访问
   - 写入：相邻线程写入连续地址，合并访问
   
   优化方案 - 使用共享内存：
   ```cuda
   __global__ void transpose_optimized(float* out, float* in, int n) {
       __shared__ float tile[32][33];  // +1避免bank冲突
       
       int x = blockIdx.x * 32 + threadIdx.x;
       int y = blockIdx.y * 32 + threadIdx.y;
       
       // 合并读取到共享内存
       if (x < n && y < n) {
           tile[threadIdx.y][threadIdx.x] = in[y * n + x];
       }
       __syncthreads();
       
       // 转置后的坐标
       x = blockIdx.y * 32 + threadIdx.x;
       y = blockIdx.x * 32 + threadIdx.y;
       
       // 合并写入
       if (x < n && y < n) {
           out[y * n + x] = tile[threadIdx.x][threadIdx.y];
       }
   }
   ```
   </details>

4. **OpenMP负载均衡**
   有一个粒子系统，每个粒子的计算成本与其邻居数成正比。如何用OpenMP实现动态负载均衡？
   
   <details>
   <summary>提示</summary>
   比较static、dynamic、guided调度策略的特点。
   </details>
   
   <details>
   <summary>答案</summary>
   
   ```cpp
   // 方案1：动态调度
   #pragma omp parallel for schedule(dynamic, 16)
   for (int i = 0; i < n_particles; i++) {
       int neighbors = count_neighbors(i);
       // 计算成本 ∝ neighbors
       compute_forces(i, neighbors);
   }
   
   // 方案2：预排序+静态调度
   // 先统计每个粒子的邻居数
   std::vector<std::pair<int, int>> workload;
   for (int i = 0; i < n_particles; i++) {
       workload.push_back({count_neighbors(i), i});
   }
   
   // 按工作量排序，交错分配
   std::sort(workload.begin(), workload.end());
   
   #pragma omp parallel
   {
       int tid = omp_get_thread_num();
       int nthreads = omp_get_num_threads();
       
       // 循环分配确保负载均衡
       for (int i = tid; i < n_particles; i += nthreads) {
           int particle_id = workload[i].second;
           compute_forces(particle_id, workload[i].first);
       }
   }
   
   // 方案3：任务并行
   #pragma omp parallel
   #pragma omp single
   {
       for (int i = 0; i < n_particles; i++) {
           #pragma omp task
           {
               int neighbors = count_neighbors(i);
               compute_forces(i, neighbors);
           }
       }
   }
   ```
   </details>

### 挑战题

5. **多GPU域分解优化**
   设计一个3D流体仿真的多GPU域分解方案，考虑：
   - 最小化通信量
   - 支持非均匀网格
   - 动态负载均衡
   - 容错性
   
   <details>
   <summary>提示</summary>
   考虑空间填充曲线、递归二分、图分割等方法。
   </details>
   
   <details>
   <summary>答案</summary>
   
   综合方案设计：
   
   1. **初始分解** - 使用递归坐标二分(RCB)：
   ```cpp
   struct Domain {
       float3 min, max;
       int workload;
       int gpu_id;
   };
   
   void recursive_bisection(Domain& domain, int level, int gpu_start, int gpu_count) {
       if (gpu_count == 1) {
           domain.gpu_id = gpu_start;
           return;
       }
       
       // 选择最长维度切分
       int dim = longest_dimension(domain);
       
       // 按工作量平衡找切分点
       float split = find_median_workload(domain, dim);
       
       // 递归分解子域
       Domain left, right;
       split_domain(domain, dim, split, left, right);
       
       int left_gpus = gpu_count / 2;
       recursive_bisection(left, level+1, gpu_start, left_gpus);
       recursive_bisection(right, level+1, gpu_start+left_gpus, gpu_count-left_gpus);
   }
   ```
   
   2. **动态负载均衡** - 基于扩散的方法：
   ```cpp
   class LoadBalancer {
       struct GPUMetrics {
           float computation_time;
           float communication_time;
           int particle_count;
       };
       
       void rebalance() {
           // 收集所有GPU的负载信息
           GPUMetrics metrics[num_gpus];
           gather_metrics(metrics);
           
           // 计算平均负载
           float avg_time = compute_average_time(metrics);
           
           // 扩散算法：邻居间迁移粒子
           for (int iter = 0; iter < max_iterations; iter++) {
               bool converged = true;
               
               for (int gpu = 0; gpu < num_gpus; gpu++) {
                   if (metrics[gpu].computation_time > avg_time * 1.1) {
                       // 找到最空闲的邻居
                       int target = find_least_loaded_neighbor(gpu);
                       
                       // 计算迁移量
                       int migrate_count = estimate_migration(gpu, target, metrics);
                       
                       if (migrate_count > 0) {
                           migrate_particles(gpu, target, migrate_count);
                           converged = false;
                       }
                   }
               }
               
               if (converged) break;
           }
       }
   };
   ```
   
   3. **非均匀网格支持** - 使用八叉树：
   ```cpp
   class AdaptiveGrid {
       struct OctreeNode {
           float3 center, half_size;
           OctreeNode* children[8];
           int particle_count;
           int gpu_owner;
           
           bool should_refine() {
               return particle_count > threshold && !is_leaf();
           }
           
           void assign_to_gpu(int* gpu_loads) {
               if (is_leaf()) {
                   // 分配给负载最小的GPU
                   gpu_owner = min_element(gpu_loads, gpu_loads + num_gpus) - gpu_loads;
                   gpu_loads[gpu_owner] += particle_count;
               } else {
                   for (auto child : children) {
                       child->assign_to_gpu(gpu_loads);
                   }
               }
           }
       };
   };
   ```
   
   4. **容错性** - 检查点和冗余：
   ```cpp
   class FaultTolerantDomain {
       void checkpoint() {
           // 每个GPU保存邻居的ghost数据
           for (int neighbor : neighbors) {
               backup_ghost_data[neighbor] = get_ghost_data(neighbor);
           }
           
           // 异步写入持久存储
           async_write_checkpoint(local_data, backup_ghost_data);
       }
       
       void recover_from_failure(int failed_gpu) {
           // 邻居GPU接管失败GPU的区域
           std::vector<int> neighbors = get_neighbors(failed_gpu);
           
           for (int i = 0; i < neighbors.size(); i++) {
               int gpu = neighbors[i];
               // 从备份恢复ghost区域的数据
               restore_from_backup(gpu, failed_gpu, i);
               
               // 扩展计算域
               extend_domain(gpu, failed_gpu, 1.0 / neighbors.size());
           }
           
           // 重新平衡负载
           rebalance_after_failure();
       }
   };
   ```
   </details>

6. **混合精度优化**
   设计一个自适应混合精度物理仿真框架，在保证精度的同时最大化性能。
   
   <details>
   <summary>提示</summary>
   考虑误差传播、数值稳定性、硬件支持（Tensor Cores）。
   </details>
   
   <details>
   <summary>答案</summary>
   
   ```cpp
   template<typename HighPrec = double, typename LowPrec = float>
   class MixedPrecisionSolver {
       // 误差估计器
       class ErrorEstimator {
           HighPrec reference_solution;
           LowPrec approximate_solution;
           
           float estimate_error() {
               return norm(reference_solution - HighPrec(approximate_solution)) 
                      / norm(reference_solution);
           }
           
           bool needs_refinement(float threshold = 1e-4) {
               return estimate_error() > threshold;
           }
       };
       
       // 自适应精度选择
       class PrecisionSelector {
           struct RegionInfo {
               float3 min, max;
               float error_estimate;
               bool use_high_precision;
           };
           
           std::vector<RegionInfo> regions;
           
           void analyze_regions() {
               for (auto& region : regions) {
                   // 基于物理量梯度判断
                   float gradient = compute_gradient(region);
                   
                   // 基于条件数判断
                   float condition = estimate_condition_number(region);
                   
                   // 决策
                   region.use_high_precision = 
                       gradient > gradient_threshold ||
                       condition > condition_threshold ||
                       region.error_estimate > error_threshold;
               }
           }
       };
       
       // 混合精度矩阵运算（利用Tensor Cores）
       void mixed_precision_gemm(float* C, const half* A, const half* B, 
                                int M, int N, int K) {
           // 使用Tensor Cores进行半精度计算
           cublasGemmEx(handle,
               CUBLAS_OP_N, CUBLAS_OP_N,
               M, N, K,
               &alpha,
               A, CUDA_R_16F, M,
               B, CUDA_R_16F, K,
               &beta,
               C, CUDA_R_32F, M,
               CUDA_R_32F,
               CUBLAS_GEMM_DEFAULT_TENSOR_OP);
       }
       
       // 迭代细化求解器
       void iterative_refinement_solve(Matrix<HighPrec>& A, 
                                     Vector<HighPrec>& b,
                                     Vector<HighPrec>& x) {
           // 转换为低精度
           Matrix<LowPrec> A_low = convert<LowPrec>(A);
           Vector<LowPrec> b_low = convert<LowPrec>(b);
           Vector<LowPrec> x_low = convert<LowPrec>(x);
           
           // 低精度LU分解
           LU<LowPrec> lu(A_low);
           
           for (int iter = 0; iter < max_iterations; iter++) {
               // 计算残差（高精度）
               Vector<HighPrec> r = b - A * x;
               
               // 检查收敛
               if (norm(r) < tolerance) break;
               
               // 低精度求解修正量
               Vector<LowPrec> r_low = convert<LowPrec>(r);
               Vector<LowPrec> dx_low = lu.solve(r_low);
               
               // 高精度更新
               x += convert<HighPrec>(dx_low);
           }
       }
       
       // 自适应时间步长with混合精度
       void adaptive_timestep() {
           float dt_high = dt;
           float dt_low = dt;
           
           // 用两种精度各走一步
           State<HighPrec> state_high = advance<HighPrec>(state, dt_high);
           State<LowPrec> state_low = advance<LowPrec>(state, dt_low);
           
           // 估计误差
           float error = estimate_error(state_high, state_low);
           
           if (error > tolerance) {
               // 需要更高精度或更小步长
               dt *= 0.5;
               use_high_precision = true;
           } else if (error < 0.1 * tolerance) {
               // 可以用更低精度或更大步长
               dt *= 1.5;
               use_high_precision = false;
           }
       }
   };
   ```
   </details>

7. **SIMD优化的SPH实现**
   实现一个高度优化的SPH邻居搜索和力计算，充分利用AVX-512指令集。
   
   <details>
   <summary>提示</summary>
   考虑数据布局、掩码操作、向量化的距离计算。
   </details>
   
   <details>
   <summary>答案</summary>
   
   ```cpp
   class SimdSPH {
       // 16路SIMD的粒子数据（AVX-512）
       struct ParticleBlock {
           __m512 x[BLOCK_SIZE/16];
           __m512 y[BLOCK_SIZE/16];
           __m512 z[BLOCK_SIZE/16];
           __m512 vx[BLOCK_SIZE/16];
           __m512 vy[BLOCK_SIZE/16];
           __m512 vz[BLOCK_SIZE/16];
           __m512 density[BLOCK_SIZE/16];
           __m512 pressure[BLOCK_SIZE/16];
       };
       
       // 向量化的邻居搜索
       void find_neighbors_simd(int particle_id, std::vector<int>& neighbors) {
           __m512 px = _mm512_set1_ps(particles.x[particle_id]);
           __m512 py = _mm512_set1_ps(particles.y[particle_id]);
           __m512 pz = _mm512_set1_ps(particles.z[particle_id]);
           __m512 h2 = _mm512_set1_ps(h * h);
           
           // 获取粒子所在的网格单元
           int3 cell = get_cell(particle_id);
           
           // 遍历27个邻居单元
           for (int dz = -1; dz <= 1; dz++) {
               for (int dy = -1; dy <= 1; dy++) {
                   for (int dx = -1; dx <= 1; dx++) {
                       int3 ncell = cell + make_int3(dx, dy, dz);
                       int start = cell_start[ncell];
                       int end = cell_end[ncell];
                       
                       // SIMD处理16个粒子
                       for (int i = start; i < end; i += 16) {
                           // 加载16个粒子位置
                           __m512 qx = _mm512_load_ps(&particles.x[i]);
                           __m512 qy = _mm512_load_ps(&particles.y[i]);
                           __m512 qz = _mm512_load_ps(&particles.z[i]);
                           
                           // 计算距离平方
                           __m512 dx = _mm512_sub_ps(px, qx);
                           __m512 dy = _mm512_sub_ps(py, qy);
                           __m512 dz = _mm512_sub_ps(pz, qz);
                           
                           __m512 r2 = _mm512_fmadd_ps(dx, dx,
                                       _mm512_fmadd_ps(dy, dy,
                                       _mm512_mul_ps(dz, dz)));
                           
                           // 生成掩码：r2 < h2
                           __mmask16 mask = _mm512_cmp_ps_mask(r2, h2, _CMP_LT_OQ);
                           
                           // 压缩存储邻居索引
                           if (mask != 0) {
                               int indices[16];
                               for (int j = 0; j < 16; j++) {
                                   indices[j] = i + j;
                               }
                               
                               // 使用compress存储有效邻居
                               __m512i idx = _mm512_load_epi32(indices);
                               __m512i compressed = _mm512_maskz_compress_epi32(mask, idx);
                               
                               // 存储到邻居列表
                               int count = _mm_popcnt_u32(mask);
                               _mm512_storeu_epi32(&neighbors[neighbors.size()], compressed);
                               neighbors.resize(neighbors.size() + count);
                           }
                       }
                   }
               }
           }
       }
       
       // 向量化的密度计算
       void compute_density_simd() {
           const __m512 factor = _mm512_set1_ps(mass * poly6_constant);
           
           #pragma omp parallel for
           for (int i = 0; i < particle_count; i += 16) {
               __m512 density_sum = _mm512_setzero_ps();
               
               // 对每个粒子的邻居进行向量化处理
               for (int j = 0; j < 16; j++) {
                   int pid = i + j;
                   if (pid >= particle_count) break;
                   
                   __m512 px = _mm512_set1_ps(particles.x[pid]);
                   __m512 py = _mm512_set1_ps(particles.y[pid]);
                   __m512 pz = _mm512_set1_ps(particles.z[pid]);
                   
                   // 处理邻居（已对齐到16的倍数）
                   int neighbor_count = neighbor_lists[pid].size();
                   for (int k = 0; k < neighbor_count; k += 16) {
                       // 加载16个邻居
                       __m512i indices = _mm512_load_epi32(&neighbor_lists[pid][k]);
                       
                       // Gather邻居位置
                       __m512 qx = _mm512_i32gather_ps(indices, particles.x, 4);
                       __m512 qy = _mm512_i32gather_ps(indices, particles.y, 4);
                       __m512 qz = _mm512_i32gather_ps(indices, particles.z, 4);
                       
                       // 计算核函数值
                       __m512 dx = _mm512_sub_ps(px, qx);
                       __m512 dy = _mm512_sub_ps(py, qy);
                       __m512 dz = _mm512_sub_ps(pz, qz);
                       
                       __m512 r2 = _mm512_fmadd_ps(dx, dx,
                                   _mm512_fmadd_ps(dy, dy,
                                   _mm512_mul_ps(dz, dz)));
                       
                       // W = (h² - r²)³（简化的poly6核）
                       __m512 h2 = _mm512_set1_ps(h * h);
                       __m512 diff = _mm512_sub_ps(h2, r2);
                       __m512 diff2 = _mm512_mul_ps(diff, diff);
                       __m512 w = _mm512_mul_ps(diff2, diff);
                       
                       // 累加密度贡献
                       density_sum = _mm512_fmadd_ps(factor, w, density_sum);
                   }
               }
               
               // 水平归约并存储
               _mm512_store_ps(&particles.density[i], density_sum);
           }
       }
       
       // 向量化的压力计算
       void compute_forces_simd() {
           #pragma omp parallel for schedule(dynamic, 16)
           for (int i = 0; i < particle_count; i++) {
               __m512 fx = _mm512_setzero_ps();
               __m512 fy = _mm512_setzero_ps(); 
               __m512 fz = _mm512_setzero_ps();
               
               __m512 pi = _mm512_set1_ps(particles.pressure[i]);
               __m512 di = _mm512_set1_ps(particles.density[i]);
               
               // 处理所有邻居
               process_neighbors_forces_simd(i, fx, fy, fz, pi, di);
               
               // 归约力到标量
               particles.fx[i] = _mm512_reduce_add_ps(fx);
               particles.fy[i] = _mm512_reduce_add_ps(fy);
               particles.fz[i] = _mm512_reduce_add_ps(fz);
           }
       }
   };
   ```
   </details>

8. **能耗感知的性能优化**
   设计一个能够在性能和能耗之间自动平衡的物理引擎调度器。
   
   <details>
   <summary>提示</summary>
   考虑DVFS、任务调度、精度调整、硬件特性。
   </details>
   
   <details>
   <summary>答案</summary>
   
   ```cpp
   class EnergyAwareScheduler {
       // 能耗模型
       struct PowerModel {
           // 动态功耗：P = C × V² × f
           float capacitance;
           float voltage;
           float frequency;
           
           // 静态功耗
           float static_power;
           
           float compute_power() {
               return capacitance * voltage * voltage * frequency + static_power;
           }
           
           // 能效：Performance per Watt
           float efficiency(float performance) {
               return performance / compute_power();
           }
       };
       
       // DVFS控制器
       class DVFSController {
           struct FrequencyLevel {
               float frequency;
               float voltage;
               float performance_scale;
           };
           
           std::vector<FrequencyLevel> levels;
           int current_level;
           
           void adjust_frequency(float target_fps, float current_fps, float power_budget) {
               // PID控制器
               float error = target_fps - current_fps;
               float power = measure_power();
               
               if (power > power_budget) {
                   // 超过功耗预算，降频
                   decrease_frequency();
               } else if (error > 0 && power < 0.9 * power_budget) {
                   // 性能不足且有功耗余量，升频
                   increase_frequency();
               }
               
               // 考虑温度限制
               float temp = read_temperature();
               if (temp > THERMAL_LIMIT) {
                   current_level = std::min(current_level, get_safe_level(temp));
               }
           }
       };
       
       // 任务能效调度
       class TaskScheduler {
           struct Task {
               std::function<void()> work;
               float expected_time;
               float expected_energy;
               int preferred_core;  // 大核/小核偏好
           };
           
           // 大小核调度
           void schedule_heterogeneous(std::vector<Task>& tasks) {
               // 按能效比排序
               std::sort(tasks.begin(), tasks.end(), 
                   [](const Task& a, const Task& b) {
                       return a.expected_time / a.expected_energy > 
                              b.expected_time / b.expected_energy;
                   });
               
               // 能效高的任务分配到小核
               for (auto& task : tasks) {
                   if (task.expected_energy < ENERGY_THRESHOLD) {
                       task.preferred_core = LITTLE_CORE;
                   } else {
                       task.preferred_core = BIG_CORE;
                   }
               }
               
               // 动态迁移
               monitor_and_migrate(tasks);
           }
           
           void monitor_and_migrate(std::vector<Task>& tasks) {
               for (auto& task : tasks) {
                   float actual_energy = measure_task_energy(task);
                   
                   // 如果实际能耗与预期差异大，迁移到合适的核
                   if (std::abs(actual_energy - task.expected_energy) > THRESHOLD) {
                       migrate_task(task);
                   }
               }
           }
       };
       
       // 自适应精度控制
       class AdaptivePrecision {
           enum PrecisionLevel {
               FP64,     // 最高精度，最高能耗
               FP32,     // 标准精度
               FP16,     // 半精度
               INT8      // 量化，最低能耗
           };
           
           PrecisionLevel select_precision(float error_tolerance, float power_budget) {
               struct Config {
                   PrecisionLevel level;
                   float relative_power;
                   float max_error;
               };
               
               Config configs[] = {
                   {FP64, 1.0,   1e-15},
                   {FP32, 0.5,   1e-7},
                   {FP16, 0.25,  1e-3},
                   {INT8, 0.1,   1e-2}
               };
               
               // 选择满足精度要求的最低能耗配置
               for (auto& config : configs) {
                   if (config.max_error <= error_tolerance &&
                       config.relative_power <= power_budget) {
                       return config.level;
                   }
               }
               
               return FP32;  // 默认
           }
       };
       
       // 能耗感知的负载均衡
       class PowerAwareLoadBalancer {
           void balance_with_power_constraint(float total_power_budget) {
               int num_gpus = get_gpu_count();
               std::vector<float> gpu_powers(num_gpus);
               std::vector<float> gpu_performances(num_gpus);
               
               // 测量每个GPU的功耗和性能
               for (int i = 0; i < num_gpus; i++) {
                   gpu_powers[i] = measure_gpu_power(i);
                   gpu_performances[i] = measure_gpu_performance(i);
               }
               
               // 优化问题：最大化总性能，约束总功耗
               // maximize: Σ performance[i] × active[i]
               // subject to: Σ power[i] × active[i] ≤ power_budget
               
               // 贪心算法：按性能功耗比排序
               std::vector<int> indices(num_gpus);
               std::iota(indices.begin(), indices.end(), 0);
               
               std::sort(indices.begin(), indices.end(),
                   [&](int a, int b) {
                       return gpu_performances[a] / gpu_powers[a] >
                              gpu_performances[b] / gpu_powers[b];
                   });
               
               // 选择GPU直到达到功耗预算
               float used_power = 0;
               std::vector<bool> active(num_gpus, false);
               
               for (int idx : indices) {
                   if (used_power + gpu_powers[idx] <= total_power_budget) {
                       active[idx] = true;
                       used_power += gpu_powers[idx];
                   }
               }
               
               // 重新分配工作负载
               redistribute_work(active);
           }
       };
       
       // 主调度循环
       void energy_aware_schedule() {
           while (running) {
               // 读取系统状态
               SystemState state = read_system_state();
               
               // 预测下一帧的计算需求
               WorkloadPrediction prediction = predict_workload();
               
               // 优化能效
               OptimalConfig config = optimize_energy_efficiency(
                   state, prediction, constraints);
               
               // 应用配置
               apply_configuration(config);
               
               // 执行仿真步
               execute_simulation_step();
               
               // 记录能耗数据用于future优化
               record_energy_metrics();
           }
       }
   };
   ```
   </details>

## 常见陷阱与错误

1. **伪共享**：多线程访问同一缓存行的不同数据导致性能严重下降。使用padding或__cacheline_aligned避免。

2. **内存带宽饱和**：盲目增加线程数不一定提升性能，要考虑内存带宽限制。

3. **GPU占用率陷阱**：100%占用率不等于最佳性能，要平衡占用率和寄存器使用。

4. **分支发散**：GPU中的条件分支会导致warp内线程执行不同路径，严重影响性能。

5. **原子操作竞争**：大量线程竞争同一原子变量会造成串行化，考虑使用归约树。

6. **NUMA远程访问**：未考虑NUMA拓扑的内存分配会导致远程访问，性能下降50%以上。

7. **过度优化**：不要优化非热点代码，先profile再优化。

8. **缓存抖动**：工作集略大于缓存容量时性能断崖式下跌，考虑分块或流式处理。

## 最佳实践检查清单

### 设计阶段
- [ ] 识别计算密集vs内存密集的部分
- [ ] 选择合适的并行粒度
- [ ] 设计缓存友好的数据结构
- [ ] 考虑NUMA和GPU内存层次
- [ ] 预估通信开销

### 实现阶段
- [ ] 数据对齐到缓存行边界
- [ ] 使用SoA布局提高向量化效率
- [ ] 避免false sharing
- [ ] 最小化同步点
- [ ] 合并内存访问（GPU）

### 优化阶段
- [ ] Profile识别真正的瓶颈
- [ ] 测量并优化缓存命中率
- [ ] 检查向量化效率
- [ ] 验证负载均衡
- [ ] 监控功耗和温度

### 调试阶段
- [ ] 使用内存检查工具（valgrind、cuda-memcheck）
- [ ] 检查数据竞争（ThreadSanitizer）
- [ ] 验证数值精度
- [ ] 测试极端情况（空数据、超大规模）
- [ ] 检查内存泄漏

### 部署阶段
- [ ] 针对目标硬件调优
- [ ] 设置合理的默认参数
- [ ] 提供性能调节选项
- [ ] 记录性能特征
- [ ] 监控生产环境表现
