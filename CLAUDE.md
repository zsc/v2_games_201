（交流可以用英文，本文档中文，保留这句）

## 项目目标
编写一份高级物理引擎中文教程markdown，要包含大量的习题和参考答案（答案默认折叠）。能数学算尽量数学算，合适时提及相关函数名但不写代码。
内容范围1:1严格参考 old.md。挑战型练习题可以适度拔高。
内容在 https://forum.taichi-lang.cn/t/topic/272 的各讲讲义的基础上，进行扩展补充相关材料。
文件组织是 index.md + chapter1.md + ...


## Audience
verteran programmer and AI scientists

## 章节结构要求
每个章节应包含：
1. **开篇段落**：简要介绍本章内容和学习目标
2. **本章小结**：总结关键概念和公式
3. **练习题**：
   - 每章包含6-8道练习题
   - 50%基础题（帮助熟悉材料）
   - 50%挑战题（包括开放性思考题）
   - 每题提供提示（Hint）
   - 答案默认折叠，不包含代码
4. **常见陷阱与错误** (Gotchas)：每章包含该主题的常见错误和调试技巧
5. **最佳实践检查清单**：每章末尾提供设计审查要点

## 术语中英对照表

| 英文 | 中文 |
|------|------|
| Lagrangian | 拉格朗日 |
| Eulerian | 欧拉 |
| Mass-spring systems | 弹簧质点系统 |
| Explicit/implicit time integrators | 显式/隐式时间积分器 |
| Smoothed particle hydrodynamics | 光滑粒子流体动力学 |
| Position-based fluids | 基于位置的流体 |
| Voxelization | 体素化 |
| Neighborhood search | 邻居搜索 |
| Weak form | 弱形式 |
| Hexahedron grid | 六面体网格 |
| Tetrahedra mesh | 四面体网格 |
| Invertible | 可逆 |
| Topology optimization | 拓扑优化 |
| Semi-Lagrangian advection | 半拉格朗日输送 |
| Chorin-Style projection | Chorin式压力投影 |
| Staggered grid | Staggered网格 |
| Nullspaces | 零空间 |
| Krylov-subspace solvers | Krylov子空间求解器 |
| Matrix-free | 无矩阵 |
| Poisson equations | 泊松方程 |
| Preconditioning | 预条件 |
| Multigrid methods | 多重网格方法 |
| Advection schemes | 输送格式 |
| Signed-distance fields | 有符号距离场 |
| Level sets | 等势面 |
| Free-surface | 自由表面 |
| Path tracing | 路径追踪 |
| Sphere tracing | 球面追踪 |
| Motion blur | 运动模糊 |
| Marching cubes | 行军立方体 |
| Digital differential analyzer | 数字微分分析器 |
| Voxel rendering | 体素渲染 |
| Volumetric rendering | 体积渲染 |
| Particle-grid transfers | 粒子-网格传输 |
| Particle-in-Cell | 粒子元胞法 |
| Fluid Implicit Particles | 流体隐粒子 |
| Affine Particle-in-Cell | 仿射粒子元胞法 |
| Material Point Method | 物质点法 |
| Constitutive models | 本构模型 |
| Moving Least Squares MPM | 移动最小二乘物质点法 |
| Lagrangian forces | 拉格朗日力 |
| Numeric fracture | 数值断裂 |
| Continuum Damage Mechanics | 连续介质伤害力学 |
| Processor microarchitectures | 处理器微结构 |
| Memory hierarchy | 内存层级 |
| Sparse data structures | 稀疏数据结构 |
| Differentiable programming | 可微编程 |
| Checkpointing | Checkpointing技巧 |

