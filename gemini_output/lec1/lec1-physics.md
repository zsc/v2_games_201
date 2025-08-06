Okay, I have extracted the text from lec1-physics.pdf. Due to the limitations of OCR, there might be some inaccuracies, especially with equations and formatting.

Here is the transcribed text:

```text
GAMES 201
Advanced Physics Engines 2020: A Hands-on Tutorial
高级物理引擎实战2020
(基于太极编程语言)
Yuanming Hu
胡渊鸣
MIT CSAIL
麻省理工学院 计算机科学与人工智能实验室
Inflow speed: 0.5
Inflow speed: 1.0
Wheel density: 10
Wheel density: 40
Wheel
Taichi
Programming Language
高级物理引擎实战2020
♦ 课程目标:自己动手打造影视级物理引擎
◆ 适合人群:0-99岁的计算机图形学爱好者
预备知识:高等数学、Python或任何一门程序设计语言
♦ 课程安排:每周一北京时间晚上20:30-21:30 共10节课
♦ 课程内容:Taichi语言基础 刚体 液体 烟雾弹塑性体 PIC/FLIP法 Krylov-子空
间求解器 预条件 无矩阵法 多重网格弱形式与有限元 隐式积分器 辛积分器
拓扑优化带符号距离场 自由表面追踪物质点法 大规模物理效果渲染 现代处
理器微架构 内存层级并行编程 GPU编程 稀疏数据结构 可微编程...
2
Physics engine
◆ “A physics engine is computer software that provides an
approximate simulation of certain physical systems, such as
rigid body dynamics (including collision detection), soft body
dynamics, and fluid dynamics, of use in the domains of
computer graphics, video games and film."
https://en.wikipedia.org/wiki/Physics engine
3
TL; DR:
Simulate the world in your computer!
4
Applications
CAD/CAE
◆ Visual effects (films)
VR/AR
◆ Training robots
Games!
https://en.wikipedia.org/wiki/Computer-aided engineering
SALOME 2.2.2-[plasticityPost.hdf#1]
Ele Edit Ver Voussion Selection Representation fods Preferences ndow p
Mesh Lin
Fanles
onNodes
Groups
0.768
Son DEPL
onNodes
0.142057, INCONNUE
0205714, INCONNUE
0420571, INCONNUE
-0.571429, INCONNUE
0714206, INCONNUE
0.640
0857143, INCONNUE
1, INCONNUE
Sabitip
Def Shape 7
on Cells
-0.142057, INCONNUE
Son EQUIELNO SIGM
0.512
0205714, INCONNUE
0420571, INCONNUE
ho Surfices
0571429, INCONNUE
-071426, INCONNUE
0857143, INCONNUE
0.384
-1, INCONNUE
ho Surfices 2
-CutPines
Cut Planes:1
CutPlanes 2
Cut Planes
CutPanes 4
0.256
Obanct Browser
Post-Pro
https://gym.openai.com/envs/#box2d
-OX
SALOME
https://www.math.ucla.edu/~jteran/papers/SSCTS13.pdf
5
https://www.youtube.com/watch?v=P_796gF0Cq4&t
The Incredible machines (1993)
ADG
SCORE
7
8
1
ENE
000000
BONUS
BONUS 2
謝
6
2-2-0
7
Angry Birds (2009)
"
11
8
15
3
50
TNT
https://www.theguardian.com/artanddesign/2016/feb/23/how-we-made-angry-birds
7
Phun (Algodoo) 2009
File
?
Algodoo for Education v2.0.1
2
..
土
T
mg
N
mg
8
mg
P
Plot - Spring
X
3 m/s
spring
لاع
rrig
http://www.algodoo.com/what-is-it/#
1 m
ㅓ
9
Besieged (2015)
: The Legend of Zelda: Breath of the Wild (2017)
https://www.youtube.com/watch?v=LamMQ47ccdc
12:30 PM
Ori and the Will of the Wisps (2020)
2009-2020
初中二年级~博士三年级
My Own Physical Simulation Story
12
My first physics engine (2009)
Applications Places System
1.672
graphics.h [Game] Code::Blocks 10.05
Eile Edit View Search Project Build Debug Tools Plugins Settings Help
Screen
Dra
Management
main.cpp game.h graphics.h
Projects Symbols
28
*((Uint32
▼Workspace
29
}
30
Game
31
►Sources
32
void DrawLine(i
Headers
33
for (double
base.h
34
}
big.h
35
36
inline void Dra
game.h
37
if (x < 0)
graphics.
38
if (y < 0)
input.h
39
if (x+w>
object.h
40
if (y+h>
41
for (int i
42
for (in
43
Dra
44
}
45
46
inline void Dra
47
//
if (x < 6
48
//
if (y < 6
Logs & others
Code::Blocks Search results Be
Checking for existence: /home/huyuanming/Documen
Executing: xterm -T Ganee /usr/bin/cb_console
Documents/Gane/WOG/.)
/home/huyuanming/Documents/Game/WOG/graphics.h
WebQQ... graphic... Music P...
Sat Jan 15, 8:49 PM
huyuanming
UTF-8
Line 4, Column 17
Insert
Read/Write default
Game
huyuan...
data
Orz...
2008 π...
Game
Untitle...
13
My first physics engine (2009)
Applications Places System
22
Sun Jan 16, 6:28 PM
huyuanming
WebQQ2.0....
Untitled.gif
huyuanming@... Gnuplot (wind...
base.h (Game]
Game
Untitled window
14
My own rigid body simulator (~2011)
Settings
Iterations
Velocity Iterations 39
Position Iterations 5
Step Iterations
30
World Properties
Friction
Pause
Gravity
9.8
Information
Frames Per Second 20
Objects
19
Shapes
259
Forces
0
(Middle school hobby project)
◆ Still works today after 9 years! (Demo)
https://github.com/yuanming-hu/FLAT
15
Settings
-Iterations
Velocity Iterations
40
Position Iterations
5
Step Iterations
30
World Properties
Friction
ي
Pause
Gravity
9.8
-Information
Frames Per Second 60
Objects
65
Shapes
69
Forces
126
Constraints
D
Settings
-Iterations
Velocity Iterations
40
Position Iterations
5
Step Iterations
30
-World Properties
Friction
Pause
Gravity
0.8
-Information
Frames Per Second 60
Objects
33
Shapes
33
Forces
D
Constraints
0
My rigid body game (2014)
http://yuanming-hu.github.io/gear/
16
My fluid simulator (2016)
GPU-based Fluid Simulation
Visualization
Control
Parameters
Warm Starting
II PAUSE
RESET
RK2 Advection
Jacobi Iterations
Simulation Method
PIC/FLIP
Simulation Resolution
128x128
Initial State
Dam Break (Left)
Settings above will be
applied after RESETTING.
Visual Particle Size
10
Jacobi Damping (Param. for the Damped
Jacobi pressure solver)
Frame Time Step
0.67
0.030
2.5
Substeps (subdivisions of frame timestep)
FLIP Blending (smaller = more viscous)
http://yuanming-hu.github.io/fluid/
17
10
0.80
My smoke simulator (2017)
18
My continuum + cutting simulator
(MLS-MPM+CPIC) SIGGRAPH 2018
米
Inflow speed: 0.5
Inflow speed: 1.0
Wheel density: 1.0
Wheel density: 4.0
Inflow speed: 0.5
Wheel density: 4.0
000
19
+
+
+
https://github.com/yuanming-hu/taichi_mpm
Solving FEM linear elasticity on 1,040,875,347 voxels (SIGGRAPH Asia 2018)
20
https://github.com/yuanming-hu/spgrid_topo_opt
Differentiable Simulation (ICRA 2018, ICLR 2019)
Center Activation
Iteration 60
Iteration 180
Height field fluid simulation.
Optimize the initial height field so that it forms "Taichi" after 256 time steps.
21
Now (2020)
♦ MIT EECS PhD student (3rd year)
• I still work on simulations
...
...
...
but spend most of the time crafting a compiler
so that people can more easily write simulators
The Taichi Programming Language (BDFL)
GAMES 201 (Lecturer)
22
Keywords of this course
23
Discretizing the force. To reach the discrete force, Eq. 12 requires
the derivative of q. Differentiating Eq. 13 gives
9α,β(x, t) =
ӘРТ (x – x)
дхв
n
X
n
PM¯¹(x)(x)P(x – x)δαά· (14)
n
A Moving Least Squares Material Point Method with Displacement Discontinuity
and Two-Way Rigid Body Coupling, Hu, Fang, Ge, Qu, Zhu, Pradhana, Jiang
PHI
Multigrid V-Cycle: Solving PHI in PDE f(PHI) = F
F
R
R1
R2
R
Gauss Seidel
Compute Residuals
Repeat Until Convergence
Restrict
Gauss Seidel
← Compute Residuals
Restrict
Gauss Seidel
← Compute Residuals
Restrict
Gauss Seidel
Set R = 0
Interpolate
Gauss Seidel
Correct
Interpolate
Gauss Seidel
Correct
Interpolate
Correct
PHI
Efficient Solvers
https://en.wikipedia.org/wiki/Multigrid_method
Productivity
68
for p in x: # grid to particle (G2P)
69
base = (x[p] * inv_dx - 0.5).cast(int)
70
fx = x[p] * inv_dx - base.cast(float)
71
W =
[0.5 * (1.5 - fx) ** 2, 0.75
-
(fx - 1.0) ** 2, 0.5 * (fx - 0.5) ** 2]
72
new_v = ti.Vector.zero(ti.f32, 2)
73
new_C = ti.Matrix.zero(ti.f32, 2, 2)
74
for i, jin ti.static(ti.ndrange(3, 3)): # loop over 3x3 grid node neighborhood
75
dpos = ti. Vector((i, j)).cast(float)
-
fx
76
g_v = grid_v [base + ti. Vector([i, j])]
77
weight = w[i][0] * w[j][1]
78
new_v += weight * g_v
79
80
new_C += 4 * inv_dx * weight * ti.outer_product(g_v, dpos)
v[p], C[p] = new_v, new_C
81
x[p] += dt * v [p] # advection
https://github.com/taichi-dev/taichi/blob/master/examples/mpm99.py
524
// Loop start
525
#ifdef MLSMPM
Performance
526 #define LOOP(node_id)
527
{
\
\
528
_m128 dpos = _mm_sub_ps(rela_pos, grid_pos_offset_[node_id]);
\
529
m128 g
=
530
grid_cache
\
531
.linear [grid_cache.kernel_linearized (node_id) + grid_cache_offset] \
532
.v;
\
533
_m128 weight =
534
_mm_set1_ps(kernels [node_id / 9] [node_id / 3 % 3] [node_id % 3]);
535
_ m128 affine_prod = _mm_fmadd_ps(
\
536
affine [2], broadcast(dpos, 2),
\
537
_mm_fmadd_ps(affine [1], broadcast(dpos, 1),
538
_mm_fmadd_ps(affine [0], broadcast(dpos, 0), mass_v)));
\
539
540
m128 contrib = _mm_blend_ps (mass_, affine_prod, 0x7);
_m128 delta = _mm_mul_ps (weight, contrib);
541
g = _mm_add_ps(g, delta);
542
grid_cache
543
.linear[grid_cache.kernel_linearized (node_id) + grid_cache_offset]
544
.v = g;
545
}
\
\
\
https://github.com/yuanming-hu/taichi_mpm/blob/3bb90fbe4c901aafc048dbb2d8d8aa388226d011/src/transfer.cpp#L524-L545
Main Memory
35.8 GB/s
256 cyc latency
L3 cache 2M/core
134.4 GB/s
42 cyc latency
L2 cache 256KB
268.8 GB/s
12 cyc latency
L1 data cache 32KB
403.2 GB/s
4 cyc latency
L2 Unified TLB (STLB)
Hardware architecture
(Part of) the Memory Hierarchy
4 KB/2MB pages - 1536 entries
1G pages - 16 entries
L1 Data TLB
4 KB pages - 64 entries
2/4 MB pages - 32 entries
1G pages - 4 entries
*
*
*
Figures are not drawn to scale.
Instruction caches are omitted.
Main memory BW is shared by all
cores.
CPU core
closer to CPU,
smaller capacity,
lower latency,
higher bandwidth.
Integer Physical Registers
8 bytes per entry, 180
entries
1 cyc latency
Vector Physical Registers
32 byte entries, 168
Execution Engine
4.20 GHz
entries
1 cyc latency
134.4 G FLOP/s, i.e. 806.4 GB/s bandwidth requirement
// "Hybrid Eulerian-
Lagrangian Grid"
auto &block = root
.hash(ij, 8)
.dense(ij, {16,8})
.bitmasked();
// Child 1: grid nodes
ock.dense(ij, {4, 8})
.place(grid_vx)
.place(grid_vy);
// Child 2: particles
block.dynamic(i, 32)
.place (part_x)
.place (part_y)
.place(part_vx)
.place (part_vy);
Data Structures
hashed
({i,j},{8,8})
dense({i,j}, {16,8})
dense
.bitmasked()
dynamic
({i,j},{4,8}) (k, 32)
place
float32
i=[256,512)
|j=[256,512)
particle
position x
float32
i=[256,512)
particle
position y
j=[0,256)
float32
particle
i=[0,256)
velocity x
|j=[0,256)
float32
particle
velocity y
...
place float32 grid velocity y
float32 grid velocity y
(c) "HPB", "SPVDB", “HLEG”: We can easily design new data structures with customized features.
Taichi: A Language for High-Performance Computation on Spatially Sparse Data Structures SIGGRAPH Asia 2019
Yuanming Hu, Tzu-Mao Li, Luke Anderson, Jonathan Ragan-Kelley, Fredo Durand
Single-GPU
update grid
halo gather & send
halo
Naïve
Multi-GPU G2P2G
MGSP
{
update grid
halo
G2P2G
G2P2G
update
partitions
receive & reduce halo
non-halo
G2P2G
update
partitions
GPU 1 Block
GPU 2 Block
Halo Regions
non-halo
G2P2G
update
halo
partitions tagging
gather & send
halo
receive & reduce halo
Fig. 9. Instruction pipelines. The additional operations in the multi-GPU
MPM compared to the single-GPU MPM are displayed in red. By masking
these halo-region-related data transfers with the execution of the G2P2G
kernel, one can achieve more optimized scaling results with multi-GPUs.
A Massively Parallel and Scalable Multi-GPU Material Point Method,
Advected
Transferred and Advected
Assigned, Transferred,
and Advected
Xinlei Wang*, Yuxing Qiu* (equal contributions), Stuart Slattery, Yu Fang, Minchen Li, Song-Chun Zhu, Yixin Zhu, Min Tang, Dinesh Manocha, Chenfanfu Jiang (SIGGRAPH 2020)
ӘХ
diffmpm.py
liquid.py
rigid_body.py mass_spring.py
water surface
water_renderer.py volume_renderer.py smoke.py
wave.py
electric.py
billiards.py
Differentiability
DiffTaichi: Differentiable Programming for Physical Simulation, ICLR 2020
Yuanming Hu, Luke Anderson, Tzu-Mao Li, Qi Sun, Nathan Carr, Jonathan Ragan-Kelley, Fredo Durand
Introduction to
☑ Taichi
Programming Language
32
Homework
◆五次简单的编码练习(Coding exercises)
●提供标准答案,但是建议看答案之前先自行尝试
·无需提交
◆ 三个开放项目 (Open projects)
●开放项目可以1-3人组队(玩得开心!)
●评奖
▶ 获奖作品会得到课上点评
▶ 获奖同学会得到神秘礼物
◆作业自愿完成
33
「纸上得来终觉浅,绝知此事要躬行。」
———《冬夜读书示子聿》陆游 1199 CE
34
Homework 0
1.Install Taichi (with Python 3.6+: pip install taichi)
2.Have fun with the Taichi examples!
• python3 -m taichi example mpm88
3.Write your own Taichi programs
• Practice makes perfect
• Don't forget to share them with your classmates on the Taichi forum!
• We will pick outstanding programs on the Taichi forum and
demonstrate the results to the class
• Finish early if you would like to be picked (we collect 48h before the
next lecture)
35
Teaching Assistants
◆理论组负责在论坛答疑、筛选优秀开放项目:
• 夏一鸣 ymxia@iphy.ac.cn
• 史雨宸 syc0412@mail.ustc.edu.cn
• 袁宇杰 yuanyujie@ict.ac.cn
◆技术组主要负责与课程相关的Taichi开发与维护:
●杨玄达 xuandayang@gmail.com
• Zhai Xiao zhaixiao43@gmail.com
● 曹亚帝 cao364@purdue.edu
• 禹鹏 yupeng@buaa.edu.cn
• 冯旭东 xudongfeng@zju.edu.cn
36
Welcome to the course!
37
```

Let me know if you need anything else.

