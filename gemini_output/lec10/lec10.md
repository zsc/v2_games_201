Okay, I have processed the file. Here is the transcribed text from lec10.pdf:

```
GAMES 201
Advanced Physics Engines 2020: A Hands-on Tutorial
高级物理引擎实战2020
(基于太极编程语言)
第十讲:总结
Yuanming Hu
胡渊鸣
MIT CSAIL
Taichi
麻省理工学院 计算机科学与人工智能实验室
Programming Language
Inflow speed: 0.5
Inflow speed: 1.0
Wheel density: 1.0
Wheel density: 4.0
Wheel density: 4.0
Homework 2 获奖作品
• Wimaxs & g1n0st: Taichi小游戏合集
BillXu2000:整合mpm lagrangian_forces和mpm_99
JYLeeLYJ: OOP Taichi and Eulerian Fluid Engine Design
hejob: 2D FVM Compressible CFD Solver for multiblock structured mesh
Timber: MPM 3D Liquid
huahuo: PBF-driven game (Water gun)
Iqxu: 3D PBF
Jerry Yan: Eulerian Fluid -- Jacobi method and CG method
linyaodong: 3D SPH with Grid-based Neighborhood Search
robs: 2D Free Surface Fluid Simulation (Eluerian/PIC/FLIP/APIC)
limengfan: MPM 的三种边界条件
♦课程结束后欢迎继续提交!(但是可能得等比较久才能收到纪念品...)
Capture
Inflow
0.5
1.5
2.5
Slip Wall
Outflow
ChinaVR Taichi竞赛
◆一分钟完成报名:https://www.wjx.top/jq/88136866.aspx
◆特等奖1名、一等奖3名、二等奖6名、三等奖和优秀奖若干
♦ 竞赛网站:https://www.chinavr.info/creative.html
中国计算机学会
https://www.chinavr.info/
中国计算机学会
中国图象图形学学会
CCF
China Computer Federation
ESP
CHINA SOCIETY OF IMAGE AND GRAPHICS
中国仿真学会
China Simulation Federation
恭喜各位!
♦ 收获了新技能:Taichi编程语言 刚体 液体 烟雾弹塑性体 PIC/FLIP法
Krylov-子空间求解器 预条件 无矩阵法 多重网格 弱形式与有限元 隐式积分
器拓扑优化物质点法 现代处理器微架构 内存层级并行编程 GPU编程稀
疏数据结构 可微编程...
♦ 课程难度不小,能坚持到最后的都是勇者!
太极人专用点赞
Overview
◆ What I learned from building physics engines and Taichi
• Simplicity is good!(“简单”好!)
◆物理引擎的未来
◆我自己在课程中的收获
◆ Acknowledgments
Simplicity is good
♦考试:证明我会某种方法做研究、工程:我能找到好方法
◆用哪种方法不重要,重要的是系统性地解决问题
◆一天只有24个小时
●“最聪明的人都是一样聪明的”
●用简单可行的方法解决大部分简单的问题,才有时间去解决真正复杂
的问题
◆有影响力的创新往往伴随着复杂性的减少
●“不要给他人带来麻烦”
Complexity is bad. Why?
◆ No one else understands your ideas (limited impact)
Slow to implement (won't fail quickly)
Error-prone and hard to analyze
◆ Hard to implement efficiently (by hands and by compilers)
• Sparse data structures (as C++ libraries)
• As a result: slow physical solvers
Why do people still go complicated?
Why do people still go complicated?
Because it's way easier to be complicated than to be simple.
Why do people still go complicated?
To prove smartness
• Maybe a good idea at school...
...
but almost never a good idea when you need to interact with people
• In fact if you look at the bigger picture it's often better not to be too smart
To show off the amount of work
• so that you have more pages in your paper-not a great idea either
To hide potential errors
• "Nobody can prove I'm wrong because nobody understands what I'm doing, so I'm right."
●“科学的表述在于其是可证伪的”- 卡尔·波普尔
To differentiate from prior work (and have the paper accepted)
• "Simple ideas have been done by previous people, so I have to do a more complex one."
Being lazy to simplify
Audience being too nice to tolerate
Use a simple solution != being an idiot
◆The smartest people know the importance of
keeping things simple.
"Everything should be made as simple as possible, but not simpler."
Albert Einstein
cppcon+
Morgan Stanley
Father of C++
Make simple things simple!
What does "simple" mean?
Everything should be made
as simple as possible,
but not simpler.
MAKE SIMPLE TASKS SIMPLE!
Bjarne Stroustrup
Albert Erutein
Make the code directly express intent
-
Simple code is easier to understand
Simple code is often beautiful
Simple code is often fast
Not all code can be simple
•
Hide the complexity behind a simple interface
Stroustrup - Simple Cppcon'14
8
cppcon
MAKE SIMPLE TASKS SIMPLE!
Bjarne Stroustrup
Morgan Stanley
Being too clever is not clever
BWK (1976)
- "debugging is twice as hard as writing the program in the first place. So, if
you're as clever as you can be when you write it, how can you debug it?"
Design for clarity
-
-
Essential complexity
.
.
What is needed to defeat to solve the problem
Try to find the simplest solution
Accidental complexity
introduced by the expressions of the solution
Minimize it by using C++ well
Don't over-abstract
.
Base generalizations on concrete solutions
Hide complexity
.
Good interfaces are essential
Stroustrup-Simple - Cppcon'14
7
"There are two ways of constructing a software design: One way is to make it
so simple that there are obviously no deficiencies, and the other way is to
make it so complicated that there are no obvious deficiencies. The first
method is far more difficult."
C.A.R. Hoare
“有两种软件设计的方法:第一种是使得软件足够简单以至于显然没有错误,
第二种是使得软件足够复杂以至于没有显然的错误。第一种方法难得多。”
——Tony Hoare (1980年图灵奖得主,快速排序发明人)
How to Have a Bad Career by David Patterson (Turing Award 2017)
Gazte
David Patterson
How to Have a Bad Career
Bad Career Move #2: Let Complexity Be Your Guide
(Confuse Thine Enemies)
• Best compliment:
"Its so complicated, I can't understand the ideas"
• Easier to claim credit for subsequent good ideas
If no one understands, how can they contradict your claim?
• It's easier to be complicated
Also: to publish it must be different; N+1st incremental change
• If it were not unsimple then how could distinguished
colleagues in departments around the world be positively
appreciative of both your extraordinary intellectual grasp of
the nuances of issues as well as the depth of your
contribution?
The story of ChainQueen [ICRA 2019]
ChainQueen aims to be a differentiable MPM solver
But we didn't have a tailored AutoDiff system at that time
• I had to manually derive the gradients of MPM
• Fortunately we have MLS-MPM, which is a lot simpler than
traditional MPM
→ Thanks to the simplicity of MLS-MPM, we had a ICRA paper
finished in 25 days (otherwise it won't be doable)
• DiffTaichi further simplifies differentiable MPM - no manual
gradient derivation anymore!
ChainQueen's supplemental material...
(J, G2P) For x, we have
n+1
=
x + Atv+1
vn+1
=
(85)
(86)
Cn+1
i
4
-1
=
ΔΣΝ(x – x) (xi - x)
(87)
P
=
N(xi
-
mpv +
4
Δ2
pp
T + mp Cm (xi - x)
(88)
p
P
4
m
=
∑N(xi - x) mp
Gp := (-VPF +mpC)
(89)
(90)
←
(91)
ƏL
дха
ƏL
дх
n+1
Ldx+1
+
Lv+1
+
2+1
+
+1 +1 +1 д
OL dp JLdm
+
дриддт д
(92)
α
(93)
VB β
ƏL
ρβα
-
ΣΣ
OL
θχία
LN(x – x)
θχία
+ΣΣ+1
ί β
+ΣΣΔ
i
ί β
ί β
4
+mp am
i
pß
L [N(x – x)
n
β
N(x – x)
θχία
(mpvp3 + [Gp(xi – x)]s) – N(xi - x)GpBa
(94)
LN(x-x)
2n+1
pβγ
θχία
Vi (Xiy - Xpy)}
(95)
-
(96)
(97)
Another example: Taichi and MGPCG
PHI
Multigrid V-Cycle: Solving PHI in PDE f(PHI) = F
F
R
R1
R2
R
PHI
Gauss Seidel
Compute esid
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
MGPCG (multigrid preconditioned conjugate gradients) is hard to implement.
With Taichi, a programming system that aims to be simple, I can implement it in 80 minutes and 300 LoC(according to git log).
After GAMES 201, now everyone knows how to write a multigrid solver!
很少有人能正确实现-> 我能轻易正确实现-> 大家都能轻易正确实现(大家好才是真的好)
Recall: our ambitious goal is to
...
Simulate the world in your computer!
19
物理引擎的未来
♦易用性
●如何快速、正确设置参数?
◆健壮性 (robustness)
●手抖一下模拟炸了
◆完备性
●多材料耦合
◆真实感
obj crowd_sim
Crowd Solver crowdsolver1
Start Frame 1
Particle Behavior Animation Behavior Substeps
Orientation Constraint
Min Speed 0.25
Speed Constraint
Min Angle 0
Max Angle 90
Drag 1
Spin Drag 1
Update Particles
Up Constraint
Source Attribute
Avoidance Look At Terrain Attributes
Agent Group
Weight 1
Particle Scale Multipl... 1
Force Scale 1.5
Anticipation Time 10
Neighbor Distance
10
Max Neighbors 100
Agent Avoidance
+
H
物理引擎的未来
♦易用性
●如何快速、正确设置参数?
◆健壮性(robustness)
●手抖一下模拟炸了
◆完备性
●多材料耦合
◆真实感
物理引擎的未来
♦易用性
●如何快速、正确设置参数?
◆健壮性(robustness)
●手抖一下模拟炸了
◆完备性
●多材料耦合
◆真实感
物理引擎的未来
♦易用性
●如何快速、正确设置参数?
◆健壮性(robustness)
●手抖一下模拟炸了
◆完备性
●多材料耦合
◆真实感
Orthotropic
物理引擎的未来
◆性能和引擎制造的生产力是关键
●高性能计算、领域特定编程语言和编译器(Taichi)必不可少
• 自动adaptive grids/多GPU/多node simulation
• Simplicity, simplicity, simplicity
Learning-for-Simulation
◆ Simulation-for-Learning
我自己在课程中的收获
◆同学们的作业(每天醒来第一件事...)
◆世界上第一次用中文开设的物理引擎课?
●近千名同学参加,近五十份引擎项目提交
● Taichi的贡献:复制,粘贴,运行!
●同学们通过GAMES 201和Taichi快速掌握了曾今的我折腾了好久才弄懂的知识
●第一次普及了MGPCG, MLS-MPM, APIC等实用算法
♦讲课还是挺不容易的
●周末备课像打仗(“一周只工作五天原来是有道理的!”)
●自己懂了!=能流利的讲出来
● LaTeX真难用(但是公式排版效果还是好)
◆直播技术升级了。学会了正确使用麦克风。
TONOR
Pors
POST
POPS
00
200
200
POVS
1 Cardioid pickup
0
330
30
300
60
270
240
90
120
2
210
180
↑
150
Built-in quality
chipset
TONOR
鸣谢
GAMES社区:刘利刚教授(USTC),闫令琪教授(UCSB),周晓巍教授(ZJU); 陈凌昊、
张远青、董峻廷同学
Professors who (directly or indirectly) taught me so much on physical simulation:
●蒋陈凡夫教授(UPenn)
●朱博教授(Dartmouth)
• Prof. Eftychios Sifakis (UWisconsin-Madison)
• Prof. Joseph M. Teran (UCLA)
Prof. Ron Fedkiw (Stanford)
Ph.D. advisors at MIT:
Prof. Fredo Durand
Prof. Bill Freeman
鸣谢
◆每一位来上课的同学
●特别是交了作业的同学:-)
♦ Taichi论坛管理员:woclass
♦ Taichi核心开发组:archibate, k-ye, xumingkuan, Rullec,
FantasyVR, TH3CHARLie, rexwangcc, KLozes, Eydcao,
...
♦ Taichi中文社区: StephenArk30, archibate, ArkhamWJZ, isdanni,
rexwangcc, liaopeiyuan, zhiyaluo
♦ 后勤部:Linan Zhou
感谢助教同学们的鼎力相助!
●曹亚帝
●禹鹏
●杨玄达
●冯旭东
◎翟骁
●夏一鸣
●史雨宸
●袁宇杰
•PLY dumper,
• compiler bugs,
•ti.field,
• Export video,
•v0.6.7->v0.6.27,
....
26
Bo Zhu @ Dartmouth
• Research Interests: Physics-based Simulation, Computational Fluid
Dynamics, Scientific Machine Learning
• Grabbing my research and teaching materials at:
http://www.dartmouth.edu/~boolzhu/
• If you are interested in pursuing a PhD/MS/Intern/Visitor at
Dartmouth, please reach out via: bo.zhu@dartmouth.edu
Codimension-0
2+1 = 2n+1 + δх
position correction
C
s=
Ax-SVC
||C||
Iterative Neural Projection
DARTMOUTH
热烈祝贺GAMES 201结课!!!附宾夕法
尼亚大学丁IANG GROUP招生宣传
Machine
Learning
Numerical
Optimization
High
Performance
Computing
Physics-Based Modeling
Xc
Scientific
Computing
Π
and Simulation for :
1. Animation and Visual
Applications
2. Computational
Physics
3. Robotics
t = 183 s
CHICAGOO
GEFORCE GTX
GEFORCE GTX
GEFORCE GTX
GEFORCE GTX
欢迎PhD申请、访问学者/联合培养研究生、大
三暑假的本科生、远程合作者、以及相信基于物理
的模拟与数字孪生将促进第四次工业革命的你
4. Visual Scene
Understanding
5. And more for you to
开拓…!
https://www.seas.upenn.edu/~cffjiang/
More resources for graduate study
Toshiya Hachisuka's Graduate Study Survival Guide
◆ Fredo Durand's advices
• http://people.csail.mit.edu/fredo/student.html
• How to Write a Bad Paper
Bill Freeman's advices
• How to do research
课程预告:GAMES 102 - 几何建模与处理
·主讲:刘利刚,中国科学技术大学
• http://staff.ustc.edu.cn/~Igliu
·预期时间:2020年10-12月
·课程定位:几何建模与处理的基础入门
・内容概要:
・离散微分几何基础理论
•几何建模:重建、网格化、编辑、变形等
•几何处理:去噪、光顺、简化、压缩、参数化、修复等
・形状理解和分析基础
・前序课程:
・计算机图形学、C++编程
谢谢大家!
谢谢大家!
「孔子曰:三人行,则必有我师。是故弟子不必不如师,师
不必贤于弟子,闻道有先后,术业有专攻,如是而已。」
————韩愈《师说》
谢谢大家!
「孔子曰:三人行,则必有我师。是故弟子不必不如师,师
不必贤于弟子,闻道有先后,术业有专攻,如是而已。」
————韩愈《师说》
「师父领进门,修行在个人。」
谢谢大家!
「孔子曰:三人行,则必有我师。是故弟子不必不如师,师
不必贤于弟子,闻道有先后,术业有专攻,如是而已。」
————韩愈《师说》
「师父领进门,修行在个人。」
期待十年以后,你的物理引擎课程!
```

