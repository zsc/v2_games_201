Okay, I have extracted the text from the PDF. Due to the nature of the tool, I had to OCR the document. Here is the transcribed text:

GAMES 201
Advanced Physics Engines 2020: A Hands-on Tutorial
高级物理引擎实战2020 (基于太极编程语言)
第四讲:欧拉视角下的流体模拟
Yuanming Hu 胡渊鸣
MIT CSAIL
麻省理工学院 计算机科学与人工智能实验室
☑ Taichi
Programming Language
Inflow speed: 0.5
Inflow speed: 1.0
Wheel density: 1.0
Wheel density: 4.0
Wheel density: 4.0

Homework 0
新作品
3D explicit FEM by yucrazing
BarnsleyFern 巴恩斯利蕨 by ElementMo
• MultiBarnsleyFern by shadowalker
3D vortex ring leapfrogging by citadel
Nodes: 20
Elements: 24

醉酒的蝴蝶飞不出花花的世界的原因找到了!
by 中科院物理所
Homework0 计算流体力学视角的流体求解器
by zwang
有了流体力学的理论之后,我们就可以用它进行模拟,从而康康流体流速的具体分布情况。
模拟流体的算法非常多,像有限差分(FDM),有限元(FEM),光滑粒子动力学(SPH),物质
点法(MPM)等等。这里的模拟是用了格子玻尔兹曼方法(LBM)。
层流(非静止画面)
卡门涡街
湍流
注:自上而下,雷诺数增大,图中颜色越深,表示流速越小。
● 从模拟的结果可以看出,雷诺数小的时候,流体形成稳定的层流,这对应着实验中风速较小
•
• 的情况。此时,只有来流自身的不稳定所造成的小波动。当雷诺数变大时,流体经过圆柱,
• 开始变得不稳定起来,形成周期交替出现的涡,即卡门涡街。当雷诺数进一步增大时,流体
• 的运动更加复杂,开始形成湍流。实验中风速较大的情况就对应着这两种情形,纸条的波动
会非常大。

Homework 1:
隐式时间积分器
◆ 在物理模拟中,显式时间积分器容易实现,但是通常数值上较为不稳定,对时间步长 (dt)更为敏感。而隐式时间积分器
较难实现(通常需要求解线性系统、多次迭代等操作),但是允许较大的时间步长。
◆ 本次作业中,大家可以选择自己最喜欢的模拟器,实现隐式时间积分。从性能(模拟同样时长的物理系统所需的运行
时间)、准确度(如能量守恒)、数值稳定性等角度进行对比。
* 可选模拟器:
* implicit mass-spring/FEM, PCI-SPH, DF-SPH, MPS, Eulerian fluid (with pressure projection)...
* 对比:
* 显式 v.s. 隐式(如semi-implicit Euler v.s. backward Euler)
* 或隐式 v.s. 别的隐式(如Jacobi iteraiton v.s. conjugate gradients)
♦ 建议组队
• 分工合作(GitHub)
● 相互验算隐式时间积分器的公式
提交格式(论坛提交GitHub链接和一些关键图表)
• README.md 中加入性能(如每帧计算时间)、准确度
4
(机械能守恒)、稳定性的分析(允许的dt)

Logistics
Taichi v0.6.12
♦论坛功能升级(感谢@woclass)
◆第五讲,6月29日:多体问题与涡方法 客座讲师:张心欣[知乎]
• 多体问题以及他们与柏松方程的联系
• 涡方法的乐趣
● 从直观的角度引导同学认识几种不同的快速求和方法
• 涡方法 Demo
♦7月6日,空一周,实现开放作业1
◆7月11日,开放作业1截止,点评
5

课程纪念品
Homework 0已经邮寄~
获奖同学请查收
import taichi as ti
Taichi
Programming Language

Lagrangian View
https://pixabay.com/photos/paper-boat-coloured-colored-2770974/
Sensors that move passively
with the simulated material
"What are my position and velocity?"
7

Eulerian View
Today's topic
"What is the material velocity passing by?"
Still sensors that never moves
https://www.peakpx.com/15/gray-wooden-pillar-lot-on-body-of-water
8

