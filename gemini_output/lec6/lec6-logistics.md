I have transcribed the text from the PDF. Since the PDF contains images, the transcription may not be perfect.

Page 1:
GAMES 201
Advanced Physics Engines 2020: A Hands-on Tutorial
高级物理引擎实战2020 (基于太极编程语言)
第六讲:线性弹性有限元与拓扑优化
Yuanming Hu 胡渊鸣 MIT CSAIL
麻省理工学院 计算机科学与人工智能实验室 Taichi Programming Language

Page 2:
Homework 1 获奖选手
十余个用Taichi实现的(半)隐式时间积分器
https://zhuanlan.zhihu.com/p/158962220
2

Page 3:
作业安排
Homework 1截止
●没有赶上deadline/没有获奖的同学不用灰心
• 加把劲一个月后就可以作为作业2提交
Homework 2为最终作业:实现一个自己满意的物理模拟器
● 选项1:可交互的2D物理模拟游戏
●选项2:优化性能,提高粒子数/网格精度
●选项3:实现高精度格式(Advection-reflection等)
●选项4:...
●可以基于自己的或别人的Homework 1(评分标准是原创部分)
●建议组队:-)
●北京时间8月15日23:59 截止
3

Page 4:
课程安排
◆第六讲,7月13日 线性弹性有限元与拓扑优化
◆第七讲,7月20日 混合欧拉-拉格朗日视角(1)
◆第八讲,7月27日 混合欧拉-拉格朗日视角(2)
◆第九讲,8月3日 高性能计算与物理引擎
◆8月10日空一次,大家完善自己的物理引擎(开放作业2)
◆第十讲,8月17日 总结
4

