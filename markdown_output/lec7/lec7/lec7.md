# lec7


## Page 1

⾼高级物理理引擎实战2020 Advanced Physics Engines 2020: A Hands-on Tutorial （基于太极编程语⾔言） Yuanming Hu 胡渊鸣 MIT CSAIL 麻省理理⼯工学院 计算机科学与⼈人⼯工智能实验室 1 GAMES 201 第六讲：混合欧拉-拉格朗⽇日视⻆角(1)

![page1_img0.png](images/page1_img0.png)

![page1_img1.png](images/page1_img1.png)

### Tables


## Page 2

Schedule ✦The following two lectures will be focused on hybrid Lagrangian-Eulerian methods ✦July 20: ๏Overview ๏Particle advection schemes: ‣ Particle-in-cell (PIC), Afﬁne PIC (APIC), Polynomial PIC (Poly PIC) ‣ Fluid Implicit Particles (FLIP) ๏Material Point Method basics ✦July 27: ๏Moving Least Squares MPM (MLS-MPM), theory and implementation ๏Constitutive models in MPM ๏Lagrangian forces in MPM ๏Implicit MPM ๏Advanced Taichi features 2

### Tables


## Page 3

Lagrangian v.s. Eulerian: Two Views of Continuums 3

### Tables


## Page 4

Lagrangian View Sensors that move passively with the simulated material “What are my position and velocity?” 4 https://pixabay.com/photos/paper-boat-coloured-colored-2770974/

![page4_img0.png](images/page4_img0.png)

### Tables


## Page 5

Eulerian View https://www.peakpx.com/15/gray-wooden-pillar-lot-on-body-of-water Still sensors that never moves “What is the material velocity passing by?” 5

![page5_img0.png](images/page5_img0.png)

### Tables


## Page 6

Deformable Body Simulation Lagrangian representation Eulerian representation [Macklin et al. 2014, Uniﬁed Particle Physics for Real-Time Applications] [Fedkiw 2001, Visual Simulation of Smoke] 6

![page6_img0.png](images/page6_img0.png)

![page6_img1.png](images/page6_img1.png)

![page6_img2.png](images/page6_img2.png)

### Tables


## Page 7

7 Which is better?

### Tables


## Page 8

Key factors to consider… I.e., deﬁne “better”: ✦Conservation of physical quantities ๏Momentum ๏Angular momentum ๏Volume (incompressibility) ๏Energy (low dissipation) ๏… ✦Performance (parallelism & locality on modern hardware) ✦Complexity ✦… 8

### Tables


## Page 9

Hybrid Eulerian-Lagrangian Schemes 9

### Tables


## Page 10

Motivation 10 ✦Recall that a ﬂuid solver usually has two components: ๏Advection (evolving the ﬁelds) ๏Projection (enforcing incompressibility) ✦Eulerian grids are really good at projection: ๏Easy to discretize ๏Efﬁcient neighbor look-up ๏Easy to precondition (geometric multigrid) ✦But Eulerian grids are bad at advection… ๏Dissipative: loss of energy and geometry

![page10_img0.png](images/page10_img0.png)

![page10_img1.png](images/page10_img1.png)

### Tables


## Page 11

Motivation 11 ✦Lagrangian particles are good at advection ๏Simply move their coordinates :-) ๏More conservative (lower dissipation) ✦But projection on particles can be tricky: ๏Tricky to discretize ๏Need complex data structures for neighbor look-up ✦Can we somehow smartly combine Lagrangian particles and Eulerian grids?

### Tables


## Page 12

Eulerian Grids (3) Grid to Particle transfer G2P (1) Particle to Grid transfer P2G (4) Particle operations: ๏Move particles ๏Update material ๏… Lagrangian Particles (2) Grid operations: ๏Pressure projection ๏Boundary conditions ๏… (which stores most of the information) (often auxiliary) 12

![page12_img0.png](images/page12_img0.png)

![page12_img1.png](images/page12_img1.png)

### Tables

Table 1:

| (3
(2) Grid operations:
Grid to Part
Pressure projection G
๏
Boundary conditions
๏
…
๏
(1
Particle to G
P2
Eulerian
Grids
(often auxiliary)
1 | )
icle transfer
2P
(4) Particle operations:
Move particles
๏
)
rid transfer Update material
๏
G
…
๏
Lagrangian
Particles
(which stores most
of the information)
2 |

Table 2:

| (3
Grid to Part
G | )
icle transfer
2P |

Table 3:

| (1
Particle to G
P2 | )
rid transfer
G |

Table 4:

| Eulerian
Grids
(often auxiliary) |  |


## Page 13

Particle-in-cell Particle to Grid, P2G Grid to Particle, G2P 13 Harlow, F.H. (1964) The Particle-in-Cell Computing Method for Fluid Dynamics. Methods in Computational Physics, 3, 319-343. velocity, temperature, force, …

![page13_img0.png](images/page13_img0.png)

![page13_img1.png](images/page13_img1.png)

![page13_img2.png](images/page13_img2.png)

![page13_img3.png](images/page13_img3.png)

![page13_img4.png](images/page13_img4.png)

![page13_img5.png](images/page13_img5.png)

![page13_img6.png](images/page13_img6.png)

![page13_img7.png](images/page13_img7.png)

![page13_img8.png](images/page13_img8.png)

![page13_img9.png](images/page13_img9.png)

![page13_img10.png](images/page13_img10.png)

![page13_img11.png](images/page13_img11.png)

![page13_img12.png](images/page13_img12.png)

![page13_img13.png](images/page13_img13.png)

![page13_img14.png](images/page13_img14.png)

![page13_img15.png](images/page13_img15.png)

![page13_img16.png](images/page13_img16.png)

![page13_img17.png](images/page13_img17.png)

![page13_img18.png](images/page13_img18.png)

![page13_img19.png](images/page13_img19.png)

![page13_img20.png](images/page13_img20.png)

![page13_img21.png](images/page13_img21.png)

![page13_img22.png](images/page13_img22.png)

![page13_img23.png](images/page13_img23.png)

![page13_img24.png](images/page13_img24.png)

![page13_img25.png](images/page13_img25.png)

![page13_img26.png](images/page13_img26.png)

![page13_img27.png](images/page13_img27.png)

![page13_img28.png](images/page13_img28.png)

![page13_img29.png](images/page13_img29.png)

![page13_img30.png](images/page13_img30.png)

![page13_img31.png](images/page13_img31.png)

![page13_img32.png](images/page13_img32.png)

![page13_img33.png](images/page13_img33.png)

![page13_img34.png](images/page13_img34.png)

![page13_img35.png](images/page13_img35.png)

![page13_img36.png](images/page13_img36.png)

### Tables


## Page 14

The particle does not treat neighbors equally Closer = more importance 14

![page14_img0.png](images/page14_img0.png)

![page14_img1.png](images/page14_img1.png)

### Tables


## Page 15

B-Spline Kernels N(x) 15 https://www.seas.upenn.edu/~cffjiang/research/mpmcourse/mpmcourse.pdf Cubic Quadratic Linear

![page15_img0.png](images/page15_img0.png)

![page15_img1.png](images/page15_img1.png)

![page15_img2.png](images/page15_img2.png)

![page15_img3.png](images/page15_img3.png)

### Tables


## Page 16

PIC P2G code (transfer velocity) 16 assuming cell-centered grid (0, 0) (2.5dx, 0.5dx)

![page16_img0.png](images/page16_img0.png)

### Tables


## Page 17

PIC grid normalization code 17

![page17_img0.png](images/page17_img0.png)

### Tables


## Page 18

PIC G2P code (gather velocity) 18

![page18_img0.png](images/page18_img0.png)

### Tables


## Page 19

Combing PIC and grid-based Poisson solver ✦Simulation cycle ๏P2G: scatter velocity from particles to grid ๏normalize velocity ๏Pressure projection ๏G2P: gather divergence-free velocity from grid to particles ‣ Move particles according to the divergence-free velocity ﬁeld ‣ Use RK2/3/4 if necessary ✦Demo: http://yuanming-hu.github.io/ﬂuid/ (FLIP blending=0) ✦Does it look good? 19

### Tables


## Page 20

Let’s run a very simple PIC simulation… 20

### Tables


## Page 21

PIC G2P: information loss (Assume we have only 1 particle) Figure from Fu et al 2017, A Polynomial Particle-In-Cell Method (SIGGRAPH Asia 2017) Problem: 18 Do Fs on grid, 2 Do Fs on particle 21

![page21_img0.png](images/page21_img0.png)

### Tables


## Page 22

Reducing particle-in-cell dissipation ✦Two solutions: 1.Transfer more information: APIC, Poly PIC 2.Transfer the delta: FLIP (later in this lecture) 22

### Tables


## Page 23

Afﬁne Particle-in-cell (APIC) Jiang et al., SIGGRAPH & JCP 2016 23 Figure from Fu et al 2017, A Polynomial Particle-In-Cell Method (SIGGRAPH Asia 2017) v0 v1 C00 C11 C10 C01

![page23_img0.png](images/page23_img0.png)

### Tables


## Page 24

Homework (highly recommended!) ✦Watch Bilibili: https://www.bilibili.com/video/BV1St411G7nm/

![page24_img0.png](images/page24_img0.png)

### Tables


## Page 25

APIC conserves angular momentum! 25 https://www.seas.upenn.edu/~cffjiang/research/apic/tech-doc.pdf An angular momentum conserving afﬁne-particle-in-cell method, Jiang et al., JCP 2017

![page25_img0.png](images/page25_img0.png)

![page25_img1.png](images/page25_img1.png)

### Tables


## Page 26

APIC P2G, G2P 26

![page26_img0.png](images/page26_img0.png)

![page26_img1.png](images/page26_img1.png)

### Tables


## Page 27

Let’s run it… 27

### Tables


## Page 28

Poly PIC A Polynomial Particle-In-Cell Method, Fu et al. 2017 18 modes=9 nodes X 2 Do Fs per node: Lossless transfer! Figures from Fu et al 2017, A Polynomial Particle-In-Cell Method (SIGGRAPH Asia 2017) 28

![page28_img0.png](images/page28_img0.png)

![page28_img1.png](images/page28_img1.png)

### Tables


## Page 29

Fluid implicit particles (FLIP) 29 ✦Idea: don’t gather the physical quantity. Gather the delta of the physical quantities before/after grid operation. ๏grid op = pressure projection in incompressible ﬂuid simulation ๏grid op = internal force computation in solid simulation (MPM) ✦Note: some VFX people may use “FLIP” for “a ﬂuid solver using Chorin-Style pressure projection with FLIP advection”, but FLIP itself is just an advection scheme. Yongning Zhu and Robert Bridson. 2005. Animating Sand As a Fluid. SIGGRAPH 2005 BRACKBILL, J. U., AND RUPPEL, H. M. 1986. FLIP: a method for adaptively zoned, particle-in-cell calculations of ﬂuid ﬂows in two dimensions. JCP

### Tables


## Page 30

Combining PIC and FLIP ✦PIC is too dissipative, yet FLIP is too noisy: can we interpolate between the two methods? ✦FLIP0.99 = 0.99 * FLIP + 0.01 * PIC ✦Demo: http://yuanming-hu.github.io/ﬂuid/ (adjust FLIP blending) 30 vt+1 p = gather(vt+1 i ) vt+1 p = vt p + gather(vt+1 i −vt i)

### Tables


## Page 31

Figure from Qu et al., Efﬁcient and Conservative Fluids Using Bidirectional Mapping, SIGGRAPH 2019 ✦Suggestion: start with APIC. PIC is almost never used in graphics ✦APIC is easy to implement ✦APIC is storage friendly (no backup velocity ﬁeld compared with FLIP) ✦APIC preserves angular momentum ✦APIC is stable ✦APIC leads to nice visual results ✦APIC is the basis of MLS-MPM! PIC/FLIP/APIC/Poly PIC: which one to use? 31 Bimocq

![page31_img0.png](images/page31_img0.png)

### Tables


## Page 32

Material Point Method (MPM) 32

### Tables


## Page 33

Material Point Method (MPM) ✦Hybrid Lagrangian-Eulerian simulation scheme ๏Not just “advection” ๏Particle carries a lot more than velocity ✦Very hot research topic: at least 5 papers at SIGGRAPH 2020 ✦Invented by Sulsky & Schreyer in 1996 ✦First introduced to graphics in 2013 ๏A material point method for snow simulation (Stomakhin et al.) 33

### Tables


## Page 34

The Material Point Method as of 2018 Particles (Constitutive models) Snow [Stomakhin et al. 2013], Foam [Ram et al. 2015, Yue et al. 2015] Sand [Klar et al. 2015, Pradhana et al 2017] Grid SPGrid [Setaluri et al. 2014], Open VDB [Museth 2013] Multiple Grids [Pradhana et al. 2017] Particle to Grid (P2G) Grid to Particle (G2P) Transfer (Particle-in-Cell, PIC) Afﬁne PIC, APIC [Jiang et al. 2016] Polynomial PIC, Poly PIC [Fu et al. 2017] High-performance GIMP [Gao et al. 2017] Moving Least Squares [Hu et al. 2018] Compatible PIC [Hu et al. 2018] … 34

### Tables

Table 1:

|  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |


## Page 35

[Stomakhin et al. 2013, A material point method for snow simulation] The Material Point Method (MPM) 35

![page35_img0.png](images/page35_img0.png)

### Tables


## Page 36

The Material Point Method (MPM) ✦MPM is popular because of … ๏Automatic coupling of different materials (liquids, solids, granular materials etc.) ๏Automatic (self-)collision handling ๏Automatic fracture ๏Capable of simulating large deformations ✦Hybrid Lagrangian-Eulerian: both a grid and particles are used ๏An Eulerian grid is used for collision handling and momentum update ๏Lagrangian particles are used for state tracking such as advection and deformation 36

### Tables


## Page 37

The Material Point Method (MPM) Stomaching et al., A material point method for snow simulation, SIGGRAPH 2013 37

![page37_img0.png](images/page37_img0.png)

![page37_img1.png](images/page37_img1.png)

### Tables


## Page 38

Classical MPM in graphics 38 Stomaching et al., A material point method for snow simulation, SIGGRAPH 2013 https://www.math.ucla.edu/~jteran/papers/SSCTS13.pdf

![page38_img0.png](images/page38_img0.png)

![page38_img1.png](images/page38_img1.png)

### Tables


## Page 39

Moving Least Squares MPM ✦Based on APIC ✦Halves the required FLOPs (2x faster!) ✦Much easier to implement than traditional MPM ๏88 lines of code using Taichi ✦Demos: ๏ti example mpm88 ๏ti example mpm99 ๏ti example mpm128 39 A Moving Least Squares Material Point Method with Displacement Discontinuity and Two-Way Rigid Body Coupling Hu et al, SIGGRAPH 2018 http://taichi.graphics/wp-content/uploads/2019/03/mls-mpm-cpic.pdf

### Tables


## Page 40

Implementing MLS-MPM 40

### Tables


## Page 41

The MLS-MPM Simulation Cycle Particle 2 grid (P2G) Grid Op Grid 2 particle (G2P) For each particle: Update particles using (afﬁne) velocity; Scatter mass & momentum to nearby 3x3x3 nodes. For each grid node: divide momentum by mass to get velocity; apply gravity and boundary conditions. For each particle: Gather velocity/afﬁne velocity from 3x3x3 nodes. MLS- MPM Cycle Bandwidth-saving version in Hu2019Taichi 41

### Tables


## Page 42

MLS-MPM in 88 Lines of Taichi Code 42

![page42_img0.png](images/page42_img0.png)

### Tables


## Page 43

MLS-MPM in 88 Lines of Taichi Code 43

![page43_img0.png](images/page43_img0.png)

### Tables


## Page 44

MLS-MPM in 88 Lines of Taichi Code 44

![page44_img0.png](images/page44_img0.png)

### Tables


## Page 45

MLS-MPM in 88 Lines of Taichi Code 45 C00 C11 C10 C01

![page45_img0.png](images/page45_img0.png)

![page45_img1.png](images/page45_img1.png)

### Tables


## Page 46

MLS-MPM in 88 Lines of Taichi Code 46

![page46_img0.png](images/page46_img0.png)

### Tables


## Page 47

Recap ✦Hybrid Eulerian-Lagrangian schemes ๏Use particles to track material (position, velocity, deformation) ๏Use grids to compute force ﬁelds (Chorin projection/Cauchy stress) ✦Reducing dissipation of PIC: ๏APIC/Poly PIC (more modes) ๏FLIP (gather delta) ✦Material Point Method ๏Use particles to store deformation information 47

### Tables


## Page 48

MPM Courses/Paper list ✦The material point method for simulating continuum materials ๏Chenfanfu Jiang, Craig Schroeder, Joseph Teran, Alexey Stomakhin, and Andrew Selle ๏In ACM SIGGRAPH 2016 Courses (SIGGRAPH ’16) ✦On hybrid Lagrangian-Eulerian simulation methods: practical notes and high-performance aspects ๏Yuanming Hu, Xinxin Zhang, Ming Gao, and Chenfanfu Jiang ๏In ACM SIGGRAPH 2019 Courses (SIGGRAPH ’19) ✦MPM in computer graphics by Chenfanfu Jiang 48

### Tables


## Page 49

The end Questions are welcome! 49

### Tables

