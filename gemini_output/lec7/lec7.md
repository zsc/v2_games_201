Okay, I have transcribed the text from lec7.pdf. Due to the limitations of the tool, I was only able to extract the text through OCR, which might not be perfect, but I have done my best to provide an accurate transcription.
```text
GAMES 201
Advanced Physics Engines 2020: A Hands-on Tutorial
高级物理引擎实战2020 (基于太极编程语言)
第六讲:混合欧拉-拉格朗日视角(1)
Yuanming Hu 胡渊鸣 MIT CSAIL
麻省理工学院 计算机科学与人工智能实验室 Taichi
Programming Language
Inflow speed: 0.5
Inflow speed: 1.0
Wheel density: 1.0
Wheel density: 4.0
Wheel density: 4.0

Schedule
The following two lectures will be focused on hybrid Lagrangian-Eulerian methods
July 20:
• Overview
• Particle advection schemes:
Particle-in-cell (PIC), Affine PIC (APIC), Polynomial PIC (PolyPIC)
►
►
Fluid Implicit Particles (FLIP)
• Material Point Method basics
July 27:
• Moving Least Squares MPM (MLS-MPM), theory and implementation
• Constitutive models in MPM
• Lagrangian forces in MPM
• Implicit MPM
• Advanced Taichi features 2

Lagrangian v.s. Eulerian:
Two Views of Continuums
3

Lagrangian View https://pixabay.com/photos/paper-boat-coloured-colored-2770974/
Sensors that move passively
with the simulated material
"What are my position and velocity?"
4

Eulerian View "What is the material velocity passing by?"
Still sensors that never moves
https://www.peakpx.com/15/gray-wooden-pillar-lot-on-body-of-water
5

Deformable Body Simulation
Lagrangian representation
[Macklin et al. 2014,
Unified Particle Physics for Real-Time Applications]
6
W
Eulerian representation
[Fedkiw 2001, Visual Simulation of Smoke]

Which is better?
7

Key factors to consider...
I.e., define “better":
◆ Conservation of physical quantities
• Momentum
• Angular momentum
• Volume (incompressibility)
• Energy (low dissipation)
...
Performance (parallelism & locality on modern hardware)
+ Complexity
8

Hybrid Eulerian-Lagrangian Schemes
9

Motivation
◆ Recall that a fluid solver usually has two components:
• Advection (evolving the fields)
• Projection (enforcing incompressibility)
◆ Eulerian grids are really good at projection:
• Easy to discretize
• Efficient neighbor look-up
• Easy to precondition (geometric multigrid)
But Eulerian grids are bad at advection...
• Dissipative: loss of energy and geometry
10

Motivation
Lagrangian particles are good at advection
• Simply move their coordinates :-)
• More conservative (lower dissipation)
◆ But projection on particles can be tricky:
• Tricky to discretize
• Need complex data structures for neighbor look-up
Can we somehow smartly combine Lagrangian particles and
Eulerian grids?
11

(3)
(2) Grid operations:
• Pressure projection
• Boundary conditions
Grid to Particle transfer
G2P
•
U03
U13
U23
U33
U02
U12
U22
U32
U01
U11
U22
U32
400
U10
U20
U30
(1)
Particle to Grid transfer
Eulerian
Grids
(often auxiliary)
P2G
11
(4) Particle operations:
• Move particles
• Update material
Lagrangian
Particles
(which stores most
of the information)

Particle-in-cell
Particle to Grid, P2G
velocity,
temperature,
force,
Harlow, F.H. (1964) The Particle-in-Cell Computing Method for Fluid Dynamics.
Methods in Computational Physics, 3, 319-343.
13
Grid to Particle, G2P

The particle does not treat neighbors equally
Ñ(x)
(x)N
Particle to Grid, P2G
Closer = more importance
14

B-Spline Kernels N(x)
Linear
N(x) =
{
10
1-|x| 0 ≤ |x| < 1
|x| > L
Quadratic
-1x12
2
0 < x <
1
(-x)²½≤x<
N(x) =
1/3
0
>>
< |x|
N(x)
3
{}
Cubic
11x|3-1x|2+3 0 ≤ x < 1
N(x) = (2-|x|)3
0
1 ≤ |x| < 2
2 ≤ |x|
:
15 https://www.seas.upenn.edu/~cffjiang/research/mpmcourse/mpmcourse.pdf

PIC P2G code (transfer velocity)
assuming
cell-centered grid
for p in x:
base = (x[p] * inv_dx - 0.5).cast(int)
fx = x[p] * inv_dx - base.cast(float)
# Quadratic B-spline
W =
[0.5 * (1.5 fx) ** 2, 0.75 - (fx 1) ** 2, 0.5 * (fx 0.5) ** 2]
for i in ti.static(range(3)):
for j in ti.static(range(3)):
offset = ti. Vector([i, j])
weight = w[i][0] * w[j][1]
-
grid_v[base + offset] += weight * v[p]
grid_m[base + offset] += weight
16
(0, 0)
(2.5dx, 0.5dx)

PIC grid normalization code
54
for i, j in grid_m:
55
if grid_m[i, j] > 0:
56
57
inv_m = 1 / grid_m[i, j]
grid_v[i, j] = inv_m * grid_v[i, j]
17

PIC G2P code (gather velocity)
59
for p in x:
60
61
62
63
64
65
66
67
68
69
70
71
72
73
base = (x[p] * inv_dx 0.5).cast(int)
fx = x[p] * inv_dx - base.cast(float)
# Quadratic B-spline
W =
[
0.5 * (1.5
-
fx) ** 2, 0.75
-
(fx
1.0) ** 2, 0.5 * (fx
-
0.5) ** 2
]
new_v = ti. Vector.zero(ti.f32, 2)
for i in ti.static(range(3)):
for j in ti.static(range(3)):
weight = w[i][0] * w[j][1]
new_v += weight * grid_v [base + ti. Vector([i, j])]
x[p] = clamp_pos(x[p] + v[p] * dt)
v [p] = new_v
18

Combing PIC and grid-based Poisson
solver
◆ Simulation cycle
• P2G: scatter velocity from particles to grid
• normalize velocity
• Pressure projection
• G2P: gather divergence-free velocity from grid to particles
Move particles according to the divergence-free velocity field
Use RK2/3/4 if necessary
♦ Demo: http://yuanming-hu.github.io/fluid/ (FLIP blending=0)
◆ Does it look good?
19

Let's run a very simple PIC
simulation...
20

PIC G2P: information loss
(Assume we have only 1 particle)
constant
•
n
p
On
P
α = 1
α = 2
11 = 0, 12 = 0
11 = 0, 12 = 0
Problem: 18 DoFs on grid, 2 DoFs on particle
Figure from Fu et al 2017, A Polynomial Particle-In-Cell Method (SIGGRAPH Asia 2017)
21

Reducing particle-in-cell dissipation
◆Two solutions:
1.Transfer more information: APIC, PolyPIC
2.Transfer the delta: FLIP (later in this lecture)
22

Affine Particle-in-cell (APIC)
Jiang et al., SIGGRAPH & JCP 2016
Vo
V1
Coo
C10
C01
C11
constant
linear
On
X
Ap
On
p
On
p
n
On
X
p
n
α = 1
α = 2
α = 1
11 = 0, 12 = 0
11 = 0, i2 = 0
11 = 1, 12 = 0
α = 2
1 = 1, 12 = 0
α = 1
11 = 0, і2 = 1
α = 2
11 = 012 = 1
Figure from Fu et al 2017, A Polynomial Particle-In-Cell Method (SIGGRAPH Asia 2017)
23

Homework (highly recommended!)
APIC
affine particle in cell
Chenfanfu Jiang
Craig Schroeder
Andrew Selle
Joseph Teran
Alexey Stomakhin
University of California, Los Angeles
Walt Disney Animation Studios
ACMSIGGRAPH
Watch Bilibili: https://www.bilibili.com/video/BV1St411G7nm/

APIC conserves angular momentum!
5.3.1 Particle to grid
G,n
tot
Proposition 5.4. Angular momentum is conserved during the APIC transfer from particles to the grid. Lor = Lor under transfer 3.
Proof. The angular momentum on the grid after transferring from particles is
G,n
tot
=
=ΣΣx mp + B (D)-1 (x - 2))
=
=
P
P
3=
=
=
Σmp
X
P
P
Σmp
xxn
X
p
P
P
X
P
p
-1
T
P
P
T
: € -
mp B
P
T
: €
T
)T
: €
Σ xxmpv+mp B (D)-(x) - x(x)
= ∑xxmpv + mp (B(D)-1D): €
P
P
= Σx × mpv + mp(B): €
p
= LPn
tot
P
-1
T
: €
5.3.2 Grid to particle
Proposition 5.5. Angular momentum is conserved during the APIC transfer from the grid to particles. Lon+1 = LGn+1 una
Proof. As before, the manipulation (vuT)T : ∈ = u x v is used to convert the permutation tensor into a cross product. Th
on the particles after transferring from the grid is
LP,n+1
tot
= x+1 x mpv+1 + Σmp(Bn+1) : €
P
=Σ+1
P
P
P
T
+1
mp
2+1(x – x)T
: €
P
i
(xx) x +1
X
p
P
i
Σ
+1
X
P
= x+1 x mpv+1+mp
X
X
= 2+1 x mpon+1 - Σmp x + 1 + 2mp wax + 1
X
P
P
= Σα+1 xmpv+1 - Σmpxx+1+Σ Σmpw+1
P
X
P
Vi
i
P
= Σx+1 x mpv+1 - mpxxv+1 + x x vn+1
p
P
i
= ∑(x+1 – x) × mpv+1 + xxmon+1
p
P
P
i
G,n+1
= ΣΔτυ+1 x mpvn+1 + Lor
P
= LG,n+1
https://www.seas.upenn.edu/~cffjiang/research/apic/tech-doc.pdf 25
An angular momentum conserving affine-particle-in-cell method,
Jiang et al., JCP 2017

APIC P2G, G2P
for p in x:
base = (x[p] * inv_dx 0.5).cast(int)
fx = x[p] * inv_dx base.cast(float)
# Quadratic B-spline
w = [
0.5 * (1.5 fx) ** 2, 0.75
(fx
1.0) ** 2, 0.5 * (fx 0.5) ** 2
]
0.5) ** 2]
for p in x:
base = (x[p] * inv_dx 0.5).cast(int)
fx = x[p] * inv_dx base.cast(float)
# Quadratic B-spline
W = [0.5 * (1.5 fx) ** 2, 0.75 (fx 1) ** 2, 0.5 * (fx
affine = C[p]
for i in ti.static (range(3)):
for j in ti.static (range(3)):
offset = ti. Vector([i, j])
dpos = (offset.cast(float)
weight = w[i] [0] * w[j] [1]
fx) * dx
grid_v [base + offset] += weight * (v[p] + affine @ dpos)
grid_m[base + offset] += weight
26
new_v = ti. Vector.zero(ti.f32, 2)
new_C = ti.Matrix.zero(ti.f32, 2, 2)
for i in ti.static(range(3)):
for j in ti.static(range(3)):
dposti. Vector((i, j)).cast(float)
g_v = grid_v [base + ti. Vector([i, j])]
weight = w[i] [0] * w[j][1]
new_v += weight * g_v
fx
new_C += 4 * weight * g_v.outer_product (dpos) * inv_dx
x[p] = clamp_pos(x[p] + new_v * dt)
v [p] = new_v
C [p] = new_C

Let's run it...
27

PolyPIC
A Polynomial Particle-In-Cell Method, Fu et al. 2017
18 modes=9 nodes X 2 DoFs per node: Lossless transfer!
constant
linear
FLIP
APIC PolyPIC
n
α = 1
110, 120
bilinear
α = 1
11 = 1, 12 = 1
α = 2
11 = 0, 12 = 0
n
n
t = Os
p
α = 1
11 = 1,120
α = 2
1 = 1,12 = 0
α =
11 = 0, 12 = 1
α = 2
1012= 1
quadratic
t = 1s
xn
t = 2.6s
α = 2
α = 1
2
α =
α = 2
11 = 1, 12 =
11 = 2, 12 = 0
11 = 2, 12 =
11 = 0,12 = 2
11 = 0,12 = 2
P
t = 6.3s
α = 1
α = 2
α = 1
11 = 2, 12 = 1
21 = 2, 12 = 1
11 = 1, 12 = 2
α = 2
21 = 1, 12 = 2
α = 1
α2
스
11 = 2, 12 = 2
11 = 2, 12 = 2
t = 16.7s
Figures from Fu et al 2017, A Polynomial Particle-In-Cell Method (SIGGRAPH Asia 2017)
28

Fluid implicit particles (FLIP)
Idea: don't gather the physical quantity. Gather the delta of the
physical quantities before/after grid operation.
• grid op = pressure projection in incompressible fluid
simulation
• grid op = internal force computation in solid simulation
(MPM)
Note: some VFX people may use "FLIP" for “a fluid solver using Chorin-Style pressure projection
with FLIP advection", but FLIP itself is just an advection scheme.
BRACKBILL, J. U., AND RUPPEL, H. M. 1986. FLIP: a method for adaptively zoned, particle-in-cell calculations of fluid flows in two dimensions. JCP
Yongning Zhu and Robert Bridson. 2005. Animating Sand As a Fluid. SIGGRAPH 2005
29

Combining PIC and FLIP
t
υ
+1 = gather(+1) +1 = v + gather (v+1 – v)
υπ
p
t
υ
p
p
◆ PIC is too dissipative, yet FLIP is too noisy: can we interpolate
between the two methods?
FLIP0.99 = 0.99 * FLIP + 0.01 * PIC
♦ Demo: http://yuanming-hu.github.io/fluid/ (adjust FLIP
blending)
30

PIC/FLIP/APIC/PolyPIC: which one to use?
Suggestion: start with APIC. PIC is almost never used in graphics
APIC is easy to implement
APIC is storage friendly (no backup velocity field compared with FLIP)
APIC preserves angular momentum
APIC is stable
APIC leads to nice visual results
SemiLag
BFECC
MacCormack
FLIP
SSSS
APIC is the basis of MLS-MPM!
APIC
PolyPIC
MC+R
Bimocq
SS
Figure from Qu et al., Efficient and Conservative Fluids Using Bidirectional Mapping, SIGGRAPH 2019
31

Material Point Method (MPM)
32

Material Point Method (MPM)
→ Hybrid Lagrangian-Eulerian simulation scheme
• Not just "advection"
• Particle carries a lot more than velocity
◆ Very hot research topic: at least 5 papers at SIGGRAPH 2020
Invented by Sulsky & Schreyer in 1996
First introduced to graphics in 2013
• A material point method for snow simulation (Stomakhin et al.)
33

The Material Point Method as of 2018
Particle to Grid (P2G)
Grid to Particle (G2P)
Transfer (Particle-in-Cell, PIC)
Affine PIC, APIC [Jiang et al. 2016]
Polynomial PIC, PolyPIC [Fu et al. 2017]
High-performance GIMP [Gao et al. 2017]
Moving Least Squares [Hu et al. 2018]
Compatible PIC [Hu et al. 2018]
Particles (Constitutive models)
Snow [Stomakhin et al. 2013],
Foam [Ram et al. 2015, Yue et al. 2015]
Sand [Klar et al. 2015, Pradhana et al 2017]
34
...
Grid
SPGrid [Setaluri et al. 2014],
OpenVDB [Museth 2013]
Multiple Grids [Pradhana et al. 2017]

The Material Point Method (MPM)
Figure 9: Walking character. A character stepping through snow produces an interesting structure. Disney.
[Stomakhin et al. 2013, A material point method for snow simulation]
35

The Material Point Method (MPM)
MPM is popular because of
...
• Automatic coupling of different materials (liquids, solids, granular
materials etc.)
• Automatic (self-)collision handling
• Automatic fracture
• Capable of simulating large deformations
→ Hybrid Lagrangian-Eulerian: both a grid and particles are used
• An Eulerian grid is used for collision handling and momentum update
• Lagrangian particles are used for state tracking such as advection and
deformation
36

The Material Point Method (MPM)
Particle Domain (Lagrangian)
1
Particle
states
only once
3
2
Particle
volumes
7
8
9
10
Updated deformation Particle velocities Collided particles Updated positions
4
5
gradients
6
Grid velocities
and mass
Grid View (Eulerian)
Material Point
Grid forces
Grid velocities
(RHS)
Collided
grid velocities
Implicitly
solved
velocities
Method Overview
Stomaching et al., A material point method for snow simulation, SIGGRAPH 2013
37

Classical MPM in graphics
1. Rasterize particle data to the grid. The first step is to trans-
fer mass from particles to the grid. The mass is transferred us-
ing the weighting functions m² = ∑p mpwip. Velocity also
should be transferred to the grid, but weighting with wip does
not result in momentum conservation. Instead, we use nor-
malized weights for velocity v = ∑pumpwp/m². This
contrasts with most incompressible FLIP implementations.
2. Compute particle volumes and densities. First timestep
only. Our force discretization requires a notion of a particle's
volume in the initial configuration. We can estimate a cell's
density as mi mi/h³, which we can weight back to the particle
as pp = Σi miwip/h³. We can now estimate a particle's
volume as Vp = mp/pp.
0
3. Compute grid forces using equation (6) with 2i = xi.
4. Update velocities on grid to v using equation (10).
5. Grid-based body collisions on vi as described in Section 8.
6. Solve the linear system in equation (9) for semi-implicit inte-
gration. For explicit time integration, simply let v+1 = v.
T
7. Update deformation gradient. The deformation gradient for
each particle is updated as Fn+1 = (I + ∆t∇v+1)Fm,
where we have computed ∇vn+1 = Συ+1(ρ).
Section 7 gives a detailed description of the update rule for
elastic and plastic parts of F.
8. Update particle velocities.
are vn+1 = (1 − a)PIC +
n+1
=
n
n+1
Vi Wip
n+1
i
Our new particle velocities
FLIPp, where the PIC part is
n+1
n+1
UFLIPP = v +
UPICP Συ+wi and the FLIP part is v
Σ(+1 – ). We typically used a = 0.95.
9. Particle-based body collisions on vn+1 as detailed in Sec-
tion 8.
n
10. Update particle positions using x+1 = x + ∆tv+1.
Stomaching et al., A material point method for
snow simulation, SIGGRAPH 2013
https://www.math.ucla.edu/~jteran/papers/SSCTS13.pdf
38

Moving Least Squares MPM
A Moving Least Squares Material Point Method with Displacement Discontinuity and Two-Way Rigid Body Coupling
Hu et al, SIGGRAPH 2018 http://taichi.graphics/wp-content/uploads/2019/03/mls-mpm-cpic.pdf
Based on APIC
◆ Halves the required FLOPs (2x faster!)
→ Much easier to implement than traditional MPM
• 88 lines of code using Taichi
Demos:
• ti example mpm88
⚫ ti example mpm99
• ti example mpm128
39

Implementing MLS-MPM
40

The MLS-MPM Simulation Cycle
Bandwidth-saving version in Hu2019Taichi
For each particle:
Update particles using (affine) velocity;
Scatter mass & momentum to nearby 3x3x3 nodes.
MLS- Particle 2 grid (P2G)
MPM Grid Op
Cycle
Grid 2 particle (G2P)
For each grid node:
divide momentum by mass to get velocity;
apply gravity and boundary conditions.
For each particle:
Gather velocity/affine velocity from 3x3x3 nodes.
41

MLS-MPM in 88 Lines of Taichi Code
1 import taichi as ti
2
import random
3
ti.init(arch=ti.gpu)
4
5
dim = 2
6
n_particles = 8192
7
n_grid = 128
8
dx = 1 / n_grid
9
inv_dx = 1 / dx
10
dt = 2.0e-4
11
p_vol = (dx * 0.5)**2
12
p_rho = 1
13
p_mass = p_vol * p_rho
14
E = 400
15
16
x = ti. Vector(dim, dt=ti.f32, shape=n_particles)
17
v = ti. Vector(dim, dt=ti.f32, shape=n_particles)
18
C = ti.Matrix(dim, dim, dt=ti.f32, shape=n_particles)
19
J = ti.var(dt=ti.f32, shape=n_particles)
20
grid_v = ti. Vector(dim, dt=ti.f32, shape=(n_grid, n_grid))
21 grid_m = ti.var(dt=ti.f32, 48hape=(n_grid, n_grid))

MLS-MPM in 88 Lines of Taichi Code
23
@ti.kernel
24
def substep():
25
for p in x:
26
27
fx = x[p] * inv_dx
-
base.cast(float)
28
W =
-
29
30
31
base = (x[p] * inv_dx - 0.5).cast(int)
[0.5 * (1.5 fx) ** 2, 0.75- (fx – 1) ** 2, 0.5 * (fx 0.5) ** 2]
stress = -dt * p_vol * (J[p]
affine = ti. Matrix(([stress, 0], [0, stress]]) + p_mass * C[p]
for i in ti.static(range(3)):
-
-
-
1) * 4 * inv_dx * inv_dx * E
32
for j in ti.static(range(3)):
33
offset = ti. Vector([i, j])
34
dpos = (offset.cast(float)
-
fx) * dx
35
weight = w[i][0] * w[j][1]
36
grid_v [base + offset] += weight * (p_mass * v[p] + affine @ dpos)
37
grid_m[base + offset] += weight * p_mass
43

MLS-MPM in 88 Lines of Taichi Code
39
for i, jin grid_m:
40
if grid_m[i, j] > 0:
41
bound = 3
42
43
44
45
46
inv_m = 1 / grid_m[i, j]
grid_v[i, j] = inv_m * grid_v[i, j]
grid_v[i, j][1] -= dt * 9.8
if i < bound and grid_v[i, j][0] < 0:
grid_v[i, j][0] = 0
47
if i > n_grid
-
bound and grid_v[i, j] [0] > 0:
48
49
50
grid_v[i, j][0] = 0
if j < bound and grid_v[i, j][1] < 0:
grid_v[i, j][1] = 0
51
if j > n_grid
-
bound and grid_v[i, j] [1] > 0:
52
grid_v[i, j][1] = 0
44

MLS-MPM in 88 Lines of Taichi Code
54
55
56
57
W =
[
58
-
-
for p in x:
base = (x[p] * inv_dx - 0.5).cast(int)
fx = x[p] * inv_dx base.cast(float)
0.5 * (1.5 fx) ** 2, 0.75 (fx -
1.0) ** 2, 0.5 * (fx - 0.5) ** 2
59
]
60
new_v = ti.Vector.zero(ti.f32, 2)
61
new_C = ti.Matrix.zero(ti.f32, 2, 2)
62
for i in ti.static(range(3)):
63
for j in ti.static(range(3)):
64
65
66
67
68
69
70
71
72
dpos = ti. Vector((i, j)).cast(float)
g_v = grid_v [base + ti. Vector([i, j])]
weight = w[i][0] * w[j] [1]
new_v += weight * g_v
fx
new_C += 4 * weight * g_v.outer_product(dpos) * inv_dx
v [p] = new_v
x[p] += dt * v[p]
J[p] *= 1 + dt * new_C.trace()
C[p] = new_C
linear
Coo
C10
C01
C11
On
p
p-
P
α = 1
11 = 1,12 = 0
α = 2
1 = 1,12 = 0
α =
11 = 0, 12 = 1
α = 2
1012=1
A
45

MLS-MPM in 88 Lines of Taichi Code
74
for i in range(n_particles):
75
x[i]
=
[random.random() * 0.4 + 0.2, random.random() * 0.4 + 0.2]
76
v[i] = [0, -1]
77
J[i] = 1
78
79
80
81
gui = ti.GUI("MPM88", (512,512))
for frame in range(20000):
for s in range(50):
82
grid_v.fill([0, 0])
83
grid_m.fill(0)
84
substep()
85
86
87
gui.clear(0x112F41)
gui.circles(x.to_numpy(), radius=1.5, color=0x068587)
gui.show()
46

Recap
Hybrid Eulerian-Lagrangian schemes
• Use particles to track material (position, velocity, deformation)
• Use grids to compute force fields (Chorin projection/Cauchy
stress)
Reducing dissipation of PIC:
• APIC/PolyPIC (more modes)
• FLIP (gather delta)
Material Point Method
Use particles to store deformation information
47

MPM Courses/Paper list
The material point method for simulating continuum materials
• Chenfanfu Jiang, Craig Schroeder, Joseph Teran, Alexey Stomakhin, and Andrew
Selle
• In ACM SIGGRAPH 2016 Courses (SIGGRAPH '16)
On hybrid Lagrangian-Eulerian simulation methods: practical notes and
high-performance aspects
• Yuanming Hu, Xinxin Zhang, Ming Gao, and Chenfanfu Jiang
In ACM SIGGRAPH 2019 Courses (SIGGRAPH '19)
MPM in computer graphics by Chenfanfu Jiang
48

The end
Questions are welcome!
49
```
