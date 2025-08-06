Okay, I have extracted the text from the PDF. Due to the limitations of the tool, I had to perform OCR to extract the text, so there may be some errors. Here is the transcription:

Eulerian Fluid Simulation
Yuanming Hu
Overview
Grid
Advection
Projection
Solving
large-scale linear
systems
Eulerian Fluid Simulation
GAMES 201 Lecture 4
Yuanming Hu
MIT CSAIL
June 22, 2020

Today's topic: Eulerian fluid simulation
Eulerian Fluid
Simulation
Yuanming Hu
Overview
Grid
Advection
Projection
Solving
large-scale linear
systems
• Eulerian representation uses still sensors in space, usually arranged in a
regular grid/triangular mesh.
• A little bit of math - but not too much.
• This course: intuitive derivation - instead of finite volume/finite difference.
Recommended book
A great introduction to Eulerian fluid simulation:
Fluid simulation for computer graphics¹ by Robert Bridson.
1R. Bridson (2015). Fluid simulation for computer graphics. CRC press.

8
Table of Contents
Eulerian Fluid
Simulation
Yuanming Hu
1 Overview
Overview
Grid
Advection
Projection
Solving
large-scale linear
2 Grid
3 Advection
systems
4 Projection
5 Solving large-scale linear systems

Eulerian Fluid
Simulation
Yuanming Hu
Overview
Grid
Material Derivatives: Lagrangian v.s. Eulerian
D
Dt
д
:= + u.
dt
E.g.,
Advection
Projection
Solving
large-scale linear
systems
DT
Dt
=
Әт
dt
Dux дих
dt
Dt
=
+u∇T
+ u · ux
u: material (fluid) velocity. Many other names: Advective/Lagrangian/particle
derivative.
Intuitively, change of physical quantity on a piece of material
1 change due to time + (Eulerian).
② change due to material movement u 7.
=

Eulerian Fluid
Simulation
Yuanming Hu
(Incompressible) Navier-Stokes equations
Overview
Grid
Advection
Projection
Solving
large-scale linear
systems
Du
2
ρ
Dt
=
-∇p + µ∇²u + pg, ∇ • u = 0
or
Du
1
=
∇p + v∇²u + g
72u+g,
•
· u = 0
Dt
ρ
μ: dynamic viscosity; v = : kinematic viscosity.
Variants of the N-S equations
The Navier-Stokes equations have many variants. Here we show a version that is
the most friendly to fluids simulation in computer graphics. In graphics we usually
drop viscosity except for highly viscous materials (e.g., honey).

Eulerian Fluid
Simulation
Yuanming Hu
Overview
Grid
Operator splitting [More details]
1
=--∇p+g
ρ
Du
Dt
Split the equations above into three parts:
Advection
Projection
Solving
large-scale linear
systems
• u
0
Du
Da
=
0,
= 0 (advection)
(1)
Dt
Dt
ди
=
g (external forces, optional)
(2)
dt
ди
1
-
ps.t.
0 (projection)
(3)
dt
ρ
a: any physical property (temperature, color, smoke density etc.)

Eulerian Fluid
Simulation
Yuanming Hu
Eulerian fluid simulation cycle
Time discretization with splitting: for each time step,
1 Advection: "move” the fluid field. Solve u* using ut
Overview
Grid
Advection
Projection
Solving
large-scale linear
systems
Du
Da
= 0,
=
0
Dt
Dt
② External forces (optional): evaluate u** using u*
ди
dt
at = g (external forces, optional)
● Projection: make velocity field ut+1 divergence-free based on u**
ди
1
ps.t.
•
V. ut+1 = 0
dt
ρ

8
Table of Contents
Eulerian Fluid
Simulation
Yuanming Hu
1 Overview
Overview
Grid
Advection
Projection
Solving
large-scale linear
2 Grid
3 Advection
systems
4 Projection
5 Solving large-scale linear systems

Eulerian Fluid
Simulation
Yuanming Hu
Overview
Grid
Advection
Projection
Solving
large-scale linear
systems
Spatial discretization using cell-centered grids
Figure: ux, uy, p are all stored at the center (orange) of cells.
u
=
=
n, m 3,3
ti.var(ti.f32, shape=(n, m)) # x-component of velocity
V
=
ti.var(ti.f32, shape=(n, m)) # y-component of velocity
p
= ti.var(ti.f32, shape=(n, m)) # pressure

Eulerian Fluid
Simulation
Yuanming Hu
Overview
Grid
Advection
Projection
Solving
large-scale linear
systems
Spatial discretization using staggered grids
Figure: Red: ux; Green: uy; Blue: p.
n, m = 3,3
u
=
ti.var(ti.f32, shape=(n+1, m)) # x-component of velocity
V
=
ti.var(ti.f32, shape=(n, m+1)) # y-component of velocity
=
p ti.var(ti.f32, shape=(n, m)) # pressure

8
Eulerian Fluid
Simulation
Yuanming Hu
Overview
Grid
Advection
Projection
Solving
large-scale linear
systems
Bilinear interpolation
(X1,Y2) (X2,Y2)
(x,y)
(x1,y1)
(X2,Y1)
=0 +0
+
。
+0
Figure: Bilinear interpolation: value at (x, y) is a weighted average of the four corners.
Source: Wikepedia

8
Table of Contents
Eulerian Fluid
Simulation
Yuanming Hu
1 Overview
Overview
Grid
Advection
Projection
Solving
large-scale linear
2 Grid
3 Advection
systems
4 Projection
5 Solving large-scale linear systems

8
Eulerian Fluid
Simulation
Yuanming Hu
Overview
Grid
Advection
Projection
Solving
large-scale linear
systems
Advection schemes
A trade-off between numerical viscosity, stability, performance and complexity:
• Semi-Lagrangian advection²
• MacCormack/BFECC³
• "BiMocq2"4
• Particle advection (PIC/FLIP/APIC/PolyPIC, later in this course)
...
2R. Courant, E. Isaacson, and M. Rees (1952). "On the solution of nonlinear hyperbolic
differential equations by finite differences". In: Communications on pure and applied mathematics
5.3, pp. 243-255; J. Stam (1999). "Stable fluids". In: Proceedings of the 26th annual conference
on Computer graphics and interactive techniques, pp. 121–128.
3B. Kim et al. (2005). Flowfixer: Using BFECC for fluid simulation. Tech. rep. Georgia
Institute of Technology; A. Selle et al. (2008). "An unconditionally stable MacCormack method".
In: Journal of Scientific Computing 35.2-3, pp. 350-371.
4Z. Qu et al. (2019). "Efficient and conservative fluids using bidirectional mapping". In: ACM
Transactions on Graphics (TOG) 38.4, pp. 1–12.

Eulerian Fluid
Simulation
Yuanming Hu
Overview
Grid
Advection
Projection
Solving
large-scale linear
systems
Semi-Lagrangian advection
-
-
P
u(p) At
Figure: What should be the field value at p now based on the field and velocity at the
previous time step? Well, just let reverse the simulation...
@ti.func
def semi_lagrangian(x, new_x, dt):
for I in ti.grouped(x):
new_x[I]
=
sample_bilinear(x, backtrace(I, dt))

Eulerian Fluid
Simulation
Yuanming Hu
Overview
Grid
Advection
Projection
Solving
large-scale linear
systems
What if
...
p
Figure: The real trajectory of material parcels can be complex... Red: a naive estimation
of last position; Light gree: the true previous position.

Overview
Grid
Advection
Eulerian Fluid
Simulation
Yuanming Hu
Going back in time (Demo)
Initial value problem (ODE): simply use explicit time integration schemes, e.g.,
• Forward Euler (“RK1")
-=
p dt * velocity(p)
• Explicit Midpoint ("RK2")
Projection
p_mid
=
-
P 0.5 * dt * velocity(p)
Solving
large-scale linear
systems
-=
P dt * velocity(p_mid)
• RK3
v1 =
p1
=
-
v2
=
p2
=
-
v3
=
velocity(p)
P 0.5 * dt * v1
velocity(p1)
p 0.75 * dt * v2
velocity(p2)
-=
P
dt * (2 / 9 * v1 + 1 / 3 * v2 + 4 / 9 * v3)

Eulerian Fluid
Simulation
Yuanming Hu
Overview
Grid
Advection
Projection
Solving
large-scale linear
systems
BFECC and MacCormack advection schemes
BFECC: Back and Forth Error Compensation and Correction
• x* = SL(x,
**
x
t)
= SL(x*, -∆t)
• Estimate the error x
final
error 1
=
(x** - x)
• Apply the error x = x* + x
error
Be careful: need to prevent overshooting.
Demo!
@ti.func
def maccormack(x, dt):
semi_lagrangian(x, new_x, dt)
semi_lagrangian(new_x, new_x_aux, -dt)
for I in ti.grouped(x):
new_x[I]
=
new_x[I] + 0.5 * (x[I]
new_x_aux[I])

8
Table of Contents
Eulerian Fluid
Simulation
Yuanming Hu
1 Overview
Overview
Grid
Advection
Projection
Solving
large-scale linear
2 Grid
3 Advection
systems
4 Projection
5 Solving large-scale linear systems

ди
dt
Eulerian Fluid
Simulation
Yuanming Hu
Overview
Chorin-style projection
How to ensure the velocity field is divergence free after projection?
1
ρ
ps.t.
u = 0 (projection)
Grid
Expand (using finite difference in time):
Advection
Projection
At
Solving
u*
*
u
ps.t.
•
u = 0
(4)
large-scale linear
ρ
systems
At
*
*
u
u
ps.t.
· u
0
(5)
ρ
At
*
· u
=
· (u
p)
(6)
ρ
At
0
· u -
∇∇p
(7)
ρ
∇∇p
Pu
(8)
At

8
Eulerian Fluid
Simulation
Yuanming Hu
Overview
Grid
Advection
Projection
Solving
large-scale linear
systems
Poisson's equation
(From the previous slide)
p=
Pu
(9)
t
...
which is Poisson's equation
∇∇p = f or Ap = f.
(10)
2
is the Laplace operator. If f = 0, the equation is called
Laplace's equation.

Eulerian Fluid
Simulation
Yuanming Hu
Overview
Grid
Advection
Projection
Solving
large-scale linear
systems
Spatial discretization (2D)
Recall the equation for p:
Discretize on a 2D grid:
ρ
∇∇p=u
At
1
(Ap) ij = (∇∇p) i,j = 2 (-4Pinj+
X
(11)
22(-4pij + Pi+1,j + Pi−1,j + Pi,j−1 + Pi,j+1) (12)
(น+1,j
bj = () = (+15 -
bi,j
At
Again, a linear system:
u
i,j
ρ
AtAx
X
+ + - )
uj u u)
Anmxnmpnm = bnm
n, m: numbers of cells along the x- and y-axis.
i,j+1
i,j
(13)

Eulerian Fluid
Simulation
Yuanming Hu
Overview
Grid
Spatial discretization (2D): ∇ · u
•
uvi,j+1
uxi+1,j
Advection
Projection
Solving
large-scale linear
systems
uxi,j
uvi,j
1
Δα
X
X
(u) i,j = (u+1,j - uj + uj+1 – uj
i,j

Eulerian Fluid
Simulation
Yuanming Hu
Spatial discretization (2D): ∇ · ∇ p
•
Overview
Grid
Advection
Projection
Solving
large-scale linear
systems
Pi,j+1
Pi-1,j
Pi,j Pi+1,j
Pi,j-1
1
(∇∇p)ij =
△22 (-4Pi,j+Pi+1, j + Pi-1,j + Pi.j−1 + Pi,j+1)
Question: How to handle Dirichlet and Neumann boundaries?

8
Table of Contents
Eulerian Fluid
Simulation
Yuanming Hu
1 Overview
Overview
Grid
Advection
Projection
2 Grid
Solving
large-scale linear
systems
3 Advection
4 Projection
5 Solving large-scale linear systems

8
Eulerian Fluid
Simulation
Yuanming Hu
Overview
Grid
Advection
Projection
Solving
large-scale linear
systems
The top 10 algorithms from the 20th century
• 1946: The Metropolis Algorithm for Monte Carlo.
• 1947: Simplex Method for Linear Programming.
• 1950: Krylov Subspace Iteration Method.
• 1951: The Decompositional Approach to Matrix Computations.
• 1957: The Fortran Optimizing Compiler.
• 1959: QR Algorithm for Computing Eigenvalues.
• 1962: Quicksort Algorithms for Sorting.
• 1965: Fast Fourier Transform.
• 1977: Integer Relation Detection.
• 1987: Fast Multipole Method.

8
Eulerian Fluid
Simulation
Yuanming Hu
Overview
Grid
Advection
Projection
Solving
large-scale linear
systems
The top 10 algorithms from the 20th century
• 1946: The Metropolis Algorithm for Monte Carlo.
• 1947: Simplex Method for Linear Programming.
• 1950: Krylov Subspace Iteration Method.
• 1951: The Decompositional Approach to Matrix Computations.
• 1957: The Fortran Optimizing Compiler.
• 1959: QR Algorithm for Computing Eigenvalues.
• 1962: Quicksort Algorithms for Sorting.
• 1965: Fast Fourier Transform.
• 1977: Integer Relation Detection.
• 1987: Fast Multipole Method.

Eulerian Fluid
Simulation
Yuanming Hu
Solving large-scale linear systems
Many physics engines boil down to a (huge) linear system solve:
Overview
Grid
Advection
Projection
Solving
Ax = b
How to solve it:
• Direct solvers (e.g., PARDISO)
large-scale linear
systems
• Iterative solvers:
• Gauss-Seidel
• (Damped) Jacobi
• (Preconditioned) Krylov-subspace solvers (e.g., conjugate gradients)
Good numeric solvers are usually composed of different solvers: e.g.,
multigrid-preconditioned conjugate gradients with damped Jacobi smoothing and
PARDISO at the bottom multigrid level.

Matrix storage
Eulerian Fluid
Simulation
Yuanming Hu
Overview
Grid
Advection
Projection
Solving
large-scale linear
systems
What's special about A: often sparse, symmetric & positive-definite (SPD).
How to store A? Options:
1 As a dense matrix (e.g., float A[1024] [1024] doesn't scale but works)
2 As a sparse matrix (various sparse matrix formats: CSR, COO, ...)
• Don't store it at all (aka. Matrix-free, often the ultimate solution...)
Modern computer architecture: memory bandwidth is expensive but FLOPs are
free. So compute matrix entries on-the-fly (instead of fetching values from
memory) can sometimes be good to performance.

Eulerian Fluid
Simulation
Yuanming Hu
Overview
Grid
Advection
Projection
Solving
large-scale linear
systems
Krylov-subspace solvers
Krylov-subspace solvers are among most efficient linear system solvers. The most
well-known version: conjugate gradients (CG).
Less frequently used (in graphics):
• Conjugate residuals (CR)
• Generalized minimal residual method (GMRES)
• Biconjugate gradient stabilized (BiCGStab)
...
Recommended book
An Introduction to the Conjugate Gradient Method Without the Agonizing Pain5
by Jonathan Richard Shewchuk.
5J. R. Shewchuk et al. (1994). An introduction to the conjugate gradient method without the
agonizing pain.

Eulerian Fluid
Simulation
Yuanming Hu
Overview
Grid
Advection
Projection
Solving
large-scale linear
systems
Conjugate gradients
ro = b
Po = ro
k=0
while True:
ακ =
Axo
T
rkrk
T
Pk Apk
Xk+1 = xk + akPk
rk+1 = rk akApk
if ||rk+1|| is sufficiently small, break
T
rk+1rk+1
Bk
=
rT
rerk
Pk+1 = rk+1 + ẞkPk
k = k+1
return xk+1
Residual v.s. Error
Small residual r does not mean small error e, especially when A is poorly
conditioned (i.e., with a huge condition number (next slide)).

Eigenvalues and condition numbers
Eulerian Fluid
Simulation
Yuanming Hu
Recall that if
Overview
Grid
Advection
Ax = λχ,
then A is an eigenvalue of A and x is an eigenvector of A.
Projection
Solving
large-scale linear
systems
The condition number K of SPD matrix A:
κ(A) = Amax/min
In general: a smaller condition numbers means faster convergence.
(Note that condition numbers have many different definitions.)

Iterative solver trick: Warm starting
Eulerian Fluid
Simulation
Yuanming Hu
Overview
Grid
Advection
Projection
Solving
large-scale linear
systems
If you start with an initial guess that is close to the solution, very likely fewer
iterations are needed.
"Warm starting": use the p from the last frame as the initial guess of the current
frame.
Online demo
In practice works well for (damped) Jacobi/Gauss-Seidel/CG, but for MGPCG
(later in this lecture) it doesn't work well.

8
Eulerian Fluid
Simulation
Yuanming Hu
Overview
Grid
Preconditioning
Find an approximate operator M that is close to A but easier to invert. Then,
Advection
Ax = b
M¯¹Ax = M¯¹b
Projection
Solving
large-scale linear
systems
Intuition: M¯¹A may have a smaller condition number (closer to identity) or
better eigenvalue clustering than A itself.
Question: why not directly let M = A?

8
Eulerian Fluid
Simulation
Yuanming Hu
Overview
Grid
Advection
Projection
Solving
large-scale linear
systems
Common preconditioners
• Jacobi (diagonal) preconditioner M = diag(A)
• Poisson preconditioner
• (Incomplete) Cholesky decomposition
• Multigrid: M = very complex linear operator that almost inverts A...
• Geometric multigrid
• Algebraic multigrid
• Fast multipole method (FMM)

Eulerian Fluid
Simulation
Yuanming Hu
(Geometric) Multigrid methods
PHI
Multigrid V-Cycle: Solving PHI in PDE f(PHI) = F
F
R
R1
R2
R
Overview
Grid
Advection
Projection
Solving
large-scale linear
systems
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
Figure: Multigrid V-cycle (source: Wikipedia)
PHI

8
Eulerian Fluid
Simulation
Yuanming Hu
Overview
Grid
The Multigrid design space
1 Restriction/prolongation
● Cycle (V/W/F cycle)
• Smoothers ((red-black) Gauss-Seidel, Jacobi, damped Jacobi, etc.)
Advection
Projection
Solving
large-scale linear
systems
⑤ Bottom level solver (Brute-force Jacobi or direct solvers)
4 Number of levels (e.g., coarsen until the bottom level has 50K voxels)
⑥ Number of pre/post iterations (usually, 2-5)
⑦ Coarsening and boundary handling (e.g., Galerkin coarsening, semi-algebraic
8
multigrid)
...

Eulerian Fluid
Simulation
Yuanming Hu
Overview
Grid
Advection
Projection
Solving
large-scale linear
systems
Multigrid preconditioned conjugate gradients (MGPCG)
When solving Poisson's equation in graphics, people ususally use geometric
multigrid as the preconditioner for conjugate gradients.
Recommended reading
If you want to learn more about multigrid (and linear solvers in general):
• A multigrid tutorial6.
• A seminal and easy-to-understand multigrid paper in graphics: A parallel
multigrid Poisson solver for fluids simulation on large grids
Taichi demo: 2D/3D multigrd: ti example mgpcg_advanced
7
“W. L. Briggs, V. E. Henson, and S. F. McCormick (2000). A multigrid tutorial. SIAM.
7
A. McAdams, E. Sifakis, and J. Teran (2010). “A Parallel Multigrid Poisson Solver for Fluids
Simulation on Large Grids.". In: Symposium on Computer Animation, pp. 65–73.

Summary: Eulerian fluid simulation in graphics
Eulerian Fluid
Simulation
Yuanming Hu
For each time step,
• Advection: “move” the fluid field. Solve for u* using ut
Overview
Grid
Advection
Projection
Solving
large-scale linear
systems
Du
Da
=
0,
=
0
Dt
Dt
Key: Use a advection scheme with low numerical viscosity (e.g.,
MacCormack/BFECC/Particle advection)
• Projection: make velocity field ut+1 divergence-free based on u*
ди
1
ps.t.
V. ut+1 = 0
dt
ρ
Key: Use a fast linear solver (e.g., MGPCG).

Eulerian Fluid
Simulation
Yuanming Hu
Combining advection with reflection: IVOCK
SIGGRAPH 2015
Overview
Grid
Advection
Projection
Solving
large-scale linear
systems
◉◉◉
Original velocity field
Velocity field after
advection
Velocity field after
pressure projection
Figure: IVOCK: Restoring Missing Vortices in Advection-Projection Fluid Solvers
8X. Zhang, R. Bridson, and C. Greif (2015). "Restoring the missing vorticity in
advection-projection fluid solvers". In: ACM Transactions on Graphics (TOG) 34.4, pp. 1–8.

Eulerian Fluid
Simulation
Yuanming Hu
Advection-Reflection solver
SIGGRAPH 2018
1
u
1/2
1/2
1
Overview
Grid
Advection
Projection
Solving
large-scale linear
systems
u
u1
0
u
1
u
0
u
1
u1/2 u1
û1/2
Fig. 2. A geometric interpretation of our method. Left: In a standard advec-
tion-projection solver, projection to the divergence-free subspace causes
kinetic energy loss (red). Middle: Our reflection solver uses an energy-
preserving reflection (yellow) halfway through the advection step, dramati-
cally reducing the energy loss caused by the final projection. Our method
has effectively identical computational cost to an advection-projection solver
with half the time step (right), but loses less energy.
Figure: Source: An Advection-Reflection Solver for Detail-Preserving Fluid Simulation
9J. Zehnder, R. Narain, and B. Thomaszewski (2018). "An advection-reflection solver for
detail-preserving fluid simulation". In: ACM Transactions on Graphics (TOG) 37.4, pp. 1–8.

8
Possible extensions
Eulerian Fluid
Simulation
Yuanming Hu
Overview
Grid
Advection
Projection
Solving
large-scale linear
systems
• Going to 3D
• Accurate boundary conditions and fluid-solid coupling10
• Two phase fluid simulation11
• Handling free surfaces (level sets)12
• Vortex methods13
A well-implemented Eulerian fluid solver counts as Homework 1 :-)
10C. Batty, F. Bertails, and R. Bridson (2007). "A fast variational framework for accurate
solid-fluid coupling". In: ACM Transactions on Graphics (TOG) 26.3, 100–es.
11R. Ando, N. Thuerey, and C. Wojtan (2015). “A stream function solver for liquid
simulations". In: ACM Transactions on Graphics (TOG) 34.4, pp. 1–9.
12S. Osher, R. Fedkiw, and K Piechor (2004). "Level set methods and dynamic implicit
surfaces". In: Appl. Mech. Rev. 57.3, B15–B15.
13X. Zhang and R. Bridson (2014). "A PPPM fast summation method for fluids and beyond”.
In: ACM Transactions on Graphics (TOG) 33.6, pp. 1–11.

