I have extracted the text from the PDF. Due to the limitations of the tool, I can only provide the text recognized via OCR, and some formatting may be lost. The OCR results are as follows:

Finite Elements and Topology Optimization
GAMES 201 Lecture 6
Yuanming Hu
MIT CSAIL
July 13, 2020

Table of Contents
1 FEM Overview
2 Discretizing Poisson's equation
3 Discretizing linear elasticity
4 Topology optimization

Finite element method
Finite element method (FEM) belongs to the family of Galerkin methods. In FEM, continuous PDEs are converted to discrete (linear) systems.
Typical steps:
1 Convert strong-form PDEs to weak forms, using a test function w.
2 Integrate by parts to redistribute gradient operators.
3 Use the divergence theorem to simplify equations and enforce Neumann boundary conditions (BCs).
4 Discretization (build the stiffness matrix and right-hand side).
5 Solve the (discrete) linear system.
Understanding FEM is important for many other discretization methods, including the Material Point Method (MPM) later in this course.

2D Poisson's equation
u(x) = f(x), x ∈ ∂Ω
∇ · ∇u = 0
∇u(x) · n = g(x), x ∈ ∂Ω
Application: pressure projection in fluid simulations.

Weak formulation
Arbitrary 2D test function w(x):
∇ · ∇u = 0 <=> ∀w, ∫∫ w(∇ · ∇u)dA = 0
Intuitively:
1 =>: trivial
2 <=: if ∇ · ∇u(x) != 0, we can always construct a test function w(x) s.t. ...

Getting rid of second-order derivatives (I)
We want to get rid of ∇ · ∇ in
∇ · ∇v = 0.
Recall integration by parts, or derivative of products:
∇w · ∇u + w∇ · ∇u = ∇ · (w∇u)
Since ∇ · ∇v = 0, we have
∇w · ∇u = ∇ · (w∇u)
To summarize,
∇ · ∇u = 0 <=> ∀w, ∫∫ ∇w · ∇udA - ∫ ∇ · (w∇u)dA.

Getting rid of second-order derivatives (II)
∫∫ ∇w · ∇udA = ∫∫ ∇ · (w∇u)dA
Divergence theorem applied to the RHS:
∫∫ ∇w · ∇udA = ∮ w∇u · dn

Discretization (I) Basis functions
u(x) = ∑ ujφj(x)
Figure: Compusing 1D piece-wise linear functions using 1D Basis functions.
In this course we focus on linear basis functions, which means the fields are exactly linear/bilinear/trilinear interpolated versions of the degrees of freedoms uj. (Higher order basis functions are rarely used in graphics.)

Discretization (I) Basis functions
Figure: 2D basis functions on a triangular mesh. Source: https://www.comsol.com/multiphysics/finite-element-method

Discretization (1)
u(x) = f(x), x ∈ ∂Ω
∇ · ∇u = 0
∇u(x) · n = g(x), x ∈ ∂Ω
Figure: We use rectangular (quadrilateral) finite elements. Use your imagination to visualize the basis functions :-)

Discretization (II)
Now we represent u(x) as
u(x) = ∑ ujφj(x),
Recall that we wan to solve for u s.t.
∫∫ ∇w · ∇udA = ∮ w∇u · dn,
Simple substitution gives
∀w, ∫∫ ∇w · ∇(∑ ujφj)dA = ∮ w∇u · dn.
It's sufficient to use only the basis function φi as test functions w:
∀i, ∫∫ ∇φi · ∇(∑ ujφj)dA = ∮ φi∇u · dn.

Discretization (III)
∀i, ∫∫ ∇φi · ∇(∑ ujφj)dA = ∮ φi∇u · dn.
Extract ∑ uj out of ∫∫:
∀i, ∑ (∫∫ ∇φi · ∇φj dA)uj = ∮ φi∇u · dn
In matrix form...
Ku = f
K: "stiffness" matrix; u: degree of freedoms/solution vector; f: load vector

Discretization (IV)
Now we need to compute Kij = ∫∫ ∇φi · ∇φj dA. Here we are using a simply basis function so it's not hard to compute analytically. (In more difficult cases people use Gaussian quadrature).

9-point Laplacian stencil
[ -1 -1 -1
-1 8 -1
-1 -1 -1 ]

Recall the 5-point Laplacian stencil we obtained using finite difference in previous lectures:
[ 0 1 0
1 -4 1
0 1 0 ]

Boundary Conditions
Recall that our linear system is
∀i, ∫∫ ∇φi · ∇φj dA uj = ∮ φi∇u · dn
Two types of boundary conditions
1 Dirichlet boundary conditions u(x) = f(x), x ∈ ∂Ω.
Easy: directly set corresponding ui = f(xi).
2 Neumann boundary conditions ∇u(x) · n = g(x), x ∈ ∂Ω
Plug g into the RHS of the equation, which yields non-zero entries in f.

Linear elasticity FEM
Cauchy momentum equation:
ρ Dv/Dt = ∇ · σ + g
v: velocity
ρ: density
σ: Cauchy stress tensor (symmetric 2/3D “matrix")
g: body force (e.g., gravity)
Quasistatic state (v = 0), constant density, no gravity:
∇ · σ = 0
Degree of freedom: displacement u. Note that σ = σ(u) (more on this later).
Infinitesimal deformation: Lagrangian/Eulerian classification does not make sense.

Index notation
Spacial axis x, y, z, . . . are uniformly represented as xα, xβ, xγ, ...
Comma "," means spatial derivatives. For example, σαβ,γ = ∂σαβ/∂xγ
Vector notation v.s. index notation:
ρ Dv/Dt = ∇ · σ + g <=> ρ Dvα/Dt = ∑ σaβ,β + gα
(σαβ,β stands for ∂σαβ/∂xβ .)
(In this lecture we do not use implicit summation for clarity.)

Discretize Cauchy momentum equation using FEM
More difficult compared to Poisson's problem: a) scalar v.s. vector; b) direct v.s. extra linear mapping.
Weak form with test function: w(x) : R2 → R2:
∑ σaβ,β wα = 0,
Integration by parts:
∑ σaβ,β wα + ∑ σaβ wα,β - ∑ (σaβ wα),β => ∑ σaβ wα,β - ∑ (σaβ wα),β.
Divergence theorem:
∀α∀w, ∫∫ ∑ σaβ wα,β dA = ∮ ∑ (σaβ wα) dnβ.

Discretization
∀α∀w, ∫∫ ∑ σaβ wα,β dA = ∮ ∑ (σaβ wα) dnβ.
Replace w and u with their discrete versions:
wα(x) = ∑ wiα φiα(x), uα(x) = ∑ ujα φjα(x)
∀α∀i, ∫∫ ∑ [σ(u(x))]αβ φiα(x) dA = ∮ ∑ (σaβ φiα) dnβ.

Relating σ to u
From infinitesimal strain theory:
Strain tensor:
e = 1/2 (∇u + (∇u)T)
Cauchy stress tensor:
σ = λtr(e)I + 2µe
... or in index notation:
eαβ = 1/2 (uα,β + uβ,α)
σαβ = λδαβ ∑ eαα + 2µeαβ
δαβ = { 1, if α = β
{ 0, if α != β
In one word: here σ is a linear function of u.

Building the linear system
σ is a linear function of u.
∀α∀i, ∫∫ ∑ [σ(u(x))]αβ φiα(x) dA = ∮ ∑ (σaβ φiα) dnβ.
Again,
Ku = f
Stencil size: away from the boundary, how many non-zero entries are there per row in the K matrix? 3^2 × 2 = 18 in 2D; 3^3 × 3 = 81 in 3D.
What does K look like?

The 8 × 8 Ke matrix
k=[ 1/2-nu/6 1/8+nu/8 -1/4-nu/12 -1/8+3*nu/8 ...
-1/4+nu/12 -1/8-nu/8 nu/6 1/8-3*nu/8];
[ k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8)
k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3)
k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2)
k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5)
k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4)
k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7)
k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6)
k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1)];
After some transforms (e.g. strain-displacement matrix B, stress-train matrix E).
(Source: A 99 line topology optimization code written in Matlab¹)
¹O. Sigmund (2001). "A 99 line topology optimization code written in Matlab”. In: Structural and multidisciplinary optimization 21.2, pp. 120–127.

Topology optimization
The minimal compliance topology optimization problem can be formulated as:
min L(ρ) = uTK(ρ)u
s.t. K(ρ)u = f
∑ ρe <= cV,
ρe ∈ [ρmin, 1]
L: measure of deformation energy, or the loss function
c: volume fraction (e.g., 0.3)
ρe: material occupancy (0 = empty, 1 = filled) of cell e.
V: total volume

Topology optimization (Demo)
Check out the supplementary material for more details. Keywords: Solid Isotropic Material with Penalization (SIMP), Optimility Criterion (OC)

Narrow-band TopOpt on a sparsely populated grid
Optimizing 1,040, 875,347 FEM voxels². [Bilibili]
²H. Liu et al. (2018). "Narrow-band topology optimization on a sparsely populated grid". In: ACM Transactions on Graphics (TOG) 37.6, pp. 1–14.

