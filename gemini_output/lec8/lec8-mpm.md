Okay, I have transcribed the text from lec8-mpm.pdf. Due to the limitations of the tool, I had to rely on OCR, which might not be perfect, but I have done my best to capture all the information. Here is the transcribed text:

Material Point Methods: A Hands-on Tutorial
GAMES 201 Lecture 8
Yuanming Hu
MIT CSAIL
July 27, 2020

Table of Contents
1 Overview
2 Moving Least Squares MPM
3 Constitutive models in MPM
4 Lagrangian forces in MPM

A little bit of MPM theory (in graphics)
Just like FEM, MPM belongs to the family of Galerkin methods. There are no
elements in MPM, so MPM ∈ Element-free Galerkin (EFG).
• MPM particles correspond to FEM quadrature points, instead of elements.
MPM typically uses one-point quadrature rule.
• MPM equations are derived using weak formulation.

Moving Least Squares MPM (MLS-MPM)
TL; DR: use MLS shape function in MPM.
• Originally proposed in SIGGRAPH 2018¹.
• Further improved in the SIGGRAPH Asia 2019 Taichi paper² to save memory
bandwidth.
• Faster and easier to implement than classical B-spline MPM.
• Reason for simplicity and performance: MPM almost always uses the APIC
transfer scheme, and MLS-MPM reuses APIC as much as possible.
1Y. Hu et al. (2018). "A moving least squares material point method with displacement
discontinuity and two-way rigid body coupling". In: ACM Transactions on Graphics (TOG) 37.4,
pp. 1-14.
2Y. Hu et al. (2019). "Taichi: a language for high-performance computation on spatially
sparse data structures". In: ACM Transactions on Graphics (TOG) 38.6, pp. 1–16.

Notations
In this lecture,
• Scalars are non-bold. E.g., mi and Vo.
p
• Vectors/matrices are bold lower-/upper-case letters respectively. E.g., vp and
Cp.
• Subscript i for grid nodes; p for particles. E.g., vi and vp.
• Superscripts are for time steps, e.g. x and xn+1

Recap: Affine Particle-in-Cell³ for incompressible fluids
1 Particle to grid (P2G)
• (mv)n+1 = ∑p Wip [Mpvn + mp Cm (xi - x)] (Grid momentum)
• mi
n+1
= ∑p mp Wip (Grid mass)
2 Grid operations
• vn+1
= (mv)n+1/mn+1 (Grid velocity)
i
• Apply Chorin-style pressure projection: vn+1 = Project(vn+1)
3 Grid to particle (G2P)
• vn+1 = ∑ Wipvn+1 (Particle velocity)
p
• Cn+1
p
•
2
4
=
2
Σi i Wipvn+1 (xi - x)T (Particle velocity gradient)
n
xn+1 = x + tv +1 (Particle position)
p
3C. Jiang, C. Schroeder, and J. Teran (2017). "An angular momentum conserving
affine-particle-in-cell method". In: Journal of Computational Physics 338, pp. 137-164.

☑ (Explicit) Moving Least Squares MPM (MLS-MPM)
1 Particle to grid (P2G)
• Fn+1 = (I+AC)F,... (Deformation update)
p
•
n
(mv) n+1 = p Wip{mpv + [mp Cm - VP (Fn+1)(n+1) T] (xi - x)}
n
p
X
p
p
(Grid momentum)
• mn+1
= ∑p mp Wip (Grid mass)
2 Grid operations
n+1
= (mv)n+1/m7+1 (Grid velocity)
•
n+1
= BC(+1) (Grid boundary condition. BC is the boundary condition
Vi
operator.)
③ Grid to particle (G2P)
• v vn+1 = 2; Wipvn+1 (Particle velocity)
V
Cn+1 = 2 Wipvn+1 (xi - x) T (Particle velocity gradient)
•
p
• xn+1
p
V
4
= X
p
x + Atvn+1 (Particle position)
Note that in classical B-spline MPM, deformation update usually happens after
G2P.

Deformation update
Deformation gradients evolve because v =
dvn
дх
x=x≠0.
(Local velocity field is not constant, so the material keeps deforming.)
Evaluating new deformation gradients:
p
Fn+1 = (I + Atv) Fn.
p
In MLS-MPM, APIC CD is used as an approximation of V.
Therefore in MLS-MPM we have
Fn+1 = (I + AtCn)Fn.

P2G: Computing internal forces
Recall that
n
(mv) = p Wip{mpv + [mpCpVP(Fn+1) (Fn+1)T] (x – x)} (Grid
i
n
p
2
momentum).
Two momentum terms:
n
• APIC: Wip [mpv + mpCn(xi-x)]
• Particle elastic force (impulse):
Atfip = - Wip VOP(Fn+1) (Fn+1) T] (xi - x)
4At
p
Assuming hyperelastic materials. Deriving f; using potential energy gradients:
U = Σνρ (Fp)
p
JU
fi
=
дхі
Ψp: elastic energy density of particle p; U: total elastic potential energy.
Vo: particle initial volume.

P2G: Computing nodal force fi
Assume we move forward τ → 0, and then compute deformed grid node location
Xi = xi + TVi, Cp = ∑i Wipi(xi - xp) T, updated F = (I + TCp)Fp:
2
4
2
X
году (F)
VO
p
V
JU
fi
=
=
дхі
p
=
Σ
P
V dp (Fp)
τ
dvi
p
=
Σ
p
τ
p
VO
=
di
p
V (F) ӘЕ ӘСр
ΣP(F) TF2
PPP(FD) TFT τF
τ
p
4
Δx2
X
•
p
T
VP(F) F Wip(Xi - Xp)
p

Grid operations: enforcing boundary conditions (BC)
BC in MPM should be applied on the grid. For all grid nodes i within the
boundary:
vn+1 =
BCsticky (n+1) = 0
vn+1
.n+1
=
BCslip (n+1)
= vn+1 = n(nTvn+1)
vn+1 = BCseparate (n+1)
i
= vn+1 – n. min(n+1,0)
(n: surface normal)
Extras:
• Adding gravity vn+1 + = Atg
i
② Moving collision object
3 Coulomb Friction

Summary: benefits of MLS-MPM
Why is MLS-MPM (SIGGRAPH 2018) easier and faster than classical B-spline
MPM (SIGGRAPH 2013)?
① Directly reuse APIC Cp as an approximation of Vv for deformation gradient
update. No need to evaluate ∇ Wip (Fewer FLOPs)
② Easy to move deformation update from G2P to P2G, because we only need
Cp for deformation update. (Fewer bytes to fetch from main memory)
● In P2G, APIC and MLS-MPM momentum contribution can be fused, since
they are both "MLS”. (Fewer FLOPs)
MLS-MPM is consistent with the weak formulation of the Cauchy momemtum
equation. See the original MLS-MPM paper for a correctness proof.
4Y. Hu et al. (2018). "A moving least squares material point method with displacement
discontinuity and two-way rigid body coupling". In: ACM Transactions on Graphics (TOG) 37.4,
pp. 1-14.

Constitutive Models
Common constitutive models in MPM:
1 Elastic objects: NeoHookean & Corotated
② Fluid: Equation-of-States (EOS)
• Elastoplastic objects (snow, sand etc.): Yield criteria: ad-hoc boxing5,
Cam-clay6, Drucker-prager7, NACC,
...
Two critical aspects of a constitutive model in MPM:
① (Elastic/plastic) deformation update
② (PK1) stress evaluation
5A. Stomakhin et al. (2013). “A material point method for snow simulation". In: ACM
Transactions on Graphics (TOG) 32.4, pp. 1–10.
6J. Gaume et al. (2018). "Dynamic anticrack propagation in snow". In: Nature
communications 9.1, pp. 1-10.
7G. Klár et al. (2016). "Drucker-prager elastoplasticity for sand animation". In: ACM
Transactions on Graphics (TOG) 35.4, pp. 1–12.

Constitutive models for elastic solids
Deformation update: simply Fn+1 = (I + AtC)Fn.
PK1 stresses of hyperelastic material models:
• Neo-Hookean:
•
•
p
2
(F) = [(FTF) ii - 1] – µ log(J) + log²(J).
μ
2
дψ
-
P(F) = = μ(F – F-T) + λ log(J)F-T
• (Fixed) Corotated:
2
• ψ(F) = μ Σ₁(σ₁ − 1)² + (J − 1)². σ₁ are singular values of F.
• P(F) = = 2μ(F – R) + λ(J − 1).JF-T
ду
Cauchy stress σ =
PFT is usually unused in MPM.
More details: check out the SIGGRAPH 2016 MPM course.
°C. Jiang et al. (2016). "The material point method for simulating continuum materials". In:
ACM SIGGRAPH 2016 Courses, pp. 1–52.

Constitutive models weakly compressible fluids⁹
Volume ratio Jp = V/V = det (Fr).
The simplest equation of state: p = K(1 – J), Cauchy stress σ = −pI. K: bulk
modulus.
Computing det (Fm) can be numerically unstable.
11
Recall that for F2x2, det (F) = F00F11 - F01F10. The "-" opeartion leads to
catastrophic cancellation. Same for F3×3 (Question: why doesn't this happen
to NeoHookean/corotated materials?)
Deformation update: instead of maintaining Fp, directly maintain J = det (Fm):
Fn+1 = (I+AtCn)Fn
⇒ det(Fn+1) = det (I+AtCp) det (Fm)
⇒ J+1 = (1 + Attr(C)) J
p
p
A. P. Tampubolon et al. (2017). "Multi-species simulation of porous sand and water
mixtures". In: ACM Transactions on Graphics (TOG) 36.4, pp. 1–11.

Simulating weakly compressible fluids (lazy solution)
Setting u to zero in (Fixed) corotated model.
(Recap) In corotated materials:
2
• ψ(F) = μ Σί(σε – 1)² + (J − 1)2. σ₁ are singular values of F.
•
дψ
ƏF
P(F) = = 2μ(F – R) + 2(J − 1)JF-T

Recap: Singular value decomposition (SVD)
Theorem
(Existence of singular value decompositions) Every real matrix Mnxm can be
decomposed into Mnxm = UnxnnxmV where U and V are orthonormal
matrices, and ∑ is a diagonal matrix.
mxm'
Diagonal entries σ₁ = ∑ii are called singular values.
To learn more about linear algebra: check out Gilbert Strang's MIT OCW.

SVD: Intuition

2 × 2 and 3 × 3 SVD10 in Taichi
Note that SVD is not unique. We additionally require
• det(U) = det(V) = 1.
• Σii are sorted in decreasing order.
• Only the singular value with smallest magnitude can be negative.
Example
U, sig, V
=
ti.svd(M) %2 sig is an NxN diagonal matrix.
10A. McAdams et al. (2011). Computing the singular value decomposition of 3x3 matrices with
minimal branching and elementary floating point operations. Tech. rep. University of
Wisconsin-Madison Department of Computer Sciences.

Simulating elastoplastic solids
In hyperelastic settings:
Fp
=
n
F p,elastic Fp,plastic, p = (Fp,elastic),
i.e., the potential energy penalizes elastic deformation only.
Example
"Box" yield criterion11: deformation udpate:
• Evolve Fn+1 = (I+AC) Fr, elastic
● SVD: Fn+1 = UÊVT
p
● Clamping: Eii = max(min(∑ii, 1 + 0s), 1 – 0c) (forget about too large
deformations)
4 Reconstruct: Fn+1 = UEVT; move clamped parts to Fn+1
p,elastic
p,plastic
11A. Stomakhin et al. (2013). "A material point method for snow simulation". In: ACM
Transactions on Graphics (TOG) 32.4, pp. 1–10.

Lagrangian forces in MPM12
TL; DR: Treat MPM particles as FEM vertices, and use FEM potential energy
model. A triangular mesh is needed.
Benefits:
• (Compared to FEM): Self-collision is handled on the grid;
• (Compared to MPM): Numerical fracture is avoided due to the mesh
connectivity.
• Can easily couple MPM and FEM.
Easy to implement in Taichi using AutoDiff: ti example mpm_lagrangian_forces
12C. Jiang et al. (2015). "The affine particle-in-cell method". In: ACM Transactions on
Graphics (TOG) 34.4, pp. 1–10.

Introducing Taichi “field"
Upgrading Taichi: pip install --upgrade taichi==0.6.22
Use "field" instead of "tensor" since Taichi v0.6.22
• The name "tensor” is deprecated. Always use "field" instead.
• ti.var is deprecated. Use ti.field instead.
• Argument dt is deprecated. Use dtype instead.
Declaring fields in Taichi
# particle_x
=
ti.Vector(3, dt=ti.f32, shape=1024)
particle_x
=
ti.Vector.field (3, dtype=ti.f32, shape=1024)
=
# density
density
=
=
num_springs ti.field(dtype=ti.i32, shape=())
particle_F = ti. Matrix.field (3, 3, dtype=ti.f32, shape=1024)
ti.field(dtype=ti.f32, shape=(256, 256))
ti.var(dtype=ti.f32, shape=(256, 256))

Fields := global variables in Taichi
Distinguishing global fields from local variables
Global variables are always declared with "field". Local variables are always
declared without "field":
X =
ti.Vector.field (3, dtype=ti.f32, shape=(128, 512)) # global
@ti.kernel
def foo():
a
=
ti.Vector([0.2, 0.4]) # local
The word "field" refers to
...
1 a component of a (database) record. For example, mass and volume properties
of a particle array.
② a (physical) quantity that is assigned to every point in space. E.g., "velocity
fields", and "magnetic fields". High-dimensional arrays of
scalars/vectors/matrices are exactly "fields" sampled at discrete grid points.

