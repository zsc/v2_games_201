Okay, I have processed the file and extracted the text. Due to the presence of mathematical equations and diagrams, the transcription might not be perfect, but I've done my best to capture the content accurately.

Here is the transcribed text:

```
Lagrangian
Simulation
Approaches
Yuanming Hu
Mass-spring
systems
Time integration
Lagrangian fluid
simulation:
Smoothed
particle
hydrodynamics
Lagrangian Simulation Approaches
Mass-Spring Systems and Smoothed Particle Hydrodynamics
Yuanming Hu
MIT CSAIL
June 8, 2020
```

```
Lagrangian
Simulation
Approaches
Yuanming Hu
Mass-spring
systems
Time integration
Lagrangian fluid
simulation:
Smoothed
particle
Table of Contents
1 Mass-spring systems
2 Time integration
hydrodynamics
3 Lagrangian fluid simulation: Smoothed particle hydrodynamics
```

```
Lagrangian
Simulation
Approaches
Yuanming Hu
Mass-spring Systems
Mass-spring
systems
Time integration
Lagrangian fluid
simulation:
Demo!
Smoothed
particle
hydrodynamics
```

```
Mass-spring
systems
Time integration
Lagrangian
Simulation
Approaches
Yuanming Hu
Mass-spring systems
fij
f₁ fi = fij
=
-k(||xi - xj||2 - lij) (xi - xj) (Hooke's Law)
j≠i
Lagrangian fluid
j
simulation:
Smoothed
particle
hydrodynamics
dvi
dt
дхі
dt
1
=
f; (Newton's second law of motion)
Mi
=
Vi
k: spring stiffness; lij: spring rest length between particle i and particle j;
mi: the mass of particle i. (x – xj): direction vector from particle i to particle i;
means normalization.
```

```
Lagrangian
Simulation
Approaches
Yuanming Hu
Mass-spring
systems
Time integration
Lagrangian fluid
simulation:
Smoothed
particle
Table of Contents
1 Mass-spring systems
2 Time integration
hydrodynamics
3 Lagrangian fluid simulation: Smoothed particle hydrodynamics
```

```
Lagrangian
Simulation
Time integration
Approaches
1 Forward Euler (explicit)
Yuanming Hu
Mass-spring
Vt+1
=
systems
Time integration
Xt+1
=
Lagrangian fluid
simulation:
Smoothed
particle
f+
t
Vt + At
m
xt + Atvt
hydrodynamics
② Semi-implicit Euler (aka. symplectic Euler, explicit)
=
Vt+1
Xt+1
=
ft
Vt + At
m
xt + ∆tvt+1
3 Backward Euler (often with Newton's method, implicit)
```

```
Lagrangian
Simulation
Approaches
Yuanming Hu
Mass-spring
systems
Implementing a mass-spring system with symplectic Euler
Time integration
Lagrangian fluid
1 Vt
tft
• Compute new velocity using Vt+1 = vt + th
m
simulation:
Smoothed
② Collision with ground
particle
hydrodynamics
3 Compute new position using Xt+1 = Xt + tvt+1
```

```
Lagrangian
Simulation
Approaches
Yuanming Hu
Mass-spring
systems
Time integration
Lagrangian fluid
simulation:
Smoothed
particle
hydrodynamics
Implementing a mass-spring system with symplectic Euler
Showcase mass_spring.py
@ti.kernel
def substep():
n = num_particles [None]
# Compute force and new velocity
for i in range(n):
v[i] *= ti.exp(-dt * damping [None]) # damping
total_force = ti. Vector (gravity) * particle_mass
for j in range(n):
if rest_length[i, j] != 0:
x_ij = x[i]
x[j]
total_force += -spring_stiffness [None] * (x_ij. norm() rest_length[i, j]) * x_ij.
normalized()
v[i] += dt * total_force / particle_mass
# Collide with ground
for i in range(n):
if x[i].y < bottom_y:
x[i].y = bottom_y
v[i].y = 0
# Compute new position
for i in range(n):
x[i] += v[i] * dt
```

```
Lagrangian
Simulation
Approaches
Yuanming Hu
Mass-spring
systems
Time integration
Lagrangian fluid
simulation:
Smoothed
particle
hydrodynamics
Explicit v.s. implicit time integrators
Explicit (forward Euler, symplectic Euler, RK, ...):
• Future depends only on past
• Easy to implement
• Easy to explode:
• Bad for stiff materials
-
At≤cm (c~1)
Implicit (backward Euler, middle-point, ...):
• Future depends on both future and past
• Chicken-egg problem: need to solve a system of (linear) equations
• In general harder to implement
• Each step is more expensive but time steps are larger
• Sometimes brings you benefits
...
•
but sometimes not
• Numerical damping and locking
```

```
Lagrangian
Simulation
Approaches
Mass-spring systems
Implicit time integration:
Yuanming Hu
Mass-spring
Xt+1
=
xt + ∆tvt+1
systems
=
Vt+1 vt+ ∆tM¯¹f(xt+1)
(1)
(2)
Time integration
Lagrangian fluid
simulation:
Eliminate Xt+1:
Smoothed
particle
hydrodynamics
Vt+1
vt + AtM¯¹f(xt + Atvt+1)
(3)
Linearize (one step of Newton's method):
Vt+1
= vt + AtM-1 f(xt) + (xt)tvt+1
[f(x1) + x(x) tvt+1]
(4)
```

```
Lagrangian
Simulation
Approaches
Yuanming Hu
After linearization
Linearize:
df
Vt+AtM-1 f(x) +x (xt)∆tvt+1
Mass-spring
systems
Time integration
Vt+1
=
vt
Lagrangian fluid
simulation:
Smoothed
particle
hydrodynamics
Clean up:
2
Ι – ΔΜ-1
A nice linear system!
1 df
dx(xt) Vt+1
(5)
=
vt + AtM¯¹f(xt)
(6)
```

```
Lagrangian
Simulation
Approaches
Yuanming Hu
Solving the system
Mass-spring
systems
Time integration
Lagrangian fluid
simulation:
Smoothed
particle
hydrodynamics
How to solve it?
2
Ι – Δ²M-1
df
дх
=
(xt) Vt+1 vt + AtM¯¹f(xt)
1
(7)
• Jacobi/Gauss-Seidel iterations (easy to implement!)
• Conjugate gradients (later in this course)
```

```
Lagrangian
Simulation
Approaches
Yuanming Hu
Solving the system
Mass-spring
systems
Time integration
Lagrangian fluid
simulation:
A =
Smoothed
particle
hydrodynamics
b
=
Avt+1
=
b
2
-1 дf
I – AM-1(xt)
дх
1
vt + ∆tM¯¹f(xt)
```

```
Lagrangian
Simulation
Approaches
Yuanming Hu
Solving linear systems with Jacobi iterations (Demo!)
Mass-spring
systems
Time integration
Lagrangian fluid
simulation:
Smoothed
particle
hydrodynamics
A
=
X
=
ti.var(dt=ti.f32, shape=(n, n))
ti.var(dt=ti.f32, shape=n)
=
new_x ti.var(dt=ti.f32, shape=n)
b = ti.var(dt=ti.f32, shape=n)
@ti.kernel
def iterate():
for i in range(n):
r =
b[i]
for j in range(n):
if i != j:
r
=-
A[i, j] * x[j]
new_x[i] = r / A[i, i]
for i in range(n):
x[i]
=
new_x[i]
```

```
Lagrangian
Simulation
Approaches
Yuanming Hu
Mass-spring
systems
Time integration
Lagrangian fluid
simulation:
Smoothed
particle
hydrodynamics
Unifying explicit and implicit integrators
df
[1-1 I-βAM-1(xt) (xt) Vt+1
=
vt + ∆tM¯¹f(xt)
1 β = 0: forward/semi-implicit Euler (explicit)
2 β = 1/2: middle-point (implicit)
③ β = 1: backward Euler (implicit)
```

```
Lagrangian
Simulation
Approaches
Yuanming Hu
Mass-spring
systems
Time integration
Lagrangian fluid
simulation:
Smoothed
particle
hydrodynamics
Solve faster
What if we have millions of mass points and springs?
• Sparse matrices
• Conjugate gradients
• Preconditioning
• Use position-based dynamics¹
A different (yet much faster) approach: Fast mass-spring system solver2
1M. Müller et al. (2007). "Position based dynamics". In: Journal of Visual Communication
and Image Representation 18.2, pp. 109–118.
2T. Liu et al. (2013). "Fast simulation of mass-spring systems". In: ACM Transactions on
Graphics (TOG) 32.6, pp. 1–7.
```

```
Lagrangian
Simulation
Approaches
Yuanming Hu
Mass-spring
systems
Time integration
Lagrangian fluid
simulation:
Smoothed
particle
Table of Contents
1 Mass-spring systems
2 Time integration
hydrodynamics
3 Lagrangian fluid simulation: Smoothed particle hydrodynamics
```

```
Lagrangian
Simulation
Approaches
Yuanming Hu
Mass-spring
systems
Time integration
Lagrangian fluid
simulation:
Smoothed
particle
hydrodynamics
Smoothed particle hydrodynamics
High-level idea: use particles carrying samples of physical quantities, and a kernel
function W, to approximate continuous fields: (A can be almost any spatially
varying physical attributes: density, pressure, etc. Derivatives: different story)
A(x) = ∑ A W(||x − xj||2, h)
i
sh
Mi
i
Pi
X
W(|r-r, h)
i
Figure: SPH particles and their kernel (source: Wikipedia)
```

```
Lagrangian
Simulation
Approaches
Yuanming Hu
Mass-spring
systems
Time integration
Lagrangian fluid
simulation:
Smoothed
particle
Smoothed particle hydrodynamics
1 Originally proposed for astrophysical problems³
2 No meshes. Very suitable for free-surface flows!
• Easy to understand intuitively: just imagine each particle is a small parcel of
water (although strictly not the case!)
hydrodynamics
3R. A. Gingold and J. J. Monaghan (1977). "Smoothed particle hydrodynamics: theory and
application to non-spherical stars". In: Monthly notices of the royal astronomical society 181.3,
pp. 375-389.
4J. J. Monaghan (1994). "Simulating free surface flows with SPH". In: Journal of
computational physics 110.2, pp. 399–406.
```

```
Lagrangian
Simulation
Approaches
Yuanming Hu
Mass-spring
systems
Time integration
Lagrangian fluid
simulation:
Smoothed
particle
hydrodynamics
Implementing SPH using the Equation of States (EOS)
Also known as Weakly Compressible SPH (WCSPH)5.
Momentum equation: (p: density; B: bulk modulus; y: constant, usually ~ 7)
1
Dr-p+g, p= B(()-1)
Dv
Dt
Mi
ρ
▽p+g,
A(x) = ∑ A W(||x - xj||2, h), pi =
i
Pi
3=
j
Extras: surface tension, viscosity; (very) nice tutorial6
mj W(||xi - xj||2, h)
Note: the WCSPH paper should have used material derivatives.
5M. Becker and M. Teschner (2007). "Weakly compressible SPH for free surface flows". In:
Proceedings of the 2007 ACM SIGGRAPH/Eurographics symposium on Computer animation.
Eurographics Association, pp. 209–217.
6
'D. Koschier et al. (2019). "Smoothed Particle Hydrodynamics Techniques for the Physics
Based Simulation of Fluids and Solids". In:
```

```
Lagrangian
Simulation
Approaches
Yuanming Hu
Gradients in SPH
A(x) = Σ Αί
i
mi W(||x - xj||2,
Pi
Ai Aj
2
h)
+ 'x₁ W(||xi - xj||2, h)
P ρρ P
2
Mass-spring
systems
Time integration
Lagrangian fluid
simulation:
Smoothed
Ai = ρεΣ mj
particle
j
hydrodynamics
• Not really accurate...
• but at least symmetric and momentum conserving!
Now we can compute Vpi.
Extension: Laplace operator (viscosity etc.)...
```

```
Lagrangian
Simulation
Approaches
SPH Simulation Cycle
Recall:
Yuanming Hu
Mass-spring
systems
Dv
Dt
1
ρ
p+g,p=B(()-1)
Time integration
Lagrangian fluid
simulation:
Smoothed
particle
1 For each particle i, compute pi = ∑j mj W(||xi - xj||2, h)
② For each particle i, compute Vpi using the gradient operator
Xi
hydrodynamics
• Symplectic Euler step (again...):
Dv
Vt+1 = vt + At-
Dt
Xt+1 = xt + ∆tvt+1
```

```
Variants of SPH
Recent updates:
...
• Predictive-Corrective Incompressible SPH (PCI-SPH)7
Lagrangian
Simulation
Approaches
Yuanming Hu
•
Mass-spring
systems
Time integration
Lagrangian fluid
simulation:
Smoothed
particle
hydrodynamics
• Position-based fluids (PBF) Demo: ti example pbf2d
• Divergence-free SPH (DFSPHº)
•
...
Survey paper: SPH Fluids in Computer Graphics10
7
B. Solenthaler and R. Pajarola (2009). "Predictive-corrective incompressible SPH". In: ACM
SIGGRAPH 2009 papers, pp. 1–6.
8M. Macklin and M. Müller (2013). "Position based fluids". In: ACM Transactions on
Graphics (TOG) 32.4, pp. 1–12.
9J. Bender and D. Koschier (2015). "Divergence-free smoothed particle hydrodynamics". In:
Proceedings of the 14th ACM SIGGRAPH/Eurographics symposium on computer animation,
pp. 147-155.
10M. Ihmsen et al. (2014). "SPH fluids in computer graphics". In:
```

```
Lagrangian
Simulation
Approaches
Yuanming Hu
Mass-spring
systems
Time integration
Lagrangian fluid
simulation:
Smoothed
particle
hydrodynamics
Courant-Friedrichs-Lewy (CFL) condition
One upper bound of time step size:
C=
uAt
Δχ
< Cmax ~ 1
• C: CFL number (Courant number, or simple the CFL)
• t: time step
•
•
x: length interval (e.g. particle radius and grid size)
น: maximum (velocity)
Application: estimating allowed time step in (explicit) time integrations.
Typical Cmax in graphics:
• SPH: ~ 0.4
• MPM: 0.3 ~ 1
• FLIP fluid (smoke): 1 ~ 5+
```

```
Lagrangian
Simulation
Approaches
Yuanming Hu
Mass-spring
systems
Time integration
Lagrangian fluid
simulation:
Smoothed
particle
hydrodynamics
Accelerating SPH: Neighborhood search
So far, per substep complexity of SPH is O(n²). This is too costly to be practical.
In practice, people build spatial data structure such as voxel grids to accelerate
neighborhood search. This reduces time complexity to O(n).
h
h
Figure: Neighborhood search with hashing. Source: Koschier et al. 2019.
Reference: Compact hashing
```

```
Lagrangian
Simulation
Approaches
Yuanming Hu
Mass-spring
systems
Time integration
Lagrangian fluid
Other particle-based simulation methods
• Discrete element method, e.g.
11
• Moving Particle Semi-implicit (MPS)12
• Power Particles: An incompressible fluid solver based on power diagrams13
• A peridynamic perspective on spring-mass fracture14
simulation:
Smoothed
particle
hydrodynamics
...
11N. Bell, Y. Yu, and P. J. Mucha (2005). "Particle-based simulation of granular materials”.
In: Proceedings of the 2005 ACM SIGGRAPH/Eurographics symposium on Computer animation,
pp. 77-86.
12S. Koshizuka and Y. Oka (1996). "Moving-particle semi-implicit method for fragmentation
of incompressible fluid". In: Nuclear science and engineering 123.3, pp. 421–434.
13F. de Goes et al. (2015). "Power particles: an incompressible fluid solver based on power
diagrams.". In: ACM Trans. Graph. 34.4, pp. 50–1.
14J. A. Levine et al. (2014). "A peridynamic perspective on spring-mass fracture". In:
Proceedings of the ACM SIGGRAPH/Eurographics Symposium on Computer Animation.
Citeseer, pp. 47–55.
```

```
8 Exporting your results taichi v0.6.8
Lagrangian
Simulation
Approaches
Yuanming Hu
Mass-spring
systems
Time integration
Lagrangian fluid
simulation:
Smoothed
particle
hydrodynamics
Make an mp4 video out of your frames
1 Use ti.GUI.show [doc] to save the screenshots. Or simply use
ti.imwrite(img, filename) [doc].
2 ti video creates video.mp4 using frames under the current folder. To specify
frame rate, use ti video -f 24 or ti video -f 60.
3 Convert mp4 to gif and share it online: ti gif -i input.mp4.
Make sure ffmpeg works!
• Linux and OS X: with high probability you already have ffmpeg.
• Windows: install ffmpeg on your own [doc].
More information: [Documentation] Export your results.
```

Let me know if you have any other questions.

