Okay, I have transcribed the text from the PDF. Due to the limitations of the tool, I wasn't able to directly extract the text. Instead, I processed the binary content and used OCR to extract the text. Here's the transcription:

Page 1:
Basics of deformation, elasticity, and finite elements
Yuanming Hu
Deformation
Elasticity
FEM basics
Basics of deformation, elasticity, and finite elements
Yuanming Hu
MIT CSAIL
June 15, 2020

Page 2:
:
Basics of
deformation,
elasticity, and
finite elements
Yuanming Hu
Deformation
Elasticity
FEM basics
Overview
Simulating elastic materials is a lot of fun!
• Cool visual effects
• Not too hard to implement (using Taichi (Demos))
• Base of other materials (viscoelastic, elastoplastic, viscoplastic...)
Recommended reading
① The classical FEM method and discretization methodology by Eftychios
Sifakis¹
② The Material Point Method for Simulating Continuum Materials by
Chenfanfu Jiang et al.2
1E. Sifakis and J. Barbic (2012). "FEM simulation of 3D deformable solids: a practitioner's
guide to theory, discretization and model reduction". In: Acm siggraph 2012 courses, pp. 1–50.
2C. Jiang et al. (2016). "The material point method for simulating continuum materials". In:
ACM SIGGRAPH 2016 Courses, pp. 1–52.

Page 3:
Basics of
deformation,
elasticity, and
finite elements
Yuanming Hu
Deformation
Elasticity
FEM basics
Table of Contents
1 Deformation
2 Elasticity
3 FEM basics

Page 4:
Basics of
deformation,
elasticity, and
finite elements
Yuanming Hu
Deformation
Deformation map : a (vector to vector) function that relates rest material
position and deformed material position.
Deformation
Elasticity
Xdeformed = (Xrest)
FEM basics
Deformation gradient F
F :=
dxdeformed
dxrest
Deformation gradients are translational invariant
Φ1
=
$(Xrest) and $2 = (Xrest) + c have the same deformation gradients!
Deform/rest volume ratio J = det(F)

Page 5:
Basics of
deformation,
elasticity, and
finite elements
Yuanming Hu
Deformation
Elasticity
FEM basics
Table of Contents
1 Deformation
2 Elasticity
3 FEM basics

Page 6:
Hyperelasticity
Basics of
deformation,
elasticity, and
finite elements
Yuanming Hu
Deformation
Elasticity
FEM basics
Hyperelastic materials: materials whose stress-strain relationship is defined by a
strain energy density function
ψ = y(F)
Intuitive understanding: y is a potential function that penalizes deformation.
"Stress": the material's internal elastic forces.
"Strain”: just replace it with deformation gradients F for now.
Be careful
We use y as the strain energy density function and & as the deformation map.
They are completely different.

Page 7:
Stress tensor
Basics of
deformation,
elasticity, and
finite elements
Yuanming Hu
Deformation
Elasticity
FEM basics
Stress stands for internal forces that infinitesimal material components exert on
their neighborhood.
Based on our need, we use different measures of stress
• The First Piola-Kirchhoff stress tensor (PK1): P(F) =
compute, but in rest space)
• Kirchhoff stress: τ
ƏF
(F) (easy to
• Cauchy stress tensor: σ (symmetric, because of conservation of angular
momentum)
Relationship: τ
=
PFT P Jo =
=
JoF-T Traction t = στη.
Intuition of P = JoF-T: F-T compensates for material deformation. (Note
that it's F-T instead of F-1 since we transform the normal n instead of x.)

Page 8:
Basics of
deformation,
elasticity, and
finite elements
Elastic moduli (isotropic materials)
•
Yuanming Hu
Young's modulus E =
ε
• Bulk modulus K =
Deformation
Elasticity
FEM basics
VdP
dV
• Poisson's ratio v ∈ [0.0,0.5) (Auxetics have negative Poisson's ratio)
Lamé parameters:
• Lamé's first parameter μ
• Lamé's second parameter λ (aka. shear modulus, denoted by G)
Useful conversion formula:
K =
E
3(1 – 2v)
Ev
E
λ =
μ
=
(1 + v)(1 – 2v)
2(1 + v)

Page 9:
Hyperelastic material models
Basics of
deformation,
elasticity, and
finite elements
Yuanming Hu
Deformation
Elasticity
FEM basics
Popular ones in graphics:
• Linear elasticity (small deformation only)
• Neo-Hookean:
2
• (F) = Σ[(FTF) ii - 1] – μ log(J) + log² (J).
2
дү
• P(F) = = μ(F – FT) + 2 log(J)F-T
• (Fixed) Corotated:
• ψ(F) = μ Σ₁(σε − 1)² + (J − 1)2. σ₁ are singular values of F.
ду
• P(F) = = 2μ(F – R) + 2(J − 1).JF-T
JF
More details: The Material Point Method for Simulating Continuum Materials3
3C. Jiang et al. (2016). "The material point method for simulating continuum materials". In:
ACM SIGGRAPH 2016 Courses, pp. 1–52.

Page 10:
Basics of
deformation,
elasticity, and
finite elements
Yuanming Hu
Deformation
Elasticity
FEM basics
Table of Contents
1 Deformation
2 Elasticity
3 FEM basics

Page 11:
The finite element method
Basics of
deformation,
elasticity, and
finite elements
Yuanming Hu
Finite element method: Galerkin discretization scheme that builds discrete
equations using weak formulations of continuous PDEs. (More details later in this
course.)
Deformation
Elasticity
FEM basics
1
0.5
0
1
1
0
-1:1
Figure: A solution to a discretized partial differential equation, obtained with FEM
(source: Wikipedia)

Page 12:
Linear tetrahedral (triangular) FEM
Basics of
deformation,
elasticity, and
finite elements
Yuanming Hu
Deformation
Elasticity
FEM basics
Linear tetrahedral finite elements (for elasticity) assume the deformation map ф
is affine and thereby deformation gradient F is constant within a single
tetrahedral element:
Xdeformed = Fxrest + p.
For every element e, its elastic potential energy
U(e)
=
y(F(x))x = Vey(Fe).
e
Question: how to compute Fe(x)?

Page 13:
:
Basics of
deformation,
elasticity, and
finite elements
Yuanming Hu
Deformation
Elasticity
FEM basics
Computing Fe in linear triangular finite elements (1)
Recall that
Xdeformed = FXrest + p.
In 2D triangular elements (3D would be tetrahedral elements), assuming the rest
positions of the vertices are arest, brest, Crest and deformed positions are
adeformed, bdeformed, Cdeformed. Since within an linear triangular element F is
constant, we have
adeformed
Farest + p
(1)
bdeformed
Fbrest + p
(2)
Cdeformed
FCrest + p
(3)
Let's eliminate p:
(adeformed - Cdeformed)
=
F(arest - Crest)
(4)
(bdeformed - Cdeformed)
=
F(brest - Crest)
(5)

Page 14:
:
Basics of
deformation,
elasticity, and
finite elements
Yuanming Hu
Deformation
Elasticity
FEM basics
Computing Fe in linear triangular finite elements (2)
(adeformed - Cdeformed)
=
F(arest - Crest)
(6)
(bdeformed - Cdeformed
=
F(brest - Crest)
(7)
Note that F2×2 now has four linear constraints (equations).
B
=
[arest - Crest Brest - Crest] -1
(8)
D
=
F
adeformed
DB
Cdeformed | bdeformed
Cdeformed
(9)
(10)
(B is constant through out the physical process. Therefore it should be
pre-computed.)

Page 15:
Basics of
deformation,
elasticity, and
finite elements
Yuanming Hu
Deformation
Elasticity
FEM basics
Explicit linear triangular FEM simulation
Recall the Semi-implicit Euler (aka. symplectic Euler) time integration
Vt+1,i
=
=
ft,i
Vt,i + At
Mi
Xt+1,i xt,i + Atvt+1,i
Note that xt,i and vt,i are stored on the vertices of finite elements
(triangles/tetrahedrons).
ft,i
=
JU
дхі
=
Συ(e)
=
дхі
e
ΣVe Σν
-
e
дψ(Fe) Fe
dFе дхі
e
=
-
Fe
ΣVeP(Fe) dx
Don't want to compute P(Fe)? Use Taichi's AutoDiff system.
e

Let me know if you have any other questions.

