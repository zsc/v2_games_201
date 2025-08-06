Okay, I have transcribed the text from the PDF. Due to the nature of the document (a scanned PDF), the transcription includes some OCR errors and might not be perfectly accurate. Also, some mathematical equations and symbols might not be rendered correctly in plain text.

Here's the transcription:

```
Poisson's Equation and
Fast Method
Poisson's Euqation and its Fundamental
Solution
2
$=-p
$(x-200)=0
((y)
Φ(x) = SPLYD dy
-
4711x-Y12
M₃
4 Rij
For Example: the Gravitational Problem
f(x) = ∇¢
N
f(x)=- Σ
j=1,j≠i
N
Xi
Pj v j
4π||xi - xj||2/
-
X
fi =- PjV; (xi - x₁)
j=1,j≠i 4n||xi - xj||
3
2
Given N particles and M evaluation position, direct
computation requires O(NM) time!
Introduction to Fast Summation
ZD. Complex Number. $(z) $(x)=Re(&=\log(2-23) = Re (Egglog(2-23))
Consider a source and its
potential
ο εί
Z
Apply Taylor expansion gives:
中(x)=9;log(2-Eⅰ)
P
= Pilogz - Bi
k
KkZk
Multipole Expansion:
⇒ 109(2-2))
P
=(28;)(og(z) - Σ
Q=32%.
Q(k)=-2
jk
k
Σ
3 k
P Qk
4(z)=Qloyz +/
셔水
k
NlogN Algorithm:
Compute
Q. Qk for
each level
Eval:
Tree Code
Take one step further:
If we know M-Expansion at z1 (M1={z1, Q, Q_k},
What is the M-Expansion at z2 (M2={z2, Q2, Q2_k}?
We want to obtain the coefficients from M1, not q_i's
Q
QK
Q,QK?
Z2
Pbk
$(z) = Q log(2-2) +2 + Σ
P
Recall: $(z)=Qlag z+ 2
셔
Qk
|
bk is a generalization of QR:
k
bk=-Q(z-zz) K
k
k
Recall: 2; (2;-1)*
k
十
k
=
k-ik-1
Qi (z-z) (2-1.
Rest of terms
Pbk
k
$(z)=Qlog(2-c)+/01(水
View Source as Multipole:
Q=½,
Reveal "Multipole Expansion"
$(z) = Qloq(2)-2-()
k
K
(2-c)k
From "Multipoles
こ
Compute ble with "Rest of terms"
QK
⑧ M2M Transform
struct Multipole{
};
vec2 center;
complex q[p];
//source charge is a special Multipole,
//with center = charge pos, q0 = qi, q1...q4 =0
Multipole M2M(std::vector<Multipole> &qlist)
{
Multipole res;
res.center = weightedAverageof(qlist[i].center*qlist[i].q[0]);
q[0] = sumof(qlist[i].q[0]);
for(k=1:p){
res.q[k]=0;
for(j=0:qlist.size()) {
res.q[k] += computeBk(qlist[i]);
}
}
return res;
}
Question:
If we know M-Expansion at c1 (M1={c1, Q, Q_k},
What is the polynomial at z1, so that potentials at
neighbor z can be evaluated.
Q,Qk,
P
Φ(z) = Q log (z-c) + ∑ bk
=
=(Z-K
P
=Qlog(z₁-c+z-z)+=bk
k=1
= Qlog(zi-c) + {bk
셔(라
k
(21-(+Z-8)k
M2L
Transform
Z
7
Φ(z)
+ H.O.T
P
H.O. T. = ∑ by(z-z)²
(=1
Q
P
1 bk ltk-1
|
et
k-1
struct Localpole{
};
vec2 center;
complex b[p];
Question:
If we know L-Expansion at c1 (L1={c1, B}),
What is the polynomial at c2, so that potentials at neighbor z around c2 can be evaluated.
L2L
Transform
P
P
k
Σακ (2-20) k∑ bk (2)
k=0
k=0
Honer scheme
Multipole Expansion :
Coarsening
Localpole Expansion: Interpolation
O(N) Algorithm:
M2M
Mutipole
k
=
O(N)
M2L
L2L+ M2L
P
= O(CN)
Tree code
54
2
X
1
3
6
X
2
Phi(x1) = contribution from (node1, 2, 3, 4, 5, 6)
Phi(x2) = contribution from (node1, 2, 3, 4, 5, 6)
VS
FMM
54
X
6
1
X
2
2
3
= contribution from (node1, 2,
3)
Phi(x1) = L2L frorn
+ contribution from (4,
5, 6)
Phi(x2) = L2L from
5, 6)
+ contribution from (4,
In 3D, the algorithm can be obtained via:
• Taylor expansion of the Green's function in Cartesian system
• Taylor expansion of the Green's function in Spherical
coordinates
For more details, checkout:
https://math.nyu.edu/faculty/greengar/shortcourse_fmm.pdf
Many, Many important applications
• Gravitational Force → Dark matter, cosmology...
Electrostatic
2
∇²φ = ερ
f = ∇¢
Major force
between
moleculars!!
Cancer research
Drug Design
Virus Analysis
Helmholtz
2
2
∇²¢ + k²¢ =− f
Acoustics
0.027
0.026
0.022
0.021
0.0.18
0.012
0.0 18
0.014
0.0 12
0.011
0.008
0.00%
0.006
0.004
0.002
0.000
X
Electromagnetic
P11
•
Boundary Element Method
T₂ →²=0
1
2
==font
u = 4 = g on T₁
Solve for u on 12 and q, on T1, s.t.
N
+Σ
N
dA; = 22; RdAj
اتر
PDE → Integral Equations
Matrix is Dense!
Condition number is Good!(usually converge in O(1)
iterations)
Boundary Element Method
• In 2D
• Full Domain is N^2, with Multigrid Methods, O(N^2) computation.
• Boundary element has N elements, BiCGSTAB method converges in
constant iteration, and each iteration took O(N^2) for Dense Matrix-
vector!
• In 3D
• Full Domain is N^3, with Multigrids, O(N^3)
• Boundary element has N^2 elements, in total N^4 computation!
• With FMM replacing the matrix vector multiplication operation
• O(N) in 2D
• O(N^2) in 3D.
• Semi-Analytic!
Large scale, deep
wave
Wavepacket
source
surface
(Pa) (Pa) (dB)
Pr
PT2г
BDDD
6000
4000
2000
0
-2000
-4000
Insertion loss (IL)
on far-field sphere
(drawn smaller)
Real part of
wavepacket total
field (PT2r)
Sound Pressure
-6000
-8000
187654321012
Real part of
wavepacket
prescribed field (Pr)
Vortex flow
Ui
=
a
N
-
vjwj × (xi − xj)
Σ
j=1.j≠i 4π||xi - xj||
3
2
2
∇²ψ =− ω
u = ∇ × ψ
Type of summations
N
√xp(x) = ½ mj dij
3=1
3
R
BIE: ㄥˇ庁口(応)叮」
Biot-Savart
j=1
리~
wjx dij duj
R
Compute them with same routine
• Given Routine:
• F = computeGradPhi(s); //returns gravity force using FMM from many source s
• Β.Ι.Ε:
• Let s1[i] = n[i].x*s[i],
s2[i] = n[i].y*s[i],
s3[i] = n[i].z*s[i];
• F1 = computeGradPhi(s1), F2 = computeGradPhi(s2), F3 = computeGradPhi(s3)
• Res = F1.x + F2.y + F3.z;
• Biot-Savart:
• Your homework
Other fast summation methods:
• PPPM
PPPM: Combining PDE form and summation
forms
Ui
=
N
-
vjwj × (xi − xj)
j=1,j≠i 4π||xi - xj||
fi=-∈
N
Σ
-
3
2
PjVj (xi − xj)
j=1;j≠i 4π||xi - xj||
Xi
Direct summation for
the turbulent part
3
2
2
∇²ψ =− ω
u = ∇ × ψ
2
∇²φ = ερ
f = ∇¢
Poisson's Equation
for the smooth part
PPPM
• Fast solution uses near-far decomposition to get acceleration.
Can we do similar thing on a particle-mesh setup?
Far field construction
UCO
gipbal
Find bounding box of vortex
Determine the
global domain
Assign particle values to
grid, compute B.C.
Solve Poisson system to
get global velocity field
For each gridcell, cancel the
short range contribution
27
PPPM
• Fast solution uses near-far decomposition to get acceleration.
Can we do similar thing on a particle-mesh setup?
Velocity evaluation
Ufar
Ufinal
UCO
gipbal
Find bounding box of vortex
Dete
glok
Get far field velocity
by interpolation
Add near field direct sum- on system to
mation to get final velocity velocity field
For each gridcell, cancel the
short range contribution
28
timging
10
measured performance
O(N)
103
O(N2)
102
10
10°
10-1
10-2
-3
10
Performance of the PPPM fast summation method
accuracy of PPPM summation
10
4
10
102
3
10
10
number of N-bodies
105
6
10
10
Figure 4: Performance of the PPPM fast summation. Computation
time grows linearly with the number of computational elements.
ا
average error compare to direct summation
10
0-2
10
0-3
0
0-1
10
256 vortex particle
2048 vortex particle
16384 vortex particle
less than 1%
0.5
1
1.5
2
2.5
3
Local correction range K
Figure 5: Accuracy statistics of the PPPM fast summation.
29
Local Correction
• In 3D, for a correction window of size K in each dimension, a
local matrix of size K3 × K³ can be precomputed to cancel the
local influence from grid.
• T(N) = O( c K^6 N)
30
Local Correction
31
Local Correction
The influence made by
neighbor cells.
32
Local Correction
• The matrix inverse reveals
how the center cell's value
depends linearly on its
neighbors(including itself).
Sc =
Σ
jen
ajrj
33
PPPM in few lines
• w_bar = particle_to_grid(w_p);
• dw = w_p - interpolate(w_bar);
• Psi = Poisson.Solve(w_bar);
• v_smooth = curl(Psi);
• v_p = interpolate(v_smooth) + nearSum(dw);
Summarize
• Fast Summation Methods
• FMM
• PPPM
• Equations solved by Fast Summation Methods:
• Poisson's Equation
• Laplace Equation
• Helmholtz Equation
• BEM
• Applications of Fast Summation Methods
• Electrostatic::Molecular Dynamics::Cancer, drug design research
Magnetics::Ship design
• Acoustics::Urban planning, vehicle shape design, theatre design
• Potential flow::aircraft, wave
• Vortex method::turbulent flow
```

