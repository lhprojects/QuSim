(* ::Package:: *)

ClearAll[f0, f1, f2, f3 ]
ClearAll[g0, g1, g2, g3 ]
ClearAll[h0, h1, h2, h3 ]
ClearAll[l0, l1, l2, l3 ]
ClearAll[c1,c2,c3,c4]
ClearAll[b1,b2,b3,b4]
ClearAll[a21,a32,a43,a42]

I2 = IdentityMatrix[2];



A={{0,1},{f0, 0}};
dAdx={{0, 0},{f1, 0}};
d2Ad2x={{0, 0},{f2, 0}};
d3Ad3x={{0, 0},{f3, 0}};
Generate[c_] := {{0,1}, {f0 + f1 (c h) + 1/2 f2 (c h)^2 + 1/6 f3 (c h)^3,0}};


(*
A={{1,0},{0, f0}};
dAdx={{0, 0},{0, f1}};
d2Ad2x={{0, 0},{0, f2}};
d3Ad3x={{0, 0},{0, f3}};
Generate[c_] := {{1,0},{0, f0 + f1 (c h) + 1/2 f2 (c h)^2 + 1/6 f3 (c h)^3 + 1/24 f4 (c h)^4}};
*)


A={{h0, l0},{g0, f0}};
dAdx={{h1, l1},{g1, f1}};
d2Ad2x={{h2, l2},{g2, f2}};
d3Ad3x={{h3, l3},{g3, f3}};
Generate[c_] := {{h0 + h1 (c h) + 1/2 h2 (c h)^2 + 1/6 h3 (c h)^3 + 1/24 h4 (c h)^4, l0 + l1 (c h) + 1/2 l2 (c h)^2 + 1/6 l3 (c h)^3 + 1/24 l4 (c h)^4},
{g0 + g1 (c h) + 1/2 g2 (c h)^2 + 1/6 g3 (c h)^3 + 1/24 g4 (c h)^4, f0 + f1 (c h) + 1/2 f2 (c h)^2 + 1/6 f3 (c h)^3 + 1/24 f4 (c h)^4}};



Print["Begin Exact"]
dydx = A;
d2yd2x = dAdx + A . A;
d3yd3x = d2Ad2x + A . dAdx + dAdx . A + (dAdx + A . A) . A;
d4yd4x = d3Ad3x + dAdx . dAdx + A . d2Ad2x + d2Ad2x . A + dAdx . dAdx +
  (d2Ad2x + A . dAdx + dAdx . A). A + (dAdx + A . A) . dAdx  + d3yd3x . A;

Print["Exact"]
Exact = I2 + h dydx + 1/2 h^2 d2yd2x + 1/6 h^3 d3yd3x + 1/24 h^4 d4yd4x;


Print["MidPoint"]
A1=Generate[1/2];
W=(I2 + 1/2 A1 h). Inverse[I2 - 1/2 A1 h];
Print["Series"]
Series[ W - Exact, {h,0,5}];
Print["check simplic"]
Simplify[Transpose[W] . {{0,1},{-1,0}} . W - {{0,1},{-1,0}}]


A1={{0,1},{a,0}};
W=(I2 + 1/2 A1 h). Inverse[I2 - 1/2 A1 h];

h=1;
Print["decay"]
{eigen1,eigen2}=Eigenvalues[W];
Plot[Eigenvalues[W], {a, 0, 10}]
Print["osc"]
Plot[Abs[Eigenvalues[W]], {a, -100, 0}]
ClearAll[h];



Print["Begin RK4 Classical"]

A1=Generate[0];
A2=Generate[1/2];
A3=Generate[1/2];
A4=Generate[1];

k1 = A1;
k2 = A2 .(I2 + h 1/2 k1);
k3 = A3 .(I2 + h 1/2 k2);
k4 = A4 .(I2 + h k3);
ClassicalRK4 = I2 + h(1/6 k1 + 1/3 k2 + 1/3 k3 + 1/6 k4);

Print["Classical RK4 Series"];
Series[ClassicalRK4, {h, 0, 4}];

Print["Classical RK4 Diff Series"];
Simplify[Series[ClassicalRK4 - Exact,{h, 0, 4}]]



Print["Begin Check simplic"];
A1={{0,1},{a1,0}};
A2={{0,1},{a2,0}};
A3={{0,1},{a3,0}};
A4={{0,1},{a4,0}};
k1 = A1;
k2 = A2 .(I2 + h 1/2 k1);
k3 = A3 .(I2 + h 1/2 k2);
k4 = A4 .(I2 + h k3);
ClassicalRK4 = I2 + h(1/6 k1 + 1/3 k2 + 1/3 k3 + 1/6 k4);
Series[Simplify[Transpose[ClassicalRK4] . {{0,1},{-1,0}} . ClassicalRK4 - {{0,1},{-1,0}}],{h,0,5}]


Print["check stability"]
h=1;
A1={{0,1},{a,0}};
A2={{0,1},{a,0}};
A3={{0,1},{a,0}};
A4={{0,1},{a,0}};

k1 = A1;
k2 = A2 .(I2 + h 1/2 k1);
k3 = A3 .(I2 + h 1/2 k2);
k4 = A4 .(I2 + h k3);
ClassicalRK4 = I2 + h(1/6 k1 + 1/3 k2 + 1/3 k3 + 1/6 k4);
Plot[Abs[Eigenvalues[ClassicalRK4]],{a,-10,0}]
ClearAll[h, a]

Print["End Classical BK4"]


Print["Begin GLO4"]
A1={{0,1},{a1, 0}};
A2={{0,1},{a2, 0}};

(*k1 = A1 .(I2 + h 1/4 k1 + h (1/4-Sqrt[3]/6) k2 );
k2 = A2 .(I2 + h (1/4+Sqrt[3]/6) k1 + 1/4 k2);*)

AA=Inverse[A1 ]- 1/4 h I2;
BB=-(1/4 - Sqrt[3]/6) h I2;
CC=-(1/4 + Sqrt[3]/6) h I2;
DD=Inverse[A2] - 1/4 h I2;
kk1 = Inverse[Inverse[BB] . AA - Inverse[DD] . CC] . (Inverse[BB]-Inverse[DD]);
kk2 = Inverse[Inverse[AA] . BB - Inverse[CC] . DD] . (Inverse[AA]-Inverse[CC]);
Print["Check Zero"]
Simplify[kk1 - A1 .(I2 + h 1/4 kk1 + h (1/4-Sqrt[3]/6) kk2 )]
Simplify[kk2 - A2 .(I2 + h (1/4+Sqrt[3]/6) kk1 + h 1/4 kk2)]
GLO4 = I2 + h(1/2 kk1 + 1/2 kk2);

Print["Check simplic"]
Simplify[Transpose[GLO4] . {{0,1},{-1,0}} . GLO4 - {{0,1},{-1,0}}]


Print["Check stability"]
A1={{0,1},{a, 0}};
A2={{0,1},{a, 0}};
h=1;
(*k1 = A1 .(I2 + h 1/4 k1 + h (1/4-Sqrt[3]/6) k2 );
k2 = A2 .(I2 + h (1/4+Sqrt[3]/6) k1 + 1/4 k2);*)

AA=Inverse[A1 ]- 1/4 h I2;
BB=-(1/4 - Sqrt[3]/6) h I2;
CC=-(1/4 + Sqrt[3]/6) h I2;
DD=Inverse[A2] - 1/4 h I2;
kk1 = Inverse[Inverse[BB] . AA - Inverse[DD] . CC] . (Inverse[BB]-Inverse[DD]);
kk2 = Inverse[Inverse[AA] . BB - Inverse[CC] . DD] . (Inverse[AA]-Inverse[CC]);
GLO4 = Simplify[I2 + h(1/2 kk1 + 1/2 kk2)];
Print["decay"]
Plot[Eigenvalues[GLO4], {a, 0, 10000}]
Print["osc"]
Plot[Abs[Eigenvalues[GLO4]], {a, -10000, 0}]
ClearAll[h];

A1=Generate[1/2-Sqrt[3]/6];
A2=Generate[1/2+Sqrt[3]/6];

k1 = A1 .(I2 + h 1/4 k1 + h (1/4-Sqrt[3]/6) k2 );
k2 = A2 .(I2 + h (1/4+Sqrt[3]/6) k1 + 1/4 k2);

AA=Inverse[A1 ]- 1/4 h I2;
BB=-(1/4 - Sqrt[3]/6) h I2;
CC=-(1/4 + Sqrt[3]/6) h I2;
DD=Inverse[A2] - 1/4 h I2;
kk1 = Inverse[Inverse[BB] . AA - Inverse[DD] . CC] . (Inverse[BB]-Inverse[DD]);
kk2 = Inverse[Inverse[AA] . BB - Inverse[CC] . DD] . (Inverse[AA]-Inverse[CC]);

(*
Simplify[kk1 - A1 .(I2 + h 1/4 kk1 + h (1/4-Sqrt[3]/6) kk2 )]
Simplify[kk2 - A2 .(I2 + h (1/4+Sqrt[3]/6) kk1 + h 1/4 kk2)]
*)

GLO4 = I2 + h(1/2 kk1 + 1/2 kk2);

Print["GLO4 Diff Series"];
(*Series[GLO4- Exact, {h, 0, 4}]*)










