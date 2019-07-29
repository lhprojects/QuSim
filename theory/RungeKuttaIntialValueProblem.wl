(* ::Package:: *)

Clear[a1];
Clear[a2];
Clear[a3];
Clear[a4];

A1={{0,1},{a1,0}} dx;
A2={{0,1},{a2,0}} dx;
A3={{0,1},{a3,0}} dx;
A4={{0,1},{a4,0}} dx;

I2 = IdentityMatrix[2];
Sym = {{0,1},{-1,0}};


Print["Midpoint"]
K = A1 . Inverse[I2 - 1/2 A1];
T=I2+K;
Print[T]
Print["Check correctness"]
Print[Simplify[K - A1 . (I2 + 1/2 K)]]
Print["check sympl"]
Print[Simplify[Transpose[T] . Sym . T - Sym]]


Print["GLO4"]
a11 = 1/4;
a22 = 1/4;
a12 = 1/4 - Sqrt[3]/6;
a21 = 1/4 + Sqrt[3]/6;
K1 = A1 . Inverse[I2 - a11 A1 - a22 A2 +(a11 a22 - a12 a21) A2 . A1] . (I2 + (a12 - a22)A2);
K2 = A2 . Inverse[I2 - a11 A1 - a22 A2 +(a11 a22 - a12 a21) A1 . A2] . (I2 + (a21 - a11)A1);
T=I2+1/2(K1+K2);
Print[T]
Print["Check correctness"]
Print[Simplify[K1 - A1 . (I2 + a11 K1 + a12 K2)]]
Print[Simplify[K2 - A2 . (I2 + a21 K1 + a22 K2)]]
Print["check sympl"]
Print[Simplify[Transpose[T] . Sym . T - Sym]]



