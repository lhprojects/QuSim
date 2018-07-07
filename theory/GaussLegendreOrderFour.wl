(* ::Package:: *)

A =  { {1/4, 1/4 - 1/6 Sqrt[3]},  {1/4 + 1/6 Sqrt[3], 1/4}   };
(* {k1, k2} = I K(1 + h A {k1, k2})*)
IdentityMatrix[2]/(I K)- h A;
Inverse[IdentityMatrix[2]/(I K)- h A]
{k1, k2} = Inverse[IdentityMatrix[2]/(I K)- h A] . {1, 1}
Simplify[k1];
Simplify[k2];
W = 1 + 1/2 k1 h + 1/2 k2 h;
Simplify[W]

WW[h_]:=(-12-6 I h +h^2)/(-12+6 I h+h^2);
Series[WW[h],{h, 0, 5}]
Series[Exp[I h],{h, 0, 5}]




