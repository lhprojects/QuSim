(* ::Package:: *)


A =  { {5/36, 2/9 - 1/15 Sqrt[15], 5/36 - 1/30 Sqrt[15]},
{5/36 + 1/24 Sqrt[15], 2/9 , 5/36 - 1/24 Sqrt[15]},
  {5/36 + 1/30 Sqrt[15], 2/9 + 1/15 Sqrt[15], 5/36}   };
(* {k1, k2} = I K(1 + h A {k1, k2})*)
IdentityMatrix[3]/(I K)- h A;
Inverse[IdentityMatrix[3]/(I K)- h A];
{k1, k2, k3} = Inverse[IdentityMatrix[3]/(I K)- h A] . {1, 1, 1};
Simplify[k1];
Simplify[k2];
Simplify[k3];
W = 1 + 5/18 k1 h + 4/9 k2 h + 5/18 k3 h;
Simplify[W]
WW[h_] := (120 I-60 h-12 I h^2+h^3)/(120 I+60 h-12 I h^2-h^3);
Series[WW[h], {h , 0, 8}]
Series[Exp[I h], {h , 0, 8}]





