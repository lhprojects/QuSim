(* ::Package:: *)

A =  { {0, 0, 0, 0},
{1/2, 0, 0, 0},
{0, 1/2, 0, 0},
{0, 0, 1, 0} };
  
(* {k1, k2} = I K(1 + h A {k1, k2})*)
IdentityMatrix[4]/(I K)- h A;
Inverse[IdentityMatrix[4]/(I K)- h A];
{k1, k2, k3, k4} = Inverse[IdentityMatrix[4]/(I K)- h A] . {1, 1, 1, 1};
Simplify[k1];
Simplify[k2];
Simplify[k3];
Simplify[k4];
W = 1 + 1/6 k1 h + 1/3 k2 h + 1/3 k3 h + 1/6 k4 h;
Simplify[W]
Series[W,{h, 0, 10}]
Series[Exp[I h K],{h, 0, 10}]
Plot[Abs[1+I h-h^2/2-1/6 I  h^3+ (h^4)/24], {h,-3,3}]








