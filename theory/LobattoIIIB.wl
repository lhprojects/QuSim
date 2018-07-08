(* ::Package:: *)

A =  { {1/6 , -1/6, 0},
{1/6, 1/3, 0},
  {1/6, 5/6, 0}
  };
(* {k1, k2} = I K(1 + h A {k1, k2})*)
IdentityMatrix[3]/(I K)- h A;
Inverse[IdentityMatrix[3]/(I K)- h A];
{k1, k2, k3} = Inverse[IdentityMatrix[3]/(I K)- h A] . {1, 1, 1};
Simplify[k1];
Simplify[k2];
Simplify[k3];
W = 1 + 1/6 k1 h + 2/3 k2 h + 1/6 k3 h;
Simplify[W]

