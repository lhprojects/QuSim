(* ::Package:: *)

A =  { {5/36, 2/9 - 1/15 Sqrt[15], 5/36 - 1/30 Sqrt[15]},
{5/36 + 1/24 Sqrt[15], 2/9 , 5/36 - 1/24 Sqrt[15]},
  {5/36 + 1/30 Sqrt[15], 2/9 + 1/15 Sqrt[15], 5/36}   };
(* {k1, k2} = I (-K)(1 + h A {k1, k2})*)
Print["{k1,k2,k3}"]
{k1, k2, k3} = Inverse[IdentityMatrix[3]/(I (-K))- h A] . {1, 1, 1}
GL6[h] = 1 + 5/18 k1 h + 4/9 k2 h + 5/18 k3 h;

Print["GL6"]
Simplify[GL6[h]]

Print["GL6 Series"]
Series[GL6[h], {h, 0, 8}]

Print["Exp[-I K h]"]
Series[Exp[-I K h], {h , 0, 8}]














