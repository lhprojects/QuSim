(* ::Package:: *)

A =  { {1/4, 1/4 - 1/6 Sqrt[3]},  {1/4 + 1/6 Sqrt[3], 1/4}   };
(* {k1, k2} = I (-K)(1 + h A {k1, k2})*)
Print["{k1, k2}:"]
{k1, k2} = Inverse[IdentityMatrix[2]/(-I K)- h A] . {1, 1}
GL4[h_] = 1 + 1/2 k1 h + 1/2 k2 h;

Print["GL4:"]
Simplify[GL4[h]]

Print["GL4:"]
Series[GL4[h],{h, 0, 5}]

Print["Series Exp[-I h]:"]
Series[Exp[I (-K) h],{h, 0, 5}]

Simplify[
(1 - 1/12 (3 I-Sqrt[3])h K)(1 - 1/12 (3 I+Sqrt[3])h K) / ((1 - 1/12 (-3 I-Sqrt[3])h K)(1 - 1/12 (-3 I+Sqrt[3])h K))
]


















