(* ::Package:: *)

slns=k /. Solve[k^2 -2 k + 1 == 1/12 (h^2)(10)(k^2 + 10 k + 1 ) , k];
slns
Plot[Abs[slns],{h, -10, 10}]











