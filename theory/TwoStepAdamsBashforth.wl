(* ::Package:: *)

slns=k /. Solve[k^2 == (1 + 3/2 I h)k - 1/2 I h ,k];
Plot[Abs[slns],{h,-1,1}]


