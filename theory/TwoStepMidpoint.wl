(* ::Package:: *)

slns=k /. Solve[k^2 == I h k + 1 ,k];
slns
Plot[Abs[slns],{h,-4,4}]





