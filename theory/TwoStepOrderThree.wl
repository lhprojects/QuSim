(* ::Package:: *)

slns = k /. Solve[k^2 == (-4 + 4 I h)k  + 5 + 2 I h ,k ];
Series[slns,{h, 0,4}]
Series[Exp[I h],{h, 0,4}]
slns = Abs[slns];
Plot[slns,{h,-1,1}]




