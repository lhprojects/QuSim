(* ::Package:: *)


slns=k /. Solve[k^2 == 5/12 I h k^2 + (1 + 2/3 I h)k - 1/12 I h, k];

Plot[Abs[slns],{h,-1,1}]
Series[slns,{h,0,4}]
Series[Exp[h],{h,0,4}]






