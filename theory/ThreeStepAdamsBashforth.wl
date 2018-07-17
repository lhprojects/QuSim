(* ::Package:: *)

slns = k /. Solve[k^3 == (1 + 23/12 I h)k^2 - 4/3 I h k + 5/12 I h ,k ];
{a,b,c} = Abs[slns];
Plot[{a,b,c},{h,-1,1}]

