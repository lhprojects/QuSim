(* ::Package:: *)


inv[x_] := 1/(1-x);


valuePairs = Table[{x, inv[x]}, {x, -2+1/20, 2, 1/10}];

Interpolation[valuePairs]

ExactValues = Table[inv[x], {x, -2+1/20, 2, 1/10}];
Mat=Table[Table[x^i,{i, 0, 4*10-1}], {x, -2+1/20, 2, 1/10}];
Coef=Inverse[Mat] . ExactValues;
(* numerical *)
roots = x /. NSolve[Table[x^i, {i, 0, 4*10-1}] . Coef == 0, x];
norm := For[i=1;t=1,  i < 40,   i=i+1,
t=t*(-roots[[i]]);
If[i==39, Return[t]]
];

inter2[x_] := For[i=1;t=1,  i < 40,   i=i+1,
t=t*(x-roots[[i]]);
If[i==39, Return[t]]
];
inter[x_] := inter2[x]/Abs[norm]


Plot[{inv[x], Abs[inter[x]]}, {x,-1, 1}, PlotRange->6]




