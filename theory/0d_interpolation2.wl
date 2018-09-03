(* ::Package:: *)


inv[x_] := 1/(1-x);
Dinv[x_] := 1/(1-x)^2;

ORDER=8;
POINTS=ORDER/2;
x0=-2;
x1=+2;
valuePairs = Table[{x, inv[x]}, {x, x0+(x1-x0)/, x1, (x1 - x0)/POINTS}];
Interpolation[valuePairs];

Mat=Table[Table[x^i,{i, 0, ORDER}], {x, -2+1/20, 2, 1/10}];
ExactValues = Table[inv[x], {x, -2+1/20, 2, 1/10}];

Mat1=Table[Table[i x^(i-1),{i, 0, ORDER}], {x, -2+1/20, 2, 1/10}];
ExactValues1 = Table[Dinv[x], {x, -2+1/20, 2, 1/10}];
TotalMat = Join[Mat, Mat1];
TotalExactValues = Join[ExactValues, ExactValues1];

Coef=Inverse[TotalMat] . TotalExactValues;


ExactValues
ExactValues1
(* numerical *)
roots = x /. NSolve[Table[x^i, {i, 0, 8*10-1}] . Coef == 0, x]
norm := For[i=1;t=1,  i <= ORDER,   i=i+1,
t=t*(-roots[[i]]);
If[i==ORDER-1, Return[t]]
];

inter2[x_] := For[i=1;t=1,  i <= ORDER,   i=i+1,
t=t*(x-roots[[i]]);
If[i==ORDER-1, Return[t]]
];
inter[x_] := inter2[x]/Abs[norm];


Plot[{inv[x], Abs[inter[x]]}, {x,-0.5, 0.5}, PlotRange->6]




