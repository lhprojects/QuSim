(* ::Package:: *)

Exact[x_]:=1/(1-x);
O2[x_]:= 1+ x+x^2;
O3[x_]:= 1+ x+x^2 + x^3;

O1General[x_,c_]:=  c + c^2 x;
O2General[x_,c_]:=  c + c^2 x + c^3 x^2;
O3General[x_,c_]:=  c + c^2 x + c^3 x^2 +c^4 x^3;

O1Dn[x_,n_] := For[i=0;t=1,i < n, ++i,
t=O1General[x/n,t];
If[i==n-1, Return[t]]
];

O2Dn[x_,n_] := For[i=0;t=1,i < n, ++i,
t=O2General[x/n,t];
If[i==n-1, Return[t]]
];

O3Dn[x_,n_] := For[i=0;t=1,i < n, ++i,
t=O3General[x/n,t];
If[i==n-1, Return[t]]
];



Fnl=1+I
WhatFuck[f_, x_, n_] := Re[f[x Fnl,n]]
WhatFuck2[f_, x_] := Re[f[x Fnl]]
rng=5;

Plot[{ 
WhatFuck2[Exact, x],
WhatFuck[O3Dn, x, 4],
WhatFuck[O3Dn, x, 16]
},
 {x,-rng,rng},PlotRange->{-10,10}, PlotLegends->"Expressions"]

Plot[{ 
WhatFuck2[Exact, x],
WhatFuck[O2Dn, x, 4],
WhatFuck[O2Dn, x, 16],
},
 {x,-rng,rng},PlotRange->{-5,5}, PlotLegends->"Expressions"]
 
 Plot[{ 
WhatFuck2[Exact, x],
WhatFuck[O1Dn, x, 4],
WhatFuck[O1Dn, x, 16],
WhatFuck[O1Dn, x, 1024]
},
 {x,-rng,rng},PlotRange->{-10,10}, PlotLegends->"Expressions"]

