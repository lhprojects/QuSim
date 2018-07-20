(* ::Package:: *)

data={
{1.5, 1.4 10^-6},
{2.0, 1.2 10^-8},
{3.0, 5.0 10^-12},
{4.0, 7.9 10^-15}
};

CC[p_] := Exp[-2 NIntegrate[ Sqrt[2(p Exp[-0.1 x^2] - 0.50)],{x, -Sqrt[-Log[0.50/p]/0.1], +Sqrt[-Log[0.50/p]/0.1]}]]
WKB[p_] := CC[p] / (1 + 0.25 CC[p]);
wkb={
{1.5,   WKB[1.5]},
{2,   WKB[2]},
{3.0, WKB[3.0]},
{4.0, WKB[4.0]}
};

ListLogPlot[{data, wkb}]


e=2;

data={
{2.0, 0.57},
{3.0, 0.078},
{4.0, 0.0089},
{5.0, 0.0013}
};


CC[p_] := Exp[-2 NIntegrate[ Sqrt[2(p Exp[-x^2] - e)],{x, -Sqrt[-Log[e/p]], +Sqrt[-Log[e/p]]}]]
WKB[p_] := CC[p] / (1 + 0.25 CC[p]);
wkb={
{2,   WKB[2]},
{3,   WKB[3]},
{4,   WKB[4]},
{5,   WKB[5]}
};

ListLogPlot[{data, wkb}]




e=4;

data={
{4.0, 0.55},
{5.0, 0.134},
{6.0, 0.0241},
{7.0, 0.0047},
{8.0, 0.0010}
};

CC[p_] := Exp[-2 NIntegrate[ Sqrt[2(p Exp[-x^2] - e)],{x, -Sqrt[-Log[e/p]], +Sqrt[-Log[e/p]]}]]
WKB[p_] := CC[p] / (1 + 0.25 CC[p]);
wkb={
{4,   WKB[4]},
{5,   WKB[5]},
{6,   WKB[6]},
{7,   WKB[7]},
{8,   WKB[8]}
};

ListLogPlot[{data, wkb}]







