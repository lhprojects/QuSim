(* ::Package:: *)

(*GlobalPrecision=MachinePrecision*)
GlobalPrecision=MachinePrecision;
epsilon=N[1/10, GlobalPrecision];
hbar=1;
(* particle *)
mass=1;
PotWidth=1;

EnCont = 1/2;
k0Cont = Sqrt[2 EnCont/mass]/hbar;

(* range *)
M=200;
x0=-20;
x1=20;
dx=(x1-x0)/(M+1);
center=N[(x1+x0)/2,GlobalPrecision];

dk = 2 Pi/(x1 - x0);
k0 = N[Floor[k0Cont (x1-x0)/(2 Pi)] 2 Pi /(x1 - x0), GlobalPrecision];
E0 = N[1/2 mass hbar^2 k0^2, GlobalPrecision];
Print["dk ", N[dk]]
Print["E0 ", N[E0]]

MinK=2. Pi/(x1 - x0);
MaxK=2. Pi/(x1 - x0) M;
ErrBoundary = Exp[-1/2 epsilon /E0 k0 (x1 - x0)];
Print["ErrBound decay ", ErrBoundary]

MyFourier[x_] := Fourier[x, FourierParameters->{0,-1}]; (*  -   *)
MyInverseFourier[x_] := InverseFourier[x, FourierParameters->{0,-1}]; (*  +   *)


G0k[en_] := Table[1/(1/2 mass hbar^2 (2 Pi If[ki < M/2, ki, M - ki] /(x1 - x0))^2 - en - I epsilon), {ki,0,M}];
ListPlot[Abs[Take[G0k,M+1]],PlotRange->Full, PlotLegends->"G0 (k) abs"]
ListPlot[Arg[Take[G0k,M+1]],PlotRange->Full, PlotLegends->"G0 (k) arg"]






Psi0 = N[Table[Exp[(ix dx + x0) k0 I], {ix, 0, M}], GlobalPrecision];
(* potential *)
Potential[x_] := N[-Exp[-x^2], GlobalPrecision];
PotUnit := Table[Potential[i dx + x0],{i,0,M}];
ListPlot[PotUnit, PlotRange->Full, PlotLegends->"V"]
ListPlot[Abs[MyFourier[PotUnit Psi0]], PlotRange->Full, PlotLegends->"V Psi0 (k)"]

G0Vx[en_, psix_] := MyInverseFourier[G0k[en] MyFourier[-PotUnit psix]];

O1x[en_, M_, psix_] := 1/M G0Vx[en, psix];
O2x[en_, M_, psix_] := 1/M G0Vx[en, psix + 1/M G0Vx[en, psix]];
O3x[en_, M_, psix_] := 1/M G0Vx[en, psix + 1/M G0Vx[en, psix + 1/M G0Vx[en, psix]]];

O3Generalx[en_, MM_, M_, psix_] := If[ MM == 0,
   O3x[en, M, psix],
   1/M G0Vx[en, MM-1, M, psix + 1/M G0Vx[en, MM-1, M, psix + 1/M G0Vx[en, MM - 1, M, psix]]]
 ];


ListPlot[Abs[O3Generalx[E0, Psi0]], PlotLegends->"Psi"]













