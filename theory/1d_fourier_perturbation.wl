(* ::Package:: *)

(*GlobalPrecision=MachinePrecision*)
GlobalPrecision=MachinePrecision;
epsilon=N[1/20, GlobalPrecision];
hbar=1;
(* particle *)
mass=1;
PotWidth=1;

EnCont = 1/2;
k0Cont = Sqrt[2 EnCont/mass]/hbar;

(* range *)
M=400;
x0=-80;
x1=80;
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
ListPlot[{Re[G0k[E0]],Im[G0k[E0]]}, PlotRange->Full, PlotLegends->"G0 (k)"]






(* ::InheritFromParent:: *)
(**)


(* ::InheritFromParent:: *)
(**)


Psi0 = N[Table[Exp[(ix dx + x0) k0Cont I], {ix, 0, M}], GlobalPrecision];

(* potential *)
Potential[x_] := N[Exp[-((x-center)/PotWidth)^2],GlobalPrecision];
PotUnit := Table[Potential[i dx + x0],{i,0,M}];
ListPlot[PotUnit, PlotRange->Full, PlotLegends->"V (x)"]

ListPlot[{Re[MyFourier[PotUnit Psi0]], Im[MyFourier[PotUnit Psi0]]}, PlotRange->Full, PlotLegends->"V Psi0 (k)"]


G0Vx[psix_] := MyInverseFourier[G0k[E0] MyFourier[-PotUnit psix]];


O1x[psix_] := G0Vx[psix];
O2x[psix_] := G0Vx[psix + G0Vx[psix]];
O3x[psix_] := G0Vx[psix + G0Vx[psix + G0Vx[psix]]];


ListPlot[Abs[O1x[1/10 Psi0]]]
ListPlot[Arg[O1x[1/10 Psi0]]]

a=Join[Take[G0k[E0] MyFourier[-PotUnit Psi0], 200], Table[0, M+1-200]];
ListPlot[Abs[MyInverseFourier[a]], PlotRange->Full]











