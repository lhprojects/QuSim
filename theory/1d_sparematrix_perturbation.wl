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
M=800;
x1=200;
x0=0;
dx=(x1-x0)/(M+1);
center=N[(x1+x0)/2,GlobalPrecision];
K[i]:= 2 Pi i /(x1-x0);

k0 = N[Floor[k0Cont (x1-x0)/(2 Pi)] 2 Pi /(x1 - x0), GlobalPrecision];
En= N[1/2 mass hbar^2 k0^2, GlobalPrecision];


PotCut= 2. PotWidth;
MinK=2. Pi/(x1 - x0)
MaxK=2. Pi/(x1 - x0) M
ErrBoundary = Exp[-1/2 epsilon /En k0 (x1 - x0)]
ErrDecay = 1/2 epsilon /En k0 PotCut
ErrCut = Exp[-(PotCut/PotWidth)^2]


MyFourier[x_] := Fourier[x, FourierParameters->{0,-1}]; (*  -   *)
MyInverseFourier[x_] := InverseFourier[x, FourierParameters->{0,-1}]; (*  +   *)


(* potential *)
Potential[x_] := N[Exp[-((x-center)/PotWidth)^2],GlobalPrecision];
PotUnit := Table[Potential[i dx + x0],{i,0,M}];



Psi0 = N[Table[Exp[(ix dx + x0) k0 I], {ix, 0, M}], GlobalPrecision];
RebIdx[i_]:= If[i<0,  M+1+i, If[i>M,i-M-1, i]]+1;
 
T1=Table[{RebIdx[i], RebIdx[i+1]}-> 1, {i,0,M}];
T2=Table[{RebIdx[i], RebIdx[i-1]}-> 1, {i,0,M}];
T3=Table[{RebIdx[i], RebIdx[i+0]}-> -2, {i,0,M}];
Tt=Join[T1,T2,T3];
T=N[-1/2 mass hbar^2/dx^2, GlobalPrecision] SparseArray[Tt]

V=SparseArray[Table[{RebIdx[i], RebIdx[i]}-> PotUnit[[i+1]], {i,0,M}]];
EM=(En + epsilon I) SparseArray[Table[{RebIdx[i], RebIdx[i]}-> 1, {i,0,M}]]


ASolver = LinearSolve[T-EM];
psi = ASolver[V . Psi0];
ListPlot[Abs[psi], PlotRange->Full]
ListPlot[Arg[psi], PlotRange->Full]





