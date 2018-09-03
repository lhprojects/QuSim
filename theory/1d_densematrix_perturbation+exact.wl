(* ::Package:: *)

ClearAll[i]
ClearAll[Rn]
(* constants *)
(*GlobalPrecision=MachinePrecision*)
GlobalPrecision=80;
hbar=1;
(* particle *)
mass=1;
En=1/2;
k=Sqrt[2 mass En]/hbar;
(* range *)
M=100;
x1=20;
x0=0;
dx=(x1-x0)/M;
center=(x1+x0)/2;
(* potential *)
Potential[x_, p_] := p Exp[-(x-center)^2];

Pot[p_] := N[Table[Potential[i dx, p],{i,0,M}],GlobalPrecision];
V[p_] := DiagonalMatrix[Pot[p]];
G0 = -Table[Table[N[Exp[ Abs[i-j] dx k I], GlobalPrecision],{i, 0, M}],{j, 0, M}] (mass/(hbar^2 I k));
T1[p_] := -G0 . V[p]  dx;


ListLinePlot[Re[Pot[1]], PlotRange->Full]

Psi0 = N[Table[Exp[k i dx I], {i, 0, M}], GlobalPrecision];

Plot[{ Re[Eigenvalues[T1[p], 2]], Im[Eigenvalues[T1[p], 2]] },{p,0,1}]


EigenVec=Eigenvalues[T1[1], 1];
EigenVal=Eigenvectors[T1[1], 1];


EigenT=Eigenvectors[T1[1],1];
ListLinePlot[Re[EigenT[[1]]]]
ListLinePlot[Im[EigenT[[1]]]]

ORDER=4;
Psin[p_] := Table[MatrixPower[T1[p], i] . Psi0, {i, 1, ORDER}];
Sn[p_] := Table[Sum[Psin[p][[j]] ,{j, 1, i}], {i, 1, ORDER}];
Rn[p_] := Abs[Map[First, Sn[p]]]^2;


POINTS=8;
Exact = {
0.00586814,
0.0312010,
0.0887220,
0.1866818,
0.3198141,
0.4677157,
0.6060175,
0.7193307
};
PsiInf[p_] := (Inverse[IdentityMatrix[M+1] - T1[p]] - IdentityMatrix[M+1]) . Psi0;
RInf[p_] := Abs[First[PsiInf[p]]]^2;

RnList=Table[    Table[Rn[0.1 i][[o]],{i,1,POINTS}],               {o,1,ORDER}];
RnList=Append[RnList,   Table[RInf[0.1 i], {i,1,POINTS}]         ];
RnList;
RnNormList=Table[Table[RnList[[o]][[i]]/Exact[[i]], {i,1,POINTS}], {o,1,ORDER+1}];
RnNormList=Append[RnNormList,Table[1,{i,1,POINTS}]];
RnNormList;

ListPlot[Pot[1], PlotRange->Full]
ListLogPlot[RnNormList,PlotLegends->{"O1","O2","O3", "O4","Inf","Exact"},PlotRange->Full, PlotMarkers->{Automatic, Large}]




(* ::InheritFromParent:: *)
(**)


(* ::InheritFromParent:: *)
(**)


(* ::InheritFromParent:: *)
(**)


(* ::InheritFromParent:: *)
(**)
