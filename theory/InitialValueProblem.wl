(* ::Package:: *)

ClearAll[psi, x]
En=1/2;
k=Sqrt[2 En];
initialCondition={psi[-10]==1, Derivative[1][psi][-10]==-k I};

sln=NDSolveValue[{-1/2 D[psi[x],{x,2}] + 10 Exp[-x^2] psi[x] == En psi[x], initialCondition},
 psi[x], {x,-10,10}, PrecisionGoal->35,WorkingPrecision->50,MaxSteps->Infinity,MaxStepSize->0.01,
  Method->{"ExplicitRungeKutta","DifferenceOrder"->9}];
dsln=D[sln, x];

x=10;
B=(sln-dsln/(I k))/2;
N[1/Abs[B]^2, 25]
Plot[{Exp[-x^2], Abs[sln]},{x,-10,10}]




0.1355284179587569045304922`25.
8.666893019856854462914111 10^-9




