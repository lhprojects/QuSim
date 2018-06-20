(* ::Package:: *)

Psi[x_] := Sqrt[alpha / Sqrt[Pi]] Exp[I k x] Exp[- alpha^2 x^2 / 2];
PsiConj[x_] := Sqrt[alpha / Sqrt[Pi]] Exp[-I k x] Exp[- alpha^2 x^2 / 2]; 
D[Psi[x],x]
Integrate[Psi[x] PsiConj[x],{x,-Infinity,Infinity}]
Integrate[D[D[Psi[x],x],x ] PsiConj[x],{x,-Infinity,Infinity}]















