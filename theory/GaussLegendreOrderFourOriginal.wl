(* ::Package:: *)


A =  { {0, 0, 0, 0},
{1/2, 0, 0, 0},
{0, 1/2, 0, 0},
{0, 0, 1, 0} };
  
(* {k1, k2} = I K(1 + h A {k1, k2})*)
IdentityMatrix[4]/(I K)- h A;
Inverse[IdentityMatrix[4]/(I K)- h A];
{k1, k2, k3, k4} = Inverse[IdentityMatrix[4]/(I K)- h A] . {1, 1, 1, 1};
Simplify[k1];
Simplify[k2];
Simplify[k3];
Simplify[k4];
W = 1 + 1/6 k1 h + 1/3 k2 h + 1/3 k3 h + 1/6 k4 h;
Simplify[W]
Series[W,{h, 0, 10}]
Series[Exp[I h K],{h, 0, 10}]






(* ::Output:: *)
(*SeriesData[h, 0, {1, Complex[0, 1], Rational[-1, 2], Complex[0, Rational[-1, 6]], Rational[1, 24], Complex[0, Rational[1, 120]], Rational[-1, 720], Complex[0, Rational[-1, 4800]], Rational[1, 28800], Complex[0, Rational[7, 864000]], Rational[-1, 432000], Complex[0, Rational[-11, 17280000]], Rational[1, 6480000], Complex[0, Rational[17, 518400000]], Rational[-13, 2073600000], Complex[0, Rational[-71, 62208000000]]}, 0, 16, 1]*)
