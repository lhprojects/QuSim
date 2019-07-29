(* ::Package:: *)


GL2[h] = (1-1/2 I h K)/(1+1/2 I h K);
Print["GL2:"]
GL2[h]

Print["GL2:"]
Series[GL2[h],{h, 0, 5}]

Print["Series Exp[-I h]:"]
Series[Exp[I (-K) h],{h, 0, 5}]





















