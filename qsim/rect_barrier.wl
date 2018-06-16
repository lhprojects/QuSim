(* ::Package:: *)

(* system parameter*)
hbar=1;

(* particle parameter*)
k=1;
m=1;

(* rectangular barrier*)
V0=1;
a=1;
en=k^2 hbar^2/(2 m);

k0=Sqrt[2 m en/hbar^2];
k1B=Sqrt[2 m (en-V0)/hbar^2];
k1A=Sqrt[2 m (V0-en)/hbar^2];

T=If[V0==en, 1/(1 + m a^2 V0 / (2 hbar^2)),
If[V0 > en, 1/(1 + V0^2 Sinh[k1A a]^2 / (4 en (V0 - en))),
1/(1 + V0^2 Sin[k1B a]^2/(4 en (en - V0)))]
];

Print["V0: ", V0]
Print["en: ", en]
Print["en/V0: ", en/V0]
Print["T: ", N[T]]









