(* ::Package:: *)

(* ::Title:: *)
(*MST Package*)


(* ::Section::Closed:: *)
(*Create Package*)


BeginPackage[MST`$MasterFunction<>"`MST`MST`", {MST`$MasterFunction<>"`", MST`$MasterFunction<>"`MST`RenormalizedAngularMomentum`", "SpinWeightedSpheroidalHarmonics`"}];

Begin["`Private`"];


(* ::Section::Closed:: *)
(*Utility functions*)


(* ::Subsection::Closed:: *)
(*Continued Fraction*)


(* ::Text:: *)
(*Continued fraction with automatic convergence checking.*)
(*FIXME: There is a potentially better algorithm in Numerical recipes which handles cases where successive terms have large magnitude differences.*)


CF[a_, b_, {n_, n0_}] := Module[{A, B, ak, bk, res = Indeterminate, j = n0},
  A[n0 - 2] = 1;
  B[n0 - 2] = 0;
  ak[k_] := ak[k] = (a /. n -> k);
  bk[k_] := bk[k] = (b /. n -> k);
  A[n0 - 1] = 0(*bk[n0-1]*);
  B[n0 - 1] = 1;
  A[k_] := A[k] = bk[k] A[k - 1] + ak[k] A[k - 2];
  B[k_] := B[k] = bk[k] B[k - 1] + ak[k] B[k - 2];
  While[res =!= (res = A[j]/B[j]), j++];
  Clear[A, B, ak, bk];
  res
];


(* ::Section::Closed:: *)
(*Master function dependent settings*)


(* ::Subsection::Closed:: *)
(*Parameters for the Hypergeometric functions*)


Switch[MST`$MasterFunction,
"ReggeWheeler",
  (* Parameters for Hypergeometric2F1 *)
  aF[s_, \[Nu]_, \[Tau]_, \[Epsilon]_] := \[Nu]+s+1-I \[Epsilon];
  bF[s_, \[Nu]_, \[Tau]_, \[Epsilon]_] := -\[Nu]+s-I \[Epsilon];
  cF[s_, \[Nu]_, \[Tau]_, \[Epsilon]_] := 1-2 I \[Epsilon];

  (* Parameters for HypergeometricU *)
  aU[s_, \[Nu]_, \[Tau]_, \[Epsilon]_] := \[Nu] + 1 - I \[Epsilon];,
 
"Teukolsky",
  (* Parameters for Hypergeometric2F1 *)
  aF[s_, \[Nu]_, \[Tau]_, \[Epsilon]_] := \[Nu]+1-I \[Tau];
  bF[s_, \[Nu]_, \[Tau]_, \[Epsilon]_] := -\[Nu]-I \[Tau];
  cF[s_, \[Nu]_, \[Tau]_, \[Epsilon]_] := 1-s-I(\[Epsilon]+\[Tau]);

  (* Parameters for HypergeometricU *)
  aU[s_, \[Nu]_, \[Tau]_, \[Epsilon]_] := \[Nu] + s + 1 - I \[Epsilon];,

_, Abort[]
];



(* ::Subsection::Closed:: *)
(*Radial solutions*)


Switch[MST`$MasterFunction,
"ReggeWheeler",
  fIn[q_, \[Epsilon]_, \[Kappa]_, \[Tau]_, \[Nu]_, \[Lambda]_, s_, m_, n_] := Pochhammer[-\[Nu]+s-I \[Epsilon],-n]Pochhammer[\[Nu]+s-I \[Epsilon]+1,n]Pochhammer[\[Nu]+I \[Epsilon]+1,n]/Pochhammer[\[Nu]-I \[Epsilon]+1,n](-1)^n fn[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, n];
  prefacIn[s_, \[Epsilon]_, \[Tau]_, \[Kappa]_, x_] := (1-x)^(s+1) (-x)^(-I \[Epsilon]) E^(I \[Epsilon] x);
  fUp[q_, \[Epsilon]_, \[Kappa]_, \[Tau]_, \[Nu]_, \[Lambda]_, s_, m_, n_] := (-1)^n Pochhammer[\[Nu] + 1 + s - I \[Epsilon], n]/Pochhammer[\[Nu] + 1 - s + I \[Epsilon], n] fn[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, n];
  prefacUp[s_, \[Epsilon]_, \[Kappa]_, \[Tau]_, \[Nu]_, zhat_] := 2^\[Nu] E^(-\[Pi] \[Epsilon]) E^(-I \[Pi] (\[Nu]+1)) E^(I zhat) zhat^(\[Nu]+I (\[Epsilon]+\[Tau])/2) (zhat-\[Epsilon] \[Kappa])^(-I (\[Epsilon]+\[Tau])/2) zhat;,
"Teukolsky", 
  fIn[q_, \[Epsilon]_, \[Kappa]_, \[Tau]_, \[Nu]_, \[Lambda]_, s_, m_, n_] := fn[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, n];
  prefacIn[s_, \[Epsilon]_, \[Tau]_, \[Kappa]_, x_] := (-x)^(-s - I (\[Epsilon] + \[Tau])/2) (1 - x)^(I (\[Epsilon] - \[Tau])/2) E^(I \[Epsilon] \[Kappa] x);
  fUp[q_, \[Epsilon]_, \[Kappa]_, \[Tau]_, \[Nu]_, \[Lambda]_, s_, m_, n_] := (-1)^n Pochhammer[\[Nu] + 1 + s - I \[Epsilon], n]/Pochhammer[\[Nu] + 1 - s + I \[Epsilon], n] fn[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, n];
  prefacUp[s_, \[Epsilon]_, \[Kappa]_, \[Tau]_, \[Nu]_, zhat_] := 2^\[Nu] E^(-\[Pi] \[Epsilon]) E^(-I \[Pi] (\[Nu]+1)) E^(I zhat) zhat^(\[Nu]+I (\[Epsilon]+\[Tau])/2) (zhat-\[Epsilon] \[Kappa])^(-I (\[Epsilon]+\[Tau])/2) E^(-I \[Pi] s) (zhat-\[Epsilon] \[Kappa])^(-s);
];



(* ::Subsection::Closed:: *)
(*Asymptotic amplitudes*)


Switch[MST`$MasterFunction,
"ReggeWheeler",
  prefacInTrans[s_, \[Epsilon]_, \[Tau]_, \[Kappa]_] := E^(I \[Epsilon]);
  prefacUpTrans[s_, \[Epsilon]_, \[Tau]_, \[Kappa]_] := (2I)^s Exp[I \[Epsilon] Log[\[Epsilon]]];
  prefacAplus[s_, \[Epsilon]_, \[Tau]_, \[Kappa]_, \[Nu]_] := (2^(-1 - I \[Epsilon]) E^(-((\[Pi] \[Epsilon])/2) + 1/2 I \[Pi] (1 + \[Nu])) Gamma[\[Nu]+I \[Epsilon]+1])/Gamma[\[Nu]-I \[Epsilon]+1];
  prefacInInc[s_, \[Epsilon]_, \[Tau]_, \[Kappa]_, \[Nu]_, K\[Nu]1_, K\[Nu]2_] := (K\[Nu]1 - I E^(-I \[Pi] \[Nu]) Sin[\[Pi] (\[Nu] + I \[Epsilon])] / Sin[\[Pi] (\[Nu] - I \[Epsilon])] K\[Nu]2) E^(-I*\[Epsilon]*Log[\[Epsilon]]);,
"Teukolsky", 
  prefacInTrans[s_, \[Epsilon]_, \[Tau]_, \[Kappa]_] := 4^s \[Kappa]^(2 s) E^(I (\[Epsilon] + \[Tau]) \[Kappa] (1/2 + Log[\[Kappa]]/(1 + \[Kappa])));
  prefacUpTrans[s_, \[Epsilon]_, \[Tau]_, \[Kappa]_] := (\[Epsilon]/2)^(-1 - 2 s) Exp[I \[Epsilon] (Log[\[Epsilon]] - (1 - \[Kappa])/2)];
  prefacAplus[s_, \[Epsilon]_, \[Tau]_, \[Kappa]_, \[Nu]_] := (2^(-1 + s - I \[Epsilon]) E^(-((\[Pi] \[Epsilon])/2) + 1/2 I \[Pi] (1 - s + \[Nu])) Gamma[1 - s + I \[Epsilon] + \[Nu]])/Gamma[1 + s - I \[Epsilon] + \[Nu]];
  prefacInInc[s_, \[Epsilon]_, \[Tau]_, \[Kappa]_, \[Nu]_, K\[Nu]1_, K\[Nu]2_] := (\[Epsilon]/2)^-1 (K\[Nu]1 - I E^(-I \[Pi] \[Nu]) Sin[\[Pi] (\[Nu] - s + I \[Epsilon])] / Sin[\[Pi] (\[Nu] + s - I \[Epsilon])] K\[Nu]2) Exp[-I \[Epsilon] (Log[\[Epsilon]] - (1 - \[Kappa])/2)];,
_, Abort[];
];


(* ::Subsection::Closed:: *)
(*Radial equation*)


Switch[MST`$MasterFunction,
  "ReggeWheeler",
  d2R[s_, l_, m_, q_, \[Epsilon]_, \[Lambda]_, r_, R_] := -1/(1-2/r)2/r^2 Derivative[1][R][r]+1/(1-2/r)(l (l+1)/r^2+2(1-s^2)/r^3)R[r]-(\[Epsilon]/2)^2/(1-2/r)^2 R[r];,
  "Teukolsky",
  d2R[s_, l_, m_, q_, \[Epsilon]_, \[Lambda]_, r_, R_] := (-(-\[Lambda] + 2 I r s \[Epsilon] + (-2 I (-1 + r) s (-q m + (q^2 + r^2) \[Epsilon]/2) + (-q m + (q^2 + r^2) \[Epsilon]/2)^2)/(q^2 - 2 r + r^2)) R[r] - (-2 + 2 r) (1 + s) Derivative[1][R][r])/(q^2 - 2 r + r^2);,
  _, Abort[]
];


(* ::Section::Closed:: *)
(*Hypergeometric functions*)


(* ::Text:: *)
(*All recurrence relations for the hypergeometric functions below can be derived from equations provided by DLMF*)


(* ::Subsection::Closed:: *)
(*Hypergeometric2F1*)


H2F1Exact[n_, s_, \[Nu]_, \[Tau]_, \[Epsilon]_, x_] :=
 Module[{a = aF[s, \[Nu], \[Tau], \[Epsilon]], b = bF[s, \[Nu], \[Tau], \[Epsilon]], c = cF[s, \[Nu], \[Tau], \[Epsilon]]},
  Hypergeometric2F1[n + a, b-n, c, x]
];

H2F1Up[n_, s_, \[Nu]_, \[Tau]_, \[Epsilon]_, x_] :=
 Module[{a = aF[s, \[Nu], \[Tau], \[Epsilon]], b = bF[s, \[Nu], \[Tau], \[Epsilon]], c = cF[s, \[Nu], \[Tau], \[Epsilon]]},
  1/((3-a+b-2 n) (1+b-c-n) (-1+a+n)){-(1-a+b-2 n) (1+b-n) (-1+a-c+n) H2F1[-2+n], -(-2+a-b+2 n) (-2+2 a-2 b+2 a b+c-a c-b c+4 n-2 a n+2 b n-2 n^2+3 x-4 a x+a^2 x+4 b x-2 a b x+b^2 x-8 n x+4 a n x-4 b n x+4 n^2 x) H2F1[-1+n]}
];

H2F1Down[n_, s_, \[Nu]_, \[Tau]_, \[Epsilon]_, x_] :=
 Module[{a = aF[s, \[Nu], \[Tau], \[Epsilon]], b = bF[s, \[Nu], \[Tau], \[Epsilon]], c = cF[s, \[Nu], \[Tau], \[Epsilon]]},
  1/((-3-a+b-2 n) (-1+b-n) (1+a-c+n)){(-2-a+b-2 n) (-2-2 a+2 b+2 a b+c-a c-b c-4 n-2 a n+2 b n-2 n^2+3 x+4 a x+a^2 x-4 b x-2 a b x+b^2 x+8 n x+4 a n x-4 b n x+4 n^2 x) H2F1[1+n], -(-1-a+b-2 n) (-1+b-c-n) (1+a+n) H2F1[2+n]}
];

dH2F1Exact[n_, s_, \[Nu]_, \[Tau]_, \[Epsilon]_, x_] :=
 Module[{a = aF[s, \[Nu], \[Tau], \[Epsilon]], b = bF[s, \[Nu], \[Tau], \[Epsilon]], c = cF[s, \[Nu], \[Tau], \[Epsilon]]},
  (n+a)(-n+b)/c Hypergeometric2F1[n + a + 1, -n + b + 1, c + 1, x]
];

dH2F1Up[n_, s_, \[Nu]_, \[Tau]_, \[Epsilon]_, x_] :=
 Module[{a = aF[s, \[Nu], \[Tau], \[Epsilon]], b = bF[s, \[Nu], \[Tau], \[Epsilon]], c = cF[s, \[Nu], \[Tau], \[Epsilon]]},
  1/((1+b-c-n) (-1+a+n)){-(((1-a+b-2 n) (1+b-n) (-1+a-c+n) dH2F1[-2+n])/(3-a+b-2 n) ), 1/(3-a+b-2 n)  (2-a+b-2 n) (-2+2 a-2 b+2 a b+c-a c-b c+4 n-2 a n+2 b n-2 n^2+3 x-4 a x+a^2 x+4 b x-2 a b x+b^2 x-8 n x+4 a n x-4 b n x+4 n^2 x) dH2F1[-1+n],(1-a+b-2 n) (2-a+b-2 n) H2F1[-1+n]}
];

dH2F1Down[n_, s_, \[Nu]_, \[Tau]_, \[Epsilon]_, x_] :=
 Module[{a = aF[s, \[Nu], \[Tau], \[Epsilon]], b = bF[s, \[Nu], \[Tau], \[Epsilon]], c = cF[s, \[Nu], \[Tau], \[Epsilon]]},
  1/((-1+b-n) (1+a-c+n)){1/(-3-a+b-2 n)  (-2-a+b-2 n) (-2-2 a+2 b+2 a b+c-a c-b c-4 n-2 a n+2 b n-2 n^2+3 x+4 a x+a^2 x-4 b x-2 a b x+b^2 x+8 n x+4 a n x-4 b n x+4 n^2 x) dH2F1[1+n], -(((-1-a+b-2 n) (-1+b-c-n) (1+a+n) dH2F1[2+n])/(-3-a+b-2 n)),(-2-a+b-2 n) (-1-a+b-2 n) H2F1[1+n]}
];


(* ::Subsection::Closed:: *)
(*HypergeometricU*)


HUExact[n_, s_, \[Nu]_, \[Epsilon]_, zhat_] :=
 Module[{a = aU[s, \[Nu], \[Tau], \[Epsilon]], b = 2 \[Nu] + 2, c = -2 I zhat},
  (c)^n HypergeometricU[n+a,2n+b,c]
];

HUUp[n_, s_, \[Nu]_, \[Epsilon]_, zhat_] :=
 Module[{a = aU[s, \[Nu], \[Tau], \[Epsilon]], b = 2 \[Nu] + 2, c = -2 I zhat},
  1/((-1+a+n) (-4+b+2 n) ) {(-2-a+b+n) (-2+b+2 n) HU[-2+n],(-3+b+2 n) (8+(b+2 n)^2+2 (a+n) c-(b+2 n) (6+c)) HU[-1+n]/c}
];

HUDown[n_, s_, \[Nu]_, \[Epsilon]_, zhat_] :=
 Module[{a = aU[s, \[Nu], \[Tau], \[Epsilon]], b = 2 \[Nu] + 2, c = -2 I zhat},
  1/((-a+b+n) (2+b+2 n)){-(((1+b+2 n) (b^2+4 n (1+n)+b (2+4 n-c)+2 a c) HU[1+n])/ c),(1+a+n) (b+2 n) HU[2+n]}
];

dHUExact[n_, s_, \[Nu]_, \[Epsilon]_, zhat_] :=
 Module[{a = aU[s, \[Nu], \[Tau], \[Epsilon]], b = 2 \[Nu] + 2, c = -2 I zhat},
  (-2 I) (c^(-1+n) n HypergeometricU[a+n,b+2 n,c]-c^n (a+n) HypergeometricU[1+a+n,1+b+2 n,c])
];

dHUUp[n_, s_, \[Nu]_, \[Epsilon]_, zhat_] :=
 Module[{a = aU[s, \[Nu], \[Tau], \[Epsilon]], b = 2 \[Nu] + 2, c = -2 I zhat},
  1/(-1+a+n) {((-2-a+b+n) (-2+b+2 n) dHU[-2+n])/(-4+b+2 n),((-3+b+2 n) (8+b^2+4 (-3+n) n+b (-6+4 n-c)+2 a c) dHU[-1+n])/( (-4+b+2 n) c),(2 I (-3+b+2 n) (-2+b+2 n) HU[-1+n])/c^2}
];

dHUDown[n_, s_, \[Nu]_, \[Epsilon]_, zhat_] :=
 Module[{a = aU[s, \[Nu], \[Tau], \[Epsilon]], b = 2 \[Nu] + 2, c = -2 I zhat},
  1/((a-b-n) (2+b+2 n) c^2){(1+b+2 n) c (b^2+4 n (1+n)+b (2+4 n-c)+2 a c) dHU[1+n],(b+2 n) (-(1+a+n) c^2 dHU[2+n]),(b+2 n)(2 I (1+b+2 n) (2+b+2 n) HU[1+n])}
];


(* ::Section::Closed:: *)
(*MST Series Coefficients*)


(* ::Subsection::Closed:: *)
(*Recurrence formula coefficients*)


(* ::Text:: *)
(*\[Alpha], \[Beta], \[Gamma] defined in Eq. (124) from Sasaki & Tagoshi.*)


\[Alpha][q_, \[Epsilon]_, \[Kappa]_, \[Tau]_, \[Nu]_, \[Lambda]_, s_, m_, n_] :=
 \[Alpha][q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, n] =
  (I \[Epsilon] \[Kappa] (n + \[Nu] + 1 + s + I \[Epsilon]) (n + \[Nu] + 1 + s - I \[Epsilon]) (n + \[Nu] + 1 + I \[Tau]))/((n + \[Nu] + 1) (2 n + 2 \[Nu] + 3));

\[Beta][q_, \[Epsilon]_, \[Kappa]_, \[Tau]_, \[Nu]_, \[Lambda]_, s_, m_, n_] :=
 \[Beta][q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, n] =
  -\[Lambda] - s (s + 1) + (n + \[Nu]) (n + \[Nu] + 1) + \[Epsilon]^2 + \[Epsilon] (\[Epsilon] - m q) + (\[Epsilon] (\[Epsilon] - m q) (s^2 + \[Epsilon]^2))/((n + \[Nu]) (n + \[Nu] + 1));

\[Gamma][q_, \[Epsilon]_, \[Kappa]_, \[Tau]_, \[Nu]_, \[Lambda]_, s_, m_, n_] :=
 \[Gamma][q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, n] =
  -((I \[Epsilon] \[Kappa] (n + \[Nu] - s + I \[Epsilon]) (n + \[Nu] - s - I \[Epsilon]) (n + \[Nu] - I \[Tau]))/((n + \[Nu]) (2 n + 2 \[Nu] - 1)));


(* ::Subsection::Closed:: *)
(*MST series coefficients*)


(* ::Text:: *)
(*fn are the MST coefficients as defined by Sasaki and Tagoshi.*)
(*Sasaki and Tagoshi denote the ingoing MST coefficients as an and the upgoing MST coefficients as fn*)
(*an and fn turn out to be equivalent. We shall therefore only use fn to denote MST coefficients.*)
(*Note: The fn defined in Sasaki and Tagoshi are equivalent to anT, as defined by Casals and Ottewill.*)


fn[q_, \[Epsilon]_, \[Kappa]_, \[Tau]_, \[Nu]_, \[Lambda]_, s_, m_, 0] = 1;

fn[q_, \[Epsilon]_, \[Kappa]_, \[Tau]_, \[Nu]_, \[Lambda]_, s_, m_, nf_] :=
 fn[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, nf] =
 Module[{\[Alpha]n, \[Beta]n, \[Gamma]n, i, n, ret},
  \[Alpha]n[n_] := \[Alpha][q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, n];
  \[Beta]n[n_] := \[Beta][q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, n];
  \[Gamma]n[n_] := \[Gamma][q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, n];

  If[nf > 0,
    ret = fn[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, nf - 1] CF[-\[Alpha]n[i - 1] \[Gamma]n[i], \[Beta]n[i], {i, nf}]/\[Alpha]n[nf - 1];
  ,
    ret = fn[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, nf + 1] CF[-\[Alpha]n[2 nf - i] \[Gamma]n[2 nf - i + 1], \[Beta]n[2 nf - i], {i, nf}]/\[Gamma]n[nf + 1];
  ];
  Clear[\[Alpha]n, \[Beta]n, \[Gamma]n];
  ret
];


(* ::Section::Closed:: *)
(*Coefficients of the Coulomb-series representation*)


(* ::Text:: *)
(*Sasaki & Tagoshi Eqs. (157), (158) and (165) (with r = 0), for the Teukolsky master function. Shared by*)
(*the asymptotic amplitudes and by the large-radius representation of the "In" solution. They must be called*)
(*within an Internal`InheritedBlock[{alpha, beta, gamma, fn}, ...] so that the memoised coefficients stay local.*)


(* sum term[n] from n0 in direction dir until the sum stops changing *)
sumUntil[term_, n0_, dir_] := Module[{res = 0, k = n0}, While[res != (res += term[k]), k += dir]; res];

(* K_nu, ST Eq. (165) with r = 0, CO (3.32) *)
KCoefficient[s_Integer, m_Integer, q_, \[Epsilon]_, \[Kappa]_, \[Tau]_, \[Nu]_, \[Lambda]_] :=
 Module[{\[Epsilon]p = 1/2 (\[Tau] + \[Epsilon])},
  ((2^-\[Nu]) (E^(I \[Epsilon] \[Kappa])) ((\[Epsilon] \[Kappa])^(s - \[Nu])) Gamma[1 - s - 2 I \[Epsilon]p] Gamma[2 + 2 \[Nu]])/(Gamma[1 - s + I \[Epsilon] + \[Nu]] Gamma[1 + s + I \[Epsilon] + \[Nu]] Gamma[1 + \[Nu] + I \[Tau]]) *
   sumUntil[((-1)^# Gamma[1 + # + s + I \[Epsilon] + \[Nu]] Gamma[1 + # + 2 \[Nu]] Gamma[1 + # + \[Nu] + I \[Tau]])/(#! Gamma[1 + # - s - I \[Epsilon] + \[Nu]] Gamma[1 + # + \[Nu] - I \[Tau]]) fn[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, #] &, 0, 1] /
   sumUntil[(((-1)^#) Pochhammer[1 + s - I \[Epsilon] + \[Nu], #])/((-#)! Pochhammer[1 - s + I \[Epsilon] + \[Nu], #] Pochhammer[2 + 2 \[Nu], #]) fn[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, #] &, 0, -1]
 ];

(* A_+^nu, ST Eq. (157), CO (3.38) and (3.41) *)
AplusCoefficient[s_Integer, m_Integer, q_, \[Epsilon]_, \[Kappa]_, \[Tau]_, \[Nu]_, \[Lambda]_] :=
 prefacAplus[s, \[Epsilon], \[Tau], \[Kappa], \[Nu]] (sumUntil[fn[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, #] &, 0, 1] + sumUntil[fn[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, #] &, -1, -1]);

(* A_-^nu, ST Eq. (158), CO (3.19) *)
AminusCoefficient[s_Integer, m_Integer, q_, \[Epsilon]_, \[Kappa]_, \[Tau]_, \[Nu]_, \[Lambda]_] :=
 2^(-s - 1 + I \[Epsilon]) E^(-\[Pi] \[Epsilon] / 2 - I \[Pi] (\[Nu]+1+s) / 2) (sumUntil[(-1)^# Pochhammer[\[Nu] + 1 + s - I \[Epsilon], #]/Pochhammer[\[Nu] + 1 - s + I \[Epsilon], #] fn[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, #] &, 0, 1] + sumUntil[(-1)^# Pochhammer[\[Nu] + 1 + s - I \[Epsilon], #]/Pochhammer[\[Nu] + 1 - s + I \[Epsilon], #] fn[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, #] &, -1, -1]);


(* ::Section::Closed:: *)
(*Asymptotic amplitudes*)


Switch[MST`$MasterFunction,
"ReggeWheeler",
Amplitudes[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, {wp_, prec_, acc_}] :=
 Module[{\[Kappa], \[Tau], \[Epsilon]p, \[Omega], K\[Nu], K\[Nu]1, K\[Nu]2, Aminus, Aplus, D1, D12, D2, D22, InTrans, UpTrans, InInc, UpInc, InRef, UpRef, n, fSumUp, fSumDown, fSumK\[Nu]1Up, fSumK\[Nu]1Down, fSumK\[Nu]2Up, fSumK\[Nu]2Down, fSumAminusUp, fSumAminusDown, fSumD1Up, fSumD1Down, fSumD12Up, fSumD12Down, termf, termK\[Nu]1Up, termK\[Nu]1Down, termK\[Nu]2Up, termK\[Nu]2Down, termAminus, termD1, termD12},
 Internal`InheritedBlock[{\[Alpha], \[Beta], \[Gamma], fn},
  \[Kappa] = Sqrt[1 - q^2];
  \[Tau] = (\[Epsilon] - m q)/\[Kappa];
  \[Epsilon]p = 1/2 (\[Tau] + \[Epsilon]);
  \[Omega] = \[Epsilon] / 2;

  (* All of the formulae are taken from Sasaki & Tagoshi, Living Rev. Relativity 6:6 (ST)
     and Casals & Ottewill, Phys. Rev. D 92, 124055 (CO) *)

  (* There are three formally infinite sums which must be computed, but which may be numerically 
     truncated after a finite number of terms. We determine how many terms to include by summing
     until the result doesn't change. *)

  (* Sum MST series coefficients with no extra factors *)
  termf[n_] := termf[n] = fIn[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, n];
  fSumUp = fSumDown = 0;

  n = 0;
  While[fSumUp != (fSumUp += termf[n]), n++];

  n = -1;
  While[fSumDown != (fSumDown += termf[n]), n--];

  (* Sums appearing in ST Eq. (165) with r=0. We evaluate these with 1: \[Nu] and 2:-\[Nu]-1 *)
  termK\[Nu]1Up[n_] := termK\[Nu]1Up[n] = ((-1)^(2 n) Gamma[-n+s-I \[Epsilon]-\[Nu]] Gamma[1+n-s+I \[Epsilon]+\[Nu]] Gamma[1+n+s+I \[Epsilon]+\[Nu]] Gamma[1+n+2 \[Nu]] Pochhammer[1+I \[Epsilon]+\[Nu],n])/(n! Gamma[1+n-s-I \[Epsilon]+\[Nu]] Gamma[1-s+I \[Epsilon]+\[Nu]] Gamma[1+s+I \[Epsilon]+\[Nu]] Pochhammer[1-I \[Epsilon]+\[Nu],n]) fn[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, n];
  termK\[Nu]1Down[n_] := termK\[Nu]1Down[n] = ((-1)^n Gamma[1+n-I \[Epsilon]+\[Nu]] Gamma[1+n+s-I \[Epsilon]+\[Nu]] Gamma[1+I \[Epsilon]+\[Nu]] Pochhammer[1+I \[Epsilon]+\[Nu],n])/((-n)! Gamma[1+n+I \[Epsilon]+\[Nu]] Gamma[1+n-s+I \[Epsilon]+\[Nu]] Gamma[2+n+2 \[Nu]] Pochhammer[1-I \[Epsilon]+\[Nu],n]) fn[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, n];
  fSumK\[Nu]1Up = fSumK\[Nu]1Down = 0;

  n = 0;
  While[fSumK\[Nu]1Up != (fSumK\[Nu]1Up += termK\[Nu]1Up[n]), n++];

  n = 0;
  While[fSumK\[Nu]1Down != (fSumK\[Nu]1Down += termK\[Nu]1Down[n]), n--];

  termK\[Nu]2Up[n_] := termK\[Nu]2Up[n] = ((-1)^(2 n) Gamma[-n+s-I \[Epsilon]-(-1-\[Nu])] Gamma[1+n-s+I \[Epsilon]+(-1-\[Nu])] Gamma[1+n+s+I \[Epsilon]+(-1-\[Nu])] Gamma[1+n+2 (-1-\[Nu])] Pochhammer[1+I \[Epsilon]+(-1-\[Nu]),n])/(n! Gamma[1+n-s-I \[Epsilon]+(-1-\[Nu])] Gamma[1-s+I \[Epsilon]+(-1-\[Nu])] Gamma[1+s+I \[Epsilon]+(-1-\[Nu])] Pochhammer[1-I \[Epsilon]+(-1-\[Nu]),n])fn[q, \[Epsilon], \[Kappa], \[Tau], (-1-\[Nu]), \[Lambda], s, m, n];
  termK\[Nu]2Down[n_] := termK\[Nu]2Down[n] = ((-1)^n Gamma[1+n-I \[Epsilon]+(-1-\[Nu])] Gamma[1+n+s-I \[Epsilon]+(-1-\[Nu])] Gamma[1+I \[Epsilon]+(-1-\[Nu])] Pochhammer[1+I \[Epsilon]+(-1-\[Nu]),n])/((-n)! Gamma[1+n+I \[Epsilon]+(-1-\[Nu])] Gamma[1+n-s+I \[Epsilon]+(-1-\[Nu])] Gamma[2+n+2 (-1-\[Nu])] Pochhammer[1-I \[Epsilon]+(-1-\[Nu]),n])fn[q, \[Epsilon], \[Kappa], \[Tau], (-1-\[Nu]), \[Lambda], s, m, n];
  fSumK\[Nu]2Up = fSumK\[Nu]2Down = 0;

  n = 0;
  While[fSumK\[Nu]2Up != (fSumK\[Nu]2Up += termK\[Nu]2Up[n]), n++];

  n = 0;
  While[fSumK\[Nu]2Down != (fSumK\[Nu]2Down += termK\[Nu]2Down[n]), n--];

  (* Sum appearing in ST (158), CO (3.19) *)
  termAminus[n_] := termAminus[n] = (-1)^n Pochhammer[\[Nu] + 1 + s - I \[Epsilon], n]/Pochhammer[\[Nu] + 1 - s + I \[Epsilon], n] fn[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, n];
  fSumAminusUp = fSumAminusDown = 0;

  n = 0;
  While[fSumAminusUp != (fSumAminusUp += termAminus[n]), n++];

  n = -1;
  While[fSumAminusDown != (fSumAminusDown += termAminus[n]), n--];

  (* In transmission coefficient: Btrans in ST (167) and CO (3.12) *)
  InTrans = prefacInTrans[s, \[Epsilon], \[Tau], \[Kappa]] (fSumUp+fSumDown);

  (* A-: ST (158), CO (3.19) *)
  Aminus = 2^(-s - 1 + I \[Epsilon]) E^(-\[Pi] \[Epsilon] / 2 - I \[Pi] (\[Nu]+1+s) / 2) (fSumAminusUp+fSumAminusDown);

  (* Up Transmission coefficient: Ctrans in ST (170), CO (3.20) *)
  UpTrans = prefacUpTrans[s, \[Epsilon], \[Tau], \[Kappa]] Aminus;

  (* K\[Nu]: ST (165), CO (3.32) *)
  K\[Nu]1 = ((2^-\[Nu])( E^(I \[Epsilon])) \[Epsilon]^(-1-\[Nu]) Gamma[1-s-2 I \[Epsilon]] Gamma[1+s-I \[Epsilon]+\[Nu]])/(Gamma[-I \[Epsilon]-\[Nu]] Gamma[1+I \[Epsilon]+\[Nu]] Gamma[1-s+I \[Epsilon]+\[Nu]]) fSumK\[Nu]1Up / fSumK\[Nu]1Down;
  K\[Nu]2 = ((2^-(-1-\[Nu]))( E^(I \[Epsilon])) \[Epsilon]^(-1-(-1-\[Nu])) Gamma[1-s-2 I \[Epsilon]] Gamma[1+s-I \[Epsilon]+(-1-\[Nu])])/(Gamma[-I \[Epsilon]-(-1-\[Nu])] Gamma[1+I \[Epsilon]+(-1-\[Nu])] Gamma[1-s+I \[Epsilon]+(-1-\[Nu])]) fSumK\[Nu]2Up / fSumK\[Nu]2Down;

  (* In reflection coefficient: Bref in ST (169), CO (3.37) *)
  InRef = (Gamma[1-2 I \[Epsilon]] Gamma[-I \[Epsilon]-\[Nu]] Gamma[1-I \[Epsilon]+\[Nu]])/(Gamma[1-s-2 I \[Epsilon]] Gamma[s-I \[Epsilon]-\[Nu]] Gamma[1+s-I \[Epsilon]+\[Nu]]) UpTrans (K\[Nu]1 + I E^(I \[Pi] \[Nu]) K\[Nu]2);

  (* A+: ST (157), CO (3.38) and (3.41) *)
  Aplus = prefacAplus[s, \[Epsilon], \[Tau], \[Kappa], \[Nu]] (fSumUp+fSumDown);

  (* In incidence coefficient: Binc from ST (168), CO (3.36) and (3.39) *)
  InInc = (Gamma[1-2 I \[Epsilon]] Gamma[-I \[Epsilon]-\[Nu]] Gamma[1-I \[Epsilon]+\[Nu]])/(Gamma[1-s-2 I \[Epsilon]] Gamma[s-I \[Epsilon]-\[Nu]] Gamma[1+s-I \[Epsilon]+\[Nu]]) prefacInInc[s, \[Epsilon], \[Tau], \[Kappa], \[Nu], K\[Nu]1, K\[Nu]2] Aplus;

  (* Compute Up incidence and reflection coefficients from other coefficients. The incidence follows
     from the Wronskian at the two ends (any frequency); the reflection uses the conjugation symmetry of
     the Regge-Wheeler equation, which holds for real frequencies only, so for a complex frequency it is
     Indeterminate (with a message) rather than wrong. *)
  UpInc = UpTrans/InTrans InInc;
  If[Im[\[Epsilon]] == 0,
    UpRef = -Conjugate[InRef/InTrans] UpTrans;,
    With[{sym = $radialFunctionSymbol}, Message[sym::upref, \[Epsilon]/2]];
    UpRef = Indeterminate;
  ];

  (* Clear local symbols with DownValues to avoid memory leaks *)
  Clear[termf, termK\[Nu]1Up, termK\[Nu]1Down, termK\[Nu]2Up, termK\[Nu]2Down, termAminus, termD1, termD12];

  (* Return results as an Association *)
  <| "In" -> <| "Incidence" -> InInc, "Transmission" -> InTrans, "Reflection" -> InRef|>,
     "Up" -> <| "Incidence" -> UpInc, "Transmission" -> UpTrans, "Reflection" -> UpRef |>
   |>
]],

"Teukolsky",
Amplitudes[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, {wp_, prec_, acc_}] :=
 Module[{\[Kappa], \[Tau], \[Epsilon]p, \[Omega], K\[Nu], K\[Nu]1, K\[Nu]2, Aminus, Aplus, D1, D12, D2, D22, InTrans, UpTrans, InInc, UpInc, InRef, UpRef, n, fSumUp, fSumDown, fSumD1Up, fSumD1Down, fSumD12Up, fSumD12Down, termf, termD1, termD12},
 Internal`InheritedBlock[{\[Alpha], \[Beta], \[Gamma], fn},
  \[Kappa] = Sqrt[1 - q^2];
  \[Tau] = (\[Epsilon] - m q)/\[Kappa];
  \[Epsilon]p = 1/2 (\[Tau] + \[Epsilon]);
  \[Omega] = \[Epsilon] / 2;

  (* At the superradiant bound frequency omega = m Omega_H (\[Epsilon]p = 0) the two horizon solutions
     Delta^-s Exp[+-i k r_*] coincide and the formulae below are singular (1/Sin[2 Pi I \[Epsilon]p] and
     Gamma[1 - s - 2 I \[Epsilon]p] factors). The transmitted flux vanishes there but the transmission
     amplitudes do not: they and, for s <= 0, the "In" amplitudes have finite limits, whereas the "Up"
     horizon coefficients diverge (for s >= 1 so do the "In" amplitudes, the hypergeometric series having
     c = 1 - s - 2 I \[Epsilon]p at a pole). Those limits are taken by the caller from neighbouring
     frequencies (see TeukolskyRadial); here every amplitude is returned as Indeterminate, so that a
     radial function is never silently normalised by a transmission amplitude of 1. *)
  If[\[Epsilon]p == 0,
    Return[<| "In" -> <| "Incidence" -> Indeterminate, "Transmission" -> Indeterminate, "Reflection" -> Indeterminate|>,
              "Up" -> <| "Incidence" -> Indeterminate, "Transmission" -> Indeterminate, "Reflection" -> Indeterminate|>|>];
  ];

  (* All of the formulae are taken from Sasaki & Tagoshi, Living Rev. Relativity 6:6 (ST)
     and Casals & Ottewill, Phys. Rev. D 92, 124055 (CO) *)

  (* There are three formally infinite sums which must be computed, but which may be numerically 
     truncated after a finite number of terms. We determine how many terms to include by summing
     until the result doesn't change. *)

  (* Sum MST series coefficients with no extra factors *)
  termf[n_] := termf[n] = fIn[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, n];
  fSumUp = fSumDown = 0;

  n = 0;
  While[fSumUp != (fSumUp += termf[n]), n++];

  n = -1;
  While[fSumDown != (fSumDown += termf[n]), n--];



  (* Sums appearing in "up" incidence coefficient *)
  termD1[n_] := termD1[n] = (Gamma[1+\[Nu]+n+s+I \[Epsilon]] Gamma[1+\[Nu]+n+I \[Tau]])/(Gamma[1+\[Nu]+n-s-I \[Epsilon]] Gamma[1+\[Nu]+n-I \[Tau]]) fn[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, n];
  fSumD1Up = fSumD1Down = 0;

  n = 0;
  While[fSumD1Up != (fSumD1Up += termD1[n]), n++];

  n = -1;
  While[fSumD1Down != (fSumD1Down += termD1[n]), n--];

  termD12[n_] := termD12[n] = (Gamma[1+(-1-\[Nu])+n+s+I \[Epsilon]] Gamma[1+(-1-\[Nu])+n+I \[Tau]])/(Gamma[1+(-1-\[Nu])+n-s-I \[Epsilon]] Gamma[1+(-1-\[Nu])+n-I \[Tau]]) fn[q, \[Epsilon], \[Kappa], \[Tau], (-1-\[Nu]), \[Lambda], s, m, n];
  fSumD12Up = fSumD12Down = 0;

  n = 0;
  While[fSumD12Up != (fSumD12Up += termD12[n]), n++];

  n = -1;
  While[fSumD12Down != (fSumD12Down += termD12[n]), n--];
  
  (* In transmission coefficient: Btrans in ST (167) and CO (3.12) *)
  InTrans = prefacInTrans[s, \[Epsilon], \[Tau], \[Kappa]] (fSumUp+fSumDown);

  (* A-: ST (158), CO (3.19) *)
  Aminus = AminusCoefficient[s, m, q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda]];

  (* Up Transmission coefficient: Ctrans in ST (170), CO (3.20) *)
  UpTrans = prefacUpTrans[s, \[Epsilon], \[Tau], \[Kappa]] Aminus;

  (* K\[Nu]: ST (165), CO (3.32), for \[Nu] and for -\[Nu]-1 *)
  K\[Nu]1 = KCoefficient[s, m, q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda]];
  K\[Nu]2 = KCoefficient[s, m, q, \[Epsilon], \[Kappa], \[Tau], -1 - \[Nu], \[Lambda]];

  (* In reflection coefficient: Bref in ST (169), CO (3.37) *)
  InRef = UpTrans (K\[Nu]1 + I E^(I \[Pi] \[Nu]) K\[Nu]2);

  (* D2 *)
  D2 = -Exp[(I \[Kappa] (\[Epsilon]+\[Tau]) (1+\[Kappa]+2 Log[\[Kappa]]))/(2 (1+\[Kappa]))] (2\[Kappa])^(2 s) ( Sin[\[Pi] (\[Nu]-I \[Epsilon])] Sin[\[Pi] (\[Nu]-I \[Tau])])/(Sin[2 \[Pi] \[Nu]] Sin[\[Pi] I (\[Epsilon]+\[Tau])]) (fSumUp+fSumDown);
  D22 = -Exp[(I \[Kappa] (\[Epsilon]+\[Tau]) (1+\[Kappa]+2 Log[\[Kappa]]))/(2 (1+\[Kappa]))] (2\[Kappa])^(2 s) ( Sin[\[Pi] ((-1-\[Nu])-I \[Epsilon])] Sin[\[Pi] ((-1-\[Nu])-I \[Tau])])/(Sin[2 \[Pi] (-1-\[Nu])] Sin[\[Pi] I (\[Epsilon]+\[Tau])]) (fSumUp+fSumDown);

  (* Up reflection coefficient *)
  UpRef = Exp[-\[Pi] \[Epsilon]-I \[Pi] s]/Sin[2\[Pi] \[Nu]] ((Exp[-I \[Pi] \[Nu]]Sin[\[Pi](\[Nu]-s+I \[Epsilon])])/K\[Nu]1 D2-I Sin[\[Pi](\[Nu]+s-I \[Epsilon])]/K\[Nu]2 D22);

  (* A+: ST (157), CO (3.38) and (3.41) *)
  Aplus = AplusCoefficient[s, m, q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda]];

  (* In incidence coefficient: Binc from ST (168), CO (3.36) and (3.39) *)
  InInc = prefacInInc[s, \[Epsilon], \[Tau], \[Kappa], \[Nu], K\[Nu]1, K\[Nu]2] Aplus;

  (* D1 *)
  D1 = Exp[-((I \[Kappa] (\[Epsilon]+\[Tau]) (1+\[Kappa]+2 Log[\[Kappa]]))/(2 (1+\[Kappa])))] ( Sin[\[Pi] (\[Nu]+I \[Epsilon])] Sin[\[Pi] (\[Nu]+I \[Tau])] Gamma[1-s-I (\[Epsilon]+\[Tau])])/(Sin[2 \[Pi] \[Nu]] Sin[\[Pi] I (\[Epsilon]+\[Tau])]  Gamma[1+s+I \[Epsilon]+I \[Tau]]) (fSumD1Up+fSumD1Down);
  D12 = Exp[-((I \[Kappa] (\[Epsilon]+\[Tau]) (1+\[Kappa]+2 Log[\[Kappa]]))/(2 (1+\[Kappa])))] ( Sin[\[Pi] ((-1-\[Nu])+I \[Epsilon])] Sin[\[Pi] ((-1-\[Nu])+I \[Tau])] Gamma[1-s-I (\[Epsilon]+\[Tau])])/(Sin[2 \[Pi] (-1-\[Nu])] Sin[\[Pi] I (\[Epsilon]+\[Tau])]  Gamma[1+s+I \[Epsilon]+I \[Tau]]) (fSumD12Up+fSumD12Down);

  (* Up incidence coefficient *)
  UpInc = Exp[-\[Pi] \[Epsilon]-I \[Pi] s]/Sin[2\[Pi] \[Nu]] ((Exp[-I \[Pi] \[Nu]]Sin[\[Pi](\[Nu]-s+I \[Epsilon])])/K\[Nu]1 D1-I Sin[\[Pi](\[Nu]+s-I \[Epsilon])]/K\[Nu]2 D12);

  (* Clear local symbols with DownValues to avoid memory leaks *)
  Clear[termf, termD1, termD12];

  (* Return results as an Association *)
  <| "In" -> <| "Incidence" -> InInc, "Transmission" -> InTrans, "Reflection" -> InRef|>,
     "Up" -> <| "Incidence" -> UpInc, "Transmission" -> UpTrans, "Reflection" -> UpRef |>
   |>
]]
];


(* ::Section::Closed:: *)
(*Radial "In" solution*)


(* ::Text:: *)
(*Ingoing MST Teukolsky Radial Function: Throwe B.1, Sasaki & Tagoshi Eqs. (116) and (120)*)
(*Ingoing MST Regge Wheeler Radial Function: Casals & Ottewill Eqs. (3.1) and (3.4a)*)
(*Ingoing Regge Wheeler MST coefficients: Casals & Ottewill Eqs. (3.8) & (3.9)*)


SetAttributes[MSTRadialIn, {NumericFunction}];


(* ::Subsection::Closed:: *)
(*Radial function*)


(* Sum term[n] from n0 in direction dir until the partial sum stops changing or the term falls below the
   goals; term[n] is a number, or a list {value, derivative} for a combined evaluation, in which case every
   component must have converged. *)
sumSeries[term_, n0_, dir_, prec_, acc_] :=
 Module[{res = 0, old, t, n = n0},
  While[True,
    t = term[n]; old = res; res = old + t;
    If[!(Or @@ Thread[Flatten[{res}] != Flatten[{old}]] && Or @@ Thread[Abs[Flatten[{t}]] > 10^-acc + Abs[Flatten[{res}]] 10^-prec]), Break[]];
    n += dir;
  ];
  res
 ];

(* The hypergeometric series for the "In" solution (Sasaki & Tagoshi Eq. (116)) at r: the value (deriv 0),
   the first derivative (deriv 1), or both from a single summation (deriv All), which shares the
   coefficients and the hypergeometric functions between the two. *)
mstRadialInSeriesCore[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_, {wp_, prec_, acc_}, deriv_][r_?NumericQ] :=
 Module[{\[Kappa], \[Tau], rp, x, dxdr, prefac, dprefac, term, res},
 Block[{H2F1, dH2F1},
 Internal`InheritedBlock[{\[Alpha], \[Beta], \[Gamma], fn},
  \[Kappa] = Sqrt[1 - q^2];
  \[Tau] = (\[Epsilon] - m q)/\[Kappa];
  rp = 1 + \[Kappa];
  x = (rp - r)/(2 \[Kappa]);
  dxdr = - 1/(2\[Kappa]);

  H2F1[n : (0 | 1)] := H2F1[n] = H2F1Exact[n, s, \[Nu], \[Tau], \[Epsilon], x];

  H2F1[n_Integer] := H2F1[n] =
   Module[{t1, t2, res},
    {t1, t2} = If[n>0, H2F1Up[n, s, \[Nu], \[Tau], \[Epsilon], x], H2F1Down[n, s, \[Nu], \[Tau], \[Epsilon], x]];
    res = t1 + t2;
    If[Max[Abs[{t1, t2}/res]] > 2.,
      res = H2F1Exact[n, s, \[Nu], \[Tau], \[Epsilon], x];
    ];
    res
  ];

  dH2F1[n : (0 | 1)] := dH2F1[n] = dH2F1Exact[n, s, \[Nu], \[Tau], \[Epsilon], x];

  dH2F1[n_Integer] := dH2F1[n] =
   Module[{t1, t2, t3, res},
    {t1, t2, t3} = If[n>0, dH2F1Up[n, s, \[Nu], \[Tau], \[Epsilon], x], dH2F1Down[n, s, \[Nu], \[Tau], \[Epsilon], x]];
    res = t1 + t2 + t3;
    If[Max[Abs[{t1, t2, t3}/res]] > 2.,
      res = dH2F1Exact[n, s, \[Nu], \[Tau], \[Epsilon], x];
    ];
    res
  ];
 
  prefac = prefacIn[s, \[Epsilon], \[Tau], \[Kappa], x]/norm;
  If[deriv =!= 0, dprefac = Derivative[0,0,0,0,1][prefacIn][s, \[Epsilon], \[Tau], \[Kappa], x]/norm];

  Switch[deriv,
    0,   term[n_] := term[n] = fIn[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, n] prefac H2F1[n];,
    1,   term[n_] := term[n] = fIn[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, n] dxdr (dprefac H2F1[n] + prefac dH2F1[n]);,
    All, term[n_] := term[n] = fIn[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, n] {prefac H2F1[n], dxdr (dprefac H2F1[n] + prefac dH2F1[n])};
  ];
  res = sumSeries[term, 0, 1, prec, acc] + sumSeries[term, -1, -1, prec, acc];
  Clear[term];
  res
]]];

mstRadialInSeries[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_, {wp_, prec_, acc_}][r_?NumericQ] :=
  mstRadialInSeriesCore[s, l, m, q, \[Epsilon], \[Nu], \[Lambda], norm, {wp, prec, acc}, 0][r];

Derivative[1][mstRadialInSeries[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_, {wp_, prec_, acc_}]][r_?NumericQ] :=
  mstRadialInSeriesCore[s, l, m, q, \[Epsilon], \[Nu], \[Lambda], norm, {wp, prec, acc}, 1][r];


(* ::Section::Closed:: *)
(*Radial "Up" solution*)


(* ::Text:: *)
(*Upgoing MST Teukolsky Radial Function: Throwe B.5, Sasaki & Tagoshi Eqs. (153) and (159)*)
(*Upgoing MST Regge Wheeler Radial Function: Casals & Ottewill Eq. (3.15)*)
(*Regge Wheeler MST coefficients: Casals & Ottewill Eqs. (3.8) & (3.9)*)


SetAttributes[MSTRadialUp, {NumericFunction}];


(* ::Subsection::Closed:: *)
(*Radial function*)


(* The series of Tricomi functions for the "Up" solution (Sasaki & Tagoshi Eqs. (153), (159)) at r: the
   value (deriv 0), the first derivative (deriv 1), or both from a single summation (deriv All). *)
mstRadialUpSeriesCore[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_, {wp_, prec_, acc_}, deriv_][r_?NumericQ] :=
 Module[{\[Kappa], \[Tau], \[Epsilon]p, rm, z, zm, zhat, dzhatdr, prefac, dprefac, term, res},
 Block[{HU, dHU},
 Internal`InheritedBlock[{\[Alpha], \[Beta], \[Gamma], fn},
  \[Kappa] = Sqrt[1 - q^2];
  \[Tau] = (\[Epsilon] - m q)/\[Kappa];
  \[Epsilon]p = 1/2 (\[Tau]+\[Epsilon]);
  rm = 1 - Sqrt[1 - q^2];
  z = \[Epsilon] r / 2;
  zm = \[Epsilon] rm / 2;
  zhat = z - zm;
  dzhatdr = \[Epsilon] / 2;
 
  HU[n : (0 | 1)] := HU[n] = HUExact[n, s, \[Nu], \[Epsilon], zhat];
 
  HU[n_Integer] := HU[n] =
   Module[{t1, t2, res},
    {t1, t2} = If[n>0, HUUp[n, s, \[Nu], \[Epsilon], zhat], HUDown[n, s, \[Nu], \[Epsilon], zhat]];
    res = t1 + t2;
    If[Max[Abs[{t1, t2}/res]] > 2.,
      res = HUExact[n, s, \[Nu], \[Epsilon], zhat];
    ];
    res
  ];
 
  dHU[n : (0 | 1)] := dHU[n] = dHUExact[n, s, \[Nu], \[Epsilon], zhat];
  
  dHU[n_Integer] := dHU[n] =
   Module[{t1, t2, t3, res},
    {t1, t2, t3} = If[n>0, dHUUp[n, s, \[Nu], \[Epsilon], zhat], dHUDown[n, s, \[Nu], \[Epsilon], zhat]];
    res = t1 + t2 + t3 ;
    If[Max[Abs[{t1, t2, t3}/res]] > 2.,
      res = dHUExact[n, s, \[Nu], \[Epsilon], zhat];
    ];
    res
  ];

  prefac = prefacUp[s, \[Epsilon], \[Kappa], \[Tau], \[Nu], zhat]/norm;
  If[deriv =!= 0, dprefac = Derivative[0,0,0,0,0,1][prefacUp][s, \[Epsilon], \[Kappa], \[Tau], \[Nu], zhat]/norm];

  Switch[deriv,
    0,   term[n_] := term[n] = fUp[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, n] prefac HU[n];,
    1,   term[n_] := term[n] = fUp[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, n] dzhatdr (dprefac HU[n] + prefac dHU[n]);,
    All, term[n_] := term[n] = fUp[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, n] {prefac HU[n], dzhatdr (dprefac HU[n] + prefac dHU[n])};
  ];
  res = sumSeries[term, 0, 1, prec, acc] + sumSeries[term, -1, -1, prec, acc];
  Clear[term];
  res
]]];

mstRadialUpSeries[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_, {wp_, prec_, acc_}][r_?NumericQ] :=
  mstRadialUpSeriesCore[s, l, m, q, \[Epsilon], \[Nu], \[Lambda], norm, {wp, prec, acc}, 0][r];

Derivative[1][mstRadialUpSeries[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_, {wp_, prec_, acc_}]][r_?NumericQ] :=
  mstRadialUpSeriesCore[s, l, m, q, \[Epsilon], \[Nu], \[Lambda], norm, {wp, prec, acc}, 1][r];


(* ::Section::Closed:: *)
(*Radial "In" solution at large radius: series of Coulomb wave functions*)


(* ::Text:: *)
(*At large omega (r - r+) the hypergeometric series for the "In" solution suffers catastrophic cancellation*)
(*(roughly 0.8 omega (r - r+) digits are lost), whereas the series of Coulomb wave functions do not degrade*)
(*with radius. Sasaki & Tagoshi Eq. (166) gives R_in = K_nu R_C^nu + K_{-nu-1} R_C^{-nu-1}, and by their*)
(*Eq. (152) R_C^nu = R_+^nu + R_-^nu, with R_-^nu the outgoing series of Tricomi functions used for the "Up"*)
(*solution (Eqs. (153), (159)) and R_+^nu its incoming partner. Since R_+^{-nu-1} and R_+^nu are both purely*)
(*incoming at infinity they are proportional, with ratio A_+^{-nu-1}/A_+^nu (Eqs. (155), (157)), and likewise*)
(*for R_-, so that R_in = [K_nu + K_{-nu-1} A_+^{-nu-1}/A_+^nu] R_+^nu + [K_nu + K_{-nu-1} A_-^{-nu-1}/A_-^nu] R_-^nu.*)
(*Only available for the Teukolsky master function.*)


(* R_+^nu of ST Eq. (153) and its r-derivative, up to an overall sign: the series of Tricomi functions
   U(n+nu+1-s+i eps, 2n+2nu+2, +2i zhat), obtained from the "Up" series (U(n+nu+1+s-i eps, 2n+2nu+2, -2i zhat),
   with its recurrences) by the substitution (s, eps, zhat) -> (-s, -eps, -zhat), times the factors relating
   the two (DLMF 33.2.7: R_+ is built on the Coulomb function H^-, R_- on H^+). *)
(* zz is the formal variable of the prefactor, differentiated symbolically below; a package symbol rather
   than a Module local, which a message issued during the evaluation could keep alive as a leaked symbol *)
mstRadialPlusSeries[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, {wp_, prec_, acc_}, deriv_][r_?NumericQ] :=
 Module[{\[Kappa], \[Tau], rm, zhat, \[Eta], Q, dQ, G, term, res},
 Block[{HU, dHU},
 Internal`InheritedBlock[{\[Alpha], \[Beta], \[Gamma], fn},
  \[Kappa] = Sqrt[1 - q^2];
  \[Tau] = (\[Epsilon] - m q)/\[Kappa];
  rm = 1 - \[Kappa];
  zhat = \[Epsilon] (r - rm)/2;
  \[Eta] = -I s - \[Epsilon];
  (* n-independent prefactor: prefacUp times the factors relating the H^- series to the H^+ series;
     equal to minus the prefactor of ST Eq. (153) *)
  Q = prefacUp[s, \[Epsilon], \[Kappa], \[Tau], \[Nu], zz] Exp[-2 I (zz - \[Eta] Log[2 zz] - \[Nu] \[Pi]/2)] (-2 I zz)^(-\[Nu] - 1 - s + I \[Epsilon]) (2 I zz)^(\[Nu] + 1 - s + I \[Epsilon]);
  If[deriv =!= 0, dQ = D[Q, zz] /. zz -> zhat];
  Q = Q /. zz -> zhat;
  G[n_] := Gamma[n + \[Nu] + 1 - s + I \[Epsilon]]/Gamma[n + \[Nu] + 1 + s - I \[Epsilon]];

  (* (2 i zhat)^n U(n+nu+1-s+i eps, 2n+2nu+2, 2 i zhat) via the "Up" recurrences with flipped signs *)
  HU[n : (0 | 1)] := HU[n] = HUExact[n, -s, \[Nu], -\[Epsilon], -zhat];
  HU[n_Integer] := HU[n] =
   Module[{t1, t2, res},
    {t1, t2} = If[n > 0, HUUp[n, -s, \[Nu], -\[Epsilon], -zhat], HUDown[n, -s, \[Nu], -\[Epsilon], -zhat]];
    res = t1 + t2;
    If[Max[Abs[{t1, t2}/res]] > 2., res = HUExact[n, -s, \[Nu], -\[Epsilon], -zhat]];
    res
  ];
  (* derivative with respect to the recurrences' own argument, -zhat (the zhat derivative is its negative) *)
  dHU[n : (0 | 1)] := dHU[n] = dHUExact[n, -s, \[Nu], -\[Epsilon], -zhat];
  dHU[n_Integer] := dHU[n] =
   Module[{t1, t2, t3, res},
    {t1, t2, t3} = If[n > 0, dHUUp[n, -s, \[Nu], -\[Epsilon], -zhat], dHUDown[n, -s, \[Nu], -\[Epsilon], -zhat]];
    res = t1 + t2 + t3;
    If[Max[Abs[{t1, t2, t3}/res]] > 2., res = dHUExact[n, -s, \[Nu], -\[Epsilon], -zhat]];
    res
  ];

  Switch[deriv,
    0,   term[n_] := term[n] = (-1)^n fUp[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, n] G[n] Q HU[n];,
    1,   term[n_] := term[n] = (-1)^n fUp[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, n] G[n] (dQ HU[n] - Q dHU[n]) \[Epsilon]/2;,
    All, term[n_] := term[n] = (-1)^n fUp[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, n] G[n] {Q HU[n], (dQ HU[n] - Q dHU[n]) \[Epsilon]/2};
  ];
  res = sumSeries[term, 0, 1, prec, acc] + sumSeries[term, -1, -1, prec, acc];
  Clear[term];
  -res
]]];


(* Coefficients {c+, c-} of R_in = c+ R_+^nu + c- R_-^nu (ST Eqs. (152), (155)-(158), (165), (166)) *)
$inConnectionCache = <||>;   (* bounded cache of the connection coefficients, keyed by the parameters *)

inConnectionCoefficients[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_] :=
 Module[{key = {s, l, m, q, \[Epsilon], \[Nu], \[Lambda]}, res},
  res = Lookup[$inConnectionCache, Key[key], None];
  If[res =!= None, Return[res]];
  res = inConnectionCoefficientsCompute[s, l, m, q, \[Epsilon], \[Nu], \[Lambda]];
  If[Length[$inConnectionCache] >= 50, $inConnectionCache = <||>];
  $inConnectionCache[key] = res
 ];

inConnectionCoefficientsCompute[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_] :=
 Module[{\[Kappa], \[Tau], K1, K2},
 Internal`InheritedBlock[{\[Alpha], \[Beta], \[Gamma], fn},
  \[Kappa] = Sqrt[1 - q^2];
  \[Tau] = (\[Epsilon] - m q)/\[Kappa];
  K1 = KCoefficient[s, m, q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda]];
  K2 = KCoefficient[s, m, q, \[Epsilon], \[Kappa], \[Tau], -1 - \[Nu], \[Lambda]];
  {K1 + K2 AplusCoefficient[s, m, q, \[Epsilon], \[Kappa], \[Tau], -1 - \[Nu], \[Lambda]]/AplusCoefficient[s, m, q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda]],
   K1 + K2 AminusCoefficient[s, m, q, \[Epsilon], \[Kappa], \[Tau], -1 - \[Nu], \[Lambda]]/AminusCoefficient[s, m, q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda]]}
]];


mstRadialInCoulomb[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_, {wp_, prec_, acc_}, deriv_][r_?NumericQ] :=
 Module[{cp, cm},
  {cp, cm} = inConnectionCoefficients[s, l, m, q, \[Epsilon], \[Nu], \[Lambda]];
  (cp mstRadialPlusSeries[s, l, m, q, \[Epsilon], \[Nu], \[Lambda], {wp, prec, acc}, deriv][r]
   + cm mstRadialUpSeriesCore[s, l, m, q, \[Epsilon], \[Nu], \[Lambda], 1, {wp, prec, acc}, deriv][r])/norm
 ];


(* ::Section::Closed:: *)
(*Radial "In" solution at large radius: hypergeometric series in 1/(1-x)*)


(* ::Text:: *)
(*Sasaki & Tagoshi Eqs. (137) and (138): R_in = R_0^nu + R_0^{-nu-1}, with R_0^nu a series of hypergeometric*)
(*functions of argument 1/(1-x) = eps kappa/zhat, which tends to zero at large radius, so that this*)
(*representation converges best where the series in x (Eq. (120)) is worst. Same normalisation as Eq. (116).*)


mstRadialInLargeRadiusSeries[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_, {wp_, prec_, acc_}, deriv_][r_?NumericQ] :=
 Module[{\[Kappa], \[Tau], \[Epsilon]p, rp, x, xx, R0, res},
 Internal`InheritedBlock[{\[Alpha], \[Beta], \[Gamma], fn},
  \[Kappa] = Sqrt[1 - q^2];
  \[Tau] = (\[Epsilon] - m q)/\[Kappa];
  \[Epsilon]p = 1/2 (\[Tau] + \[Epsilon]);
  rp = 1 + \[Kappa];
  x = (rp - r)/(2 \[Kappa]);
  (* R_0^nu of ST Eq. (138), or its r-derivative, as a function of the symbolic xx (dx/dr = -1/(2 kappa)) *)
  R0[nu_] := Module[{pref, coef, g, term, res},
    pref = E^(I \[Epsilon] \[Kappa] xx) (-xx)^(-s - I \[Epsilon]p) (1 - xx)^(I \[Epsilon]p + nu);
    coef[n_] := Gamma[1 - s - I \[Epsilon] - I \[Tau]] Gamma[2 n + 2 nu + 1]/(Gamma[n + nu + 1 - I \[Tau]] Gamma[n + nu + 1 - s - I \[Epsilon]]) fn[q, \[Epsilon], \[Kappa], \[Tau], nu, \[Lambda], s, m, n];
    g[n_] := (1 - xx)^n Hypergeometric2F1[-n - nu - I \[Tau], -n - nu - s - I \[Epsilon], -2 n - 2 nu, 1/(1 - xx)];
    Switch[deriv,
      0,   term[n_] := term[n] = coef[n] (pref g[n] /. xx -> x);,
      1,   term[n_] := term[n] = coef[n] (-1/(2 \[Kappa])) (D[pref g[n], xx] /. xx -> x);,
      All, term[n_] := term[n] = coef[n] {pref g[n] /. xx -> x, (-1/(2 \[Kappa])) (D[pref g[n], xx] /. xx -> x)};
    ];
    res = sumSeries[term, 0, 1, prec, acc] + sumSeries[term, -1, -1, prec, acc];
    Clear[term];
    res
  ];
  (R0[\[Nu]] + R0[-1 - \[Nu]])/norm
]];


(* ::Section::Closed:: *)
(*Evaluation with precision padding*)


(* ::Text:: *)
(*The MST series lose digits to cancellation: the "Up" (Coulomb) series lose an amount that grows with omega,*)
(*the hypergeometric "In" series additionally about 0.8 omega (r - r+) digits. The loss is arithmetic, not a*)
(*sensitivity to the inputs, so it is recovered by summing the series with the inputs padded to a higher working*)
(*precision. The deficit is measured from the tracked precision of a first evaluation (for machine-precision*)
(*input the first evaluation is done in arbitrary precision a few digits above machine precision) and the series*)
(*are re-summed with the inputs padded by the deficit; the result is returned at the precision of the input.*)


$MSTRepresentationThreshold = 2;
(* representation of the "In" solution beyond the threshold: "Coulomb" (ST Eq. (166), series of Coulomb wave
   functions) or "Hypergeometric" (ST Eq. (138), series of hypergeometric functions in 1/(1-x)) *)
$MSTInLargeRadiusRepresentation = "Coulomb";
$masterFunction = MST`$MasterFunction;   (* captured at load time; MST`$MasterFunction is only set while the package loads *)
$radialFunctionSymbol = Symbol[MST`$MasterFunction <> "`" <> MST`$MasterFunction <> "RadialFunction"];   (* carries the messages *)

With[{sym = $radialFunctionSymbol},
  sym::prec = "The MST series for the `1` radial function at r = `2` could only be evaluated to a precision of `3` (`4` requested).";
  sym::upref = "The \"Up\" reflection amplitude is only available for real frequencies (\[Omega] = `1` given) and is Indeterminate.";
];

(* The large-radius representations were derived and validated for real frequencies; for a complex
   frequency the hypergeometric series is used at every radius. *)
mstInRepresentation[q_, \[Epsilon]_, r_] :=
 If[$masterFunction === "Teukolsky" && Im[\[Epsilon]] == 0 && Abs[\[Epsilon]] (r - (1 + Sqrt[1 - q^2]))/2 > $MSTRepresentationThreshold, $MSTInLargeRadiusRepresentation, "Series"];

mstInCore[rep_] := Switch[rep, "Coulomb", mstRadialInCoulomb, "Hypergeometric", mstRadialInLargeRadiusSeries, _, mstRadialInSeriesCore];

(* The eigenvalue and the renormalized angular momentum at the precision pp of a padded evaluation. The MST
   series amplify an error in nu (or lambda) by roughly the number of digits they lose to cancellation, so
   padding the digits of a nu known to fewer digits than pp (which SetPrecision would do) gives a wrong
   result whose tracked precision does not show it. Instead the eigenvalue is recomputed at pp and nu is
   recomputed from it, padded until it carries pp digits, on the same representative as the given nu. The
   results are cached per parameter set and precision. For a master function other than Teukolsky the
   eigenvalue is used as given. *)
$refinedParameterCache = <||>;

refinedParameters[s_, l_, m_, q_, \[Epsilon]_, \[Lambda]_, \[Nu]_, pp_] :=
 Module[{key = {s, l, m, q, \[Epsilon], \[Lambda], \[Nu], pp}, \[Lambda]p = \[Lambda], \[Nu]p = \[Nu], res},
  If[Precision[\[Lambda]] >= pp && Precision[\[Nu]] >= pp, Return[SetPrecision[{\[Lambda], \[Nu]}, pp]]];
  res = Lookup[$refinedParameterCache, Key[key], None];
  If[res =!= None, Return[res]];
  If[$masterFunction === "Teukolsky" && Precision[\[Lambda]] < pp,
    (* the eigenvalue code compares against a machine-number tolerance, which underflows at a few
       hundred digits with a harmless General::munfl *)
    \[Lambda]p = Quiet[SpinWeightedSpheroidalEigenvalue[s, l, m, SetPrecision[q \[Epsilon]/2, pp]], General::munfl];
  ];
  If[Precision[\[Nu]] < pp,
    \[Nu]p = paddedNu[s, l, m, q, \[Epsilon], \[Lambda]p, \[Nu], pp];
  ];
  res = SetPrecision[{\[Lambda]p, \[Nu]p}, pp];
  If[Length[$refinedParameterCache] >= 50, $refinedParameterCache = <||>];
  $refinedParameterCache[key] = res
 ];

(* nu to pp digits, from the renormalized angular momentum computed with padded inputs, on the
   representative (among +-nu + k) of the given nu *)
paddedNu[s_, l_, m_, q_, \[Epsilon]_, \[Lambda]_, \[Nu]_, pp_] :=
 Module[{p = pp, res, tries = 0, ramAt},
  ramAt[p1_] := RenormalizedAngularMomentum[s, l, m, SetPrecision[q, p1], SetPrecision[\[Epsilon]/2, p1], SetPrecision[\[Lambda], p1]];
  res = ramAt[p];
  While[tries < 4 && (!NumericQ[res] || Precision[res] < pp - 1),
    p = If[NumericQ[res], p + Ceiling[pp - Precision[res]] + 3, 2 p];
    res = ramAt[p];
    tries++;
  ];
  If[!NumericQ[res], Return[\[Nu]]];
  nearestRepresentative[res, \[Nu]]
 ];

(* the member of {+-nu + k, k integer} closest to nu0 *)
nearestRepresentative[\[Nu]_, \[Nu]0_] :=
 Module[{cands},
  cands = Flatten[{\[Nu] + Round[Re[\[Nu]0 - \[Nu]]], -\[Nu] + Round[Re[\[Nu]0 + \[Nu]]]}];
  First[MinimalBy[cands, Abs[# - \[Nu]0] &]]
 ];

(* Evaluate core (a series definition taking a derivative order) at r with the working precision padded
   until the result carries the precision of the input. The loss is arithmetic cancellation, so it is
   measured from the tracked precision of a first evaluation and the series are re-summed with the inputs
   padded by the deficit; a non-numeric result (a precision-zero intermediate, 1/0) is retried at twice
   the precision. The eigenvalue and nu are refined to the padded precision (see refinedParameters).
   $lastPaddingLoss records the digits lost at the last working precision used, $lastPaddingPrecision
   that precision; maxTries limits the retries and p0 sets the initial working precision. *)
$lastPaddingLoss = 0;
$lastPaddingPrecision = 0;

(* Extra working precision for a mode {s, l, m, q, epsilon}, set by the master package when an independent
   consistency check (the Wronskian) shows that the padded evaluations are wrong although their tracked
   precision looks fine: for some modes the recurrence for the MST coefficients passes through a region
   where its two solutions are hard to separate, and below a certain working precision it silently yields
   the wrong one, which no precision tracking reveals. *)
$modePadding = <||>;
modePadding[s_, l_, m_, q_, \[Epsilon]_] := Lookup[$modePadding, Key[{s, l, m, q, \[Epsilon]}], 0];
setModePadding[s_, l_, m_, q_, \[Epsilon]_, extra_] := (If[Length[$modePadding] >= 50, $modePadding = <||>]; $modePadding[{s, l, m, q, \[Epsilon]}] = extra);

allNumericQ[x_] := VectorQ[Flatten[{x}], NumericQ];   (* a number, or a list of numbers ({value, derivative}) *)

mstPaddedEvaluation[core_, {s_, l_, m_, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_}, {wp_, prec_, acc_}, deriv_, r_, maxTries_:4, p0_:Automatic] :=
 Module[{target, p, res, deficit, tries = 0, eval, numericQ},
  numericQ[x_] := allNumericQ[x] && !AllTrue[Flatten[{x}], # == 0 &];
  If[wp === MachinePrecision,
    target = $MachinePrecision; p = $MachinePrecision + 4;,
    target = wp; p = wp;
  ];
  p += modePadding[s, l, m, q, \[Epsilon]];
  If[NumericQ[p0], p = Max[p, p0]];
  p = Ceiling[p];   (* an integer, so that 10^-prec below stays exact (10^-318. would underflow) *)
  eval[pp_] := Module[{params = SetPrecision[{q, \[Epsilon], norm}, pp], \[Lambda]p, \[Nu]p, rr = SetPrecision[r, pp], f, precgoal},
    {\[Lambda]p, \[Nu]p} = refinedParameters[s, l, m, q, \[Epsilon], \[Lambda], \[Nu], pp];
    precgoal = If[pp > target, pp - 2, prec];
    f = core[s, l, m, params[[1]], params[[2]], \[Nu]p, \[Lambda]p, params[[3]], {pp, precgoal, acc}, deriv];
    (* a precision-zero intermediate at too low a working precision is retried below, not reported *)
    Quiet[f[rr], {Power::infy, Infinity::indet, Divide::infy}]
  ];
  res = eval[p];
  While[tries < maxTries && (!numericQ[res] || (deficit = target - Precision[res]) > 1),
    p = If[numericQ[res], p + Ceiling[deficit] + 3, 2 p];
    res = eval[p];
    tries++;
  ];
  $lastPaddingPrecision = p;
  $lastPaddingLoss = If[numericQ[res], p - Precision[res], Infinity];
  If[maxTries > 1 && (!allNumericQ[res] || (numericQ[res] && Precision[res] < target - 1)),
    With[{sym = $radialFunctionSymbol}, Message[sym::prec, If[core === mstRadialUpSeriesCore, "Up", "In"], r, If[allNumericQ[res], Precision[res], res], target]];
  ];
  Which[
    !allNumericQ[res], res,
    wp === MachinePrecision, N[res],
    Precision[res] > target, SetPrecision[res, target],
    True, res
  ]
];

(* "In" solution: the representation is chosen by the radius (mstInRepresentation), except that a
   representation which loses many more digits than the hypergeometric series would (a mode deep under
   its potential barrier, where the Coulomb-type solutions are exponentially larger than the "In"
   solution) is abandoned for the series: the loss of the first evaluation of a mode is measured and
   cached, and compared with the roughly 0.8 omega (r - r+) digits that the series loses. *)
$inRepresentationCache = <||>;

mstRadialInEvaluate[params:{s_, l_, m_, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_}, goals_, deriv_, r_] :=
 Module[{rep = mstInRepresentation[q, \[Epsilon], r], key = {s, l, m, q, \[Epsilon], \[Nu], \[Lambda]}, seriesLoss, cached, res},
  If[rep === "Series", Return[mstPaddedEvaluation[mstRadialInSeriesCore, params, goals, deriv, r]]];
  seriesLoss = 0.8 Abs[\[Epsilon]] (r - (1 + Sqrt[1 - q^2]))/2 + 8;
  cached = Lookup[$inRepresentationCache, Key[key], None];
  If[cached =!= None && cached > seriesLoss, Return[mstPaddedEvaluation[mstRadialInSeriesCore, params, goals, deriv, r]]];
  If[cached =!= None, Return[mstPaddedEvaluation[mstInCore[rep], params, goals, deriv, r]]];
  (* first evaluation of this mode: a single pass, to measure the loss *)
  res = mstPaddedEvaluation[mstInCore[rep], params, goals, deriv, r, 1];
  If[Length[$inRepresentationCache] >= 50, $inRepresentationCache = <||>];
  $inRepresentationCache[key] = $lastPaddingLoss;
  If[$lastPaddingLoss > seriesLoss, Return[mstPaddedEvaluation[mstRadialInSeriesCore, params, goals, deriv, r]]];
  If[!allNumericQ[res] || Precision[res] < If[goals[[1]] === MachinePrecision, $MachinePrecision, goals[[1]]] - 1,
    res = mstPaddedEvaluation[mstInCore[rep], params, goals, deriv, r, 4, $lastPaddingPrecision + Ceiling[$lastPaddingLoss] + 3];
  ];
  res
 ];

(* cores taking a derivative order, wrapping the series definitions above *)
(* Public evaluation: value, first derivative, or both from one summation (f[r, {0, 1}]) *)
MSTRadialIn[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_, {wp_, prec_, acc_}][r_?NumericQ] :=
 mstRadialInEvaluate[{s, l, m, q, \[Epsilon], \[Nu], \[Lambda], norm}, {wp, prec, acc}, 0, r];

Derivative[1][MSTRadialIn[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_, {wp_, prec_, acc_}]][r_?NumericQ] :=
 mstRadialInEvaluate[{s, l, m, q, \[Epsilon], \[Nu], \[Lambda], norm}, {wp, prec, acc}, 1, r];

MSTRadialIn[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_, {wp_, prec_, acc_}][r_?NumericQ, {0, 1}] :=
 mstRadialInEvaluate[{s, l, m, q, \[Epsilon], \[Nu], \[Lambda], norm}, {wp, prec, acc}, All, r];

MSTRadialUp[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_, {wp_, prec_, acc_}][r_?NumericQ] :=
 mstPaddedEvaluation[mstRadialUpSeriesCore, {s, l, m, q, \[Epsilon], \[Nu], \[Lambda], norm}, {wp, prec, acc}, 0, r];

Derivative[1][MSTRadialUp[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_, {wp_, prec_, acc_}]][r_?NumericQ] :=
 mstPaddedEvaluation[mstRadialUpSeriesCore, {s, l, m, q, \[Epsilon], \[Nu], \[Lambda], norm}, {wp, prec, acc}, 1, r];

MSTRadialUp[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_, {wp_, prec_, acc_}][r_?NumericQ, {0, 1}] :=
 mstPaddedEvaluation[mstRadialUpSeriesCore, {s, l, m, q, \[Epsilon], \[Nu], \[Lambda], norm}, {wp, prec, acc}, All, r];


(* ::Section::Closed:: *)
(*Second and higher derivatives*)


Derivative[n_Integer?Positive][(MSTR:MSTRadialIn|MSTRadialUp)[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_, {wp_, prec_, acc_}]][r0_?NumericQ] :=
 Module[{Rderivs, R, r, i, res},
  pderivs = D[R[r_], {r_, i_}] :> D[d2R[s, l, m, q, \[Epsilon], \[Lambda], r, R], {r, i - 2}] /; i >= 2;
  Do[Derivative[i][R][r] = Collect[D[Derivative[i - 1][R][r], r] /. pderivs,{R'[r], R[r]}, Simplify];, {i, 2, n}];
  res = Derivative[n][R][r] /. {
    R'[r] -> MSTR[s, l, m, q, \[Epsilon], \[Nu], \[Lambda], norm, {wp, prec, acc}]'[r0],
    R[r] -> MSTR[s, l, m, q, \[Epsilon], \[Nu], \[Lambda], norm, {wp, prec, acc}][r0], r -> r0, \[Epsilon]L -> \[Epsilon], qL -> q, \[Lambda]L -> \[Lambda]};
  Remove[R, r];
  res
];


(* ::Section::Closed:: *)
(*Listable functions*)


MSTRadialIn[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_, {wp_, prec_, acc_}][r:{_?NumericQ..}] :=
  Map[MSTRadialIn[s, l, m, q, \[Epsilon], \[Nu], \[Lambda], norm, {wp, prec, acc}], r];


MSTRadialUp[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_, {wp_, prec_, acc_}][r:{_?NumericQ..}] :=
  Map[MSTRadialUp[s, l, m, q, \[Epsilon], \[Nu], \[Lambda], norm, {wp, prec, acc}], r];


Derivative[n_][MSTRadialIn[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_, {wp_, prec_, acc_}]][r:{_?NumericQ..}] :=
  Map[Derivative[n][MSTRadialIn[s, l, m, q, \[Epsilon], \[Nu], \[Lambda], norm, {wp, prec, acc}]], r];


Derivative[n_][MSTRadialUp[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_, {wp_, prec_, acc_}]][r:{_?NumericQ..}] :=
  Map[Derivative[n][MSTRadialUp[s, l, m, q, \[Epsilon], \[Nu], \[Lambda], norm, {wp, prec, acc}]], r];


(* ::Section::Closed:: *)
(*End package*)


End[];
EndPackage[];
