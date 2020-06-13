(* ::Package:: *)

(* ::Title:: *)
(*MST Package*)


(* ::Section::Closed:: *)
(*Create Package*)


BeginPackage[MST`$MasterFunction<>"`MST`MST`"];

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
  prefacInInc[s_, \[Epsilon]_, \[Tau]_, \[Kappa]_, \[Nu]_, K\[Nu]1_, K\[Nu]2_] := (K\[Nu]1 - I E^(-I \[Pi] \[Nu]) Sin[\[Pi] (\[Nu] + I \[Epsilon])] / Sin[\[Pi] (\[Nu] - I \[Epsilon])] K\[Nu]2);,
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
  
  ret
];


(* ::Section::Closed:: *)
(*Asymptotic amplitudes*)


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
  termK\[Nu]1Up[n_] := termK\[Nu]1Up[n] = ((-1)^n Gamma[1 + n + s + I \[Epsilon] + \[Nu]] Gamma[1 + n + 2 \[Nu]] Gamma[1 + n + \[Nu] + I \[Tau]])/(n! Gamma[1 + n - s - I \[Epsilon] + \[Nu]] Gamma[1 + n + \[Nu] - I \[Tau]]) fn[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, n];
  termK\[Nu]1Down[n_] := termK\[Nu]1Down[n] = (((-1)^n) Pochhammer[1 + s - I \[Epsilon] + \[Nu], n])/((-n)! Pochhammer[1 - s + I \[Epsilon] + \[Nu], n] Pochhammer[2 + 2 \[Nu], n]) fn[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, n];
  fSumK\[Nu]1Up = fSumK\[Nu]1Down = 0;

  n = 0;
  While[fSumK\[Nu]1Up != (fSumK\[Nu]1Up += termK\[Nu]1Up[n]), n++];

  n = 0;
  While[fSumK\[Nu]1Down != (fSumK\[Nu]1Down += termK\[Nu]1Down[n]), n--];

  termK\[Nu]2Up[n_] := termK\[Nu]2Up[n] = ((-1)^n Gamma[1 + n + s + I \[Epsilon] + (-1-\[Nu])] Gamma[1 + n + 2 (-1-\[Nu])] Gamma[1 + n + (-1-\[Nu]) + I \[Tau]])/(n! Gamma[1 + n - s - I \[Epsilon] + (-1-\[Nu])] Gamma[1 + n + (-1-\[Nu]) - I \[Tau]]) fn[q, \[Epsilon], \[Kappa], \[Tau], (-1-\[Nu]), \[Lambda], s, m, n];
  termK\[Nu]2Down[n_] := termK\[Nu]2Down[n] = (((-1)^n) Pochhammer[1 + s - I \[Epsilon] + (-1-\[Nu]), n])/((-n)! Pochhammer[1 - s + I \[Epsilon] + (-1-\[Nu]), n] Pochhammer[2 + 2 (-1-\[Nu]), n]) fn[q, \[Epsilon], \[Kappa], \[Tau], (-1-\[Nu]), \[Lambda], s, m, n];
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
  Aminus = 2^(-s - 1 + I \[Epsilon]) E^(-\[Pi] \[Epsilon] / 2 - I \[Pi] (\[Nu]+1+s) / 2) (fSumAminusUp+fSumAminusDown);

  (* Up Transmission coefficient: Ctrans in ST (170), CO (3.20) *)
  UpTrans = prefacUpTrans[s, \[Epsilon], \[Tau], \[Kappa]] Aminus;

  (* K\[Nu]: ST (165), CO (3.32) *)
  K\[Nu]1 = ((2^-\[Nu]) (E^(I \[Epsilon] \[Kappa])) ((\[Epsilon] \[Kappa])^(s - \[Nu])) Gamma[1 - s - 2 I \[Epsilon]p] Gamma[2 + 2 \[Nu]])/(Gamma[1 - s + I \[Epsilon] + \[Nu]] Gamma[1 + s + I \[Epsilon] + \[Nu]] Gamma[1 + \[Nu] + I \[Tau]]) fSumK\[Nu]1Up / fSumK\[Nu]1Down;
  K\[Nu]2 = ((2^-(-1-\[Nu])) (E^(I \[Epsilon] \[Kappa])) ((\[Epsilon] \[Kappa])^(s - (-1-\[Nu]))) Gamma[1 - s - 2 I \[Epsilon]p] Gamma[2 + 2 (-1-\[Nu])])/(Gamma[1 - s + I \[Epsilon] + (-1-\[Nu])] Gamma[1 + s + I \[Epsilon] + (-1-\[Nu])] Gamma[1 + (-1-\[Nu]) + I \[Tau]]) fSumK\[Nu]2Up / fSumK\[Nu]2Down;

  (* In reflection coefficient: Bref in ST (169), CO (3.37) *)
  InRef = UpTrans (K\[Nu]1 + I E^(I \[Pi] \[Nu]) K\[Nu]2);

  (* D2 *)
  D2 = -Exp[(I \[Kappa] (\[Epsilon]+\[Tau]) (1+\[Kappa]+2 Log[\[Kappa]]))/(2 (1+\[Kappa]))] (2\[Kappa])^(2 s) ( Sin[\[Pi] (\[Nu]-I \[Epsilon])] Sin[\[Pi] (\[Nu]-I \[Tau])])/(Sin[2 \[Pi] \[Nu]] Sin[\[Pi] I (\[Epsilon]+\[Tau])]) (fSumUp+fSumDown);
  D22 = -Exp[(I \[Kappa] (\[Epsilon]+\[Tau]) (1+\[Kappa]+2 Log[\[Kappa]]))/(2 (1+\[Kappa]))] (2\[Kappa])^(2 s) ( Sin[\[Pi] ((-1-\[Nu])-I \[Epsilon])] Sin[\[Pi] ((-1-\[Nu])-I \[Tau])])/(Sin[2 \[Pi] (-1-\[Nu])] Sin[\[Pi] I (\[Epsilon]+\[Tau])]) (fSumUp+fSumDown);

  (* Up reflection coefficient *)
  UpRef = Exp[-\[Pi] \[Epsilon]-I \[Pi] s]/Sin[2\[Pi] \[Nu]] ((Exp[-I \[Pi] \[Nu]]Sin[\[Pi](\[Nu]-s+I \[Epsilon])])/K\[Nu]1 D2-I Sin[\[Pi](\[Nu]+s-I \[Epsilon])]/K\[Nu]2 D22);

  (* A+: ST (157), CO (3.38) and (3.41) *)
  Aplus = prefacAplus[s, \[Epsilon], \[Tau], \[Kappa], \[Nu]] (fSumUp+fSumDown);

  (* In incidence coefficient: Binc from ST (168), CO (3.36) and (3.39) *)
  InInc = prefacInInc[s, \[Epsilon], \[Tau], \[Kappa], \[Nu], K\[Nu]1, K\[Nu]2] Aplus;

  (* D1 *)
  D1 = Exp[-((I \[Kappa] (\[Epsilon]+\[Tau]) (1+\[Kappa]+2 Log[\[Kappa]]))/(2 (1+\[Kappa])))] ( Sin[\[Pi] (\[Nu]+I \[Epsilon])] Sin[\[Pi] (\[Nu]+I \[Tau])] Gamma[1-s-I (\[Epsilon]+\[Tau])])/(Sin[2 \[Pi] \[Nu]] Sin[\[Pi] I (\[Epsilon]+\[Tau])]  Gamma[1+s+I \[Epsilon]+I \[Tau]]) (fSumD1Up+fSumD1Down);
  D12 = Exp[-((I \[Kappa] (\[Epsilon]+\[Tau]) (1+\[Kappa]+2 Log[\[Kappa]]))/(2 (1+\[Kappa])))] ( Sin[\[Pi] ((-1-\[Nu])+I \[Epsilon])] Sin[\[Pi] ((-1-\[Nu])+I \[Tau])] Gamma[1-s-I (\[Epsilon]+\[Tau])])/(Sin[2 \[Pi] (-1-\[Nu])] Sin[\[Pi] I (\[Epsilon]+\[Tau])]  Gamma[1+s+I \[Epsilon]+I \[Tau]]) (fSumD12Up+fSumD12Down);

  (* Up incidence coefficient *)
  UpInc = Exp[-\[Pi] \[Epsilon]-I \[Pi] s]/Sin[2\[Pi] \[Nu]] ((Exp[-I \[Pi] \[Nu]]Sin[\[Pi](\[Nu]-s+I \[Epsilon])])/K\[Nu]1 D1-I Sin[\[Pi](\[Nu]+s-I \[Epsilon])]/K\[Nu]2 D12);

  (* Return results as an Association *)
  <| "In" -> <| "Incidence" -> InInc, "Transmission" -> InTrans, "Reflection" -> InRef|>,
     "Up" -> <| "Incidence" -> UpInc, "Transmission" -> UpTrans, "Reflection" -> UpRef |>
   |>
]];


(* ::Section::Closed:: *)
(*Radial "In" solution*)


(* ::Text:: *)
(*Ingoing MST Teukolsky Radial Function: Throwe B.1, Sasaki & Tagoshi Eqs. (116) and (120)*)
(*Ingoing MST Regge Wheeler Radial Function: Casals & Ottewill Eqs. (3.1) and (3.4a)*)
(*Ingoing Regge Wheeler MST coefficients: Casals & Ottewill Eqs. (3.8) & (3.9)*)


SetAttributes[MSTRadialIn, {NumericFunction}];


(* ::Subsection::Closed:: *)
(*Radial function*)


MSTRadialIn[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_, {wp_, prec_, acc_}][r_?NumericQ] :=
 Module[{\[Kappa], \[Tau], rp, x, resUp, nUp, resDown, nDown, term, prefac},
 Block[{H2F1},
 Internal`InheritedBlock[{\[Alpha], \[Beta], \[Gamma], fn},
  \[Kappa] = Sqrt[1 - q^2];
  \[Tau] = (\[Epsilon] - m q)/\[Kappa];
  rp = 1 + \[Kappa];
  x = (rp - r)/(2 \[Kappa]);

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
 
  prefac = prefacIn[s, \[Epsilon], \[Tau], \[Kappa], x]/norm;
  term[n_] := term[n] = prefac fIn[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, n]H2F1[n];
  resUp = resDown = 0;

  nUp = 0;
  While[resUp != (resUp += term[nUp]) && (Abs[term[nUp]] > 10^-acc + Abs[resUp] 10^-prec), nUp++];

  nDown = -1;
  While[resDown != (resDown += term[nDown]) && (Abs[term[nDown]] > 10^-acc + Abs[resDown] 10^-prec), nDown--];

  resUp + resDown
]]];


(* ::Subsection::Closed:: *)
(*First derivative*)


Derivative[1][MSTRadialIn[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_, {wp_, prec_, acc_}]][r_?NumericQ] :=
 Module[{\[Kappa], \[Tau], rp, x, dxdr, prefac, dprefac, resUp, nUp, resDown, nDown, term},
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
 
  prefac = prefacIn[s, \[Epsilon], \[Tau], \[Kappa], x] dxdr/norm;
  dprefac = Derivative[0,0,0,0,1][prefacIn][s, \[Epsilon], \[Tau], \[Kappa], x] dxdr/norm;

  term[n_] := term[n] = fIn[q,\[Epsilon],\[Kappa],\[Tau],\[Nu],\[Lambda],s,m,n](dprefac H2F1[n] + prefac dH2F1[n]);
  resUp = resDown = 0;

  nUp = 0;
  While[resUp != (resUp+= term[nUp]) && (Abs[term[nUp]] > 10^-acc + Abs[resUp] 10^-prec), nUp++];

  nDown = -1;
  While[resDown != (resDown+= term[nDown]) && (Abs[term[nDown]] > 10^-acc + Abs[resDown] 10^-prec), nDown--];

  (resUp+resDown)
]]];


(* ::Section::Closed:: *)
(*Radial "Up" solution*)


(* ::Text:: *)
(*Upgoing MST Teukolsky Radial Function: Throwe B.5, Sasaki & Tagoshi Eqs. (153) and (159)*)
(*Upgoing MST Regge Wheeler Radial Function: Casals & Ottewill Eq. (3.15)*)
(*Regge Wheeler MST coefficients: Casals & Ottewill Eqs. (3.8) & (3.9)*)


SetAttributes[MSTRadialUp, {NumericFunction}];


(* ::Subsection::Closed:: *)
(*Radial function*)


MSTRadialUp[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_, {wp_, prec_, acc_}][r_?NumericQ] :=
 Module[{\[Kappa], \[Tau], \[Epsilon]p, rm, z, zm, zhat, resUp, nUp, resDown, nDown, term, prefac},
 Block[{HU},
 Internal`InheritedBlock[{\[Alpha], \[Beta], \[Gamma], fn},
  \[Kappa] = Sqrt[1 - q^2];
  \[Tau] = (\[Epsilon] - m q)/\[Kappa];
  \[Epsilon]p = 1/2 (\[Tau]+\[Epsilon]);
  rm = 1 - Sqrt[1 - q^2];
  z = \[Epsilon] r / 2;
  zm = \[Epsilon] rm / 2;
  zhat = z - zm;
 
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

  prefac = prefacUp[s, \[Epsilon], \[Kappa], \[Tau], \[Nu], zhat]/norm;
  term[n_] := term[n] = prefac fUp[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, n] HU[n];
  resDown = resUp = 0;
  nUp = 0;
  While[resUp != (resUp += term[nUp]) && (Abs[term[nUp]] > 10^-acc + Abs[resUp] 10^-prec), nUp++];
  
  nDown = -1;
  While[resDown != (resDown += term[nDown]) && (Abs[term[nDown]] > 10^-acc + Abs[resDown] 10^-prec), nDown--];
  
  resUp+resDown
]]];


(* ::Subsection::Closed:: *)
(*First derivative*)


Derivative[1][MSTRadialUp[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_, {wp_, prec_, acc_}]][r_?NumericQ] :=
 Module[{\[Kappa], \[Tau], \[Epsilon]p, rm, z, zm, zhat, dzhatdr, prefac, dprefac, resUp, nUp, resDown, nDown, term},
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

  prefac = prefacUp[s, \[Epsilon], \[Kappa], \[Tau], \[Nu], zhat] dzhatdr/norm;
  dprefac = Derivative[0,0,0,0,0,1][prefacUp][s, \[Epsilon], \[Kappa], \[Tau], \[Nu], zhat] dzhatdr/norm;

  term[n_] := term[n] = fUp[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, n] (dprefac HU[n] + prefac dHU[n]);
  resDown = resUp = 0;
  nUp = 0;
  While[resUp != (resUp += term[nUp]) && (Abs[term[nUp]] > 10^-acc + Abs[resUp] 10^-prec), nUp++];
  
  nDown = -1;
  While[resDown != (resDown += term[nDown]) && (Abs[term[nDown]] > 10^-acc + Abs[resDown] 10^-prec), nDown--];
  
  (resUp+resDown)
]]];


(* ::Section::Closed:: *)
(*Second and higher derivatives*)


Derivative[n_Integer?Positive][(MSTR:MSTRadialIn|MSTRadialUp)[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_, {wp_, prec_, acc_}]][r0_?NumericQ] :=
 Module[{Rderivs, R, r, i},
  pderivs = D[R[r_], {r_, i_}] :> D[d2R[s, l, m, q, \[Epsilon], \[Lambda], r, R], {r, i - 2}] /; i >= 2;
  Do[Derivative[i][R][r] = Collect[D[Derivative[i - 1][R][r], r] /. pderivs,{R'[r], R[r]}, Simplify];, {i, 2, n}];
  Derivative[n][R][r] /. {
    R'[r] -> MSTR[s, l, m, q, \[Epsilon], \[Nu], \[Lambda], norm, {wp, prec, acc}]'[r0],
    R[r] -> MSTR[s, l, m, q, \[Epsilon], \[Nu], \[Lambda], norm, {wp, prec, acc}][r0], r -> r0, \[Epsilon]L -> \[Epsilon], qL -> q, \[Lambda]L -> \[Lambda]}
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
