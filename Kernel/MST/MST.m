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
  (* bounded: convergent fractions need tens to a few hundred terms (1000 is ample); on terms without correct digits (a
     Cos[2 Pi nu] of no precision in the estimate of the monodromy method, at 24 digits for s = 1, l = 8,
     omega = -2.218 - 0.180 I) the convergents never repeat and the loop ran until the kernel died *)
  While[res =!= (res = A[j]/B[j]), If[j - n0 > 1000, res = Indeterminate; Break[]]; j++];
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


(* ::Text:: *)
(*The hypergeometric functions of the "In" series are the regularised ones, 2F1(a, b; c; x)/Gamma(c), which are entire in c.*)
(*The "In" solution and its unscaled asymptotic amplitudes are therefore those of Sasaki & Tagoshi divided by Gamma(c),*)
(*c = 1 - s - 2 I epsilon_+, and stay finite at the degeneracies 2 I epsilon_+ = n >= 1 - s, where c is a non-positive*)
(*integer and the unregularised series has a pole. There the regularised series is the solution of larger exponent*)
(*at the horizon, whose transmission amplitude vanishes. The recurrences below are linear and homogeneous in the*)
(*hypergeometric functions, so they hold unchanged for the regularised ones.*)


(* ::Subsection::Closed:: *)
(*Hypergeometric2F1*)


(* At a degeneracy 2 I epsilon_+ = n >= 1 - s the parameter c is a non-positive integer up to the rounding of
   the inexact inputs, e.g. -1 + 6 10^-58 I with two digits of precision at 60-digit working precision.
   The regularised functions are entire in c, so such a c is taken to be the integer itself: Mathematica
   evaluates them stably for an exact integer c (and for a genuinely small offset carrying full precision)
   but returns Indeterminate for an offset that is zero to working precision, in particular for c + 1 in the
   derivative. *)
(* the exponent is an integer so that the power is exact: a machine-number power underflows beyond 308 digits *)
cSnapped[c_] := With[{n = Round[Re[c]], p = Precision[c]}, If[n <= 0 && Abs[c - n] <= 10^(3 - Floor[If[p === MachinePrecision, $MachinePrecision, p]]), n, c]];

H2F1Exact[n_, s_, \[Nu]_, \[Tau]_, \[Epsilon]_, x_] :=
 Module[{a = aF[s, \[Nu], \[Tau], \[Epsilon]], b = bF[s, \[Nu], \[Tau], \[Epsilon]], c = cSnapped[cF[s, \[Nu], \[Tau], \[Epsilon]]]},
  Hypergeometric2F1Regularized[n + a, b-n, c, x]
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
 Module[{a = aF[s, \[Nu], \[Tau], \[Epsilon]], b = bF[s, \[Nu], \[Tau], \[Epsilon]], c = cSnapped[cF[s, \[Nu], \[Tau], \[Epsilon]]]},
  (n+a)(-n+b) Hypergeometric2F1Regularized[n + a + 1, -n + b + 1, c + 1, x]
];

(* The second Kummer solution of the same hypergeometric equation, x^(1-c) F(a-c+1, b-c+1; 2-c; x), which
   carries the outgoing horizon exponent, normalised by Gamma[a-c+1] Gamma[b-c+1]/(Gamma[a] Gamma[b]) so
   that it satisfies the same contiguous relations in a and b as F(a, b; c; x) itself, and therefore the
   recurrences H2F1Up and H2F1Down: with the same coefficients as the "In" series it then gives the second
   solution of the radial equation at the horizon (see mstRadialUpHorizon). *)
outNormalisation[n_, a_, b_, c_] := Gamma[n + a - c + 1] Gamma[b - n - c + 1]/(Gamma[n + a] Gamma[b - n]);

H2F1OutExact[n_, s_, \[Nu]_, \[Tau]_, \[Epsilon]_, x_] :=
 Module[{a = aF[s, \[Nu], \[Tau], \[Epsilon]], b = bF[s, \[Nu], \[Tau], \[Epsilon]], c = cF[s, \[Nu], \[Tau], \[Epsilon]]},
  outNormalisation[n, a, b, c] (-x)^(1 - c) Hypergeometric2F1[n + a - c + 1, b - n - c + 1, 2 - c, x]
];

dH2F1OutExact[n_, s_, \[Nu]_, \[Tau]_, \[Epsilon]_, x_] :=
 Module[{a = aF[s, \[Nu], \[Tau], \[Epsilon]], b = bF[s, \[Nu], \[Tau], \[Epsilon]], c = cF[s, \[Nu], \[Tau], \[Epsilon]], ap, bp, cp},
  {ap, bp, cp} = {n + a - c + 1, b - n - c + 1, 2 - c};
  outNormalisation[n, a, b, c] (-(1 - c) (-x)^(-c) Hypergeometric2F1[ap, bp, cp, x] + (-x)^(1 - c) ap bp/cp Hypergeometric2F1[ap + 1, bp + 1, cp + 1, x])
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


(* The Coulomb-type series are derived for Re epsilon > 0, where the argument c = -2 I zhat of the
   confluent hypergeometric functions has -Pi < Arg[c] < 0. At a purely imaginary frequency with
   Im epsilon < 0 the argument is negative real, on the branch cut of HypergeometricU, and the value
   continuous with Re epsilon > 0 is the limit from below the cut, whereas Mathematica's principal value
   is the limit from above. We evaluate U just below the cut, with an imaginary part far below the
   precision goal (raising the precision of the argument so that the shift is representable): the
   connection formula for the two sides of the cut (DLMF 13.2.12) suffers from cancellation between its
   terms and loses up to fifteen digits at machine precision for |Im a| of a few. Frequencies with
   Re epsilon < 0 are handled by the conjugation symmetry (see MSTRadialIn), so the negative real axis
   is the only part of the cut that is reached. *)
(* The side of the cut depends on the series: for the outgoing series R_- ("Up") the argument -2 I zhat is
   negative real at epsilon = -I t, and the value continuous with Re epsilon > 0 is the limit from below; the
   incoming series R_+ (the Coulomb-type "In" representation, mstRadialPlusSeries) reaches the same function
   with the same epsilon and zhat but at epsilon = +I t, where the continuous value is the limit from above,
   Mathematica's principal value. mstRadialPlusSeries therefore evaluates with $uLowerSide = False. *)
$uLowerSide = True;
hypergeometricU[a_, b_, c_] /; $uLowerSide && Im[c] == 0 && Re[c] < 0 :=
  With[{p = Precision[{a, b, c}]},
    Which[
      p === MachinePrecision, HypergeometricU[a, b, c - I 10^-26 Abs[c]],
      p === Infinity, HypergeometricU[a, b, c - I 10^-60 Abs[c]],
      (* the arguments are evaluated at their lowest common precision, so all three are raised, and the
         result is set back to that precision *)
      True, SetPrecision[HypergeometricU[SetPrecision[a, p + 20], SetPrecision[b, p + 20], SetPrecision[c, p + 20] - I 10^(-Floor[p] - 10) Abs[c]], p]]];
hypergeometricU[a_, b_, c_] := HypergeometricU[a, b, c];

HUExact[n_, s_, \[Nu]_, \[Epsilon]_, zhat_] :=
 Module[{a = aU[s, \[Nu], \[Tau], \[Epsilon]], b = 2 \[Nu] + 2, c = -2 I zhat},
  (c)^n hypergeometricU[n+a,2n+b,c]
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
  (-2 I) (c^(-1+n) n hypergeometricU[a+n,b+2 n,c]-c^n (a+n) hypergeometricU[1+a+n,1+b+2 n,c])
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
(* bounded: a non-convergent sum (garbage parameters at an absurd frequency) returns Indeterminate instead of
   running forever; convergent sums need at most a few hundred terms *)
sumUntil[term_, n0_, dir_] := Module[{res = 0, k = n0, t}, While[res != (res += (t = term[k])), If[!NumericQ[t] || Abs[k - n0] > 5000, Return[Indeterminate, Module]]; k += dir]; res];

(* K_nu / Gamma(1 - s - 2 I epsilon_+), ST Eq. (165) with r = 0, CO (3.32): the connection coefficient of the
   regularised "In" series (see the note on the hypergeometric functions), finite at 2 I epsilon_+ = n >= 1 - s *)
KCoefficient[s_Integer, m_Integer, q_, \[Epsilon]_, \[Kappa]_, \[Tau]_, \[Nu]_, \[Lambda]_] :=
 Module[{},
  ((2^-\[Nu]) (E^(I \[Epsilon] \[Kappa])) ((\[Epsilon] \[Kappa])^(s - \[Nu])) Gamma[2 + 2 \[Nu]])/(Gamma[1 - s + I \[Epsilon] + \[Nu]] Gamma[1 + s + I \[Epsilon] + \[Nu]] Gamma[1 + \[Nu] + I \[Tau]]) *
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


(* Re epsilon < 0: the conjugate partner's amplitudes, conjugated (see MSTRadialIn; the powers of epsilon in
   the amplitude formulae are on the wrong side of their cuts there) *)
Amplitudes[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, {wp_, prec_, acc_}] /; Re[\[Epsilon]] < 0 :=
  Map[Conjugate, Amplitudes[s, l, -m, q, -Conjugate[\[Epsilon]], Conjugate[\[Nu]], Conjugate[\[Lambda]], {wp, prec, acc}], {2}];

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

  (* In transmission coefficient: Btrans in ST (167) and CO (3.12), divided by Gamma[c] = Gamma[1 - 2 I \[Epsilon]]
     like the "In" series (see the note on the hypergeometric functions) *)
  InTrans = prefacInTrans[s, \[Epsilon], \[Tau], \[Kappa]] (fSumUp+fSumDown) / Gamma[1 - 2 I \[Epsilon]];

  (* A-: ST (158), CO (3.19) *)
  Aminus = 2^(-s - 1 + I \[Epsilon]) E^(-\[Pi] \[Epsilon] / 2 - I \[Pi] (\[Nu]+1+s) / 2) (fSumAminusUp+fSumAminusDown);

  (* Up Transmission coefficient: Ctrans in ST (170), CO (3.20) *)
  UpTrans = prefacUpTrans[s, \[Epsilon], \[Tau], \[Kappa]] Aminus;

  (* K\[Nu]: ST (165), CO (3.32) *)
  K\[Nu]1 = ((2^-\[Nu])( E^(I \[Epsilon])) \[Epsilon]^(-1-\[Nu]) Gamma[1-s-2 I \[Epsilon]] Gamma[1+s-I \[Epsilon]+\[Nu]])/(Gamma[-I \[Epsilon]-\[Nu]] Gamma[1+I \[Epsilon]+\[Nu]] Gamma[1-s+I \[Epsilon]+\[Nu]]) fSumK\[Nu]1Up / fSumK\[Nu]1Down;
  K\[Nu]2 = ((2^-(-1-\[Nu]))( E^(I \[Epsilon])) \[Epsilon]^(-1-(-1-\[Nu])) Gamma[1-s-2 I \[Epsilon]] Gamma[1+s-I \[Epsilon]+(-1-\[Nu])])/(Gamma[-I \[Epsilon]-(-1-\[Nu])] Gamma[1+I \[Epsilon]+(-1-\[Nu])] Gamma[1-s+I \[Epsilon]+(-1-\[Nu])]) fSumK\[Nu]2Up / fSumK\[Nu]2Down;

  (* In reflection coefficient: Bref in ST (169), CO (3.37), divided by Gamma[c] = Gamma[1-2 I \[Epsilon]] like the "In" series *)
  InRef = (Gamma[-I \[Epsilon]-\[Nu]] Gamma[1-I \[Epsilon]+\[Nu]])/(Gamma[1-s-2 I \[Epsilon]] Gamma[s-I \[Epsilon]-\[Nu]] Gamma[1+s-I \[Epsilon]+\[Nu]]) UpTrans (K\[Nu]1 + I E^(I \[Pi] \[Nu]) K\[Nu]2);

  (* A+: ST (157), CO (3.38) and (3.41) *)
  Aplus = prefacAplus[s, \[Epsilon], \[Tau], \[Kappa], \[Nu]] (fSumUp+fSumDown);

  (* In incidence coefficient: Binc from ST (168), CO (3.36) and (3.39), divided by Gamma[c] = Gamma[1-2 I \[Epsilon]] *)
  InInc = (Gamma[-I \[Epsilon]-\[Nu]] Gamma[1-I \[Epsilon]+\[Nu]])/(Gamma[1-s-2 I \[Epsilon]] Gamma[s-I \[Epsilon]-\[Nu]] Gamma[1+s-I \[Epsilon]+\[Nu]]) prefacInInc[s, \[Epsilon], \[Tau], \[Kappa], \[Nu], K\[Nu]1, K\[Nu]2] Aplus;

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
 Module[{\[Kappa], \[Tau], \[Omega], x, xs, nDeg, sRef, sInc, K\[Nu], K\[Nu]1, K\[Nu]2, Aminus, Aplus, D1, D12, D2, D22, InTrans, UpTrans, InInc, UpInc, InRef, UpRef, n, fSumUp, fSumDown, fSumD1Up, fSumD1Down, fSumD12Up, fSumD12Down, termf, termD1, termD12},
 Internal`InheritedBlock[{\[Alpha], \[Beta], \[Gamma], fn},
  \[Kappa] = Sqrt[1 - q^2];
  \[Tau] = (\[Epsilon] - m q)/\[Kappa];
  \[Omega] = \[Epsilon] / 2;

  (* Degeneracies x = 2 I epsilon_+ = n, an integer: for real frequencies the superradiant bound frequency
     omega = m Omega_H (n = 0, k = 0), on the imaginary axis every n. There the horizon basis Exp[i k r_*],
     Delta^-s Exp[-i k r_*] is resonant (its Frobenius exponents at r_+ differ by the integer s + n; the two
     solutions coincide for s = n = 0) and the MST formulae are singular in two explicit scalar factors: the
     "In" amplitudes carry Gamma[1 - s - x], divided out below and in the "In" series (see the note on the
     hypergeometric functions), so that they are finite, with a transmission amplitude that vanishes for
     n >= 1 - s (the regularised "In" solution is then the horizon solution of larger exponent and is
     normalised to unit incidence by the caller); the "Up" horizon coefficients carry
     1/(Gamma[1 - s - x] Sin[Pi x]) (reflection) and 1/(Gamma[1 + s + x] Sin[Pi x]) (incidence). By
     Gamma[x] Gamma[1 - x] = Pi/Sin[Pi x] these have the finite limits (-1)^s (s + n - 1)!/Pi and
     (-1)^(s + 1) (-s - n - 1)!/Pi when the factorial's argument is non-negative, and are poles otherwise:
     the coefficients of a finite solution along the basis function that has ceased to be independent,
     returned as Indeterminate. x is snapped to n within the tolerance of cSnapped, so that a frequency
     numerically at a degeneracy is treated exactly. *)
  x = I (\[Epsilon] + \[Tau]);
  nDeg = With[{n = Round[Re[x]], p = Precision[x]}, If[Abs[x - n] <= 10^(3 - Floor[If[p === MachinePrecision, $MachinePrecision, p]]), n, None]];
  xs = If[nDeg === None, x, nDeg];
  {sRef, sInc} = If[nDeg === None,
    {1/(Gamma[1 - s - x] Sin[\[Pi] x]), 1/(Gamma[1 + s + x] Sin[\[Pi] x])},
    {If[s + nDeg - 1 >= 0, (-1)^s (s + nDeg - 1)!/\[Pi], Indeterminate], If[-s - nDeg - 1 >= 0, (-1)^(s + 1) (-s - nDeg - 1)!/\[Pi], Indeterminate]}];

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
  
  (* In transmission coefficient: Btrans in ST (167) and CO (3.12), divided by Gamma[c] like the "In" series
     (see the note on the hypergeometric functions); it vanishes at 2 I \[Epsilon]_+ = n >= 1 - s *)
  InTrans = prefacInTrans[s, \[Epsilon], \[Tau], \[Kappa]] (fSumUp+fSumDown) / Gamma[1 - s - xs];

  (* A-: ST (158), CO (3.19) *)
  Aminus = AminusCoefficient[s, m, q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda]];

  (* Up Transmission coefficient: Ctrans in ST (170), CO (3.20) *)
  UpTrans = prefacUpTrans[s, \[Epsilon], \[Tau], \[Kappa]] Aminus;

  (* K\[Nu]/Gamma[c]: ST (165), CO (3.32), for \[Nu] and for -\[Nu]-1, with the factor Gamma[c] = Gamma[1 - s - 2 I \[Epsilon]_+]
     removed. The "In" amplitudes built from them are then divided by Gamma[c] like the "In" series; the "Up"
     horizon coefficients, in which K\[Nu] appears in the denominator, are corrected below. *)
  K\[Nu]1 = KCoefficient[s, m, q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda]];
  K\[Nu]2 = KCoefficient[s, m, q, \[Epsilon], \[Kappa], \[Tau], -1 - \[Nu], \[Lambda]];

  (* In reflection coefficient: Bref in ST (169), CO (3.37), divided by Gamma[c] *)
  InRef = UpTrans (K\[Nu]1 + I E^(I \[Pi] \[Nu]) K\[Nu]2);

  (* D2 Sin[Pi x]: the common factor 1/Sin[Pi x] of D2 and D22 is in sRef *)
  D2 = -Exp[(I \[Kappa] (\[Epsilon]+\[Tau]) (1+\[Kappa]+2 Log[\[Kappa]]))/(2 (1+\[Kappa]))] (2\[Kappa])^(2 s) ( Sin[\[Pi] (\[Nu]-I \[Epsilon])] Sin[\[Pi] (\[Nu]-I \[Tau])])/Sin[2 \[Pi] \[Nu]] (fSumUp+fSumDown);
  D22 = -Exp[(I \[Kappa] (\[Epsilon]+\[Tau]) (1+\[Kappa]+2 Log[\[Kappa]]))/(2 (1+\[Kappa]))] (2\[Kappa])^(2 s) ( Sin[\[Pi] ((-1-\[Nu])-I \[Epsilon])] Sin[\[Pi] ((-1-\[Nu])-I \[Tau])])/Sin[2 \[Pi] (-1-\[Nu])] (fSumUp+fSumDown);

  (* Up reflection coefficient; 1/K\[Nu] = (1/Gamma[c]) / (K\[Nu]/Gamma[c]), and 1/(Gamma[c] Sin[Pi x]) = sRef *)
  UpRef = Exp[-\[Pi] \[Epsilon]-I \[Pi] s]/Sin[2\[Pi] \[Nu]] sRef ((Exp[-I \[Pi] \[Nu]]Sin[\[Pi](\[Nu]-s+I \[Epsilon])])/K\[Nu]1 D2-I Sin[\[Pi](\[Nu]+s-I \[Epsilon])]/K\[Nu]2 D22);

  (* A+: ST (157), CO (3.38) and (3.41) *)
  Aplus = AplusCoefficient[s, m, q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda]];

  (* In incidence coefficient: Binc from ST (168), CO (3.36) and (3.39), divided by Gamma[c] *)
  InInc = prefacInInc[s, \[Epsilon], \[Tau], \[Kappa], \[Nu], K\[Nu]1, K\[Nu]2] Aplus;

  (* D1 Gamma[1 + s + x] Sin[Pi x]/Gamma[c]: the factor Gamma[c] = Gamma[1 - s - x] of D1 cancels against the
     1/K\[Nu] in the "Up" incidence, and the common factor 1/(Gamma[1 + s + x] Sin[Pi x]) of D1 and D12 is in sInc *)
  D1 = Exp[-((I \[Kappa] (\[Epsilon]+\[Tau]) (1+\[Kappa]+2 Log[\[Kappa]]))/(2 (1+\[Kappa])))] ( Sin[\[Pi] (\[Nu]+I \[Epsilon])] Sin[\[Pi] (\[Nu]+I \[Tau])])/Sin[2 \[Pi] \[Nu]] (fSumD1Up+fSumD1Down);
  D12 = Exp[-((I \[Kappa] (\[Epsilon]+\[Tau]) (1+\[Kappa]+2 Log[\[Kappa]]))/(2 (1+\[Kappa])))] ( Sin[\[Pi] ((-1-\[Nu])+I \[Epsilon])] Sin[\[Pi] ((-1-\[Nu])+I \[Tau])])/Sin[2 \[Pi] (-1-\[Nu])] (fSumD12Up+fSumD12Down);

  (* Up incidence coefficient; (D1/Gamma[c]) / (K\[Nu]/Gamma[c]) = D1/K\[Nu] *)
  UpInc = Exp[-\[Pi] \[Epsilon]-I \[Pi] s]/Sin[2\[Pi] \[Nu]] sInc ((Exp[-I \[Pi] \[Nu]]Sin[\[Pi](\[Nu]-s+I \[Epsilon])])/K\[Nu]1 D1-I Sin[\[Pi](\[Nu]+s-I \[Epsilon])]/K\[Nu]2 D12);

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
(* A non-numeric term (Indeterminate, ComplexInfinity) ends the summation with Indeterminate, which the
   padded evaluation then retries or reports; the convergence test below would never be met by it and the
   loop would not terminate. *)
sumSeries[term_, n0_, dir_, prec_, acc_] :=
 Module[{res = 0, old, t, n = n0},
  While[True,
    t = term[n];
    If[!AllTrue[Flatten[{t}], NumericQ], Return[Indeterminate]];
    old = res; res = old + t;
    If[!(Or @@ Thread[Flatten[{res}] != Flatten[{old}]] && Or @@ Thread[Abs[Flatten[{t}]] > 10^-acc + Abs[Flatten[{res}]] 10^-prec]), Break[]];
    n += dir;
  ];
  res
 ];

(* The hypergeometric series for the "In" solution (Sasaki & Tagoshi Eq. (116)) at r: the value (deriv 0),
   the first derivative (deriv 1), or both from a single summation (deriv All), which shares the
   coefficients and the hypergeometric functions between the two. *)
(* The series of hypergeometric functions about the horizon, built on the regularised F(a, b; c; x) (the
   "In" solution, exact = H2F1Exact) or on the normalised second Kummer solution (the outgoing horizon
   solution, exact = H2F1OutExact); both share the coefficients and the recurrences. *)
mstRadialInSeriesCore[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_, {wp_, prec_, acc_}, deriv_][r_?NumericQ] :=
  mstRadialHorizonSeriesCore[H2F1Exact, dH2F1Exact][s, l, m, q, \[Epsilon], \[Nu], \[Lambda], norm, {wp, prec, acc}, deriv][r];

mstRadialOutSeriesCore[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_, {wp_, prec_, acc_}, deriv_][r_?NumericQ] :=
  mstRadialHorizonSeriesCore[H2F1OutExact, dH2F1OutExact][s, l, m, q, \[Epsilon], \[Nu], \[Lambda], norm, {wp, prec, acc}, deriv][r];

mstRadialHorizonSeriesCore[exact_, dexact_][s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_, {wp_, prec_, acc_}, deriv_][r_?NumericQ] :=
 Module[{\[Kappa], \[Tau], rp, x, dxdr, prefac, dprefac, term, res},
 Block[{H2F1, dH2F1},
 Internal`InheritedBlock[{\[Alpha], \[Beta], \[Gamma], fn},
  \[Kappa] = Sqrt[1 - q^2];
  \[Tau] = (\[Epsilon] - m q)/\[Kappa];
  rp = 1 + \[Kappa];
  x = (rp - r)/(2 \[Kappa]);
  dxdr = - 1/(2\[Kappa]);

  H2F1[n : (0 | 1)] := H2F1[n] = exact[n, s, \[Nu], \[Tau], \[Epsilon], x];

  H2F1[n_Integer] := H2F1[n] =
   Module[{t1, t2, res},
    {t1, t2} = If[n>0, H2F1Up[n, s, \[Nu], \[Tau], \[Epsilon], x], H2F1Down[n, s, \[Nu], \[Tau], \[Epsilon], x]];
    res = t1 + t2;
    If[Max[Abs[{t1, t2}/res]] > 2.,
      res = exact[n, s, \[Nu], \[Tau], \[Epsilon], x];
    ];
    res
  ];

  dH2F1[n : (0 | 1)] := dH2F1[n] = dexact[n, s, \[Nu], \[Tau], \[Epsilon], x];

  dH2F1[n_Integer] := dH2F1[n] =
   Module[{t1, t2, t3, res},
    {t1, t2, t3} = If[n>0, dH2F1Up[n, s, \[Nu], \[Tau], \[Epsilon], x], dH2F1Down[n, s, \[Nu], \[Tau], \[Epsilon], x]];
    res = t1 + t2 + t3;
    If[Max[Abs[{t1, t2, t3}/res]] > 2.,
      res = dexact[n, s, \[Nu], \[Tau], \[Epsilon], x];
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


(* ::Subsection::Closed:: *)
(*"Up" solution near the horizon: the horizon basis of hypergeometric series*)


(* ::Text:: *)
(*The Coulomb-type series is an expansion about infinity and converges slowly near the horizon (hundreds of*)
(*Tricomi functions at r+ + 0.05, thousands at r+ + 0.01, where an evaluation took minutes at 40 digits).*)
(*There the "Up" solution is written in the horizon basis instead, mirroring the Coulomb-type representation*)
(*of the "In" solution at large radius: R_up = c1 R_in + c2 R_out, with R_in the regularised "In" series and*)
(*R_out the series built on the normalised second Kummer solution (mstRadialOutSeriesCore), whose horizon*)
(*behaviour is (-x)^(I epsilon_+), the outgoing one. The coefficients follow from the horizon amplitudes:*)
(*c1 = C^ref/B^trans (the unscaled amplitudes of the package, both relative to the same regularised "In"*)
(*normalisation), and c2 = C^inc Exp[I (epsilon + tau) kappa (1/2 + Log[kappa]/(1 + kappa))]/Sum[a_n g_n],*)
(*where the phase is that of the "In" transmission prefactor (prefacInTrans without 4^s kappa^(2 s)) and*)
(*Sum[a_n g_n] is the coefficient of (-x)^(I epsilon_+) in R_out, g_n the normalisation of the second*)
(*solution (outNormalisation); both were verified against the Coulomb-type series to 33 digits. At the*)
(*degeneracies 2 I epsilon_+ = n the second Kummer solution has a logarithm and one of C^inc, C^ref is*)
(*Indeterminate, so the Coulomb-type series is kept there. Only for the Teukolsky master function.*)


(* r - r+ below which the "Up" solution uses the horizon basis: at r+ + 1 the two representations cost about
   the same (0.1 s at 40 digits), at r+ + 0.3 the Coulomb-type series costs ten times more and at r+ + 0.05 two
   hundred times more, while the horizon series, whose loss grows with the radius like the "In" series', is
   still at full precision *)
$MSTUpHorizonThreshold = 1;

mstUpRepresentation[s_, m_, q_, \[Epsilon]_, r_] :=
 Module[{\[Kappa] = Sqrt[1 - q^2], x, n, p},
  If[$masterFunction =!= "Teukolsky" || !(r - (1 + \[Kappa]) < $MSTUpHorizonThreshold), Return["Coulomb"]];
  (* a degeneracy 2 I epsilon_+ = n, within the tolerance the amplitudes snap with *)
  x = I (\[Epsilon] + (\[Epsilon] - m q)/\[Kappa]); n = Round[Re[x]]; p = Precision[x];
  If[Abs[x - n] <= 10^(3 - Floor[If[p === MachinePrecision, $MachinePrecision, p]]), "Coulomb", "Horizon"]
 ];

$upHorizonCache = <||>;   (* bounded cache of {c1, c2}, keyed by the parameters and the goals *)

(* the goals are part of the key: the coefficients come from the amplitude formulae, whose sums are truncated
   according to them, so coefficients computed with looser goals must not be reused for stricter ones *)
upHorizonCoefficients[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, goals_] :=
 Module[{key = {s, l, m, q, \[Epsilon], \[Nu], \[Lambda], goals}, res},
  res = Lookup[$upHorizonCache, Key[key], None];
  If[res =!= None, Return[res]];
  res = upHorizonCoefficientsCompute[s, l, m, q, \[Epsilon], \[Nu], \[Lambda], goals];
  If[Length[$upHorizonCache] >= 50, $upHorizonCache = <||>];
  $upHorizonCache[key] = res
 ];

upHorizonCoefficientsCompute[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, goals_] :=
 Module[{\[Kappa], \[Tau], a, b, c, amps, sumG},
 Internal`InheritedBlock[{\[Alpha], \[Beta], \[Gamma], fn},
  \[Kappa] = Sqrt[1 - q^2];
  \[Tau] = (\[Epsilon] - m q)/\[Kappa];
  amps = Amplitudes[s, l, m, q, \[Epsilon], \[Nu], \[Lambda], goals];
  If[!(AssociationQ[amps] && NumericQ[amps["Up"]["Reflection"]] && NumericQ[amps["Up"]["Incidence"]] && NumericQ[amps["In"]["Transmission"]] && amps["In"]["Transmission"] != 0),
    Return[$Failed]];
  {a, b, c} = {aF[s, \[Nu], \[Tau], \[Epsilon]], bF[s, \[Nu], \[Tau], \[Epsilon]], cF[s, \[Nu], \[Tau], \[Epsilon]]};
  sumG = sumUntil[fn[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, #] outNormalisation[#, a, b, c] &, 0, 1] + sumUntil[fn[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, #] outNormalisation[#, a, b, c] &, -1, -1];
  {amps["Up"]["Reflection"]/amps["In"]["Transmission"], amps["Up"]["Incidence"] Exp[I (\[Epsilon] + \[Tau]) \[Kappa] (1/2 + Log[\[Kappa]]/(1 + \[Kappa]))]/sumG}
 ]];

mstRadialUpHorizon[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_, {wp_, prec_, acc_}, deriv_][r_?NumericQ] :=
 Module[{cs = upHorizonCoefficients[s, l, m, q, \[Epsilon], \[Nu], \[Lambda], {wp, prec, acc}]},
  If[cs === $Failed, Return[Indeterminate]];
  (cs[[1]] mstRadialInSeriesCore[s, l, m, q, \[Epsilon], \[Nu], \[Lambda], 1, {wp, prec, acc}, deriv][r]
   + cs[[2]] mstRadialOutSeriesCore[s, l, m, q, \[Epsilon], \[Nu], \[Lambda], 1, {wp, prec, acc}, deriv][r])/norm
 ];


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
 Block[{HU, dHU, $uLowerSide = False},   (* principal side of the cut, see hypergeometricU *)
 Internal`InheritedBlock[{\[Alpha], \[Beta], \[Gamma], fn},
  \[Kappa] = Sqrt[1 - q^2];
  \[Tau] = (\[Epsilon] - m q)/\[Kappa];
  rm = 1 - \[Kappa];
  zhat = \[Epsilon] (r - rm)/2;
  \[Eta] = -I s - \[Epsilon];
  (* n-independent prefactor: prefacUp times the factors relating the H^- series to the H^+ series;
     equal to minus the prefactor of ST Eq. (153) *)
  Q = prefacUp[s, \[Epsilon], \[Kappa], \[Tau], \[Nu], zz] Exp[-2 I (zz - \[Eta] Log[2 zz] - \[Nu] \[Pi]/2)] (-2 I zz)^(-\[Nu] - 1 - s + I \[Epsilon]) (2 I zz)^(\[Nu] + 1 - s + I \[Epsilon]);
  (* at a purely imaginary frequency -2 I zhat is negative real and the power (-2 I zz)^(...) must be taken
     on the side of its cut continuous with Re epsilon > 0, i.e. with Arg = -Pi rather than the principal +Pi
     (see hypergeometricU) *)
  If[Re[zhat] == 0 && Im[zhat] < 0, Q = Q Exp[-2 Pi I (-\[Nu] - 1 - s + I \[Epsilon])]];
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
(*representation converges best where the series in x (Eq. (120)) is worst. Same normalisation as the series in x,*)
(*i.e. Eq. (116) divided by Gamma[1 - s - 2 I epsilon_+] (see the note on the hypergeometric functions).*)


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
    coef[n_] := Gamma[2 n + 2 nu + 1]/(Gamma[n + nu + 1 - I \[Tau]] Gamma[n + nu + 1 - s - I \[Epsilon]]) fn[q, \[Epsilon], \[Kappa], \[Tau], nu, \[Lambda], s, m, n];
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


(* The Coulomb-type representation is used where the hypergeometric series would lose more than about
   0.8 * threshold digits to cancellation. One pass of the Coulomb-type series costs two to four times one
   pass of the hypergeometric series at the same working precision (its Tricomi functions against the
   series' Gauss functions), independently of the radius, while the series' cost grows with the padding its
   loss requires; measured at 32 digits for s = -2, 0, 2 and omega = 0.1 to 2, the two cost the same at
   omega (r - r+) between 10 and 35, and below that the series is up to four times faster. The threshold
   used to be 2, chosen for accuracy alone. *)
$MSTRepresentationThreshold = 20;
(* representation of the "In" solution beyond the threshold: "Coulomb" (ST Eq. (166), series of Coulomb wave
   functions) or "Hypergeometric" (ST Eq. (138), series of hypergeometric functions in 1/(1-x)) *)
$MSTInLargeRadiusRepresentation = "Coulomb";
$masterFunction = MST`$MasterFunction;   (* captured at load time; MST`$MasterFunction is only set while the package loads *)
$radialFunctionSymbol = Symbol[MST`$MasterFunction <> "`" <> MST`$MasterFunction <> "RadialFunction"];   (* carries the messages *)

With[{sym = $radialFunctionSymbol},
  sym::prec = "The MST series for the `1` radial function at r = `2` could only be evaluated to a precision of `3` (`4` requested).";
  sym::upref = "The \"Up\" reflection amplitude is only available for real frequencies (\[Omega] = `1` given) and is Indeterminate.";
];

(* The Coulomb-type representation is valid at complex frequencies too, now that its series are taken
   on the right side of their branch cuts (see hypergeometricU and the conjugation rules of MSTRadialIn). *)
mstInRepresentation[q_, \[Epsilon]_, r_] :=
 If[$masterFunction === "Teukolsky" && Abs[\[Epsilon]] (r - (1 + Sqrt[1 - q^2]))/2 > $MSTRepresentationThreshold, $MSTInLargeRadiusRepresentation, "Series"];

mstInCore[rep_] := Switch[rep, "Coulomb", mstRadialInCoulomb, "Hypergeometric", mstRadialInLargeRadiusSeries, _, mstRadialInSeriesCore];

(* The eigenvalue and the renormalized angular momentum at the precision pp of a padded evaluation. The MST
   series amplify an error in nu (or lambda) by roughly the number of digits they lose to cancellation, so
   padding the digits of a nu known to fewer digits than pp (which SetPrecision would do) gives a wrong
   result whose tracked precision does not show it. Instead the eigenvalue is recomputed at pp and nu is
   recomputed from it, padded until it carries pp digits, on the same representative as the given nu. The
   results are cached per parameter set and precision. For a master function other than Teukolsky the
   eigenvalue is used as given. *)
$refinedParameterCache = <||>;

eigenvalueAt[s_, l_, m_, q_, \[Epsilon]_, \[Lambda]_, pp_] :=
 Module[{res = $Failed},
  Do[
    res = Quiet[Check[SpinWeightedSpheroidalEigenvalue[s, l, m, SetPrecision[q, p1] SetPrecision[\[Epsilon], p1]/2], $Failed,
      {FindRoot::cvmit, SpinWeightedSpheroidalEigenvalue::findroot}], {General::munfl, FindRoot::cvmit, FindRoot::lstol, SpinWeightedSpheroidalEigenvalue::findroot}];
    If[NumericQ[res], Return[SetPrecision[res, pp], Module]],
    {p1, {pp, pp + 16, pp + 40}}];
  \[Lambda]
 ];

refinedParameters[s_, l_, m_, q_, \[Epsilon]_, \[Lambda]_, \[Nu]_, pp_] :=
 Module[{key = {s, l, m, q, \[Epsilon], \[Lambda], \[Nu], pp}, \[Lambda]p = \[Lambda], \[Nu]p = \[Nu], res},
  If[Precision[\[Lambda]] >= pp && Precision[\[Nu]] >= pp, Return[SetPrecision[{\[Lambda], \[Nu]}, pp]]];
  res = Lookup[$refinedParameterCache, Key[key], None];
  If[res =!= None, Return[res]];
  If[$masterFunction === "Teukolsky" && Precision[\[Lambda]] < pp,
    (* the eigenvalue code compares against a machine-number tolerance, which underflows at a few
       hundred digits with a harmless General::munfl *)
    (* At some padded precisions the eigenvalue solver does not converge (FindRoot::cvmit, ::findroot) and
       returns a value accurate only to machine precision but labelled with pp digits, which the precision
       tracking downstream cannot see (s = -2, l = 7, m = 6, a omega = -0.1 came out at 1e-6.5). Such a value is
       recomputed at a higher precision; if that fails too, the input eigenvalue is kept with its own precision,
       so that the evaluation built on it sees the shortfall instead. The convergence notices are not shown. *)
    \[Lambda]p = eigenvalueAt[s, l, m, q, \[Epsilon], \[Lambda], pp];
  ];
  If[Precision[\[Nu]] < pp,
    \[Nu]p = paddedNu[s, l, m, q, \[Epsilon], \[Lambda]p, \[Nu], pp];
  ];
  (* rounded down to pp; a component whose refinement fell short keeps its precision, so that the evaluations
     built on it see the shortfall (and retry or report it) instead of a value dressed up as pp digits *)
  res = If[Precision[#] >= pp, SetPrecision[#, pp], #] & /@ {\[Lambda]p, \[Nu]p};
  If[Length[$refinedParameterCache] >= 50, $refinedParameterCache = <||>];
  $refinedParameterCache[key] = res
 ];

(* nu to pp digits, from the renormalized angular momentum computed with padded inputs, on the
   representative (among +-nu + k) of the given nu *)
paddedNu[s_, l_, m_, q_, \[Epsilon]_, \[Lambda]_, \[Nu]_, pp_] :=
 Module[{p = pp, res, tries = 0, ramAt, lamAt},
  (* a failure at a precision the retries will raise is not reported; a final failure returns the given nu,
     whose precision then shows the shortfall *)
  (* the eigenvalue is recomputed at p1 rather than having its precision raised, which would leave it accurate
     only to the precision it was computed at, a deficit that nu amplifies (ten digits at omega = 5 for l = 4);
     where the eigenvalue solver fails at such a precision (hundreds of digits at l = 36) the raised value is
     used, as before *)
  lamAt[p1_] := If[$masterFunction =!= "Teukolsky" || Precision[\[Lambda]] >= p1, SetPrecision[\[Lambda], p1],
    Quiet[Check[SpinWeightedSpheroidalEigenvalue[s, l, m, SetPrecision[q, p1] SetPrecision[\[Epsilon], p1]/2], SetPrecision[\[Lambda], p1], {FindRoot::cvmit, SpinWeightedSpheroidalEigenvalue::findroot}]]];
  ramAt[p1_] := Quiet[RenormalizedAngularMomentum[s, l, m, SetPrecision[q, p1], SetPrecision[\[Epsilon], p1]/2, lamAt[p1]], RenormalizedAngularMomentum::conv];
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
(* The padding is keyed by the parameters the series are evaluated with: for Re epsilon < 0 those of the
   conjugate partner (-m, -Conjugate[epsilon]), to which MSTRadialIn, MSTRadialUp and Amplitudes map such a
   mode (see conjugateRules), so that a padding set by the Wronskian check of TeukolskyRadial under the
   original key is seen by the evaluations. *)
paddingKey[s_, l_, m_, q_, \[Epsilon]_] := If[Re[\[Epsilon]] < 0, {s, l, -m, q, -Conjugate[\[Epsilon]]}, {s, l, m, q, \[Epsilon]}];
modePadding[s_, l_, m_, q_, \[Epsilon]_] := Lookup[$modePadding, Key[paddingKey[s, l, m, q, \[Epsilon]]], 0];
setModePadding[s_, l_, m_, q_, \[Epsilon]_, extra_] := (If[Length[$modePadding] >= 50, $modePadding = <||>]; $modePadding[paddingKey[s, l, m, q, \[Epsilon]]] = extra);

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
    With[{sym = $radialFunctionSymbol}, Message[sym::prec, If[MemberQ[{mstRadialUpSeriesCore, mstRadialUpHorizon}, core], "Up", "In"], r, If[allNumericQ[res], Precision[res], res], target]];
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

(* Pre-padding from the loss of earlier evaluations. mstPaddedEvaluation starts at the working precision
   and, once it has measured the loss, repeats the evaluation with the padding that loss requires, so an
   evaluation costs two passes, the first of them wasted, whenever the loss exceeds a digit, which is
   always. The loss is predictable: for the hypergeometric "In" series it is about 0.8 |epsilon| (r - r+)/2
   plus a mode-dependent part, and for the Coulomb-type "In" representation and the "Up" series it does
   not depend on the radius. The mode-dependent part measured at every evaluation is cached per mode and
   core, and the next evaluation starts at the predicted padding, with the retry loop still there should
   the prediction fall short. *)
$lossCache = <||>;

seriesLossEstimate[q_, \[Epsilon]_, r_] := 0.8 Abs[\[Epsilon]] (r - (1 + Sqrt[1 - q^2]))/2;

prepaddedEvaluation[core_, params:{s_, l_, m_, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_}, goals:{wp_, prec_, acc_}, deriv_, r_] :=
 Module[{key = {core, s, l, m, q, \[Epsilon], \[Nu], \[Lambda]}, rdep, base, p0, res},
  rdep = If[MemberQ[{mstRadialInSeriesCore, mstRadialUpHorizon}, core], seriesLossEstimate[q, \[Epsilon], r], 0];
  base = Lookup[$lossCache, Key[key], None];
  (* four digits of margin: significance arithmetic overstates the precision of the summed series by a
     digit or so, and the margin keeps the result beyond the requested precision as before *)
  p0 = If[base === None, Automatic, If[wp === MachinePrecision, $MachinePrecision, wp] + base + rdep + 4];
  res = mstPaddedEvaluation[core, params, goals, deriv, r, 4, p0];
  If[NumericQ[$lastPaddingLoss],
    If[Length[$lossCache] >= 50, $lossCache = <||>];
    $lossCache[key] = Max[$lastPaddingLoss - rdep, 0]];
  res
 ];

mstRadialInEvaluate[params:{s_, l_, m_, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_}, goals_, deriv_, r_] :=
 Module[{rep = mstInRepresentation[q, \[Epsilon], r], key = {s, l, m, q, \[Epsilon], \[Nu], \[Lambda]}, seriesLoss, cached, res},
  If[rep === "Series", Return[prepaddedEvaluation[mstRadialInSeriesCore, params, goals, deriv, r]]];
  seriesLoss = seriesLossEstimate[q, \[Epsilon], r] + 8;
  cached = Lookup[$inRepresentationCache, Key[key], None];
  If[cached =!= None && cached > seriesLoss, Return[prepaddedEvaluation[mstRadialInSeriesCore, params, goals, deriv, r]]];
  If[cached =!= None, Return[prepaddedEvaluation[mstInCore[rep], params, goals, deriv, r]]];
  (* first evaluation of this mode: a single pass, to measure the loss *)
  res = mstPaddedEvaluation[mstInCore[rep], params, goals, deriv, r, 1];
  If[Length[$inRepresentationCache] >= 50, $inRepresentationCache = <||>];
  $inRepresentationCache[key] = $lastPaddingLoss;
  If[$lastPaddingLoss > seriesLoss, Return[prepaddedEvaluation[mstRadialInSeriesCore, params, goals, deriv, r]]];
  (* the precision reached is judged from the tracked loss, not from res: at machine precision
     mstPaddedEvaluation returns N[res], whose MachinePrecision says nothing about the digits it carries *)
  If[!allNumericQ[res] || $lastPaddingPrecision - $lastPaddingLoss < If[goals[[1]] === MachinePrecision, $MachinePrecision, goals[[1]]] - 1,
    res = prepaddedEvaluation[mstInCore[rep], params, goals, deriv, r];
  ];
  res
 ];

(* Frequencies with Re epsilon < 0: the radial equation and its boundary conditions are invariant under
   complex conjugation combined with (m, epsilon) -> (-m, -Conjugate[epsilon]), which maps the solutions
   normalised to unit transmission onto each other, so R[m, epsilon] = Conjugate[R[-m, -Conjugate[epsilon]]]
   with conjugate eigenvalue, nu and amplitudes. The Coulomb-type series and the amplitude formulae are
   derived for Re epsilon > 0 (their powers and logarithms of zhat and epsilon cross branch cuts otherwise,
   which gave "Up" solutions violating this symmetry by O(1) for Re omega < 0), so every evaluation with
   Re epsilon < 0 is mapped to its partner. *)
conjugateRules[f_Symbol] := (
  f[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_, {wp_, prec_, acc_}][r_?NumericQ] /; Re[\[Epsilon]] < 0 :=
    Conjugate[f[s, l, -m, q, -Conjugate[\[Epsilon]], Conjugate[\[Nu]], Conjugate[\[Lambda]], Conjugate[norm], {wp, prec, acc}][r]];
  Derivative[1][f[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_, {wp_, prec_, acc_}]][r_?NumericQ] /; Re[\[Epsilon]] < 0 :=
    Conjugate[Derivative[1][f[s, l, -m, q, -Conjugate[\[Epsilon]], Conjugate[\[Nu]], Conjugate[\[Lambda]], Conjugate[norm], {wp, prec, acc}]][r]];
  f[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_, {wp_, prec_, acc_}][r_?NumericQ, {0, 1}] /; Re[\[Epsilon]] < 0 :=
    Conjugate[f[s, l, -m, q, -Conjugate[\[Epsilon]], Conjugate[\[Nu]], Conjugate[\[Lambda]], Conjugate[norm], {wp, prec, acc}][r, {0, 1}]];
);
conjugateRules /@ {MSTRadialIn, MSTRadialUp};

(* cores taking a derivative order, wrapping the series definitions above *)
(* Public evaluation: value, first derivative, or both from one summation (f[r, {0, 1}]) *)
MSTRadialIn[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_, {wp_, prec_, acc_}][r_?NumericQ] :=
 mstRadialInEvaluate[{s, l, m, q, \[Epsilon], \[Nu], \[Lambda], norm}, {wp, prec, acc}, 0, r];

Derivative[1][MSTRadialIn[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_, {wp_, prec_, acc_}]][r_?NumericQ] :=
 mstRadialInEvaluate[{s, l, m, q, \[Epsilon], \[Nu], \[Lambda], norm}, {wp, prec, acc}, 1, r];

MSTRadialIn[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_, {wp_, prec_, acc_}][r_?NumericQ, {0, 1}] :=
 mstRadialInEvaluate[{s, l, m, q, \[Epsilon], \[Nu], \[Lambda], norm}, {wp, prec, acc}, All, r];

(* "Up" solution: the horizon basis near the horizon (mstUpRepresentation), the Coulomb-type series
   elsewhere and wherever the horizon representation is not available *)
mstRadialUpEvaluate[params:{s_, l_, m_, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_}, goals_, deriv_, r_] :=
 Module[{res},
  If[mstUpRepresentation[s, m, q, \[Epsilon], r] === "Horizon",
    res = prepaddedEvaluation[mstRadialUpHorizon, params, goals, deriv, r];
    If[allNumericQ[res], Return[res]]];
  prepaddedEvaluation[mstRadialUpSeriesCore, params, goals, deriv, r]
 ];

MSTRadialUp[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_, {wp_, prec_, acc_}][r_?NumericQ] :=
 mstRadialUpEvaluate[{s, l, m, q, \[Epsilon], \[Nu], \[Lambda], norm}, {wp, prec, acc}, 0, r];

Derivative[1][MSTRadialUp[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_, {wp_, prec_, acc_}]][r_?NumericQ] :=
 mstRadialUpEvaluate[{s, l, m, q, \[Epsilon], \[Nu], \[Lambda], norm}, {wp, prec, acc}, 1, r];

MSTRadialUp[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_, {wp_, prec_, acc_}][r_?NumericQ, {0, 1}] :=
 mstRadialUpEvaluate[{s, l, m, q, \[Epsilon], \[Nu], \[Lambda], norm}, {wp, prec, acc}, All, r];


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
