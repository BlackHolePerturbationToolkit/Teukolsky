(* ::Package:: *)

BeginPackage["Teukolsky`MST`"];

(*MSTRadialIn::usage = "MSTRadialIn[s, l, m, q, \[Epsilon], \[Nu], \[Lambda]][r]";
MSTRadialUp::usage = "MSTRadialUp[s, l, m, q, \[Epsilon], \[Nu], \[Lambda]][r]";*)

Begin["`Private`"];

(******************************************************************************)
(***************************** Utility functions ******************************)
(******************************************************************************)

(* Continued fraction with automatic convergence checking. *)
(* FIXME: There is a potentially better algorithm in Numerical recipes which handles
          cases where successive terms have large magnitude differences *)
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


(******************************************************************************)
(************************** Hypergeometric functions **************************)
(******************************************************************************)

H2F1Exact[n_, s_, \[Nu]_, \[Tau]_, \[Epsilon]_, x_] := 
  Hypergeometric2F1[n + \[Nu] + 1 - I \[Tau], -n - \[Nu] - I \[Tau], 1 - s - I (\[Epsilon] + \[Tau]), x];

H2F1Up[n_, s_, \[Nu]_, \[Tau]_, \[Epsilon]_, x_] :=
  (n + \[Nu])/((n + \[Nu] - (s + I \[Epsilon])) ((n + \[Nu]) - I \[Tau])) {
    -((((n + \[Nu]) - 1 + (s + I \[Epsilon])) ((n + \[Nu]) - 1 + I \[Tau]) )/((n + \[Nu]) - 1) ) H2F1[n - 2],
    (2 (n + \[Nu]) - 1) (1 - 2 x + (I (s + I \[Epsilon]) \[Tau])/(((n + \[Nu]) - 1) (n + \[Nu]))) H2F1[n - 1]
  };

H2F1Down[n_, s_, \[Nu]_, \[Tau]_, \[Epsilon]_, x_] :=
  (1 + n + \[Nu])/((1 + n + s + I \[Epsilon] + \[Nu]) (1 + n + \[Nu] + I \[Tau])) {
    (2 (n + \[Nu]) + 3) (1 - 2 x + (I (s + I \[Epsilon]) \[Tau])/(((n + \[Nu]) + 1) ((n + \[Nu]) + 2))) H2F1[n + 1],
    -((((n + \[Nu]) + 2 - (s + I \[Epsilon])) ((n + \[Nu]) + 2 - I \[Tau]))/(2 + n + \[Nu]) ) H2F1[n + 2]
  };

(* FIXME: There's probably a better way to write this using identities for derivaties of 2F1 *)
dH2F1Exact[n_, s_, \[Nu]_, \[Tau]_, \[Epsilon]_, x_] := 
  Derivative[0, 0, 0, 1][Hypergeometric2F1][n + \[Nu] + 1 - I \[Tau], -n - \[Nu] - I \[Tau], 1 - s - I (\[Epsilon] + \[Tau]), x];

dH2F1Up[n_, s_, \[Nu]_, \[Tau]_, \[Epsilon]_, x_] :=
  (n + \[Nu])/((n + \[Nu] - (s + I \[Epsilon])) ((n + \[Nu]) - I \[Tau])) {
   -(((-1 + n + s + I \[Epsilon] + \[Nu]) (-1 + n + \[Nu] + I \[Tau]) dH2F1[n - 2])/(-1 + n + \[Nu])),
   (-1 + 2 (n + \[Nu])) (-2 H2F1[n - 1]),
   (-1 + 2 (n + \[Nu])) (1 - 2 x + (I (s + I \[Epsilon]) \[Tau])/((-1 + n + \[Nu]) (n + \[Nu]))) dH2F1[n - 1]};

dH2F1Down[n_, s_, \[Nu]_, \[Tau]_, \[Epsilon]_, x_] :=
  (1 + n + \[Nu])/((1 + n + s + I \[Epsilon] + \[Nu]) (1 + n + \[Nu] + I \[Tau])) {
    (3 + 2 (n + \[Nu])) (-2 H2F1[n + 1]),
    (3 + 2 (n + \[Nu])) (1 - 2 x + (I (s + I \[Epsilon]) \[Tau])/((1 + n + \[Nu]) (2 + n + \[Nu]))) dH2F1[n + 1],
    -(((2 + n - s - I \[Epsilon] + \[Nu]) (2 + n + \[Nu] - I \[Tau]) dH2F1[n + 2])/(2 + n + \[Nu]))};

HUExact[n_, s_, \[Nu]_, \[Epsilon]_, zhat_] :=
  E^(I zhat) (-2 I zhat)^(n + \[Nu] + 1) HypergeometricU[n + \[Nu] + 1 + s - I \[Epsilon], 2 (n + \[Nu]) + 2, -2 I zhat];

HUUp[n_, s_, \[Nu]_, \[Epsilon]_, zhat_] :=
  1/(((n + \[Nu]) - 1) ((n + \[Nu]) + s - I \[Epsilon])) {
    (n + \[Nu]) ((n + \[Nu]) - 1 - (s - I \[Epsilon])) HU[n - 2],
    (I (2 (n + \[Nu]) - 1) ((n + \[Nu])^2 - (n + \[Nu]) - zhat (\[Epsilon] + I s)))/zhat HU[n - 1]
  };

HUDown[n_, s_, \[Nu]_, \[Epsilon]_, zhat_] :=
  {-((I (2 (n + \[Nu]) + 3) ((n + \[Nu])^2 + 3 (n + \[Nu]) + 2 - zhat (\[Epsilon] + I s)))/(((n + \[Nu]) + 2) zhat ((n + \[Nu]) + 1 - (s - I \[Epsilon])))) HU[n + 1],
   (((n + \[Nu]) + 1) ((n + \[Nu]) + 2 + (s - I \[Epsilon])))/(((n + \[Nu]) + 2) ((n + \[Nu]) + 1 - (s - I \[Epsilon]))) HU[n + 2]
  };

dHUExact[n_, s_, \[Nu]_, \[Epsilon]_, zhat_] :=
 2^(1 + n + \[Nu]) E^(I zhat) (-I zhat)^(n + \[Nu]) (-I (1 + n + I zhat + \[Nu]) HypergeometricU[n + \[Nu] + 1 + s - I \[Epsilon], 2 (n + \[Nu]) + 2, -2 I zhat] - 2 zhat Derivative[0,0,1][HypergeometricU][n + \[Nu] + 1 + s - I \[Epsilon], 2 (n + \[Nu]) + 2, -2 I zhat]);

dHUUp[n_, s_, \[Nu]_, \[Epsilon]_, zhat_] :=
 1/(((n + \[Nu]) - 1) ((n + \[Nu]) + s - I \[Epsilon])) {
   (n + \[Nu]) ((n + \[Nu]) - 1 - (s - I \[Epsilon])) dHU[n - 2],
   I (-1 + 2 n + 2 \[Nu]) (n - n^2 + \[Nu] - 2 n \[Nu] - \[Nu]^2) HU[n - 1] / zhat^2,
   I (-1 + 2 n + 2 \[Nu]) (n^2 - I s zhat - zhat \[Epsilon] - \[Nu] + \[Nu]^2 + n (-1 + 2 \[Nu])) dHU[n - 1] / zhat
 };

dHUDown[n_, s_, \[Nu]_, \[Epsilon]_, zhat_] :=
 1 / ((2 + n + \[Nu]) (1 + n - s + I \[Epsilon] + \[Nu])) {
   -I (3 + 2 (n + \[Nu])) / zhat^2 (-(2 + n^2 + 3 \[Nu] + \[Nu]^2 + n (3 + 2 \[Nu])) HU[n + 1]),
   -I (3 + 2 (n + \[Nu])) / zhat (2 + n^2 - I s zhat - zhat \[Epsilon] + 3 \[Nu] + \[Nu]^2 + n (3 + 2 \[Nu])) dHU[n + 1],
   (1 + n + \[Nu]) (2 + n + s - I \[Epsilon] + \[Nu]) dHU[n + 2]};

(******************************************************************************)
(************************** MST series coefficients ***************************)
(******************************************************************************)

\[Alpha][q_, \[Epsilon]_, \[Kappa]_, \[Tau]_, \[Nu]_, \[Lambda]_, s_, m_, n_] :=
 \[Alpha][q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, n] =
  (I \[Epsilon] \[Kappa] (n + \[Nu] + 1 + s + I \[Epsilon]) (n + \[Nu] + 1 + s - I \[Epsilon]) (n + \[Nu] + 1 + I \[Tau]))/((n + \[Nu] + 1) (2 n + 2 \[Nu] + 3));

\[Beta][q_, \[Epsilon]_, \[Kappa]_, \[Tau]_, \[Nu]_, \[Lambda]_, s_, m_, n_] :=
 \[Beta][q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, n] =
  -\[Lambda] - s (s + 1) + (n + \[Nu]) (n + \[Nu] + 1) + \[Epsilon]^2 + \[Epsilon] (\[Epsilon] - m q) + (\[Epsilon] (\[Epsilon] - m q) (s^2 + \[Epsilon]^2))/((n + \[Nu]) (n + \[Nu] + 1));

\[Gamma][q_, \[Epsilon]_, \[Kappa]_, \[Tau]_, \[Nu]_, \[Lambda]_, s_, m_, n_] :=
 \[Gamma][q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, n] =
  -((I \[Epsilon] \[Kappa] (n + \[Nu] - s + I \[Epsilon]) (n + \[Nu] - s - I \[Epsilon]) (n + \[Nu] - I \[Tau]))/((n + \[Nu]) (2 n + 2 \[Nu] - 1)));

f[q_, \[Epsilon]_, \[Kappa]_, \[Tau]_, \[Nu]_, \[Lambda]_, s_, m_, 0] = 1;

f[q_, \[Epsilon]_, \[Kappa]_, \[Tau]_, \[Nu]_, \[Lambda]_, s_, m_, nf_] :=
 f[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, nf] =
 Module[{\[Alpha]n, \[Beta]n, \[Gamma]n, i, n, ret},
  \[Alpha]n[n_] := \[Alpha][q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, n];
  \[Beta]n[n_] := \[Beta][q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, n];
  \[Gamma]n[n_] := \[Gamma][q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, n];

  If[nf > 0,
    ret = f[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, nf - 1] CF[-\[Alpha]n[i - 1] \[Gamma]n[i], \[Beta]n[i], {i, nf}]/\[Alpha]n[nf - 1];
  ,
    ret = f[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, nf + 1] CF[-\[Alpha]n[2 nf - i] \[Gamma]n[2 nf - i + 1], \[Beta]n[2 nf - i], {i, nf}]/\[Gamma]n[nf + 1];
  ];
  
  ret
];

(******************************************************************************)
(*************************** Asymptotic amplitudes ****************************)
(******************************************************************************)

Amplitudes[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_] :=
 Module[{\[Kappa], \[Tau], \[Epsilon]p, \[Omega], K\[Nu], K\[Nu]1, K\[Nu]2, Aplus, Btrans, Ctrans, Binc, Bref, n, fSumUp, fSumDown, fSumK\[Nu]1Up, fSumK\[Nu]1Down, fSumK\[Nu]2Up, fSumK\[Nu]2Down, fSumCUp, fSumCDown},
 Internal`InheritedBlock[{\[Alpha], \[Beta], \[Gamma], f},
  \[Kappa] = Sqrt[1 - q^2];
  \[Tau] = (\[Epsilon] - m q)/\[Kappa];
  \[Epsilon]p = 1/2 (\[Tau] + \[Epsilon]);
  \[Omega] = \[Epsilon] / 2;

  (* All of the formulae are taken from Sasaki & Tagoshi, Living Rev. Relativity 6:6 *)

  (* There are three formally infinite sums which must be computed, but which may be numerically 
     truncated after a finite number of terms. We determine how many terms to include by summing
     until the result doesn't change. *)

  (* Sum MST series coefficients with no extra factors*)
  fSumUp = fSumDown = 0;

  n = 0;
  While[fSumUp != (fSumUp += 
    f[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, n]),
    n++;
  ];

  n = -1;
  While[fSumDown != (fSumDown += 
    f[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, n]), 
    n--;
  ];
  
  (* Sums appearing in ST Eq. (165). We evaluate these with 1: \[Nu] and 2:-\[Nu]-1 *)
  fSumK\[Nu]1Up = fSumK\[Nu]1Down = 0;

  n = 0;
  While[fSumK\[Nu]1Up != (fSumK\[Nu]1Up += 
    ((-1)^n Gamma[1 + n + s + I \[Epsilon] + \[Nu]] Gamma[1 + n + 2 \[Nu]] Gamma[1 + n + \[Nu] + I \[Tau]])/(n! Gamma[1 + n - s - I \[Epsilon] + \[Nu]] Gamma[1 + n + \[Nu] - I \[Tau]]) f[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, n]),
    n++;
  ];

  n = 0;
  While[fSumK\[Nu]1Down != (fSumK\[Nu]1Down += 
    (((-1)^n) Pochhammer[1 + s - I \[Epsilon] + \[Nu], n])/((-n)! Pochhammer[1 - s + I \[Epsilon] + \[Nu], n] Pochhammer[2 + 2 \[Nu], n]) f[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, n]), 
    n--;
  ];

  fSumK\[Nu]2Up = fSumK\[Nu]2Down = 0;

  n = 0;
  While[fSumK\[Nu]2Up != (fSumK\[Nu]2Up += 
    ((-1)^n Gamma[1 + n + s + I \[Epsilon] + (-1-\[Nu])] Gamma[1 + n + 2 (-1-\[Nu])] Gamma[1 + n + (-1-\[Nu]) + I \[Tau]])/(n! Gamma[1 + n - s - I \[Epsilon] + (-1-\[Nu])] Gamma[1 + n + (-1-\[Nu]) - I \[Tau]]) f[q, \[Epsilon], \[Kappa], \[Tau], (-1-\[Nu]), \[Lambda], s, m, n]),
    n++;
  ];

  n = 0;
  While[fSumK\[Nu]2Down != (fSumK\[Nu]2Down += 
    (((-1)^n) Pochhammer[1 + s - I \[Epsilon] + (-1-\[Nu]), n])/((-n)! Pochhammer[1 - s + I \[Epsilon] + (-1-\[Nu]), n] Pochhammer[2 + 2 (-1-\[Nu]), n]) f[q, \[Epsilon], \[Kappa], \[Tau], (-1-\[Nu]), \[Lambda], s, m, n]), 
    n--;
  ];

  (* Sum appearing in ST Eq. (158) *)
  fSumCUp = fSumCDown = 0;

  n = 0;
  While[fSumCUp != (fSumCUp += 
    (-1)^n Pochhammer[\[Nu] + 1 + s - I \[Epsilon], n]/Pochhammer[\[Nu] + 1 - s + I \[Epsilon], n] f[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, n]),
    n++;
  ];

  n = -1;
  While[fSumCDown != (fSumCDown += 
    (-1)^n Pochhammer[\[Nu] + 1 + s - I \[Epsilon], n]/Pochhammer[\[Nu] + 1 - s + I \[Epsilon], n] f[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, n]), 
    n--;
  ];
  

  (* K\[Nu]: ST Eq. (165) *)
  K\[Nu]1 = ((2^-\[Nu]) (E^(I \[Epsilon] \[Kappa])) ((\[Epsilon] \[Kappa])^(s - \[Nu])) Gamma[1 - s - 2 I \[Epsilon]p] Gamma[2 + 2 \[Nu]])/(Gamma[1 - s + I \[Epsilon] + \[Nu]] Gamma[1 + s + I \[Epsilon] + \[Nu]] Gamma[1 + \[Nu] + I \[Tau]]) fSumK\[Nu]1Up / fSumK\[Nu]1Down;

  K\[Nu]2 = ((2^-(-1-\[Nu])) (E^(I \[Epsilon] \[Kappa])) ((\[Epsilon] \[Kappa])^(s - (-1-\[Nu]))) Gamma[1 - s - 2 I \[Epsilon]p] Gamma[2 + 2 (-1-\[Nu])])/(Gamma[1 - s + I \[Epsilon] + (-1-\[Nu])] Gamma[1 + s + I \[Epsilon] + (-1-\[Nu])] Gamma[1 + (-1-\[Nu]) + I \[Tau]]) fSumK\[Nu]2Up / fSumK\[Nu]2Down;

  (* A+: ST Eq. (157) *)
  Aplus = (2^(-1 + s - I \[Epsilon]) E^(-((\[Pi] \[Epsilon])/2) + 1/2 I \[Pi] (1 - s + \[Nu])) Gamma[1 - s + I \[Epsilon] + \[Nu]])/Gamma[1 + s - I \[Epsilon] + \[Nu]] (fSumUp+fSumDown);

  (* Btrans: ST Eq. (167) *)
  Btrans = 4^s \[Kappa]^(2 s) E^(I (\[Epsilon] + \[Tau]) \[Kappa] (1/2 + Log[\[Kappa]]/(1 + \[Kappa]))) (fSumUp+fSumDown);

  (* Binc: ST Eq. (168) *)
  Binc = \[Omega]^-1 (K\[Nu]1 - I E^(-I \[Pi] \[Nu]) Sin[\[Pi] (\[Nu] - s + I \[Epsilon])] / Sin[\[Pi] (\[Nu] + s - I \[Epsilon])] K\[Nu]2) Exp[-I \[Epsilon] (Log[\[Epsilon]] - (1 - \[Kappa])/2)] Aplus;

  (* Bref: ST Eq. (169) *)
  Bref = Ctrans (K\[Nu]1 + I E^(I \[Pi] \[Nu]) K\[Nu]2);

  (* Ctrans: ST Eqs. (158) and (170) *)
  Ctrans = (\[Epsilon]/2)^(-1 - 2 s) Exp[I \[Epsilon] (Log[\[Epsilon]] - (1 - \[Kappa])/2)] 2^(-1 - s + I \[Epsilon]) Exp[-\[Pi] (\[Epsilon] + I (\[Nu] + 1 + s))/2] (fSumCUp+fSumCDown);

  (* Return results as an Association *)
  <| "In" -> <| "Incidence" -> Binc, "Transmission" -> Btrans, "Reflection" -> Bref|>,
     "Up" -> <| "Transmission" -> Ctrans |>
   |>
]];

(******************************************************************************)
(****************************** Radial functions ******************************)
(******************************************************************************)

(* Throwe B.1, Sasaki & Tagoshi Eqs. (116) and (120) *)
SetAttributes[MSTRadialIn, {NumericFunction}];

MSTRadialIn[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_][r_?NumericQ] :=
 Module[{\[Kappa], \[Tau], rp, x, resUp, nUp, resDown, nDown},
 Block[{H2F1},
 Internal`InheritedBlock[{\[Alpha], \[Beta], \[Gamma], f},
  \[Kappa] = Sqrt[1 - q^2];
  \[Tau] = (\[Epsilon] - m q)/\[Kappa];
  rp = 1 + Sqrt[1 - q^2];
  x = (r - rp)/(2 \[Kappa]);

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
 
  resUp = resDown = 0;

  nUp = 0;
  While[resUp != (resUp += 
    f[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, nUp] H2F1[nUp]),
    nUp++;
  ];

  nDown = -1;
  While[resDown != (resDown += 
    f[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, nDown] H2F1[nDown]), 
    nDown--;
  ];

  E^(I \[Epsilon] \[Kappa] x) (-x)^(-s - I (\[Epsilon] + \[Tau])/2) (1 - x)^(I (\[Epsilon] - \[Tau])/2) (resUp + resDown) / norm
]]];

Derivative[1][MSTRadialIn[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_]][r_?NumericQ] :=
 Module[{\[Kappa], \[Tau], rp, x, dxdr, resUp, nUp, resDown, nDown},
 Block[{H2F1, dH2F1},
 Internal`InheritedBlock[{\[Alpha], \[Beta], \[Gamma], f},
  \[Kappa] = Sqrt[1 - q^2];
  \[Tau] = (\[Epsilon] - m q)/\[Kappa];
  rp = 1 + Sqrt[1 - q^2];
  x = (r - rp)/(2 \[Kappa]);
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
 
  resUp = resDown = 0;

  nUp = 0;
  While[resUp != (resUp += 
    f[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, nUp] ((-2 s (-1 + x) + I (\[Epsilon] - 2 x \[Epsilon] \[Kappa] + 2 x^2 \[Epsilon] \[Kappa] + \[Tau] - 2 x \[Tau])) H2F1[nUp] + 2 (-1 + x) x dH2F1[nUp])),
    nUp++;
  ];

  nDown = -1;
  While[resDown != (resDown += 
    f[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, nDown] ((-2 s (-1 + x) + I (\[Epsilon] - 2 x \[Epsilon] \[Kappa] + 2 x^2 \[Epsilon] \[Kappa] + \[Tau] - 2 x \[Tau])) H2F1[nDown] + 2 (-1 + x) x dH2F1[nDown])),
    nDown--;
  ];

  1/2 E^(I x \[Epsilon] \[Kappa]) (1 - x)^(1/2 I (2 I + \[Epsilon] - \[Tau])) (-x)^(-1 - s - 1/2 I (\[Epsilon] + \[Tau])) (resUp + resDown) dxdr / norm
]]];

Derivative[n_Integer?Positive][MSTRadialIn[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_]][r0_?NumericQ] :=
 Module[{d2R, Rderivs, R, r, i},
 (*FIXME: Add appropriate factors of M here *)
  d2R = (-(-\[Lambda] + 2 I r s \[Epsilon] + (-2 I (-1 + r) s (-q m + (q^2 + r^2) \[Epsilon]/2) + (-q m + (q^2 + r^2) \[Epsilon]/2)^2)/(q^2 - 2 r + r^2)) R[r] - (-2 + 2 r) (1 + s) Derivative[1][R][r])/(q^2 - 2 r + r^2);

  pderivs = D[R[r_], {r_, i_}] :> D[d2R, {r, i - 2}] /; i >= 2;
  Do[Derivative[i][R][r] = Simplify[D[Derivative[i - 1][R][r], r] /. pderivs];, {i, 2, n}];
  Derivative[n][R][r] /. {
    R'[r] -> MSTRadialIn[s, l, m, q, \[Epsilon], \[Nu], \[Lambda], norm]'[r0],
    R[r] -> MSTRadialIn[s, l, m, q, \[Epsilon], \[Nu], \[Lambda], norm][r0], r -> r0}
];

(* Throwe B.5, Sasaki & Tagoshi (153) and (159) *)
SetAttributes[MSTRadialUp, {NumericFunction}];

MSTRadialUp[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_][r_?NumericQ]  :=
 Module[{\[Kappa], \[Tau], \[Epsilon]p, rm, z, zm, zhat, resUp, nUp, resDown, nDown},
 Block[{HU},
 Internal`InheritedBlock[{\[Alpha], \[Beta], \[Gamma], f},
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

  resDown = resUp = 0;
  nUp = 0;

  While[resUp != (resUp +=
    I^nUp Pochhammer[\[Nu] + 1 + s - I \[Epsilon], nUp]/Pochhammer[\[Nu] + 1 - s + I \[Epsilon], nUp] f[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, nUp](2zhat)^nUp HU[nUp]/(E^(I zhat) (-2 I zhat)^(nUp + \[Nu] + 1) )),
    nUp++;
  ];
  
  nDown = -1;
  While[resDown != (resDown +=
    I^nDown Pochhammer[\[Nu] + 1 + s - I \[Epsilon], nDown]/Pochhammer[\[Nu] + 1 - s + I \[Epsilon], nDown] f[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, nDown](2zhat)^nDown HU[nDown]/(E^(I zhat) (-2 I zhat)^(nDown + \[Nu] + 1) )),
    nDown--;
  ];

  2^\[Nu] E^(-\[Pi] \[Epsilon]) E^(-I \[Pi](\[Nu]+1+s)) E^(I zhat) zhat^(\[Nu]+I \[Epsilon]p) (zhat-\[Epsilon] \[Kappa])^(-s-I \[Epsilon]p) (resUp + resDown) / norm
]]];

Derivative[1][MSTRadialUp[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_]][r_?NumericQ] :=
 Module[{\[Kappa], \[Tau], \[Epsilon]p, rm, z, zm, zhat, dzhatdr, resUp, nUp, resDown, nDown},
 Block[{HU, dHU},
 Internal`InheritedBlock[{\[Alpha], \[Beta], \[Gamma], f},
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
    res = t1 + t2 + t3;
    If[Max[Abs[{t1, t2, t3}/res]] > 2.,
      res = dHUExact[n, s, \[Nu], \[Epsilon], zhat];
    ];
    res
  ];

  resDown = resUp = 0;
  nUp = 0;

  While[resUp != (resUp +=
    I^nUp (-I zhat)^(1 - nUp - \[Nu]) zhat^(-3 + nUp + I \[Epsilon]p + \[Nu]) f[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, nUp] Pochhammer[1 + s - I \[Epsilon] + \[Nu], nUp] (((1 + s) zhat + I \[Epsilon] (I + \[Epsilon]p) \[Kappa]) HU[nUp] + zhat (-zhat + \[Epsilon] \[Kappa]) dHU[nUp]) / (2 Pochhammer[1 - s + I \[Epsilon] + \[Nu], nUp])),
    nUp++;
  ];
  
  nDown = -1;
  While[resDown != (resDown +=
    I^nDown (-I zhat)^(1 - nDown - \[Nu]) zhat^(-3 + nDown + I \[Epsilon]p + \[Nu]) f[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, nDown] Pochhammer[1 + s - I \[Epsilon] + \[Nu], nDown] (((1 + s) zhat + I \[Epsilon] (I + \[Epsilon]p) \[Kappa]) HU[nDown] + zhat (-zhat + \[Epsilon] \[Kappa]) dHU[nDown]) / (2 Pochhammer[1 - s + I \[Epsilon] + \[Nu], nDown])),
    nDown--;
  ];

  E^(-\[Pi] (\[Epsilon] + I (1 + s + \[Nu]))) (zhat - \[Epsilon] \[Kappa])^(-1 - s - I \[Epsilon]p) (resUp + resDown) dzhatdr / norm
]]];

Derivative[n_Integer?Positive][MSTRadialUp[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_]][r0_?NumericQ] :=
 Module[{d2R, Rderivs, R, r, i},
  d2R = (-(-\[Lambda] + 2 I r s \[Epsilon] + (-2 I (-1 + r) s (-q m + (q^2 + r^2) \[Epsilon]/2) + (-q m + (q^2 + r^2) \[Epsilon]/2)^2)/(q^2 - 2 r + r^2)) R[r] - (-2 + 2 r) (1 + s) Derivative[1][R][r])/(q^2 - 2 r + r^2);

  pderivs = D[R[r_], {r_, i_}] :> D[d2R, {r, i - 2}] /; i >= 2;
  Do[Derivative[i][R][r] = Simplify[D[Derivative[i - 1][R][r], r] /. pderivs];, {i, 2, n}];
  Derivative[n][R][r] /. {
    R'[r] -> MSTRadialUp[s, l, m, q, \[Epsilon], \[Nu], \[Lambda], norm]'[r0],
    R[r] -> MSTRadialUp[s, l, m, q, \[Epsilon], \[Nu], \[Lambda], norm][r0], r -> r0}
];

End[];
EndPackage[];
