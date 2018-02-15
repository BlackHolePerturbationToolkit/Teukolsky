(* ::Package:: *)

BeginPackage["Teukolsky`MST`"];

MSTRadialIn::usage = "MSTRadialIn[M, a, \[Omega], s, l, m, \[Nu], \[Lambda], r]";
MSTRadialUp::usage = "MSTRadialUp[M, a, \[Omega], s, l, m, \[Nu], \[Lambda], r]";

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
(****************************** Radial functions ******************************)
(******************************************************************************)

(* Throwe B.1, Sasaki & Tagoshi Eqs. (116) and (120) *)
MSTRadialIn[M_, a_, \[Omega]_, s_Integer, l_Integer, m_Integer, \[Nu]_, \[Lambda]_, r_?InexactNumberQ] := 
 Module[{q, \[Epsilon], \[Kappa], \[Tau], rp, z, zp, x, resUp, nUp, resDown, nDown},
 Block[{H2F1},
 Internal`InheritedBlock[{\[Alpha], \[Beta], \[Gamma], f},
  q = a/M;
  \[Epsilon] = 2 M \[Omega];
  \[Kappa] = Sqrt[1 - q^2];
  \[Tau] = (\[Epsilon] - m q)/\[Kappa];
  rp = M + Sqrt[M^2 - a^2];
  z = \[Omega] r;
  zp = \[Omega] rp;
  x = (zp - z)/(\[Epsilon] \[Kappa]);

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

  E^(I \[Epsilon] \[Kappa] x) (-x)^(-s - I (\[Epsilon] + \[Tau])/2) (1 - x)^(I (\[Epsilon] - \[Tau])/2) (resUp + resDown)
]]];

(* Throwe B.5, Sasaki & Tagoshi (153) and (159) *)
MSTRadialUp[M_, a_, \[Omega]_, s_Integer, l_Integer, m_Integer, \[Nu]_, \[Lambda]_, r_?InexactNumberQ] := 
 Module[{q, \[Epsilon], \[Kappa], \[Tau], \[Epsilon]p, rm, z, zm, zhat, resUp, nUp, resDown, nDown},
 Block[{HU},
 Internal`InheritedBlock[{\[Alpha], \[Beta], \[Gamma], f},
  q = a/M;
  \[Epsilon] = 2 M \[Omega];
  \[Kappa] = Sqrt[1 - q^2];
  \[Tau] = (\[Epsilon] - m q)/\[Kappa];
  \[Epsilon]p = 1/2 (\[Tau]+\[Epsilon]);
  rm = M - Sqrt[M^2 - a^2];
  z = \[Omega] r;
  zm = \[Omega] rm;
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

  2^\[Nu] E^(-\[Pi] \[Epsilon]) E^(-I \[Pi](\[Nu]+1+s)) E^(I zhat) zhat^(\[Nu]+I \[Epsilon]p) (zhat-\[Epsilon] \[Kappa])^(-s-I \[Epsilon]p) (resUp + resDown)
]]];

End[];
EndPackage[];
