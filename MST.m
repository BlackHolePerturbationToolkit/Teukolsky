(* ::Package:: *)

BeginPackage["Teukolsky`MST`"];

(*MSTRadialIn::usage = "MSTRadialIn[s, l, m, q, \[Epsilon], \[Nu], \[Lambda]][r]";
MSTRadialUp::usage = "MSTRadialUp[s, l, m, q, \[Epsilon], \[Nu], \[Lambda]][r]";*)

Begin["`Private`"];

$MasterFunction = "ReggeWheeler";

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
 Module[{a, b, c},
  Switch[$MasterFunction,
    "ReggeWheeler",
    a=\[Nu]+s+1-I \[Epsilon]; b=-\[Nu]+s-I \[Epsilon]; c=1-2 I \[Epsilon];,
    "Teukolsky",
    a=\[Nu]+1-I \[Tau]; b=-\[Nu]-I \[Tau]; c=1-s-I(\[Epsilon]+\[Tau]);,
    _, Abort[]
  ];
  Hypergeometric2F1[n + a, -n +b,c, x]
];

H2F1Up[n_, s_, \[Nu]_, \[Tau]_, \[Epsilon]_, x_] :=
 Module[{a, b, c},
  Switch[$MasterFunction,
    "ReggeWheeler",
    a=\[Nu]+s+1-I \[Epsilon]; b=-\[Nu]+s-I \[Epsilon]; c=1-2 I \[Epsilon];,
    "Teukolsky",
    a=\[Nu]+1-I \[Tau]; b=-\[Nu]-I \[Tau]; c=1-s-I(\[Epsilon]+\[Tau]);,
    _, Abort[]
  ];
  1/((3-a+b-2 n) (1+b-c-n) (-1+a+n)){-(1-a+b-2 n) (1+b-n) (-1+a-c+n) H2F1[-2+n], -(-2+a-b+2 n) (-2+2 a-2 b+2 a b+c-a c-b c+4 n-2 a n+2 b n-2 n^2+3 x-4 a x+a^2 x+4 b x-2 a b x+b^2 x-8 n x+4 a n x-4 b n x+4 n^2 x) H2F1[-1+n]}
]

H2F1Down[n_, s_, \[Nu]_, \[Tau]_, \[Epsilon]_, x_] :=
 Module[{a, b, c},
  Switch[$MasterFunction,
    "ReggeWheeler",
    a=\[Nu]+s+1-I \[Epsilon]; b=-\[Nu]+s-I \[Epsilon]; c=1-2 I \[Epsilon];,
    "Teukolsky",
    a=\[Nu]+1-I \[Tau]; b=-\[Nu]-I \[Tau]; c=1-s-I(\[Epsilon]+\[Tau]);,
    _, Abort[]
  ];
  1/((-3-a+b-2 n) (-1+b-n) (1+a-c+n)){(-2-a+b-2 n) (-2-2 a+2 b+2 a b+c-a c-b c-4 n-2 a n+2 b n-2 n^2+3 x+4 a x+a^2 x-4 b x-2 a b x+b^2 x+8 n x+4 a n x-4 b n x+4 n^2 x) H2F1[1+n], -(-1-a+b-2 n) (-1+b-c-n) (1+a+n) H2F1[2+n]}
];

dH2F1Exact[n_, s_, \[Nu]_, \[Tau]_, \[Epsilon]_, x_] :=
 Module[{a, b, c},
  Switch[$MasterFunction,
    "ReggeWheeler",
    a=\[Nu]+s+1-I \[Epsilon]; b=-\[Nu]+s-I \[Epsilon]; c=1-2 I \[Epsilon];,
    "Teukolsky",
    a=\[Nu]+1-I \[Tau]; b=-\[Nu]-I \[Tau]; c=1-s-I(\[Epsilon]+\[Tau]);,
    _, Abort[]
  ];
  (n+a)(-n+b)/c Hypergeometric2F1[n + a + 1, -n + b + 1, c + 1, x]
];

dH2F1Up[n_, s_, \[Nu]_, \[Tau]_, \[Epsilon]_, x_] :=
 Module[{a, b, c},
  Switch[$MasterFunction,
    "ReggeWheeler",
    a=\[Nu]+s+1-I \[Epsilon]; b=-\[Nu]+s-I \[Epsilon]; c=1-2 I \[Epsilon];,
    "Teukolsky",
    a=\[Nu]+1-I \[Tau]; b=-\[Nu]-I \[Tau]; c=1-s-I(\[Epsilon]+\[Tau]);,
    _, Abort[]
  ];
  1/((1+b-c-n) (-1+a+n)){-(((1-a+b-2 n) (1+b-n) (-1+a-c+n) dH2F1[-2+n])/(3-a+b-2 n) ), 1/(3-a+b-2 n)  (2-a+b-2 n) (-2+2 a-2 b+2 a b+c-a c-b c+4 n-2 a n+2 b n-2 n^2+3 x-4 a x+a^2 x+4 b x-2 a b x+b^2 x-8 n x+4 a n x-4 b n x+4 n^2 x) dH2F1[-1+n],(1-a+b-2 n) (2-a+b-2 n) H2F1[-1+n]}
];

dH2F1Down[n_, s_, \[Nu]_, \[Tau]_, \[Epsilon]_, x_] :=
 Module[{a, b, c},
  Switch[$MasterFunction,
    "ReggeWheeler",
    a=\[Nu]+s+1-I \[Epsilon]; b=-\[Nu]+s-I \[Epsilon];c=1-2 I \[Epsilon];,
    "Teukolsky",
    a=\[Nu]+1-I \[Tau]; b=-\[Nu]-I \[Tau]; c=1-s-I(\[Epsilon]+\[Tau]);,
    _, Abort[]
  ];
  1/((-1+b-n) (1+a-c+n)){1/(-3-a+b-2 n)  (-2-a+b-2 n) (-2-2 a+2 b+2 a b+c-a c-b c-4 n-2 a n+2 b n-2 n^2+3 x+4 a x+a^2 x-4 b x-2 a b x+b^2 x+8 n x+4 a n x-4 b n x+4 n^2 x) dH2F1[1+n], -(((-1-a+b-2 n) (-1+b-c-n) (1+a+n) dH2F1[2+n])/(-3-a+b-2 n)),(-2-a+b-2 n) (-1-a+b-2 n) H2F1[1+n]}
];

HUExact[n_, s_,\[Nu]_,\[Epsilon]_,zhat_] := Module[{a,b,c},Switch[$MasterFunction,"ReggeWheeler",a=\[Nu]+1-I \[Epsilon];,"Teukolsky",
a=\[Nu]+s+1-I \[Epsilon];,_,Abort[]];b=2\[Nu]+2;c=-2 I zhat;
  (c)^n HypergeometricU[n+a,2n+b,c]];

HUUp[n_,s_,\[Nu]_,\[Epsilon]_,zhat_]:=Module[{a,b,c},Switch[$MasterFunction,"ReggeWheeler",a=\[Nu]+1-I \[Epsilon];,"Teukolsky",
a=\[Nu]+s+1-I \[Epsilon];,_,Abort[]];b=2\[Nu]+2;c=-2 I zhat;1/((-1+a+n) (-4+b+2 n) ) {(-2-a+b+n) (-2+b+2 n) HU[-2+n],(-3+b+2 n) (8+(b+2 n)^2+2 (a+n) c-(b+2 n) (6+c)) HU[-1+n]/c}];

HUDown[n_,s_,\[Nu]_,\[Epsilon]_,zhat_]:=Module[{a,b,c},Switch[$MasterFunction,"ReggeWheeler",a=\[Nu]+1-I \[Epsilon];,"Teukolsky",
a=\[Nu]+s+1-I \[Epsilon];,_,Abort[]];b=2\[Nu]+2;c=-2 I zhat;1/((-a+b+n) (2+b+2 n)){-(((1+b+2 n) (b^2+4 n (1+n)+b (2+4 n-c)+2 a c) HU[1+n])/ c),(1+a+n) (b+2 n) HU[2+n]}];

dHUExact[n_, s_,\[Nu]_,\[Epsilon]_,zhat_] := Module[{a,b,c},Switch[$MasterFunction,"ReggeWheeler",a=\[Nu]+1-I \[Epsilon];,"Teukolsky",
a=\[Nu]+s+1-I \[Epsilon];,_,Abort[]];b=2\[Nu]+2;c=-2 I zhat;
(-2 I) (c^(-1+n) n HypergeometricU[a+n,b+2 n,c]-c^n (a+n) HypergeometricU[1+a+n,1+b+2 n,c])];

dHUUp[n_,s_,\[Nu]_,\[Epsilon]_,zhat_]:=Module[{a,b,c},Switch[$MasterFunction,"ReggeWheeler",a=\[Nu]+1-I \[Epsilon];,"Teukolsky",
a=\[Nu]+s+1-I \[Epsilon];,_,Abort[]];b=2\[Nu]+2;c=-2 I zhat;1/(-1+a+n) {((-2-a+b+n) (-2+b+2 n) dHU[-2+n])/(-4+b+2 n),((-3+b+2 n) (8+b^2+4 (-3+n) n+b (-6+4 n-c)+2 a c) dHU[-1+n])/( (-4+b+2 n) c),(2 I (-3+b+2 n) (-2+b+2 n) HU[-1+n])/c^2}];

dHUDown[n_,s_,\[Nu]_,\[Epsilon]_,zhat_]:=Module[{a,b,c},Switch[$MasterFunction,"ReggeWheeler",a=\[Nu]+1-I \[Epsilon];,"Teukolsky",
a=\[Nu]+s+1-I \[Epsilon];,_,Abort[]];b=2\[Nu]+2;c=-2 I zhat;1/((a-b-n) (2+b+2 n) c^2){(1+b+2 n) c (b^2+4 n (1+n)+b (2+4 n-c)+2 a c) dHU[1+n],(b+2 n) (-(1+a+n) c^2 dHU[2+n]),(b+2 n)(2 I (1+b+2 n) (2+b+2 n) HU[1+n])}];
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

fT[q_, \[Epsilon]_, \[Kappa]_, \[Tau]_, \[Nu]_, \[Lambda]_, s_, m_, 0] = 1;

fT[q_, \[Epsilon]_, \[Kappa]_, \[Tau]_, \[Nu]_, \[Lambda]_, s_, m_, nf_] :=
 fT[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, nf] =
 Module[{\[Alpha]n, \[Beta]n, \[Gamma]n, i, n, ret},
  \[Alpha]n[n_] := \[Alpha][q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, n];
  \[Beta]n[n_] := \[Beta][q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, n];
  \[Gamma]n[n_] := \[Gamma][q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, n];

  If[nf > 0,
    ret = fT[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, nf - 1] CF[-\[Alpha]n[i - 1] \[Gamma]n[i], \[Beta]n[i], {i, nf}]/\[Alpha]n[nf - 1];
  ,
    ret = fT[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, nf + 1] CF[-\[Alpha]n[2 nf - i] \[Gamma]n[2 nf - i + 1], \[Beta]n[2 nf - i], {i, nf}]/\[Gamma]n[nf + 1];
  ];
  
  ret
];

f[q_, \[Epsilon]_, \[Kappa]_, \[Tau]_, \[Nu]_, \[Lambda]_, s_, m_, nf_] :=
  Switch[$MasterFunction,
    "ReggeWheeler",
    Pochhammer[-\[Nu]+s-I \[Epsilon],-nf]Pochhammer[\[Nu]+1+s-I \[Epsilon],nf](-1)^nf Pochhammer[\[Nu]+I \[Epsilon]+1,nf]/Pochhammer[\[Nu]-I \[Epsilon]+1,nf],
    "Teukolsky",
    1
   ]fT[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, nf];
   
an[q_, \[Epsilon]_, \[Kappa]_, \[Tau]_, \[Nu]_, \[Lambda]_, s_, m_, nf_] :=(-1)^nf Pochhammer[\[Nu]+I \[Epsilon]+1,nf]/Pochhammer[\[Nu]-I \[Epsilon]+1,nf]fT[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, nf];
(******************************************************************************)
(*************************** Asymptotic amplitudes ****************************)
(******************************************************************************)
(*Does this section need a RW option?*)

Amplitudes[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_] :=
 Module[{\[Kappa], \[Tau], \[Epsilon]p, \[Omega], K\[Nu], Aplus, Btrans, Ctrans, Binc, Bref, n, fSumUp, fSumDown, fSumK\[Nu]1Up, fSumK\[Nu]1Down, fSumK\[Nu]2Up, fSumK\[Nu]2Down, fSumCUp, fSumCDown},
 Internal`InheritedBlock[{\[Alpha], \[Beta], \[Gamma], fT},
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
MSTRadialIn[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_][r_?InexactNumberQ] :=
 Module[{\[Kappa], \[Tau], rp, z, zp, x, resUp, nUp, resDown, nDown},
 Block[{H2F1},
 Internal`InheritedBlock[{\[Alpha], \[Beta], \[Gamma], fT},
  \[Kappa] = Sqrt[1 - q^2];
  \[Tau] = (\[Epsilon] - m q)/\[Kappa];
  rp = 1 + Sqrt[1 - q^2];
  z = \[Epsilon] r / 2;
  zp = \[Epsilon] rp / 2;
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

  Switch[$MasterFunction,
    "ReggeWheeler",
    (1-x)^(s+1) (-x)^(-I \[Epsilon]) E^(I \[Epsilon] x),
    "Teukolsky",
    (-x)^(-s - I (\[Epsilon] + \[Tau])/2) (1 - x)^(I (\[Epsilon] - \[Tau])/2) E^(I \[Epsilon] \[Kappa] x),
    _, Abort[]
  ](resUp + resDown) / norm
]]];


Derivative[1][MSTRadialIn[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_]][r_?InexactNumberQ] :=
 Module[{\[Kappa], \[Tau], rp, z, zp, x, dxdr, resUp, nUp, resDown, nDown},
 Block[{H2F1, dH2F1},
 Internal`InheritedBlock[{\[Alpha], \[Beta], \[Gamma], fT},
  \[Kappa] = Sqrt[1 - q^2];
  \[Tau] = (\[Epsilon] - m q)/\[Kappa];
  rp = 1 + Sqrt[1 - q^2];
  z = \[Epsilon] r / 2;
  zp = \[Epsilon] rp / 2;
  x = (zp - z)/(\[Epsilon] \[Kappa]);
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
While[resUp!=(resUp+=f[q,\[Epsilon],\[Kappa],\[Tau],\[Nu],\[Lambda],s,m,nUp] Switch[$MasterFunction,
"ReggeWheeler",-(((I (-I (1+s) x+(-1+x)^2 \[Epsilon])) H2F1[nUp])/x)+(1-x) dH2F1[nUp],
"Teukolsky",(-2 s (-1+x)+I (\[Epsilon]-2 x \[Epsilon] \[Kappa]+2 x^2 \[Epsilon] \[Kappa]+\[Tau]-2 x \[Tau])) H2F1[nUp]+2 (-1+x) x dH2F1[nUp],
_,Abort[]]),nUp++;];

  nDown = -1;
While[resDown!=(resDown+=f[q,\[Epsilon],\[Kappa],\[Tau],\[Nu],\[Lambda],s,m,nDown] Switch[$MasterFunction,
"ReggeWheeler",-(((I (-I (1+s) x+(-1+x)^2 \[Epsilon])) H2F1[nDown])/x)+(1-x) dH2F1[nDown],
"Teukolsky",(-2 s (-1+x)+I (\[Epsilon]-2 x \[Epsilon] \[Kappa]+2 x^2 \[Epsilon] \[Kappa]+\[Tau]-2 x \[Tau])) H2F1[nDown]+2 (-1+x) x dH2F1[nDown],
_,Abort[]]),nDown--;];

Switch[$MasterFunction,
"ReggeWheeler",E^(I x \[Epsilon]) (1-x)^s (-x)^(-I \[Epsilon]) ,
"Teukolsky",((E^(I x \[Epsilon] \[Kappa])) ((1-x)^(1/2 I (2 I+\[Epsilon]-\[Tau]))) ((-x)^(-1-s-1/2 I (\[Epsilon]+\[Tau]))) )/2,
_,Abort[]](resUp+resDown) dxdr/norm
]]];


Derivative[n_Integer?Positive][MSTRadialIn[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_]][r0_?InexactNumberQ] :=
 Module[{d2R, Rderivs, R, r, i},
 (*FIXME: Add appropriate factors of M here *)
  d2R = Switch[$MasterFunction,"ReggeWheeler",-1/(1-2/r)2/r^2 Derivative[1][R][r]+1/(1-2/r)(l (l+1)/r^2-2(1-s^2)/r^3)R[r]-(\[Epsilon]/2)^2/(1-2/r)^2 R[r],"Teukolsky",(-(-\[Lambda] + 2 I r s \[Epsilon] + (-2 I (-1 + r) s (-q m + (q^2 + r^2) \[Epsilon]/2) + (-q m + (q^2 + r^2) \[Epsilon]/2)^2)/(q^2 - 2 r + r^2)) R[r] - (-2 + 2 r) (1 + s) Derivative[1][R][r])/(q^2 - 2 r + r^2),_,Abort[]];

  pderivs = D[R[r_], {r_, i_}] :> D[d2R, {r, i - 2}] /; i >= 2;
  Do[Derivative[i][R][r] = Simplify[D[Derivative[i - 1][R][r], r] /. pderivs];, {i, 2, n}];
  Derivative[n][R][r] /. {
    R'[r] -> MSTRadialIn[s, l, m, q, \[Epsilon], \[Nu], \[Lambda], norm]'[r0],
    R[r] -> MSTRadialIn[s, l, m, q, \[Epsilon], \[Nu], \[Lambda], norm][r0], r -> r0}
];

(* Throwe B.5, Sasaki & Tagoshi (153) and (159) *)
SetAttributes[MSTRadialUp, {NumericFunction}];

MSTRadialUp[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_][r_?InexactNumberQ]   := 
 Module[{\[Kappa], \[Tau], \[Epsilon]p, rm, z, zm, zhat, resUp, nUp, resDown, nDown},
 Block[{HU},
 Internal`InheritedBlock[{\[Alpha], \[Beta], \[Gamma], fT},
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
    I^nUp Pochhammer[\[Nu] + 1 + s - I \[Epsilon], nUp]/Pochhammer[\[Nu] + 1 - s + I \[Epsilon], nUp](2zhat)^nUp /((*E^(I zhat)*) (-2 I zhat)^(nUp (*+ \[Nu] + 1*)) )Switch[$MasterFunction,"ReggeWheeler",Pochhammer[\[Nu]+1 -I \[Epsilon],nUp]/Pochhammer[\[Nu]+1+I \[Epsilon],nUp] an[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, nUp](-1)^nUp HU[nUp],"Teukolsky",fT[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, nUp]HU[nUp],_,Abort[]]),
    nUp++;
  ];
  
  nDown = -1;
  While[resDown != (resDown +=
    I^nDown Pochhammer[\[Nu] + 1 + s - I \[Epsilon], nDown]/Pochhammer[\[Nu] + 1 - s + I \[Epsilon], nDown](2zhat)^nDown/((*E^(I zhat)*) (-2 I zhat)^(nDown (*+ \[Nu] + 1*)) )Switch[$MasterFunction,"ReggeWheeler",Pochhammer[\[Nu]+1 -I \[Epsilon],nDown]/Pochhammer[\[Nu]+1+I \[Epsilon],nDown] an[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, nDown](-1)^nDown HU[nDown],"Teukolsky",fT[q, \[Epsilon], \[Kappa], \[Tau], \[Nu], \[Lambda], s, m, nDown]HU[nDown],_,Abort[]]),
    nDown--;
  ];

2^\[Nu] E^(-\[Pi] \[Epsilon]) E^(-I \[Pi] (\[Nu]+1)) E^(I zhat) zhat^(\[Nu]+I \[Epsilon]p) (zhat-\[Epsilon] \[Kappa])^(-I \[Epsilon]p) Switch[$MasterFunction,
"ReggeWheeler",zhat ,
"Teukolsky",E^(-I \[Pi] s) (zhat-\[Epsilon] \[Kappa])^(-s) ,
_,Abort[]](resUp+resDown)/norm
]]];


Derivative[1][MSTRadialUp[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_]][r_?InexactNumberQ] :=
 Module[{\[Kappa], \[Tau], \[Epsilon]p, rm, z, zm, zhat, dzhatdr, resUp, nUp, resDown, nDown},
 Block[{HU, dHU},
 Internal`InheritedBlock[{\[Alpha], \[Beta], \[Gamma], fT},
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

  resDown = resUp = 0;
  nUp = 0;

While[resUp!=(resUp+=Switch[$MasterFunction,
"ReggeWheeler",((-I)^nUp Pochhammer[\[Nu]+1+s-I \[Epsilon],nUp] Pochhammer[\[Nu]+1-I \[Epsilon],nUp] an[q,\[Epsilon],\[Kappa],\[Tau],\[Nu],\[Lambda],s,m,nUp] (2 zhat)^nUp ((1+I \[Epsilon]p+zhat (I+(I \[Epsilon]p)/(-zhat+\[Epsilon] \[Kappa]))+\[Nu]) HU[nUp]+zhat dHU[nUp]))/((Pochhammer[\[Nu]+1-s+I \[Epsilon],nUp] Pochhammer[\[Nu]+1+I \[Epsilon],nUp]) (-2 I zhat)^nUp),
"Teukolsky",(I^nUp Pochhammer[\[Nu]+1+s-I \[Epsilon],nUp] (2 zhat)^nUp f[q,\[Epsilon],\[Kappa],\[Tau],\[Nu],\[Lambda],s,m,nUp] (((-s zhat+I (zhat^2-\[Epsilon] (zhat+\[Epsilon]p) \[Kappa])+(zhat-\[Epsilon] \[Kappa]) \[Nu]) HU[nUp])/(zhat (zhat-\[Epsilon] \[Kappa]))+dHU[nUp]))/(Pochhammer[\[Nu]+1-s+I \[Epsilon],nUp] (-2 I zhat)^nUp),
_,Abort[]]), 
nUp++;
  ];
  
  nDown = -1;
 While[resDown!=(resDown+=Switch[$MasterFunction,
"ReggeWheeler",((-I)^nUp Pochhammer[\[Nu]+1+s-I \[Epsilon],nUp] Pochhammer[\[Nu]+1-I \[Epsilon],nUp] an[q,\[Epsilon],\[Kappa],\[Tau],\[Nu],\[Lambda],s,m,nUp] (2 zhat)^nUp ((1+I \[Epsilon]p+zhat (I+(I \[Epsilon]p)/(-zhat+\[Epsilon] \[Kappa]))+\[Nu]) HU[nUp]+zhat dHU[nUp]))/((Pochhammer[\[Nu]+1-s+I \[Epsilon],nUp] Pochhammer[\[Nu]+1+I \[Epsilon],nUp]) (-2 I zhat)^nUp),
"ReggeWheeler",((-I)^nDown Pochhammer[\[Nu]+1+s-I \[Epsilon],nDown] Pochhammer[\[Nu]+1-I \[Epsilon],nDown] an[q,\[Epsilon],\[Kappa],\[Tau],\[Nu],\[Lambda],s,m,nDown] (2 zhat)^nDown ((1+I \[Epsilon]p+zhat (I+(I \[Epsilon]p)/(-zhat+\[Epsilon] \[Kappa]))+\[Nu]) HU[nDown]+zhat dHU[nDown]))/((Pochhammer[\[Nu]+1-s+I \[Epsilon],nDown] Pochhammer[\[Nu]+1+I \[Epsilon],nDown]) (-2 I zhat)^nDown),
"Teukolsky",(I^nDown Pochhammer[\[Nu]+1+s-I \[Epsilon],nDown] (2 zhat)^nDown f[q,\[Epsilon],\[Kappa],\[Tau],\[Nu],\[Lambda],s,m,nDown] (((-s zhat+I (zhat^2-\[Epsilon] (zhat+\[Epsilon]p) \[Kappa])+(zhat-\[Epsilon] \[Kappa]) \[Nu]) HU[nDown])/(zhat (zhat-\[Epsilon] \[Kappa]))+dHU[nDown]))/(Pochhammer[\[Nu]+1-s+I \[Epsilon],nDown] (-2 I zhat)^nDown),
_,Abort[]]),
 nDown--;
  ];

2^\[Nu] E^(-\[Pi] \[Epsilon]) E^(-I \[Pi] (\[Nu]+1)) E^(I zhat) zhat^(\[Nu]+I \[Epsilon]p) (zhat-\[Epsilon] \[Kappa])^(-I \[Epsilon]p) Switch[$MasterFunction,
"ReggeWheeler",1,
"Teukolsky",E^(-I \[Pi] s) (zhat-\[Epsilon] \[Kappa])^(-s),
_,Abort[]] (resUp+resDown) dzhatdr/norm
]]];



Derivative[n_Integer?Positive][MSTRadialUp[s_Integer, l_Integer, m_Integer, q_, \[Epsilon]_, \[Nu]_, \[Lambda]_, norm_]][r0_?InexactNumberQ] :=
 Module[{d2R, Rderivs, R, r, i},
  d2R = Switch[$MasterFunction,"ReggeWheeler",-1/(1-2/r)2/r^2 Derivative[1][R][r]+1/(1-2/r)(l (l+1)/r^2-2(1-s^2)/r^3)R[r]-(\[Epsilon]/2)^2/(1-2/r)^2 R[r],"Teukolsky",(-(-\[Lambda] + 2 I r s \[Epsilon] + (-2 I (-1 + r) s (-q m + (q^2 + r^2) \[Epsilon]/2) + (-q m + (q^2 + r^2) \[Epsilon]/2)^2)/(q^2 - 2 r + r^2)) R[r] - (-2 + 2 r) (1 + s) Derivative[1][R][r])/(q^2 - 2 r + r^2),_,Abort[]];

  pderivs = D[R[r_], {r_, i_}] :> D[d2R, {r, i - 2}] /; i >= 2;
  Do[Derivative[i][R][r] = Simplify[D[Derivative[i - 1][R][r], r] /. pderivs];, {i, 2, n}];
  Derivative[n][R][r] /. {
    R'[r] -> MSTRadialUp[s, l, m, q, \[Epsilon], \[Nu], \[Lambda], norm]'[r0],
    R[r] -> MSTRadialUp[s, l, m, q, \[Epsilon], \[Nu], \[Lambda], norm][r0], r -> r0}
];


End[];
EndPackage[];
