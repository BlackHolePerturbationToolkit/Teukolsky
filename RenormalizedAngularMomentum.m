(* ::Package:: *)

BeginPackage["Teukolsky`RenormalizedAngularMomentum`",
  {"SpinWeightedSpheroidalHarmonics`"}
];

RenormalizedAngularMomentum::usage =
 "RenormalizedAngularMomentum[s, l, m, q, \[Epsilon], \[Lambda]] gives the renormalized angular momentum \[Nu].\n" <>
 "RenormalizedAngularMomentum[s, l, m, q, \[Epsilon]] gives the renormalized angular momentum \[Nu].";

Begin["`Private`"];

(**********************************************************)
(* Internal functions                                     *)
(**********************************************************)

(* Read in data for series coefficients from hdf5 files *)

(* If the h5mma is not found, then just use Mathematica's built-in HDF5 support *)
$h5mma = If[Quiet[Get["h5mma`"], {Needs::nocont, Get::noopen}]===$Failed, False, True];
If[$h5mma, SetOptions[ImportHDF5, Turbo->True]];

ReadHDF5[file_String, opts_:"Datasets"] :=
  If[$h5mma, ImportHDF5[file, opts], Import[file, opts]];

coeffsFile = FileNameJoin[{FileNameDrop[FindFile["Teukolsky`"], -2], "nu_coeffs.h5"}];

modesAvailable := modesAvailable =
  Sort[Flatten[StringCases[ReadHDF5[coeffsFile, "Datasets"],
    "s" ~~ s : NumberString ~~ "l" ~~ l : NumberString ~~ "m" ~~ m : NumberString :> {ToExpression[s], ToExpression[l], ToExpression[m]}], 1]];

coeffs[s_Integer, l_Integer, m_Integer] := coeffs[s, l, m] =
 Module[{},
  If[!MemberQ[modesAvailable, {s, l, m}], Return[$Failed]];
  ReadHDF5[coeffsFile, {"Datasets", "s"<>ToString[s]<>"l"<>ToString[l]<>"m"<>ToString[m]}]
]

(* Compute an approximation for nu from its series expansion *)
Cos2\[Pi]\[Nu]Series[a_, \[Omega]_, s_, l_, m_] :=
 Module[{data, \[Omega]Exp, qExp, q, Cos2\[Pi]\[Nu]},
  data = coeffs[s, l, m];
  If[data == $Failed, Return[$Failed]];
  
  (* Exponents of \[Omega] and q *)
  \[Omega]Exp = Range[0, Length[data] - 1];
  qExp = Range[0, Length[First[data]] - 1];
  
  q = a;
  
  Cos2\[Pi]\[Nu] = 1 + \[Omega]^\[Omega]Exp.(data.q^qExp);
  Cos2\[Pi]\[Nu]
]

(* Find \[Nu] using root finding with an initial guess *)
\[Nu]RootFind[a_, \[Omega]_, \[Lambda]_, s_, l_, m_, Cos2\[Pi]\[Nu]_?NumericQ] :=
 Module[{q, \[Epsilon], \[Kappa], \[Tau], \[Alpha]\[Gamma], \[Beta], R, L, \[Nu], \[Nu]0, \[Nu]1, \[Nu]2, \[Nu]3, \[Nu]i, i, precision, ap, \[Omega]p, \[Lambda]p, Cos2\[Pi]\[Nu]p, nmax=25(*FIXME*), res},
  precision = Precision[{a, \[Omega], \[Lambda]}];
  {ap, \[Omega]p, \[Lambda]p, Cos2\[Pi]\[Nu]p} = SetPrecision[{a, \[Omega], \[Lambda], Cos2\[Pi]\[Nu]}, precision];

  q = ap;
  \[Epsilon] = 2 \[Omega]p;
  \[Kappa]=Sqrt[1-q^2];
  \[Tau]=(\[Epsilon]-m q)/\[Kappa];

  \[Alpha]\[Gamma][n_, \[Nu]_?InexactNumberQ] := \[Epsilon]^2 \[Kappa]^2 (n + \[Nu]) (2 + n + \[Nu]) ((1 + n + \[Nu] - s)^2 + \[Epsilon]^2) ((1 + n + \[Nu] + s)^2 + \[Epsilon]^2) (-1 + 2 n + 2 \[Nu]) (5 + 2 n + 2 \[Nu]) ((1 + n + \[Nu])^2 + \[Tau]^2);
  \[Beta][n_, \[Nu]_?InexactNumberQ] := (2 n + 2 \[Nu] + 3) (2 n + 2 \[Nu] - 1) ((-\[Lambda] - s (s + 1) + (n + \[Nu]) (n + \[Nu] + 1) + \[Epsilon]^2 + \[Epsilon] (\[Epsilon] - m q)) ((n + \[Nu]) (n + \[Nu] + 1)) + (\[Epsilon] (\[Epsilon] - m q) (s^2 + \[Epsilon]^2)));

  R[n_, \[Nu]_] := ContinuedFractionK[-\[Alpha]\[Gamma][i-1, \[Nu]], \[Beta][i, \[Nu]], {i, n, n+nmax}];
  L[n_, \[Nu]_] := ContinuedFractionK[-\[Alpha]\[Gamma][2n-i, \[Nu]], \[Beta][2n-i, \[Nu]], {i, n, n+nmax}];

  \[Nu]0 = l - ArcCos[Cos2\[Pi]\[Nu]p] / (2*\[Pi]);

  (* There are three possible cases: *)
  Check[Which[
    -1 <= Cos2\[Pi]\[Nu] <= 1, (* \[Nu] Real *)
      \[Nu] /. FindRoot[Re[\[Beta][0, \[Nu]] + R[1, \[Nu]] + L[-1, \[Nu]]] == 0, {\[Nu], \[Nu]0, 9/10 \[Nu]0, 11/10 \[Nu]0}, WorkingPrecision -> precision]
    ,
    Cos2\[Pi]\[Nu] < -1, (* \[Nu] = 1/2 + I \[Nu]i *)
      1/2 + I \[Nu]i /. FindRoot[Re[(\[Beta][0, \[Nu]] + R[1, \[Nu]]+ L[-1, \[Nu]] /. \[Nu] -> 1/2 + I \[Nu]i)] == 0, {\[Nu]i, Im[\[Nu]0], 9/10 Im[\[Nu]0], 11/10 Im[\[Nu]0]}, WorkingPrecision -> precision]
    ,
    Cos2\[Pi]\[Nu] > 1, (* \[Nu] = I \[Nu]i *)
      I \[Nu]i /. FindRoot[Re[\[Beta][0, \[Nu]] + R[1, \[Nu]] + L[-1, \[Nu]] /. \[Nu] -> I \[Nu]i] == 0, {\[Nu]i, Im[\[Nu]0], 11/10 Im[\[Nu]0], 9/10 Im[\[Nu]0]}, WorkingPrecision -> precision]
    ,
    True,
      $Failed
    ],
    $Failed
  ]
];

\[Nu]RCHMonodromy[a_, \[Omega]_, \[Lambda]_, s_, l_, m_, Npmax_] :=
 Module[{q, \[Epsilon], \[Kappa], \[Tau], \[Gamma]CH, \[Delta]CH, \[Epsilon]CH, \[Alpha]CH, qCH, \[Mu]1C, \[Mu]2C, a1, a2, a1sum, a2sum, Pochhammerp1m2, Pochhammerm1p2, Cos2\[Pi]\[Nu], nmax, Cos2\[Pi]\[Nu]precision, \[Nu]},
  q = a;
  \[Epsilon] = 2 \[Omega];
  \[Kappa] = Sqrt[1-q^2];
  \[Tau] = (\[Epsilon]-m q)/\[Kappa];

  (* Compute parameters in confluent Heun equation *)
  \[Gamma]CH = 1-s-I \[Epsilon]-I \[Tau];
  \[Delta]CH = 1+s+I \[Epsilon]-I \[Tau];
  \[Epsilon]CH = 2I \[Epsilon] \[Kappa];
  \[Alpha]CH = 2 \[Epsilon] \[Kappa] (I-I s-\[Epsilon]+\[Tau]);
  qCH = -(-s (1+s)+\[Epsilon]^2+I (-1+2 s) \[Epsilon] \[Kappa]-\[Lambda]-\[Tau] (I+\[Tau]));

  \[Mu]1C = \[Alpha]CH/\[Epsilon]CH-(\[Gamma]CH+\[Delta]CH);
  \[Mu]2C=-(\[Alpha]CH/\[Epsilon]CH);

  (* Recurrence relations *)
  a1[-1] = 0;
  a1[0] = 1;
  a1[n_] := a1[n] = ((\[Alpha]CH-(-1+n+\[Delta]CH) \[Epsilon]CH) (\[Alpha]CH-(-2+n+\[Gamma]CH+\[Delta]CH) \[Epsilon]CH) )/(n \[Epsilon]CH) a1[n-2]-1/(n \[Epsilon]CH^2) (\[Alpha]CH^2+\[Alpha]CH \[Epsilon]CH (1-2 n-\[Gamma]CH-\[Delta]CH+\[Epsilon]CH)+\[Epsilon]CH^2 (n^2-qCH+n (-1+\[Gamma]CH+\[Delta]CH-\[Epsilon]CH)+\[Epsilon]CH-\[Delta]CH \[Epsilon]CH)) a1[n-1];

  a2[-1] = 0;
  a2[0] = 1;
  a2[n_] := a2[n] = -(((\[Alpha]CH+(-2+n) \[Epsilon]CH) (\[Alpha]CH+(-1+n-\[Gamma]CH) \[Epsilon]CH) )/(n \[Epsilon]CH))a2[n-2]+1/(n \[Epsilon]CH^2) (\[Alpha]CH^2+(n^2-qCH+\[Gamma]CH+\[Delta]CH-n (1+\[Gamma]CH+\[Delta]CH-\[Epsilon]CH)-\[Epsilon]CH) \[Epsilon]CH^2+\[Alpha]CH \[Epsilon]CH (-1+2 n-\[Gamma]CH-\[Delta]CH+\[Epsilon]CH)) a2[n-1];

  Pochhammerp1m2[0] = 1;
  Pochhammerp1m2[i_] := Pochhammerp1m2[i] = (-\[Mu]2C+\[Mu]1C+i-1)Pochhammerp1m2[i-1];
  Pochhammerm1p2[0] = 1;
  Pochhammerm1p2[i_] := Pochhammerm1p2[i] = (\[Mu]2C-\[Mu]1C+i-1)Pochhammerm1p2[i-1];

  a1sum[n_] := Gamma[-\[Mu]2C+\[Mu]1C] Sum[a1[j]Pochhammerp1m2[n-j], {j, 0, Ceiling[n/2]}]; 
  a2sum[n_] := Gamma[\[Mu]2C-\[Mu]1C] Sum[(-1)^j a2[j]Pochhammerm1p2[n-j], {j, 0, Ceiling[n/2]}];

  (* Compute \[Nu], with error estimate (precision of output) based on the assumption that Cos[2\[Pi]\[Nu]] should be real. *)
  Cos2\[Pi]\[Nu][nmax_] := Cos2\[Pi]\[Nu][nmax] = Cos[\[Pi](\[Mu]1C-\[Mu]2C)]+(2\[Pi]^2)/(a1sum[nmax] a2sum[nmax]) (-1)^(nmax-1) a1[nmax]a2[nmax];
  If[IntegerQ[Npmax],
    nmax = Npmax;
    If[Precision[Cos2\[Pi]\[Nu][nmax]] == 0, Return[$Failed]];
    Cos2\[Pi]\[Nu]precision = -RealExponent[Im[Cos2\[Pi]\[Nu][nmax]/Re[Cos2\[Pi]\[Nu][nmax]]]];
  ,
    (* FIXME: we should be able tor predict nmax based on the convergence for large nmax and the loss of precision in a1 and a2 *)
    nmax = 2 Ceiling[E^ProductLog[Precision[{a, \[Omega], \[Lambda]}] Log[100]]];
    If[Precision[Cos2\[Pi]\[Nu][nmax]] == 0, Return[$Failed]];
    Cos2\[Pi]\[Nu]precision = 0;
    (* Increase nmax by 10% until the precision of the result decreases *)
    While[Cos2\[Pi]\[Nu]precision < (Cos2\[Pi]\[Nu]precision = -RealExponent[Im[Cos2\[Pi]\[Nu][nmax]/Re[Cos2\[Pi]\[Nu][nmax]]]]),
      nmax = Round[11/10 nmax];
      If[Precision[Cos2\[Pi]\[Nu][nmax]] == 0, Return[$Failed]];
    ];
    nmax = Round[10/11 nmax];
  ];
    
  If[Precision[Cos2\[Pi]\[Nu][nmax]]=!=MachinePrecision,
    Cos2\[Pi]\[Nu][nmax] = N[Cos2\[Pi]\[Nu][nmax], Max[-RealExponent[Im[Cos2\[Pi]\[Nu][nmax]]/Re[Cos2\[Pi]\[Nu][nmax]]],0]];
  ];

  \[Nu] = Which[
    Re[Cos2\[Pi]\[Nu][nmax]]<-1, 
      1/2-Im[ArcCos[Re[Cos2\[Pi]\[Nu][nmax]]]/(2\[Pi])]I,
    -1<=Re[Cos2\[Pi]\[Nu][nmax]]<1,
      l-ArcCos[Re[Cos2\[Pi]\[Nu][nmax]]]/(2\[Pi]),
    Re[Cos2\[Pi]\[Nu][nmax]]>1,
      -I Im[ArcCos[Re[Cos2\[Pi]\[Nu][nmax]]]/(2\[Pi])],
    True,
      $Failed
  ];
  
  \[Nu]
];

(**********************************************************)
(* RenormalizedAngularMomentum                            *)
(**********************************************************)

SyntaxInformation[RenormalizedAngularMomentum] =
 {"ArgumentsPattern" -> {_, _, _, _, _, ___}};
Options[RenormalizedAngularMomentum] = {Method -> Automatic};
SetAttributes[RenormalizedAngularMomentum, {NumericFunction}];

RenormalizedAngularMomentum[s_, l_, m_, q_, \[Epsilon]_, \[Lambda]_, OptionsPattern[RenormalizedAngularMomentum]] /; l < Abs[s] := 0;

RenormalizedAngularMomentum[s_, l_, m_, q_, (0|0.), \[Lambda]_, OptionsPattern[RenormalizedAngularMomentum]] := l(l-1);

RenormalizedAngularMomentum[s_Integer, l_Integer, m_Integer, q_?NumericQ, \[Epsilon]_?NumericQ, \[Lambda]_?NumericQ,
 Method -> {"FindRoot", "InitialGuess" -> \[Nu]_}] /; InexactNumberQ[q] || InexactNumberQ[\[Epsilon]] || InexactNumberQ[\[Lambda]] :=
  \[Nu]RootFind[q, \[Epsilon]/2, \[Lambda], s, l, m, Cos[2 \[Pi] \[Nu]]];

RenormalizedAngularMomentum[s_Integer, l_Integer, m_Integer, q_?NumericQ, \[Epsilon]_?NumericQ, \[Lambda]_?NumericQ, Method -> ("FindRoot"|{"FindRoot"})] /; InexactNumberQ[q] || InexactNumberQ[\[Epsilon]] || InexactNumberQ[\[Lambda]] :=
 Module[{\[Nu], Cos2\[Pi]\[Nu]},
  Cos2\[Pi]\[Nu] = Cos2\[Pi]\[Nu]Series[N[q], N[\[Epsilon]/2], s, l, m];
  If[Cos2\[Pi]\[Nu] == $Failed, Return[$Failed]];
  \[Nu] = \[Nu]RootFind[q, \[Epsilon]/2, \[Lambda], s, l, m, Cos2\[Pi]\[Nu]];

  \[Nu]
];

RenormalizedAngularMomentum[s_Integer, l_Integer, m_Integer, q_?NumericQ, \[Epsilon]_?NumericQ, \[Lambda]_?NumericQ,
 Method -> ("Monodromy"|{"Monodromy"})] /; InexactNumberQ[q] || InexactNumberQ[\[Epsilon]] || InexactNumberQ[\[Lambda]] :=
  RenormalizedAngularMomentum[s, l, m, q, \[Epsilon], \[Lambda], Method -> {"Monodromy", "nmax" -> Automatic}];

RenormalizedAngularMomentum[s_Integer, l_Integer, m_Integer, q_?NumericQ, \[Epsilon]_?NumericQ, \[Lambda]_?NumericQ, 
 Method -> {"Monodromy", "nmax" -> nmax_}] /; InexactNumberQ[q] || InexactNumberQ[\[Epsilon]] || InexactNumberQ[\[Lambda]] :=
  \[Nu]RCHMonodromy[q, \[Epsilon]/2, \[Lambda], s, l, m, nmax];

RenormalizedAngularMomentum[s_Integer, l_Integer, m_Integer, q_?NumericQ, \[Epsilon]_?NumericQ, \[Lambda]_?NumericQ, 
 Method -> ("Series"|{"Series"})] /; InexactNumberQ[q] || InexactNumberQ[\[Epsilon]] || InexactNumberQ[\[Lambda]] :=
  l - ArcCos[Cos2\[Pi]\[Nu]Series[q, \[Epsilon]/2, s, l, m]] / (2*\[Pi]);

RenormalizedAngularMomentum[s_Integer, l_Integer, m_Integer, q_?NumericQ, \[Epsilon]_?NumericQ,
 Method -> ("Series"|{"Series"})] /; InexactNumberQ[q] || InexactNumberQ[\[Epsilon]] || InexactNumberQ[\[Lambda]] :=
  l - ArcCos[Cos2\[Pi]\[Nu]Series[q, \[Epsilon]/2, s, l, m]] / (2*\[Pi]);

RenormalizedAngularMomentum[s_Integer, l_Integer, m_Integer, q_?NumericQ, \[Epsilon]_?NumericQ, \[Lambda]_?NumericQ] /;
 InexactNumberQ[q] || InexactNumberQ[\[Epsilon]] || InexactNumberQ[\[Lambda]] :=
   RenormalizedAngularMomentum[s, l, m, q, \[Epsilon], \[Lambda], Method -> "FindRoot"];

RenormalizedAngularMomentum[s_Integer, l_Integer, m_Integer, q_?NumericQ, \[Epsilon]_?NumericQ, opts:OptionsPattern[RenormalizedAngularMomentum]] :=
  RenormalizedAngularMomentum[s, l, m, q, \[Epsilon], SpinWeightedSpheroidalEigenvalue[s, l, m, q \[Epsilon]/2], opts];

RenormalizedAngularMomentum /: N[RenormalizedAngularMomentum[s_Integer, l_Integer, m_Integer, q_?NumericQ, \[Epsilon]_?NumericQ, \[Lambda]_?NumericQ], Nopts:OptionsPattern[N]] :=
  RenormalizedAngularMomentum[s, l, m, N[q, Nopts], N[\[Epsilon], Nopts], N[\[Lambda], Nopts]];

End[];
EndPackage[];

(* Add h5mma to $ContextPath since Get[] inside `Private` does not do so. *)
If[SimulationTools`ReadHDF5`Private`$h5mma,
  If[!MemberQ[$ContextPath, "h5mma`"], AppendTo[$ContextPath, "h5mma`"]];
];
