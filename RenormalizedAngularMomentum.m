(* ::Package:: *)

BeginPackage["Teukolsky`RenormalizedAngularMomentum`",
  {"SpinWeightedSpheroidalHarmonics`"}
];

RenormalizedAngularMomentum::usage = "RenormalizedAngularMomentum[M, a, \[Omega], s, l, m] gives the renormalized angular momentum \[Nu].";

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
Cos2\[Pi]\[Nu]Series[M_, a_, \[Omega]_, s_, l_, m_] :=
 Module[{data, \[Omega]Exp, qExp, q, Cos2\[Pi]\[Nu]},
  data = coeffs[s, l, m];
  If[data == $Failed, Return[$Failed]];
  
  (* Exponents of \[Omega] and q *)
  \[Omega]Exp = Range[0, Length[data] - 1];
  qExp = Range[0, Length[First[data]] - 1];
  
  q = a / M;
  
  Cos2\[Pi]\[Nu] = 1 + \[Omega]^\[Omega]Exp.(data.q^qExp);
  Cos2\[Pi]\[Nu]
]

(* Find \[Nu] using root finding with an initial guess *)
\[Nu]RootFind[M_, a_, \[Omega]_, s_, l_, m_, Cos2\[Pi]\[Nu]_] :=
  \[Nu]RootFind[M, a, \[Omega], s, l, m, SpinWeightedSpheroidalEigenvalue[s, l, m, a \[Omega]], Cos2\[Pi]\[Nu]];

\[Nu]RootFind[M_, a_, \[Omega]_, s_, l_, m_, \[Lambda]_, Cos2\[Pi]\[Nu]_] :=
 Module[{q, \[Epsilon], \[Kappa], \[Tau], \[Alpha], \[Beta], \[Gamma], R, L, \[Nu], \[Nu]0, \[Nu]1, \[Nu]2, \[Nu]3, \[Nu]i, i, precision, Mp, ap, \[Omega]p, \[Lambda]p, Cos2\[Pi]\[Nu]p, nmax=25(*FIXME*), res},
  precision = Precision[{M, a, \[Omega], \[Lambda]}];
  {Mp, ap, \[Omega]p, \[Lambda]p, Cos2\[Pi]\[Nu]p} = SetPrecision[{M, a, \[Omega], \[Lambda], Cos2\[Pi]\[Nu]}, precision];

  q = ap/Mp;
  \[Epsilon] = 2 Mp \[Omega]p;
  \[Kappa]=Sqrt[1-q^2];
  \[Tau]=(\[Epsilon]-m q)/\[Kappa];

  \[Alpha][n_, \[Nu]_?InexactNumberQ] := (I \[Epsilon] \[Kappa](n+\[Nu]+1+s+I \[Epsilon])(n+\[Nu]+1+s-I \[Epsilon])(n+\[Nu]+1+I \[Tau]))/((n+\[Nu]+1)(2n+2\[Nu]+3));
  \[Beta][n_, \[Nu]_?InexactNumberQ] := -\[Lambda]p-s(s+1)+(n+\[Nu])(n+\[Nu]+1)+\[Epsilon]^2+\[Epsilon](\[Epsilon]-m q)+(\[Epsilon](\[Epsilon]-m q)(s^2+\[Epsilon]^2))/((n+\[Nu])(n+\[Nu]+1));
  \[Gamma][n_, \[Nu]_?InexactNumberQ] := -((I \[Epsilon] \[Kappa](n+\[Nu]-s+I \[Epsilon])(n+\[Nu]-s-I \[Epsilon])(n+\[Nu]-I \[Tau]))/((n+\[Nu])(2n+2\[Nu]-1)));

  R[n_, \[Nu]_] := ContinuedFractionK[-\[Alpha][i-1, \[Nu]] \[Gamma][i, \[Nu]],\[Beta][i, \[Nu]], {i, n, n+nmax}]/\[Alpha][n-1, \[Nu]];
  L[n_, \[Nu]_] := ContinuedFractionK[-\[Alpha][2n-i, \[Nu]] \[Gamma][2n-i+1, \[Nu]], \[Beta][2n-i, \[Nu]], {i, n, n+nmax}]/\[Gamma][n+1, \[Nu]];

  \[Nu]0 = l - ArcCos[Cos2\[Pi]\[Nu]p] / (2*\[Pi]);

  (* There are three possible cases: *)
  Check[Which[
    -1 <= Cos2\[Pi]\[Nu] <= 1,
      \[Nu] /. FindRoot[Re[\[Beta][0, \[Nu]]+\[Alpha][0, \[Nu]] R[1, \[Nu]]+\[Gamma][0, \[Nu]] L[-1, \[Nu]]] == 0, {\[Nu], \[Nu]0, 9/10 \[Nu]0, 11/10 \[Nu]0}, WorkingPrecision -> precision]
    ,
    Cos2\[Pi]\[Nu] < 1,
      \[Nu]0 = Im[\[Nu]0];
      Check[
        1/2 + I \[Nu]i /. FindRoot[Re[(\[Beta][0, \[Nu]]+\[Alpha][0, \[Nu]] R[1, \[Nu]]+\[Gamma][0, \[Nu]] L[-1, \[Nu]] /. \[Nu] -> 1/2 + I \[Nu]i)] == 0, {\[Nu]i, 11/10 \[Nu]0, 9/10 \[Nu]0}, Method->"Brent", WorkingPrecision -> precision],
        1/2 + I \[Nu]i /. FindRoot[Re[(\[Beta][0, \[Nu]]+\[Alpha][0, \[Nu]] R[1, \[Nu]]+\[Gamma][0, \[Nu]] L[-1, \[Nu]] /. \[Nu] -> 1/2 + I \[Nu]i)] == 0, {\[Nu]i, \[Nu]0, 11/10 \[Nu]0, 9/10 \[Nu]0}, WorkingPrecision -> precision],
        FindRoot::bbrac]
    ,
    Cos2\[Pi]\[Nu] > 1,
      I \[Nu]i /. FindRoot[Re[\[Beta][0, \[Nu]]+\[Alpha][0, \[Nu]] R[1, \[Nu]]+\[Gamma][0, \[Nu]] L[-1, \[Nu]] /. \[Nu] -> I \[Nu]i] == 0, {\[Nu]i, \[Nu]0, 9/10 \[Nu]0, 11/10 \[Nu]0}, WorkingPrecision -> precision]
    ,
    True,
      $Failed
    ],
    $Failed
  ]
];

(**********************************************************)
(* RenormalizedAngularMomentum                            *)
(**********************************************************)

SyntaxInformation[RenormalizedAngularMomentum] =
 {"ArgumentsPattern" -> {_, _, _, _, _, _, ___}};
Options[RenormalizedAngularMomentum] = {Method -> Automatic};
SetAttributes[RenormalizedAngularMomentum, {NumericFunction}];

RenormalizedAngularMomentum[M_, a_, \[Omega]_, s_, l_, m_] /; l < Abs[s] := 0; (* FIXME: unphysical cases *)

RenormalizedAngularMomentum[M_, (0|0.), (0|0.), 0, 0, 0] := 0; (* FIXME: special cases *)

RenormalizedAngularMomentum[M_?NumericQ, a_?NumericQ, \[Omega]_?NumericQ, s_Integer, l_Integer, m_Integer] /; InexactNumberQ[a] || InexactNumberQ[\[Omega]] :=
 Module[{\[Nu], Cos2\[Pi]\[Nu]},
  Cos2\[Pi]\[Nu] = Cos2\[Pi]\[Nu]Series[N[M], N[a], N[\[Omega]], s, l, m]; (* FIXME: catch failure *)
  \[Nu] = \[Nu]RootFind[M, a, \[Omega], s, l, m, Cos2\[Pi]\[Nu]];

  \[Nu]
];

RenormalizedAngularMomentum /: N[RenormalizedAngularMomentum[M_?NumericQ, a_?NumericQ, \[Omega]_?NumericQ, s_Integer, l_Integer, m_Integer], Nopts___] :=
  RenormalizedAngularMomentum[M, N[a, Nopts], N[\[Omega], Nopts], s, l, m];

End[];
EndPackage[];

(* Add h5mma to $ContextPath since Get[] inside `Private` does not do so. *)
If[SimulationTools`ReadHDF5`Private`$h5mma,
  If[!MemberQ[$ContextPath, "h5mma`"], AppendTo[$ContextPath, "h5mma`"]];
];
