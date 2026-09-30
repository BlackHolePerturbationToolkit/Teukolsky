(* ::Package:: *)

BeginPackage[MST`$MasterFunction<>"`MST`RenormalizedAngularMomentum`",
  {MST`$MasterFunction<>"`", "SpinWeightedSpheroidalHarmonics`"}
];

ClearAttributes[RenormalizedAngularMomentum, {Protected, ReadProtected}];

RenormalizedAngularMomentum::usage =
 "RenormalizedAngularMomentum[s, l, m, a, \[Omega], \[Lambda]] gives the renormalized angular momentum \[Nu].\n" <>
 "RenormalizedAngularMomentum[s, l, m, a, \[Omega]] gives the renormalized angular momentum \[Nu].\n" <>
 "Of the equivalent values \[PlusMinus]\[Nu] + k, the one returned is continuous with \[Nu] = l at \[Omega] = 0 (l - ArcCos[Cos[2\[Pi]\[Nu]]]/(2\[Pi])), for real and complex \[Omega].";

(* Messages *)
RenormalizedAngularMomentum::precision = "Method \"Monodromy\" currently only works reliably with arbitrary precision input parameters.";
RenormalizedAngularMomentum::conv = "The monodromy method did not converge for \[Omega] = `1`.";

Begin["`Private`"];

(**********************************************************)
(* Internal functions                                     *)
(**********************************************************)

(* Compute an approximation for nu from its series expansion - cos of Eq. (173) in Sasaki-Tagoshi *)
Cos2\[Pi]\[Nu]Series[a_, \[Omega]_, s_, l_, m_] :=
 Module[{Cos2\[Pi]\[Nu]},
  If[l==0,
    Cos2\[Pi]\[Nu] = 1 - (8 ((-11 + 15 l + 15 l^2)^2 \[Pi]^2) \[Omega]^4)/((1 - 2 l)^2 (1 + 2 l)^2 (3 + 2 l)^2);
    ,
    Cos2\[Pi]\[Nu] = 1 - (8 (\[Pi]^2 (30 l^3 + 15 l^4 + 3 s^2 (-2 + s^2) + l (-11 + 6 s^2) + l^2 (4 + 6 s^2))^2) \[Omega]^4)/((1 - 2 l)^2 l^2 (1 + l)^2 (1 + 2 l)^2 (3 + 2 l)^2);
  ];
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
  res = Check[Which[
    -1 <= Cos2\[Pi]\[Nu] <= 1, (* \[Nu] Real *)
      (* FIXME: The order of the arguments assumes \[Nu]0 is positive *)
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
  ];
  Clear[\[Alpha]\[Gamma], \[Beta], R, L];
  res
];

(* Estimate precision of \[Nu] based on how well the continued fraction equation is satisfied. *)
\[Nu]precision[Cos2\[Pi]\[Nu]_, q_, \[Epsilon]_Complex, \[Kappa]_, \[Tau]_, s_, \[Lambda]_, m_] :=
 Module[{\[Alpha]\[Gamma], \[Beta], R, L, \[Nu]0, prec},
  \[Alpha]\[Gamma][n_, \[Nu]_?InexactNumberQ] := \[Epsilon]^2 \[Kappa]^2 (n + \[Nu]) (2 + n + \[Nu]) ((1 + n + \[Nu] - s)^2 + \[Epsilon]^2) ((1 + n + \[Nu] + s)^2 + \[Epsilon]^2) (-1 + 2 n + 2 \[Nu]) (5 + 2 n + 2 \[Nu]) ((1 + n + \[Nu])^2 + \[Tau]^2);
  \[Beta][n_, \[Nu]_?InexactNumberQ] := (2 n + 2 \[Nu] + 3) (2 n + 2 \[Nu] - 1) ((-\[Lambda] - s (s + 1) + (n + \[Nu]) (n + \[Nu] + 1) + \[Epsilon]^2 + \[Epsilon] (\[Epsilon] - m q)) ((n + \[Nu]) (n + \[Nu] + 1)) + (\[Epsilon] (\[Epsilon] - m q) (s^2 + \[Epsilon]^2)));
  R[n_, \[Nu]_] := Module[{i}, Teukolsky`MST`MST`Private`CF[-\[Alpha]\[Gamma][i-1, \[Nu]], \[Beta][i, \[Nu]], {i, n}]];
  L[n_, \[Nu]_] := Module[{i}, Teukolsky`MST`MST`Private`CF[-\[Alpha]\[Gamma][2n-i, \[Nu]], \[Beta][2n-i, \[Nu]], {i, n}]];
  prec = With[{\[Nu] = ArcCos[Cos2\[Pi]\[Nu]]/(2\[Pi])}, -RealExponent[\[Beta][0, \[Nu]] + R[1, \[Nu]] + L[-1, \[Nu]]]];
  Clear[\[Alpha]\[Gamma], \[Beta], R, L];
  prec
];

(* Estimate precision of \[Nu] based on the complex part of Cos[2 \[Pi] \[Nu]]. This is only valid
   for real frequencies where Cos[2 \[Pi] \[Nu]] is expected to be purely real. *)
\[Nu]precision[Cos2\[Pi]\[Nu]_, q_, \[Epsilon]_, \[Kappa]_, \[Tau]_, s_, \[Lambda]_, m_] := -RealExponent[Im[Cos2\[Pi]\[Nu]]/Re[Cos2\[Pi]\[Nu]]];

(* Find \[Nu] using monodromy of confluent Heun equation *)
\[Nu]RCHMonodromy[a_, \[Omega]_, \[Lambda]_, s_, l_, m_, Npmax_] :=
 Module[{q, \[Epsilon], \[Kappa], \[Tau], \[Gamma]CH, \[Delta]CH, \[Epsilon]CH, \[Alpha]CH\[Epsilon]CH, qCH, \[Mu]1C, \[Mu]2C, a1, a2, a1sum, a2sum, Cos2\[Pi]\[Nu], extend, nmax, nmin, precision, \[Nu], iterations = 0},
  q = a;
  \[Epsilon] = 2 \[Omega];
  \[Kappa] = Sqrt[1-q^2];
  \[Tau] = (\[Epsilon]-m q)/\[Kappa];

  (* Compute parameters in confluent Heun equation *)
  \[Gamma]CH = 1-s-I \[Epsilon]-I \[Tau];
  \[Delta]CH = 1+s+I \[Epsilon]-I \[Tau];
  \[Epsilon]CH = 2I \[Epsilon] \[Kappa];
  \[Alpha]CH\[Epsilon]CH =(1-s+I(\[Epsilon]-\[Tau]));
  qCH = -(-s (1+s)+\[Epsilon]^2+I (-1+2 s) \[Epsilon] \[Kappa]-\[Lambda]-\[Tau] (I+\[Tau]));

  \[Mu]1C = \[Alpha]CH\[Epsilon]CH-(\[Gamma]CH+\[Delta]CH);
  \[Mu]2C = -\[Alpha]CH\[Epsilon]CH;

  (* Recurrence relations *)
  a1[-1] = 0;
  a1[0] = 1;
  a1[n_] := a1[n] = (((\[Alpha]CH\[Epsilon]CH-(-1+n+\[Delta]CH)) (\[Alpha]CH\[Epsilon]CH-(-2+n+\[Gamma]CH+\[Delta]CH))\[Epsilon]CH) a1[n-2])/n-((\[Alpha]CH\[Epsilon]CH^2+\[Alpha]CH\[Epsilon]CH (1-2 n-\[Gamma]CH-\[Delta]CH+\[Epsilon]CH)+(n^2-qCH+n (-1+\[Gamma]CH+\[Delta]CH-\[Epsilon]CH)+\[Epsilon]CH-\[Delta]CH \[Epsilon]CH)) a1[n-1])/n ;

  a2[-1] = 0;
  a2[0] = 1;
  a2[n_] := a2[n] = -((((\[Alpha]CH\[Epsilon]CH+(-2+n)) (\[Alpha]CH\[Epsilon]CH+(-1+n-\[Gamma]CH))\[Epsilon]CH) a2[n-2])/n)+((\[Alpha]CH\[Epsilon]CH^2+(n^2-qCH+\[Gamma]CH+\[Delta]CH-n (1+\[Gamma]CH+\[Delta]CH-\[Epsilon]CH)-\[Epsilon]CH)+\[Alpha]CH\[Epsilon]CH(-1+2 n-\[Gamma]CH-\[Delta]CH+\[Epsilon]CH)) a2[n-1])/n;

  (* The sums are Gamma[mu1 - mu2] Sum[a1[j] Pochhammer[mu1 - mu2, n - j]] and its mirror image, written
     with Gamma[mu1 - mu2] Pochhammer[mu1 - mu2, k] = Gamma[mu1 - mu2 + k]: at the degeneracies
     2 I epsilon = integer (omega = -I n/(4M)) the difference mu1 - mu2 = 2 I epsilon - 2 s is an integer,
     Gamma[mu1 - mu2] or Gamma[mu2 - mu1] has a pole and the Pochhammer symbol the cancelling zero, so the
     combined form is finite there, provided n - j + |mu1 - mu2| > 0 for every term (nmax > 2 |mu1 - mu2|,
     enforced below), and nu is evaluated directly instead of from neighbouring frequencies. *)
  a1sum[n_] := Sum[a1[j] Gamma[\[Mu]1C-\[Mu]2C+n-j], {j, 0, Ceiling[n/2]}];
  a2sum[n_] := Sum[(-1)^j a2[j] Gamma[\[Mu]2C-\[Mu]1C+n-j], {j, 0, Ceiling[n/2]}];
  nmin = 2 Ceiling[Abs[\[Mu]1C-\[Mu]2C]] + 4;

  (* Fill the memoised tables up to n in increasing order, so that the recursion depth stays at one
     whatever nmax is (evaluating a1[nmax] directly recursed nmax deep and exceeded $RecursionLimit) *)
  extend[n_] := Do[a1[i]; a2[i], {i, 1, n}];

  (* Compute \[Nu]. The memoised tables are cleared on every path, including failure. *)
  Cos2\[Pi]\[Nu][nmax_] := Cos2\[Pi]\[Nu][nmax] = (extend[nmax]; Cos[\[Pi](\[Mu]1C-\[Mu]2C)]+(2\[Pi]^2)/(a1sum[nmax] a2sum[nmax]) (-1)^(nmax-1) a1[nmax]a2[nmax]);
  \[Nu] = Catch[
  If[IntegerQ[Npmax],
    nmax = Max[Npmax, nmin];
    If[Precision[Cos2\[Pi]\[Nu][nmax]] == 0, Throw[$Failed, \[Nu]RCHMonodromy]];
  ,
    (* FIXME: we should be able to predict nmax based on the convergence for large nmax and the loss of precision in a1 and a2 *)
    nmax = Max[2 Ceiling[E^ProductLog[Precision[{a, \[Omega], \[Lambda]}] Log[100]]], nmin];
    If[Precision[Cos2\[Pi]\[Nu][nmax]] == 0, Throw[$Failed, \[Nu]RCHMonodromy]];

    (* Increase nmax by 10% until the precision of the result decreases; a bounded number of times, so
       that a non-convergent case fails instead of exhausting memory *)
    precision = -Infinity;
    While[precision < (precision = \[Nu]precision[Cos2\[Pi]\[Nu][nmax], q, \[Epsilon], \[Kappa], \[Tau], s, \[Lambda], m]),
      nmax = Round[11/10 nmax];
      If[++iterations > 25, Message[RenormalizedAngularMomentum::conv, \[Omega]]; Throw[$Failed, \[Nu]RCHMonodromy]];
      If[Precision[Cos2\[Pi]\[Nu][nmax]] == 0, Throw[$Failed, \[Nu]RCHMonodromy]];
    ];
    nmax = Round[10/11 nmax];
  ];
    
  If[Precision[Cos2\[Pi]\[Nu][nmax]]=!=MachinePrecision,
    Cos2\[Pi]\[Nu][nmax] = N[Cos2\[Pi]\[Nu][nmax], Max[\[Nu]precision[Cos2\[Pi]\[Nu][nmax], q, \[Epsilon], \[Kappa], \[Tau], s, \[Lambda], m],0]];
  ];

  (* The representative of nu (the class {+-nu + k}): the one continuous with nu = l + O(epsilon^2) at
     small frequency, l - ArcCos[Cos[2 Pi nu]]/(2 Pi), for real and complex frequencies alike (the principal
     value ArcCos[...]/(2 Pi), returned before for complex frequencies, is near 0 at small frequency, where
     Pochhammer[2 nu + 2, n] in the K_nu sums of the amplitudes has an exact pole). *)
  Which[
    Im[\[Omega]] != 0,
      l - ArcCos[Cos2\[Pi]\[Nu][nmax]]/(2\[Pi]),
    Re[Cos2\[Pi]\[Nu][nmax]]<-1, 
      1/2-Im[ArcCos[Re[Cos2\[Pi]\[Nu][nmax]]]/(2\[Pi])]I,
    -1<=Re[Cos2\[Pi]\[Nu][nmax]]<=1,
      l-ArcCos[Re[Cos2\[Pi]\[Nu][nmax]]]/(2\[Pi]),
    Re[Cos2\[Pi]\[Nu][nmax]]>1,
      -I Im[ArcCos[Re[Cos2\[Pi]\[Nu][nmax]]]/(2\[Pi])],
    True,
      $Failed
  ], \[Nu]RCHMonodromy];
  (* Remove rather than Clear: a message issued during the evaluation (1/0 near a degeneracy) can keep a
     reference to a memoised table, which would then survive as a leaked symbol *)
  Remove[a1, a2, a1sum, a2sum, Cos2\[Pi]\[Nu], extend];
  \[Nu]
];

(**********************************************************)
(* RenormalizedAngularMomentum                            *)
(**********************************************************)

SyntaxInformation[RenormalizedAngularMomentum] =
 {"ArgumentsPattern" -> {_, _, _, _, _, ___}};
Options[RenormalizedAngularMomentum] = {Method -> "Monodromy"};
SetAttributes[RenormalizedAngularMomentum, {NumericFunction}];

RenormalizedAngularMomentum[s_, l_, m_, a_, \[Omega]_, \[Lambda]_, OptionsPattern[RenormalizedAngularMomentum]] /; l < Abs[s] := 0;

RenormalizedAngularMomentum[s_, l_, m_, a_, \[Omega]_?PossibleZeroQ, \[Lambda]_, OptionsPattern[RenormalizedAngularMomentum]] := l;

RenormalizedAngularMomentum[s_Integer, l_Integer, m_Integer, a_?NumericQ, \[Omega]_?NumericQ, \[Lambda]_?NumericQ,
 Method -> {"FindRoot", "InitialGuess" -> \[Nu]_}] /; InexactNumberQ[a] || InexactNumberQ[\[Omega]] || InexactNumberQ[\[Lambda]] :=
  \[Nu]RootFind[a, \[Omega], \[Lambda], s, l, m, Cos[2 \[Pi] \[Nu]]];

RenormalizedAngularMomentum[s_Integer, l_Integer, m_Integer, a_?NumericQ, \[Omega]_?NumericQ, \[Lambda]_?NumericQ, Method -> ("FindRoot"|{"FindRoot"})] /; InexactNumberQ[a] || InexactNumberQ[\[Omega]] || InexactNumberQ[\[Lambda]] :=
 Module[{\[Nu], Cos2\[Pi]\[Nu]},
  Cos2\[Pi]\[Nu] = Cos2\[Pi]\[Nu]Series[N[a], N[\[Omega]], s, l, m];
  If[Cos2\[Pi]\[Nu] == $Failed, Return[$Failed]];
  \[Nu] = \[Nu]RootFind[a, \[Omega], \[Lambda], s, l, m, Cos2\[Pi]\[Nu]];

  \[Nu]
];

RenormalizedAngularMomentum[s_Integer, l_Integer, m_Integer, a_?NumericQ, \[Omega]_?NumericQ, \[Lambda]_?NumericQ,
 Method -> ("Monodromy"|{"Monodromy"})] /; InexactNumberQ[a] || InexactNumberQ[\[Omega]] || InexactNumberQ[\[Lambda]] :=
  RenormalizedAngularMomentum[s, l, m, a, \[Omega], \[Lambda], Method -> {"Monodromy", "nmax" -> Automatic}];

RenormalizedAngularMomentum[s_Integer, l_Integer, m_Integer, a_?NumericQ, \[Omega]_?NumericQ, \[Lambda]_?NumericQ, 
 Method -> {"Monodromy", "nmax" -> nmax_}] /; InexactNumberQ[a] || InexactNumberQ[\[Omega]] || InexactNumberQ[\[Lambda]] :=
 Module[{\[Nu]},
  If[AnyTrue[{a, \[Omega], \[Lambda]}, MachineNumberQ],
    (* Print a warning if run with machine precision input *)
    Message[RenormalizedAngularMomentum::precision];
  ];
  \[Nu] = \[Nu]RCHMonodromy[a, \[Omega], \[Lambda], s, l, m, nmax]
];

RenormalizedAngularMomentum[s_Integer, l_Integer, m_Integer, a_?NumericQ, \[Omega]_?NumericQ, \[Lambda]_?NumericQ, 
 Method -> ("Series"|{"Series"})] /; InexactNumberQ[a] || InexactNumberQ[\[Omega]] || InexactNumberQ[\[Lambda]] :=
  l - ArcCos[Cos2\[Pi]\[Nu]Series[a, \[Omega], s, l, m]] / (2*\[Pi]);

RenormalizedAngularMomentum[s_Integer, l_Integer, m_Integer, a_?NumericQ, \[Omega]_?NumericQ,
 Method -> ("Series"|{"Series"})] /; InexactNumberQ[a] || InexactNumberQ[\[Omega]] || InexactNumberQ[\[Lambda]] :=
  l - ArcCos[Cos2\[Pi]\[Nu]Series[a, \[Omega], s, l, m]] / (2*\[Pi]);

RenormalizedAngularMomentum[s_Integer, l_Integer, m_Integer, a_?NumericQ, \[Omega]_?NumericQ, \[Lambda]_?NumericQ] /;
 InexactNumberQ[a] || InexactNumberQ[\[Omega]] || InexactNumberQ[\[Lambda]] :=
   RenormalizedAngularMomentum[s, l, m, a, \[Omega], \[Lambda], Method -> "Monodromy"];

RenormalizedAngularMomentum[s_Integer, l_Integer, m_Integer, a_?NumericQ, \[Omega]_?NumericQ, opts:OptionsPattern[RenormalizedAngularMomentum]] :=
  RenormalizedAngularMomentum[s, l, m, a, \[Omega], SpinWeightedSpheroidalEigenvalue[s, l, m, a \[Omega]], opts];

(* With exact arguments the five-argument form rewrites itself in terms of the (unevaluated, exact)
   spheroidal eigenvalue, so N must accept a non-numeric eigenvalue and apply itself to it. *)
N[RenormalizedAngularMomentum[s_Integer, l_Integer, m_Integer, a_?NumericQ, \[Omega]_?NumericQ, \[Lambda]_, opts:OptionsPattern[RenormalizedAngularMomentum]], p_:MachinePrecision] ^:=
  RenormalizedAngularMomentum[s, l, m, N[a, p], N[\[Omega], p], N[\[Lambda], p], opts];

SetAttributes[RenormalizedAngularMomentum, {Protected, ReadProtected}];

End[];
EndPackage[];
