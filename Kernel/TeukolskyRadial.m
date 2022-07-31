(* ::Package:: *)

(* ::Title:: *)
(*TeukolskyRadial*)


(* ::Section::Closed:: *)
(*Create Package*)


(* ::Subsection::Closed:: *)
(*BeginPackage*)


BeginPackage["Teukolsky`TeukolskyRadial`",
  {
  "Teukolsky`SasakiNakamura`",
   "Teukolsky`NumericalIntegration`",
   "Teukolsky`MST`RenormalizedAngularMomentum`",
   "Teukolsky`MST`MST`",
   "SpinWeightedSpheroidalHarmonics`"
  }
];


(* ::Subsection::Closed:: *)
(*Unprotect symbols*)


ClearAttributes[{TeukolskyRadial, TeukolskyRadialFunction}, {Protected, ReadProtected}];


(* ::Subsection::Closed:: *)
(*Usage messages*)


TeukolskyRadial::usage = "TeukolskyRadial[s, l, m, a, \[Omega]] computes homogeneous solutions to the radial Teukolsky equation."
TeukolskyRadialFunction::usage = "TeukolskyRadialFunction[s, l, m, a, \[Omega], assoc] is an object representing a homogeneous solution to the radial Teukolsky equation."


(* ::Subsection::Closed:: *)
(*Error Messages*)


TeukolskyRadial::precw = "The precision of `1`=`2` is less than WorkingPrecision (`3`).";
TeukolskyRadial::optx = "Unknown options in `1`";
TeukolskyRadial::dm = "Option `1` is not valid with BoundaryConditions \[RightArrow] `2`.";
TeukolskyRadial::sopt = "Option `1` not supported for static (\[Omega]=0) modes.";
TeukolskyRadial::hc = "Method HeunC is only supported with Mathematica version 12.1 and later.";
TeukolskyRadial::hcopt = "Option `1` not supported for HeunC method.";
TeukolskyRadialFunction::dmval = "Radius `1` lies outside the computational domain. Results may be incorrect.";


(* ::Subsection::Closed:: *)
(*Begin Private section*)


Begin["`Private`"];


(* ::Section::Closed:: *)
(*Utility Functions*)


(* ::Subsection::Closed:: *)
(*Horizon Locations*)


rp[a_,M_] := M+Sqrt[M^2-a^2];
rm[a_,M_] := M-Sqrt[M^2-a^2];


(* ::Subsection::Closed:: *)
(*Hyperboloidal Transformation Functions*)


f[r_] := 1-2/r;
rs[r_,a_]:=r+2/(rp[a,1]-rm[a,1]) (rp[a,1] Log[(r-rp[a,1])/2]-rm[a,1] Log[(r-rm[a,1])/2]);
\[CapitalDelta][r_,a_]:=r^2+a^2-2r;
\[Phi]Reg[r_,a_]:=a /(rp[a,1]-rm[a,1]) Log[(r-rp[a,1])/(r-rm[a,1])];


(* ::Section::Closed:: *)
(*TeukolskyRadial*)


(* ::Subsection::Closed:: *)
(*Numerical Integration Method*)


Options[TeukolskyRadialNumericalIntegration] = Join[
  {"Domain" -> All},
  FilterRules[Options[NDSolve], Except[WorkingPrecision|AccuracyGoal|PrecisionGoal]]];


domainQ[domain_] := MatchQ[domain, {_?NumericQ, _?NumericQ} | (_?NumericQ) | All];


TeukolskyRadialNumericalIntegration[s_Integer, l_Integer, m_Integer, a_, \[Omega]_, BCs_, {wp_, prec_, acc_}, opts:OptionsPattern[]] :=
 Module[{\[Lambda], TRF, ndsolveopts, solFuncs, domains},
  (* Compute the eigenvalue *)
  \[Lambda] = SpinWeightedSpheroidalEigenvalue[s, l, m, a \[Omega]];
  
  (* Function to construct a single TeukolskyRadialFunction *)
  TRF[bc_, sf_, domain_, ndsolveopts___] :=
   Module[{solutionFunction, bcdir},
    solutionFunction = sf[domain];
    bcdir = bc /. {"In" -> -1, "Up" -> +1};
    TeukolskyRadialFunction[s, l, m, a, \[Omega],
     Association["s" -> s, "l" -> l, "m" -> m, "a" -> a, "\[Omega]" -> \[Omega], "Eigenvalue" -> \[Lambda],
      "Method" -> {"NumericalIntegration", ndsolveopts},
      "BoundaryConditions" -> bc, "Amplitudes" -> solutionFunction[[2]],
      "Domain" -> If[domain === All, {rp[a, 1], \[Infinity]}, First[solutionFunction[[1]]["Domain"]]],
      "RadialFunction" -> Function[{Global`r}, Evaluate[Global`r^-1 \[CapitalDelta][Global`r,a]^-s Exp[bcdir I \[Omega] rs[Global`r,a]] Exp[I m \[Phi]Reg[Global`r,a]] solutionFunction[[1]][Global`r]]]
     ]
    ]
   ];

  (* Domain over which the numerical solution can be evaluated *)
  domains = OptionValue["Domain"];
  If[ListQ[BCs],
    If[domains === All, domains = Thread[BCs -> All]];
    If[!MatchQ[domains, (List|Association)[Rule["In"|"Up",_?domainQ]..]],
      Message[TeukolskyRadial::dm, "Domain" -> domains, BCs];
      Return[$Failed];
    ];
    domains = Lookup[domains, BCs, None]; 
    If[!AllTrue[domains, domainQ],
      Message[TeukolskyRadial::dm, "Domain" -> OptionValue["Domain"], BCs];
      Return[$Failed];
    ];
  ,
    If[!domainQ[domains],
      Message[TeukolskyRadial::dm, "Domain" -> domains, BCs];
      Return[$Failed];
    ];
  ];
  
  (* Solution functions for the specified boundary conditions *)
  ndsolveopts = Sequence@@FilterRules[{opts}, Options[NDSolve]];
  solFuncs =
   <|"In" :> Teukolsky`NumericalIntegration`Private`psi[s, \[Lambda], l, m, a, \[Omega], "In", WorkingPrecision -> wp, PrecisionGoal -> prec, AccuracyGoal -> acc, ndsolveopts],
     "Up" :> Teukolsky`NumericalIntegration`Private`psi[s, \[Lambda], l, m, a, \[Omega], "Up", WorkingPrecision -> wp, PrecisionGoal -> prec, AccuracyGoal -> acc, ndsolveopts]|>;
  solFuncs = Lookup[solFuncs, BCs];

  If[ListQ[BCs],
    Return[Association[MapThread[#1 -> TRF[#1, #2, #3, ndsolveopts]&, {BCs, solFuncs, domains}]]],
    Return[TRF[BCs, solFuncs, domains, ndsolveopts]]
  ];
];


(* ::Subsection::Closed:: *)
(*Sasaki-Nakamura Method*)


Options[TeukolskyRadialSasakiNakamura] = Join[
  {"Domain" -> None},
  FilterRules[Options[NDSolve], Except[WorkingPrecision|AccuracyGoal|PrecisionGoal]]];


domainQ[domain_] := MatchQ[domain, {_?NumericQ, _?NumericQ} | (_?NumericQ) | All];


TeukolskyRadialSasakiNakamura[s_Integer, l_Integer, m_Integer, a_, \[Omega]_, BCs_, {wp_, prec_, acc_}, opts:OptionsPattern[]] :=
 Module[{\[Lambda], TRF, norms, ndsolveopts, solFuncs, domains},
  (* Compute the eigenvalue *)
  \[Lambda] = SpinWeightedSpheroidalEigenvalue[s, l, m, a \[Omega]];

  (* Function to construct a single TeukolskyRadialFunction *)
  TRF[bc_, ns_, sf_, domain_, ndsolveopts___] :=
   Module[{solutionFunction},
    If[sf === $Failed, Return[$Failed]];
    solutionFunction = sf[domain];
    TeukolskyRadialFunction[s, l, m, a, \[Omega],
     Association["s" -> s, "l" -> l, "m" -> m, "a" -> a, "\[Omega]" -> \[Omega], "Eigenvalue" -> \[Lambda],
      "Method" -> {"SasakiNakamura", ndsolveopts},
      "BoundaryConditions" -> bc, "Amplitudes" -> ns,
      "Domain" -> If[domain === All, {rp[a, 1], \[Infinity]}, First[solutionFunction["Domain"]]],
      "RadialFunction" -> solutionFunction
     ]
    ]
   ];

  (* Domain over which the numerical solution can be evaluated *)
  domains = OptionValue["Domain"];
  If[ListQ[BCs],
    If[!MatchQ[domains, (List|Association)[Rule["In"|"Up",_?domainQ]..]],
      Message[TeukolskyRadial::dm, "Domain" -> domains, BCs];
      Return[$Failed];
    ];
    domains = Lookup[domains, BCs, None]; 
    If[!AllTrue[domains, domainQ],
      Message[TeukolskyRadial::dm, "Domain" -> OptionValue["Domain"], BCs];
      Return[$Failed];
    ];
  ,
    If[!domainQ[domains],
      Message[TeukolskyRadial::dm, "Domain" -> domains, BCs];
      Return[$Failed];
    ];
  ];

  (* Asymptotic normalizations such that we have unit transmission coefficient *)
  norms = <|"In" -> <|"Transmission" -> 1|>, "Up" -> <|"Transmission" -> 1|>|>;
  norms = Lookup[norms, BCs];

  (* Solution functions for the specified boundary conditions *)
  ndsolveopts = Sequence@@FilterRules[{opts}, Options[NDSolve]];
  solFuncs =
   <|"In" :> $Failed,
     "Up" :> Teukolsky`SasakiNakamura`Private`TeukolskyRadialUp[s, \[Lambda], m, a, \[Omega](*, WorkingPrecision -> wp, PrecisionGoal -> prec, AccuracyGoal -> acc, ndsolveopts*)]
     |>;
  solFuncs = Lookup[solFuncs, BCs];

  If[ListQ[BCs],
    Return[Association[MapThread[#1 -> TRF[#1, #2, #3, #4, ndsolveopts]&, {BCs, norms, solFuncs, domains}]]],
    Return[TRF[BCs, norms, solFuncs, domains, ndsolveopts]]
  ];
];


(* ::Subsection::Closed:: *)
(*MST Method*)


Options[TeukolskyRadialMST] = {
  "RenormalizedAngularMomentum" -> "Monodromy"};


TeukolskyRadialMST[s_Integer, l_Integer, m_Integer, a_, \[Omega]_, BCs_, {wp_, prec_, acc_}, opts:OptionsPattern[]] :=
 Module[{\[Lambda], \[Nu], norms, solFuncs, TRF},
  (* Compute the eigenvalue and renormalized angular momentum *)
  \[Lambda] = SpinWeightedSpheroidalEigenvalue[s, l, m, a \[Omega]];
  \[Nu] = RenormalizedAngularMomentum[s, l, m, a, \[Omega], \[Lambda], Method -> OptionValue["RenormalizedAngularMomentum"]];

  (* Function to construct a TeukolskyRadialFunction *)
  TRF[bc_, ns_, sf_] :=
    TeukolskyRadialFunction[s, l, m, a, \[Omega],
     Association["s" -> s, "l" -> l, "m" -> m, "a" -> a, "\[Omega]" -> \[Omega], "Eigenvalue" -> \[Lambda],
      "Method" -> {"MST", "RenormalizedAngularMomentum" -> \[Nu]},
      "BoundaryConditions" -> bc, "Amplitudes" -> ns,
      "Domain" -> {rp[a, 1], \[Infinity]}, "RadialFunction" -> sf
     ]
    ];

  (* Compute the asymptotic normalisations *)
  norms = Teukolsky`MST`MST`Private`Amplitudes[s, l, m, a, 2\[Omega], \[Nu], \[Lambda], {wp, prec, acc}];

  (* Solution functions for the specified boundary conditions *)
  solFuncs =
    <|"In" :> Teukolsky`MST`MST`Private`MSTRadialIn[s,l,m,a,2\[Omega],\[Nu],\[Lambda],norms["In"]["Transmission"], {wp, prec, acc}],
      "Up" :> Teukolsky`MST`MST`Private`MSTRadialUp[s,l,m,a,2\[Omega],\[Nu],\[Lambda],norms["Up"]["Transmission"], {wp, prec, acc}]|>;
  solFuncs = Lookup[solFuncs, BCs];

  (* Select normalisation coefficients for the specified boundary conditions and rescale
     to give unit transmission coefficient. *)
  norms = norms/norms[[All, "Transmission"]];
  norms = Lookup[norms, BCs];

  If[ListQ[BCs],
    Return[Association[MapThread[#1 -> TRF[#1, #2, #3]&, {BCs, norms, solFuncs}]]],
    Return[TRF[BCs, norms, solFuncs]]
  ];
];


(* ::Subsection::Closed:: *)
(*HeunC*)


Options[TeukolskyRadialHeunC] = {};


TeukolskyRadialHeunC[s_Integer, l_Integer, m_Integer, a_, \[Omega]_, BCs_, {wp_, prec_, acc_}, opts:OptionsPattern[]] :=
 Module[{\[Lambda], \[Nu], norms, solFuncs, TRF},
  (* The HeunC method is only supported on version 12.1 and newer *)
  If[$VersionNumber < 12.1,
    Message[TeukolskyRadial::hc];
    Return[$Failed]
  ];
 
  (* Compute the eigenvalue and renormalized angular momentum *)
  \[Lambda] = SpinWeightedSpheroidalEigenvalue[s, l, m, a \[Omega]];

  (* Function to construct a TeukolskyRadialFunction *)
  TRF[bc_, ns_, sf_] :=
    If[sf === $Failed, $Failed,
      TeukolskyRadialFunction[s, l, m, a, \[Omega],
        Association["s" -> s, "l" -> l, "m" -> m, "a" -> a, "\[Omega]" -> \[Omega], "Eigenvalue" -> \[Lambda],
          "Method" -> {"HeunC"},
          "BoundaryConditions" -> bc, "Amplitudes" -> ns,
          "Domain" -> {rp[a, 1], \[Infinity]}, "RadialFunction" -> sf
        ]
      ]
    ];

  (* Compute the asymptotic normalisations *)
  norms = <|"In" -> <|"Transmission" -> 1|>, "Up" -> <|"Transmission" -> 1|>|>;

  (* Solution functions for the specified boundary conditions *)
  solFuncs =
    <|"In" :> Function[{r}, 2^((-I a m-Sqrt[1-a^2] s+2 I \[Omega])/Sqrt[1-a^2]) (1-a^2)^((I a m)/(2 (1+Sqrt[1-a^2]))-s-I \[Omega]) E^(1/2 I (a m-2 r \[Omega])) ((-1-Sqrt[1-a^2]+r)/Sqrt[1-a^2])^((I a m)/(2 Sqrt[1-a^2])-s-I (1+1/Sqrt[1-a^2]) \[Omega]) ((-1+Sqrt[1-a^2]+r)/Sqrt[1-a^2])^((I (a m+2 (-1+Sqrt[1-a^2]) \[Omega]))/(2 Sqrt[1-a^2])) HeunC[s+s^2+\[Lambda]-(a m-2 \[Omega])^2/(-1+a^2)-4 \[Omega]^2+(I (-a m-4 (-1+s) \[Omega]+2 a^2 (-1+2 s) \[Omega]))/Sqrt[1-a^2],-4 \[Omega] (-I Sqrt[1-a^2]+a m+I Sqrt[1-a^2] s-2 \[Omega]+2 Sqrt[1-a^2] \[Omega]),1-s+(I (a m-2 \[Omega]))/Sqrt[1-a^2]-2 I \[Omega],1+s+(I (a m-2 \[Omega]))/Sqrt[1-a^2]+2 I \[Omega],4 I Sqrt[1-a^2] \[Omega],(1+Sqrt[1-a^2]-r)/(2 Sqrt[1-a^2])]],
      "Up" :> $Failed |>;
  solFuncs = Lookup[solFuncs, BCs];

  (* Select normalisation coefficients for the specified boundary conditions and rescale
     to give unit transmission coefficient. *)
  norms = norms[[All, {"Transmission"}]]/norms[[All, "Transmission"]];
  norms = Lookup[norms, BCs];

  If[ListQ[BCs],
    Return[Association[MapThread[#1 -> TRF[#1, #2, #3]&, {BCs, norms, solFuncs}]]],
    Return[TRF[BCs, norms, solFuncs]]
  ];
];


(* ::Subsection::Closed:: *)
(*Static modes*)


TeukolskyRadialStatic[s_Integer, l_Integer, m_Integer, a_, \[Omega]_, BCs_] :=
 Module[{\[Lambda], norms, normIn, normUp, solFuncs, TRF},
  (* Compute the eigenvalue *)
  \[Lambda] = SpinWeightedSpheroidalEigenvalue[s, l, m, a \[Omega]];

  (* Function to construct a TeukolskyRadialFunction *)
  TRF[bc_, ns_, sf_] :=
    TeukolskyRadialFunction[s, l, m, a, \[Omega],
     Association["s" -> s, "l" -> l, "m" -> m, "a" -> a, "\[Omega]" -> \[Omega], "Eigenvalue" -> \[Lambda],
      "Method" -> {"Static"},
      "BoundaryConditions" -> bc, "Amplitudes" -> ns,
      "Domain" -> {rp[a, 1], \[Infinity]}, "RadialFunction" -> sf
     ]
    ];

  (* Compute the asymptotic normalisations *)
  norms = <|"In" -> <|"Transmission" -> 1|>, "Up" -> <|"Transmission" -> 1|>|>;

  (* Solution functions for the specified boundary conditions *)
  With[{\[Tau] = -((m a)/Sqrt[1-a^2]), \[Kappa] = Sqrt[1 - a^2]},
    normIn = (2 \[Kappa])^(-2s) Exp[-I \[Tau]/2 \[Kappa] (1+2Log[\[Kappa]]/(1+\[Kappa]))] If[s >= 1 && \[Tau] == 0, 1, Gamma[1 - s - I \[Tau]]];
    normUp = (2 \[Kappa])^(-s-l-1);
    solFuncs =
      <|"In" :> Function[{r}, normIn (-(1 + \[Kappa] - r)/(2 \[Kappa]))^(-s - I \[Tau]/2) (1 - (1 + \[Kappa] - r)/(2 \[Kappa]))^(-I \[Tau]/2) Hypergeometric2F1Regularized[-l - I \[Tau], l + 1 - I \[Tau], 1 - s - I \[Tau], (1 + \[Kappa] - r)/(2 \[Kappa])]],
      "Up" :> Function[{r}, normUp (-(1 + \[Kappa] - r)/(2 \[Kappa]))^(-s - (I \[Tau])/2) (1 - (1 + \[Kappa] - r)/(2 \[Kappa]))^((I \[Tau])/2 - l - 1) Hypergeometric2F1[l + 1 - I \[Tau], l + 1 - s, 2 l + 2, 1/(1 - (1 + \[Kappa] - r)/(2 \[Kappa]))]]
       |>;
  ];
  solFuncs = Lookup[solFuncs, BCs];

  (* Select normalisation coefficients for the specified boundary conditions and rescale
     to give unit transmission coefficient. *)
  norms = norms/norms[[All, "Transmission"]];
  norms = Lookup[norms, BCs];

  If[ListQ[BCs],
    Return[Association[MapThread[#1 -> TRF[#1, #2, #3]&, {BCs, norms, solFuncs}]]],
    Return[TRF[BCs, norms, solFuncs]]
  ];
];


(* ::Subsection::Closed:: *)
(*TeukolskyRadial*)


SyntaxInformation[TeukolskyRadial] =
 {"ArgumentsPattern" -> {_, _, _, _, _, OptionsPattern[]}};


Options[TeukolskyRadial] = {
  Method -> Automatic,
  "BoundaryConditions" -> {"In", "Up"},
  WorkingPrecision -> Automatic,
  PrecisionGoal -> Automatic,
  AccuracyGoal -> Automatic
};


(* ::Subsubsection::Closed:: *)
(*Static modes*)


TeukolskyRadial[s_Integer, l_Integer, m_, a_, \[Omega]_, opts:OptionsPattern[]] /; \[Omega] == 0 :=
 Module[{BCs, wp, prec, acc},
  (* Determine which boundary conditions the homogeneous solution(s) should satisfy *)
  BCs = OptionValue["BoundaryConditions"];
  If[!MatchQ[BCs, "In"|"Up"|{("In"|"Up")..}], 
    Message[TeukolskyRadial::optx, "BoundaryConditions" -> BCs];
    Return[$Failed];
  ];

  (* Options are not supported for static modes *)
  Do[
    If[OptionValue[opt] =!= Automatic, Message[TeukolskyRadial::sopt, opt]];,
    {opt, {Method, WorkingPrecision, PrecisionGoal, AccuracyGoal}}
  ];

  (* Call the chosen implementation *)
  TeukolskyRadialStatic[s, l, m, a, \[Omega], BCs]
]


(* ::Subsubsection::Closed:: *)
(*Non-static modes*)


TeukolskyRadial[s_Integer, l_Integer, m_Integer, a_, \[Omega]_, opts:OptionsPattern[]] /; InexactNumberQ[a] || InexactNumberQ[\[Omega]] :=
 Module[{TRF, subopts, BCs, wp, prec, acc},
  (* Extract suboptions from Method to be passed on. *)
  If[ListQ[OptionValue[Method]],
    subopts = Rest[OptionValue[Method]];,
    subopts = {};
  ];

  (* Determine which boundary conditions the homogeneous solution(s) should satisfy *)
  BCs = OptionValue["BoundaryConditions"];
  If[!MatchQ[BCs, "In"|"Up"|{("In"|"Up")..}], 
    Message[TeukolskyRadial::optx, "BoundaryConditions" -> BCs];
    Return[$Failed];
  ];

  (* Options associated with precision and accuracy *)
  {wp, prec, acc} = OptionValue[{WorkingPrecision, PrecisionGoal, AccuracyGoal}];
  If[wp === Automatic, wp = Precision[{a, \[Omega]}]];
  If[prec === Automatic, prec = wp / 2];
  If[acc === Automatic, acc = Infinity];
  If[Precision[a] < wp, Message[TeukolskyRadial::precw, "a", a, wp]];
  If[Precision[\[Omega]] < wp, Message[TeukolskyRadial::precw, "\[Omega]", \[Omega], wp]];

  (* Decide which implementation to use *)
  Switch[OptionValue[Method],
    Automatic,
      TRF = TeukolskyRadialMST,
    "MST" | {"MST", OptionsPattern[TeukolskyRadialMST]},
      TRF = TeukolskyRadialMST,
    "NumericalIntegration" | {"NumericalIntegration", OptionsPattern[TeukolskyRadialNumericalIntegration]},
      TRF = TeukolskyRadialNumericalIntegration;,
    "HeunC",
      (* Some options are not supported for the HeunC method modes *)
      Do[
        If[OptionValue[opt] =!= Automatic, Message[TeukolskyRadial::hcopt, opt]];,
        {opt, {WorkingPrecision, PrecisionGoal, AccuracyGoal}}
      ];
      TRF = TeukolskyRadialHeunC;,
    _,
      Message[TeukolskyRadial::optx, Method -> OptionValue[Method]];
      Return[$Failed];
  ];

  (* Check only supported sub-options have been specified *)  
  If[subopts =!= (subopts = FilterRules[subopts, Options[TRF]]),
    Message[TeukolskyRadial::optx, Method -> OptionValue[Method]];
  ];

  (* Call the chosen implementation *)
  TRF[s, l, m, a, \[Omega], BCs, {wp, prec, acc}, Sequence@@subopts]
];


(* ::Section::Closed:: *)
(*TeukolskyRadialFunction*)


(* ::Subsection::Closed:: *)
(*Output format*)


(* ::Subsubsection::Closed:: *)
(*Icons*)


icons = <|
 "In" -> Graphics[{
         Line[{{0,1/2},{1/2,1},{1,1/2},{1/2,0},{0,1/2}}],
         Line[{{3/4,1/4},{1/2,1/2}}],
         {Arrowheads[0.2],Arrow[Line[{{1/2,1/2},{1/4,3/4}}]]},
         {Arrowheads[0.2],Arrow[Line[{{1/2,1/2},{3/4,3/4}}]]}},
         Background -> White,
         ImageSize -> Dynamic[{Automatic, 3.5 CurrentValue["FontCapHeight"]/AbsoluteCurrentValue[Magnification]}]],
 "Up" -> Graphics[{
         Line[{{0,1/2},{1/2,1},{1,1/2},{1/2,0},{0,1/2}}],
         Line[{{1/4,1/4},{1/2,1/2}}],
         {Arrowheads[0.2],Arrow[Line[{{1/2,1/2},{1/4,3/4}}]]},
         {Arrowheads[0.2],Arrow[Line[{{1/2,1/2},{3/4,3/4}}]]}},
         Background -> White,
         ImageSize -> Dynamic[{Automatic, 3.5 CurrentValue["FontCapHeight"]/AbsoluteCurrentValue[Magnification]}]]
|>;


(* ::Subsubsection::Closed:: *)
(*Formatting of TeukolskyRadialFunction*)


TeukolskyRadialFunction /:
 MakeBoxes[trf:TeukolskyRadialFunction[s_, l_, m_, a_, \[Omega]_, assoc_], form:(StandardForm|TraditionalForm)] :=
 Module[{summary, extended},
  summary = {Row[{BoxForm`SummaryItem[{"s: ", s}], "  ",
                  BoxForm`SummaryItem[{"l: ", l}], "  ",
                  BoxForm`SummaryItem[{"m: ", m}], "  ",
                  BoxForm`SummaryItem[{"a: ", a}], "  ",
                  BoxForm`SummaryItem[{"\[Omega]: ", \[Omega]}]}],
             BoxForm`SummaryItem[{"Domain: ", assoc["Domain"]}],
             BoxForm`SummaryItem[{"Boundary Conditions: " , assoc["BoundaryConditions"]}]};
  extended = {BoxForm`SummaryItem[{"Eigenvalue: ", assoc["Eigenvalue"]}],
              BoxForm`SummaryItem[{"Transmission Amplitude: ", assoc["Amplitudes", "Transmission"]}],
              BoxForm`SummaryItem[{"Incidence Amplitude: ", Lookup[assoc["Amplitudes"], "Incidence", Missing]}],
              BoxForm`SummaryItem[{"Reflection Amplitude: ", Lookup[assoc["Amplitudes"], "Reflection", Missing]}],
              BoxForm`SummaryItem[{"Method: ", First[assoc["Method"]]}],
              BoxForm`SummaryItem[{"Method options: ",Column[Rest[assoc["Method"]]]}]};
  BoxForm`ArrangeSummaryBox[
    TeukolskyRadialFunction,
    trf,
    Lookup[icons, assoc["BoundaryConditions"], None],
    summary,
    extended,
    form]
];


(* ::Subsection::Closed:: *)
(*Accessing attributes*)


TeukolskyRadialFunction[s_, l_, m_, a_, \[Omega]_, assoc_][y_String] /; !MemberQ[{"RadialFunction"}, y] :=
  assoc[y];


(* ::Subsection::Closed:: *)
(*Numerical evaluation*)


SetAttributes[TeukolskyRadialFunction, {NumericFunction}];


outsideDomainQ[r_, rmin_, rmax_] := Min[r]<rmin || Max[r]>rmax;


TeukolskyRadialFunction[s_, l_, m_, a_, \[Omega]_, assoc_][r:(_?NumericQ|{_?NumericQ..})] :=
 Module[{rmin, rmax},
  {rmin, rmax} = assoc["Domain"];
  If[outsideDomainQ[r, rmin, rmax],
    Message[TeukolskyRadialFunction::dmval, #]& /@ Select[Flatten[{r}], outsideDomainQ[#, rmin, rmax]&];
  ];
  Quiet[assoc["RadialFunction"][r], InterpolatingFunction::dmval]
 ];


Derivative[n_][TeukolskyRadialFunction[s_, l_, m_, a_, \[Omega]_, assoc_]][r:(_?NumericQ|{_?NumericQ..})] :=
 Module[{rmin, rmax},
  {rmin, rmax} = assoc["Domain"];
  If[outsideDomainQ[r, rmin, rmax],
    Message[TeukolskyRadialFunction::dmval, #]& /@ Select[Flatten[{r}], outsideDomainQ[#, rmin, rmax]&];
  ];
  Quiet[Derivative[n][assoc["RadialFunction"]][r], InterpolatingFunction::dmval]
 ];


(* ::Section::Closed:: *)
(*End Package*)


(* ::Subsection::Closed:: *)
(*Protect symbols*)


SetAttributes[{TeukolskyRadial, TeukolskyRadialFunction}, {Protected, ReadProtected}];


(* ::Subsection::Closed:: *)
(*End*)


End[];
EndPackage[];
