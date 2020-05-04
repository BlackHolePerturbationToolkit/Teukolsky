(* ::Package:: *)

(* ::Title:: *)
(*TeukolskyRadial*)


(* ::Section::Closed:: *)
(*Create Package*)


(* ::Subsection::Closed:: *)
(*BeginPackage*)


BeginPackage["Teukolsky`TeukolskyRadial`",
  {"SpinWeightedSpheroidalHarmonics`",
   "Teukolsky`MST`RenormalizedAngularMomentum`",
   "Teukolsky`MST`MST`"}
];


(* ::Subsection::Closed:: *)
(*Usage messages*)


TeukolskyRadial::usage = "TeukolskyRadial[s,l,m,a,\[Omega]] computes solutions to the radial Teukolsky equation."
TeukolskyRadialFunction::usage = "TeukolskyRadialFunction[s, l, m, a, \[Omega], assoc] an object representing solutions to the Teukolsky equation."


(* ::Subsection:: *)
(*Error Messages*)


(* ::Subsection::Closed:: *)
(*Begin Private section*)


Begin["`Private`"];


(* ::Section::Closed:: *)
(*TeukolskyRadial*)


(* ::Subsection::Closed:: *)
(*TeukolskyRadial*)


Options[TeukolskyRadial] = {Method -> {"MST", "RenormalizedAngularMomentum" -> "Monodromy"}, "BoundaryConditions" -> {"In", "Up"}};

TeukolskyRadial[s_Integer, l_Integer, m_Integer, a_, \[Omega]_, OptionsPattern[]] := Module[{assoc, \[Lambda], \[Nu], rmin, rmax, solFuncs, method, norms},
	\[Lambda] = SpinWeightedSpheroidalEigenvalue[s,l,m,a \[Omega]];
	
	Switch[OptionValue[Method],
		"MST"|{"MST", ___},
			\[Nu] = RenormalizedAngularMomentum[s, l, m, a, \[Omega], \[Lambda], Method->"RenormalizedAngularMomentum"/.OptionValue[Method][[2;;]] ];
			method = {"MST", "RenormalizedAngularMomentum" -> \[Nu]};
			norms = Teukolsky`MST`MST`Private`Amplitudes[s,l,m,a,2\[Omega],\[Nu],\[Lambda]];
			solFuncs = OptionValue["BoundaryConditions"] /.
			  {"In" -> Teukolsky`MST`MST`Private`MSTRadialIn[s,l,m,a,2\[Omega],\[Nu],\[Lambda],norms["In"]["Transmission"]],
			   "Up" -> Teukolsky`MST`MST`Private`MSTRadialUp[s,l,m,a,2\[Omega],\[Nu],\[Lambda],norms["Up"]["Transmission"]]};
			norms = norms/norms[[All, "Transmission"]];,
		{"SasakiNakamura", "rmin" -> _, "rmax" -> _},
			method = {"SasakiNakamura", rmin, rmax} /. OptionValue[Method];
			solFuncs = $Failed;
	];
	
	
	
	assoc = Association[
		"s" -> s,
		"l" -> l,
		"m" -> m,
		"\[Omega]" -> \[Omega],
		"Method" -> method,
		"BoundaryConditions" -> OptionValue["BoundaryConditions"],
		"Eigenvalue" -> \[Lambda],
		"Amplitudes" -> norms,
		"SolutionFunctions" -> solFuncs
		];

	TeukolskyRadialFunction[s,l,m,a,\[Omega],assoc]

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
    form,
    "Interpretable" -> Automatic]
];


(* ::Subsection::Closed:: *)
(*Accessing attributes*)


TeukolskyRadialFunction[s_,l_,m_,a_,\[Omega]_,assoc_][y:("In"|"Up")] := Module[{assocNew=assoc},
	assocNew["SolutionFunctions"] = First[Pick[assoc["SolutionFunctions"], assoc["BoundaryConditions"], y]];
	assocNew["Amplitudes"] = First[Pick[assoc["Amplitudes"], assoc["BoundaryConditions"], y]];
	assocNew["BoundaryConditions"] = y;
	TeukolskyRadialFunction[s, l, m, a, \[Omega], assocNew]
];

TeukolskyRadialFunction[s_,l_,m_,a_,\[Omega]_,assoc_][y_String] /; !MemberQ[{"SolutionFunctions"},y]:= assoc[y];


(* ::Subsection::Closed:: *)
(*Numerical evaluation*)


TeukolskyRadialFunction[s_,l_,m_,a_,\[Omega]_,assoc_][r_?NumericQ] := Module[{},
	If[
		Head[assoc["BoundaryConditions"]] === List,
		Return[Association[MapThread[#1 -> #2[r]/assoc["Amplitudes"][#1]["Transmission"] &, {assoc["BoundaryConditions"], assoc["SolutionFunctions"]}]]], 
		Return[assoc["SolutionFunctions"][r]/assoc["Amplitudes"]["Transmission"]]
	];	
];

Derivative[n_][TeukolskyRadialFunction[s_,l_,m_,a_,\[Omega]_,assoc_]][r_?NumericQ] := Module[{},
	If[
		Head[assoc["BoundaryConditions"]] === List,
		Return[Association[MapThread[#1 -> Derivative[n][#2][r]/assoc["Amplitudes"][#1]["Transmission"] &, {assoc["BoundaryConditions"], assoc["SolutionFunctions"]}]]], 
		Return[Derivative[n][assoc["SolutionFunctions"]][r]/assoc["Amplitudes"]["Transmission"]]
	];	
];


(* ::Section::Closed:: *)
(*End Package*)


End[];
EndPackage[];
