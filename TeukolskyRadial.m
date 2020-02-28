(* ::Package:: *)

(* ::Chapter:: *)
(*TeukolskyRadial: calculate the radial solutions to the Teukolsky equation*)


BeginPackage["Teukolsky`TeukolskyRadial`",
  {"SpinWeightedSpheroidalHarmonics`",
   "Teukolsky`MST`RenormalizedAngularMomentum`",
   "Teukolsky`MST`MST`"}
];

TeukolskyRadial::usage = "TeukolskyRadial[s,l,m,a,\[Omega]] computes solutions to the radial Teukolsky equation."
TeukolskyRadialFunction::usage = "TeukolskyRadialFunction[s, l, m, a, \[Omega], assoc] an object representing solutions to the Teukolsky equation."

Begin["`Private`"];


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


Format[TeukolskyRadialFunction[s_,l_,m_,a_,\[Omega]_,assoc_]] := "TeukolskyRadialFunction["<>ToString[s]<>","<>ToString[l]<>","<>ToString[m]<>","<>ToString[a]<>","<>ToString[\[Omega]]<>",<<>>]";

TeukolskyRadialFunction[s_,l_,m_,a_,\[Omega]_,assoc_][y:("In"|"Up")] := Module[{assocNew=assoc},
	assocNew["SolutionFunctions"] = First[Pick[assoc["SolutionFunctions"], assoc["BoundaryConditions"], y]];
	assocNew["Amplitudes"] = First[Pick[assoc["Amplitudes"], assoc["BoundaryConditions"], y]];
	assocNew["BoundaryConditions"] = y;
	TeukolskyRadialFunction[s, l, m, a, \[Omega], assocNew]
];
TeukolskyRadialFunction[s_,l_,m_,a_,\[Omega]_,assoc_][y_String] /; !MemberQ[{"SolutionFunctions"},y]:= assoc[y];
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


End[];
EndPackage[];
