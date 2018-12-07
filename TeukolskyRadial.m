(* ::Package:: *)

(* ::Chapter:: *)
(*TeukolskyRadial: calculate the radial solutions to the Teukolsky equation*)


BeginPackage["Teukolsky`TeukolskyRadial`",
  {"SpinWeightedSpheroidalHarmonics`",
   "Teukolsky`RenormalizedAngularMomentum`",
   "Teukolsky`MST`"}
];

TeukolskyRadial::usage = "TeukolskyRadial[s,l,m,a,\[Omega]] computes solutions to the radial Teukolsky equation."
TeukolskyRadialFunction::usage = "TeukolskyRadialFunction[s, l, m, a, \[Omega], assoc] an object representing solutions to the Teukolsky equation."

Begin["`Private`"];


Options[TeukolskyRadial] = {Method -> {"MST", "RenormalizedAngularMomentum" -> "Monodromy"}, "BoundaryConditions" -> {"In", "Up"}};

TeukolskyRadial[s_Integer, l_Integer, m_Integer, a_, \[Omega]_, OptionsPattern[]] := Module[{assoc, \[Lambda], \[Nu], rmin, rmax, solFuncs, method},
	\[Lambda] = SpinWeightedSpheroidalEigenvalue[s,l,m,a \[Omega]];
	
	Switch[OptionValue[Method],
		"MST"|{"MST", ___},
			\[Nu] = RenormalizedAngularMomentum[s, l, m, a, 2 \[Omega], \[Lambda], Method->"RenormalizedAngularMomentum"/.OptionValue[Method][[2;;]] ];
			method = {"MST", "RenormalizedAngularMomentum" -> \[Nu]};
			solFuncs = OptionValue["BoundaryConditions"] /. {"In" -> Teukolsky`MST`Private`MSTRadialIn[s,l,m,a,2\[Omega],\[Lambda],\[Nu]], "Up" -> Teukolsky`MST`Private`MSTRadialUp[s,l,m,a,2\[Omega],\[Lambda],\[Nu]]};,
		{"SasakiNakamura","rmin" -> _, "rmax" -> _},
			method = {"SasakiNakamura", rmin, rmax};
			solFuncs = $Failed;
	];
	
	
	
	assoc = Association[
		"Method" -> method,
		"BoundaryConditions" -> OptionValue["BoundaryConditions"],
		"Eigenvalue" -> \[Lambda],
		"SolutionFunctions" -> solFuncs
		];

	TeukolskyRadialFunction[s,l,m,a,\[Omega],assoc]

]


Format[TeukolskyRadialFunction[s_,l_,m_,a_,\[Omega]_,assoc_]] := "TeukolskyRadialFunction["<>ToString[s]<>","<>ToString[l]<>","<>ToString[m]<>","<>ToString[a]<>","<>ToString[m]<>","<>ToString[\[Omega]]<>",<<>>]";

TeukolskyRadialFunction[s_,l_,m_,a_,\[Omega]_,assoc_][y_String/;(y=="In"||y=="Up")] := Module[{assocNew=assoc},
	assocNew["SolutionFunctions"] = {Association[MapThread[#1 -> #2 &, {assoc["BoundaryConditions"], assoc["SolutionFunctions"]}]][y]};
	assocNew["BoundaryConditions"]={y};
	TeukolskyRadialFunction[s, l, m, a, \[Omega], assocNew]
];
TeukolskyRadialFunction[s_,l_,m_,a_,\[Omega]_,assoc_][y_String] /; !MemberQ[{"SolutionFunctions"},y]:= assoc[y];
TeukolskyRadialFunction[s_,l_,m_,a_,\[Omega]_,assoc_][r_?NumericQ] := Module[{},
	If[
		Length[assoc["BoundaryConditions"]]==1, 
		Return[assoc["SolutionFunctions"][[1]][r]],
		Return[Association[MapThread[#1 -> #2[r] &, {assoc["BoundaryConditions"], assoc["SolutionFunctions"]}]]]
	];	
]


End[];
EndPackage[];
