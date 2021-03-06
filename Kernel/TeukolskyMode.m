(* ::Package:: *)

(* ::Title:: *)
(*TeukolskyMode*)


(* ::Section::Closed:: *)
(*Create Package*)


(* ::Subsection::Closed:: *)
(*BeginPackage*)


BeginPackage["Teukolsky`TeukolskyMode`",
	{"Teukolsky`TeukolskySource`",
	 "Teukolsky`TeukolskyRadial`",
	 "Teukolsky`ConvolveSource`",
	 "KerrGeodesics`KerrGeoOrbit`",
	 "KerrGeodesics`OrbitalFrequencies`",
	 "SpinWeightedSpheroidalHarmonics`"}
];


(* ::Subsection::Closed:: *)
(*Unprotect symbols*)


ClearAttributes[{TeukolskyMode, TeukolskyPointParticleMode}, {Protected}];


(* ::Subsection::Closed:: *)
(*Usage messages*)


TeukolskyMode::usage = "TeukolskyMode[assoc] is an object which represents a Teukolsky mode.";
TeukolskyPointParticleMode::usage = "TeukolskyPointParticleMode[s, l, m, n, k, orbit] produces a "<>
 "TeukolskyMode representing a solution to the radial Teukolsky equation with a point particle source.";


(* ::Subsection::Closed:: *)
(*Error Messages*)


TeukolskyPointParticleMode::params = "Parameters s=`1`, e=`2`, x=`3` are not currently supported.";


(* ::Subsection::Closed:: *)
(*Begin Private section*)


Begin["`Private`"];


(* ::Section::Closed:: *)
(*TeukolskyPointParticleMode*)


SyntaxInformation[TeukolskyPointParticleMode] =
 {"ArgumentsPattern" -> {_, _, _, _, _, _, OptionsPattern[]}};


Options[TeukolskyPointParticleMode] = {};


TeukolskyPointParticleMode[s_Integer, l_Integer, m_Integer, n_Integer, k_Integer, orbit_KerrGeoOrbitFunction, opts:OptionsPattern[]] /; AllTrue[orbit["Frequencies"], InexactNumberQ] :=
 Module[{source, assoc, R, S, \[Omega], \[CapitalOmega]r, \[CapitalOmega]\[Phi], \[CapitalOmega]\[Theta], Z, a, \[Lambda]},
  If[(s != -2  && s != -1 && s != 0) || orbit["e"] != 0 || Abs[orbit["x"]] != 1,
    Message[TeukolskyPointParticleMode::params, s, orbit["e"], orbit["Inclination"]];
    Return[$Failed];
  ];

  (*{\[CapitalOmega]r, \[CapitalOmega]\[Theta], \[CapitalOmega]\[Phi]} = orbit["Frequencies"];*) (*This gives Mino frequencies, need BL frequencies*)
  {\[CapitalOmega]r, \[CapitalOmega]\[Theta], \[CapitalOmega]\[Phi]} = Values[KerrGeoFrequencies[orbit["a"], orbit["p"], orbit["e"], orbit["Inclination"]]];
  \[Omega] = m \[CapitalOmega]\[Phi] + n \[CapitalOmega]r + k \[CapitalOmega]\[Theta];
  a = orbit["a"];

  source = Teukolsky`TeukolskySource`Private`TeukolskyPointParticleSource[s, orbit];

  R = TeukolskyRadial[s, l, m, a, \[Omega]];  
  S = SpinWeightedSpheroidalHarmonicS[s, l, m, a \[Omega]];
  Z = Teukolsky`ConvolveSource`Private`ConvolveSource[R, S, source];

  assoc = <| "s" -> s, 
		     "l" -> l,
		     "m" -> m,
		     "a" -> a,
  		   "\[Omega]" -> \[Omega],
		     "Eigenvalue" -> R["In"]["Eigenvalue"],
 		    "Type" -> {"PointParticleCircular", "Radius" -> orbit["p"]},
		     "RadialFunctions" -> R,
		     "AngularFunction" -> S,
		     "Amplitudes" -> Z
		    |>;

  TeukolskyMode[assoc]
]


(* ::Section::Closed:: *)
(*TeukolskyMode*)


(* ::Subsection::Closed:: *)
(*Output format*)


TeukolskyMode /:
 MakeBoxes[tm:TeukolskyMode[assoc_], form:(StandardForm|TraditionalForm)] :=
 Module[{summary, extended},
  summary = {Row[{BoxForm`SummaryItem[{"s: ", assoc["s"]}], "  ",
                  BoxForm`SummaryItem[{"l: ", assoc["l"]}], "  ",
                  BoxForm`SummaryItem[{"m: ", assoc["m"]}], "  ",
                  BoxForm`SummaryItem[{"a: ", assoc["a"]}], "  ",
                  BoxForm`SummaryItem[{"\[Omega]: ", assoc["\[Omega]"]}]}],
             BoxForm`SummaryItem[{"Type: ", First[assoc["Type"]]}]};
  extended = {BoxForm`SummaryItem[{"Eigenvalue: ", assoc["Eigenvalue"]}],
              BoxForm`SummaryItem[{"Amplitude at \[ScriptCapitalI]: ", assoc["Amplitudes"]["\[ScriptCapitalI]"]}],
              BoxForm`SummaryItem[{"Amplitude at \[ScriptCapitalH]: ", assoc["Amplitudes"]["\[ScriptCapitalH]"]}],
              BoxForm`SummaryItem[{"Type details: ", Column[Rest[assoc["Type"]]]}]};
  BoxForm`ArrangeSummaryBox[
    TeukolskyMode,
    tm,
    None,
    summary,
    extended,
    form
  ]
];


(* ::Subsection::Closed:: *)
(*Accessing attributes*)


TeukolskyMode[assoc_]["EnergyFlux"] := EnergyFlux[TeukolskyMode[assoc]];


TeukolskyMode[assoc_]["Fluxes"] := <|"Energy" -> TeukolskyMode[assoc]["EnergyFlux"], "AngularMomentum" -> TeukolskyMode[assoc]["AngularMomentumFlux"]|>;


TeukolskyMode[assoc_]["AngularMomentumFlux"] := AngularMomentumFlux[TeukolskyMode[assoc]];


TeukolskyMode[assoc_][string_] := assoc[string];


(* ::Section::Closed:: *)
(*Fluxes*)


(* ::Subsection::Closed:: *)
(*Energy Flux*)


EnergyFlux[mode_TeukolskyMode] :=
 Module[{M = 1, s, l, m, a, \[Omega], \[Lambda], Z, rh, \[CapitalOmega]h, \[Kappa], \[Epsilon], AbsCSq, \[Alpha], p, FluxInf, FluxHor},
  a = mode["a"];
  s = mode["s"];
  l = mode["l"];
  m = mode["m"];
  \[Omega] = mode["\[Omega]"];
  \[Lambda] = mode["Eigenvalue"];
  Z = mode["Amplitudes"];

  rh = M + Sqrt[M^2-a^2];
  \[CapitalOmega]h = a/(2 M rh);
  \[Kappa] = \[Omega] - m \[CapitalOmega]h;
  \[Epsilon] = Sqrt[M^2-a^2]/(4 M rh);

  FluxInf = Abs[Z["\[ScriptCapitalI]"]]^2 \[Omega]^(2(1-Abs[s]))/(4 \[Pi]);
  FluxHor = Switch[s,
			-2,
			  AbsCSq = ((\[Lambda]+2)^2 + 4 a m \[Omega] - 4a^2 \[Omega]^2)(\[Lambda]^2+36 m a \[Omega] - 36 a^2 \[Omega]^2) + (2\[Lambda]+3)(96 a^2 \[Omega]^2 - 48 m a \[Omega]) + 144 \[Omega]^2 (M^2-a^2);
              \[Alpha] = (256(2M rh)^5 \[Kappa](\[Kappa]^2+4\[Epsilon]^2)(\[Kappa]^2+16\[Epsilon]^2)\[Omega]^3)/AbsCSq;
              \[Alpha] Abs[Z["\[ScriptCapitalH]"]]^2/(4 \[Pi] \[Omega]^2),
			-1,
			  p = \[Lambda]^2 + 4*a*\[Omega]*(m - a*\[Omega]);
			  \[Omega] Abs[Z["\[ScriptCapitalH]"]]^2 (2 M rh \[Kappa]) 4 ((2 M rh \[Kappa])^2+(M^2-a^2))/ (p \[Pi]),
			0,
			  (* The rh^2 factor vs arXiv:1003.1860 Eq. (55) is needed as \[Psi] = r R*)
			  1/(2 \[Pi] rh) \[Omega](\[Omega]-m \[CapitalOmega]h) Abs[Z["\[ScriptCapitalH]"]]^2*rh^2
			];

  <| "\[ScriptCapitalI]" -> FluxInf, "\[ScriptCapitalH]" -> FluxHor |>
];


(* ::Subsection::Closed:: *)
(*Angular Momentum Flux*)


AngularMomentumFlux[mode_TeukolskyMode] := EnergyFlux[mode] mode["m"]/mode["\[Omega]"];


(* ::Section::Closed:: *)
(*End Package*)


(* ::Subsection::Closed:: *)
(*Protect symbols*)


SetAttributes[{TeukolskyMode, TeukolskyPointParticleMode}, {Protected, ReadProtected}];


(* ::Subsection::Closed:: *)
(*End*)


End[];
EndPackage[];
