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
	 "KerrGeodesics`",
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


TeukolskyPointParticleMode::mode = "Mode with n=`1`, k=`2` not defined for `3` orbit.";


TeukolskyPointParticleMode::noret = "This package does not compute the retarded solution `1`. If you want the extended homogeneous solutions at this radius use [\"ExtendedHomogeneous\" -> \"\[ScriptCapitalH]|\[ScriptCapitalI]\"]`2`[`3`] instead.";


TeukolskyPointParticleMode::type = "Only bound orbits supported but orbit type is: `1`.";


TeukolskyPointParticleMode::mino = "Orbit parametrization \"`1`\" is not currently supported. Use \"Mino\" parametrization instead.";


(* ::Subsection::Closed:: *)
(*Begin Private section*)


Begin["`Private`"];


(* ::Section::Closed:: *)
(*TeukolskyPointParticleMode*)


SyntaxInformation[TeukolskyPointParticleMode] =
 {"ArgumentsPattern" -> {_, _, _, _., _., _, OptionsPattern[]}};


Options[TeukolskyPointParticleMode] = {"Domain" -> Automatic};


TeukolskyPointParticleMode[s_Integer, l_Integer, m_Integer, n_Integer, k_Integer, orbit_KerrGeoOrbitFunction, opts:OptionsPattern[]] /; AllTrue[orbit["Frequencies"], InexactNumberQ] :=
 Module[{source, assoc, domain, Ruser, R, S, \[Omega], \[CapitalOmega]r, \[CapitalOmega]\[Phi], \[CapitalOmega]\[Theta], Z, \[Lambda], rmin, rmax, a, p, e, x},
  {a, p, e, x} = orbit /@ {"a", "p", "e", "Inclination"};  

  If[orbit["Parametrization"] =!= "Mino",
    Message[TeukolskyPointParticleMode::mino, orbit["Parametrization"]];
    Return[$Failed]];
  If[!MatchQ[orbit["Type"], {"Bound", ___} | {"Unbound" | "MarginallyBound", "Circular", __}],
    Message[TeukolskyPointParticleMode::type, orbit["Type"]];
    Return[$Failed]];
  If[{e, Abs[x]} == {0, 1} && (n != 0 || k != 0),
    Message[TeukolskyPointParticleMode::mode, n, k, "circular"];
    Return[$Failed]];
  If[e == 0 && n != 0,
    Message[TeukolskyPointParticleMode::mode, n, k, "spherical"];
    Return[$Failed]];
  If[Abs[x] == 1 && k != 0,
    Message[TeukolskyPointParticleMode::mode, n, k, "eccentric"];
    Return[$Failed]];
  If[!MemberQ[{-2,-1,0,1,2},s], 
    Message[TeukolskyPointParticleMode::params,s];
    Return[$Failed];
    ];

  (*{\[CapitalOmega]r, \[CapitalOmega]\[Theta], \[CapitalOmega]\[Phi]} = orbit["Frequencies"];*) (*This gives Mino frequencies, need BL frequencies*)
  {\[CapitalOmega]r, \[CapitalOmega]\[Theta], \[CapitalOmega]\[Phi]} = {"\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(r\)]\)","\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(\[Theta]\)]\)","\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(\[Phi]\)]\)"}/.KerrGeoFrequencies[orbit["a"], orbit["p"], orbit["e"], orbit["Inclination"]];
  \[Omega] = m \[CapitalOmega]\[Phi] + n \[CapitalOmega]r + k \[CapitalOmega]\[Theta];

  source = Teukolsky`TeukolskySource`Private`TeukolskyPointParticleSource[s, orbit];

  domain = OptionValue["Domain"];
  If[MatchQ[domain, {_?NumericQ, _?NumericQ}],
    R = Ruser = TeukolskyRadial[s, l, m, a, \[Omega], Method ->
      {"NumericalIntegration", "Domain"-> {"In" -> domain, "Up" -> domain}}];
  ,
    R = Ruser = TeukolskyRadial[s, l, m, a, \[Omega]];
  ];

  If[e != 0,
    rmin = p/(1+e);
    rmax = p/(1-e);
    R = TeukolskyRadial[s, l, m, a, \[Omega], Method->{"NumericalIntegration","Domain"-> {"In"->{rmin,rmax}, "Up"->{rmin,rmax}}},
        "Amplitudes" -> <|"In"-> R["In"]["UnscaledAmplitudes"], "Up"-> R["Up"]["UnscaledAmplitudes"]|>,
        "RenormalizedAngularMomentum"-> R["In"]["RenormalizedAngularMomentum"], "Eigenvalue" -> R["In"]["Eigenvalue"]];
  ,
    rmin = rmax = p;
  ];
  
  S = SpinWeightedSpheroidalHarmonicS[s, l, m, a \[Omega]];
  Z = Teukolsky`ConvolveSource`Private`ConvolveSource[l, m, n, k, R, S, source];

  assoc = <| "s" -> s, 
		     "l" -> l,
		     "m" -> m,
		     "n" -> n,
		     "k" -> k,
		     "a" -> a,
  		   "\[Omega]" -> \[Omega],
		     "Eigenvalue" -> R["In"]["Eigenvalue"],
 		    "Type" ->
               Which[
               {e, Abs[x]} == {0, 1}, {"PointParticleCircular", "Radius" -> p},
               e == 0, {"PointParticleSpherical", "Radius" -> p, "Inclination" -> x},
               Abs[x] == 1, {"PointParticleEccentric", "Semi-latus Rectum" -> p, "Eccentricity" -> e},
               True, {"PointParticleGeneric", "Semi-latus Rectum" -> p, "Eccentricity" -> e , "Inclination" -> x}],
 		    "rmin" -> rmin,
 		    "rmax" -> rmax,
 		    "Domain" -> domain,
		     "RadialFunctions" -> Ruser,
		     "AngularFunction" -> S,
		     "Amplitudes" -> Z
		    |>;

  TeukolskyMode[assoc]
]


(* ::Subsection::Closed:: *)
(*Special cases*)


circularOrbitQ[orbit_KerrGeoOrbitFunction] := orbit["e"] == 0 && Abs[orbit["Inclination"]] == 1;


sphericalOrbitQ[orbit_KerrGeoOrbitFunction] := orbit["e"] == 0;


eccentricOrbitQ[orbit_KerrGeoOrbitFunction] := Abs[orbit["Inclination"]] == 1;


TeukolskyPointParticleMode[s_Integer, l_Integer, m_Integer, orbit_KerrGeoOrbitFunction, opts:OptionsPattern[]] /; circularOrbitQ[orbit] :=
  TeukolskyPointParticleMode[s, l, m, 0, 0, orbit, opts]


TeukolskyPointParticleMode[s_Integer, l_Integer, m_Integer, k_Integer, orbit_KerrGeoOrbitFunction, opts:OptionsPattern[]] /; sphericalOrbitQ[orbit] :=
  TeukolskyPointParticleMode[s, l, m, 0, k, orbit, opts]


TeukolskyPointParticleMode[s_Integer, l_Integer, m_Integer, n_Integer, orbit_KerrGeoOrbitFunction, opts:OptionsPattern[]] /; eccentricOrbitQ[orbit] :=
  TeukolskyPointParticleMode[s, l, m, n, 0, orbit, opts]


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


TeukolskyMode[assoc_][key_String] /; KeyExistsQ[assoc, key] := assoc[key];


Keys[m_TeukolskyMode] ^:= Join[Keys[m[[1]]], {"Fluxes", "EnergyFlux", "AngularMomentumFlux", "ExtendedHomogeneous" -> "\[ScriptCapitalH]", "ExtendedHomogeneous" -> "\[ScriptCapitalI]"}];


(* ::Subsection::Closed:: *)
(*Numerical Evaluation*)


TeukolskyMode[assoc_][r:(_?NumericQ|{_?NumericQ..})] :=
  Piecewise[{{assoc["Amplitudes"]["\[ScriptCapitalH]"]assoc["RadialFunctions"]["In"][r], r < assoc["rmin"]},
             {assoc["Amplitudes"]["\[ScriptCapitalI]"]assoc["RadialFunctions"]["Up"][r], r > assoc["rmax"]}},
            Message[TeukolskyPointParticleMode::noret, If[assoc["rmin"]==assoc["rmax"], "at the particle", "inside the libration region"], "", r]; Indeterminate
  ];


Derivative[n_][TeukolskyMode[assoc_]][r:(_?NumericQ|{_?NumericQ..})] :=
  Piecewise[{{assoc["Amplitudes"]["\[ScriptCapitalH]"]Derivative[n][assoc["RadialFunctions"]["In"]][r], r < assoc["rmin"]},
             {assoc["Amplitudes"]["\[ScriptCapitalI]"]Derivative[n][assoc["RadialFunctions"]["Up"]][r], r > assoc["rmax"]}},
            Message[TeukolskyPointParticleMode::noret, If[assoc["rmin"]==assoc["rmax"], "at the particle", "inside the libration region"], StringRepeat["'", n], r]; Indeterminate
  ];


(* ::Subsection::Closed:: *)
(*Numerical Evaluation using extended homogeneous solutions*)


TeukolskyMode[assoc_]["ExtendedHomogeneous" -> "\[ScriptCapitalH]"][r:(_?NumericQ|{_?NumericQ..})] :=
  assoc["Amplitudes"]["\[ScriptCapitalH]"]assoc["RadialFunctions"]["In"][r];


Derivative[n_][TeukolskyMode[assoc_]["ExtendedHomogeneous" -> "\[ScriptCapitalH]"]][r:(_?NumericQ|{_?NumericQ..})] :=
  assoc["Amplitudes"]["\[ScriptCapitalH]"]Derivative[n][assoc["RadialFunctions"]["In"]][r];


TeukolskyMode[assoc_]["ExtendedHomogeneous" -> "\[ScriptCapitalI]"][r:(_?NumericQ|{_?NumericQ..})] :=
  assoc["Amplitudes"]["\[ScriptCapitalI]"]assoc["RadialFunctions"]["Up"][r];


Derivative[n_][TeukolskyMode[assoc_]["ExtendedHomogeneous" -> "\[ScriptCapitalI]"]][r:(_?NumericQ|{_?NumericQ..})] :=
  assoc["Amplitudes"]["\[ScriptCapitalI]"]Derivative[n][assoc["RadialFunctions"]["Up"]][r];


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

  If[\[Omega] == 0, Return[<| "\[ScriptCapitalI]" -> 0, "\[ScriptCapitalH]" -> 0 |>]];

  rh = M + Sqrt[M^2-a^2];
  \[CapitalOmega]h = a/(2 M rh);
  \[Kappa] = \[Omega] - m \[CapitalOmega]h;
  \[Epsilon] = Sqrt[M^2-a^2]/(4 M rh);  	

  FluxInf = 
  Switch[s,
  -2, Abs[Z["\[ScriptCapitalI]"]]^2 \[Omega]^(2(1-Abs[s]))/(4 \[Pi]),
  -1, 2 Abs[Z["\[ScriptCapitalI]"]]^2 \[Omega]^(2(1-Abs[s]))/(4 \[Pi]),
  0, Abs[Z["\[ScriptCapitalI]"]]^2 \[Omega]^(2(1-Abs[s]))/(4 \[Pi]),
  1, AbsCSq = (\[Lambda]+2)^2+ 4 a \[Omega](m-a \[Omega]);
  2 (4 \[Omega]^4 Abs[Z["\[ScriptCapitalI]"]]^2)/(AbsCSq) \[Omega]^(2(1-Abs[s]))/(4 \[Pi]),
  2, AbsCSq = (4+\[Lambda])^2 (6+\[Lambda])^2+144 M^2 \[Omega]^2+8 a (4+\[Lambda]) (-4+5 (6+\[Lambda])) \[Omega] (m-a \[Omega])+48 a^2 \[Omega]^2 (2 (4+\[Lambda])+3 (m-a \[Omega])^2);
  (16 \[Omega]^8 Abs[Z["\[ScriptCapitalI]"]]^2)/(AbsCSq) \[Omega]^(2(1-Abs[s]))/(4 \[Pi])
  ];
                
  
  
  (*Abs[Z["\[ScriptCapitalI]"]]^2 \[Omega]^(2(1-Abs[s]))/(4 \[Pi]);*)
  FluxHor = Switch[s,
			-2,
			  AbsCSq = ((\[Lambda]+2)^2 + 4 a m \[Omega] - 4a^2 \[Omega]^2)(\[Lambda]^2+36 m a \[Omega] - 36 a^2 \[Omega]^2) + (2\[Lambda]+3)(96 a^2 \[Omega]^2 - 48 m a \[Omega]) + 144 \[Omega]^2 (M^2-a^2);
              \[Alpha] = (256(2M rh)^5 \[Kappa](\[Kappa]^2+4\[Epsilon]^2)(\[Kappa]^2+16\[Epsilon]^2)\[Omega]^3)/AbsCSq;
              \[Alpha] Abs[Z["\[ScriptCapitalH]"]]^2/(4 \[Pi] \[Omega]^2),
			-1,
			  p = \[Lambda]^2 + 4*a*\[Omega]*(m - a*\[Omega]);
			  2 \[Omega] Abs[Z["\[ScriptCapitalH]"]]^2 (2 M rh \[Kappa]) 4 ((2 M rh \[Kappa])^2+(M^2-a^2))/ (p \[Pi]),
			0,
			  (* The rh^2 factor vs arXiv:1003.1860 Eq. (55) is needed as \[Psi] = r R*)
			  1/(2 \[Pi] rh) \[Omega](\[Omega]-m \[CapitalOmega]h) Abs[Z["\[ScriptCapitalH]"]]^2*rh^2,
			 1,
			 2 (\[Omega] Abs[Z["\[ScriptCapitalH]"]]^2)/(32 \[Pi] \[Kappa] rh),
			 2,
			  (\[Omega] Abs[Z["\[ScriptCapitalH]"]]^2)/(512 \[Pi] rh^3 \[Kappa] (\[Kappa]^2+4 \[Epsilon]^2))
			];

  <| "\[ScriptCapitalI]" -> FluxInf, "\[ScriptCapitalH]" -> FluxHor |>
];


(* ::Subsection::Closed:: *)
(*Angular Momentum Flux*)


AngularMomentumFlux[mode_TeukolskyMode] := If[mode["\[Omega]"]!=0,EnergyFlux[mode] mode["m"]/mode["\[Omega]"], <| "\[ScriptCapitalI]" -> 0, "\[ScriptCapitalH]" -> 0 |>];


(* ::Section::Closed:: *)
(*End Package*)


(* ::Subsection::Closed:: *)
(*Protect symbols*)


SetAttributes[{TeukolskyMode, TeukolskyPointParticleMode}, {Protected, ReadProtected}];


(* ::Subsection::Closed:: *)
(*End*)


End[];
EndPackage[];
