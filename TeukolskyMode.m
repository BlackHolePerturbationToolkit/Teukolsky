(* ::Package:: *)

(* ::Chapter:: *)
(*TeukolskyMode*)


BeginPackage["Teukolsky`TeukolskyMode`",
	{"Teukolsky`TeukolskySource`",
	 "Teukolsky`TeukolskyRadial`",
	 "Teukolsky`ConvolveSource`",
	 "KerrGeodesics`KerrGeoOrbit`",
	 "KerrGeodesics`OrbitalFrequencies`",
	 "SpinWeightedSpheroidalHarmonics`"}
];

TeukolskyModeObject::usage = "TeukolskyModeObject[assoc] an object which contains a Teukolsky mode"

TeukolskyPointParticleMode::usage = "TeukolskyPointParticleMode[s, l, m, n, k, orbit] Solve the Teukolsky equation with a point particle source"

Begin["`Private`"];


TeukolskyPointParticleMode[s_Integer, l_Integer, m_Integer, n_Integer, k_Integer, orbit_KerrGeoOrbitFunction] := Module[{source},
	source = TeukolskyPointParticleSource[s, orbit];
	TeukolskyPointParticleMode[s,l,m,n,k,orbit,source]
]

TeukolskyPointParticleMode[s_Integer, l_Integer, m_Integer, n_Integer, k_Integer, orbit_KerrGeoOrbitFunction, source_TeukolskySourceObject]:=Module[{assoc, R, S, \[Omega], \[CapitalOmega]r, \[CapitalOmega]\[Phi], \[CapitalOmega]\[Theta], Z, a, Fluxes},

(*{\[CapitalOmega]r, \[CapitalOmega]\[Theta], \[CapitalOmega]\[Phi]} = orbit["Frequencies"];*) (*This gives Mino frequencies, need BL frequencies*)

{\[CapitalOmega]r, \[CapitalOmega]\[Theta], \[CapitalOmega]\[Phi]} = Values[KerrGeoFrequencies[orbit["a"], orbit["p"], orbit["e"], orbit["Inclination"]]];

\[Omega] = m \[CapitalOmega]\[Phi] + n \[CapitalOmega]r + k \[CapitalOmega]\[Theta];

a = orbit[[1]]; (*FIXME: should be able to get parameters from the orbit via assocication*)

R = TeukolskyRadial[s, l, m, a, \[Omega]];  

S = SpinWeightedSpheroidalHarmonicS[s, l, m, a \[Omega]];

Z = Teukolsky`ConvolveSource`Private`ConvolveSource[R, S, source];

(*For flux calculation we need to check the normalization*)
Fluxes = <|"FluxInf" -> Abs[Z["ZInf"]]^2/(4 \[Pi] \[Omega]^2)|>;

assoc = <| "s" -> s, 
		   "l" -> l,
		   "m" -> m,
		   "n" -> n,
		   "k" -> k,
		   "\[Omega]" -> \[Omega],
		    "Type" -> "PointParticle", 
		    "Radial" -> R,
		    "Angular" -> S,
		    "Amplitudes" -> Z,
		    "Fluxes" -> Fluxes
		    |>;

TeukolskyModeObject[assoc]
]


Format[TeukolskyModeObject[assoc_]] := "TeukolskyModeObject["<>ToString[assoc["s"]]<>","<>ToString[assoc["l"]]<>","<>ToString[assoc["m"]]<>","<>ToString[assoc["n"]]<>","<>ToString[assoc["k"]]<>",<<>>]";

TeukolskyModeObject[assoc_][string_] := assoc[string]


End[];
EndPackage[];
