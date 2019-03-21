(* ::Package:: *)

(* ::Chapter:: *)
(*TeukolskyMode*)


BeginPackage["Teukolsky`TeukolskyMode`",
	{"Teukolsky`TeukolskySource`",
	 "Teukolsky`TeukolskyRadial`",
	 "Teukolsky`ConvolveSource`",
	 "KerrGeodesics`",
	 "SpinWeightedSpheroidalHarmonics`"}
];

TeukolskyModeObject::usage = "TeukolskyModeObject[assoc] an object which contains a Teukolsky mode"

TeukolskyPointParticleMode::usage = "TeukolskyPointParticleMode[s, l, m, n, k, orbit] Solve the Teukolsky equation with a point particle source"

Begin["`Private`"];


TeukolskyPointParticleMode[s_Integer, l_Integer, m_Integer, n_Integer, k_Integer, orbit_KerrGeoOrbitFunction] := Module[{source},
	source = TeukolskyPointParticleSource[s, orbit];
	TeukolskyPointParticleMode[s,l,m,n,k,orbit,source]
]

TeukolskyPointParticleMode[s_Integer, l_Integer, m_Integer, n_Integer, k_Integer, orbit_KerrGeoOrbitFunction, source_TeukolskySourceObject]:=Module[{assoc, R, S, \[Omega], \[CapitalOmega]r, \[CapitalOmega]\[Phi], \[CapitalOmega]\[Theta], Z, a},

{\[CapitalOmega]r, \[CapitalOmega]\[Theta], \[CapitalOmega]\[Phi]} = orbit["Frequencies"];
\[Omega] = m \[CapitalOmega]\[Phi] + n \[CapitalOmega]r + k \[CapitalOmega]\[Theta];

a = orbit[[1]]; (*FIXME: should be able to get parameters from the orbit via assocication*)

R = TeukolskyRadial[s, l, m, a, \[Omega]];  

S = SpinWeightedSpheroidalHarmonicS[s, l, m, a \[Omega]];

Z = Teukolsky`ConvolveSource`Private`ConvolveSource[R, S, source];

assoc = <| "s" -> s, 
		   "l" -> l,
		   "m" -> m,
		   "n" -> n,
		   "k" -> k,
		    "Type" -> "PointParticle", 
		    "Radial" -> R,
		    "Angular" -> S,
		    "Amplitudes" -> Z
		    |>;

TeukolskyModeObject[assoc]
]


Format[TeukolskyModeObject[assoc_]] := "TeukolskyModeObject["<>ToString[assoc["s"]]<>","<>ToString[assoc["l"]]<>","<>ToString[assoc["m"]]<>","<>ToString[assoc["n"]]<>","<>ToString[assoc["k"]]<>",<<>>]";

TeukolskyModeObject[assoc_][string_] := assoc[string]


End[];
EndPackage[];
