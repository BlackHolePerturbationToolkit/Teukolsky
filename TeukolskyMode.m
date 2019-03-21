(* ::Package:: *)

(* ::Chapter:: *)
(*TeukolskyMode*)


BeginPackage["Teukolsky`TeukolskyMode`",
	{"Teukolsky`TeukolskySource`",
	 "Teukolsky`TeukolskyRadial`",
	 "KerrGeodesics`"}
];

TeukolskyModeObject::usage = "TeukolskyModeObject[assoc] an object which contains a Teukolsky mode"

TeukolskyPointParticleMode::usage = "TeukolskyPointParticleMode[s, l, m, n, k, orbit] Solve the Teukolsky equation with a point particle source"

Begin["`Private`"];


TeukolskyPointParticleMode[s_, l_, m_, n_, k_, orbit_]:=Module[{assoc, source, radial, \[Omega], \[CapitalOmega]r, \[CapitalOmega]\[Phi], \[CapitalOmega]\[Theta]},

{\[CapitalOmega]r, \[CapitalOmega]\[Theta], \[CapitalOmega]\[Phi]} = orbit["Frequencies"];
\[Omega] = m \[CapitalOmega]\[Phi] + n \[CapitalOmega]r + k \[CapitalOmega]\[Theta];
radial = TeukolskyRadial[s, l, m, orbit[[1]], \[Omega]];  (*FIXME: should be able to get parameters from the orbit via assocication*)

source = TeukolskyPointParticleSource[s, orbit];

assoc = <| "s" -> s, 
		    "modeType" -> "PointParticle", 
		    "source" -> source, 
		    "radial" -> radial,
		    "orbit" -> orbit
		    |>;

TeukolskyModeObject[assoc]
]


Format[TeukolskyModeObject[assoc_]] := "TeukolskyModeObject[<<>>]";

TeukolskyModeObject[assoc_][string_] := assoc[string]


End[];
EndPackage[];
