(* ::Package:: *)

(* ::Chapter:: *)
(*TeukolskySource*)


BeginPackage["Teukolsky`TeukolskySource`"];

TeukolskySourceObject::usage = "TeukolskySourceObject[assoc] an object which contains a Teukolsky source"

TeukolskyPointParticleSource::usage = "TeukolskyPointParticleSource[s, orbit] Point particle source for the Teukolsky equation"

Begin["`Private`"];


TeukolskyPointParticleSource[s_, orbit_]:=Module[{assoc},

    If[(s==-2 || s==0) && orbit["e"] == 0 && Abs[orbit["Inclination"]] == 1,
         Return[TeukolskyPointParticleSourceCircular[s,orbit]],
         Print["No point-particle source yet available for those parameters"];
      ]
	
]


TeukolskyPointParticleSourceCircular[0, orbit_] := Module[{assoc, \[Alpha], gtt, gt\[Phi], \[CapitalDelta], ut, a, r0, E0, Lz, \[CapitalSigma], \[Rho]},
	a = orbit["a"];
	r0 = orbit["p"];
	
	
	(*BEGIN FIXME: this should come from KerrGeoOrbit but currently we don't have code to compute the four-velocity in the KerrGeodesics package*)
	E0 = orbit["Energy"];
	Lz = orbit["AngularMomentum"];
	
	\[CapitalSigma] = r0^2;
	\[Rho] = -1/r0;

	\[CapitalDelta] = r0^2 - 2r0 + a^2;
	gtt = -(1/\[CapitalSigma])((r0^2+a^2)^2/\[CapitalDelta]-a^2);
	gt\[Phi] = -((2 r0 a)/(\[CapitalSigma] \[CapitalDelta]));
	ut = gt\[Phi] Lz - gtt E0;
	(*END FIXME*)

	\[Alpha] = -((4 \[Pi] r0)/(ut \[CapitalDelta]));

	assoc = <|  "s" -> 0,
				"SourceType" -> "PointParticle",
				"Orbit" -> orbit,
				"\[Alpha]" -> \[Alpha]
			|>;

	TeukolskySourceObject[assoc]

]


TeukolskyPointParticleSourceCircular[-2, orbit_] := Module[{assoc, a,r0, E0, Lz, \[CapitalSigma], Cnn, Cnmb, Cmbmb, ut, \[Rho], \[CapitalDelta], gtt, gt\[Phi]},
	a = orbit["a"];
	r0 = orbit["p"];
	
	(*BEGIN FIXME: this should come from KerrGeoOrbit but currently we don't have code to compute the four-velocity in the KerrGeodesics package*)
	
	E0 = orbit["Energy"];
	Lz = orbit["AngularMomentum"];
	
	\[CapitalSigma] = r0^2;
	\[Rho] = -1/r0;
	
	\[CapitalDelta] = r0^2 - 2r0 + a^2;
	gtt = -(1/\[CapitalSigma])((r0^2+a^2)^2/\[CapitalDelta]-a^2);
	gt\[Phi] = -((2 r0 a)/(\[CapitalSigma] \[CapitalDelta]));
	ut = gt\[Phi] Lz - gtt E0;
	(*END FIXME*)
	
	Cnn  = 1/(4 \[CapitalSigma]^3 ut) (E0(r0^2+a^2)-a Lz)^2;
	Cnmb = \[Rho]/(2Sqrt[2]\[CapitalSigma]^2 ut) (E0(r0^2+a^2)-a Lz)(I (a E0 - Lz));
	Cmbmb = \[Rho]^2/(2 \[CapitalSigma] ut) (I (a E0 - Lz))^2;

	assoc = <| "s" -> -2, 
			   "SourceType" -> "PointParticle",
			   "Orbit" -> orbit,
			   "Cnn" -> Cnn,
			   "Cnmb" -> Cnmb,
			   "Cmbmb" -> Cmbmb
			   |>;
	
	
	
	TeukolskySourceObject[assoc]

]


Format[TeukolskySourceObject[assoc_]] := "TeukolskySourceObject[<<>>]";

TeukolskySourceObject[assoc_][string_] := assoc[string]


End[];
EndPackage[];
