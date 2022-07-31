(* ::Package:: *)

(* ::Title:: *)
(*TeukolskySource*)


(* ::Section::Closed:: *)
(*Create Package*)


(* ::Subsection::Closed:: *)
(*BeginPackage*)


BeginPackage["Teukolsky`TeukolskySource`", {"KerrGeodesics`", "KerrGeodesics`OrbitalFrequencies`"}];


(* ::Subsection:: *)
(*Usage messages*)


(* ::Subsection::Closed:: *)
(*Begin Private section*)


Begin["`Private`"];


(* ::Section::Closed:: *)
(*TeukolskyPointParticleSource*)


TeukolskyPointParticleSource[s_, orbit_] :=
 Module[{e, x},
  {e, x} = {orbit["e"], orbit["Inclination"]};

  Which[
  {e, Abs[x]} == {0, 1},
    Return[TeukolskyPointParticleSourceCircular[s,orbit]],
  e == 0,
    Return[TeukolskyPointParticleSourceSpherical[s,orbit]],
  Abs[x] == 1,
    Return[TeukolskyPointParticleSourceEccentric[s,orbit]];,
  True,
    Return[TeukolskyPointParticleSourceGeneric[s,orbit]];
  ];

  $Failed
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


TeukolskyPointParticleSourceSpherical[0, orbit_] := Module[{assoc, \[Alpha], gtt, gt\[Phi], \[CapitalDelta], ut, dtd\[Lambda], a, r0, e, x, \[ScriptCapitalE], \[ScriptCapitalL], \[CapitalUpsilon], \[CapitalSigma], \[Rho],\[CapitalOmega]r,\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi]},
	r0 =orbit["p"];
	a= orbit["a"];
	x = orbit["Inclination"];
	e = orbit["e"];
	{\[CapitalOmega]r,\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi]} = {"\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(r\)]\)", "\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(\[Theta]\)]\)", "\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(\[Phi]\)]\)"} /. KerrGeoFrequencies[a,r0,e,x];
	\[ScriptCapitalE]=orbit["Energy"];
	\[ScriptCapitalL]=orbit["AngularMomentum"];
	
	\[CapitalSigma] = Function[\[Theta],r0^2+a^2 Cos[\[Theta]]^2];
	\[CapitalDelta] = r0^2-2 r0 + a^2;	

  gtt = Function[\[Theta], -((a^4+2 r0^4+a^2 r0 (2 +3 r0)+a^2 (a^2+r0 (-2 +r0)) Cos[2 \[Theta]])/((a^2+r0 (-2 +r0)) (a^2+2 r0^2+a^2 Cos[2 \[Theta]])))];
  gt\[Phi] = Function[\[Theta], -((4 a  r0)/((a^2+r0 (-2 +r0)) (a^2+2 r0^2+a^2 Cos[2 \[Theta]])))];
  ut = Function[\[Theta],gt\[Phi][\[Theta]]\[ScriptCapitalL] - gtt[\[Theta]]\[ScriptCapitalE]];
  dtd\[Lambda] = Function[\[Theta],(\[ScriptCapitalE]((r0^2+a^2)^2/\[CapitalDelta] - a^2 Sin[\[Theta]]^2)+a \[ScriptCapitalL](1-(r0^2+a^2)/\[CapitalDelta]))];
  
  \[Alpha] = Function[\[Theta],-(8 \[CapitalOmega]\[Theta] r0 dtd\[Lambda][\[Theta]])/(\[CapitalDelta] ut[\[Theta]])];

  assoc = <|  "s" -> 0,
        "SourceType" -> "PointParticle",
        "Orbit" -> orbit,
        "\[Alpha]" -> \[Alpha]
      |>;

  TeukolskySourceObject[assoc]
]


TeukolskyPointParticleSourceEccentric[0, orbit_] := Module[{assoc, a, r0, e, x, \[ScriptCapitalE], \[ScriptCapitalL]},
	a= orbit["a"];
	
	assoc = <|  "s" -> 0,
        "SourceType" -> "PointParticle",
        "Orbit" -> orbit,
        "\[Alpha]" -> Function[{r,\[Theta]},-4 Pi(r^2+a^2 Cos[\[Theta]]^2)]
        |>;
        
     TeukolskySourceObject[assoc]
]


TeukolskyPointParticleSourceGeneric[0, orbit_] := Module[{assoc, a, r0, e, x, \[ScriptCapitalE], \[ScriptCapitalL]},
	a= orbit["a"];
	
	assoc = <|  "s" -> 0,
        "SourceType" -> "PointParticle",
        "Orbit" -> orbit,
        "\[Alpha]1" -> Function[r,r^2],
        "\[Alpha]2" -> Function[\[Theta],-4 Pi],
        "\[Alpha]3" -> Function[r,1],
        "\[Alpha]4" -> Function[\[Theta],-4 Pi a^2 Cos[\[Theta]]^2]
        |>;

  TeukolskySourceObject[assoc]
]


TeukolskyPointParticleSourceCircular[s:(-1|+1), orbit_] := Module[{assoc, a, r0, E0, \[CapitalOmega], Lz, \[CapitalSigma], c, S, B, Ar, Ati, ut, \[Rho], \[CapitalDelta], \[CapitalDelta]p, gtt, gt\[Phi]},
  a = orbit["a"];
  r0 = orbit["p"];

  (*BEGIN FIXME: this should come from KerrGeoOrbit but currently we don't have code to compute the four-velocity in the KerrGeodesics package*)

  E0 = orbit["Energy"];
  Lz = orbit["AngularMomentum"];

  \[CapitalSigma] = r0^2;
  \[Rho] = -1/r0;

  \[CapitalDelta] = r0^2 - 2r0 + a^2;
  \[CapitalDelta]p = 2r0 - 2;
  gtt = -(1/\[CapitalSigma])((r0^2+a^2)^2/\[CapitalDelta]-a^2);
  gt\[Phi] = -((2 r0 a)/(\[CapitalSigma] \[CapitalDelta]));
  ut = gt\[Phi] Lz - gtt E0;
  \[CapitalOmega] = 1/(Sqrt[r0^3]+a);
  (*END FIXME*)
  
  (* Define source terms: arXiv:2008.12703 Eq. (45). We have an extra factor of 1/2
     for s=-1 because we want \[Zeta]^2Subscript[\[Phi], 2], not 2\[Zeta]^2Subscript[\[Phi], 2], cf. Eq. (20a).  *)
  S = (4*\[Pi])/(Sqrt[2]r0)*Which[s==-1, 1/2, s==+1, 1];
  B = \[CapitalDelta]((r0^2+a^2)\[CapitalOmega] - a);
  Ar = r0(r0((r0^2+a^2)\[CapitalOmega]^2-1)+2 (1-a \[CapitalOmega])^2);
  Ati = r0 \[CapitalDelta] \[CapitalOmega];
  c = - \[CapitalDelta](1-a \[CapitalOmega]);

  assoc = <| "s" -> s, 
         "SourceType" -> "PointParticle",
         "Orbit" -> orbit,
         "\[ScriptCapitalS]" -> S,
         "B" -> B,
         "\!\(\*SuperscriptBox[\(A\), \((r)\)]\)" -> Ar,
         "\!\(\*SuperscriptBox[OverscriptBox[\(A\), \(~\)], \((i)\)]\)" -> Ati,
         "C" -> c 
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


TeukolskyPointParticleSourceSpherical[-2, orbit_]:=Module[
{a,r0,e,x,\[ScriptCapitalE],\[ScriptCapitalL],\[ScriptCapitalQ],\[Theta],\[CapitalDelta],\[CapitalSigma],\[Rho],\[Rho]b,ut,\[CapitalTheta],\[Lambda]0,CnnPlus,CnmPlus,CmmPlus,assoc},
(*Orbital Params*)
a   = orbit["a"];
r0  = orbit["p"];
e   = orbit["Eccentricity"];
x   = orbit["Inclination"];

(*Constants of motion & Frequencies*)
\[ScriptCapitalE] = orbit["Energy"];
\[ScriptCapitalL] = orbit["AngularMomentum"];
\[ScriptCapitalQ] = orbit["CarterConstant"];

(*Definitions and reparameterisation*)
\[Theta]=orbit["Trajectory"][[3]];
\[CapitalTheta] = Derivative[1][\[Theta]];

\[CapitalDelta]=r0^2-2r0+a^2;
\[CapitalSigma][\[Lambda]0_]:=r0^2+a^2 Cos[\[Theta][\[Lambda]0]]^2;
\[Rho][\[Lambda]0_]:=-1/(r0-I a Cos[\[Theta][\[Lambda]0]]);
\[Rho]b[\[Lambda]0_]:=-1/(r0+I a Cos[\[Theta][\[Lambda]0]]);

ut[\[Lambda]0_]:=1/\[CapitalSigma][\[Lambda]0] (\[ScriptCapitalE]((r0^2+a^2)^2/\[CapitalDelta]-a^2 Sin[\[Theta][\[Lambda]0]]^2)+a \[ScriptCapitalL](1-(r0^2+a^2)/\[CapitalDelta]));

(*Stress energy projected onto tetrad (\[PlusMinus]\[CapitalTheta] for travelling up or down)*)
CnnPlus=Function[{\[Lambda]0},Evaluate[1/(4 \[CapitalSigma][\[Lambda]0]^3 ut[\[Lambda]0]) (\[ScriptCapitalE](r0^2+a^2)-a \[ScriptCapitalL])^2]];
CnmPlus=Function[{\[Lambda]0},Evaluate[\[Rho][\[Lambda]0]/(2Sqrt[2] \[CapitalSigma][\[Lambda]0]^2 ut[\[Lambda]0]) (\[ScriptCapitalE](r0^2+a^2)-a \[ScriptCapitalL])(I Sin[\[Theta][\[Lambda]0]](a \[ScriptCapitalE]-\[ScriptCapitalL]/Sin[\[Theta][\[Lambda]0]]^2)+ \[CapitalTheta][\[Lambda]0])]];
CmmPlus=Function[{\[Lambda]0},Evaluate[\[Rho][\[Lambda]0]^2/(2\[CapitalSigma][\[Lambda]0] ut[\[Lambda]0]) (I Sin[\[Theta][\[Lambda]0]](a \[ScriptCapitalE]-\[ScriptCapitalL]/Sin[\[Theta][\[Lambda]0]]^2)+ \[CapitalTheta][\[Lambda]0])^2]];
(*Rewriting source*)
assoc = <|"s"->-2, "SourceType"->"PointParticle","Orbit"->orbit,"Cnn+"->CnnPlus,"Cnm+"->CnmPlus,"Cmm+"->CmmPlus|>;
Clear[a,r0,e,x,\[ScriptCapitalE],\[ScriptCapitalL],\[ScriptCapitalQ],\[Theta],\[CapitalDelta],\[CapitalSigma],\[Rho],\[Rho]b,ut,\[CapitalTheta],CnnPlus,CnmPlus,CmmPlus];
assoc
]


TeukolskyPointParticleSourceEccentric[-2, orbit_] := Module[{assoc, a, E0, Lz, Qc, Cab},
  a = orbit["a"];
  E0 = orbit["Energy"];
  Lz = orbit["AngularMomentum"];
  Qc = orbit["CarterConstant"];
  Cab[r_,\[Theta]_,sgnUr_,sgnU\[Theta]_] := Module[{\[CapitalSigma],\[Rho],ur,u\[Theta],rcomp,\[Theta]comp},
		\[CapitalSigma]=r^2+a^2Cos[\[Theta]]^2;
		\[Rho]=-1/(r-I a Cos[\[Theta]]);
		ur = Sqrt[-((r^2 - 2r + a^2)*(Qc + ((a*E0) - Lz)^2 + r^2)) + (a*Lz - E0*(a^2 + r^2))^2];
		u\[Theta] = Sqrt[Qc - a^2*(1 - E0^2)*Cos[\[Theta]]^2 - (Lz*Cos[\[Theta]]/Sin[\[Theta]])^2];
		rcomp=(E0(r^2+a^2) - a Lz + sgnUr ur)/(2\[CapitalSigma]);
		\[Theta]comp=\[Rho] (I Sin[\[Theta]](a E0 - Lz/Sin[\[Theta]]^2)+sgnU\[Theta] u\[Theta] )/Sqrt[2];
		{rcomp^2, rcomp \[Theta]comp, \[Theta]comp^2}
  ];

  assoc = <|  "s" -> -2,
        "SourceType" -> "PointParticle",
        "Orbit" -> orbit,
        "Cab" -> Cab
      |>;

  TeukolskySourceObject[assoc]
]


TeukolskyPointParticleSourceGeneric[-2, orbit_] := Module[{assoc, a, E0, Lz, Qc, Cab},
  a = orbit["a"];
  E0 = orbit["Energy"];
  Lz = orbit["AngularMomentum"];
  Qc = orbit["CarterConstant"];
  Cab[r_,\[Theta]_,sgnUr_,sgnU\[Theta]_] := Module[{\[CapitalSigma],\[Rho],ur,u\[Theta],rcomp,\[Theta]comp},
		\[CapitalSigma]=r^2+a^2Cos[\[Theta]]^2;
		\[Rho]=-1/(r-I a Cos[\[Theta]]);
		ur = Sqrt[-((r^2 - 2r + a^2)*(Qc + ((a*E0) - Lz)^2 + r^2)) + (a*Lz - E0*(a^2 + r^2))^2];
		u\[Theta] = Sqrt[Qc - a^2*(1 - E0^2)*Cos[\[Theta]]^2 - (Lz*Cos[\[Theta]]/Sin[\[Theta]])^2];
		rcomp=(E0(r^2+a^2) - a Lz + sgnUr ur)/(2\[CapitalSigma]);
		\[Theta]comp=\[Rho] (I Sin[\[Theta]](a E0 - Lz/Sin[\[Theta]]^2)+sgnU\[Theta] u\[Theta] )/Sqrt[2];
		{rcomp^2, rcomp \[Theta]comp, \[Theta]comp^2}
  ];

  assoc = <|  "s" -> -2,
        "SourceType" -> "PointParticle",
        "Orbit" -> orbit,
        "Cab" -> Cab
      |>;

  TeukolskySourceObject[assoc]
]


TeukolskySourceObject[assoc_][string_] := assoc[string]


(* ::Section::Closed:: *)
(*End Package*)


End[];
EndPackage[];
