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
  \[CapitalDelta] = r0^2 - 2r0 + a^2;
  ut = orbit["FourVelocity"][[1]][0];
  
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


TeukolskyPointParticleSourceCircular[s:-2, orbit_] :=
  TeukolskySourceObject[<|"s" -> s, "SourceType" -> "PointParticle", "Orbit" -> orbit|>]


TeukolskyPointParticleSourceSpherical[s:-2, orbit_] :=
  TeukolskySourceObject[<|"s" -> s, "SourceType" -> "PointParticle", "Orbit" -> orbit|>]


TeukolskyPointParticleSourceEccentric[s:-2, orbit_] :=
  TeukolskySourceObject[<|"s" -> s, "SourceType" -> "PointParticle", "Orbit" -> orbit|>]


TeukolskyPointParticleSourceGeneric[s:-2, orbit_] :=
  TeukolskySourceObject[<|"s" -> s, "SourceType" -> "PointParticle", "Orbit" -> orbit|>];


TeukolskySourceObject[assoc_][string_] := assoc[string]


(* ::Section::Closed:: *)
(*End Package*)


End[];
EndPackage[];
