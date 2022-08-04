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
    


TeukolskyPointParticleSourceCircular[s:0, orbit_] :=
  TeukolskySourceObject[<|"s" -> s, "SourceType" -> "PointParticle", "Orbit" -> orbit|>];


TeukolskyPointParticleSourceSpherical[s:0, orbit_] :=
  TeukolskySourceObject[<|"s" -> s, "SourceType" -> "PointParticle", "Orbit" -> orbit|>];


TeukolskyPointParticleSourceEccentric[s:0, orbit_] :=
  TeukolskySourceObject[<|"s" -> s, "SourceType" -> "PointParticle", "Orbit" -> orbit|>];


TeukolskyPointParticleSourceGeneric[s:0, orbit_] :=
  TeukolskySourceObject[<|"s" -> s, "SourceType" -> "PointParticle", "Orbit" -> orbit|>];


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
