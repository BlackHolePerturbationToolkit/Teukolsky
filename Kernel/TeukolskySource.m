(* ::Package:: *)

(* ::Title:: *)
(*TeukolskySource*)


(* ::Section::Closed:: *)
(*Create Package*)


(* ::Subsection::Closed:: *)
(*BeginPackage*)


BeginPackage["Teukolsky`TeukolskySource`",{"KerrGeodesics`OrbitalFrequencies`"}];


(* ::Subsection:: *)
(*Usage messages*)


(* ::Subsection::Closed:: *)
(*Begin Private section*)


Begin["`Private`"];


(* ::Section::Closed:: *)
(*TeukolskyPointParticleSource*)


TeukolskyPointParticleSource[s_, orbit_] :=
Module[{},
  If[(s==-2 || s==-1 || s==+1 || s==0) && orbit["e"] == 0 && Abs[orbit["Inclination"]] == 1,
    Return[TeukolskyPointParticleSourceCircular[s,orbit]]];
  If[(s==0||s==-2) && orbit["e"] == 0 && Abs[orbit["Inclination"]] != 1,
    Return[TeukolskyPointParticleSourceSpherical[s,orbit]]];
    
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
	{\[CapitalOmega]r,\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi]} = Values[Evaluate[KerrGeoFrequencies[a,r0,e,x]]];
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
{a,r0,e,x,\[ScriptCapitalE],\[ScriptCapitalL],\[ScriptCapitalQ],\[Beta],zm,zp,\[Theta],z,\[CapitalDelta],\[CapitalSigma],\[Rho],\[Rho]b,ut,\[CapitalTheta],\[Chi],CnnPlus,CnnMinus,CnmPlus,CnmMinus,CmmPlus,CmmMinus,assoc},
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
\[Beta]=a^2 (1-\[ScriptCapitalE]^2);
{zm,zp}={(\[ScriptCapitalL]^2+\[ScriptCapitalQ]+\[Beta]-\[Beta] Sqrt[(\[ScriptCapitalL]^4+(\[ScriptCapitalQ]-\[Beta])^2+2 \[ScriptCapitalL]^2 (\[ScriptCapitalQ]+\[Beta]))/\[Beta]^2])/(2 \[Beta]),(\[ScriptCapitalL]^2+\[ScriptCapitalQ]+\[Beta]+\[Beta] Sqrt[(\[ScriptCapitalL]^4+(\[ScriptCapitalQ]-\[Beta])^2+2 \[ScriptCapitalL]^2 (\[ScriptCapitalQ]+\[Beta]))/\[Beta]^2])/(2 \[Beta])};

\[Theta]=ArcCos[Sqrt[zm] Cos[\[Chi]]];

\[CapitalDelta]=r0^2-2r0+a^2;
\[CapitalSigma]=r0^2+a^2 Cos[\[Theta]]^2;
\[Rho]=-1/(r0-I a Cos[\[Theta]]);
\[Rho]b=-1/(r0+I a Cos[\[Theta]]);

ut=1/\[CapitalSigma] (\[ScriptCapitalE]((r0^2+a^2)^2/\[CapitalDelta]-a^2 Sin[\[Theta]]^2)+a \[ScriptCapitalL](1-(r0^2+a^2)/\[CapitalDelta]));
\[CapitalTheta]=Sqrt[\[ScriptCapitalQ]-Cot[\[Theta]]^2 \[ScriptCapitalL]^2-a^2 Cos[\[Theta]]^2 (1-\[ScriptCapitalE]^2)];

(*Stress energy projected onto tetrad (\[PlusMinus]\[CapitalTheta] for travelling up or down)*)
CnnPlus=1/(4 \[CapitalSigma]^3 ut) (\[ScriptCapitalE](r0^2+a^2)-a \[ScriptCapitalL])^2;
CnmPlus=\[Rho]/(2Sqrt[2] \[CapitalSigma]^2 ut) (\[ScriptCapitalE](r0^2+a^2)-a \[ScriptCapitalL])(I Sin[\[Theta]](a \[ScriptCapitalE]-\[ScriptCapitalL]/Sin[\[Theta]]^2)+ \[CapitalTheta]);
CmmPlus=\[Rho]^2/(2\[CapitalSigma] ut) (I Sin[\[Theta]](a \[ScriptCapitalE]-\[ScriptCapitalL]/Sin[\[Theta]]^2)+ \[CapitalTheta])^2;

CnnMinus=1/(4 \[CapitalSigma]^3 ut) (\[ScriptCapitalE](r0^2+a^2)-a \[ScriptCapitalL])^2;
CnmMinus=\[Rho]/(2Sqrt[2] \[CapitalSigma]^2 ut) (\[ScriptCapitalE](r0^2+a^2)-a \[ScriptCapitalL])(I Sin[\[Theta]](a \[ScriptCapitalE]-\[ScriptCapitalL]/Sin[\[Theta]]^2)- \[CapitalTheta]);
CmmMinus=\[Rho]^2/(2\[CapitalSigma] ut) (I Sin[\[Theta]](a \[ScriptCapitalE]-\[ScriptCapitalL]/Sin[\[Theta]]^2)- \[CapitalTheta])^2;
(*Rewriting source*)


assoc = <|"s"->-2, "SourceType"->"PointParticle","Orbit"->orbit,"Cnn+"->CnnPlus,"Cnn-"->CnnMinus,"Cnm+"->CnmPlus,"Cnm-"->CnmMinus,"Cmm+"->CmmPlus,"Cmm-"->CmmMinus,"\[Chi]"->\[Chi]|>
]


TeukolskySourceObject[assoc_][string_] := assoc[string]


(* ::Section::Closed:: *)
(*End Package*)


End[];
EndPackage[];
