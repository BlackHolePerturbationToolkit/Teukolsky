(* ::Package:: *)

BeginPackage["Teukolsky`"];

(* Public symbols are declared here, in the Teukolsky` context, so that packages depending on this one
   (BeginPackage["X`", {"Teukolsky`"}]) see them; their usage messages and definitions are attached
   by the sub-packages loaded below, each of which has Teukolsky` in its list of needed contexts. *)
{TeukolskyRadial, TeukolskyRadialFunction, TeukolskyMode, TeukolskyPointParticleMode, RenormalizedAngularMomentum};

Begin["`Private`"];

(* Check appropriate versions of dependencies are installed *)
versionCheck[name_, version_] :=
 Module[{paclet},
  paclet = PacletObject[name];
  If[!(PacletNewerQ[paclet, version] || paclet["Version"] == version),
    Throw[Failure["PacletNotFound",
      <|"MessageTemplate" ->"Paclet `1` version `2` or greater not found.", 
        "MessageParameters" -> {name, version}|>], Teukolsky, #1&]
  ]
];
versionCheck["SpinWeightedSpheroidalHarmonics", "1.0.0"];
versionCheck["KerrGeodesics", "0.9.0"];

End[];

EndPackage[];

Block[{MST`$MasterFunction = "Teukolsky"},
  Get["Teukolsky`MST`RenormalizedAngularMomentum`"];
  Get["Teukolsky`MST`MST`"];
];

Get["Teukolsky`SasakiNakamura`"];
Get["Teukolsky`TeukolskyRadial`"];
Get["Teukolsky`TeukolskyMode`"];
Get["Teukolsky`NumericalIntegration`"];
Get["Teukolsky`ConvolveSource`"];
Get["Teukolsky`PN`"];
