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

(* Forwarding aliases for the contexts the public symbols lived in before they were moved to Teukolsky`, so
   that fully qualified references and expressions saved with the old heads (the Documentation notebooks hold
   Teukolsky`TeukolskyRadial`TeukolskyRadialFunction and Teukolsky`TeukolskyMode`TeukolskyMode) keep
   evaluating: an expression with an old head evaluates to the same expression with the new one. Those three
   contexts export nothing any more and are taken off $ContextPath first, so that the aliases cannot shadow
   the Teukolsky` symbols when a short name is typed; the symbols are created at evaluation time, since the
   shadowing warning is issued when a symbol whose name exists in Teukolsky` is created. *)
$ContextPath = DeleteCases[$ContextPath, "Teukolsky`TeukolskyRadial`" | "Teukolsky`TeukolskyMode`" | "Teukolsky`MST`RenormalizedAngularMomentum`"];
Teukolsky`Private`alias[old_String, new_Symbol] := Quiet[With[{sym = Symbol[old]}, sym = new; Protect[sym]], General::shdw];
Teukolsky`Private`alias["Teukolsky`TeukolskyRadial`TeukolskyRadial", Teukolsky`TeukolskyRadial];
Teukolsky`Private`alias["Teukolsky`TeukolskyRadial`TeukolskyRadialFunction", Teukolsky`TeukolskyRadialFunction];
Teukolsky`Private`alias["Teukolsky`TeukolskyMode`TeukolskyMode", Teukolsky`TeukolskyMode];
Teukolsky`Private`alias["Teukolsky`TeukolskyMode`TeukolskyPointParticleMode", Teukolsky`TeukolskyPointParticleMode];
Teukolsky`Private`alias["Teukolsky`MST`RenormalizedAngularMomentum`RenormalizedAngularMomentum", Teukolsky`RenormalizedAngularMomentum];
