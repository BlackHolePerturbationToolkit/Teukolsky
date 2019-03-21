(* ::Package:: *)

(* ::Chapter:: *)
(*TeukolskySource*)


BeginPackage["Teukolsky`TeukolskySource`"];

TeukolskySourceObject::usage = "TeukolskySourceObject[assoc] an object which contains a Teukolsky source"

TeukolskyPointParticleSource::usage = "TeukolskyPointParticleSource[s, orbit] Point particle source for the Teukolsky equation"

Begin["`Private`"];


TeukolskyPointParticleSource[s_, orbit_]:=Module[{assoc, source},

assoc = <| "s" -> s, "sourceType" -> "PointParticle" |>;

TeukolskySourceObject[assoc]
]


Format[TeukolskySourceObject[assoc_]] := "TeukolskySourceObject[<<>>]";

TeukolskySourceObject[assoc_][string_] := assoc[string]


End[];
EndPackage[];
