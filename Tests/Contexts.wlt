(* Packages that depend on Teukolsky` must see its public symbols in the Teukolsky` context *)
BeginPackage["TeukolskyContextTest`", {"Teukolsky`"}];
Begin["`Private`"];
contexts = Context /@ {TeukolskyRadial, TeukolskyRadialFunction, TeukolskyMode, TeukolskyPointParticleMode, RenormalizedAngularMomentum};
End[];
EndPackage[];

VerificationTest[
  TeukolskyContextTest`Private`contexts,
  ConstantArray["Teukolsky`", 5],
  TestID -> "Public symbols are in the Teukolsky` context"
]

VerificationTest[
  StringQ[TeukolskyRadial::usage] && StringQ[RenormalizedAngularMomentum::usage],
  True,
  TestID -> "Usage messages are attached to the Teukolsky` symbols"
]


(* The names the public symbols had before they moved to Teukolsky` are aliases: fully qualified references
   and expressions saved with the old heads (as in the Documentation notebooks) keep working *)
VerificationTest[
  {ToExpression["Teukolsky`TeukolskyRadial`TeukolskyRadial"] === TeukolskyRadial,
   ToExpression["Teukolsky`TeukolskyRadial`TeukolskyRadialFunction"] === TeukolskyRadialFunction,
   ToExpression["Teukolsky`TeukolskyMode`TeukolskyMode"] === TeukolskyMode,
   ToExpression["Teukolsky`TeukolskyMode`TeukolskyPointParticleMode"] === TeukolskyPointParticleMode,
   ToExpression["Teukolsky`MST`RenormalizedAngularMomentum`RenormalizedAngularMomentum"] === RenormalizedAngularMomentum},
  {True, True, True, True, True},
  TestID -> "Old fully qualified names are aliases of the Teukolsky` symbols"
]

VerificationTest[
  Module[{R, old},
    R = TeukolskyRadial[0, 2, 2, 0.6, 0.5]["In"];
    (* the saved form of a radial function with the old head *)
    old = ToExpression[StringReplace[ToString[R, InputForm], "Teukolsky`TeukolskyRadialFunction[" -> "Teukolsky`TeukolskyRadial`TeukolskyRadialFunction["]];
    {Head[old] === TeukolskyRadialFunction, old[6.] == R[6.]}
  ],
  {True, True},
  TestID -> "A radial function saved with the old head evaluates"
]

VerificationTest[
  {Select[$ContextPath, StringMatchQ[#, "Teukolsky`TeukolskyRadial`" | "Teukolsky`TeukolskyMode`" | "Teukolsky`MST`RenormalizedAngularMomentum`"] &], Context[TeukolskyRadialPN]},
  {{}, "Teukolsky`PN`"},
  TestID -> "The old contexts are off the context path, the PN context remains"
]
