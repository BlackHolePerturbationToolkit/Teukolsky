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
  Module[{R, str, old},
    R = TeukolskyRadial[0, 2, 2, 0.6, 0.5]["In"];
    (* the saved form of a radial function with the old head; with Teukolsky` on the context path InputForm
       prints the short name, which the replacement would miss, so the head is printed fully qualified *)
    str = StringReplace[Block[{$ContextPath = {"System`"}}, ToString[R, InputForm]], "Teukolsky`TeukolskyRadialFunction[" -> "Teukolsky`TeukolskyRadial`TeukolskyRadialFunction["];
    old = ToExpression[str];
    {StringStartsQ[str, "Teukolsky`TeukolskyRadial`TeukolskyRadialFunction["], Head[old] === TeukolskyRadialFunction, old[6.] == R[6.]}
  ],
  {True, True, True},
  TestID -> "A radial function saved with the old head evaluates"
]

VerificationTest[
  {Select[$ContextPath, StringMatchQ[#, "Teukolsky`TeukolskyRadial`" | "Teukolsky`TeukolskyMode`" | "Teukolsky`MST`RenormalizedAngularMomentum`"] &], Context[TeukolskyRadialPN]},
  {{}, "Teukolsky`PN`"},
  TestID -> "The old contexts are off the context path, the PN context remains"
]

(* a reload works without messages (the old-context aliases are removed first) *)
VerificationTest[
  (Get["Teukolsky`"]; NumericQ[TeukolskyRadial[-2, 2, 2, 0.5, 0.3]["In"][10.]]),
  True,
  TestID -> "The package can be reloaded"
]

(* the MST messages are attached to the radial-function symbol of the master package, looked up when first
   needed (for ReggeWheeler it lives in a sub-context, and creating it at load time gave a shadowing copy) *)
VerificationTest[
  {Teukolsky`MST`MST`Private`radialFunctionSymbol[] === TeukolskyRadialFunction, StringQ[MessageName[TeukolskyRadialFunction, "prec"]]},
  {True, True},
  TestID -> "MST messages belong to TeukolskyRadialFunction"
]
