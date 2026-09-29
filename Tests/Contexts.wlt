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
