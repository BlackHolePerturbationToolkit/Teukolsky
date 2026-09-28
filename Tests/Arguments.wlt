(* Exact arguments *)
VerificationTest[
  TeukolskyRadial[0, 2, 2, 6/10, 1/2],
  $Failed,
  {TeukolskyRadial::exact},
  TestID -> "Exact arguments without a WorkingPrecision fail with a message"
]

VerificationTest[
  R = TeukolskyRadial[0, 2, 2, 6/10, 1/2, WorkingPrecision -> 30];
  Rn = TeukolskyRadial[0, 2, 2, N[6/10, 30], N[1/2, 30]];
  Abs[R["In"][N[10, 30]]/Rn["In"][N[10, 30]] - 1] < 10^-25 && Precision[R["In"][N[10, 30]]] >= 28,
  True,
  TestID -> "Exact arguments are evaluated at the given WorkingPrecision"
]

VerificationTest[
  Precision[N[RenormalizedAngularMomentum[0, 2, 2, 6/10, 1/10], 30]],
  30.,
  TestID -> "N applied to RenormalizedAngularMomentum with exact arguments"
]
