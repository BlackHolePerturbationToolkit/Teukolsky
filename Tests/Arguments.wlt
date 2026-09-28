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

(* Derivatives through R[r, n] *)
VerificationTest[
  R = TeukolskyRadial[-2, 2, 2, 0.6, 0.3]["In"];
  {R[10., 0] - R[10.], R[10., 1] - R'[10.], R[10., 2] - R''[10.]},
  {0., 0., 0.},
  TestID -> "R[r, n] gives the n-th derivative"
]

(* Options of TeukolskyRadial given inside Method *)
VerificationTest[
  TeukolskyRadial[-2, 2, 2, 0.6, 0.3, Method -> {"MST", "RenormalizedAngularMomentum" -> 1.8}]["In"]["Method"],
  {"MST"},
  {TeukolskyRadial::topopt},
  TestID -> "Options of TeukolskyRadial inside Method are reported as misplaced"
]
