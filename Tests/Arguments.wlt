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
  Max[Abs[{R[10., 0] - R[10.], R[10., 1] - R'[10.], R[10., 2] - R''[10.]}]] == 0,
  True,
  TestID -> "R[r, n] gives the n-th derivative"
]

(* Options of TeukolskyRadial given inside Method *)
VerificationTest[
  TeukolskyRadial[-2, 2, 2, 0.6, 0.3, Method -> {"MST", "RenormalizedAngularMomentum" -> 1.8}]["In"]["Method"],
  {"MST"},
  {TeukolskyRadial::topopt},
  TestID -> "Options of TeukolskyRadial inside Method are reported as misplaced"
]

(* Value and derivative from one summation *)
VerificationTest[
  Module[{R = TeukolskyRadial[-2, 2, 2, N[6/10, 32], N[3/10, 32], Method -> "MST"], r = N[20, 32]},
    Max[Abs[Join[R["In"][r, {0, 1}]/{R["In"][r], R["In"]'[r]}, R["Up"][r, {0, 1}]/{R["Up"][r], R["Up"]'[r]}] - 1]] < 10^-29
  ],
  True,
  TestID -> "R[r, {0, 1}] agrees with the separate MST evaluations"
]

VerificationTest[
  Module[{R = TeukolskyRadial[-2, 2, 2, 0.6, 0.3]["In"]},
    R[10., {0, 1}] == {R[10.], R'[10.]} && R[10., {0, 1, 2}] == {R[10.], R'[10.], R''[10.]}
  ],
  True,
  TestID -> "R[r, {0, 1}] and R[r, list of orders] for a numerically integrated solution"
]
