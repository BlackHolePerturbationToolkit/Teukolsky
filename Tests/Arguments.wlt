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

(* Machine precision without the asymptotic amplitudes or nu: the integration from series boundary data needs
   neither, the functions are the same, and no accuracy warning is issued (the estimate needs the amplitudes) *)
VerificationTest[
  Module[{Rd, Ra, Rn, r = 6.},
    Rd = TeukolskyRadial[-2, 2, 2, 0.6, 0.5];
    Ra = TeukolskyRadial[-2, 2, 2, 0.6, 0.5, "Amplitudes" -> False];
    Rn = TeukolskyRadial[-2, 2, 2, 0.6, 0.5, "Amplitudes" -> False, "RenormalizedAngularMomentum" -> False];
    {Max[Abs[{Ra["In"][r], Ra["Up"][r], Rn["In"][r], Rn["Up"][r]}/{Rd["In"][r], Rd["Up"][r], Rd["In"][r], Rd["Up"][r]} - 1]] < 10^-12,
     Ra["In"]["Amplitudes"], Rn["In"]["RenormalizedAngularMomentum"]}
  ],
  {True, <|"Transmission" -> 1|>, Indeterminate},
  TestID -> "Machine precision without amplitudes or nu gives the same functions and no messages"
]

(* A single boundary condition builds only that solution, and no accuracy warning is attempted *)
VerificationTest[
  Module[{Rd, Ri, Ru, r = 6.},
    Rd = TeukolskyRadial[0, 2, 0, 0.6, 0.5];
    Ri = TeukolskyRadial[0, 2, 0, 0.6, 0.5, "BoundaryConditions" -> "In"];
    Ru = TeukolskyRadial[0, 2, 0, 0.6, 0.5, "BoundaryConditions" -> "Up", "Amplitudes" -> False];
    {Head[Ri], Head[Ru], Abs[Ri[r]/Rd["In"][r] - 1] < 10^-12, Abs[Ru[r]/Rd["Up"][r] - 1] < 10^-12}
  ],
  {TeukolskyRadialFunction, TeukolskyRadialFunction, True, True},
  TestID -> "A single boundary condition at machine precision, with and without amplitudes"
]


(* A frequency far outside the range of the methods fails quickly (it used to run for minutes) *)
VerificationTest[
  Module[{t, res},
    {t, res} = AbsoluteTiming[Quiet[TeukolskyRadial[-2, 2, 2, 0.6, 270.]]];
    {res, t < 60}
  ],
  {$Failed, True},
  TestID -> "An absurd frequency fails fast"
]

(* At omega = 50 the amplitude formulae cancel to zeros of no precision at 32 digits; they used to be taken for
   exact zeros (and reported as an overflow), and are now padded until they carry the working precision *)
VerificationTest[
  Module[{R, b},
    R = TeukolskyRadial[-2, 2, 2, 0.6, 50.];
    b = R["In"]["Amplitudes"]["Incidence"];
    {NumericQ[R["In"][6.]], NumericQ[R["Up"][6.]], NumericQ[b] && b != 0}
  ],
  {True, True, True},
  TestID -> "Amplitudes that cancel to zero at large frequency are padded"
]

(* an invalid "WronskianCheck" value is rejected rather than silently disabling the check *)
VerificationTest[
  TeukolskyRadial[-2, 2, 2, 0.6, 0.5, "WronskianCheck" -> "True"],
  $Failed,
  {TeukolskyRadial::optx},
  TestID -> "Invalid WronskianCheck option fails with a message"
]

(* Unit-incidence normalisation is reserved for the "In" solution at a degeneracy: a zero transmission of the
   "Up" solution, or of the "In" solution away from a degeneracy, is an overflow of the amplitude formulae *)
VerificationTest[
  Module[{key = Teukolsky`TeukolskyRadial`Private`normalisationKey, ns = <|"Incidence" -> 1., "Transmission" -> 0., "Reflection" -> 1.|>},
    {key[ns, "In"], key[ns, "Up"], Block[{Teukolsky`TeukolskyRadial`Private`$degenerateIn = True}, {key[ns, "In"], key[ns, "Up"]}]}
  ],
  {"Transmission", "Transmission", {"Incidence", "Transmission"}},
  TestID -> "Normalisation key: unit incidence only for In at a degeneracy"
]

(* the same validation in the static path *)
VerificationTest[
  {TeukolskyRadial[-2, 2, 2, 0.6, 0, "WronskianCheck" -> "True"], Head[TeukolskyRadial[-2, 2, 2, 0.6, 0, "WronskianCheck" -> True]]},
  {$Failed, Association},
  {TeukolskyRadial::optx, TeukolskyRadial::sopt},
  TestID -> "WronskianCheck is validated and reported as unsupported for static modes"
]

(* a forced Wronskian check needs the amplitudes *)
VerificationTest[
  TeukolskyRadial[-2, 2, 2, N[3/5, 32], N[1/2, 32], "Amplitudes" -> False, "WronskianCheck" -> True],
  $Failed,
  {TeukolskyRadial::opti},
  TestID -> "A forced Wronskian check without amplitudes is refused"
]

(* the Wronskian check also applies to numerically integrated solutions, inside their domains; the HeunC
   method has no check and refuses it *)
VerificationTest[
  {Head[TeukolskyRadial[-2, 2, 2, 0.6, 0.5, Method -> "NumericalIntegration", "WronskianCheck" -> True]],
   Head[TeukolskyRadial[-2, 2, 2, N[3/5, 24], N[1/2, 24], Method -> {"NumericalIntegration", "Domain" -> {"In" -> {3, 8}, "Up" -> {5, 12}}}, "WronskianCheck" -> True]]},
  {Association, Association},
  TestID -> "A forced Wronskian check with numerical integration runs"
]
