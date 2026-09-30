(* Mathematica Test File *)

(****************************************************************)
(* ZeroFrequency                                                *)
(****************************************************************)
VerificationTest[
    RenormalizedAngularMomentum[-2, 2, 0, 0, 0, 4]
    ,
    2
    ,
    TestID->"ZeroFrequency"
]


(****************************************************************)
(* ZeroFrequencyMachinePrecision                                *)
(****************************************************************)
VerificationTest[
    RenormalizedAngularMomentum[-2, 2, 0, 0, 0., 4]
    ,
    2
    ,
    TestID->"ZeroFrequency"
]


(****************************************************************)
(* ZeroFrequencyMonodromy                                       *)
(****************************************************************)
VerificationTest[
    RenormalizedAngularMomentum[-2, 2, 0, 0, 0, 4, Method -> "Monodromy"]
    ,
    2
    ,
    TestID->"ZeroFrequencyMonodromy"
]




(****************************************************************)
(* Representative continuous across the branch points            *)
(****************************************************************)
(* For s = -2, l = m = 2, a = 3/5 the monodromy cosine drops below -1 between omega = 0.40 and 0.42 and nu turns
   complex; the representative l - ArcCos[Cos[2 Pi nu]]/(2 Pi), used on every branch, goes continuously from
   real values near l - 1/2 to l - 1/2 + i y (the 1/2 + i y returned before differs by the integer l - 1) *)
VerificationTest[
  Module[{nus},
    nus = Table[RenormalizedAngularMomentum[-2, 2, 2, N[3/5, 32], N[w, 32]], {w, 30/100, 60/100, 2/100}];
    {Max[Abs[Differences[nus]]] < 1/5, Abs[Re[nus[[-1]]] - 3/2] < 10^-25, Abs[Re[RenormalizedAngularMomentum[0, 2, 0, N[3/5, 32], N[1, 32]]] - 2] < 10^-25}
  ],
  {True, True, True},
  TestID -> "Renormalized angular momentum is continuous where it turns complex, on both branches"
]

(* the representative is a convention: the radial functions and amplitudes are the same for the equivalent value *)
VerificationTest[
  Module[{q = N[3/5, 32], w = N[1/2, 32], nu, Rn, Ro, r = N[6, 32]},
    nu = RenormalizedAngularMomentum[-2, 2, 2, q, w];
    Rn = TeukolskyRadial[-2, 2, 2, q, w, Method -> "MST"];
    Ro = TeukolskyRadial[-2, 2, 2, q, w, Method -> "MST", "RenormalizedAngularMomentum" -> 3 - nu];
    N[Max[Abs[{Rn["In"][r]/Ro["In"][r], Rn["Up"][r]/Ro["Up"][r], Rn["In"]["Amplitudes"]["Incidence"]/Ro["In"]["Amplitudes"]["Incidence"]} - 1]]] < 10^-28
  ],
  True,
  TestID -> "Equivalent representatives of nu give the same solutions"
]


(****************************************************************)
(* Failure is reported                                          *)
(****************************************************************)
(* at omega = 270 the monodromy recurrences overflow: $Failed, with the convergence message *)
VerificationTest[
  RenormalizedAngularMomentum[-2, 2, 2, N[3/5, 32], N[270, 32]],
  $Failed,
  {RenormalizedAngularMomentum::conv},
  TestID -> "A failed monodromy evaluation returns $Failed with a message"
]
