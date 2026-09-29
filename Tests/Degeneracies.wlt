(* Degeneracies of the MST method on the negative imaginary axis, omega = -I sigma (M = 1) *)

(* 2 I epsilon = 4 sigma an integer: the monodromy method for nu degenerates; nu is evaluated from
   neighbouring frequencies (it used to run until the kernel died) *)
VerificationTest[
  Module[{nu, nb},
    nu = Quiet[RenormalizedAngularMomentum[0, 2, 0, SetPrecision[3/5, 40], SetPrecision[-I/4, 40]], RenormalizedAngularMomentum::degenerate];
    nb = (RenormalizedAngularMomentum[0, 2, 0, SetPrecision[3/5, 40], SetPrecision[-I/4 (1 + 10^-8), 40]] + RenormalizedAngularMomentum[0, 2, 0, SetPrecision[3/5, 40], SetPrecision[-I/4 (1 - 10^-8), 40]])/2;
    Abs[nu/nb - 1] < 10^-12
  ],
  True,
  TestID -> "Renormalized angular momentum at a monodromy degeneracy"
]

VerificationTest[
  RenormalizedAngularMomentum[0, 2, 0, SetPrecision[3/5, 40], SetPrecision[-I/2, 40]],
  _?NumericQ,
  {RenormalizedAngularMomentum::degenerate},
  SameTest -> MatchQ,
  TestID -> "Monodromy degeneracy at sigma = 1/2 is reported"
]

(* 2 I epsilon_+ an integer: for a = 3/5, m = 0 that is 4.5 sigma. With n = 2 I epsilon_+ >= 1 - s the
   hypergeometric series cannot represent the "In" solution *)
VerificationTest[
  TeukolskyRadial[0, 2, 0, SetPrecision[3/5, 40], SetPrecision[-2 I, 40]],
  $Failed,
  {TeukolskyRadial::insing},
  TestID -> "In solution not representable at 2 I epsilon_+ = 9, s = 0"
]

(* n = 2 < 1 - s for s = -2: the amplitudes are the limit from neighbouring frequencies and the functions are fine *)
VerificationTest[
  Module[{R, Rp},
    R = Quiet[TeukolskyRadial[-2, 2, 0, SetPrecision[3/5, 40], SetPrecision[-4 I/9, 40]], TeukolskyRadial::degenerate];
    Rp = TeukolskyRadial[-2, 2, 0, SetPrecision[3/5, 40], SetPrecision[-4 I/9 (1 + 10^-8), 40]];
    {Abs[R["In"][N[6, 40]]/Rp["In"][N[6, 40]] - 1] < 10^-6, Abs[R["In"]["Amplitudes"]["Incidence"]/Rp["In"]["Amplitudes"]["Incidence"] - 1] < 10^-6, R["Up"]["Amplitudes"]["Reflection"]}
  ],
  {True, True, Indeterminate},
  TestID -> "Amplitudes at 2 I epsilon_+ = 2, s = -2, are the limit of their neighbours"
]
