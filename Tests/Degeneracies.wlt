(* Degeneracies of the MST method on the negative imaginary axis, omega = -I sigma (M = 1) *)

(* 2 I epsilon = 4 sigma an integer: Gamma[mu1 - mu2] in the monodromy formula has a pole (the method used to
   run until the kernel died, then extrapolated from neighbouring frequencies); with Gamma[mu1 - mu2]
   Pochhammer[mu1 - mu2, k] combined into Gamma[mu1 - mu2 + k] it is evaluated directly, to the precision
   of any other frequency, and agrees with the neighbours' average to their O(h^2) *)
VerificationTest[
  Module[{nu, nb},
    nu = RenormalizedAngularMomentum[0, 2, 0, SetPrecision[3/5, 40], SetPrecision[-I/4, 40]];
    nb = (RenormalizedAngularMomentum[0, 2, 0, SetPrecision[3/5, 40], SetPrecision[-I/4 (1 + 10^-8), 40]] + RenormalizedAngularMomentum[0, 2, 0, SetPrecision[3/5, 40], SetPrecision[-I/4 (1 - 10^-8), 40]])/2;
    {Abs[nu/nb - 1] < 10^-15, Precision[nu] > 28}
  ],
  {True, True},
  TestID -> "Renormalized angular momentum at a monodromy degeneracy"
]

VerificationTest[
  Module[{nu40, nu60},
    nu40 = RenormalizedAngularMomentum[-2, 2, 2, SetPrecision[3/5, 40], SetPrecision[-I/2, 40]];
    nu60 = RenormalizedAngularMomentum[-2, 2, 2, SetPrecision[3/5, 60], SetPrecision[-I/2, 60]];
    {Abs[nu40/nu60 - 1] < 10^-28, Precision[nu60] > 45}
  ],
  {True, True},
  TestID -> "Monodromy degeneracy at sigma = 1/2, s = -2: full precision, no message"
]

(* 2 I epsilon_+ an integer: for a = 3/5, m = 0 that is 4.5 sigma. With n = 2 I epsilon_+ >= 1 - s the
   transmission amplitude of the "In" solution vanishes and it is normalised to unit incidence: the limit of
   the unit-incidence "In" solution of the neighbours *)
(* omega = -2 I is also a degeneracy of the monodromy method (2 I epsilon = 8), evaluated directly *)
VerificationTest[
  Module[{R, Rp, r = N[6, 40]},
    R = Quiet[TeukolskyRadial[0, 2, 0, SetPrecision[3/5, 40], SetPrecision[-2 I, 40]], {TeukolskyRadial::degenerate, TeukolskyRadial::innorm}];
    Rp = TeukolskyRadial[0, 2, 0, SetPrecision[3/5, 40], SetPrecision[-2 I (1 + 10^-8), 40]];
    {Abs[R["In"]["Amplitudes"]["Incidence"] - 1] < 10^-30, R["In"]["Amplitudes"]["Transmission"],
     Abs[R["In"][r]/(Rp["In"][r]/Rp["In"]["Amplitudes"]["Incidence"]) - 1] < 10^-6,
     Abs[R["In"]["Amplitudes"]["Reflection"]/(Rp["In"]["Amplitudes"]["Reflection"]/Rp["In"]["Amplitudes"]["Incidence"]) - 1] < 10^-6}
  ],
  {True, 0, True, True},
  TestID -> "In solution normalised to unit incidence at 2 I epsilon_+ = 9, s = 0"
]

VerificationTest[
  TeukolskyRadial[0, 2, 0, SetPrecision[3/5, 40], SetPrecision[-2 I, 40]],
  _Association,
  {TeukolskyRadial::degenerate, TeukolskyRadial::innorm},
  SameTest -> MatchQ,
  TestID -> "Messages at 2 I epsilon_+ = 9, s = 0"
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
