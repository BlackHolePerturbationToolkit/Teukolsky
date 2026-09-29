(* The superradiant bound frequency omega = m Omega_H: for a = 3/5, Omega_H = 1/6, so l = m = 3 at omega = 1/2 *)

VerificationTest[
  R = TeukolskyRadial[0, 3, 3, N[3/5, 50], N[1/2, 50]];
  W = With[{r = N[6, 50]}, (r^2 - 2 r + 9/25) (R["In"][r] R["Up"]'[r] - R["In"]'[r] R["Up"][r])];
  Abs[W/(2 I 1/2 R["In"]["Amplitudes"]["Incidence"]) - 1] < 10^-30,
  True,
  {TeukolskyRadial::superradiant},
  TestID -> "Wronskian at the superradiant bound frequency, s = 0"
]

VerificationTest[
  R = Quiet[TeukolskyRadial[0, 3, 3, N[3/5, 50], N[1/2, 50]]];
  Rp = TeukolskyRadial[0, 3, 3, N[3/5, 50], N[1/2 + 10^-7, 50]];
  Rm = TeukolskyRadial[0, 3, 3, N[3/5, 50], N[1/2 - 10^-7, 50]];
  {Abs[R["In"][N[6, 50]]/((Rp["In"][N[6, 50]] + Rm["In"][N[6, 50]])/2) - 1] < 10^-10,
   Abs[R["Up"][N[6, 50]]/((Rp["Up"][N[6, 50]] + Rm["Up"][N[6, 50]])/2) - 1] < 10^-10,
   Abs[R["In"]["Amplitudes"]["Incidence"]/((Rp["In"]["Amplitudes"]["Incidence"] + Rm["In"]["Amplitudes"]["Incidence"])/2) - 1] < 10^-10,
   R["Up"]["Amplitudes"]["Incidence"], R["Up"]["Amplitudes"]["Reflection"]},
  {True, True, True, Indeterminate, Indeterminate},
  TestID -> "Radial functions and amplitudes at the superradiant bound frequency are the limit of their neighbours, s = 0"
]

VerificationTest[
  R = Quiet[TeukolskyRadial[-2, 3, 3, N[3/5, 50], N[1/2, 50]]];
  W = With[{r = N[6, 50]}, (r^2 - 2 r + 9/25)^(-1) (R["In"][r] R["Up"]'[r] - R["In"]'[r] R["Up"][r])];
  {Abs[W/(2 I 1/2 R["In"]["Amplitudes"]["Incidence"]) - 1] < 10^-30, NumericQ[R["Up"]["Amplitudes"]["Incidence"]], R["Up"]["Amplitudes"]["Reflection"]},
  {True, True, Indeterminate},
  TestID -> "Wronskian and Up amplitudes at the superradiant bound frequency, s = -2"
]

VerificationTest[
  R = Quiet[TeukolskyRadial[0, 3, 3, 0.6, 0.5]];
  Rp = TeukolskyRadial[0, 3, 3, 0.6, 0.5 + 10^-7];
  Abs[R["In"][6.]/Rp["In"][6.] - 1] < 10^-5 && Abs[R["In"]["Amplitudes"]["Incidence"]/Rp["In"]["Amplitudes"]["Incidence"] - 1] < 10^-5,
  True,
  TestID -> "Superradiant bound frequency at machine precision"
]

(* For s >= 1 the "In" solution is the smaller-exponent member of the resonant horizon basis and has no
   unit-transmission limit at the bound frequency: refused rather than evaluated (which used to hang) *)
VerificationTest[
  TeukolskyRadial[2, 3, 3, N[3/5, 30], N[1/2, 30]],
  $Failed,
  {TeukolskyRadial::insing},
  TestID -> "In solution refused at the superradiant bound frequency, s = 2"
]

VerificationTest[
  R = Quiet[TeukolskyRadial[2, 3, 3, N[3/5, 30], N[1/2, 30], "BoundaryConditions" -> "Up"]];
  Rp = TeukolskyRadial[2, 3, 3, N[3/5, 30], N[1/2 + 10^-7, 30]];
  Rm = TeukolskyRadial[2, 3, 3, N[3/5, 30], N[1/2 - 10^-7, 30]];
  {Abs[R[N[6, 30]]/((Rp["Up"][N[6, 30]] + Rm["Up"][N[6, 30]])/2) - 1] < 10^-8,
   NumericQ[R["Amplitudes"]["Reflection"]], R["Amplitudes"]["Incidence"]},
  {True, True, Indeterminate},
  TestID -> "Up solution at the superradiant bound frequency, s = 2"
]
