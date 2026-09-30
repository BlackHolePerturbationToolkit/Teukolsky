(* The MST "Up" solution near the horizon in the horizon basis of hypergeometric series (the "In" series and
   the series on the normalised second Kummer solution), with the connection coefficients from the horizon
   amplitudes; the Coulomb-type series, an expansion about infinity, is kept beyond r+ + 1 and at the
   degeneracies 2 I epsilon_+ = n. References: the Coulomb-type series forced on where it is still affordable. *)

horizonVsCoulomb[s_, l_, m_, a_, om_, dr_] :=
 Module[{wp = 32, R, r, vh, dvh, vc, dvc},
  R = TeukolskyRadial[s, l, m, N[a, wp], N[om, wp], Method -> "MST"];
  r = N[1 + Sqrt[1 - a^2] + dr, wp];
  {vh, dvh} = {R["Up"][r], R["Up"]'[r]};
  {vc, dvc} = Block[{Teukolsky`MST`MST`Private`$MSTUpHorizonThreshold = 0}, {R["Up"][r], R["Up"]'[r]}];
  N[Max[Abs[{vh/vc - 1, dvh/dvc - 1}]]]
 ];

VerificationTest[
  Table[horizonVsCoulomb[s, 2, 2, 3/5, 1/2, 3/10] < 10^-28, {s, {-2, 0, 2}}],
  {True, True, True},
  TestID -> "Up in the horizon basis agrees with the Coulomb-type series at r+ + 3/10, value and derivative"
]

wronskianErr[s_, l_, m_, a_, om_, dr_, wp_] :=
 Module[{R, q, w, r, Wex},
  {q, w} = If[wp === MachinePrecision, {N[a], N[om]}, {N[a, wp], N[om, wp]}];
  R = TeukolskyRadial[s, l, m, q, w, Method -> "MST"];
  r = N[1 + Sqrt[1 - a^2] + dr, wp];
  Wex = 2 I w R["In"]["Amplitudes"]["Incidence"] R["Up"]["Amplitudes"]["Transmission"];
  N[Abs[(r^2 - 2 r + q^2)^(s + 1) (R["In"][r] R["Up"]'[r] - R["In"]'[r] R["Up"][r])/Wex - 1]]
 ];

VerificationTest[
  {Table[wronskianErr[s, 2, 2, 3/5, 1/2, 1/100, 32] < 10^-28, {s, {-2, 0, 2}}], wronskianErr[-2, 2, -2, 3/5, -3/10 - I/5, 1/20, 32] < 10^-28, wronskianErr[-2, 2, 2, 3/5, 1/2, 1/100, MachinePrecision] < 10^-12},
  {{True, True, True}, True, True},
  TestID -> "Wronskian identity of the MST solutions next to the horizon, real and complex omega, 32 digits and machine precision"
]

VerificationTest[
  Module[{R, t, v},
    R = TeukolskyRadial[-2, 2, 2, N[3/5, 32], N[1/2, 32], Method -> "MST"];
    {t, v} = AbsoluteTiming[R["Up"][N[1 + 4/5 + 1/100, 32]]];
    {NumericQ[v], t < 10}
  ],
  {True, True},
  TestID -> "Up at r+ + 1/100 evaluates in seconds rather than minutes"
]

(* at the superradiant bound frequency the horizon amplitudes are Indeterminate and the Coulomb-type series is used *)
VerificationTest[
  Module[{R, r = N[1 + 4/5 + 3/10, 32], q = N[3/5, 32], w = N[1/2, 32], Wex},
    R = Quiet[TeukolskyRadial[0, 3, 3, q, w, Method -> "MST"], TeukolskyRadial::superradiant];
    Wex = 2 I w R["In"]["Amplitudes"]["Incidence"] R["Up"]["Amplitudes"]["Transmission"];
    {Teukolsky`MST`MST`Private`mstUpRepresentation[0, 3, q, 2 w, r], N[Abs[(r^2 - 2 r + q^2) (R["In"][r] R["Up"]'[r] - R["In"]'[r] R["Up"][r])/Wex - 1]] < 10^-28}
  ],
  {"Coulomb", True},
  TestID -> "The Coulomb-type series is kept at the superradiant bound frequency"
]
