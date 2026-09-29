(* Machine-precision radial functions: numerical integration from series boundary data (near-horizon power
   series for "In", large-r asymptotic series for "Up"), the "Up" solution of negative spin integrated at the
   flipped spin and mapped back with the Teukolsky-Starobinsky identity. References: MST at 32 digits. *)

seriesErr[s_, l_, m_, a_, om_, radii_] :=
 Module[{R, Rref},
  R = Quiet[TeukolskyRadial[s, l, m, N[a], N[om]]];
  Rref = Quiet[TeukolskyRadial[s, l, m, N[a, 32], N[om, 32], Method -> "MST"]];
  Max[Table[Abs[{R["In"][N[r]], R["In"]'[N[r]], R["Up"][N[r]], R["Up"]'[N[r]]}/{Rref["In"][N[r, 32]], Rref["In"]'[N[r, 32]], Rref["Up"][N[r, 32]], Rref["Up"]'[N[r, 32]]} - 1], {r, radii}]]
 ];

VerificationTest[
  Table[seriesErr[s, 2, 2, 3/5, 1/2, {1 + Sqrt[1 - 9/25] + 1/5, 3, 6, 20}] < 10^-11, {s, {-2, -1, 0, 1, 2}}],
  {True, True, True, True, True},
  TestID -> "Series boundary data: all spins at a = 3/5, omega = 1/2, to 1e-11 down to r+ + 1/5"
]

VerificationTest[
  {seriesErr[-2, 2, 2, 3/5, 1/100, {3, 20, 50}] < 10^-11, seriesErr[-2, 6, 6, 3/5, 1/100, {3, 20}] < 10^-11},
  {True, True},
  TestID -> "Series boundary data at omega = 1/100 (the large-r series starts near r = 1300)"
]

VerificationTest[
  seriesErr[-2, 2, 2, 99/100, 1/2, {1 + Sqrt[1 - 9801/10000] + 1/10, 3, 6}] < 10^-11,
  True,
  TestID -> "Series boundary data at a = 99/100 (the horizon series is evaluated at r+ + kappa)"
]

VerificationTest[
  seriesErr[-2, 2, 2, 3/5, 1/2 - I/5, {3, 6, 20}] < 10^-10,
  True,
  TestID -> "Series boundary data at a complex frequency"
]

(* the default method reports itself, and MST boundary data remain available *)
VerificationTest[
  Module[{R, RM},
    R = TeukolskyRadial[-2, 2, 2, 0.6, 0.5];
    RM = TeukolskyRadial[-2, 2, 2, 0.6, 0.5, Method -> {"NumericalIntegration", "BoundaryMethod" -> "MST"}];
    {R["Up"]["Method"], RM["Up"]["Method"], Abs[R["Up"][6.]/RM["Up"][6.] - 1] < 10^-11}
  ],
  {{"NumericalIntegration"}, {"NumericalIntegration", "BoundaryMethod" -> "MST"}, True},
  TestID -> "BoundaryMethod option"
]

(* the Teukolsky-Starobinsky map of the unit-transmission "Up" solution of spin +s is (2 I omega)^(2s) times
   that of spin -s: checked through the flipped-spin solution itself, beyond the series radius, where it is the
   map of the series *)
VerificationTest[
  Module[{R, Rref},
    R = TeukolskyRadial[-2, 2, 2, 0.6, 0.5];
    Rref = TeukolskyRadial[-2, 2, 2, N[3/5, 32], N[1/2, 32], Method -> "MST"];
    Abs[R["Up"][2000.]/Rref["Up"][N[2000, 32]] - 1] < 10^-12
  ],
  True,
  TestID -> "Flipped-spin Up solution beyond the series radius"
]

(* second derivatives through the radial equation, for the pieced-together Up function *)
VerificationTest[
  Module[{R, Rref},
    R = TeukolskyRadial[-2, 2, 2, 0.6, 0.5];
    Rref = TeukolskyRadial[-2, 2, 2, N[3/5, 32], N[1/2, 32], Method -> "MST"];
    {Abs[R["Up"][2.5, 2]/Rref["Up"][N[5/2, 32], 2] - 1] < 10^-10, Abs[R["Up"][6., 2]/Rref["Up"][N[6, 32], 2] - 1] < 10^-10, Abs[R["In"][6., 2]/Rref["In"][N[6, 32], 2] - 1] < 10^-10}
  ],
  {True, True, True},
  TestID -> "Second derivatives of the series-started solutions"
]

(* The "In" solution of negative spin at large radius: at omega = 2 the reflection is 1e-10 of the incidence, so
   the outward integration at spin -2 loses digits like r^4 (1e-8 at r = 50); beyond the radius where the large-r
   series converge (about 20 here) the solution is now Binc R_ingoing + Bref R_up from the series and the MST
   amplitudes *)
VerificationTest[
  Module[{R, Rref},
    R = TeukolskyRadial[-2, 2, 2, 0.6, 2.];
    Rref = TeukolskyRadial[-2, 2, 2, N[3/5, 32], N[2, 32], Method -> "MST"];
    Table[Abs[R["In"][N[r]]/Rref["In"][N[r, 32]] - 1] < 10^-12 && Abs[R["In"]'[N[r]]/Rref["In"]'[N[r, 32]] - 1] < 10^-12, {r, {5/2, 6, 50, 100, 300}}]
  ],
  {True, True, True, True, True},
  TestID -> "In solution of spin -2 at omega = 2 up to r = 300"
]
