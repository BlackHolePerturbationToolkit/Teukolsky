(* MST solutions at complex frequencies: the conjugation symmetry R[s, l, m, a, omega] = Conjugate[R[s, l, -m, a, -Conjugate[omega]]]
   of the radial equation (violated by O(1) for Re omega < 0 before the Coulomb-type series were restricted to Re omega > 0),
   the Wronskian identity and the continuity of the "Up" solution on the negative imaginary axis (where the argument of
   the Tricomi functions lies on their branch cut), and the Coulomb-type representation of the "In" solution at large radius. *)

conjErr[s_, l_, m_, a_, om_, radii_] :=
 Module[{R, Rc, wp = 32, fs},
  R = TeukolskyRadial[s, l, m, N[a, wp], N[om, wp], Method -> "MST"];
  Rc = TeukolskyRadial[s, l, -m, N[a, wp], N[-Conjugate[om], wp], Method -> "MST"];
  fs = Max[Table[Abs[{R["In"][N[r, wp]], R["In"]'[N[r, wp]], R["Up"][N[r, wp]], R["Up"]'[N[r, wp]]}/Conjugate[{Rc["In"][N[r, wp]], Rc["In"]'[N[r, wp]], Rc["Up"][N[r, wp]], Rc["Up"]'[N[r, wp]]}] - 1], {r, radii}]];
  N[{fs, Abs[R["In"]["Amplitudes"]["Incidence"]/Conjugate[Rc["In"]["Amplitudes"]["Incidence"]] - 1], Abs[R["Up"]["Amplitudes"]["Reflection"]/Conjugate[Rc["Up"]["Amplitudes"]["Reflection"]] - 1]}]
 ];

VerificationTest[
  Max[conjErr[-2, 3, -3, 3/5, -4/5 - 9 I/10, {6, 12, 30}]] < 10^-28,
  True,
  TestID -> "Conjugation symmetry of the MST solutions and amplitudes at Re omega < 0, complex omega"
]

VerificationTest[
  {Max[conjErr[-2, 2, -2, 3/5, -3/2, {6, 12, 30}]] < 10^-28, Max[conjErr[2, 2, 0, 3/5, -1/2, {6, 12, 30}]] < 10^-28},
  {True, True},
  TestID -> "Conjugation symmetry at negative real omega (the Coulomb-type In representation beyond the switch radius)"
]

wronskianErr[s_, l_, m_, a_, om_, radii_] :=
 Module[{R, Wex, wp = 32, q, w},
  {q, w} = N[{a, om}, wp];
  R = TeukolskyRadial[s, l, m, q, w, Method -> "MST"];
  Wex = 2 I w R["In"]["Amplitudes"]["Incidence"] R["Up"]["Amplitudes"]["Transmission"];
  N[Max[Table[Abs[(r^2 - 2 r + q^2)^(s + 1) (R["In"][r] R["Up"]'[r] - R["In"]'[r] R["Up"][r])/Wex - 1], {r, N[radii, wp]}]]]
 ];

VerificationTest[
  Table[wronskianErr[s, 2, 2, 3/5, -3 I/10, {6, 12, 30}] < 10^-28, {s, {-2, 0, 2}}],
  {True, True, True},
  TestID -> "Wronskian identity on the negative imaginary axis (Tricomi functions on their branch cut)"
]

VerificationTest[
  Module[{wp = 32, q, h = 10^-4, r = 12, up, upm, upp},
    q = N[3/5, wp];
    up = TeukolskyRadial[-2, 2, 2, q, N[-3 I/10, wp], Method -> "MST"]["Up"];
    upm = TeukolskyRadial[-2, 2, 2, q, N[-h - 3 I/10, wp], Method -> "MST"]["Up"];
    upp = TeukolskyRadial[-2, 2, 2, q, N[h - 3 I/10, wp], Method -> "MST"]["Up"];
    N[Abs[(upm[N[r, wp]] + upp[N[r, wp]])/(2 up[N[r, wp]]) - 1]] < 10^-4 && N[Abs[(upm'[N[r, wp]] + upp'[N[r, wp]])/(2 up'[N[r, wp]]) - 1]] < 10^-4
  ],
  True,
  TestID -> "Continuity of the Up solution across Re omega = 0"
]

coulombErr[s_, l_, m_, a_, om_, radii_] :=
 Module[{R, wp = 32, def, ser},
  R = TeukolskyRadial[s, l, m, N[a, wp], N[om, wp], Method -> "MST"];
  def = R["In"][N[radii, wp]];
  ser = Block[{Teukolsky`MST`MST`Private`$MSTRepresentationThreshold = Infinity}, R["In"][N[radii, wp]]];
  N[Max[Abs[def/ser - 1]]]
 ];

VerificationTest[
  {Teukolsky`MST`MST`Private`mstInRepresentation[N[3/5, 32], N[2 (1 - I/2), 32], N[30, 32]],
   coulombErr[-2, 2, 2, 3/5, 1 - I/2, {30, 60}] < 10^-28, coulombErr[0, 2, 0, 3/5, -7 I/10, {30, 60}] < 10^-28, coulombErr[-2, 2, -2, 3/5, -1 - I/2, {30, 60}] < 10^-28},
  {"Coulomb", True, True, True},
  TestID -> "Coulomb-type In representation at large radius agrees with the hypergeometric series at complex omega"
]

VerificationTest[
  Module[{R32, RMP},
    R32 = TeukolskyRadial[-2, 2, 2, N[3/5, 32], N[3/2 - 2 I, 32], Method -> "MST", "WronskianCheck" -> True];
    RMP = TeukolskyRadial[-2, 2, 2, 0.6, 1.5 - 2. I, "WronskianCheck" -> True];
    {NumericQ[R32["In"][N[6, 32]]], NumericQ[RMP["Up"][6.]], Abs[RMP["Up"][6.]/R32["Up"][N[6, 32]] - 1] < 10^-12}
  ],
  {True, True, True},
  TestID -> "Wronskian check at complex omega passes without messages"
]
