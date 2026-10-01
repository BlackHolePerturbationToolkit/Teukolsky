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

(* The same symmetry at machine precision with the default method, where the amplitudes come from the MST
   package and the functions from the numerical integration: the "In" reflection amplitude at negative real
   frequency was wrong by an order-one factor (the powers and logarithms of epsilon in the amplitude formulae
   on the wrong side of their cuts) until the MST package took the amplitudes at Re omega < 0 as the conjugates
   of those at (-m, -Conjugate[omega]) *)
mpSymmetryErr[s_, l_, m_, a_, om_] :=
 Module[{R, Rc, keys = {"Incidence", "Transmission", "Reflection"}, amps, fs},
  R = TeukolskyRadial[s, l, m, a, om];
  Rc = TeukolskyRadial[s, l, -m, a, -om];
  amps = Max[Table[Abs[Conjugate[R[bc]["Amplitudes"][k]]/Rc[bc]["Amplitudes"][k] - 1], {bc, {"In", "Up"}}, {k, keys}]];
  fs = Max[Table[Abs[Conjugate[R[bc][r]]/Rc[bc][r] - 1], {bc, {"In", "Up"}}, {r, {6., 20.}}]];
  {amps, fs}
 ];

VerificationTest[
  Max[mpSymmetryErr[-2, 2, 2, 0.6, 0.3]] < 10^-12 && Max[mpSymmetryErr[0, 2, 1, 0.6, 0.5]] < 10^-12,
  True,
  TestID -> "Amplitudes and functions at negative real omega are the conjugates of those at (-m, omega), machine precision"
]

(* The padding a failed Wronskian check sets for a mode must reach the series, which at Re omega < 0 are
   evaluated at the conjugate partner (-m, -Conjugate[omega]): the padding is keyed by the partner *)
VerificationTest[
  Module[{q = N[3/5, 32], set = Teukolsky`MST`MST`Private`setModePadding, get = Teukolsky`MST`MST`Private`modePadding},
    set[-2, 2, 2, q, N[-1, 32], 17];
    {get[-2, 2, -2, q, N[1, 32]], get[-2, 2, 2, q, N[-1, 32]]}
  ],
  {17, 17},
  TestID -> "Mode padding set at a negative frequency is read at the conjugate partner"
]

(* on the positive imaginary axis the incoming Coulomb-type series is evaluated on the principal side of the
   cut of its Tricomi functions (the Up series on the other side at the negative imaginary axis) *)
VerificationTest[
  {Max[coulombErr[0, 2, 0, 3/5, 7 I/10, {30, 60}]] < 10^-28, Max[coulombErr[-2, 2, 2, 3/5, 7 I/10, {30, 60}]] < 10^-28, wronskianErr[0, 2, 0, 3/5, 7 I/10, {5/2, 10, 60}] < 10^-28},
  {True, True, True},
  TestID -> "Coulomb-type In representation and Wronskian on the positive imaginary axis"
]

(* The first evaluation of the Coulomb-type "In" representation of a mode is a single pass that measures the
   loss; at machine precision its result was accepted on the precision of N[result], always MachinePrecision,
   although the pass had kept less than a digit (1e-3 at r = 50 for this mode, right on the second call) *)
VerificationTest[
  Module[{R = TeukolskyRadial[1, 8, 6, 0.99, 0.0027 - 0.856 I], R40 = TeukolskyRadial[1, 8, 6, N[99/100, 40], N[27/10000 - 856/1000 I, 40]]},
    Abs[R["In"][50.]/R40["In"][N[50, 40]] - 1] < 10^-12
  ],
  True,
  TestID -> "First machine-precision Coulomb-type In evaluation at a complex frequency"
]

(* a monodromy estimate without correct digits on the way gave an unevaluated nu and N::precbd *)
VerificationTest[
  Module[{R = TeukolskyRadial[0, 8, -5, 0.1, 0.592 - 0.915 I], R40 = TeukolskyRadial[0, 8, -5, N[1/10, 40], N[592/1000 - 915/1000 I, 40]]},
    Abs[R["In"][100.]/R40["In"][N[100, 40]] - 1] < 10^-12
  ],
  True,
  TestID -> "Large-radius In at a complex frequency without precision messages"
]

(* nu near an integer at a small complex frequency: the downward MST coefficients can converge to the wrong
   solution with a tracked precision that claims full accuracy (1e-3 here before), so the Wronskian check runs *)
VerificationTest[
  Module[{R = TeukolskyRadial[0, 8, 0, 0., -0.01 I], R40 = TeukolskyRadial[0, 8, 0, N[0, 40], N[-1/100 I, 40]]},
    Max[Abs[R["In"][#]/R40["In"][SetPrecision[#, 40]] - 1] & /@ {4., 10., 30.}, Abs[R["Up"][#]/R40["Up"][SetPrecision[#, 40]] - 1] & /@ {4., 10., 30.}] < 10^-12
  ],
  True,
  TestID -> "Near-integer nu at a small imaginary frequency"
]

(* where the retries at machine precision cannot repair such a mode, the result is flagged instead of silent *)
VerificationTest[
  AssociationQ[TeukolskyRadial[0, 12, 0, 0.5, 0.01 - 0.02 I]],
  True,
  {TeukolskyRadial::acc},
  TestID -> "Near-integer nu beyond the machine-precision retries is flagged"
]
