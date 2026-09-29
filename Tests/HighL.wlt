(* A mode deep under its potential barrier: l = m = 40 for the circular orbit at r0 = 9 Sqrt[11]/5 with a = 3/5.
   The MST series lose about 100 digits to cancellation and amplify an error in nu by about 27 digits, so the
   result is only right if nu is recomputed at the padded precision. The logarithmic derivative of the "In"
   solution at r0 is real (the mode sits under a 24-digit barrier); reference from a 320-digit evaluation. *)
VerificationTest[
  Module[{aa = 6/10, r0 = 9 Sqrt[11]/5, Om, R, v},
    Om = 1/(r0^(3/2) + aa);
    R = TeukolskyRadial[0, 40, 40, N[aa, 80], N[40 Om, 80]]["In"];
    v = R'[N[r0, 80]]/R[N[r0, 80]];
    {Abs[Re[v]/7.28778025417069399615869085427167465392 - 1] < 10^-30, Abs[Im[v]] < 10^-30, Precision[v] > 70}
  ],
  {True, True, True},
  TestID -> "High-l mode under the potential barrier at 80 digits"
]

(* A mode for which the recurrence for the MST coefficients silently yields the wrong solution below about
   300 digits: l = 36, m = 2, omega = 3 at a = 3/5. The Wronskian check catches it and raises the working
   precision (80 -> 160 -> 320; about 90 s). Reference from a checked evaluation (Wronskian identity to 1e-78). *)
VerificationTest[
  Module[{aa = 6/10, r0 = 9 Sqrt[11]/5, R, r, G},
    R = TeukolskyRadial[0, 36, 2, N[aa, 80], N[3, 80]];
    r = N[r0, 80];
    G = -R["In"][r] R["Up"][r]/(2 I 3 R["In"]["Amplitudes"]["Incidence"]);
    (* G has a small imaginary part (2e-14 relative), so the real part is compared *)
    Abs[Re[G]/0.0034919644133859045660219331138202791577813 - 1] < 10^-25
  ],
  True,
  TestID -> "Wronskian check raises the working precision for l = 36, m = 2, omega = 3"
]
