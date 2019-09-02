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


