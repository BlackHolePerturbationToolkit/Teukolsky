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
(* NearZeroFrequency                                            *)
(****************************************************************)
VerificationTest[
    RenormalizedAngularMomentum[-2, 2, 0, 0, 0.00001, 4,  Method -> "FindRoot"]
    ,
    1.9999999999490476
    ,
    TestID->"NearZeroFrequencyFindRoot"
]



(****************************************************************)
(* NearZeroFrequencyFindRootWithGuess                           *)
(****************************************************************)
VerificationTest[
    RenormalizedAngularMomentum[-2, 2, 0, 0, 0.00001, 4,  Method -> {"FindRoot", "InitialGuess" -> l (l - 1)}]
    ,
    1.9999999999490476
    ,
    TestID->"NearZeroFrequencyFindRootWithGuess"
]


