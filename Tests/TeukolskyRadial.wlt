(* Mathematica Test File *)

{\[Psi]In, \[Psi]Up} = Values[TeukolskyRadial[2, 2, 2, 0.5, 0.1]];

(****************************************************************)
(* TeukolskyRadial                                              *)
(****************************************************************)
VerificationTest[
    \[Psi]In
    ,
    TeukolskyRadialFunction[2, 2, 2, 0.5, 0.1, <|
      "s" -> 2, "l" -> 2, "m" -> 2, "a" -> 0.5, "\[Omega]" -> 0.1, "Eigenvalue" -> \[Lambda]_,
      "Method" -> {"MST", "RenormalizedAngularMomentum" -> \[Nu]_}, 
      "BoundaryConditions" -> "In", 
      "Amplitudes" ->
        <|(* "Incidence" -> _,  *)"Transmission" -> _(* , "Reflection" -> _  *)|>,
        "Domain" -> {_, Infinity},
      "RadialFunction" -> Teukolsky`MST`MST`Private`MSTRadialIn[2, 2, 2, 0.5, 0.2, 
         \[Nu]_, \[Lambda]_, _, {MachinePrecision, MachinePrecision/2, MachinePrecision/2}]|>
    ]
    ,
    TestID->"TeukolskyRadial",
    SameTest -> MatchQ
]

(****************************************************************)
(* InvalidKey                                                   *)
(****************************************************************)
VerificationTest[
    \[Psi]In["NotAKey"]
    ,
    Missing["KeyAbsent", "NotAKey"]
    ,
    TestID->"InvalidKey"
]


(****************************************************************)
(* BoundaryConditions                                           *)
(****************************************************************)
VerificationTest[
    \[Psi]In["BoundaryConditions"]
    ,
    "In"
    ,
    TestID->"BoundaryConditions"
]


(****************************************************************)
(* Method                                                       *)
(****************************************************************)
VerificationTest[
    \[Psi]In["Method"]
    ,
    {"MST", "RenormalizedAngularMomentum" -> _}
    ,
    TestID->"Method",
    SameTest -> MatchQ
]


(****************************************************************)
(* Eigenvalue                                                   *)
(****************************************************************)
VerificationTest[
    \[Psi]In["Eigenvalue"]
    ,
    -0.33267928615316333
    ,
    TestID->"Eigenvalue"
]


(****************************************************************)
(* SolutionFunctions                                            *)
(****************************************************************)
VerificationTest[
    \[Psi]In["RadialFunction"]
    ,
    TeukolskyRadialFunction[__]["RadialFunction"]
    ,
    TestID->"RadialFunction",
    SameTest -> MatchQ
]


(****************************************************************)
(* Numerical Evaluation                                         *)
(****************************************************************)
VerificationTest[
    \[Psi]In[10.0]
    ,
    0.8151274441805518 + 0.5569358337070693*I
    ,
    TestID->"Numerical Evaluation",
    SameTest -> withinRoundoff
]


(****************************************************************)
(* Derivative Numerical Evaluation                              *)
(****************************************************************)
VerificationTest[
    \[Psi]In'[10.0]
    ,
    0.032005837243917735 - 0.07297367925001329*I
    ,
    TestID->"Derivative Numerical Evaluation ",
    SameTest -> withinRoundoff
]


(****************************************************************)
(* Higher Derivative Numerical Evaluation                       *)
(****************************************************************)
VerificationTest[
    \[Psi]In''''[10.0]
    ,
    -0.0003139905476106219 + 0.00007805832804581858*I
    ,
    TestID->"Higher Derivative Numerical Evaluation ",
    SameTest -> withinRoundoff
]



