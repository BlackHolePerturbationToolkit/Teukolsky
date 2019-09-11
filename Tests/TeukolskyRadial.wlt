(* Mathematica Test File *)

\[Psi] = TeukolskyRadial[2, 2, 2, 0.5, 0.1];
\[Psi]In = \[Psi]["In"];

(****************************************************************)
(* TeukolskyRadial                                              *)
(****************************************************************)
VerificationTest[
    \[Psi]
    ,
    TeukolskyRadialFunction[2, 2, 2, 0.5, 0.1, <|"s" -> 2, "l" -> 2, "m" -> 2, "\[Omega]" -> 0.1,
      "Method" -> {"MST", "RenormalizedAngularMomentum" -> 1.9800007321197113}, 
      "BoundaryConditions" -> {"In", "Up"}, "Eigenvalue" -> -0.33267928615316333, 
      "Amplitudes" ->
        <|"In" -> <|"Incidence" -> -188.00403765671575 + 1.8501269745181474*I,
                    "Transmission" -> 9.504222894583032 - 1.136766701917988*I,
                    "Reflection" -> 2.0202392307641034*^6 - 1.7016381012327243*^6*I|>,
          "Up" -> <|"Transmission" -> -1805.5571508475773 - 9034.631436940852*I|>|>, 
      "SolutionFunctions" -> {Teukolsky`MST`Private`MSTRadialIn[2, 2, 2, 0.5, 0.2, 
         1.9800007321197113, -0.33267928615316333], Teukolsky`MST`Private`MSTRadialUp[2, 2, 
         2, 0.5, 0.2, 1.9800007321197113, -0.33267928615316333]}|>]
    ,
    TestID->"TeukolskyRadial",
    SameTest -> withinRoundoff
]

(****************************************************************)
(* InvalidKey                                                   *)
(****************************************************************)
VerificationTest[
    \[Psi]["NotAKey"]
    ,
    Missing["KeyAbsent", "NotAKey"]
    ,
    TestID->"InvalidKey"
]


(****************************************************************)
(* BoundaryConditions                                           *)
(****************************************************************)
VerificationTest[
    \[Psi]["BoundaryConditions"]
    ,
    {"In", "Up"}
    ,
    TestID->"BoundaryConditions"
]


(****************************************************************)
(* Method                                                       *)
(****************************************************************)
VerificationTest[
    \[Psi]["Method"]
    ,
    {"MST", "RenormalizedAngularMomentum" -> 1.9800007321197113}
    ,
    TestID->"Method"
]


(****************************************************************)
(* Eigenvalue                                                   *)
(****************************************************************)
VerificationTest[
    \[Psi]["Eigenvalue"]
    ,
    -0.33267928615316333
    ,
    TestID->"Eigenvalue"
]


(****************************************************************)
(* SolutionFunctions                                            *)
(****************************************************************)
VerificationTest[
    \[Psi]["SolutionFunctions"]
    ,
    TeukolskyRadialFunction[2, 2, 2, 0.5, 0.1, <|"s" -> 2, "l" -> 2, "m" -> 2, "\[Omega]" -> 0.1,
       "Method" -> {"MST", "RenormalizedAngularMomentum" -> 1.9800007321197113}, 
       "BoundaryConditions" -> {"In", "Up"}, "Eigenvalue" -> -0.33267928615316333, 
       "Amplitudes" ->
         <|"In" -> <|"Incidence" -> -188.00403765671575 + 1.8501269745181474*I,
                     "Transmission" -> 9.504222894583032 - 1.136766701917988*I,
                     "Reflection" -> 2.0202392307641034*^6 - 1.7016381012327243*^6*I|>,
           "Up" -> <|"Transmission" -> -1805.5571508475773 - 9034.631436940852*I|>|>, 
       "SolutionFunctions" -> {Teukolsky`MST`Private`MSTRadialIn[2, 2, 2, 0.5, 0.2, 
          1.9800007321197113, -0.33267928615316333], Teukolsky`MST`Private`MSTRadialUp[2, 2, 
          2, 0.5, 0.2, 1.9800007321197113, -0.33267928615316333]}|>]["SolutionFunctions"]
    ,
    TestID->"SolutionFunctions"
]


(****************************************************************)
(* Numerical Evaluation                                         *)
(****************************************************************)
VerificationTest[
    \[Psi][10.0]
    ,
    <|"In" -> 0.8151274455692021 + 0.5569358329985141*I,
      "Up" -> 1.999804503329915*^-6 + 0.000013378466321641552*I|>
    ,
    TestID->"Numerical Evaluation",
    SameTest -> withinRoundoff
]


(****************************************************************)
(* Derivative Numerical Evaluation                              *)
(****************************************************************)
VerificationTest[
    \[Psi]'[10.0]
    ,
    <|"In" -> 0.03200583839610331 - 0.07297367974124275*I,
      "Up" -> -2.593598129598713*^-6 - 7.151271833408181*^-6*I|>
    ,
    TestID->"Derivative Numerical Evaluation ",
    SameTest -> withinRoundoff
]


(****************************************************************)
(* Higher Derivative Numerical Evaluation                       *)
(****************************************************************)
VerificationTest[
    \[Psi]''''[10.0]
    ,
    <|"In" -> -0.0003139910712826748 + 0.00007805860280741137*I,
      "Up" -> 2.171038044040099*^-6 + 2.775091588344741*^-6*I|>
    ,
    TestID->"Higher Derivative Numerical Evaluation ",
    SameTest -> withinRoundoff
]


(****************************************************************)
(* Subcase                                                      *)
(****************************************************************)
VerificationTest[
    \[Psi]In
    ,
    TeukolskyRadialFunction[2, 2, 2, 0.5, 0.1, <|"s" -> 2, "l" -> 2, "m" -> 2, "\[Omega]" -> 0.1,
      "Method" -> {"MST", "RenormalizedAngularMomentum" -> 1.9800007321197113}, 
      "BoundaryConditions" -> "In", "Eigenvalue" -> -0.33267928615316333, 
      "Amplitudes" -> 
        <|"Incidence" -> -188.00403765671575 + 1.8501269745181474*I,
          "Transmission" -> 9.504222894583032 - 1.136766701917988*I,
          "Reflection" -> 2.0202392307641034*^6 - 1.7016381012327243*^6*I|>, 
      "SolutionFunctions" -> Teukolsky`MST`Private`MSTRadialIn[2, 2, 2, 0.5, 0.2, 
        1.9800007321197113, -0.33267928615316333]|>]
    ,
    TestID->"Subcase"
]

(****************************************************************)
(* Single Subcase Boundary Conditions                           *)
(****************************************************************)
VerificationTest[
    \[Psi]In["BoundaryConditions"]
    ,
    "In"
    ,
    TestID->"Single Subcase Boundary Conditions"
]

(****************************************************************)
(* Single Subcase Numerical Evaluation                          *)
(****************************************************************)
VerificationTest[
    \[Psi]In[10.0]
    ,
    0.8151274455692021 + 0.5569358329985141*I
    ,
    TestID->"Single Subcase Numerical Evaluation",
    SameTest -> withinRoundoff
]


(****************************************************************)
(* Single Subcase Derivative Numerical Evaluation               *)
(****************************************************************)
VerificationTest[
    \[Psi]In'[10.0]
    ,
    0.03200583839610331 - 0.07297367974124275*I
    ,
    TestID->"Single Subcase Derivative Numerical Evaluation ",
    SameTest -> withinRoundoff
]


