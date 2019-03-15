(* Mathematica Test File *)

\[Psi] = TeukolskyRadial[2, 2, 2, 0.5, 0.1];
\[Psi]In = \[Psi]["In"];

(****************************************************************)
(* TeukolskyRadial                                              *)
(****************************************************************)
VerificationTest[
    \[Psi]
    ,
    TeukolskyRadialFunction[2, 2, 2, 0.5, 0.1, 
     <|"Method" -> {"MST", "RenormalizedAngularMomentum" -> 1.9800007321197113}, 
      "BoundaryConditions" -> {"In", "Up"}, "Eigenvalue" -> -0.33267928615316333, 
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
    TeukolskyRadialFunction[2, 2, 2, 0.5, 0.1, 
      <|"Method" -> {"MST", "RenormalizedAngularMomentum" -> 1.9800007321197113}, 
       "BoundaryConditions" -> {"In", "Up"}, "Eigenvalue" -> -0.33267928615316333, 
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
    <|"In" -> 8.380259040239462 + 4.366632556855613*I, 
     "Up" -> 0.1172587510862727 - 0.04222308216793374*I|>
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
    <|"In" -> 0.22123657279829767 - 0.7299412890543474*I, 
     "Up" -> -0.05992621567089734 + 0.03634423319292853*I|>
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
    <|"In" -> -0.0028955067079097616 + 0.0010988209544549633*I, 
     "Up" -> 0.021151996439270782 - 0.02462511502507206*I|>
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
    TeukolskyRadialFunction[2, 2, 2, 0.5, 0.1, 
     <|"Method" -> {"MST", "RenormalizedAngularMomentum" -> 1.9800007321197113}, 
      "BoundaryConditions" -> "In", "Eigenvalue" -> -0.33267928615316333, 
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
    8.380259040239462 + 4.366632556855613*I
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
    0.22123657279829767 - 0.7299412890543474*I
    ,
    TestID->"Single Subcase Derivative Numerical Evaluation ",
    SameTest -> withinRoundoff
]


