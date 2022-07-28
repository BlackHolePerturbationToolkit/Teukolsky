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
        <|"Incidence" -> _,"Transmission" -> _, "Reflection" -> _|>,
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


(****************************************************************)
(* Static modes                                                 *)
(****************************************************************)
VerificationTest[
    TeukolskyRadial[2, 2, 2, 1/3, 0]["In"][10]
    ,
    ((514105 + 77004 I) 2^(-(7/2) - I/Sqrt[2] + 1/2 (2 - I/(2 Sqrt[2])))
        3^(3 + I/(2 Sqrt[2])) (9 - (2 Sqrt[2])/3)^(-2 + I/(
        2 Sqrt[2])) (1 - (3 (-9 + (2 Sqrt[2])/3))/(4 Sqrt[2]))^(-(I/(
        2 Sqrt[2]))) E^(
       1/3 I (1 + (2 Log[(2 Sqrt[2])/3])/(1 + (2 Sqrt[2])/3)))
        Gamma[-1 + I/Sqrt[2]])/((-2 I + Sqrt[2]) (2 I + Sqrt[2]) (-4 I + 
         Sqrt[2]) (4 + 27 Sqrt[2])^2 Gamma[1/2 I (2 I + Sqrt[2])])
    ,
    TestID->"Static \"in\" mode"
]

VerificationTest[
    TeukolskyRadial[2, 2, 2, 1/3, 0.0]["In"][10]
    ,
    ((514105 + 77004 I) 2^(-(7/2) - I/Sqrt[2] + 1/2 (2 - I/(2 Sqrt[2])))
        3^(3 + I/(2 Sqrt[2])) (9 - (2 Sqrt[2])/3)^(-2 + I/(
        2 Sqrt[2])) (1 - (3 (-9 + (2 Sqrt[2])/3))/(4 Sqrt[2]))^(-(I/(
        2 Sqrt[2]))) E^(
       1/3 I (1 + (2 Log[(2 Sqrt[2])/3])/(1 + (2 Sqrt[2])/3)))
        Gamma[-1 + I/Sqrt[2]])/((-2 I + Sqrt[2]) (2 I + Sqrt[2]) (-4 I + 
         Sqrt[2]) (4 + 27 Sqrt[2])^2 Gamma[1/2 I (2 I + Sqrt[2])])
    ,
    TestID->"Static \"in\" mode (machine-precision omega)"
]

VerificationTest[
    TeukolskyRadial[2, 2, 2, 1/3, 0.0``32]["In"][10]
    ,
    ((514105 + 77004 I) 2^(-(7/2) - I/Sqrt[2] + 1/2 (2 - I/(2 Sqrt[2])))
        3^(3 + I/(2 Sqrt[2])) (9 - (2 Sqrt[2])/3)^(-2 + I/(
        2 Sqrt[2])) (1 - (3 (-9 + (2 Sqrt[2])/3))/(4 Sqrt[2]))^(-(I/(
        2 Sqrt[2]))) E^(
       1/3 I (1 + (2 Log[(2 Sqrt[2])/3])/(1 + (2 Sqrt[2])/3)))
        Gamma[-1 + I/Sqrt[2]])/((-2 I + Sqrt[2]) (2 I + Sqrt[2]) (-4 I + 
         Sqrt[2]) (4 + 27 Sqrt[2])^2 Gamma[1/2 I (2 I + Sqrt[2])])
    ,
    TestID->"Static \"in\" mode (high-accuracy omega)"
]

VerificationTest[
    TeukolskyRadial[2, 2, 0, 1/3, 0]["In"][10]
    ,
    (243*(-9 + (2*Sqrt[2])/3)^2)/(256*(9 - (2*Sqrt[2])/3)^2)
    ,
    TestID->"Static m=0 \"in\" mode"
]

VerificationTest[
    TeukolskyRadial[2, 2, 0, 1/3, 0.0]["In"][10]
    ,
    (243*(-9 + (2*Sqrt[2])/3)^2)/(256*(9 - (2*Sqrt[2])/3)^2)
    ,
    TestID->"Static m=0 \"in\" mode (machine-precision omega)"
]

VerificationTest[
    TeukolskyRadial[2, 2, 0, 1/3, 0.0``32]["In"][10]
    ,
    (243*(-9 + (2*Sqrt[2])/3)^2)/(256*(9 - (2*Sqrt[2])/3)^2)
    ,
    TestID->"Static m=0 \"in\" mode (high-accuracy omega)"
]

VerificationTest[
    TeukolskyRadial[2, 2, 2, 1/3, 0]["Up"][10]
    ,
    (5 2^(-14 - I/Sqrt[2] + 1/2 (2 - I/(2 Sqrt[2]))) 3^(
       4 + I/(2 Sqrt[2])) (9 - (2 Sqrt[2])/3)^(-2 + I/(
        2 Sqrt[2])) (-4 + 27 Sqrt[2])^(-(I/Sqrt[
        2])) (4 + 27 Sqrt[2]) (1 - (3 (-9 + (2 Sqrt[2])/3))/(
         4 Sqrt[2]))^(-3 - I/(
        2 Sqrt[2])) ((514105 + 77004 I) (-4 + 27 Sqrt[2])^(I/Sqrt[2]) - 
         519841 (4 + 27 Sqrt[2])^(I/Sqrt[2])))/((-2 - I Sqrt[2]) (2 I + 
         Sqrt[2]) (-4 I + Sqrt[2]) (4 I + Sqrt[2]))
    ,
    TestID->"Static \"up\" mode"
]

VerificationTest[
    TeukolskyRadial[2, 2, 2, 1/3, 0.0]["Up"][10]
    ,
    (5 2^(-14 - I/Sqrt[2] + 1/2 (2 - I/(2 Sqrt[2]))) 3^(
       4 + I/(2 Sqrt[2])) (9 - (2 Sqrt[2])/3)^(-2 + I/(
        2 Sqrt[2])) (-4 + 27 Sqrt[2])^(-(I/Sqrt[
        2])) (4 + 27 Sqrt[2]) (1 - (3 (-9 + (2 Sqrt[2])/3))/(
         4 Sqrt[2]))^(-3 - I/(
        2 Sqrt[2])) ((514105 + 77004 I) (-4 + 27 Sqrt[2])^(I/Sqrt[2]) - 
         519841 (4 + 27 Sqrt[2])^(I/Sqrt[2])))/((-2 - I Sqrt[2]) (2 I + 
         Sqrt[2]) (-4 I + Sqrt[2]) (4 I + Sqrt[2]))
    ,
    TestID->"Static \"up\" mode (machine-precision omega)"
]

VerificationTest[
    TeukolskyRadial[2, 2, 2, 1/3, 0.0``32]["Up"][10]
    ,
    (5 2^(-14 - I/Sqrt[2] + 1/2 (2 - I/(2 Sqrt[2]))) 3^(
       4 + I/(2 Sqrt[2])) (9 - (2 Sqrt[2])/3)^(-2 + I/(
        2 Sqrt[2])) (-4 + 27 Sqrt[2])^(-(I/Sqrt[
        2])) (4 + 27 Sqrt[2]) (1 - (3 (-9 + (2 Sqrt[2])/3))/(
         4 Sqrt[2]))^(-3 - I/(
        2 Sqrt[2])) ((514105 + 77004 I) (-4 + 27 Sqrt[2])^(I/Sqrt[2]) - 
         519841 (4 + 27 Sqrt[2])^(I/Sqrt[2])))/((-2 - I Sqrt[2]) (2 I + 
         Sqrt[2]) (-4 I + Sqrt[2]) (4 I + Sqrt[2]))
    ,
    TestID->"Static \"up\" mode (high-accuracy omega)"
]


VerificationTest[
    TeukolskyRadial[2, 2, 0, 1/3, 0]["Up"][10]
    ,
    (-405*(4 + 27*Sqrt[2])*(77292*Sqrt[2] + 519841*Log[1 - (1 - (3*(-9 + (2*Sqrt[2])/3))/(4*Sqrt[2]))^(-1)]))/(524288*Sqrt[2]*(9 - (2*Sqrt[2])/3)^2*(1 - (3*(-9 + (2*Sqrt[2])/3))/(4*Sqrt[2]))^3)
    ,
    TestID->"Static m=0 \"up\" mode"
]

VerificationTest[
    TeukolskyRadial[2, 2, 0, 1/3, 0.0]["Up"][10]
    ,
    (-405*(4 + 27*Sqrt[2])*(77292*Sqrt[2] + 519841*Log[1 - (1 - (3*(-9 + (2*Sqrt[2])/3))/(4*Sqrt[2]))^(-1)]))/(524288*Sqrt[2]*(9 - (2*Sqrt[2])/3)^2*(1 - (3*(-9 + (2*Sqrt[2])/3))/(4*Sqrt[2]))^3)
    ,
    TestID->"Static m=0 \"up\" mode (machine-precision omega)"
]

VerificationTest[
    TeukolskyRadial[2, 2, 0, 1/3, 0.0``32]["Up"][10]
    ,
   (-405*(4 + 27*Sqrt[2])*(77292*Sqrt[2] + 519841*Log[1 - (1 - (3*(-9 + (2*Sqrt[2])/3))/(4*Sqrt[2]))^(-1)]))/(524288*Sqrt[2]*(9 - (2*Sqrt[2])/3)^2*(1 - (3*(-9 + (2*Sqrt[2])/3))/(4*Sqrt[2]))^3)
    ,
    TestID->"Static m=0 \"up\" mode (high-accuracy omega)"
]

(****************************************************************)
(* Fluxes                                                       *)
(****************************************************************)
VerificationTest[
    Module[{orbit = KerrGeoOrbit[0.9`32, 10.0`32, 0, 1]},
      TeukolskyPointParticleMode[-2, 2, 2, 0, 0, orbit]["Fluxes"]
    ]
    ,
    <|"Energy" -> <|"\[ScriptCapitalI]" -> 0.0000222730005511804972562,
       "\[ScriptCapitalH]" -> -5.9836792133833632312*10^-8|>,
     "AngularMomentum" -> <|"\[ScriptCapitalI]" ->
        0.00072437982117522330809,
       "\[ScriptCapitalH]" -> -1.94605862313006131615*10^-6|>|>
    ,
    TestID->"Fluxes"
]

VerificationTest[
    Module[{orbit = KerrGeoOrbit[0.1`32, 10.0`32, 0, Cos[\[Pi]/4.0`32]]},
      TeukolskyPointParticleMode[0, 2, 2, 0, 2, orbit]["Amplitudes"]
    ]
    ,
    <|"\[ScriptCapitalI]" -> 
      2.817674016780767525*10^-6 - 1.506523056130186854*10^-6 I, 
     "\[ScriptCapitalH]" -> -9.033598925366890723*10^-9 + 
       2.194950369301104043*10^-9 I|>
    ,
    TestID->"Spherical orbit amplitudes (s=0)"
]

VerificationTest[
    Module[{orbit = KerrGeoOrbit[0.1`32, 10.0`32, 0, Cos[\[Pi]/4.0`32]]},
      TeukolskyPointParticleMode[-2, 2, 2, 0, 2, orbit]["Amplitudes"]
    ]
    ,<|"ℐ" -> , "ℋ" -> |>
    <|"\[ScriptCapitalI]" -> -1.830922692793925658154841194038`19.860474179811195*10^-7 + 6.205049783176476429478414886579`19.400330311612727*10^-8 I, 
     "\[ScriptCapitalH]" -> -1.0285026158488625877011240877191`19.836131477231444*10^-7 - 5.128302017372319782302308898876`19.540942944235*10^-8 I|>
    ,
    TestID->"Spherical orbit amplitudes (s=-2)"
]
