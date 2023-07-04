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
      "RenormalizedAngularMomentum" -> \[Nu]_, 
      "Method" -> {"NumericalIntegration"}, 
      "BoundaryConditions" -> "In", 
      "Amplitudes" ->
        <|"Incidence" -> _,"Transmission" -> _, "Reflection" -> _|>,
      "UnscaledAmplitudes" ->
        <|"Incidence" -> _,"Transmission" -> _, "Reflection" -> _|>,
      "Domain" -> {_, Infinity},
      "RadialFunction" -> _Function|>
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
    {"NumericalIntegration"}
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
    SetPrecision[\[Psi]In[10.0], 6]
    ,
    0.8151274455692312 + 0.5569358329985331*I
    ,
    TestID->"Numerical Evaluation",
    SameTest -> withinRoundoff
]


(****************************************************************)
(* Derivative Numerical Evaluation                              *)
(****************************************************************)
VerificationTest[
    SetPrecision[\[Psi]In'[10.0], 6]
    ,
    0.032005838396104935 - 0.07297367974124509*I
    ,
    TestID->"Derivative Numerical Evaluation ",
    SameTest -> withinRoundoff
]


(****************************************************************)
(* Higher Derivative Numerical Evaluation                       *)
(****************************************************************)
VerificationTest[
    SetPrecision[\[Psi]In''''[10.0], 2]
    ,
    -0.00031399107128295535 + 0.00007805860280733666*I
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
      N[TeukolskyPointParticleMode[-2, 2, 2, 0, 0, orbit]["Fluxes"]]
    ]
    ,
    <|"Energy" -> <|"\[ScriptCapitalI]" -> 0.0000222730005511805,
       "\[ScriptCapitalH]" -> -5.983679213383364*^-8|>,
     "AngularMomentum" -> <|"\[ScriptCapitalI]" ->
        0.0007243798211752233,
       "\[ScriptCapitalH]" -> -1.9460586231300613*^-6|>|>
    ,
    TestID->"Fluxes",
    SameTest->withinRoundoff
]

(****************************************************************)
(* Amplitudes                                                   *)
(****************************************************************)
Module[{orbitC = KerrGeoOrbit[0, 10.0, 0, 1]},
  VerificationTest[
      SetPrecision[TeukolskyPointParticleMode[0, 2, 2, 0, 0, orbitC]["Amplitudes"], 6]
    ,
    <|"\[ScriptCapitalI]" -> -0.096869111783302688482248 + 0.034691088866233852517826 I, 
     "\[ScriptCapitalH]" -> 0.00085249088747438226386116 - 0.00021231389467651976331909 I|>
    ,
    TestID->"Circular orbit Schwarzschild amplitudes (s=0)",
    SameTest->withinRoundoff
  ];

  VerificationTest[
      SetPrecision[TeukolskyPointParticleMode[-1, 2, 2, 0, 0, orbitC]["Amplitudes"], 6]
    ,
    <|"\[ScriptCapitalI]" -> 0.0019080193702080048309563 + 0.0055230606668409233703376 I, 
     "\[ScriptCapitalH]" -> 0.00109116926835821214306336 + 0.00002231589868836575762256 I|>
    ,
    TestID->"Circular orbit Schwarzschild amplitudes (s=-1)",
    SameTest->withinRoundoff
  ];

  VerificationTest[
      SetPrecision[TeukolskyPointParticleMode[1, 2, 2, 0, 0, orbitC]["Amplitudes"], 6]
    ,
    <|"\[ScriptCapitalI]" -> -1.43101452765597547942523 - 4.14229550013065129411783 I, 
     "\[ScriptCapitalH]" -> -0.000085585743658428712077571 - 0.000369966178105824541428495 I|>
    ,
    TestID->"Circular orbit Schwarzschild amplitudes (s=1)",
    SameTest->withinRoundoff
  ];

  VerificationTest[
      SetPrecision[TeukolskyPointParticleMode[-2, 2, 2, 0, 0, orbitC]["Amplitudes"], 6]
    ,
    <|"\[ScriptCapitalI]" -> -0.0011206407919254166103293 + 0.0003057608384581628431873 I, 
     "\[ScriptCapitalH]" -> -0.0014710247263310518147472 - 0.0003861305299246959983883 I|>
    ,
    TestID->"Circular orbit Schwarzschild amplitudes (s=-2)",
    SameTest->withinRoundoff
  ];
]

Module[{orbitS = KerrGeoOrbit[0, 10.0, 0, Cos[\[Pi]/4.]]},
  VerificationTest[
      SetPrecision[TeukolskyPointParticleMode[0, 2, 2, 0, 0, orbitS]["Amplitudes"], 6]
    ,
    <|"\[ScriptCapitalI]" -> -0.070574319833484021478549 + 0.025274310416867242591122 I, 
     "\[ScriptCapitalH]" -> 0.00062108512651933021254542 - 0.00015468200783664169003028 I|>
    ,
    TestID->"Spherical orbit Schwarzschild amplitudes (s=0)",
    SameTest->withinRoundoff
  ];

  VerificationTest[
      SetPrecision[TeukolskyPointParticleMode[-2, 2, 2, 0, 0, orbitS]["Amplitudes"], 6]
    ,
    <|"\[ScriptCapitalI]" -> -0.000816446648594393652746243 + 0.0002227630955693367815631962 I, 
     "\[ScriptCapitalH]" -> -0.001071720052015030398567520 - 0.000281316706788214826313457 I|>
    ,
    TestID->"Spherical orbit Schwarzschild amplitudes (s=-2)",
    SameTest->withinRoundoff
  ];
]

Module[{orbitE = KerrGeoOrbit[0, 10.0, 0.1, 1]},
  VerificationTest[
      SetPrecision[TeukolskyPointParticleMode[0, 2, 2, 0, 0, orbitE]["Amplitudes"], 6]
    ,
    <|"\[ScriptCapitalI]" -> -0.090966715053256986219528369 + 0.032357265869968441307292536 I, 
     "\[ScriptCapitalH]" -> 0.00078064421192569808851129185 - 0.00019244803540255800765595629 I|>
    ,
    TestID->"Eccentric orbit Schwarzschild amplitudes (s=0)",
    SameTest->withinRoundoff
  ];

  VerificationTest[
      SetPrecision[TeukolskyPointParticleMode[-2, 2, 2, 0, 0, orbitE]["Amplitudes"], 6]
    ,
    <|"\[ScriptCapitalI]" -> -0.00102739509088208 + 0.00027888266315215 I, 
     "\[ScriptCapitalH]" -> -0.00135207925401707 - 0.00035093350945876 I|>
    ,
    TestID->"Eccentric orbit Schwarzschild amplitudes (s=-2)",
    SameTest->withinRoundoff
  ];
]

Module[{orbitG = KerrGeoOrbit[0, 10.0, 0.1, Cos[\[Pi]/4.]]},
  VerificationTest[
      SetPrecision[TeukolskyPointParticleMode[0, 2, 2, 0, 0, orbitG]["Amplitudes"], 6]
    ,
    <|"\[ScriptCapitalI]" -> -0.066274108683182594817608793 + 0.023573995759893462190169796 I, 
     "\[ScriptCapitalH]" -> 0.00056874098744548149271895618 - 0.00014020866870554800968163363 I|>
    ,
    TestID->"Generic orbit Schwarzschild amplitudes (s=0)",
    SameTest->withinRoundoff
  ];

  VerificationTest[
      SetPrecision[TeukolskyPointParticleMode[-2, 2, 2, 0, 0, orbitG]["Amplitudes"], 6]
    ,
    <|"\[ScriptCapitalI]" -> -0.00074851217701405 + 0.00020318090971475 I, 
     "\[ScriptCapitalH]" -> -0.00098506192486100 - 0.00025567379818942 I|>
    ,
    TestID->"Generic orbit Schwarzschild amplitudes (s=-2)",
    SameTest->withinRoundoff
  ];
]


Module[{orbitC = KerrGeoOrbit[0.1, 10.0, 0, 1]},
  VerificationTest[
      SetPrecision[TeukolskyPointParticleMode[0, 2, 2, 0, 0, orbitC]["Amplitudes"], 6]
    ,
    <|"\[ScriptCapitalI]" -> -0.096227584894766381938308 + 0.034422505898013111362993 I, 
     "\[ScriptCapitalH]" -> 0.00082969529874874947045813 - 0.00004165341379611482771148 I|>
    ,
    TestID->"Circular orbit Kerr amplitudes (s=0)",
    SameTest->withinRoundoff
  ];

  VerificationTest[
      SetPrecision[TeukolskyPointParticleMode[-1, 2, 2, 0, 0, orbitC]["Amplitudes"], 6]
    ,
    <|"\[ScriptCapitalI]" -> 0.0018866657836778435088676 + 0.0054668655210736447633660 I, 
     "\[ScriptCapitalH]" -> 0.00106346002166190597056679 + 0.00001947401101465360379365 I|>
    ,
    TestID->"Circular orbit Kerr amplitudes (s=-1)",
    SameTest->withinRoundoff
  ];

  VerificationTest[
      SetPrecision[TeukolskyPointParticleMode[1, 2, 2, 0, 0, orbitC]["Amplitudes"], 6]
    ,
    <|"\[ScriptCapitalI]" -> -1.4179854870221155479558 - 4.1088019062243623376054 I, 
     "\[ScriptCapitalH]" -> -2.446845307003628488731*10^-6 - 0.000073108439893058193988156 I|>
    ,
    TestID->"Circular orbit Kerr amplitudes (s=1)",
    SameTest->withinRoundoff
  ];

  VerificationTest[
      SetPrecision[TeukolskyPointParticleMode[-2, 2, 2, 0, 0, orbitC]["Amplitudes"], 6]
    ,
    <|"\[ScriptCapitalI]" -> -0.0011042513597484404673410 + 0.0003010513879342753094272 I, 
     "\[ScriptCapitalH]" -> -0.0014748445457123321571078 - 0.0002289511157707928466590 I|>
    ,
    TestID->"Circular orbit Kerr amplitudes (s=-2)",
    SameTest->withinRoundoff
  ];
]

Module[{orbitS = KerrGeoOrbit[0.1, 10.0, 0, Cos[\[Pi]/4.]]},
  VerificationTest[
      SetPrecision[TeukolskyPointParticleMode[0, 2, 2, 0, 0, orbitS]["Amplitudes"], 6]
    ,
    <|"\[ScriptCapitalI]" -> -0.070468528422932792881 + 0.025252522903963179701 I, 
     "\[ScriptCapitalH]" -> 0.00060463777485046264504 - 0.00003075360976614048730 I|>
    ,
    TestID->"Spherical orbit Kerr amplitudes (s=0)",
    SameTest->withinRoundoff
  ];

  VerificationTest[
      SetPrecision[TeukolskyPointParticleMode[-2, 2, 2, 0, 0, orbitS]["Amplitudes"], 6]
    ,
    <|"\[ScriptCapitalI]" -> -0.0008134090552314842569669 + 0.0002220559776171971048038 I, 
     "\[ScriptCapitalH]" -> -0.0010818095010835464640668 - 0.0001687703221397815132785 I|>
    ,
    TestID->"Spherical orbit Kerr amplitudes (s=-2)",
    SameTest->withinRoundoff
  ];

  VerificationTest[
      SetPrecision[TeukolskyPointParticleMode[0, 2, 2, 0, 2, orbitS]["Amplitudes"], 6]
      ,
      <|"\[ScriptCapitalI]" -> 2.817674016780767525*10^-6 - 1.506523056130186854*10^-6 I,
       "\[ScriptCapitalH]" -> -9.033598925366890723*10^-9 + 2.194950369301104043*10^-9 I|>
      ,
      TestID->"Spherical orbit amplitudes (s=0)",
      SameTest->withinRoundoff
  ];

  VerificationTest[
      SetPrecision[TeukolskyPointParticleMode[-2, 2, 2, 0, 2, orbitS]["Amplitudes"], 6]
      ,
      <|"\[ScriptCapitalI]" -> 1.8309226927939255*^-7 - 6.205049783176476*^-8*I,
        "\[ScriptCapitalH]" -> 1.0285026158488626*^-7 + 5.12830201737232*^-8*I|>
      ,
      TestID->"Spherical orbit amplitudes (s=-2)",
      SameTest->withinRoundoff
  ];
]

Module[{orbitE = KerrGeoOrbit[0.1, 10.0, 0.1, 1]},
  VerificationTest[
      SetPrecision[TeukolskyPointParticleMode[0, 2, 2, 0, 0, orbitE]["Amplitudes"], 6]
    ,
    <|"\[ScriptCapitalI]" -> -0.090689967567633577108184261 + 0.032215375076956083066057133 I, 
     "\[ScriptCapitalH]" -> 0.00076250722962311750779893683 - 0.00003624478217354641654074258 I|>
    ,
    TestID->"Eccentric orbit Kerr amplitudes (s=0)",
    SameTest->withinRoundoff
  ];

  VerificationTest[
      SetPrecision[TeukolskyPointParticleMode[-2, 2, 2, 0, 0, orbitE]["Amplitudes"], 6]
    ,
    <|"\[ScriptCapitalI]" -> -0.00101539073433068 + 0.00027535746450647 I, 
     "\[ScriptCapitalH]" -> -0.001359631602157010 - 0.000207192301964262 I|>
    ,
    TestID->"Eccentric orbit Kerr amplitudes (s=-2)",
    SameTest->withinRoundoff
  ];
]

Module[{orbitG = KerrGeoOrbit[0.1, 10.0, 0.1, Cos[\[Pi]/4.]]},
  VerificationTest[
      SetPrecision[TeukolskyPointParticleMode[0, 2, 2, 0, 0, orbitG]["Amplitudes"], 6]
    ,
    <|"\[ScriptCapitalI]" -> -0.066318324275833298296962 + 0.0236010362126117048938706 I, 
     "\[ScriptCapitalH]" -> 0.000554901107100978524115932 - 0.000026750779539412935880262 I|>
    ,
    TestID->"Generic orbit Kerr amplitudes (s=0)",
    SameTest->withinRoundoff
  ];

  VerificationTest[
      SetPrecision[TeukolskyPointParticleMode[-2, 2, 2, 0, 0, orbitG]["Amplitudes"], 6]
    ,
    <|"\[ScriptCapitalI]" -> -0.0007470184367075 + 0.0002028616424164 I, 
     "\[ScriptCapitalH]" -> -0.0009960463270432 - 0.0001525666277914 I|>
    ,
    TestID->"Generic orbit Kerr amplitudes (s=-2)",
    SameTest->withinRoundoff
  ];
]
