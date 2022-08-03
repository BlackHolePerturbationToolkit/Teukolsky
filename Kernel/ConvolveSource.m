(* ::Package:: *)

(* ::Title:: *)
(*ConvolveSource*)


(* ::Section::Closed:: *)
(*Create Package*)


(* ::Subsection::Closed:: *)
(*BeginPackage*)


BeginPackage["Teukolsky`ConvolveSource`", {"KerrGeodesics`", "KerrGeodesics`OrbitalFrequencies`"}];


(* ::Subsection::Closed:: *)
(*Begin Private section*)


Begin["`Private`"];


(* ::Section::Closed:: *)
(*ConvolveSource*)


(* ::Subsection::Closed:: *)
(*ConvolveSource*)


ConvolveSource[l_Integer, m_Integer, n_Integer, k_Integer, R_, S_, TS_] :=
 Module[{s, orbit, e, x},
  If[TS["SourceType"] == "PointParticle",
    orbit = TS["Orbit"];
    s = TS["s"];
    {e, x} = {orbit["e"], orbit["Inclination"]};

    Which[
    {e, Abs[x]} == {0, 1},
      Return[ConvolveSourcePointParticleCircular[s,R,S,TS]],
    e == 0,
      Return[ConvolveSourcePointParticleSpherical[s,k,R,S,TS]],
    Abs[x] == 1,
      Return[ConvolveSourcePointParticleEccentric[s,n,R,S,TS]];,
    True,
      Return[ConvolveSourcePointParticleGeneric[s,n,k,R,S,TS]];
    ]
  ,
    Return[$Failed];
  ];
]


(* ::Subsection::Closed:: *)
(*s=-2 point particle on a circular orbit*)


ConvolveSourcePointParticleCircular[s:-2, R_, SH_, TS_] :=
 Module[{a, r0, m, \[Omega], \[Lambda], \[CapitalDelta], W, Ann0, Anmb0, Ambmb0, Anmb1, Ambmb1, Ambmb2, RIn, RUp, dRIn, dRUp, d2RIn, d2RUp, CIn, COut, ZIn, ZUp,S, dS, d2S, L2dagS, L1dagL2dagS, \[Rho], \[Rho]b, K},
  a  = TS["Orbit"]["a"];
  r0 = TS["Orbit"]["p"];
  m  = R["In"]["m"];
  \[Omega]  = R["In"]["\[Omega]"];
  \[Lambda] = R["Up"]["Eigenvalue"];

  RIn = R["In"][r0];
  RUp = R["Up"][r0];
  dRIn = R["In"]'[r0];
  dRUp = R["Up"]'[r0];
  d2RUp = (-(-\[Lambda] + 2 I r0 s 2 \[Omega] + (-2 I (-1 + r0) s (-a m + (a^2 + r0^2) \[Omega]) + (-a m + (a^2 + r0^2) \[Omega])^2)/(a^2 - 2 r0 + r0^2)) RUp - (-2 + 2 r0) (1 + s) dRUp)/(a^2 - 2 r0 + r0^2);
  d2RIn = (-(-\[Lambda] + 2 I r0 s 2 \[Omega] + (-2 I (-1 + r0) s (-a m + (a^2 + r0^2) \[Omega]) + (-a m + (a^2 + r0^2) \[Omega])^2)/(a^2 - 2 r0 + r0^2)) RIn - (-2 + 2 r0) (1 + s) dRIn)/(a^2 - 2 r0 + r0^2);

  \[CapitalDelta] = r0^2 - 2r0 + a^2;
  K = (r0^2 + a^2)\[Omega] - m a;

  W = 2 I \[Omega] R["In"]["Amplitudes"]["Incidence"];

  S = SH[\[Pi]/2, 0];
  dS = Derivative[1,0][SH][\[Pi]/2,0];
  d2S = Derivative[2,0][SH][\[Pi]/2,0];

  L2dagS =  dS + (a \[Omega] - m) S;
  L1dagL2dagS = (-2+(-m+a \[Omega])^2) S + 2 (-m+a \[Omega]) dS + d2S;

  \[Rho] = -1/r0;
  \[Rho]b = -1/r0;

  Ann0 = -((2 \[Rho]^-3 \[Rho]b^-1 TS["Cnn"])/\[CapitalDelta]^2)(L1dagL2dagS + 2 I a \[Rho] L2dagS);
  Anmb0 = -((2Sqrt[2] \[Rho]^-3 TS["Cnmb"])/\[CapitalDelta])(((I K)/\[CapitalDelta]-\[Rho] -\[Rho]b)L2dagS+((I K)/\[CapitalDelta]+\[Rho]+\[Rho]b)I a S(\[Rho]-\[Rho]b));
  Ambmb0 = S \[Rho]^-3 \[Rho]b TS["Cmbmb"]((K/\[CapitalDelta])^2 + 2I \[Rho] K/\[CapitalDelta] + 2I (a m (r0-1)+a^2 \[Omega]-r0^2 \[Omega])/\[CapitalDelta]^2);
  Anmb1 = -2Sqrt[2] \[Rho]^-3 TS["Cnmb"]/\[CapitalDelta](L2dagS + I a \[Rho] (\[Rho]-\[Rho]b) S);
  Ambmb1 = 2 S \[Rho]^-3 \[Rho]b TS["Cmbmb"](\[Rho]-(I K)/\[CapitalDelta]);
  Ambmb2 = -S \[Rho]^-3 \[Rho]b TS["Cmbmb"];

  ZIn = 2 \[Pi] (RUp(Ann0 + Anmb0 + Ambmb0) - dRUp(Anmb1 + Ambmb1) + d2RUp Ambmb2)/W;
  ZUp = 2 \[Pi] (RIn(Ann0 + Anmb0 + Ambmb0) - dRIn(Anmb1 + Ambmb1) + d2RIn Ambmb2)/W;

  Clear[a, r0, m, \[Omega], \[Lambda], \[CapitalDelta], W, Ann0, Anmb0, Ambmb0, Anmb1, Ambmb1, Ambmb2, RIn, RUp, dRIn, dRUp, d2RIn, d2RUp, CIn, COut, S, dS, d2S, L2dagS, L1dagL2dagS, \[Rho], \[Rho]b, K];

  <| "\[ScriptCapitalI]" -> ZUp, "\[ScriptCapitalH]" -> ZIn |>
]


(* ::Subsection::Closed:: *)
(*s=-2 point particle on a spherical orbit*)


ConvolveSourcePointParticleSpherical[s:-2, k_Integer, R_, SH_, TS_] :=
 Module[{a, r0, m, \[Omega], \[Lambda], \[CapitalUpsilon]t, \[CapitalUpsilon]\[Theta], T\[Theta], \[Theta]\[Lambda], t\[Lambda], \[Phi]\[Lambda], T\[Lambda], \[CapitalDelta], Kf, dK\[CapitalDelta], W, RUp, RIn, dRUp, dRIn, d2RUp, d2RIn, integrand, Zin, Zup},
  a  = TS["Orbit"]["a"];
  r0 = TS["Orbit"]["p"];
  m = R["In"]["m"];
  \[Omega]  = R["In"]["\[Omega]"];
  \[Lambda] = R["Up"]["Eigenvalue"];
  {\[CapitalUpsilon]\[Theta], \[CapitalUpsilon]t} = {"\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Theta]\)]\)","\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(t\)]\)"} /. TS["Orbit"]["Frequencies"];
  T\[Theta] = (2\[Pi])/(\[CapitalUpsilon]\[Theta]/\[CapitalUpsilon]t);

  \[Theta]\[Lambda] = TS["Orbit"]["Trajectory"][[3]];
  t\[Lambda] = TS["Orbit"]["Trajectory"][[1]];
  \[Phi]\[Lambda] = TS["Orbit"]["Trajectory"][[4]];
  T\[Lambda] = Derivative[1][t\[Lambda]];

  \[CapitalDelta] = r0^2-2r0+a^2;
  Kf=(r0^2+a^2)\[Omega]-m a;
  dK\[CapitalDelta]=(2 (a m (-1+r0)+a^2 \[Omega]-r0^2 \[Omega]))/\[CapitalDelta]^2;
  
  W = 2 I \[Omega] R["In"]["Amplitudes"]["Incidence"];

  RUp = R["Up"][r0];
  RIn = R["In"][r0];
  dRUp = R["Up"]'[r0];
  dRIn = R["In"]'[r0];
  d2RUp = (-(-\[Lambda] + 2 I r0 s 2 \[Omega] + (-2 I (-1 + r0) s (-a m + (a^2 + r0^2) \[Omega]) + (-a m + (a^2 + r0^2) \[Omega])^2)/(a^2 - 2 r0 + r0^2)) RUp - (-2 + 2 r0) (1 + s) dRUp)/(a^2 - 2 r0 + r0^2);
  d2RIn = (-(-\[Lambda] + 2 I r0 s 2 \[Omega] + (-2 I (-1 + r0) s (-a m + (a^2 + r0^2) \[Omega]) + (-a m + (a^2 + r0^2) \[Omega])^2)/(a^2 - 2 r0 + r0^2)) RIn - (-2 + 2 r0) (1 + s) dRIn)/(a^2 - 2 r0 + r0^2);

  integrand[\[Lambda]0_?NumericQ, comp_]:=
   Module[{SH\[Theta], dSH\[Theta], d2SH\[Theta], \[Theta], z, \[CapitalSigma], \[Rho], \[Rho]b, Ld2S, Ld1Ld2S, \[CapitalTheta], T, t0, \[Phi]0, CnnPlus, CnnMinus, CnmPlus, CnmMinus, CmmPlus, CmmMinus, Ann0Plus, Anm0Plus, Amm0Plus, Anm1Plus, Amm1Plus, Amm2Plus, Ann0Minus, Anm0Minus, Amm0Minus, Anm1Minus, Amm1Minus, Amm2Minus, IUpPlus, IUpMinus, IInPlus, IInMinus, res},
    \[Theta] = \[Theta]\[Lambda][\[Lambda]0];
    t0 = t\[Lambda][\[Lambda]0];
    \[Phi]0 = \[Phi]\[Lambda][\[Lambda]0];
    T = T\[Lambda][\[Lambda]0];

    \[CapitalSigma] = r0^2+a^2 Cos[\[Theta]]^2;
    \[Rho] = -1/(r0-I a Cos[\[Theta]]);
    \[Rho]b = -1/(r0+I a Cos[\[Theta]]);
    
    SH\[Theta] = SH[\[Theta],0];
    dSH\[Theta] = Derivative[1,0][SH][\[Theta],0];
    d2SH\[Theta] = Derivative[2,0][SH][\[Theta],0];
    Ld2S = dSH\[Theta]+(- m Csc[\[Theta]]+a \[Omega] Sin[\[Theta]]+2 Cot[\[Theta]])SH\[Theta];
    Ld1Ld2S = d2SH\[Theta] +(- m Csc[\[Theta]]+a \[Omega] Sin[\[Theta]]+2 Cot[\[Theta]])dSH\[Theta]+(a \[Omega] Cos[\[Theta]]+ m Cot[\[Theta]]Csc[\[Theta]]-2Csc[\[Theta]]^2)SH\[Theta]+(- m Csc[\[Theta]]+a \[Omega] Sin[\[Theta]]+ Cot[\[Theta]])Ld2S;

    CnnPlus = TS["Cnn+"][\[Lambda]0];
    CnmPlus = TS["Cnm+"][\[Lambda]0];
    CmmPlus = TS["Cmm+"][\[Lambda]0];

    Ann0Plus = (-2 \[Rho]^-3 \[Rho]b^-1 CnnPlus)/\[CapitalDelta]^2 (Ld1Ld2S+2I a \[Rho] Sin[\[Theta]]Ld2S);
    Anm0Plus = -((2Sqrt[2] \[Rho]^-3 CnmPlus)/\[CapitalDelta])(((I Kf)/\[CapitalDelta]-\[Rho]-\[Rho]b)Ld2S- Kf/\[CapitalDelta] a Sin[\[Theta]]SH\[Theta](\[Rho]-\[Rho]b));
    Amm0Plus = SH\[Theta] \[Rho]^-3 \[Rho]b CmmPlus((Kf/\[CapitalDelta])^2+2I \[Rho] Kf/\[CapitalDelta]+I dK\[CapitalDelta]);
    Anm1Plus = -((2Sqrt[2] \[Rho]^-3 CnmPlus)/\[CapitalDelta])(Ld2S+I a Sin[\[Theta]](\[Rho]-\[Rho]b)SH\[Theta]);
    Amm1Plus = 2SH\[Theta] \[Rho]^-3 \[Rho]b CmmPlus(\[Rho]-(I Kf)/\[CapitalDelta]);
    Amm2Plus = -SH\[Theta] \[Rho]^-3 \[Rho]b CmmPlus;
    Ann0Minus = (-2 \[Rho]^-3 \[Rho]b^-1 CnnMinus)/\[CapitalDelta]^2 (Ld1Ld2S+2I a \[Rho] Sin[\[Theta]]Ld2S);
    Anm0Minus = -((2Sqrt[2] \[Rho]^-3 CnmMinus)/\[CapitalDelta])(((I Kf)/\[CapitalDelta]-\[Rho]-\[Rho]b)Ld2S- Kf/\[CapitalDelta] a Sin[\[Theta]]SH\[Theta](\[Rho]-\[Rho]b));
    Amm0Minus = SH\[Theta] \[Rho]^-3 \[Rho]b CmmMinus((Kf/\[CapitalDelta])^2+2I \[Rho] Kf/\[CapitalDelta]+I dK\[CapitalDelta]);
    Anm1Minus = -((2Sqrt[2] \[Rho]^-3 CnmMinus)/\[CapitalDelta])(Ld2S+I a Sin[\[Theta]](\[Rho]-\[Rho]b)SH\[Theta]);
    Amm1Minus = 2SH\[Theta] \[Rho]^-3 \[Rho]b CmmMinus(\[Rho]-(I Kf)/\[CapitalDelta]);
    Amm2Minus = -SH\[Theta]  \[Rho]^-3 \[Rho]b CmmMinus;

    IUpPlus=RUp(Ann0Plus + Anm0Plus + Amm0Plus)-dRUp(Anm1Plus+Amm1Plus)+d2RUp Amm2Plus;
    IInPlus=RIn(Ann0Plus + Anm0Plus + Amm0Plus)-dRIn(Anm1Plus+Amm1Plus)+d2RIn Amm2Plus;
    
    res = {T(E^(I(\[Omega] t0-m \[Phi]0)) IInPlus), T(E^(I(\[Omega] t0-m \[Phi]0)) IUpPlus)}[[comp]];

    Clear[SH\[Theta], dSH\[Theta], d2SH\[Theta], \[Theta], z, \[CapitalSigma], \[Rho], \[Rho]b, Ld2S, Ld1Ld2S, \[CapitalTheta], T, t0, \[Phi]0, CnnPlus, CnnMinus, CnmPlus, CnmMinus, CmmPlus, CmmMinus, Ann0Plus, Anm0Plus, Amm0Plus, Anm1Plus, Amm1Plus, Amm2Plus, Ann0Minus, Anm0Minus, Amm0Minus, Anm1Minus, Amm1Minus, Amm2Minus, IUpPlus, IUpMinus, IInPlus, IInMinus];
    res    
  ];

  Zin = ((2\[Pi])/(T\[Theta] W)) Quiet[NIntegrate[integrand[\[Lambda]0,1], {\[Lambda]0, 0, (2\[Pi])/\[CapitalUpsilon]\[Theta]}, WorkingPrecision -> Precision[\[Omega]], Method -> {"Trapezoidal", "SymbolicProcessing"->0}], NIntegrate::precw];
  Zup = ((2\[Pi])/(T\[Theta] W)) Quiet[NIntegrate[integrand[\[Lambda]0,2], {\[Lambda]0, 0, (2\[Pi])/\[CapitalUpsilon]\[Theta]}, WorkingPrecision -> Precision[\[Omega]], Method -> {"Trapezoidal", "SymbolicProcessing"->0}], NIntegrate::precw];
  
  Clear[a, r0, m, \[Omega], \[Lambda], \[CapitalUpsilon]t, \[CapitalUpsilon]\[Theta], T\[Theta], \[Theta]\[Lambda], t\[Lambda], \[Phi]\[Lambda], T\[Lambda], \[CapitalDelta], Kf, dK\[CapitalDelta], W, RUp, RIn, dRUp, dRIn, d2RUp, d2RIn];
  Remove[integrand];

  <|"\[ScriptCapitalI]" -> Zin, "\[ScriptCapitalH]" -> Zup|>
]


(* ::Subsection::Closed:: *)
(*s=-2 point particle on an eccentric orbit*)


ConvolveSourcePointParticleEccentric[s_:-2, n_Integer, R_, SH_, TS_] :=
 Module[{a, p, rpi, \[CapitalUpsilon]t, \[CapitalUpsilon]r, \[CapitalUpsilon]\[Theta], \[CapitalDelta]tr, \[CapitalDelta]\[Phi]r, qt0, qr0, q\[Theta]0, q\[Phi]0, m, \[Omega], \[Lambda], W, rq, \[Theta]0, S0, dS0, d2S0, L1, L2, L2S, L2p, L1Sp, L1L2S, integrand, \[Alpha]In, \[Alpha]Up, \[Xi], wpIn, wpUp, ZIn, ZUp},
  a = TS["Orbit"]["a"];
  p = TS["Orbit"]["p"];
  rpi = TS["Orbit"]["Trajectory"][[2]];

  {\[CapitalUpsilon]t, \[CapitalUpsilon]r, \[CapitalUpsilon]\[Theta]} = {"\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(t\)]\)", "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(r\)]\)", "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Theta]\)]\)"} /. TS["Orbit"]["Frequencies"];
  {\[CapitalDelta]tr, \[CapitalDelta]\[Phi]r} = {"\[CapitalDelta]tr", "\[CapitalDelta]\[Phi]r"} /. TS["Orbit"]["TrajectoryDeltas"];
  {qt0, qr0, q\[Theta]0, q\[Phi]0} = TS["Orbit"]["InitialPhases"];

  m = R["In"]["m"];
  \[Omega] = R["In"]["\[Omega]"];
  \[Lambda] = R["In"]["Eigenvalue"];

  W = 2 I \[Omega] R["In"]["Amplitudes"]["Incidence"];
  
  rq[qr_] := rpi[(qr-qr0)/\[CapitalUpsilon]r];
  \[Theta]0 = Pi/2;
  
  S0 = SH[\[Theta]0, 0];
  dS0 = Derivative[1,0][SH][\[Theta]0, 0];
  d2S0 = Derivative[2,0][SH][\[Theta]0, 0];
  L1 = -m/Sin[\[Theta]0] + a \[Omega] Sin[\[Theta]0] + Cos[\[Theta]0]/Sin[\[Theta]0];
  L2 = -m/Sin[\[Theta]0] + a \[Omega] Sin[\[Theta]0] + 2 Cos[\[Theta]0]/Sin[\[Theta]0];
  L2S = dS0 + L2 S0;
  L2p = m Cos[\[Theta]0]/Sin[\[Theta]0]^2 + a \[Omega] Cos[\[Theta]0] - 2/Sin[\[Theta]0]^2;
  L1Sp = d2S0 + L1 dS0;
  L1L2S = L1Sp + L2p S0 + L2 dS0 + L1 L2 S0;

  integrand[qr_, RF_]:=
   Module[{r0, R0, dR0, d2R0, \[CapitalDelta], Kt, \[Rho], \[Rho]bar, Ann0, Anmbar0, Anmbar1, Ambarmbar0, Ambarmbar1, Ambarmbar2, Cnnp1p1, Cnmbarp1p1, Cmbarmbarp1p1, Cnnm1p1, Cnmbarm1p1, Cmbarmbarm1p1, rphase, res},
    r0 = rq[qr];
    R0   = RF[r0];
    dR0  = RF'[r0];
    d2R0 = (-(-\[Lambda] + 2 I r0 s 2 \[Omega] + (-2 I (-1 + r0) s (-a m + (a^2 + r0^2) \[Omega]) + (-a m + (a^2 + r0^2) \[Omega])^2)/(a^2 - 2 r0 + r0^2)) R0 - (-2 + 2 r0) (1 + s) dR0)/(a^2 - 2 r0 + r0^2);

    \[CapitalDelta] = r0^2 + a^2 - 2 r0;
    Kt = (r0^2 + a^2) \[Omega] - m a;
    \[Rho] = -1/(r0 - I a Cos[\[Theta]0]);
    \[Rho]bar = -1/(r0 + I a Cos[\[Theta]0]);
  
    Ann0 = -\[Rho]^(-2) \[Rho]bar^(-1) (Sqrt[2] \[CapitalDelta])^(-2) (\[Rho]^(-1) L1L2S + 3 I a Sin[\[Theta]0] L1 S0 + 3 I a Cos[\[Theta]0] S0 + 2 I a Sin[\[Theta]0] dS0 - I a Sin[\[Theta]0] L2 S0 );
    Anmbar0 = \[Rho]^(-3) (Sqrt[2]\[CapitalDelta])^(-1) ( (\[Rho] + \[Rho]bar - I Kt/\[CapitalDelta]) L2S + (\[Rho] - \[Rho]bar) a Sin[\[Theta]0] Kt/\[CapitalDelta] S0 );
    Anmbar1 = -\[Rho]^(-3) (Sqrt[2]\[CapitalDelta])^(-1) ( L2S + I (\[Rho] - \[Rho]bar) a Sin[\[Theta]0] S0 );
    Ambarmbar0 = \[Rho]^(-3) \[Rho]bar S0 Kt/\[CapitalDelta]/4 ( I (2 \[Omega] r0/Kt - 2(r0 - 1)/\[CapitalDelta]) + Kt/\[CapitalDelta] + 2 I \[Rho]);
    Ambarmbar1 = -\[Rho]^(-3) \[Rho]bar S0/2 ( I Kt/\[CapitalDelta] - \[Rho] );
    Ambarmbar2 = -\[Rho]^(-3) \[Rho]bar S0/4;

    {Cnnp1p1,Cnmbarp1p1,Cmbarmbarp1p1} = TS["Cab"][r0, \[Theta]0, 1, 1];
    {Cnnm1p1,Cnmbarm1p1,Cmbarmbarm1p1} = TS["Cab"][r0, \[Theta]0, -1, 1];
    
    rphase = \[Omega] \[CapitalDelta]tr[qr] - m \[CapitalDelta]\[Phi]r[qr] + n qr;
    
    res = ((Ann0 Cnnp1p1 + Anmbar0 Cnmbarp1p1 + Ambarmbar0 Cmbarmbarp1p1) R0-(Anmbar1 Cnmbarp1p1 + Ambarmbar1 Cmbarmbarp1p1) dR0+Ambarmbar2 Cmbarmbarp1p1 d2R0)Exp[I rphase] 
        + ((Ann0 Cnnm1p1 + Anmbar0 Cnmbarm1p1 + Ambarmbar0 Cmbarmbarm1p1) R0-(Anmbar1 Cnmbarm1p1 + Ambarmbar1 Cmbarmbarm1p1) dR0+Ambarmbar2 Cmbarmbarm1p1 d2R0)Exp[-I rphase];

    Clear[r0, R0, dR0, d2R0, \[CapitalDelta], Kt, \[Rho], \[Rho]bar, Ann0, Anmbar0, Anmbar1, Ambarmbar0, Ambarmbar1, Ambarmbar2, Cnnp1p1, Cnmbarp1p1, Cmbarmbarp1p1, Cnnm1p1, Cnmbarm1p1, Cmbarmbarm1p1, rphase];
    res
  ];

  wpIn = Precision[integrand[0,R["In"]]];
  wpUp = Precision[integrand[0,R["Up"]]];

  \[Alpha]In = 1/(2\[Pi]) Quiet[NIntegrate[integrand[qr,R["In"]], {qr, 0, \[Pi]}, Method -> {"Trapezoidal", "SymbolicProcessing"->0}, WorkingPrecision -> wpIn], NIntegrate::precw];
  \[Alpha]Up = 1/(2\[Pi]) Quiet[NIntegrate[integrand[qr,R["Up"]], {qr, 0, \[Pi]}, Method -> {"Trapezoidal", "SymbolicProcessing"->0}, WorkingPrecision -> wpUp], NIntegrate::precw];
  
  \[Xi] = m(\[CapitalDelta]\[Phi]r[qr0]-q\[Phi]0) - \[Omega](\[CapitalDelta]tr[qr0]-qt0) - n qr0;
  
  ZIn = 8Pi \[Alpha]Up/W/\[CapitalUpsilon]t Exp[I \[Xi]];
  ZUp = 8Pi \[Alpha]In/W/\[CapitalUpsilon]t Exp[I \[Xi]];

  Clear[a, p, rpi, \[CapitalUpsilon]t, \[CapitalUpsilon]r, \[CapitalUpsilon]\[Theta], \[CapitalDelta]tr, \[CapitalDelta]\[Phi]r, qt0, qr0, q\[Theta]0, q\[Phi]0, m, \[Omega], \[Lambda], W, rq, \[Theta]0, S0, dS0, d2S0, L1, L2, L2S, L2p, L1Sp, L1L2S, integrand, \[Alpha]In, \[Alpha]Up, \[Xi], wpIn, wpUp];
  <| "\[ScriptCapitalI]" -> ZUp, "\[ScriptCapitalH]" -> ZIn |>
]


(* ::Subsection::Closed:: *)
(*s=-2 point particle on a generic orbit*)


ConvolveSourcePointParticleGeneric[s_:-2, n_Integer, k_Integer, R_, SH_, TS_] :=
 Module[{a, p, rpi, \[Theta]pi, \[CapitalUpsilon]t, \[CapitalUpsilon]r, \[CapitalUpsilon]\[Theta], \[CapitalDelta]tr, \[CapitalDelta]t\[Theta], \[CapitalDelta]\[Phi]r, \[CapitalDelta]\[Phi]\[Theta], qt0, qr0, q\[Theta]0, q\[Phi]0, m, \[Omega], \[Lambda], W, rq, \[Theta]q, integrand, integrandIn, integrandUp, \[Alpha]In, \[Alpha]Up, \[Xi], wpIn, wpUp, ZIn, ZUp},
  a = TS["Orbit"]["a"];
  p = TS["Orbit"]["p"];
  {rpi, \[Theta]pi}=TS["Orbit"]["Trajectory"][[2;;3]];

  {\[CapitalUpsilon]t, \[CapitalUpsilon]r, \[CapitalUpsilon]\[Theta]} = {"\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(t\)]\)", "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(r\)]\)", "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Theta]\)]\)"} /. TS["Orbit"]["Frequencies"];
  {\[CapitalDelta]tr, \[CapitalDelta]t\[Theta], \[CapitalDelta]\[Phi]r, \[CapitalDelta]\[Phi]\[Theta]} = {"\[CapitalDelta]tr", "\[CapitalDelta]t\[Theta]", "\[CapitalDelta]\[Phi]r", "\[CapitalDelta]\[Phi]\[Theta]"} /. TS["Orbit"]["TrajectoryDeltas"];
  {qt0, qr0, q\[Theta]0, q\[Phi]0} = TS["Orbit"]["InitialPhases"];

  m = R["In"]["m"];
  \[Omega] = R["In"]["\[Omega]"];
  \[Lambda] = R["In"]["Eigenvalue"];

  W = 2 I \[Omega] R["In"]["Amplitudes"]["Incidence"];

  rq[qr_] := rpi[(qr-qr0)/\[CapitalUpsilon]r];
  \[Theta]q[q\[Theta]_] := \[Theta]pi[(q\[Theta]-q\[Theta]0)/\[CapitalUpsilon]\[Theta]];

  integrand[qr_, q\[Theta]_, RF_]:=
   Module[{r0, R0, dR0, d2R0, \[Theta]0, S0, dS0, d2S0, L1, L2, L2S, L2p, L1Sp, L1L2S, \[CapitalDelta], Kt, \[Rho], \[Rho]bar, Ann0, Anmbar0, Anmbar1, Ambarmbar0, Ambarmbar1, Ambarmbar2, Cnnp1p1, Cnmbarp1p1, Cmbarmbarp1p1, Cnnp1m1, Cnmbarp1m1, Cmbarmbarp1m1, Cnnm1p1, Cnmbarm1p1, Cmbarmbarm1p1, Cnnm1m1,Cnmbarm1m1,Cmbarmbarm1m1, rphase, \[Theta]phase, res},
    r0 = rq[qr];
    R0   = RF[r0];
    dR0  = RF'[r0];
    d2R0 = (-(-\[Lambda] + 2 I r0 s 2 \[Omega] + (-2 I (-1 + r0) s (-a m + (a^2 + r0^2) \[Omega]) + (-a m + (a^2 + r0^2) \[Omega])^2)/(a^2 - 2 r0 + r0^2)) R0 - (-2 + 2 r0) (1 + s) dR0)/(a^2 - 2 r0 + r0^2);
  
    \[Theta]0 = \[Theta]q[q\[Theta]];
    S0 = SH[\[Theta]0, 0];
    dS0 = Derivative[1,0][SH][\[Theta]0, 0];
    d2S0 = Derivative[2,0][SH][\[Theta]0, 0];
    L1 = -m/Sin[\[Theta]0] + a \[Omega] Sin[\[Theta]0] + Cos[\[Theta]0]/Sin[\[Theta]0];
    L2 = -m/Sin[\[Theta]0] + a \[Omega] Sin[\[Theta]0] + 2 Cos[\[Theta]0]/Sin[\[Theta]0];
    L2S = dS0 + L2 S0;
    L2p = m Cos[\[Theta]0]/Sin[\[Theta]0]^2 + a \[Omega] Cos[\[Theta]0] - 2/Sin[\[Theta]0]^2;
    L1Sp = d2S0 + L1 dS0;
    L1L2S = L1Sp + L2p S0 + L2 dS0 + L1 L2 S0;
  
    \[CapitalDelta] = r0^2 + a^2 - 2 r0;
    Kt = (r0^2 + a^2) \[Omega] - m a;
    \[Rho] = -1/(r0 - I a Cos[\[Theta]0]);
    \[Rho]bar = -1/(r0 + I a Cos[\[Theta]0]);
  
    Ann0 = -\[Rho]^(-2) \[Rho]bar^(-1) (Sqrt[2] \[CapitalDelta])^(-2) (\[Rho]^(-1) L1L2S + 3 I a Sin[\[Theta]0] L1 S0 + 3 I a Cos[\[Theta]0] S0 + 2 I a Sin[\[Theta]0] dS0 - I a Sin[\[Theta]0] L2 S0 );
    Anmbar0 = \[Rho]^(-3) (Sqrt[2]\[CapitalDelta])^(-1) ( (\[Rho] + \[Rho]bar - I Kt/\[CapitalDelta]) L2S + (\[Rho] - \[Rho]bar) a Sin[\[Theta]0] Kt/\[CapitalDelta] S0 );
    Anmbar1 = -\[Rho]^(-3) (Sqrt[2]\[CapitalDelta])^(-1) ( L2S + I (\[Rho] - \[Rho]bar) a Sin[\[Theta]0] S0 );
    Ambarmbar0 = \[Rho]^(-3) \[Rho]bar S0 Kt/\[CapitalDelta]/4 ( I (2 \[Omega] r0/Kt - 2(r0 - 1)/\[CapitalDelta]) + Kt/\[CapitalDelta] + 2 I \[Rho]);
    Ambarmbar1 = -\[Rho]^(-3) \[Rho]bar S0/2 ( I Kt/\[CapitalDelta] - \[Rho] );
    Ambarmbar2 = -\[Rho]^(-3) \[Rho]bar S0/4;
    
    {Cnnp1p1,Cnmbarp1p1,Cmbarmbarp1p1} = TS["Cab"][r0, \[Theta]0, 1, 1];
    {Cnnp1m1,Cnmbarp1m1,Cmbarmbarp1m1} = TS["Cab"][r0, \[Theta]0, 1, -1];
    {Cnnm1p1,Cnmbarm1p1,Cmbarmbarm1p1} = TS["Cab"][r0, \[Theta]0, -1, 1];
    {Cnnm1m1,Cnmbarm1m1,Cmbarmbarm1m1} = TS["Cab"][r0, \[Theta]0, -1, -1];
    
    rphase = \[Omega] \[CapitalDelta]tr[qr] - m \[CapitalDelta]\[Phi]r[qr] + n qr;
    \[Theta]phase = \[Omega] \[CapitalDelta]t\[Theta][q\[Theta]] - m \[CapitalDelta]\[Phi]\[Theta][q\[Theta]] + k q\[Theta];
    
    res = ((Ann0*Cnnp1p1 + Anmbar0*Cnmbarp1p1 + Ambarmbar0*Cmbarmbarp1p1) R0-(Anmbar1*Cnmbarp1p1 + Ambarmbar1*Cmbarmbarp1p1) dR0+Ambarmbar2*Cmbarmbarp1p1 d2R0)Exp[I(rphase+\[Theta]phase)] 
        + ((Ann0*Cnnp1m1 + Anmbar0*Cnmbarp1m1 + Ambarmbar0*Cmbarmbarp1m1) R0-(Anmbar1*Cnmbarp1m1 + Ambarmbar1*Cmbarmbarp1m1) dR0+Ambarmbar2*Cmbarmbarp1m1 d2R0)Exp[I(rphase-\[Theta]phase)] 
        + ((Ann0*Cnnm1p1 + Anmbar0*Cnmbarm1p1 + Ambarmbar0*Cmbarmbarm1p1) R0-(Anmbar1*Cnmbarm1p1 + Ambarmbar1*Cmbarmbarm1p1) dR0+Ambarmbar2*Cmbarmbarm1p1 d2R0)Exp[-I(rphase-\[Theta]phase)] 
        + ((Ann0*Cnnm1m1 + Anmbar0*Cnmbarm1m1 + Ambarmbar0*Cmbarmbarm1m1) R0-(Anmbar1*Cnmbarm1m1 + Ambarmbar1*Cmbarmbarm1m1) dR0+Ambarmbar2*Cmbarmbarm1m1 d2R0)Exp[-I(rphase+\[Theta]phase)];
    Clear[r0, R0, dR0, d2R0, \[Theta]0, S0, dS0, d2S0, L1, L2, L2S, L2p, L1Sp, L1L2S, \[CapitalDelta], Kt, \[Rho], \[Rho]bar, Ann0, Anmbar0, Anmbar1, Ambarmbar0, Ambarmbar1, Ambarmbar2, Cnnp1p1, Cnmbarp1p1, Cmbarmbarp1p1, Cnnp1m1, Cnmbarp1m1, Cmbarmbarp1m1, Cnnm1p1, Cnmbarm1p1, Cmbarmbarm1p1, Cnnm1m1,Cnmbarm1m1,Cmbarmbarm1m1, rphase, \[Theta]phase];
    res
  ];
  
  integrandIn[qr_, q\[Theta]_] := integrand[qr,q\[Theta],R["In"]];
  integrandUp[qr_, q\[Theta]_] := integrand[qr,q\[Theta],R["Up"]];

  wpIn = Precision[integrandIn[0, 0]];
  wpUp = Precision[integrandUp[0, 0]];

  \[Alpha]In = 1/(2\[Pi])^2 Quiet[NIntegrate[integrandIn[qr, q\[Theta]], {qr, 0, \[Pi]}, {q\[Theta], 0, \[Pi]}, Method -> {"Trapezoidal", "SymbolicProcessing"->0}, WorkingPrecision -> wpIn], NIntegrate::precw];
  \[Alpha]Up = 1/(2\[Pi])^2 Quiet[NIntegrate[integrandUp[qr, q\[Theta]], {qr, 0, \[Pi]}, {q\[Theta], 0, \[Pi]}, Method -> {"Trapezoidal", "SymbolicProcessing"->0}, WorkingPrecision -> wpUp], NIntegrate::precw];
  
  \[Xi] = m(\[CapitalDelta]\[Phi]r[qr0]+\[CapitalDelta]\[Phi]\[Theta][q\[Theta]0]-q\[Phi]0) - \[Omega](\[CapitalDelta]tr[qr0]+\[CapitalDelta]t\[Theta][q\[Theta]0]-qt0) - k q\[Theta]0 - n qr0;

  ZIn = 8Pi \[Alpha]Up/W/\[CapitalUpsilon]t Exp[I \[Xi]];
  ZUp = 8Pi \[Alpha]In/W/\[CapitalUpsilon]t Exp[I \[Xi]];

  Clear[a, p, rpi, \[Theta]pi, \[CapitalUpsilon]t, \[CapitalUpsilon]r, \[CapitalUpsilon]\[Theta], \[CapitalDelta]tr, \[CapitalDelta]t\[Theta], \[CapitalDelta]\[Phi]r, \[CapitalDelta]\[Phi]\[Theta], qt0, qr0, q\[Theta]0, q\[Phi]0, m, \[Omega], \[Lambda], W, rq, \[Theta]q, integrand, integrandIn, integrandUp, \[Alpha]In, \[Alpha]Up, \[Xi], wpIn, wpUp];
  <| "\[ScriptCapitalI]" -> ZUp, "\[ScriptCapitalH]" -> ZIn |>
]


(* ::Subsection::Closed:: *)
(*s=\[PlusMinus]1 point particle on a circular orbit*)


ConvolveSourcePointParticleCircular[s:(-1|+1), R_, SH_, TS_] :=
 Module[{a, r0, \[Theta], m, \[Omega], \[CapitalDelta], \[CapitalDelta]p, W, A, B, RIn, ROut, dRIn, dROut, ZIn, ZOut, dS, S, PIn, POut, dPIn, dPOut},
  a  = TS["Orbit"]["a"];
  r0 = TS["Orbit"]["p"];
  m  = R["In"]["m"];
  \[Omega]  = R["In"]["\[Omega]"];

  \[CapitalDelta] = r0^2 - 2r0 + a^2;
  \[CapitalDelta]p = 2r0 - 2;

  S = SH[\[Pi]/2, 0];
  dS = Derivative[1,0][SH][\[Pi]/2,0];

  (* s = -1 radial functions *)
  RIn = R["In"][r0];
  ROut = R["Up"][r0];

  dRIn = R["In"]'[r0];
  dROut = R["Up"]'[r0];

  (* Convert R -> P: P+1 = \[CapitalDelta] R+1 and P-1 = R-1 *)
  Which[
    s == -1,
      PIn = RIn;
      POut = ROut;
      dPIn = dRIn;
      dPOut = dROut;,
    s == +1,
      PIn = \[CapitalDelta] RIn;
      POut = \[CapitalDelta] ROut;
      dPIn = \[CapitalDelta] dRIn + \[CapitalDelta]p RIn;
      dPOut = \[CapitalDelta] dROut + \[CapitalDelta]p ROut;
  ];

  (* Wronskian *)
  W = (PIn dPOut - POut dPIn);

  (* Define source terms: arXiv:2008.12703 Eq. (48) *)
  A = TS["\[ScriptCapitalS]"] ((m*TS["\!\(\*SuperscriptBox[\(A\), \((r)\)]\)"] + s I TS["\!\(\*SuperscriptBox[OverscriptBox[\(A\), \(~\)], \((i)\)]\)"])S + s TS["C"]dS);
  B = s I TS["\[ScriptCapitalS]"]TS["B"]*S;

  ZIn  = (POut*A - dPOut*B)/(\[CapitalDelta]*W);
  ZOut = (PIn*A - dPIn*B)/(\[CapitalDelta]*W);

  Clear[a, r0, \[Theta], m, \[Omega], \[CapitalDelta], \[CapitalDelta]p, W, A, B, RIn, ROut, dRIn, dROut, dS, S, PIn, POut, dPIn, dPOut];

  <| "\[ScriptCapitalI]" -> ZOut, "\[ScriptCapitalH]" -> ZIn |>
]


(* ::Subsection::Closed:: *)
(*s=0 point particle on a circular orbit*)


ConvolveSourcePointParticleCircular[0, R_, SH_, TS_] :=
 Module[{a, r0, \[Omega], RIn, RUp, dRIn, dRUp, \[Psi]In, \[Psi]Out, d\[Psi]In, d\[Psi]Out, W, \[Alpha], ZIn, ZOut, S},
  a  = TS["Orbit"]["a"];
  r0 = TS["Orbit"]["p"];
  \[Omega] = R["In"]["\[Omega]"];

  RIn = R["In"][r0];
  dRIn = R["In"]'[r0];
  RUp = R["Up"][r0];
  dRUp = R["Up"]'[r0];

  W = 2 I \[Omega] R["In"]["Amplitudes"]["Incidence"];

  S = SH[\[Pi]/2,0];
  \[Alpha] = (r0^2 - 2 r0 + a^2)/r0 TS["\[Alpha]"] S;

  ZIn = \[Alpha] RUp/W;
  ZOut = \[Alpha] RIn/W;

  Clear[a, r0, \[Omega], RIn, RUp, dRIn, dRUp, \[Psi]In, \[Psi]Out, d\[Psi]In, d\[Psi]Out, W, \[Alpha], S];

  <| "\[ScriptCapitalI]" -> ZOut, "\[ScriptCapitalH]" -> ZIn |>
]


(* ::Subsection::Closed:: *)
(*s=0 point particle on a spherical orbit*)


ConvolveSourcePointParticleSpherical[0, k_Integer, R_, SH_, TS_] :=
 Module[{a, r0, RIn, RUp, dRIn, dRUp, W, \[Alpha], S, tp, rp, \[Theta]p, \[Phi]p, l, m, \[Omega], \[CapitalUpsilon]\[Theta], ZIn, ZOut},
  a  = TS["Orbit"]["a"];
  r0 = TS["Orbit"]["p"];
  \[CapitalUpsilon]\[Theta] = TS["Orbit"]["Frequencies"]["\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Theta]\)]\)"];
  {tp,rp,\[Theta]p,\[Phi]p} = TS["Orbit"]["Trajectory"];
 
  \[Omega] = R["In"]["\[Omega]"];
  l = R["In"]["l"];
  m = R["In"]["m"];

  If[a!=0,
    If[MatchQ[l+m+k,_?OddQ], Return[<| "\[ScriptCapitalI]" -> 0, "\[ScriptCapitalH]" -> 0 |>]],
    If[MatchQ[l+m+k,_?OddQ]||(Abs[m+k]>l), Return[<| "\[ScriptCapitalI]" -> 0, "\[ScriptCapitalH]" -> 0 |>]];
  ];

  RIn = R["In"][r0];
  RUp = R["Up"][r0];

  dRIn = R["In"]'[r0];
  dRUp = R["Up"]'[r0];

  W = 2 I \[Omega] R["In"]["Amplitudes"]["Incidence"];

  \[Alpha] = (r0^2 - 2 r0 + a^2)/r0 Quiet[NIntegrate[TS["\[Alpha]"][\[Theta]p[\[Lambda]]]SH[\[Theta]p[\[Lambda]],0] Cos[\[Omega] tp[\[Lambda]] - m \[Phi]p[\[Lambda]]], {\[Lambda], 0, \[Pi]/(2\[CapitalUpsilon]\[Theta])},
      Method -> {"Trapezoidal", "SymbolicProcessing"->0}, WorkingPrecision -> Precision[\[Omega]]], NIntegrate::precw];

  ZIn = \[Alpha] RUp/W;
  ZOut = \[Alpha] RIn/W;

  Clear[a, r0, RIn, RUp, dRIn, dRUp, W, \[Alpha], S, tp, rp, \[Theta]p, \[Phi]p, l, m, \[Omega], \[CapitalUpsilon]\[Theta]];
  
  <| "\[ScriptCapitalI]" -> ZOut, "\[ScriptCapitalH]" -> ZIn |>
]


(* ::Subsection::Closed:: *)
(*s=0 point particle on an eccentric orbit*)


ConvolveSourcePointParticleEccentric[0, n_Integer, R_, S_, TS_] :=
 Module[{a, p, rpi, rq, \[Theta]q, qt0, qr0, q\[Theta]0, q\[Phi]0, \[Xi], W, \[Alpha]1In, \[Alpha]1Up, S0, l, m, \[Omega], \[CapitalUpsilon]t, \[CapitalUpsilon]r, \[CapitalDelta]tr, \[CapitalDelta]\[Phi]r, II1In, II1Up, ZIn, ZUp},
  a = TS["Orbit"]["a"];
  p = TS["Orbit"]["p"];
  rpi = TS["Orbit"]["Trajectory"][[2]];
  \[Theta]q = Pi/2;

  {\[CapitalUpsilon]t, \[CapitalUpsilon]r} = {"\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(t\)]\)","\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(r\)]\)"} /. TS["Orbit"]["Frequencies"];
  {\[CapitalDelta]tr, \[CapitalDelta]\[Phi]r} = {"\[CapitalDelta]tr", "\[CapitalDelta]\[Phi]r"} /. TS["Orbit"]["TrajectoryDeltas"];
  {qt0, qr0, q\[Theta]0, q\[Phi]0} = TS["Orbit"]["InitialPhases"];

  rq[qr_]:=rpi[(qr-qr0)/\[CapitalUpsilon]r];

  S0 = S[\[Theta]q, 0];
  
  \[Omega] = R["In"]["\[Omega]"];
  l = R["In"]["l"];
  m = R["In"]["m"];

  W = 2 I \[Omega] R["In"]["Amplitudes"]["Incidence"];

  II1In = Function[qr, Evaluate[2TS["\[Alpha]"][rq[qr],\[Theta]q]R["In"][rq[qr]]Cos[(\[Omega] \[CapitalDelta]tr[qr] - m \[CapitalDelta]\[Phi]r[qr] + n qr)]]];
  II1Up = Function[qr, Evaluate[2TS["\[Alpha]"][rq[qr],\[Theta]q]R["Up"][rq[qr]]Cos[(\[Omega] \[CapitalDelta]tr[qr] - m \[CapitalDelta]\[Phi]r[qr] + n qr)]]];

  \[Alpha]1In = 1/(2\[Pi]) S0 Quiet[NIntegrate[II1In[qr], {qr,0,\[Pi]}, Method -> {"Trapezoidal", "SymbolicProcessing"->0}, WorkingPrecision -> Precision[\[Omega]]], NIntegrate::precw];
  \[Alpha]1Up = 1/(2\[Pi]) S0 Quiet[NIntegrate[II1Up[qr], {qr,0,\[Pi]}, Method -> {"Trapezoidal", "SymbolicProcessing"->0}, WorkingPrecision -> Precision[\[Omega]]], NIntegrate::precw];

  \[Xi]=m(\[CapitalDelta]\[Phi]r[qr0]-q\[Phi]0) - \[Omega](\[CapitalDelta]tr[qr0]-qt0) - n qr0;

  ZIn = \[Alpha]1In/W/\[CapitalUpsilon]t Exp[I \[Xi]];
  ZUp = \[Alpha]1Up/W/\[CapitalUpsilon]t Exp[I \[Xi]];

  Clear[a, p, rpi, rq, \[Theta]q, qt0, qr0, q\[Theta]0, q\[Phi]0, \[Xi], W, \[Alpha]1In, \[Alpha]1Up, S0, l, m, \[Omega], \[CapitalUpsilon]t, \[CapitalUpsilon]r, \[CapitalDelta]tr, \[CapitalDelta]\[Phi]r, II1In, II1Up];
  <| "\[ScriptCapitalI]" -> ZUp, "\[ScriptCapitalH]" -> ZIn |>
]


(* ::Subsection::Closed:: *)
(*s=0 point particle on a generic orbit*)


ConvolveSourcePointParticleGeneric[0, n_Integer, k_Integer, R_, SH_, TS_] :=
Module[{a, p, rpi, \[Theta]pi, rq, \[Theta]q, qt0, qr0, q\[Theta]0, q\[Phi]0, \[Xi], W, \[Alpha]1In, \[Alpha]1Up, \[Alpha]2, \[Alpha]3In, \[Alpha]3Up, \[Alpha]4, S, l, m, \[Omega], \[CapitalUpsilon]t, \[CapitalUpsilon]r, \[CapitalUpsilon]\[Theta], \[CapitalDelta]tr, \[CapitalDelta]t\[Theta], \[CapitalDelta]\[Phi]r, \[CapitalDelta]\[Phi]\[Theta], II1In, II1Up, II2, II3In, II3Up, II4, ZIn, ZUp},
  a = TS["Orbit"]["a"];
  p = TS["Orbit"]["p"];
  {rpi, \[Theta]pi} = TS["Orbit"]["Trajectory"][[2;;3]];

  {\[CapitalUpsilon]t, \[CapitalUpsilon]r, \[CapitalUpsilon]\[Theta]} = {"\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(t\)]\)","\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(r\)]\)","\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Theta]\)]\)"} /. TS["Orbit"]["Frequencies"];
  {\[CapitalDelta]tr, \[CapitalDelta]t\[Theta], \[CapitalDelta]\[Phi]r, \[CapitalDelta]\[Phi]\[Theta]} = {"\[CapitalDelta]tr", "\[CapitalDelta]t\[Theta]", "\[CapitalDelta]\[Phi]r", "\[CapitalDelta]\[Phi]\[Theta]"} /. TS["Orbit"]["TrajectoryDeltas"];
  {qt0, qr0, q\[Theta]0, q\[Phi]0} = TS["Orbit"]["InitialPhases"];

  rq[qr_]:=rpi[(qr-qr0)/\[CapitalUpsilon]r];
  \[Theta]q[q\[Theta]_]:=\[Theta]pi[(q\[Theta]-q\[Theta]0)/\[CapitalUpsilon]\[Theta]];

  \[Omega] = R["In"]["\[Omega]"];
  l = R["In"]["l"];
  m = R["In"]["m"];

  W = 2 I \[Omega] R["In"]["Amplitudes"]["Incidence"];

  II1In = Function[qr, Evaluate[2TS["\[Alpha]1"][rq[qr]]R["In"][rq[qr]]Cos[(\[Omega] \[CapitalDelta]tr[qr] - m \[CapitalDelta]\[Phi]r[qr] + n qr)]]];
  II1Up = Function[qr, Evaluate[2TS["\[Alpha]1"][rq[qr]]R["Up"][rq[qr]]Cos[(\[Omega] \[CapitalDelta]tr[qr] - m \[CapitalDelta]\[Phi]r[qr] + n qr)]]];
  II2 = Function[q\[Theta], Evaluate[2TS["\[Alpha]2"][\[Theta]q[q\[Theta]]]SH[\[Theta]q[q\[Theta]], 0]Cos[(\[Omega] \[CapitalDelta]t\[Theta][q\[Theta]] - m \[CapitalDelta]\[Phi]\[Theta][q\[Theta]] + k q\[Theta])]]];
  \[Alpha]1In = 1/(2\[Pi]) Quiet[NIntegrate[II1In[qr], {qr,0,\[Pi]}, Method -> {"Trapezoidal", "SymbolicProcessing"->0}, WorkingPrecision -> Precision[\[Omega]]], NIntegrate::precw];
  \[Alpha]1Up = 1/(2\[Pi]) Quiet[NIntegrate[II1Up[qr], {qr,0,\[Pi]}, Method -> {"Trapezoidal", "SymbolicProcessing"->0}, WorkingPrecision -> Precision[\[Omega]]], NIntegrate::precw];
  \[Alpha]2 = 1/(2\[Pi]) Quiet[NIntegrate[II2[q\[Theta]], {q\[Theta],0,\[Pi]}, Method -> {"Trapezoidal", "SymbolicProcessing"->0}, WorkingPrecision -> Precision[\[Omega]]], NIntegrate::precw];

  If[a == 0,
    \[Alpha]3In = \[Alpha]3Up = \[Alpha]4 = 0;,
    II3In = Function[qr, Evaluate[2TS["\[Alpha]3"][rq[qr]]R["In"][rq[qr]]Cos[(\[Omega] \[CapitalDelta]tr[qr] - m \[CapitalDelta]\[Phi]r[qr] + n qr)]]];
    II3Up = Function[qr, Evaluate[2TS["\[Alpha]3"][rq[qr]]R["Up"][rq[qr]]Cos[(\[Omega] \[CapitalDelta]tr[qr] - m \[CapitalDelta]\[Phi]r[qr] + n qr)]]];
    II4 = Function[q\[Theta], Evaluate[2TS["\[Alpha]4"][\[Theta]q[q\[Theta]]]SH[\[Theta]q[q\[Theta]], 0]Cos[(\[Omega] \[CapitalDelta]t\[Theta][q\[Theta]] - m \[CapitalDelta]\[Phi]\[Theta][q\[Theta]] + k q\[Theta])]]];
    \[Alpha]3In = 1/(2\[Pi]) Quiet[NIntegrate[II3In[qr], {qr,0,\[Pi]}, Method -> {"Trapezoidal", "SymbolicProcessing"->0}, WorkingPrecision -> Precision[\[Omega]]], NIntegrate::precw];
    \[Alpha]3Up = 1/(2\[Pi]) Quiet[NIntegrate[II3Up[qr], {qr,0,\[Pi]}, Method -> {"Trapezoidal", "SymbolicProcessing"->0}, WorkingPrecision -> Precision[\[Omega]]], NIntegrate::precw];
    \[Alpha]4 = 1/(2\[Pi]) Quiet[NIntegrate[II4[q\[Theta]], {q\[Theta],0,\[Pi]}, Method -> {"Trapezoidal", "SymbolicProcessing"->0}, WorkingPrecision -> Precision[\[Omega]]], NIntegrate::precw];
  ];

  \[Xi]=m(\[CapitalDelta]\[Phi]r[qr0]+\[CapitalDelta]\[Phi]\[Theta][q\[Theta]0]-q\[Phi]0) - \[Omega](\[CapitalDelta]tr[qr0]+\[CapitalDelta]t\[Theta][q\[Theta]0]-qt0) - k q\[Theta]0 - n qr0;

  ZIn = (\[Alpha]1Up*\[Alpha]2 + \[Alpha]3Up*\[Alpha]4)/W/\[CapitalUpsilon]t Exp[I \[Xi]];
  ZUp = (\[Alpha]1In*\[Alpha]2 + \[Alpha]3In*\[Alpha]4)/W/\[CapitalUpsilon]t Exp[I \[Xi]];

  Clear[a, p, rpi, \[Theta]pi, rq, \[Theta]q, qt0, qr0, q\[Theta]0, q\[Phi]0, \[Xi], W, \[Alpha]1In, \[Alpha]1Up, \[Alpha]2, \[Alpha]3In, \[Alpha]3Up, \[Alpha]4, S, l, m, \[Omega], \[CapitalUpsilon]t, \[CapitalUpsilon]r, \[CapitalUpsilon]\[Theta], \[CapitalDelta]tr, \[CapitalDelta]t\[Theta], \[CapitalDelta]\[Phi]r, \[CapitalDelta]\[Phi]\[Theta], II1In, II1Up, II2, II3In, II3Up, II4];
  <| "\[ScriptCapitalI]" -> ZUp, "\[ScriptCapitalH]" -> ZIn |>
]


(* ::Section::Closed:: *)
(*Spectral Integration*)


(* ::Subsection::Closed:: *)
(*1D Fourier*)


SpectralIntegration1D[integrand_,pgTemp_,print_:False]:=
Module[{err,halfSampleInit = 2^3,pgSum,q,halfSample,sum=0,sumTerm=0,maxTerm=0,deltaQ,intSum,precisionLoss,pg,pgAdjusted,intSumCompare,errOld=0},
	
	q = 0;
	sumTerm=integrand[q];
	sum+=sumTerm;
	If[Abs[sumTerm]>maxTerm, maxTerm= Abs[sumTerm]];
	pgSum=Precision[sum];
	If[pgSum<pgTemp,pg=0.9pgSum,pg=pgTemp];
	
	halfSample = halfSampleInit;
	deltaQ = N[Pi/halfSample,1.2pgSum];

	q = Pi;
	sumTerm=integrand[q];
	sum+=sumTerm;
	If[Abs[sumTerm]>maxTerm, maxTerm= Abs[sumTerm]];

	Do[
		q = i*deltaQ;
		sumTerm=integrand[q];
		sum+=sumTerm;
		If[Abs[sumTerm]>maxTerm, maxTerm= Abs[sumTerm]];

		sumTerm=integrand[2Pi - q];
		sum+=sumTerm;
		If[Abs[sumTerm]>maxTerm, maxTerm= Abs[sumTerm]],
	{i,1,halfSample-1}
	];

	intSum = sum/(2 halfSample);

	If[Abs[intSum]>0,
		precisionLoss = Abs@RealExponent[maxTerm/intSum],
		precisionLoss = 0
	];
	
	pgAdjusted = pg - 1 - precisionLoss;
	If[pgAdjusted < pg, pgAdjusted = pg];

	intSumCompare = 0;
	err = Abs[intSum-intSumCompare];
	If[err>0, err = Abs@RealExponent[err/intSum], err = Infinity];
	While[err < pgAdjusted,
		If[print,
		Print[err];
		Print[intSum];
		];
		Do[
			q = (i+1/2)*deltaQ;
			sumTerm=integrand[q];
			sum+=sumTerm;
			If[Abs[sumTerm]>maxTerm, maxTerm= Abs[sumTerm]];

			sumTerm=integrand[2Pi - q];
			sum+=sumTerm;
			If[Abs[sumTerm]>maxTerm, maxTerm= Abs[sumTerm]],
		{i,0,halfSample-1}
		];

		halfSample *= 2;
		deltaQ *= 1/2;
		intSumCompare=intSum;
		intSum = sum/(2 halfSample);

		If[Abs[intSum]>0,
			precisionLoss = Abs@RealExponent[maxTerm/intSum],
			precisionLoss = 0
		];
		pgAdjusted = pg - 1 - precisionLoss;
		If[pgAdjusted < pg, pgAdjusted = pg];
		errOld = err;
		err = Abs[intSum-intSumCompare];
		If[err>0, err = Abs@RealExponent[err/intSum], err = Infinity];
		If[(Abs[err - errOld] < 1)&&(halfSample > 2^12), err = Infinity] (*If convergence is stagnating and you've tried more than 2000 points, then stop the integration *)
	];
	If[print,
	Print["Done"];
	Print[err];
	];

	intSum
];


(* ::Subsection::Closed:: *)
(*1D DCT*)


SpectralCosineIntegration1D[integrand_,pgTemp_]:=
Module[{err,halfSampleInit = 2^3,pgSum,q,halfSample,sum=0,sumTerm=0,maxTerm=0,deltaQ,intSum,precisionLoss,pg,pgAdjusted,intSumCompare,errOld=0},
	q = 0;
	sumTerm=integrand[q]/2;
	sum+=sumTerm;
	If[Abs[sumTerm]>maxTerm, maxTerm= Abs[sumTerm]];
	pgSum=Precision[sum];
	If[pgSum<pgTemp,pg=0.9pgSum,pg=pgTemp];
	
	halfSample = halfSampleInit;
	deltaQ = N[Pi/halfSample,2pgSum];

	q = Pi;
	sumTerm=integrand[q]/2;
	sum+=sumTerm;
	If[Abs[sumTerm]>maxTerm, maxTerm= Abs[sumTerm]];

	Do[
		q = i*deltaQ;
		sumTerm=integrand[q];
		sum+=sumTerm;
		If[Abs[sumTerm]>maxTerm, maxTerm= Abs[sumTerm]],
	{i,1,halfSample-1}
	];

	intSum = sum/(2 halfSample);

	If[Abs[intSum]>0,
		precisionLoss = Abs@RealExponent[maxTerm/intSum],
		precisionLoss = 0
	];
	
	pgAdjusted = pg - 1 - precisionLoss;
	If[pgAdjusted < pg, pgAdjusted = pg];

	intSumCompare = 0;
	err = Abs[intSum-intSumCompare];
	If[err>0, err = Abs@RealExponent[err/intSum], err = Infinity];
	While[err < pgAdjusted,
		Do[
			q = (i+1/2)*deltaQ;
			sumTerm=integrand[q];
			sum+=sumTerm;
			If[Abs[sumTerm]>maxTerm, maxTerm= Abs[sumTerm]],
		{i,0,halfSample-1}
		];

		halfSample *= 2;
		deltaQ *= 1/2;
		intSumCompare=intSum;
		intSum = sum/(2 halfSample);

		If[Abs[intSum]>0,
			precisionLoss = Abs@RealExponent[maxTerm/intSum],
			precisionLoss = 0
		];
		pgAdjusted = pg - 1 - precisionLoss;
		If[pgAdjusted < pg, pgAdjusted = pg];
		errOld = err;
		err = Abs[intSum-intSumCompare];
		If[err>0, err = Abs@RealExponent[err/intSum], err = Infinity];
		If[(Abs[err - errOld] < 1)&&(halfSample > 2^12), err = Infinity] (*If convergence is stagnating and you've tried more than 2000 points, then stop the integration *)
	];

	intSum
];


(* ::Subsection::Closed:: *)
(*2D Fourier*)


SpectralIntegration2D[TS_,pgTemp_]:=
Module[{TS2},
	TS2=Function[q\[Theta],SpectralIntegration1D[Function[qr,TS[qr,q\[Theta]]],pgTemp]];
	SpectralIntegration1D[TS2,pgTemp]
];


(* ::Subsection::Closed:: *)
(*2D DCT*)


SpectralCosineIntegration2D[TS_,pgTemp_]:=
Module[{TS2},
	TS2=Function[q\[Theta],SpectralCosineIntegration1D[Function[qr,TS[qr,q\[Theta]]],pgTemp]];
	SpectralCosineIntegration1D[TS2,pgTemp]
];


(* ::Section::Closed:: *)
(*End Package*)


End[];
EndPackage[];
