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
 Module[{a, p, \[ScriptCapitalE], \[ScriptCapitalL], r0, \[Theta]0, \[CapitalDelta], Kt, \[CapitalUpsilon]t, m, \[Omega], \[Lambda], W, RIn, dRIn, d2RIn, RUp, dRUp, d2RUp, S0, dS0, d2S0, L1, L2, L2S, L2p, L1Sp, L1L2S, rcomp, \[Theta]comp, \[Rho], \[Rho]bar, \[CapitalSigma], Ann0, Anmbar0, Anmbar1, Ambarmbar0, Ambarmbar1, Ambarmbar2, Cnnp1p1, Cnmbarp1p1, Cmbarmbarp1p1, \[Alpha]In, \[Alpha]Up, ZIn, ZUp},
  a = TS["Orbit"]["a"];
  p = TS["Orbit"]["p"];
  {\[ScriptCapitalE], \[ScriptCapitalL]} = TS["Orbit"] /@ {"Energy", "AngularMomentum"};

  \[CapitalUpsilon]t = TS["Orbit"]["Frequencies"]["\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(t\)]\)"];

  m = R["In"]["m"];
  \[Omega] = R["In"]["\[Omega]"];
  \[Lambda] = R["In"]["Eigenvalue"];

  W = 2 I \[Omega] R["In"]["Amplitudes"]["Incidence"];

  r0 = p;
  \[Theta]0 = \[Pi]/2;

  \[CapitalDelta] = r0^2-2r0+a^2;
  Kt=(r0^2+a^2)\[Omega]-m a;

  RUp = R["Up"][r0];
  RIn = R["In"][r0];
  dRUp = R["Up"]'[r0];
  dRIn = R["In"]'[r0];
  d2RUp = (-(-\[Lambda] + 2 I r0 s 2 \[Omega] + (-2 I (-1 + r0) s (-a m + (a^2 + r0^2) \[Omega]) + (-a m + (a^2 + r0^2) \[Omega])^2)/(a^2 - 2 r0 + r0^2)) RUp - (-2 + 2 r0) (1 + s) dRUp)/(a^2 - 2 r0 + r0^2);
  d2RIn = (-(-\[Lambda] + 2 I r0 s 2 \[Omega] + (-2 I (-1 + r0) s (-a m + (a^2 + r0^2) \[Omega]) + (-a m + (a^2 + r0^2) \[Omega])^2)/(a^2 - 2 r0 + r0^2)) RIn - (-2 + 2 r0) (1 + s) dRIn)/(a^2 - 2 r0 + r0^2);

  S0 = SH[\[Theta]0, 0];
  dS0 = Derivative[1,0][SH][\[Theta]0, 0];
  d2S0 = Derivative[2,0][SH][\[Theta]0, 0];
  L1 = -m/Sin[\[Theta]0] + a \[Omega] Sin[\[Theta]0] + Cos[\[Theta]0]/Sin[\[Theta]0];
  L2 = -m/Sin[\[Theta]0] + a \[Omega] Sin[\[Theta]0] + 2 Cos[\[Theta]0]/Sin[\[Theta]0];
  L2S = dS0 + L2 S0;
  L2p = m Cos[\[Theta]0]/Sin[\[Theta]0]^2 + a \[Omega] Cos[\[Theta]0] - 2/Sin[\[Theta]0]^2;
  L1Sp = d2S0 + L1 dS0;
  L1L2S = L1Sp + L2p S0 + L2 dS0 + L1 L2 S0;

  \[Rho] = -1/(r0 - I a Cos[\[Theta]0]);
  \[Rho]bar = -1/(r0 + I a Cos[\[Theta]0]);
  \[CapitalSigma] = 1/(\[Rho] \[Rho]bar);

  Ann0 = -\[Rho]^(-2) \[Rho]bar^(-1) (Sqrt[2] \[CapitalDelta])^(-2) (\[Rho]^(-1) L1L2S + 3 I a Sin[\[Theta]0] L1 S0 + 3 I a Cos[\[Theta]0] S0 + 2 I a Sin[\[Theta]0] dS0 - I a Sin[\[Theta]0] L2 S0 );
  Anmbar0 = \[Rho]^(-3) (Sqrt[2]\[CapitalDelta])^(-1) ( (\[Rho] + \[Rho]bar - I Kt/\[CapitalDelta]) L2S + (\[Rho] - \[Rho]bar) a Sin[\[Theta]0] Kt/\[CapitalDelta] S0 );
  Anmbar1 = -\[Rho]^(-3) (Sqrt[2]\[CapitalDelta])^(-1) ( L2S + I (\[Rho] - \[Rho]bar) a Sin[\[Theta]0] S0 );
  Ambarmbar0 = \[Rho]^(-3) \[Rho]bar S0 Kt/\[CapitalDelta]/4 ( I (2 \[Omega] r0/Kt - 2(r0 - 1)/\[CapitalDelta]) + Kt/\[CapitalDelta] + 2 I \[Rho]);
  Ambarmbar1 = -\[Rho]^(-3) \[Rho]bar S0/2 ( I Kt/\[CapitalDelta] - \[Rho] );
  Ambarmbar2 = -\[Rho]^(-3) \[Rho]bar S0/4;

  rcomp = (\[ScriptCapitalE](r0^2+a^2) - a \[ScriptCapitalL])/(2\[CapitalSigma]);
  \[Theta]comp = \[Rho] (I Sin[\[Theta]0](a \[ScriptCapitalE] - \[ScriptCapitalL]/Sin[\[Theta]0]^2))/Sqrt[2];
    
  {Cnnp1p1, Cnmbarp1p1, Cmbarmbarp1p1} = {rcomp^2, rcomp \[Theta]comp, \[Theta]comp^2};
    
  \[Alpha]In = ((Ann0*Cnnp1p1 + Anmbar0*Cnmbarp1p1 + Ambarmbar0*Cmbarmbarp1p1) RIn-(Anmbar1*Cnmbarp1p1 + Ambarmbar1*Cmbarmbarp1p1) dRIn+Ambarmbar2*Cmbarmbarp1p1 d2RIn);
  \[Alpha]Up = ((Ann0*Cnnp1p1 + Anmbar0*Cnmbarp1p1 + Ambarmbar0*Cmbarmbarp1p1) RUp-(Anmbar1*Cnmbarp1p1 + Ambarmbar1*Cmbarmbarp1p1) dRUp+Ambarmbar2*Cmbarmbarp1p1 d2RUp);

  ZIn = 8Pi \[Alpha]Up/W/\[CapitalUpsilon]t;
  ZUp = 8Pi \[Alpha]In/W/\[CapitalUpsilon]t;

  Clear[a, p, \[ScriptCapitalE], \[ScriptCapitalL], r0, \[Theta]0, \[CapitalDelta], Kt, \[CapitalUpsilon]t, m, \[Omega], \[Lambda], W, RIn, dRIn, d2RIn, RUp, dRUp, d2RUp, S0, dS0, d2S0, L1, L2, L2S, L2p, L1Sp, L1L2S, rcomp, \[Theta]comp, \[Rho], \[Rho]bar, \[CapitalSigma], Ann0, Anmbar0, Anmbar1, Ambarmbar0, Ambarmbar1, Ambarmbar2, Cnnp1p1, Cnmbarp1p1, Cmbarmbarp1p1, \[Alpha]In, \[Alpha]Up];
  <| "\[ScriptCapitalI]" -> ZUp, "\[ScriptCapitalH]" -> ZIn |>
]


(* ::Subsection::Closed:: *)
(*s=-2 point particle on a spherical orbit*)


ConvolveSourcePointParticleSpherical[s:-2, k_Integer, R_, SH_, TS_] :=
 Module[{a, p, \[ScriptCapitalE], \[ScriptCapitalL], \[Theta]pi, u\[Theta]pi, r0, ur0, \[CapitalDelta], Kt, \[CapitalUpsilon]t, \[CapitalUpsilon]r, \[CapitalUpsilon]\[Theta], \[CapitalDelta]t\[Theta], \[CapitalDelta]\[Phi]\[Theta], qt0, qr0, q\[Theta]0, q\[Phi]0, m, \[Omega], \[Lambda], W, \[Theta]q, u\[Theta]q, RIn, dRIn, d2RIn, RUp, dRUp, d2RUp, integrand, \[Alpha]In, \[Alpha]Up, \[Xi], wpIn, wpUp, ZIn, ZUp},
  a = TS["Orbit"]["a"];
  p = TS["Orbit"]["p"];
  {\[ScriptCapitalE], \[ScriptCapitalL]} = TS["Orbit"] /@ {"Energy", "AngularMomentum"};
  \[Theta]pi = TS["Orbit"]["Trajectory"][[3]];
  u\[Theta]pi = TS["Orbit"]["FourVelocity"][[3]];

  {\[CapitalUpsilon]t, \[CapitalUpsilon]\[Theta]} = {"\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(t\)]\)", "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Theta]\)]\)"} /. TS["Orbit"]["Frequencies"];
  {\[CapitalDelta]t\[Theta], \[CapitalDelta]\[Phi]\[Theta]} = {"\[CapitalDelta]t\[Theta]", "\[CapitalDelta]\[Phi]\[Theta]"} /. TS["Orbit"]["TrajectoryDeltas"];
  {qt0, qr0, q\[Theta]0, q\[Phi]0} = TS["Orbit"]["InitialPhases"];

  m = R["In"]["m"];
  \[Omega] = R["In"]["\[Omega]"];
  \[Lambda] = R["In"]["Eigenvalue"];

  W = 2 I \[Omega] R["In"]["Amplitudes"]["Incidence"];

  r0 = p;
  \[Theta]q[q\[Theta]_] := \[Theta]pi[(q\[Theta]-q\[Theta]0)/\[CapitalUpsilon]\[Theta]];

  (* Mino-time four-velocities *)
  ur0 = 0;
  u\[Theta]q[q\[Theta]_] := (r0^2+a^2 Cos[\[Theta]pi[(q\[Theta]-q\[Theta]0)/\[CapitalUpsilon]\[Theta]]]^2) u\[Theta]pi[(q\[Theta]-q\[Theta]0)/\[CapitalUpsilon]\[Theta]];

  \[CapitalDelta] = r0^2-2r0+a^2;
  Kt=(r0^2+a^2)\[Omega]-m a;

  RUp = R["Up"][r0];
  RIn = R["In"][r0];
  dRUp = R["Up"]'[r0];
  dRIn = R["In"]'[r0];
  d2RUp = (-(-\[Lambda] + 2 I r0 s 2 \[Omega] + (-2 I (-1 + r0) s (-a m + (a^2 + r0^2) \[Omega]) + (-a m + (a^2 + r0^2) \[Omega])^2)/(a^2 - 2 r0 + r0^2)) RUp - (-2 + 2 r0) (1 + s) dRUp)/(a^2 - 2 r0 + r0^2);
  d2RIn = (-(-\[Lambda] + 2 I r0 s 2 \[Omega] + (-2 I (-1 + r0) s (-a m + (a^2 + r0^2) \[Omega]) + (-a m + (a^2 + r0^2) \[Omega])^2)/(a^2 - 2 r0 + r0^2)) RIn - (-2 + 2 r0) (1 + s) dRIn)/(a^2 - 2 r0 + r0^2);

  integrand[q\[Theta]_, {R0_, dR0_, d2R0_}]:=
   Module[{\[Theta]0, S0, dS0, d2S0, L1, L2, L2S, L2p, L1Sp, L1L2S, u\[Theta]0, rcomp, \[Theta]comp, \[Rho], \[Rho]bar, \[CapitalSigma], Ann0, Anmbar0, Anmbar1, Ambarmbar0, Ambarmbar1, Ambarmbar2, Cnnp1p1, Cnmbarp1p1, Cmbarmbarp1p1, Cnnp1m1, Cnmbarp1m1, Cmbarmbarp1m1, \[Theta]phase, res},
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
  
    \[Rho] = -1/(r0 - I a Cos[\[Theta]0]);
    \[Rho]bar = -1/(r0 + I a Cos[\[Theta]0]);
    \[CapitalSigma] = 1/(\[Rho] \[Rho]bar);
  
    Ann0 = -\[Rho]^(-2) \[Rho]bar^(-1) (Sqrt[2] \[CapitalDelta])^(-2) (\[Rho]^(-1) L1L2S + 3 I a Sin[\[Theta]0] L1 S0 + 3 I a Cos[\[Theta]0] S0 + 2 I a Sin[\[Theta]0] dS0 - I a Sin[\[Theta]0] L2 S0 );
    Anmbar0 = \[Rho]^(-3) (Sqrt[2]\[CapitalDelta])^(-1) ( (\[Rho] + \[Rho]bar - I Kt/\[CapitalDelta]) L2S + (\[Rho] - \[Rho]bar) a Sin[\[Theta]0] Kt/\[CapitalDelta] S0 );
    Anmbar1 = -\[Rho]^(-3) (Sqrt[2]\[CapitalDelta])^(-1) ( L2S + I (\[Rho] - \[Rho]bar) a Sin[\[Theta]0] S0 );
    Ambarmbar0 = \[Rho]^(-3) \[Rho]bar S0 Kt/\[CapitalDelta]/4 ( I (2 \[Omega] r0/Kt - 2(r0 - 1)/\[CapitalDelta]) + Kt/\[CapitalDelta] + 2 I \[Rho]);
    Ambarmbar1 = -\[Rho]^(-3) \[Rho]bar S0/2 ( I Kt/\[CapitalDelta] - \[Rho] );
    Ambarmbar2 = -\[Rho]^(-3) \[Rho]bar S0/4;

    (* Save time by folding the two segments in Subscript[q, \[Theta]]\[Element][0,2\[Pi]] over to Subscript[q, \[Theta]]\[Element][0,\[Pi]] *)
    rcomp = (\[ScriptCapitalE](r0^2+a^2) - a \[ScriptCapitalL] + ur0)/(2\[CapitalSigma]);
    \[Theta]comp = \[Rho] (I Sin[\[Theta]0](a \[ScriptCapitalE] - \[ScriptCapitalL]/Sin[\[Theta]0]^2) + u\[Theta]0)/Sqrt[2];
    
    {{Cnnp1p1,Cnmbarp1p1,Cmbarmbarp1p1}, {Cnnp1m1,Cnmbarp1m1,Cmbarmbarp1m1}} =
      Table[{rcomp^2, rcomp \[Theta]comp, \[Theta]comp^2}, {u\[Theta]0, {u\[Theta]q[q\[Theta]], u\[Theta]q[2\[Pi]-q\[Theta]]}}];
    
    \[Theta]phase = \[Omega] \[CapitalDelta]t\[Theta][q\[Theta]] - m \[CapitalDelta]\[Phi]\[Theta][q\[Theta]] + k q\[Theta];
    
    res = ((Ann0*Cnnp1p1 + Anmbar0*Cnmbarp1p1 + Ambarmbar0*Cmbarmbarp1p1) R0-(Anmbar1*Cnmbarp1p1 + Ambarmbar1*Cmbarmbarp1p1) dR0+Ambarmbar2*Cmbarmbarp1p1 d2R0)Exp[I \[Theta]phase] 
        + ((Ann0*Cnnp1m1 + Anmbar0*Cnmbarp1m1 + Ambarmbar0*Cmbarmbarp1m1) R0-(Anmbar1*Cnmbarp1m1 + Ambarmbar1*Cmbarmbarp1m1) dR0+Ambarmbar2*Cmbarmbarp1m1 d2R0)Exp[-I \[Theta]phase];
    Clear[\[Theta]0, S0, dS0, d2S0, L1, L2, L2S, L2p, L1Sp, L1L2S, u\[Theta]0, rcomp, \[Theta]comp, \[Rho], \[Rho]bar, \[CapitalSigma], Ann0, Anmbar0, Anmbar1, Ambarmbar0, Ambarmbar1, Ambarmbar2, Cnnp1p1, Cnmbarp1p1, Cmbarmbarp1p1, Cnnp1m1, Cnmbarp1m1, Cmbarmbarp1m1, \[Theta]phase];
    res  
  ];

  wpIn = Precision[integrand[0, {RIn, dRIn, d2RIn}]];
  wpUp = Precision[integrand[0, {RUp, dRUp, d2RUp}]];

  \[Alpha]In = 1/(2\[Pi]) Quiet[NIntegrate[integrand[q\[Theta], {RIn, dRIn, d2RIn}], {q\[Theta], 0, \[Pi]}, Method -> {"Trapezoidal", "SymbolicProcessing"->0}, WorkingPrecision -> wpIn], NIntegrate::precw];
  \[Alpha]Up = 1/(2\[Pi]) Quiet[NIntegrate[integrand[q\[Theta], {RUp, dRUp, d2RUp}], {q\[Theta], 0, \[Pi]}, Method -> {"Trapezoidal", "SymbolicProcessing"->0}, WorkingPrecision -> wpUp], NIntegrate::precw];
  
  \[Xi] = m(\[CapitalDelta]\[Phi]\[Theta][q\[Theta]0]-q\[Phi]0) - \[Omega](\[CapitalDelta]t\[Theta][q\[Theta]0]-qt0) - k q\[Theta]0;

  ZIn = 8Pi \[Alpha]Up/W/\[CapitalUpsilon]t Exp[I \[Xi]];
  ZUp = 8Pi \[Alpha]In/W/\[CapitalUpsilon]t Exp[I \[Xi]];

  Clear[a, p, \[ScriptCapitalE], \[ScriptCapitalL], \[Theta]pi, u\[Theta]pi, r0, ur0, \[CapitalDelta], Kt, \[CapitalUpsilon]t, \[CapitalUpsilon]r, \[CapitalUpsilon]\[Theta], \[CapitalDelta]t\[Theta], \[CapitalDelta]\[Phi]\[Theta], qt0, qr0, q\[Theta]0, q\[Phi]0, m, \[Omega], \[Lambda], W, \[Theta]q, u\[Theta]q, RIn, dRIn, d2RIn, RUp, dRUp, d2RUp, integrand, \[Alpha]In, \[Alpha]Up, \[Xi], wpIn, wpUp];
  <| "\[ScriptCapitalI]" -> ZUp, "\[ScriptCapitalH]" -> ZIn |>
]


(* ::Subsection::Closed:: *)
(*s=-2 point particle on an eccentric orbit*)


ConvolveSourcePointParticleEccentric[s_:-2, n_Integer, R_, SH_, TS_] :=
 Module[{a, p, \[ScriptCapitalE], \[ScriptCapitalL], rpi, urpi, \[CapitalUpsilon]t, \[CapitalUpsilon]r, \[CapitalDelta]tr, \[CapitalDelta]\[Phi]r, qt0, qr0, q\[Theta]0, q\[Phi]0, m, \[Omega], \[Lambda], W, rq, \[Theta]0, urq, u\[Theta]0, S0, dS0, d2S0, L1, L2, L2S, L2p, L1Sp, L1L2S, integrand, \[Alpha]In, \[Alpha]Up, \[Xi], wpIn, wpUp, ZIn, ZUp},
  a = TS["Orbit"]["a"];
  p = TS["Orbit"]["p"];
  {\[ScriptCapitalE], \[ScriptCapitalL]} = TS["Orbit"] /@ {"Energy", "AngularMomentum"};
  rpi = TS["Orbit"]["Trajectory"][[2]];
  urpi = TS["Orbit"]["FourVelocity"][[2]];

  {\[CapitalUpsilon]t, \[CapitalUpsilon]r} = {"\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(t\)]\)", "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(r\)]\)"} /. TS["Orbit"]["Frequencies"];
  {\[CapitalDelta]tr, \[CapitalDelta]\[Phi]r} = {"\[CapitalDelta]tr", "\[CapitalDelta]\[Phi]r"} /. TS["Orbit"]["TrajectoryDeltas"];
  {qt0, qr0, q\[Theta]0, q\[Phi]0} = TS["Orbit"]["InitialPhases"];

  m = R["In"]["m"];
  \[Omega] = R["In"]["\[Omega]"];
  \[Lambda] = R["In"]["Eigenvalue"];

  W = 2 I \[Omega] R["In"]["Amplitudes"]["Incidence"];
  
  rq[qr_] := rpi[(qr-qr0)/\[CapitalUpsilon]r];
  \[Theta]0 = Pi/2;

  (* Mino-time four-velocity *)
  urq[qr_] := (rpi[(qr-qr0)/\[CapitalUpsilon]r]^2+a^2 Cos[\[Theta]0]^2) urpi[(qr-qr0)/\[CapitalUpsilon]r];
  u\[Theta]0 = 0;

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
   Module[{r0, R0, dR0, d2R0, ur0, rcomp, \[Theta]comp, \[CapitalDelta], Kt, \[Rho], \[Rho]bar, \[CapitalSigma], Ann0, Anmbar0, Anmbar1, Ambarmbar0, Ambarmbar1, Ambarmbar2, Cnnp1p1, Cnmbarp1p1, Cmbarmbarp1p1, Cnnm1p1, Cnmbarm1p1, Cmbarmbarm1p1, rphase, res},
    r0 = rq[qr];
    R0   = RF[r0];
    dR0  = RF'[r0];
    d2R0 = (-(-\[Lambda] + 2 I r0 s 2 \[Omega] + (-2 I (-1 + r0) s (-a m + (a^2 + r0^2) \[Omega]) + (-a m + (a^2 + r0^2) \[Omega])^2)/(a^2 - 2 r0 + r0^2)) R0 - (-2 + 2 r0) (1 + s) dR0)/(a^2 - 2 r0 + r0^2);

    \[CapitalDelta] = r0^2 + a^2 - 2 r0;
    Kt = (r0^2 + a^2) \[Omega] - m a;
    \[Rho] = -1/(r0 - I a Cos[\[Theta]0]);
    \[Rho]bar = -1/(r0 + I a Cos[\[Theta]0]);
    \[CapitalSigma] = 1/(\[Rho] \[Rho]bar);
  
    Ann0 = -\[Rho]^(-2) \[Rho]bar^(-1) (Sqrt[2] \[CapitalDelta])^(-2) (\[Rho]^(-1) L1L2S + 3 I a Sin[\[Theta]0] L1 S0 + 3 I a Cos[\[Theta]0] S0 + 2 I a Sin[\[Theta]0] dS0 - I a Sin[\[Theta]0] L2 S0 );
    Anmbar0 = \[Rho]^(-3) (Sqrt[2]\[CapitalDelta])^(-1) ( (\[Rho] + \[Rho]bar - I Kt/\[CapitalDelta]) L2S + (\[Rho] - \[Rho]bar) a Sin[\[Theta]0] Kt/\[CapitalDelta] S0 );
    Anmbar1 = -\[Rho]^(-3) (Sqrt[2]\[CapitalDelta])^(-1) ( L2S + I (\[Rho] - \[Rho]bar) a Sin[\[Theta]0] S0 );
    Ambarmbar0 = \[Rho]^(-3) \[Rho]bar S0 Kt/\[CapitalDelta]/4 ( I (2 \[Omega] r0/Kt - 2(r0 - 1)/\[CapitalDelta]) + Kt/\[CapitalDelta] + 2 I \[Rho]);
    Ambarmbar1 = -\[Rho]^(-3) \[Rho]bar S0/2 ( I Kt/\[CapitalDelta] - \[Rho] );
    Ambarmbar2 = -\[Rho]^(-3) \[Rho]bar S0/4;

    (* Save time by folding the two segments in Subscript[q, r]\[Element][0,2\[Pi]] over to Subscript[q, r]\[Element][0,\[Pi]] *)
    rcomp = (\[ScriptCapitalE](r0^2+a^2) - a \[ScriptCapitalL] + ur0)/(2\[CapitalSigma]);
    \[Theta]comp = \[Rho] (I Sin[\[Theta]0](a \[ScriptCapitalE] - \[ScriptCapitalL]/Sin[\[Theta]0]^2) + u\[Theta]0)/Sqrt[2];

    {{Cnnp1p1,Cnmbarp1p1,Cmbarmbarp1p1}, {Cnnm1p1,Cnmbarm1p1,Cmbarmbarm1p1}} =
      Table[{rcomp^2, rcomp \[Theta]comp, \[Theta]comp^2}, {ur0, {urq[qr], urq[2\[Pi]-qr]}}];
    
    rphase = \[Omega] \[CapitalDelta]tr[qr] - m \[CapitalDelta]\[Phi]r[qr] + n qr;
    
    res = ((Ann0 Cnnp1p1 + Anmbar0 Cnmbarp1p1 + Ambarmbar0 Cmbarmbarp1p1) R0-(Anmbar1 Cnmbarp1p1 + Ambarmbar1 Cmbarmbarp1p1) dR0+Ambarmbar2 Cmbarmbarp1p1 d2R0)Exp[I rphase] 
        + ((Ann0 Cnnm1p1 + Anmbar0 Cnmbarm1p1 + Ambarmbar0 Cmbarmbarm1p1) R0-(Anmbar1 Cnmbarm1p1 + Ambarmbar1 Cmbarmbarm1p1) dR0+Ambarmbar2 Cmbarmbarm1p1 d2R0)Exp[-I rphase];

    Clear[r0, R0, dR0, d2R0, ur0, rcomp, \[Theta]comp, \[CapitalDelta], Kt, \[Rho], \[Rho]bar, \[CapitalSigma], Ann0, Anmbar0, Anmbar1, Ambarmbar0, Ambarmbar1, Ambarmbar2, Cnnp1p1, Cnmbarp1p1, Cmbarmbarp1p1, Cnnm1p1, Cnmbarm1p1, Cmbarmbarm1p1, rphase];
    res
  ];

  wpIn = Precision[integrand[0,R["In"]]];
  wpUp = Precision[integrand[0,R["Up"]]];

  \[Alpha]In = 1/(2\[Pi]) Quiet[NIntegrate[integrand[qr,R["In"]], {qr, 0, \[Pi]}, Method -> {"Trapezoidal", "SymbolicProcessing"->0}, WorkingPrecision -> wpIn], NIntegrate::precw];
  \[Alpha]Up = 1/(2\[Pi]) Quiet[NIntegrate[integrand[qr,R["Up"]], {qr, 0, \[Pi]}, Method -> {"Trapezoidal", "SymbolicProcessing"->0}, WorkingPrecision -> wpUp], NIntegrate::precw];
  
  \[Xi] = m(\[CapitalDelta]\[Phi]r[qr0]-q\[Phi]0) - \[Omega](\[CapitalDelta]tr[qr0]-qt0) - n qr0;
  
  ZIn = 8Pi \[Alpha]Up/W/\[CapitalUpsilon]t Exp[I \[Xi]];
  ZUp = 8Pi \[Alpha]In/W/\[CapitalUpsilon]t Exp[I \[Xi]];

  Clear[a, p, \[ScriptCapitalE], \[ScriptCapitalL], rpi, urpi, \[CapitalUpsilon]t, \[CapitalUpsilon]r, \[CapitalDelta]tr, \[CapitalDelta]\[Phi]r, qt0, qr0, q\[Theta]0, q\[Phi]0, m, \[Omega], \[Lambda], W, rq, \[Theta]0, urq, u\[Theta]0, S0, dS0, d2S0, L1, L2, L2S, L2p, L1Sp, L1L2S, integrand, \[Alpha]In, \[Alpha]Up, \[Xi], wpIn, wpUp];
  <| "\[ScriptCapitalI]" -> ZUp, "\[ScriptCapitalH]" -> ZIn |>
]


(* ::Subsection::Closed:: *)
(*s=-2 point particle on a generic orbit*)


ConvolveSourcePointParticleGeneric[s_:-2, n_Integer, k_Integer, R_, SH_, TS_] :=
 Module[{a, p, \[ScriptCapitalE], \[ScriptCapitalL], rpi, \[Theta]pi, urpi, u\[Theta]pi, \[CapitalUpsilon]t, \[CapitalUpsilon]r, \[CapitalUpsilon]\[Theta], \[CapitalDelta]tr, \[CapitalDelta]t\[Theta], \[CapitalDelta]\[Phi]r, \[CapitalDelta]\[Phi]\[Theta], qt0, qr0, q\[Theta]0, q\[Phi]0, m, \[Omega], \[Lambda], W, rq, \[Theta]q, urq, u\[Theta]q, integrand, integrandIn, integrandUp, \[Alpha]In, \[Alpha]Up, \[Xi], wpIn, wpUp, ZIn, ZUp},
  a = TS["Orbit"]["a"];
  p = TS["Orbit"]["p"];
  {\[ScriptCapitalE], \[ScriptCapitalL]} = TS["Orbit"] /@ {"Energy", "AngularMomentum"};
  {rpi, \[Theta]pi} = TS["Orbit"]["Trajectory"][[2;;3]];
  {urpi, u\[Theta]pi} = TS["Orbit"]["FourVelocity"][[2;;3]];

  {\[CapitalUpsilon]t, \[CapitalUpsilon]r, \[CapitalUpsilon]\[Theta]} = {"\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(t\)]\)", "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(r\)]\)", "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Theta]\)]\)"} /. TS["Orbit"]["Frequencies"];
  {\[CapitalDelta]tr, \[CapitalDelta]t\[Theta], \[CapitalDelta]\[Phi]r, \[CapitalDelta]\[Phi]\[Theta]} = {"\[CapitalDelta]tr", "\[CapitalDelta]t\[Theta]", "\[CapitalDelta]\[Phi]r", "\[CapitalDelta]\[Phi]\[Theta]"} /. TS["Orbit"]["TrajectoryDeltas"];
  {qt0, qr0, q\[Theta]0, q\[Phi]0} = TS["Orbit"]["InitialPhases"];

  m = R["In"]["m"];
  \[Omega] = R["In"]["\[Omega]"];
  \[Lambda] = R["In"]["Eigenvalue"];

  W = 2 I \[Omega] R["In"]["Amplitudes"]["Incidence"];

  rq[qr_] := rpi[(qr-qr0)/\[CapitalUpsilon]r];
  \[Theta]q[q\[Theta]_] := \[Theta]pi[(q\[Theta]-q\[Theta]0)/\[CapitalUpsilon]\[Theta]];

  (* Mino-time four-velocities *)
  urq[qr_] := (rpi[(qr-qr0)/\[CapitalUpsilon]r]^2+a^2 Cos[\[Theta]pi[(qr-qr0)/\[CapitalUpsilon]r]]^2) urpi[(qr-qr0)/\[CapitalUpsilon]r];
  u\[Theta]q[q\[Theta]_] := (rpi[(q\[Theta]-q\[Theta]0)/\[CapitalUpsilon]\[Theta]]^2+a^2 Cos[\[Theta]pi[(q\[Theta]-q\[Theta]0)/\[CapitalUpsilon]\[Theta]]]^2) u\[Theta]pi[(q\[Theta]-q\[Theta]0)/\[CapitalUpsilon]\[Theta]];

  integrand[qr_, q\[Theta]_, RF_]:=
   Module[{r0, R0, dR0, d2R0, \[Theta]0, S0, dS0, d2S0, L1, L2, L2S, L2p, L1Sp, L1L2S, ur0, u\[Theta]0, rcomp, \[Theta]comp, \[CapitalDelta], Kt, \[Rho], \[Rho]bar, \[CapitalSigma], Ann0, Anmbar0, Anmbar1, Ambarmbar0, Ambarmbar1, Ambarmbar2, Cnnp1p1, Cnmbarp1p1, Cmbarmbarp1p1, Cnnp1m1, Cnmbarp1m1, Cmbarmbarp1m1, Cnnm1p1, Cnmbarm1p1, Cmbarmbarm1p1, Cnnm1m1,Cnmbarm1m1,Cmbarmbarm1m1, rphase, \[Theta]phase, res},
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
    \[CapitalSigma] = 1/(\[Rho] \[Rho]bar);
  
    Ann0 = -\[Rho]^(-2) \[Rho]bar^(-1) (Sqrt[2] \[CapitalDelta])^(-2) (\[Rho]^(-1) L1L2S + 3 I a Sin[\[Theta]0] L1 S0 + 3 I a Cos[\[Theta]0] S0 + 2 I a Sin[\[Theta]0] dS0 - I a Sin[\[Theta]0] L2 S0 );
    Anmbar0 = \[Rho]^(-3) (Sqrt[2]\[CapitalDelta])^(-1) ( (\[Rho] + \[Rho]bar - I Kt/\[CapitalDelta]) L2S + (\[Rho] - \[Rho]bar) a Sin[\[Theta]0] Kt/\[CapitalDelta] S0 );
    Anmbar1 = -\[Rho]^(-3) (Sqrt[2]\[CapitalDelta])^(-1) ( L2S + I (\[Rho] - \[Rho]bar) a Sin[\[Theta]0] S0 );
    Ambarmbar0 = \[Rho]^(-3) \[Rho]bar S0 Kt/\[CapitalDelta]/4 ( I (2 \[Omega] r0/Kt - 2(r0 - 1)/\[CapitalDelta]) + Kt/\[CapitalDelta] + 2 I \[Rho]);
    Ambarmbar1 = -\[Rho]^(-3) \[Rho]bar S0/2 ( I Kt/\[CapitalDelta] - \[Rho] );
    Ambarmbar2 = -\[Rho]^(-3) \[Rho]bar S0/4;

    (* Save time by folding the four segments in Subscript[q, \[Theta]]\[Element][0,2\[Pi]], Subscript[q, r]\[Element][0,2\[Pi]] over to Subscript[q, \[Theta]]\[Element][0,\[Pi]], Subscript[q, r]\[Element][0,\[Pi]] *)
    rcomp = (\[ScriptCapitalE](r0^2+a^2) - a \[ScriptCapitalL] + ur0)/(2\[CapitalSigma]);
    \[Theta]comp = \[Rho] (I Sin[\[Theta]0](a \[ScriptCapitalE] - \[ScriptCapitalL]/Sin[\[Theta]0]^2) + u\[Theta]0)/Sqrt[2];
    
    {{{Cnnp1p1,Cnmbarp1p1,Cmbarmbarp1p1}, {Cnnp1m1,Cnmbarp1m1,Cmbarmbarp1m1}},
      {{Cnnm1p1,Cnmbarm1p1,Cmbarmbarm1p1}, {Cnnm1m1,Cnmbarm1m1,Cmbarmbarm1m1}}} =
      Table[{rcomp^2, rcomp \[Theta]comp, \[Theta]comp^2}, {ur0, {urq[qr], urq[2\[Pi]-qr]}}, {u\[Theta]0, {u\[Theta]q[q\[Theta]], u\[Theta]q[2\[Pi]-q\[Theta]]}}];
    
    rphase = \[Omega] \[CapitalDelta]tr[qr] - m \[CapitalDelta]\[Phi]r[qr] + n qr;
    \[Theta]phase = \[Omega] \[CapitalDelta]t\[Theta][q\[Theta]] - m \[CapitalDelta]\[Phi]\[Theta][q\[Theta]] + k q\[Theta];
    
    res = ((Ann0*Cnnp1p1 + Anmbar0*Cnmbarp1p1 + Ambarmbar0*Cmbarmbarp1p1) R0-(Anmbar1*Cnmbarp1p1 + Ambarmbar1*Cmbarmbarp1p1) dR0+Ambarmbar2*Cmbarmbarp1p1 d2R0)Exp[I(rphase+\[Theta]phase)] 
        + ((Ann0*Cnnp1m1 + Anmbar0*Cnmbarp1m1 + Ambarmbar0*Cmbarmbarp1m1) R0-(Anmbar1*Cnmbarp1m1 + Ambarmbar1*Cmbarmbarp1m1) dR0+Ambarmbar2*Cmbarmbarp1m1 d2R0)Exp[I(rphase-\[Theta]phase)] 
        + ((Ann0*Cnnm1p1 + Anmbar0*Cnmbarm1p1 + Ambarmbar0*Cmbarmbarm1p1) R0-(Anmbar1*Cnmbarm1p1 + Ambarmbar1*Cmbarmbarm1p1) dR0+Ambarmbar2*Cmbarmbarm1p1 d2R0)Exp[-I(rphase-\[Theta]phase)] 
        + ((Ann0*Cnnm1m1 + Anmbar0*Cnmbarm1m1 + Ambarmbar0*Cmbarmbarm1m1) R0-(Anmbar1*Cnmbarm1m1 + Ambarmbar1*Cmbarmbarm1m1) dR0+Ambarmbar2*Cmbarmbarm1m1 d2R0)Exp[-I(rphase+\[Theta]phase)];
    Clear[r0, R0, dR0, d2R0, \[Theta]0, S0, dS0, d2S0, L1, L2, L2S, L2p, L1Sp, L1L2S, ur0, u\[Theta]0, rcomp, \[Theta]comp, \[CapitalDelta], Kt, \[Rho], \[Rho]bar, \[CapitalSigma], Ann0, Anmbar0, Anmbar1, Ambarmbar0, Ambarmbar1, Ambarmbar2, Cnnp1p1, Cnmbarp1p1, Cmbarmbarp1p1, Cnnp1m1, Cnmbarp1m1, Cmbarmbarp1m1, Cnnm1p1, Cnmbarm1p1, Cmbarmbarm1p1, Cnnm1m1,Cnmbarm1m1,Cmbarmbarm1m1, rphase, \[Theta]phase];
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

  Clear[a, p, \[ScriptCapitalE], \[ScriptCapitalL], rpi, \[Theta]pi, urpi, u\[Theta]pi, \[CapitalUpsilon]t, \[CapitalUpsilon]r, \[CapitalUpsilon]\[Theta], \[CapitalDelta]tr, \[CapitalDelta]t\[Theta], \[CapitalDelta]\[Phi]r, \[CapitalDelta]\[Phi]\[Theta], qt0, qr0, q\[Theta]0, q\[Phi]0, m, \[Omega], \[Lambda], W, rq, \[Theta]q, urq, u\[Theta]q, integrand, integrandIn, integrandUp, \[Alpha]In, \[Alpha]Up, \[Xi], wpIn, wpUp];
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
 Module[{a, r0, RIn, RUp, dRIn, dRUp, W, \[Alpha], S, tp, rp, \[Theta]p, \[Phi]p, l, m, \[Omega], \[CapitalUpsilon]\[Theta], wp, ZIn, ZOut},
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

  wp = With[{\[Lambda]=0}, Precision[TS["\[Alpha]"][\[Theta]p[\[Lambda]]]SH[\[Theta]p[\[Lambda]],0] Cos[\[Omega] tp[\[Lambda]] - m \[Phi]p[\[Lambda]]]]];
  \[Alpha] = (r0^2 - 2 r0 + a^2)/r0 Quiet[NIntegrate[TS["\[Alpha]"][\[Theta]p[\[Lambda]]]SH[\[Theta]p[\[Lambda]],0] Cos[\[Omega] tp[\[Lambda]] - m \[Phi]p[\[Lambda]]], {\[Lambda], 0, \[Pi]/(2\[CapitalUpsilon]\[Theta])},
      Method -> {"Trapezoidal", "SymbolicProcessing"->0}, WorkingPrecision -> wp], NIntegrate::precw];

  ZIn = \[Alpha] RUp/W;
  ZOut = \[Alpha] RIn/W;

  Clear[a, r0, RIn, RUp, dRIn, dRUp, W, \[Alpha], S, tp, rp, \[Theta]p, \[Phi]p, l, m, \[Omega], \[CapitalUpsilon]\[Theta], wp];
  
  <| "\[ScriptCapitalI]" -> ZOut, "\[ScriptCapitalH]" -> ZIn |>
]


(* ::Subsection::Closed:: *)
(*s=0 point particle on an eccentric orbit*)


ConvolveSourcePointParticleEccentric[0, n_Integer, R_, S_, TS_] :=
 Module[{a, p, rpi, rq, \[Theta]q, qt0, qr0, q\[Theta]0, q\[Phi]0, \[Xi], W, \[Alpha]1In, \[Alpha]1Up, S0, l, m, \[Omega], \[CapitalUpsilon]t, \[CapitalUpsilon]r, \[CapitalDelta]tr, \[CapitalDelta]\[Phi]r, II1In, II1Up, wpIn, wpUp, ZIn, ZUp},
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

  wpIn = Precision[II1In[0]];
  wpUp = Precision[II1Up[0]];

  \[Alpha]1In = 1/(2\[Pi]) S0 Quiet[NIntegrate[II1In[qr], {qr,0,\[Pi]}, Method -> {"Trapezoidal", "SymbolicProcessing"->0}, WorkingPrecision -> wpIn], NIntegrate::precw];
  \[Alpha]1Up = 1/(2\[Pi]) S0 Quiet[NIntegrate[II1Up[qr], {qr,0,\[Pi]}, Method -> {"Trapezoidal", "SymbolicProcessing"->0}, WorkingPrecision -> wpUp], NIntegrate::precw];

  \[Xi]=m(\[CapitalDelta]\[Phi]r[qr0]-q\[Phi]0) - \[Omega](\[CapitalDelta]tr[qr0]-qt0) - n qr0;

  ZIn = \[Alpha]1In/W/\[CapitalUpsilon]t Exp[I \[Xi]];
  ZUp = \[Alpha]1Up/W/\[CapitalUpsilon]t Exp[I \[Xi]];

  Clear[a, p, rpi, rq, \[Theta]q, qt0, qr0, q\[Theta]0, q\[Phi]0, \[Xi], W, \[Alpha]1In, \[Alpha]1Up, S0, l, m, \[Omega], \[CapitalUpsilon]t, \[CapitalUpsilon]r, \[CapitalDelta]tr, \[CapitalDelta]\[Phi]r, II1In, II1Up, wpIn, wpUp];
  <| "\[ScriptCapitalI]" -> ZUp, "\[ScriptCapitalH]" -> ZIn |>
]


(* ::Subsection::Closed:: *)
(*s=0 point particle on a generic orbit*)


ConvolveSourcePointParticleGeneric[0, n_Integer, k_Integer, R_, SH_, TS_] :=
Module[{a, p, rpi, \[Theta]pi, rq, \[Theta]q, qt0, qr0, q\[Theta]0, q\[Phi]0, \[Xi], W, \[Alpha]1In, \[Alpha]1Up, \[Alpha]2, \[Alpha]3In, \[Alpha]3Up, \[Alpha]4, S, l, m, \[Omega], \[CapitalUpsilon]t, \[CapitalUpsilon]r, \[CapitalUpsilon]\[Theta], \[CapitalDelta]tr, \[CapitalDelta]t\[Theta], \[CapitalDelta]\[Phi]r, \[CapitalDelta]\[Phi]\[Theta], II1In, II1Up, II2, II3In, II3Up, II4, wp12, wp34, ZIn, ZUp},
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
  wp12 = Precision[{II1In[0], II1Up[0], II2[0]}];
  \[Alpha]1In = 1/(2\[Pi]) Quiet[NIntegrate[II1In[qr], {qr,0,\[Pi]}, Method -> {"Trapezoidal", "SymbolicProcessing"->0}, WorkingPrecision -> wp12], NIntegrate::precw];
  \[Alpha]1Up = 1/(2\[Pi]) Quiet[NIntegrate[II1Up[qr], {qr,0,\[Pi]}, Method -> {"Trapezoidal", "SymbolicProcessing"->0}, WorkingPrecision -> wp12], NIntegrate::precw];
  \[Alpha]2 = 1/(2\[Pi]) Quiet[NIntegrate[II2[q\[Theta]], {q\[Theta],0,\[Pi]}, Method -> {"Trapezoidal", "SymbolicProcessing"->0}, WorkingPrecision -> wp12], NIntegrate::precw];

  If[a == 0,
    \[Alpha]3In = \[Alpha]3Up = \[Alpha]4 = 0;,
    II3In = Function[qr, Evaluate[2TS["\[Alpha]3"][rq[qr]]R["In"][rq[qr]]Cos[(\[Omega] \[CapitalDelta]tr[qr] - m \[CapitalDelta]\[Phi]r[qr] + n qr)]]];
    II3Up = Function[qr, Evaluate[2TS["\[Alpha]3"][rq[qr]]R["Up"][rq[qr]]Cos[(\[Omega] \[CapitalDelta]tr[qr] - m \[CapitalDelta]\[Phi]r[qr] + n qr)]]];
    II4 = Function[q\[Theta], Evaluate[2TS["\[Alpha]4"][\[Theta]q[q\[Theta]]]SH[\[Theta]q[q\[Theta]], 0]Cos[(\[Omega] \[CapitalDelta]t\[Theta][q\[Theta]] - m \[CapitalDelta]\[Phi]\[Theta][q\[Theta]] + k q\[Theta])]]];
    wp34 = Precision[{II3In[0], II3Up[0], II4[0]}];
    \[Alpha]3In = 1/(2\[Pi]) Quiet[NIntegrate[II3In[qr], {qr,0,\[Pi]}, Method -> {"Trapezoidal", "SymbolicProcessing"->0}, WorkingPrecision -> wp34], NIntegrate::precw];
    \[Alpha]3Up = 1/(2\[Pi]) Quiet[NIntegrate[II3Up[qr], {qr,0,\[Pi]}, Method -> {"Trapezoidal", "SymbolicProcessing"->0}, WorkingPrecision -> wp34], NIntegrate::precw];
    \[Alpha]4 = 1/(2\[Pi]) Quiet[NIntegrate[II4[q\[Theta]], {q\[Theta],0,\[Pi]}, Method -> {"Trapezoidal", "SymbolicProcessing"->0}, WorkingPrecision -> wp34], NIntegrate::precw];
  ];

  \[Xi]=m(\[CapitalDelta]\[Phi]r[qr0]+\[CapitalDelta]\[Phi]\[Theta][q\[Theta]0]-q\[Phi]0) - \[Omega](\[CapitalDelta]tr[qr0]+\[CapitalDelta]t\[Theta][q\[Theta]0]-qt0) - k q\[Theta]0 - n qr0;

  ZIn = (\[Alpha]1Up*\[Alpha]2 + \[Alpha]3Up*\[Alpha]4)/W/\[CapitalUpsilon]t Exp[I \[Xi]];
  ZUp = (\[Alpha]1In*\[Alpha]2 + \[Alpha]3In*\[Alpha]4)/W/\[CapitalUpsilon]t Exp[I \[Xi]];

  Clear[a, p, rpi, \[Theta]pi, rq, \[Theta]q, qt0, qr0, q\[Theta]0, q\[Phi]0, \[Xi], W, \[Alpha]1In, \[Alpha]1Up, \[Alpha]2, \[Alpha]3In, \[Alpha]3Up, \[Alpha]4, S, l, m, \[Omega], \[CapitalUpsilon]t, \[CapitalUpsilon]r, \[CapitalUpsilon]\[Theta], \[CapitalDelta]tr, \[CapitalDelta]t\[Theta], \[CapitalDelta]\[Phi]r, \[CapitalDelta]\[Phi]\[Theta], II1In, II1Up, II2, II3In, II3Up, II4, wp12, wp34];
  <| "\[ScriptCapitalI]" -> ZUp, "\[ScriptCapitalH]" -> ZIn |>
]


(* ::Section::Closed:: *)
(*End Package*)


End[];
EndPackage[];
