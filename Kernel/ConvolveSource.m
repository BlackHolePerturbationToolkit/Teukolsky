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


ConvolveSourcePointParticleCircular[-2, R_, SH_, TS_] :=
 Module[{a, r0, m, \[Omega], \[CapitalDelta], W, Ann0, Anmb0, Ambmb0, Anmb1, Ambmb1, Ambmb2, RIn, ROut, dRIn, dROut, CIn, COut, ZIn, ZOut,S, dS, d2S, L2dagS, L1dagL2dagS, \[Rho], \[Rho]b, K},
  a  = TS["Orbit"]["a"];
  r0 = TS["Orbit"]["p"];
  m  = R["In"]["m"];
  \[Omega]  = R["In"]["\[Omega]"];

  RIn = R["In"][r0];
  ROut = R["Up"][r0];

  dRIn = R["In"]'[r0];
  dROut = R["Up"]'[r0];

  \[CapitalDelta] = r0^2 - 2r0 + a^2;
  W = (RIn dROut - ROut dRIn)/\[CapitalDelta];
  K = (r0^2 + a^2)\[Omega] - m a;

  S = SH[\[Pi]/2, 0];
  dS = D[SH[\[Theta],0],\[Theta]]/.\[Theta]->\[Pi]/2;
  d2S = D[SH[\[Theta],0],{\[Theta],2}]/.\[Theta]->\[Pi]/2;

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

  (*FIXME, this is slow to compute the second derivative given we've already computed the R and dR*)
  CIn = RIn(Ann0 + Anmb0 + Ambmb0) - dRIn(Anmb1 + Ambmb1) + R["In"]''[r0] Ambmb2;
  COut = ROut(Ann0 + Anmb0 + Ambmb0) - dROut(Anmb1 + Ambmb1) + R["Up"]''[r0] Ambmb2;

  ZIn = 2 \[Pi] COut/W;
  ZOut = 2 \[Pi] CIn/W;

  Clear[a, r0, m, \[Omega], \[CapitalDelta], W, Ann0, Anmb0, Ambmb0, Anmb1, Ambmb1, Ambmb2, RIn, ROut, dRIn, dROut, CIn, COut, S, dS, d2S, L2dagS, L1dagL2dagS, \[Rho], \[Rho]b, K];

  <| "\[ScriptCapitalI]" -> ZOut, "\[ScriptCapitalH]" -> ZIn |>
]


(* ::Subsection::Closed:: *)
(*s=-2 point particle on a spherical orbit*)


ConvolveSourcePointParticleSpherical[s:-2, k_Integer, R_, SH_, TS_] :=
 Module[{l, m, orbit, a, r0, e, x, \[ScriptCapitalE], \[ScriptCapitalL], \[ScriptCapitalQ], \[CapitalOmega]\[Theta], \[CapitalOmega]\[Phi], T\[Theta],\[CapitalUpsilon]\[Theta], \[Omega], \[CapitalDelta], Kf, dK\[CapitalDelta], Const, \[Lambda], Rupr0, Rinr0, dRupr0, dRinr0, d2Rupr0, d2Rinr0, integrand, Zin, Zup},
  l = R["In"]["l"];
  m = R["In"]["m"];

  orbit = TS["Orbit"];
  a  = orbit["a"];
  r0 = orbit["p"];
  e  = orbit["e"];
  x  = orbit["Inclination"];

  \[ScriptCapitalE] = orbit["Energy"];
  \[ScriptCapitalL] = orbit["AngularMomentum"];
  \[ScriptCapitalQ] = orbit["CarterConstant"];

  {\[CapitalOmega]\[Theta], \[CapitalOmega]\[Phi]} = {"\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(\[Theta]\)]\)","\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(\[Phi]\)]\)"} /. KerrGeoFrequencies[a,r0,e,x];
  T\[Theta] = (2\[Pi])/\[CapitalOmega]\[Theta];
  \[Omega] = m \[CapitalOmega]\[Phi] + k \[CapitalOmega]\[Theta];

  \[CapitalDelta] = r0^2-2r0+a^2;
  Kf=(r0^2+a^2)\[Omega]-m a;
  dK\[CapitalDelta]=(2 (a m (-1+r0)+a^2 \[Omega]-r0^2 \[Omega]))/\[CapitalDelta]^2;
  
  Const = 2 I \[Omega] R["In"]["Amplitudes"]["Incidence"](* (Rin[r0]Rup'[r0] - Rin'[r0]Rup'[r0])/\[CapitalDelta]*);
  \[Lambda] = R["Up"]["Eigenvalue"];

  Rupr0 = R["Up"][r0];
  Rinr0 = R["In"][r0];
  dRupr0 = R["Up"]'[r0];
  dRinr0 = R["In"]'[r0];
  d2Rupr0 = (-(-\[Lambda] + 2 I r0 s 2 \[Omega] + (-2 I (-1 + r0) s (-a m + (a^2 + r0^2) \[Omega]) + (-a m + (a^2 + r0^2) \[Omega])^2)/(a^2 - 2 r0 + r0^2)) Rupr0 - (-2 + 2 r0) (1 + s) dRupr0)/(a^2 - 2 r0 + r0^2);
  d2Rinr0 = (-(-\[Lambda] + 2 I r0 s 2 \[Omega] + (-2 I (-1 + r0) s (-a m + (a^2 + r0^2) \[Omega]) + (-a m + (a^2 + r0^2) \[Omega])^2)/(a^2 - 2 r0 + r0^2)) Rinr0 - (-2 + 2 r0) (1 + s) dRinr0)/(a^2 - 2 r0 + r0^2);

  integrand[\[Lambda]0_?NumericQ, comp_]:=
   Module[{SH\[Theta], dSH\[Theta], d2SH\[Theta], \[Theta], z, \[CapitalSigma], \[Rho], \[Rho]b, Ld2S, Ld1Ld2S, \[CapitalTheta], T, t0, \[Phi]0, CnnPlus, CnnMinus, CnmPlus, CnmMinus, CmmPlus, CmmMinus, Ann0Plus, Anm0Plus, Amm0Plus, Anm1Plus, Amm1Plus, Amm2Plus, Ann0Minus, Anm0Minus, Amm0Minus, Anm1Minus, Amm1Minus, Amm2Minus, IUpPlus, IUpMinus, IInPlus, IInMinus, res},
    \[CapitalSigma] = r0^2+a^2 Cos[\[Theta]]^2;
    \[Rho] = -1/(r0-I a Cos[\[Theta]]);
    \[Rho]b = -1/(r0+I a Cos[\[Theta]]);
    \[Theta] = orbit["Trajectory"][[3]][\[Lambda]0];
    t0 = orbit["Trajectory"][[1]][\[Lambda]0];
    \[Phi]0 = orbit["Trajectory"][[4]][\[Lambda]0];
    
    T = Derivative[1][orbit["Trajectory"][[1]]][\[Lambda]0];
    
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

    IUpPlus=Rupr0(Ann0Plus + Anm0Plus + Amm0Plus)-dRupr0(Anm1Plus+Amm1Plus)+d2Rupr0 Amm2Plus;
    
    IInPlus=Rinr0(Ann0Plus + Anm0Plus + Amm0Plus)-dRinr0(Anm1Plus+Amm1Plus)+d2Rinr0 Amm2Plus;
    
    res = {T(E^(I(\[Omega] t0-m \[Phi]0)) IInPlus), T(E^(I(\[Omega] t0-m \[Phi]0)) IUpPlus)}[[comp]];

    Clear[SH\[Theta], dSH\[Theta], d2SH\[Theta], \[Theta], z, \[CapitalSigma], \[Rho], \[Rho]b, Ld2S, Ld1Ld2S, \[CapitalTheta], T, t0, \[Phi]0, CnnPlus, CnnMinus, CnmPlus, CnmMinus, CmmPlus, CmmMinus, Ann0Plus, Anm0Plus, Amm0Plus, Anm1Plus, Amm1Plus, Amm2Plus, Ann0Minus, Anm0Minus, Amm0Minus, Anm1Minus, Amm1Minus, Amm2Minus, IUpPlus, IUpMinus, IInPlus, IInMinus];
    res    
  ];

  \[CapitalUpsilon]\[Theta] = orbit["Frequencies"]["\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Theta]\)]\)"];
  Zin = ((2\[Pi])/(T\[Theta] Const)) Quiet[NIntegrate[integrand[\[Lambda]0,1], {\[Lambda]0, 0, (2\[Pi])/\[CapitalUpsilon]\[Theta]}, WorkingPrecision -> Precision[R["In"]["\[Omega]"]], Method -> "Trapezoidal"], NIntegrate::precw];
  Zup = ((2\[Pi])/(T\[Theta] Const)) Quiet[NIntegrate[integrand[\[Lambda]0,2], {\[Lambda]0, 0, (2\[Pi])/\[CapitalUpsilon]\[Theta]}, WorkingPrecision -> Precision[R["In"]["\[Omega]"]], Method -> "Trapezoidal"], NIntegrate::precw];
  
  Clear[l, m, orbit, a, r0, e, x, \[ScriptCapitalE], \[ScriptCapitalL], \[ScriptCapitalQ], \[CapitalOmega]\[Theta], \[CapitalOmega]\[Phi], T\[Theta], \[CapitalUpsilon]\[Theta], \[Omega], \[CapitalDelta], Kf, dK\[CapitalDelta], Const, \[Lambda], Rupr0, Rinr0, dRupr0, dRinr0, d2Rupr0, d2Rinr0];
  Remove[integrand];

  <|"\[ScriptCapitalI]" -> Zin, "\[ScriptCapitalH]" -> Zup|>
]


(* ::Subsection::Closed:: *)
(*s=-2 point particle on an eccentric orbit*)


ConvolveSourcePointParticleEccentric[s_:-2, n1_Integer, R_, SH_, TS_] :=
 Module[{a, p, x, W, \[Alpha]In, \[Alpha]Up, ZIn, ZOut, S, rq, \[Theta]q, \[Xi],l, m, \[Omega]mkn,qt0,qr0,q\[Theta]0,q\[Phi]0,\[CapitalUpsilon]t,\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],\[CapitalDelta]tr,\[CapitalDelta]t\[Theta],\[CapitalDelta]\[Phi]r,\[CapitalDelta]\[Phi]\[Theta],
 integrandIn,integrandUp,pg,RinCache,RupCache,RinPCache,RupPCache,RinPPCache,RupPPCache,\[CapitalDelta]tri,\[CapitalDelta]t\[Theta]i,rpi,zpi,\[CapitalDelta]\[Phi]ri,\[CapitalDelta]\[Phi]\[Theta]i,integrand},
  a  = TS["Orbit"]["a"];
  p  = TS["Orbit"]["p"];
  {\[CapitalUpsilon]t,\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta]} = {"\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(t\)]\)","\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(r\)]\)","\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Theta]\)]\)"}/.TS["Orbit"]["Frequencies"];
  rpi=TS["Orbit"]["Trajectory"][[2]];
  {\[CapitalDelta]tri, \[CapitalDelta]t\[Theta]i, \[CapitalDelta]\[Phi]ri, \[CapitalDelta]\[Phi]\[Theta]i} = {"\[CapitalDelta]tr", "\[CapitalDelta]t\[Theta]", "\[CapitalDelta]\[Phi]r", "\[CapitalDelta]\[Phi]\[Theta]"} /. TS["Orbit"]["TrajectoryDeltas"];
  {qt0,qr0,q\[Theta]0,q\[Phi]0}=TS["Orbit"]["InitialPhases"];

  \[CapitalDelta]tr[qr_]:=\[CapitalDelta]tr[qr]=\[CapitalDelta]tri[qr];
  rq[qr_]:=rq[qr]=rpi[(qr-qr0)/\[CapitalUpsilon]r];
  \[Theta]q=Pi/2;
  \[CapitalDelta]\[Phi]r[qr_]:=\[CapitalDelta]\[Phi]r[qr]=\[CapitalDelta]\[Phi]ri[qr];
  
  \[Omega]mkn = R["In"]["\[Omega]"];
  l = R["In"]["l"];
  m = R["In"]["m"];
  pg = Precision[\[Omega]mkn]/2; 
  
  integrand[qr_, Rp0_, Rp1_, Rp2_]:=
  Module[{rp,thp,Rt,RtP,RtPP,Slm,SlmP,SlmPP,\[CapitalDelta],Kt,\[Rho],\[Rho]bar,L1,L2,L2S,L2p,L1Sp,L1L2S,
  Ann0,Anmbar0,Anmbar1,Ambarmbar0,Ambarmbar1,Ambarmbar2,Cnnp1p1,Cnmbarp1p1,Cmbarmbarp1p1,
  Cnnm1p1,Cnmbarm1p1,Cmbarmbarm1p1,Cnnm1m1,Cnmbarm1m1,Cmbarmbarm1m1,Cnnp1m1,Cnmbarp1m1,Cmbarmbarp1m1,
  prefactorR,prefactorRp,prefactorRpp,rphase,\[Theta]phase},
    rp = rq[qr];
    thp = \[Theta]q;
    Rt = Rp0[rp];
    RtP = Rp1[rp];
    RtPP = Rp2[rp];
  
    Slm = SH[thp,0];
    SlmP = Derivative[1,0][SH][thp,0];
    SlmPP = Derivative[2,0][SH][thp,0];
  
    \[CapitalDelta] = rp*rp + a*a - 2*rp;
    Kt = (rp*rp + a*a)*\[Omega]mkn - m*a;
    \[Rho] = -1/(rp - I*a*Cos[thp]);
    \[Rho]bar = -1/(rp + I*a*Cos[thp]);
    L1 = -m/Sin[thp] + a*\[Omega]mkn*Sin[thp] + Cos[thp]/Sin[thp];
    L2 = -m/Sin[thp] + a*\[Omega]mkn*Sin[thp] + 2*Cos[thp]/Sin[thp];
    L2S = SlmP + L2*Slm;
    L2p = m*Cos[thp]/Sin[thp]^2 + a*\[Omega]mkn*Cos[thp] - 2/Sin[thp]^2;
    L1Sp = SlmPP + L1*SlmP;
    L1L2S = L1Sp + L2p*Slm + L2*SlmP + L1*L2*Slm;
  
    Ann0 = -\[Rho]^(-2)*\[Rho]bar^(-1)*(Sqrt[2]*\[CapitalDelta])^(-2)*(\[Rho]^(-1)*L1L2S + 3*I*a*Sin[thp]*L1*Slm + 3*I*a*Cos[thp]*Slm + 2*I*a*Sin[thp]*SlmP - I*a*Sin[thp]*L2*Slm );
    Anmbar0 = \[Rho]^(-3)*(Sqrt[2]\[CapitalDelta])^(-1)*( (\[Rho] + \[Rho]bar - I*Kt/\[CapitalDelta])*L2S + (\[Rho] - \[Rho]bar)*a*Sin[thp]*Kt/\[CapitalDelta]*Slm );
    Anmbar1 = -\[Rho]^(-3)*(Sqrt[2]\[CapitalDelta])^(-1)*( L2S + I*(\[Rho] - \[Rho]bar)*a*Sin[thp]*Slm );
    Ambarmbar0 = \[Rho]^(-3)*\[Rho]bar*Slm*Kt/\[CapitalDelta]/4*( I*(2*\[Omega]mkn*rp/Kt - 2(rp - 1)/\[CapitalDelta]) + Kt/\[CapitalDelta] + 2*I*\[Rho]);
    Ambarmbar1 = -\[Rho]^(-3)*\[Rho]bar*Slm/2*( I*Kt/\[CapitalDelta] - \[Rho] );
    Ambarmbar2 = -\[Rho]^(-3)*\[Rho]bar*Slm/4;
    
    {Cnnp1p1,Cnmbarp1p1,Cmbarmbarp1p1} = TS["Cab"][rp, thp, 1, 1];
    {Cnnm1p1,Cnmbarm1p1,Cmbarmbarm1p1} = TS["Cab"][rp, thp, -1, 1];
    
    rphase = \[Omega]mkn \[CapitalDelta]tr[qr] - m \[CapitalDelta]\[Phi]r[qr] + n1 qr;
    
    ((Ann0*Cnnp1p1 + Anmbar0*Cnmbarp1p1 + Ambarmbar0*Cmbarmbarp1p1) Rt-(Anmbar1*Cnmbarp1p1 + Ambarmbar1*Cmbarmbarp1p1) RtP+Ambarmbar2*Cmbarmbarp1p1 RtPP)Exp[I rphase] 
    + ((Ann0*Cnnm1p1 + Anmbar0*Cnmbarm1p1 + Ambarmbar0*Cmbarmbarm1p1) Rt-(Anmbar1*Cnmbarm1p1 + Ambarmbar1*Cmbarmbarm1p1) RtP+Ambarmbar2*Cmbarmbarm1p1 RtPP)Exp[-I rphase] 
  ];

  W = (p^2 - 2p + a^2)^(-1)(R["In"][p]R["Up"]'[p]-R["Up"][p]R["In"]'[p]);
  RinCache[r_]:= RinCache[r]=R["In"][r];
  RupCache[r_]:= RupCache[r]=R["Up"][r];
  RinPCache[r_]:= RinPCache[r]=R["In"]'[r];
  RupPCache[r_]:= RupPCache[r]=R["Up"]'[r];
  RinPPCache[r_]:= RinPPCache[r]=R["In"]''[r];
  RupPPCache[r_]:= RupPPCache[r]=R["Up"]''[r];
  
  integrandIn[qr_] := integrand[qr,RinCache,RinPCache,RinPPCache];
  integrandUp[qr_] := integrand[qr,RupCache,RupPCache,RupPPCache];

  \[Alpha]In = SpectralCosineIntegration1D[integrandIn, pg];
  \[Alpha]Up = SpectralCosineIntegration1D[integrandUp, pg];
  
  \[Xi]=m(\[CapitalDelta]\[Phi]r[qr0]-q\[Phi]0) - \[Omega]mkn(\[CapitalDelta]tr[qr0]-qt0) - n1 qr0;

  Clear[RinCache,RupCache,RinPCache,RupPCache,RinPPCache,RupPPCache];
  Clear[\[CapitalDelta]tr,\[CapitalDelta]t\[Theta],rq,\[Theta]q,\[CapitalDelta]\[Phi]r,\[CapitalDelta]\[Phi]\[Theta]];
  
  ZIn = 8Pi \[Alpha]Up/W/\[CapitalUpsilon]t Exp[I \[Xi]];
  ZOut = 8Pi \[Alpha]In/W/\[CapitalUpsilon]t Exp[I \[Xi]];
  <| "\[ScriptCapitalI]" -> ZOut, "\[ScriptCapitalH]" -> ZIn |>
]


(* ::Subsection::Closed:: *)
(*s=-2 point particle on a generic orbit*)


ConvolveSourcePointParticleGeneric[s_:-2, n1_Integer, k1_Integer, R_, SH_, TS_] :=
 Module[{a, p, x, W, \[Alpha]In, \[Alpha]Up, ZIn, ZOut, S, rq, \[Theta]q,\[Xi], l, m, k=k1,qt0,qr0,q\[Theta]0,q\[Phi]0, \[Omega]mkn,\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],\[CapitalUpsilon]t,\[CapitalDelta]tr,\[CapitalDelta]t\[Theta],\[CapitalDelta]\[Phi]r,\[CapitalDelta]\[Phi]\[Theta],
 integrandIn,integrandUp,pg,RinCache,RupCache,RinPCache,RupPCache,RinPPCache,RupPPCache,\[CapitalDelta]tri,\[CapitalDelta]t\[Theta]i,rpi,\[Theta]pi,\[CapitalDelta]\[Phi]ri,\[CapitalDelta]\[Phi]\[Theta]i,integrand},
  a  = TS["Orbit"]["a"];
  p = TS["Orbit"]["p"];
  x = TS["Orbit"]["Inclination"];
  {\[CapitalUpsilon]t,\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta]} = {"\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(t\)]\)","\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(r\)]\)","\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Theta]\)]\)"}/.TS["Orbit"]["Frequencies"];
  {rpi,\[Theta]pi}=TS["Orbit"]["Trajectory"][[2;;3]];
  {\[CapitalDelta]tri,\[CapitalDelta]t\[Theta]i,\[CapitalDelta]\[Phi]ri,\[CapitalDelta]\[Phi]\[Theta]i}=Values[TS["Orbit"]["TrajectoryDeltas"]];
  {qt0,qr0,q\[Theta]0,q\[Phi]0}=TS["Orbit"]["InitialPhases"];

  \[CapitalDelta]tr[qr_]:=\[CapitalDelta]tr[qr]=\[CapitalDelta]tri[qr];
  \[CapitalDelta]t\[Theta][q\[Theta]_]:=\[CapitalDelta]t\[Theta][q\[Theta]]=\[CapitalDelta]t\[Theta]i[q\[Theta]];
  rq[qr_]:=rq[qr]=rpi[(qr-qr0)/\[CapitalUpsilon]r];
  \[Theta]q[q\[Theta]_]:=\[Theta]q[q\[Theta]]=\[Theta]pi[(q\[Theta]-q\[Theta]0)/\[CapitalUpsilon]\[Theta]];
  \[CapitalDelta]\[Phi]r[qr_]:=\[CapitalDelta]\[Phi]r[qr]=\[CapitalDelta]\[Phi]ri[qr];
  \[CapitalDelta]\[Phi]\[Theta][q\[Theta]_]:=\[CapitalDelta]\[Phi]\[Theta][q\[Theta]]=\[CapitalDelta]\[Phi]\[Theta]i[q\[Theta]];
  
  \[Omega]mkn = R["In"]["\[Omega]"];
  l = R["In"]["l"];
  m = R["In"]["m"];
  pg = Precision[\[Omega]mkn]/2; 
  
  integrand[qr_, q\[Theta]_, Rp0_, Rp1_, Rp2_]:=
  Module[{rp,thp,Rt,RtP,RtPP,Slm,SlmP,SlmPP,\[CapitalDelta],Kt,\[Rho],\[Rho]bar,L1,L2,L2S,L2p,L1Sp,L1L2S,
  Ann0,Anmbar0,Anmbar1,Ambarmbar0,Ambarmbar1,Ambarmbar2,Cnnp1p1,Cnmbarp1p1,Cmbarmbarp1p1,
  Cnnm1p1,Cnmbarm1p1,Cmbarmbarm1p1,Cnnm1m1,Cnmbarm1m1,Cmbarmbarm1m1,Cnnp1m1,Cnmbarp1m1,Cmbarmbarp1m1,
  prefactorR,prefactorRp,prefactorRpp,rphase,\[Theta]phase},
    rp = rq[qr];
    thp = \[Theta]q[q\[Theta]];
    Rt = Rp0[rp];
    RtP = Rp1[rp];
    RtPP = Rp2[rp];
  
    Slm = SH[thp,0];
    SlmP = Derivative[1,0][SH][thp,0];
    SlmPP = Derivative[2,0][SH][thp,0];
  
    \[CapitalDelta] = rp*rp + a*a - 2*rp;
    Kt = (rp*rp + a*a)*\[Omega]mkn - m*a;
    \[Rho] = -1/(rp - I*a*Cos[thp]);
    \[Rho]bar = -1/(rp + I*a*Cos[thp]);
    L1 = -m/Sin[thp] + a*\[Omega]mkn*Sin[thp] + Cos[thp]/Sin[thp];
    L2 = -m/Sin[thp] + a*\[Omega]mkn*Sin[thp] + 2*Cos[thp]/Sin[thp];
    L2S = SlmP + L2*Slm;
    L2p = m*Cos[thp]/Sin[thp]^2 + a*\[Omega]mkn*Cos[thp] - 2/Sin[thp]^2;
    L1Sp = SlmPP + L1*SlmP;
    L1L2S = L1Sp + L2p*Slm + L2*SlmP + L1*L2*Slm;
  
    Ann0 = -\[Rho]^(-2)*\[Rho]bar^(-1)*(Sqrt[2]*\[CapitalDelta])^(-2)*(\[Rho]^(-1)*L1L2S + 3*I*a*Sin[thp]*L1*Slm + 3*I*a*Cos[thp]*Slm + 2*I*a*Sin[thp]*SlmP - I*a*Sin[thp]*L2*Slm );
    Anmbar0 = \[Rho]^(-3)*(Sqrt[2]\[CapitalDelta])^(-1)*( (\[Rho] + \[Rho]bar - I*Kt/\[CapitalDelta])*L2S + (\[Rho] - \[Rho]bar)*a*Sin[thp]*Kt/\[CapitalDelta]*Slm );
    Anmbar1 = -\[Rho]^(-3)*(Sqrt[2]\[CapitalDelta])^(-1)*( L2S + I*(\[Rho] - \[Rho]bar)*a*Sin[thp]*Slm );
    Ambarmbar0 = \[Rho]^(-3)*\[Rho]bar*Slm*Kt/\[CapitalDelta]/4*( I*(2*\[Omega]mkn*rp/Kt - 2(rp - 1)/\[CapitalDelta]) + Kt/\[CapitalDelta] + 2*I*\[Rho]);
    Ambarmbar1 = -\[Rho]^(-3)*\[Rho]bar*Slm/2*( I*Kt/\[CapitalDelta] - \[Rho] );
    Ambarmbar2 = -\[Rho]^(-3)*\[Rho]bar*Slm/4;
    
    {Cnnp1p1,Cnmbarp1p1,Cmbarmbarp1p1} = TS["Cab"][rp, thp, 1, 1];
    {Cnnp1m1,Cnmbarp1m1,Cmbarmbarp1m1} = TS["Cab"][rp, thp, 1, -1];
    {Cnnm1p1,Cnmbarm1p1,Cmbarmbarm1p1} = TS["Cab"][rp, thp, -1, 1];
    {Cnnm1m1,Cnmbarm1m1,Cmbarmbarm1m1} = TS["Cab"][rp, thp, -1, -1];
    
    rphase = \[Omega]mkn \[CapitalDelta]tr[qr] - m \[CapitalDelta]\[Phi]r[qr] + n1 qr;
    \[Theta]phase = \[Omega]mkn \[CapitalDelta]t\[Theta][q\[Theta]] - m \[CapitalDelta]\[Phi]\[Theta][q\[Theta]] + k1 q\[Theta];
    
    ((Ann0*Cnnp1p1 + Anmbar0*Cnmbarp1p1 + Ambarmbar0*Cmbarmbarp1p1) Rt-(Anmbar1*Cnmbarp1p1 + Ambarmbar1*Cmbarmbarp1p1) RtP+Ambarmbar2*Cmbarmbarp1p1 RtPP)Exp[I(rphase+\[Theta]phase)] 
    + ((Ann0*Cnnp1m1 + Anmbar0*Cnmbarp1m1 + Ambarmbar0*Cmbarmbarp1m1) Rt-(Anmbar1*Cnmbarp1m1 + Ambarmbar1*Cmbarmbarp1m1) RtP+Ambarmbar2*Cmbarmbarp1m1 RtPP)Exp[I(rphase-\[Theta]phase)] 
    + ((Ann0*Cnnm1p1 + Anmbar0*Cnmbarm1p1 + Ambarmbar0*Cmbarmbarm1p1) Rt-(Anmbar1*Cnmbarm1p1 + Ambarmbar1*Cmbarmbarm1p1) RtP+Ambarmbar2*Cmbarmbarm1p1 RtPP)Exp[-I(rphase-\[Theta]phase)] 
    + ((Ann0*Cnnm1m1 + Anmbar0*Cnmbarm1m1 + Ambarmbar0*Cmbarmbarm1m1) Rt-(Anmbar1*Cnmbarm1m1 + Ambarmbar1*Cmbarmbarm1m1) RtP+Ambarmbar2*Cmbarmbarm1m1 RtPP)Exp[-I(rphase+\[Theta]phase)]
  ];

  W = (p^2 - 2p + a^2)^(-1)(R["In"][p]R["Up"]'[p]-R["Up"][p]R["In"]'[p]);
  RinCache[r_]:= RinCache[r]=R["In"][r];
  RupCache[r_]:= RupCache[r]=R["Up"][r];
  RinPCache[r_]:= RinPCache[r]=R["In"]'[r];
  RupPCache[r_]:= RupPCache[r]=R["Up"]'[r];
  RinPPCache[r_]:= RinPPCache[r]=R["In"]''[r];
  RupPPCache[r_]:= RupPPCache[r]=R["Up"]''[r];
  
  integrandIn[qr_, q\[Theta]_] := integrand[qr,q\[Theta],RinCache,RinPCache,RinPPCache];
  integrandUp[qr_, q\[Theta]_] := integrand[qr,q\[Theta],RupCache,RupPCache,RupPPCache];

  \[Alpha]In = SpectralCosineIntegration2D[integrandIn, pg];
  \[Alpha]Up = SpectralCosineIntegration2D[integrandUp, pg];
  
  \[Xi]=m(\[CapitalDelta]\[Phi]r[qr0]+\[CapitalDelta]\[Phi]\[Theta][q\[Theta]0]-q\[Phi]0) - \[Omega]mkn(\[CapitalDelta]tr[qr0]+\[CapitalDelta]t\[Theta][q\[Theta]0]-qt0) - k1 q\[Theta]0 - n1 qr0;

  Clear[RinCache,RupCache,RinPCache,RupPCache,RinPPCache,RupPPCache];
  Clear[\[CapitalDelta]tr,\[CapitalDelta]t\[Theta],rq,\[Theta]q,\[CapitalDelta]\[Phi]r,\[CapitalDelta]\[Phi]\[Theta]];
  
  ZIn = 8Pi \[Alpha]Up/W/\[CapitalUpsilon]t Exp[I \[Xi]];
  ZOut = 8Pi \[Alpha]In/W/\[CapitalUpsilon]t Exp[I \[Xi]];
  <| "\[ScriptCapitalI]" -> ZOut, "\[ScriptCapitalH]" -> ZIn |>
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


ConvolveSourcePointParticleGeneric[0, n1_Integer, k1_Integer, R_, SH_, TS_] :=
Module[{a, p, x, \[CapitalDelta]tri,\[CapitalDelta]t\[Theta]i,rpi,\[Theta]pi,\[CapitalDelta]\[Phi]ri,\[CapitalDelta]\[Phi]\[Theta]i,rq,\[Theta]q,qt0,qr0,q\[Theta]0,q\[Phi]0,\[Xi], W, \[Alpha]1In, \[Alpha]1Up, \[Alpha]2, \[Alpha]3In, \[Alpha]3Up, \[Alpha]4, ZIn, ZOut, S, tp, rp, \[Theta]p, \[Phi]p, l, m, k=k1, \[Omega]mkn,\[CapitalUpsilon]t,\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],\[CapitalUpsilon]\[Phi],\[CapitalOmega]r,\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi],\[CapitalDelta]tr,\[CapitalDelta]t\[Theta],\[CapitalDelta]\[Phi]r,\[CapitalDelta]\[Phi]\[Theta],II1In,II1Up,II2,II3In,II3Up,II4,pg,RinCache,RupCache},
  a  = TS["Orbit"]["a"];
  p = TS["Orbit"]["p"];
  {\[CapitalUpsilon]t,\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta]} = {"\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(t\)]\)","\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(r\)]\)","\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Theta]\)]\)"}/.TS["Orbit"]["Frequencies"];
  {rpi,\[Theta]pi}=TS["Orbit"]["Trajectory"][[2;;3]];
  {\[CapitalDelta]tri,\[CapitalDelta]t\[Theta]i,\[CapitalDelta]\[Phi]ri,\[CapitalDelta]\[Phi]\[Theta]i}=Values[TS["Orbit"]["TrajectoryDeltas"]];
  {qt0,qr0,q\[Theta]0,q\[Phi]0}=TS["Orbit"]["InitialPhases"];

  \[CapitalDelta]tr[qr_]:=\[CapitalDelta]tr[qr]=\[CapitalDelta]tri[qr];
  \[CapitalDelta]t\[Theta][q\[Theta]_]:=\[CapitalDelta]t\[Theta][q\[Theta]]=\[CapitalDelta]t\[Theta]i[q\[Theta]];
  rq[qr_]:=rq[qr]=rpi[(qr-qr0)/\[CapitalUpsilon]r];
  \[Theta]q[q\[Theta]_]:=\[Theta]q[q\[Theta]]=\[Theta]pi[(q\[Theta]-q\[Theta]0)/\[CapitalUpsilon]\[Theta]];
  \[CapitalDelta]\[Phi]r[qr_]:=\[CapitalDelta]\[Phi]r[qr]=\[CapitalDelta]\[Phi]ri[qr];
  \[CapitalDelta]\[Phi]\[Theta][q\[Theta]_]:=\[CapitalDelta]\[Phi]\[Theta][q\[Theta]]=\[CapitalDelta]\[Phi]\[Theta]i[q\[Theta]];
  
  \[Omega]mkn = R["In"]["\[Omega]"];
  l = R["In"]["l"];
  m = R["In"]["m"];
  pg = Precision[\[Omega]mkn]/2;
  
  If[a!=0,
    If[MatchQ[l+m+k,_?OddQ],Return[<| "\[ScriptCapitalI]" -> 0, "\[ScriptCapitalH]" -> 0 |>]],
    If[MatchQ[l+m+k,_?OddQ]||(Abs[m+k]>l),Return[<| "\[ScriptCapitalI]" -> 0, "\[ScriptCapitalH]" -> 0 |>]]
  ];

  W = (p^2 - 2p + a^2)(R["In"][p]R["Up"]'[p]-R["Up"][p]R["In"]'[p]);
  RinCache[r_]:= RinCache[r]=R["In"][r];
  RupCache[r_]:= RupCache[r]=R["Up"][r];
  
  II1In = Function[qr, 2TS["\[Alpha]1"][rq[qr]]RinCache[rq[qr]]Cos[(\[Omega]mkn \[CapitalDelta]tr[qr] - m \[CapitalDelta]\[Phi]r[qr] + n1 qr)]];
  II1Up = Function[qr, 2TS["\[Alpha]1"][rq[qr]]RupCache[rq[qr]]Cos[(\[Omega]mkn \[CapitalDelta]tr[qr] - m \[CapitalDelta]\[Phi]r[qr] + n1 qr)]];
  II2 = Function[q\[Theta], 2TS["\[Alpha]2"][\[Theta]q[q\[Theta]]]SH[\[Theta]q[q\[Theta]], 0]Cos[(\[Omega]mkn \[CapitalDelta]t\[Theta][q\[Theta]] - m \[CapitalDelta]\[Phi]\[Theta][q\[Theta]] + k1 q\[Theta])]];
  II3In = Function[qr, 2TS["\[Alpha]3"][rq[qr]]RinCache[rq[qr]]Cos[(\[Omega]mkn \[CapitalDelta]tr[qr] - m \[CapitalDelta]\[Phi]r[qr] + n1 qr)]];
  II3Up = Function[qr, 2TS["\[Alpha]3"][rq[qr]]RupCache[rq[qr]]Cos[(\[Omega]mkn \[CapitalDelta]tr[qr] - m \[CapitalDelta]\[Phi]r[qr] + n1 qr)]];
  II4 = Function[q\[Theta], 2TS["\[Alpha]4"][\[Theta]q[q\[Theta]]]SH[\[Theta]q[q\[Theta]], 0]Cos[(\[Omega]mkn \[CapitalDelta]t\[Theta][q\[Theta]] - m \[CapitalDelta]\[Phi]\[Theta][q\[Theta]] + k1 q\[Theta])]];
  
  \[Alpha]1In = SpectralCosineIntegration1D[II1In, pg];
  \[Alpha]1Up = SpectralCosineIntegration1D[II1Up, pg];
  \[Alpha]2 = SpectralCosineIntegration1D[II2, pg];
  \[Alpha]3In = SpectralCosineIntegration1D[II3In, pg];
  \[Alpha]3Up = SpectralCosineIntegration1D[II3Up, pg];
  \[Alpha]4 = SpectralCosineIntegration1D[II4, pg];
  
  \[Xi]=m(\[CapitalDelta]\[Phi]r[qr0]+\[CapitalDelta]\[Phi]\[Theta][q\[Theta]0]-q\[Phi]0) - \[Omega]mkn(\[CapitalDelta]tr[qr0]+\[CapitalDelta]t\[Theta][q\[Theta]0]-qt0) - k1 q\[Theta]0 - n1 qr0;
  
  Clear[RinCache,RupCache];
  Clear[\[CapitalDelta]tr,\[CapitalDelta]t\[Theta],rq,\[Theta]q,\[CapitalDelta]\[Phi]r,\[CapitalDelta]\[Phi]\[Theta]];
  
  ZIn = (\[Alpha]1Up*\[Alpha]2 + \[Alpha]3Up*\[Alpha]4)/W/\[CapitalUpsilon]t Exp[I \[Xi]];
  ZOut = (\[Alpha]1In*\[Alpha]2 + \[Alpha]3In*\[Alpha]4)/W/\[CapitalUpsilon]t Exp[I \[Xi]];
  <| "\[ScriptCapitalI]" -> ZOut, "\[ScriptCapitalH]" -> ZIn |>
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
