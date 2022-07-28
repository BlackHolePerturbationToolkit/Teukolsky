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
 Module[{s, orbit},
  orbit = TS["Orbit"];
  s = TS["s"];
  
  If[
  TS["SourceType"] == "PointParticle",
  Switch[{orbit["e"],Abs[orbit["Inclination"]]},{0,1},Return[ConvolveSourcePointParticleCircular[s,R,S,TS]],{0,x_},Return[ConvolveSourcePointParticleSpherical[s,k,R,S,TS]],{e_,1},Return[ConvolveSourcePointParticleEccentric[s,n,R,S,TS]],{e_,x_},Return[ConvolveSourcePointParticleGeneric[s,n,k,R,S,TS]]]
  ];

  $Failed
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
 Module[{l, m, orbit, a, r0, e, x, \[ScriptCapitalE], \[ScriptCapitalL], \[ScriptCapitalQ], \[CapitalOmega]\[Theta], \[CapitalOmega]\[Phi], T\[Theta], \[Omega], \[CapitalDelta], Kf, dK\[CapitalDelta], \[Beta], \[Gamma], \[Delta], zm, zp, Const, \[Lambda], Rupr0, Rinr0, dRupr0, dRinr0, d2Rupr0, d2Rinr0, integrand, Zin, Zup},
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
  
  \[Beta] = a^2 (1-\[ScriptCapitalE]^2);
  \[Gamma] = \[ScriptCapitalE] ((r0^2+a^2)^2/\[CapitalDelta]-a^2)+a \[ScriptCapitalL] (1-(r0^2+a^2)/\[CapitalDelta]);
  \[Delta] = a \[ScriptCapitalE]((r0^2+a^2)/\[CapitalDelta]-1)-(a^2 \[ScriptCapitalL])/\[CapitalDelta];
  zm = (\[ScriptCapitalL]^2+\[ScriptCapitalQ]+\[Beta]-\[Beta] Sqrt[(\[ScriptCapitalL]^4+(\[ScriptCapitalQ]-\[Beta])^2+2 \[ScriptCapitalL]^2 (\[ScriptCapitalQ]+\[Beta]))/\[Beta]^2])/(2 \[Beta]);
  zp = (\[ScriptCapitalL]^2+\[ScriptCapitalQ]+\[Beta]+\[Beta] Sqrt[(\[ScriptCapitalL]^4+(\[ScriptCapitalQ]-\[Beta])^2+2 \[ScriptCapitalL]^2 (\[ScriptCapitalQ]+\[Beta]))/\[Beta]^2])/(2 \[Beta]);

  Const = 2 I \[Omega] R["In"]["Amplitudes"]["Incidence"](* (Rin[r0]Rup'[r0] - Rin'[r0]Rup'[r0])/\[CapitalDelta]*);
  \[Lambda] = R["Up"]["Eigenvalue"];

  Rupr0 = R["Up"][r0];
  Rinr0 = R["In"][r0];
  dRupr0 = R["Up"]'[r0];
  dRinr0 = R["In"]'[r0];
  d2Rupr0 = (-(-\[Lambda] + 2 I r0 s 2 \[Omega] + (-2 I (-1 + r0) s (-a m + (a^2 + r0^2) \[Omega]) + (-a m + (a^2 + r0^2) \[Omega])^2)/(a^2 - 2 r0 + r0^2)) Rupr0 - (-2 + 2 r0) (1 + s) dRupr0)/(a^2 - 2 r0 + r0^2);
  d2Rinr0 = (-(-\[Lambda] + 2 I r0 s 2 \[Omega] + (-2 I (-1 + r0) s (-a m + (a^2 + r0^2) \[Omega]) + (-a m + (a^2 + r0^2) \[Omega])^2)/(a^2 - 2 r0 + r0^2)) Rinr0 - (-2 + 2 r0) (1 + s) dRinr0)/(a^2 - 2 r0 + r0^2);

  integrand[\[Chi]_?NumericQ, comp_]:=
   Module[{SH\[Theta], dSH\[Theta], d2SH\[Theta], \[Theta], z, \[CapitalSigma], \[Rho], \[Rho]b, Ld2S, Ld1Ld2S, \[CapitalTheta], t0, \[Phi]0, CnnPlus, CnnMinus, CnmPlus, CnmMinus, CmmPlus, CmmMinus, Ann0Plus, Anm0Plus, Amm0Plus, Anm1Plus, Amm1Plus, Amm2Plus, Ann0Minus, Anm0Minus, Amm0Minus, Anm1Minus, Amm1Minus, Amm2Minus, IUpPlus, IUpMinus, IInPlus, IInMinus, res},
    \[CapitalSigma] = r0^2+a^2 Cos[\[Theta]]^2;
    \[Rho] = -1/(r0-I a Cos[\[Theta]]);
    \[Rho]b = -1/(r0+I a Cos[\[Theta]]);
    z = zm Cos[\[Chi]]^2;
    \[Theta] = ArcCos[Sqrt[zm] Cos[\[Chi]]];
    t0 = Sqrt[(zm-2 zp+zm Cos[2 \[Chi]])/(zm-zp)] (a^2 (zm-zp) \[ScriptCapitalE] EllipticE[\[Chi],zm/(zm-zp)]+(a^2 zp \[ScriptCapitalE]+\[Gamma]) EllipticF[\[Chi],zm/(zm-zp)])/(Sqrt[-\[Beta] (zm-2 zp+zm Cos[2 \[Chi]])]);
    \[Phi]0 = (Sqrt[(zm-2 zp+zm Cos[2 \[Chi]])/(zm-zp)] ((-1+zm) \[Delta] EllipticF[\[Chi],zm/(zm-zp)]-\[ScriptCapitalL] EllipticPi[zm/(-1+zm),\[Chi],zm/(zm-zp)]))/((-1+zm) Sqrt[-\[Beta] (zm-2 zp+zm Cos[2 \[Chi]])]);

    SH\[Theta] = SH[\[Theta],0];
    dSH\[Theta] = Derivative[1,0][SH][\[Theta],0];
    d2SH\[Theta] = Derivative[2,0][SH][\[Theta],0];
    Ld2S = dSH\[Theta]+(- m Csc[\[Theta]]+a \[Omega] Sin[\[Theta]]+2 Cot[\[Theta]])SH\[Theta];
    Ld1Ld2S = d2SH\[Theta] +(- m Csc[\[Theta]]+a \[Omega] Sin[\[Theta]]+2 Cot[\[Theta]])dSH\[Theta]+(a \[Omega] Cos[\[Theta]]+ m Cot[\[Theta]]Csc[\[Theta]]-2Csc[\[Theta]]^2)SH\[Theta]+(- m Csc[\[Theta]]+a \[Omega] Sin[\[Theta]]+ Cot[\[Theta]])Ld2S;

    CnnPlus = TS["Cnn+"][\[Chi]];
    CnmPlus = TS["Cnm+"][\[Chi]];
    CmmPlus = TS["Cmm+"][\[Chi]];

    CnnMinus = TS["Cnn-"][\[Chi]];
    CnmMinus = TS["Cnm-"][\[Chi]];
    CmmMinus = TS["Cmm-"][\[Chi]];

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
    IUpMinus=Rupr0(Ann0Minus + Anm0Minus + Amm0Minus)-dRupr0(Anm1Minus+Amm1Minus)+d2Rupr0 Amm2Minus;

    IInPlus=Rinr0(Ann0Plus + Anm0Plus + Amm0Plus)-dRinr0(Anm1Plus+Amm1Plus)+d2Rinr0 Amm2Plus;
    IInMinus=Rinr0(Ann0Minus + Anm0Minus + Amm0Minus)-dRinr0(Anm1Minus+Amm1Minus)+d2Rinr0 Amm2Minus;

    res = {(\[Gamma]+a^2 \[ScriptCapitalE] z)/Sqrt[\[Beta](zp-z)](E^(I(\[Omega] t0-m \[Phi]0)) IInPlus+E^(-I(\[Omega] t0-m \[Phi]0)) IInMinus),
     (\[Gamma]+a^2 \[ScriptCapitalE] z)/Sqrt[\[Beta](zp-z)](E^(I(\[Omega] t0-m \[Phi]0)) IUpPlus+E^(-I(\[Omega] t0-m \[Phi]0)) IUpMinus)}[[comp]];

    Clear[SH\[Theta], dSH\[Theta], d2SH\[Theta], \[Theta], z, \[CapitalSigma], \[Rho], \[Rho]b, Ld2S, Ld1Ld2S, \[CapitalTheta], t0, \[Phi]0, CnnPlus, CnnMinus, CnmPlus, CnmMinus, CmmPlus, CmmMinus, Ann0Plus, Anm0Plus, Amm0Plus, Anm1Plus, Amm1Plus, Amm2Plus, Ann0Minus, Anm0Minus, Amm0Minus, Anm1Minus, Amm1Minus, Amm2Minus, IUpPlus, IUpMinus, IInPlus, IInMinus];
    res
  ];

  Zin = ((2\[Pi])/(T\[Theta] Const)) NIntegrate[integrand[\[Chi]0,1], {\[Chi]0, 0, \[Pi]}, WorkingPrecision -> 0.7 Precision[R["In"]["\[Omega]"]], Method -> "Trapezoidal"];
  Zup = ((2\[Pi])/(T\[Theta] Const)) NIntegrate[integrand[\[Chi]0,2], {\[Chi]0, 0, \[Pi]}, WorkingPrecision -> 0.7 Precision[R["In"]["\[Omega]"]], Method -> "Trapezoidal"];

  Clear[l, m, orbit, a, r0, e, x, \[ScriptCapitalE], \[ScriptCapitalL], \[ScriptCapitalQ], \[CapitalOmega]\[Theta], \[CapitalOmega]\[Phi], T\[Theta], \[Omega], \[CapitalDelta], Kf, dK\[CapitalDelta], \[Beta], \[Gamma], \[Delta], zm, zp, Const, \[Lambda], Rupr0, Rinr0, dRupr0, dRinr0, d2Rupr0, d2Rinr0];
  ClearAll[integrand];

  <|"\[ScriptCapitalI]" -> Zin, "\[ScriptCapitalH]" -> Zup|>
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
 Module[{a, r0, \[Psi]In, \[Psi]Out, d\[Psi]In, d\[Psi]Out, W, \[Alpha], ZIn, ZOut, S},
  a  = TS["Orbit"]["a"];
  r0 = TS["Orbit"]["p"];

  (*\[Psi][r] = r R[r] *)
  \[Psi]In = R["In"][r0] r0;
  \[Psi]Out = R["Up"][r0] r0;

  d\[Psi]In = R["In"][r0] + r0 R["In"]'[r0];
  d\[Psi]Out = R["Up"][r0] + r0 R["Up"]'[r0];

  W = \[Psi]In d\[Psi]Out - \[Psi]Out d\[Psi]In;
 

  S = SH[\[Pi]/2,0];
    
  \[Alpha] = TS["\[Alpha]"] S;

  ZIn = \[Alpha] \[Psi]Out/W;
  ZOut = \[Alpha] \[Psi]In/W;

  Clear[a, r0, \[Psi]In, \[Psi]Out, d\[Psi]In, d\[Psi]Out, W, \[Alpha], ZIn, ZOut, S];

  <| "\[ScriptCapitalI]" -> ZOut, "\[ScriptCapitalH]" -> ZIn |>
]


(* ::Subsection::Closed:: *)
(*s=0 point particle on a spherical orbit*)


ConvolveSourcePointParticleSpherical[0, k1_Integer, R_, SH_, TS_] :=
 Module[{a, r0, x, \[Psi]In, \[Psi]Out, d\[Psi]In, d\[Psi]Out, W, \[Alpha], ZIn, ZOut, S, tp, rp, \[Theta]p, \[Phi]p, l, m, k=k1, \[Omega]mk,\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],\[CapitalUpsilon]\[Phi],\[CapitalOmega]r,\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi]},
  a  = TS["Orbit"]["a"];
  r0 = TS["Orbit"]["p"];
  x = TS["Orbit"]["Inclination"];
  {\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],\[CapitalUpsilon]\[Phi]} = {"\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(r\)]\)", "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Theta]\)]\)", "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Phi]\)]\)"} /. TS["Orbit"]["Frequencies"];
  {\[CapitalOmega]r,\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi]} = {"\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(r\)]\)", "\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(\[Theta]\)]\)", "\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(\[Phi]\)]\)"} /. KerrGeoFrequencies[TS["Orbit"]["a"],TS["Orbit"]["p"],0,TS["Orbit"]["Inclination"]];
  {tp,rp,\[Theta]p,\[Phi]p} = TS["Orbit"]["Trajectory"];
  
  \[Omega]mk = R["In"]["\[Omega]"];
  l = R["In"]["l"];
  m = R["In"]["m"];
  
  (*\[Psi][r] = r R[r] *)
  \[Psi]In = R["In"][r0] r0;
  \[Psi]Out = R["Up"][r0] r0;

  d\[Psi]In = R["In"][r0] + r0 R["In"]'[r0];
  d\[Psi]Out = R["Up"][r0] + r0 R["Up"]'[r0];

  W = \[Psi]In d\[Psi]Out - \[Psi]Out d\[Psi]In;
  
If[a!=0,
  \[Alpha] = If[MatchQ[l+m+k,_?OddQ],0,NIntegrate[TS["\[Alpha]"][\[Theta]p[\[Lambda]]]SH[\[Theta]p[\[Lambda]],0] Cos[\[Omega]mk tp[\[Lambda]] - m \[Phi]p[\[Lambda]]],{\[Lambda],0,\[Pi]/(2\[CapitalUpsilon]\[Theta])},Method->"Trapezoidal",MaxRecursion->20,WorkingPrecision->Precision[\[CapitalOmega]\[Phi]]]],
  \[Alpha] = If[MatchQ[l+m+k,_?OddQ]||(Abs[m+k]>l),0,NIntegrate[((-8 \[CapitalOmega]\[Theta] r0^2)/(r0-2))SphericalHarmonicY[l,m,\[Theta]p[\[Lambda]],0] Cos[\[Omega]mk tp[\[Lambda]] - m \[Phi]p[\[Lambda]]],{\[Lambda],0,\[Pi]/(2\[CapitalUpsilon]\[Theta])},Method->"Trapezoidal",MaxRecursion->20,WorkingPrecision->Precision[\[CapitalOmega]\[Phi]]]]];
  
  ZIn = \[Alpha] \[Psi]Out/W;
  ZOut = \[Alpha] \[Psi]In/W;

  Clear[a, r0, x, \[Psi]In, \[Psi]Out, d\[Psi]In, d\[Psi]Out, W, \[Alpha], S, tp, rp, \[Theta]p, \[Phi]p, l, m, k, \[Omega]mk,\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],\[CapitalUpsilon]\[Phi],\[CapitalOmega]r,\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi]];
  
  <| "\[ScriptCapitalI]" -> ZOut, "\[ScriptCapitalH]" -> ZIn |>
]


(* ::Section::Closed:: *)
(*End Package*)


End[];
EndPackage[];
