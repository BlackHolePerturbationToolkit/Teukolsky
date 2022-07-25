(* ::Package:: *)

(* ::Title:: *)
(*ConvolveSource*)


(* ::Section::Closed:: *)
(*Create Package*)


(* ::Subsection::Closed:: *)
(*BeginPackage*)


BeginPackage["Teukolsky`ConvolveSource`",{"KerrGeodesics`OrbitalFrequencies`"}];


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

  If[TS["SourceType"] == "PointParticle" && orbit["e"] == 0 && Abs[orbit["Inclination"]] == 1,
    Return[ConvolveSourcePointParticleCircular[s,R,S,TS]]
  ];
  If[TS["SourceType"] == "PointParticle" && orbit["e"] == 0 && Abs[orbit["Inclination"]] != 1,
    Return[ConvolveSourcePointParticleSpherical[s,k,R,S,TS]]
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

  <| "\[ScriptCapitalI]" -> ZOut, "\[ScriptCapitalH]" -> ZIn |>
]


(* ::Subsection::Closed:: *)
(*s=-2 point particle on a spherical orbit*)


ConvolveSourcePointParticleSpherical[-2, k_Integer, R_, SH_, TS_] :=
Module[
{l,m,orbit,a,r0,e,x,\[ScriptCapitalE],\[ScriptCapitalL],\[ScriptCapitalQ],\[Beta],\[Gamma],\[Delta],zm,zp,\[Theta],z,\[CapitalDelta],\[CapitalSigma],\[Rho],\[Rho]b,\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi],T\[Theta],\[Omega],Ld2S,Ld1Ld2S,Kf,dK\[CapitalDelta],ut,\[CapitalTheta],\[Chi],t0,t0\[Pi]2,t0\[Chi]m\[Pi]2,\[Phi]0,\[Phi]0\[Pi]2,\[Phi]0\[Chi]m\[Pi]2,CnnPlus,CnnMinus,CnmPlus,CnmMinus,CmmPlus,CmmMinus,Ann0Plus,Anm0Plus,Amm0Plus,Anm1Plus,Amm1Plus,Amm2Plus,Ann0Minus,Anm0Minus,Amm0Minus,Anm1Minus,Amm1Minus,Amm2Minus,Rin,Rup,Const,IUpPlus,IUpMinus,IInPlus,IInMinus,Zin,Zup,assoc, amplitudes,eval,\[Sigma]},

(*Print["Getting l, m, k, orbital params, constants"];*)
l = R["In"]["l"];
m = R["In"]["m"];

orbit = TS["Orbit"];
\[Chi] = TS["\[Chi]"];

a   = orbit["a"];
r0 = orbit["p"];
e   = orbit["e"];
x   = orbit["Inclination"];
(*Print[StringJoin["a = ",ToString[a],", r0 = ",ToString[r0],", e = ",ToString[e],", x = ",ToString[x]]];*)
(*Constants of motion & Frequencies*)
\[ScriptCapitalE] = orbit["Energy"];
\[ScriptCapitalL] = orbit["AngularMomentum"];
\[ScriptCapitalQ] = orbit["CarterConstant"];
(*Print[StringJoin["Energy = ",ToString[\[ScriptCapitalE]],", Angular Momentum = ",ToString[\[ScriptCapitalL]],", Carter constant = ",ToString[\[ScriptCapitalQ]]]];*)


(*Print["Calculating Frequencies"];*)
{\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi]} = {"\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(\[Theta]\)]\)","\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(\[Phi]\)]\)"} /. KerrGeoFrequencies[a,r0,e,x];
T\[Theta] = (2\[Pi])/\[CapitalOmega]\[Theta];
\[Omega] = m \[CapitalOmega]\[Phi] + k \[CapitalOmega]\[Theta];
\[Sigma] =.7 Precision[\[Omega]];
(*Print[StringJoin["\[CapitalOmega]\[Theta] = ",ToString[\[CapitalOmega]\[Theta]],", \[CapitalOmega]\[Phi] = ",ToString[\[CapitalOmega]\[Phi]],", \[Omega] = ",ToString[\[Omega]]]];*)


\[CapitalDelta]=r0^2-2r0+a^2;
\[CapitalSigma]=r0^2+a^2 Cos[\[Theta]]^2;
\[Rho]=-1/(r0-I a Cos[\[Theta]]);
\[Rho]b=-1/(r0+I a Cos[\[Theta]]);
\[Beta]=a^2 (1-\[ScriptCapitalE]^2);
\[Gamma] = \[ScriptCapitalE] ((r0^2+a^2)^2/\[CapitalDelta]-a^2)+a \[ScriptCapitalL] (1-(r0^2+a^2)/\[CapitalDelta]);
\[Delta] = a \[ScriptCapitalE]((r0^2+a^2)/\[CapitalDelta]-1)-(a^2 \[ScriptCapitalL])/\[CapitalDelta];
{zm,zp}={(\[ScriptCapitalL]^2+\[ScriptCapitalQ]+\[Beta]-\[Beta] Sqrt[(\[ScriptCapitalL]^4+(\[ScriptCapitalQ]-\[Beta])^2+2 \[ScriptCapitalL]^2 (\[ScriptCapitalQ]+\[Beta]))/\[Beta]^2])/(2 \[Beta]),(\[ScriptCapitalL]^2+\[ScriptCapitalQ]+\[Beta]+\[Beta] Sqrt[(\[ScriptCapitalL]^4+(\[ScriptCapitalQ]-\[Beta])^2+2 \[ScriptCapitalL]^2 (\[ScriptCapitalQ]+\[Beta]))/\[Beta]^2])/(2 \[Beta])};
z=zm Cos[\[Chi]]^2;
\[Theta]=ArcCos[Sqrt[zm] Cos[\[Chi]]];
(*t0[\[Chi]]/\[Phi]0[\[Chi]]*)
t0=Sqrt[(zm-2 zp+zm Cos[2 \[Chi]])/(zm-zp)] (a^2 (zm-zp) \[ScriptCapitalE] EllipticE[\[Chi],zm/(zm-zp)]+(a^2 zp \[ScriptCapitalE]+\[Gamma]) EllipticF[\[Chi],zm/(zm-zp)])/(Sqrt[-\[Beta] (zm-2 zp+zm Cos[2 \[Chi]])]);
\[Phi]0=(Sqrt[(zm-2 zp+zm Cos[2 \[Chi]])/(zm-zp)] ((-1+zm) \[Delta] EllipticF[\[Chi],zm/(zm-zp)]-\[ScriptCapitalL] EllipticPi[zm/(-1+zm),\[Chi],zm/(zm-zp)]))/((-1+zm) Sqrt[-\[Beta] (zm-2 zp+zm Cos[2 \[Chi]])]);
(*t0[\[Pi]/2]/\[Phi]0[\[Pi]/2]*)
t0\[Pi]2=(Sqrt[-(zp/(zm-zp))] (a^2 (zm-zp) \[ScriptCapitalE] EllipticE[zm/(zm-zp)]+(a^2 zp \[ScriptCapitalE]+\[Gamma]) EllipticK[zm/(zm-zp)]))/Sqrt[zp \[Beta]];
\[Phi]0\[Pi]2=(Sqrt[-(zp/(zm-zp))] ((-1+zm) \[Delta] EllipticK[zm/(zm-zp)]-\[ScriptCapitalL] EllipticPi[zm/(-1+zm),zm/(zm-zp)]))/((-1+zm) Sqrt[zp \[Beta]]);
(*t0[\[Chi]-\[Pi]/2]/\[Phi]0[\[Chi]-\[Pi]/2]*)
t0\[Chi]m\[Pi]2=(Sqrt[(zm-2 zp+zm Cos[2 (-(\[Pi]/2)+\[Chi])])/(zm-zp)] (-a^2 (zm-zp) \[ScriptCapitalE] EllipticE[\[Pi]/2-\[Chi],zm/(zm-zp)]-(a^2 zp \[ScriptCapitalE]+\[Gamma]) EllipticF[\[Pi]/2-\[Chi],zm/(zm-zp)]))/(Sqrt[-\[Beta] (zm-2 zp+zm Cos[2 (-(\[Pi]/2)+\[Chi])])]);
\[Phi]0\[Chi]m\[Pi]2=(Sqrt[(zm-2 zp+zm Cos[2 (-(\[Pi]/2)+\[Chi])])/(zm-zp)] (-((-1+zm) \[Delta] EllipticF[\[Pi]/2-\[Chi],zm/(zm-zp)])+\[ScriptCapitalL] EllipticPi[zm/(-1+zm),\[Pi]/2-\[Chi],zm/(zm-zp)]))/((-1+zm) Sqrt[-\[Beta] (zm-2 zp+zm Cos[2 (-(\[Pi]/2)+\[Chi])])]);

(*Print["Calculating Homogeneous Radial Solution, and SWSH"];*)
Rup = R["Up"];
Rin = R["In"];
Const =2 I \[Omega] Rin["Amplitudes"]["Incidence"](* (Rin[r0]Rup'[r0] - Rin'[r0]Rup'[r0])/\[CapitalDelta]*);
eval = Rup["Eigenvalue"];


Ld2S=Derivative[1,0][SH][\[Theta],0]+(- m Csc[\[Theta]]+a \[Omega] Sin[\[Theta]]+2 Cot[\[Theta]])SH[\[Theta],0];
Ld1Ld2S=Derivative[2,0][SH][\[Theta],0] +(- m Csc[\[Theta]]+a \[Omega] Sin[\[Theta]]+2 Cot[\[Theta]])Derivative[1,0][SH][\[Theta],0]+(a \[Omega] Cos[\[Theta]]+ m Cot[\[Theta]]Csc[\[Theta]]-2Csc[\[Theta]]^2)SH[\[Theta],0]+(- m Csc[\[Theta]]+a \[Omega] Sin[\[Theta]]+ Cot[\[Theta]])Ld2S;

Kf=(r0^2+a^2)\[Omega]-m a;
dK\[CapitalDelta]=(2 (a m (-1+r0)+a^2 \[Omega]-r0^2 \[Omega]))/\[CapitalDelta]^2;
(*Print["Retrieving Stress-energy information"];*)
CnnPlus=TS["Cnn+"];
CnmPlus=TS["Cnm+"];
CmmPlus=TS["Cmm+"];

CnnMinus=TS["Cnn-"];
CnmMinus=TS["Cnm-"];
CmmMinus=TS["Cmm-"];

(*Print["Calculating Subscript[A, abi] terms"];*)
Ann0Plus= (-2 \[Rho]^-3 \[Rho]b^-1 CnnPlus)/\[CapitalDelta]^2 (Ld1Ld2S+2I a \[Rho] Sin[\[Theta]]Ld2S);

Anm0Plus=-((2Sqrt[2] \[Rho]^-3 CnmPlus)/\[CapitalDelta])(((I Kf)/\[CapitalDelta]-\[Rho]-\[Rho]b)Ld2S- Kf/\[CapitalDelta] a Sin[\[Theta]]SH[\[Theta],0](\[Rho]-\[Rho]b));

Amm0Plus=SH[\[Theta],0]\[Rho]^-3 \[Rho]b CmmPlus((Kf/\[CapitalDelta])^2+2I \[Rho] Kf/\[CapitalDelta]+I dK\[CapitalDelta]);

Anm1Plus=-((2Sqrt[2] \[Rho]^-3 CnmPlus)/\[CapitalDelta])(Ld2S+I a Sin[\[Theta]](\[Rho]-\[Rho]b)SH[\[Theta],0]);

Amm1Plus=2SH[\[Theta],0]\[Rho]^-3 \[Rho]b CmmPlus(\[Rho]-(I Kf)/\[CapitalDelta]);

Amm2Plus=-SH[\[Theta],0] \[Rho]^-3 \[Rho]b CmmPlus;

Ann0Minus= (-2 \[Rho]^-3 \[Rho]b^-1 CnnMinus)/\[CapitalDelta]^2 (Ld1Ld2S+2I a \[Rho] Sin[\[Theta]]Ld2S);

Anm0Minus=-((2Sqrt[2] \[Rho]^-3 CnmMinus)/\[CapitalDelta])(((I Kf)/\[CapitalDelta]-\[Rho]-\[Rho]b)Ld2S- Kf/\[CapitalDelta] a Sin[\[Theta]]SH[\[Theta],0](\[Rho]-\[Rho]b));

Amm0Minus=SH[\[Theta],0]\[Rho]^-3 \[Rho]b CmmMinus((Kf/\[CapitalDelta])^2+2I \[Rho] Kf/\[CapitalDelta]+I dK\[CapitalDelta]);

Anm1Minus=-((2Sqrt[2] \[Rho]^-3 CnmMinus)/\[CapitalDelta])(Ld2S+I a Sin[\[Theta]](\[Rho]-\[Rho]b)SH[\[Theta],0]);

Amm1Minus=2SH[\[Theta],0]\[Rho]^-3 \[Rho]b CmmMinus(\[Rho]-(I Kf)/\[CapitalDelta]);

Amm2Minus=-SH[\[Theta],0] \[Rho]^-3 \[Rho]b CmmMinus;

(*Print["Preparing integrand"];*)
IUpPlus=Rup[r0](Ann0Plus + Anm0Plus + Amm0Plus)-Rup'[r0](Anm1Plus+Amm1Plus)+Rup''[r0]Amm2Plus;
IUpMinus=Rup[r0](Ann0Minus + Anm0Minus + Amm0Minus)-Rup'[r0](Anm1Minus+Amm1Minus)+Rup''[r0]Amm2Minus;

IInPlus=Rin[r0](Ann0Plus + Anm0Plus + Amm0Plus)-Rin'[r0](Anm1Plus+Amm1Plus)+Rin''[r0]Amm2Plus;
IInMinus=Rin[r0](Ann0Minus + Anm0Minus + Amm0Minus)-Rin'[r0](Anm1Minus+Amm1Minus)+Rin''[r0]Amm2Minus;

(*Print["Integrating"];*)
Zin = ((2\[Pi])/(T\[Theta] Const)) NIntegrate[((\[Gamma]+a^2 \[ScriptCapitalE] z)/Sqrt[\[Beta](zp-z)] (E^(I(\[Omega] t0-m \[Phi]0)) IInPlus+E^(-I(\[Omega] t0-m \[Phi]0)) IInMinus))/.{\[Chi]->\[Chi]0},{\[Chi]0,0,\[Pi]},WorkingPrecision->\[Sigma]];

Zup = ((2\[Pi])/(T\[Theta] Const)) NIntegrate[((\[Gamma]+a^2 \[ScriptCapitalE] z)/Sqrt[\[Beta](zp-z)] (E^(I(\[Omega] t0-m \[Phi]0)) IUpPlus+E^(-I(\[Omega] t0-m \[Phi]0)) IUpMinus))/.{\[Chi]->\[Chi]0},{\[Chi]0,0,\[Pi]},WorkingPrecision->\[Sigma]];
amplitudes = <|"\[ScriptCapitalI]"->Zin,"\[ScriptCapitalH]"->Zup|>
(*The in/up look backwards, but I'm following the naming convention from Hughes: 	arXiv:gr-qc/9910091*)
(*assoc = <|"l"->l,"m"->m,"k"->k,"\[Omega]"->\[Omega],"Eigenvalue"->eval, "Orbit"->orbit, "Amplitudes"->amplitudes|>*)

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
	(*Print[\[Alpha]];*)
  <| "\[ScriptCapitalI]" -> ZOut, "\[ScriptCapitalH]" -> ZIn |>
]


(* ::Section::Closed:: *)
(*End Package*)


End[];
EndPackage[];
