(* ::Package:: *)

(* ::Title:: *)
(*ConvolveSource*)


(* ::Section::Closed:: *)
(*Create Package*)


(* ::Subsection::Closed:: *)
(*BeginPackage*)


BeginPackage["Teukolsky`ConvolveSource`"];


(* ::Subsection::Closed:: *)
(*Begin Private section*)


Begin["`Private`"];


(* ::Section::Closed:: *)
(*ConvolveSource*)


(* ::Subsection::Closed:: *)
(*ConvolveSource*)


ConvolveSource[R_, S_, TS_] :=
 Module[{s, orbit},
  orbit = TS["Orbit"];
  s = TS["s"];

  If[TS["SourceType"] == "PointParticle" && orbit["e"] == 0 && Abs[orbit["Inclination"]] == 1,
    Return[ConvolveSourcePointParticleCircular[s,R,S,TS]]
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

  S = SH[\[Pi]/2, 0];
    
  \[Alpha] = TS["\[Alpha]"] S;

  ZIn = \[Alpha] \[Psi]Out/W;
  ZOut = \[Alpha] \[Psi]In/W;

  <| "\[ScriptCapitalI]" -> ZOut, "\[ScriptCapitalH]" -> ZIn |>
]


(* ::Section::Closed:: *)
(*End Package*)


End[];
EndPackage[];
