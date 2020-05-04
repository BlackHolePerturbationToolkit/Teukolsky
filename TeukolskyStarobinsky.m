(* ::Package:: *)

(* ::Title:: *)
(*TeukolskyStarobinsky*)


(* ::Section::Closed:: *)
(*Create Package*)


(* ::Subsection::Closed:: *)
(*BeginPackage*)


BeginPackage["Teukolsky`TeukolskyStarobinsky`"];


(* ::Subsection::Closed:: *)
(*Usage messages*)


(*TeukolskyStarobinskyUp::usage = "TeukolskyStarobinskyUp[s, l, m, a, \[Lambda], \[Omega], r, {RUp, dRUp}]";
TeukolskyStarobinskyIn::usage = "TeukolskyStarobinskyIn[s, l, m, a, \[Lambda], \[Omega], r, {RIn, dRIn}]";*)


(* ::Subsection::Closed:: *)
(*Begin Private section*)


Begin["`Private`"];


(* ::Section::Closed:: *)
(*TeukolskyStarobinsky*)


R2Func[s_,l_,m_,a_,\[Lambda]_,\[Omega]_,r1_,{R_,dR_}]:=Module[{M=1,\[CapitalDelta],K,d2R,A0,B0,r,R1,R2},

\[CapitalDelta]=r^2-2M r+a^2;
K = (r^2+a^2)\[Omega]-m a;

d2R= (-(-2+2 r) (1+s) dR-(-\[Lambda]+4 I r s \[Omega]+(-2 I (-1+r) s (-a m+(a^2+r^2) \[Omega])+(-a m+(a^2+r^2) \[Omega])^2)/(a^2-2 r+r^2)) R)/(a^2-2 r+r^2);

A0=(8I K(K^2+(r-M)^2)+(-4 I K(\[Lambda]+2)-8I \[Omega] r(r-M))\[CapitalDelta]+8 I \[Omega] \[CapitalDelta]^2) 1/\[CapitalDelta]^3;
B0=(((\[Lambda]+2-2I \[Omega] r)(\[Lambda]+6 I \[Omega] r)+12 I \[Omega] (-I K-(r-M)))\[CapitalDelta]-4 I K(-I K-(r-M))(\[Lambda]+6 I \[Omega] r)) 1/\[CapitalDelta]^3;

R2=(A0 (D[R1[r],r]-(I K)/\[CapitalDelta] R1[r])+B0 R1[r]);

{R2,D[R2,r]}//.{r->r1, R1[r]->R, R1'[r]->dR,R1''[r]->d2R}

];


(* ::Subsection::Closed:: *)
(*Transformation for Up solution*)


(*FIXME: add \[Omega]=0 case, add s\[NotEqual]-2 cases, add inhomogeneous transformation*)
TeukolskyStarobinskyUp[-2, l_, m_, a_, \[Lambda]_, \[Omega]_, r0_, {RUp_,dRUp_}]:=Module[{M=1,\[CapitalDelta],K,A0,B0,R2,R,p,\[Alpha],r,d2RUp},

R2 = R2Func[s,l,m,a,\[Lambda],\[Omega],r0,{RUp,dRUp}];

\[Alpha]=Sqrt[a^2-a m/\[Omega]];
p=\[Lambda]^2 (\[Lambda]+2)^2-8 \[Omega]^2 \[Lambda](\[Alpha]^2 (5\[Lambda]+6)-12a^2)+144\[Omega]^4 \[Alpha]^4+144\[Omega]^2 M^2;

16\[Omega]^4 p^-1 R2

]


(* ::Subsection::Closed:: *)
(*Transformation for In solution*)


TeukolskyStarobinskyIn[-2, l_, m_, a_, \[Lambda]_, \[Omega]_, r0_, {RIn_,dRIn_}]:=Module[{M=1, rp, k0, w1, q1, Q1,R2,A0,B0,r,R,\[CapitalDelta],K,d2RIn},

R2 = R2Func[s,l,m,a,\[Lambda],\[Omega],r0,{RIn,dRIn}];

rp=M+Sqrt[M^2-a^2];
k0=\[Omega]-(m a)/(2M rp);
w1=4k0 M rp;
q1=2(M^2-a^2)^(1/2);
Q1=(w1+2I q1)(w1+I q1)w1(w1-I q1);

Q1^-1 R2

]


(* ::Section::Closed:: *)
(*End Package*)


End[];
EndPackage[];
