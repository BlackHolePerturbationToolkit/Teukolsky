(* ::Package:: *)

BeginPackage["Teukolsky`TeukolskyRadialUp`"];


TeukolskyRadialUp::usage = 
"TeukolskyRadialUp[s, \[Lambda], m, a, \[Omega], r]";

SasakiNakamuraRadialUp::usage = 
"SasakiNakamuraRadialUp[\[Lambda], m, a, \[Omega], r]";

SasakiNakamuraEquation::usage =
"SasakiNakamuraEquation[\[Lambda], m, a, \[Omega], r, {X[r], X'[r], X''[r]}]";

TeuksolskyRadialEquation::usage =
"TeukolskyRadialEquation[s, \[Lambda], m, a, \[Omega], r, {R[r], R'[r], R''[r]}]";
 
Begin["`Private`"];


(**********************************************************)
(* Internal functions                                     *)
(**********************************************************)
 
(* Boundary conditions for the up solution of the Sasaki-Nakamura equation*)
(* They were derived in S. E. Gralla, A. P. Porfyriadis, N. Warburton, Phys. Rev. D 92, 064029 (2015), arXiv:1506.08496*)
OuterBCs[\[Lambda]_, m_, a_, \[Omega]_, rout_]:=Module[{M=1,recur,k,j,rp,rm,rs,X,r,c},
  recur=-((6 I a^12 (210-29 k+k^2) \[Omega]^12 c[-13+k])/(k (2 \[Lambda]+\[Lambda]^2-12 \[Omega] (-a m+I M+a^2 \[Omega]))))+(12 a^10 (-14+k) \[Omega]^11 (-a (-12+k) m+3 I (-13+k) M+a^2 (-13+k) \[Omega]) c[-12+k])/(k (2 \[Lambda]+\[Lambda]^2-12 \[Omega] (-a m+I M+a^2 \[Omega])))+(6 I a^8 \[Omega]^10 (-2 I a (729-121 k+5 k^2) m M-12 (156-25 k+k^2) M^2-4 a^3 (143-24 k+k^2) m \[Omega]+2 a^4 (156-25 k+k^2) \[Omega]^2+a^2 (k^2 (-5+2 m^2+10 I M \[Omega])+k (123-46 m^2-250 I M \[Omega])+3 (-251+87 m^2+520 I M \[Omega]))) c[-11+k])/(k (2 \[Lambda]+\[Lambda]^2-12 \[Omega] (-a m+I M+a^2 \[Omega])))+(4 a^6 \[Omega]^9 (-6 a (493-89 k+4 k^2) m M^2+12 I (132-23 k+k^2) M^3-12 a^5 (-12+k) m \[Omega]^2+6 a^6 (-12+k) \[Omega]^3-3 I a^2 M (4 k^2 (-3+m^2+2 I M \[Omega])-2 k (-135+43 m^2+91 I M \[Omega])+3 (-503+153 m^2+344 I M \[Omega]))+a^4 \[Omega] (-k^2 (-15+\[Lambda]+12 I M \[Omega])+k (-339+6 m^2+20 \[Lambda]+276 I M \[Omega])-3 (-635+25 m^2+32 \[Lambda]+528 I M \[Omega]))+a^3 m (k^2 (-12+\[Lambda]+24 I M \[Omega])+k (264-20 \[Lambda]-534 I M \[Omega])+3 (-481+m^2+32 \[Lambda]+986 I M \[Omega]))) c[-10+k])/(k (2 \[Lambda]+\[Lambda]^2-12 \[Omega] (-a m+I M+a^2 \[Omega])))-(I a^5 \[Omega]^8 (96 I (-11+k) (-9+k) m M^3+48 a^6 m \[Omega]^3+48 a M^2 (3-2 m^2+(-9+k)^2 (9-2 m^2-2 I M \[Omega])+(-9+k) (-21+4 m^2+4 I M \[Omega]))+4 a^5 \[Omega]^2 (-2970-27 k^2-30 m^2+46 \[Lambda]+264 I M \[Omega]+k (567-4 \[Lambda]-24 I M \[Omega]))-8 I a^2 m M (9 m^2-11 \[Lambda]+36 I M \[Omega]+4 (-9+k)^2 (-12+\[Lambda]+6 I M \[Omega])-2 (-9+k) (-45+\[Lambda]+30 I M \[Omega]))+4 a^4 m \[Omega] (51 k^2+k (-1011+4 \[Lambda]+48 I M \[Omega])+6 (4 m^2-9 (-92+\[Lambda]+10 I M \[Omega])))+a^3 ((-9+k)^2 (120-96 m^2+2 \[Lambda]+\[Lambda]^2-492 I M \[Omega]+32 I M \[Lambda] \[Omega]-96 M^2 \[Omega]^2)-4 (-9+6 m^4+2 \[Lambda]^2+54 I M \[Omega]+48 M^2 \[Omega]^2+\[Lambda] (4+22 I M \[Omega])-m^2 (33+8 \[Lambda]+78 I M \[Omega]))+(-9+k) (\[Lambda]^2+96 m^2 (1-I M \[Omega])+2 \[Lambda] (1-8 I M \[Omega])+12 (-22+95 I M \[Omega]+24 M^2 \[Omega]^2)))) c[-9+k])/(2 k (2 \[Lambda]+\[Lambda]^2-12 \[Omega] (-a m+I M+a^2 \[Omega])))+(a^4 \[Omega]^7 (48 I (159-36 k+2 k^2) M^3+8 a m M^2 (k (525-34 \[Lambda])+2 k^2 (-15+\[Lambda])+3 (-751+47 \[Lambda]))-4 a^5 m (51 k-2 (249+\[Lambda])) \[Omega]^2+12 a^6 (-88+9 k) \[Omega]^3+4 a^3 m (-12-15 \[Lambda]-\[Lambda]^2+m^2 (12+\[Lambda])+6 I M \[Omega]-18 I M \[Lambda] \[Omega]+2 (-8+k)^2 (-9+2 \[Lambda]+39 I M \[Omega])+(-8+k) (36+I M (-141+4 \[Lambda]) \[Omega]))+a^4 \[Omega] (-8 (-8+k)^2 (-15+2 \[Lambda]+21 I M \[Omega])+4 (9+\[Lambda]^2-3 m^2 (19+\[Lambda])-54 I M \[Omega]+2 \[Lambda] (6+5 I M \[Omega]))+(-8+k) (-264+96 m^2-\[Lambda]^2+456 I M \[Omega]+\[Lambda] (-2-16 I M \[Omega])))-I a^2 M (2 (-8+k)^2 (-108+72 m^2-\[Lambda]^2+156 I M \[Omega]+\[Lambda] (-2-8 I M \[Omega]))-(-8+k) (-432+168 m^2+\[Lambda]^2+540 I M \[Omega]+2 \[Lambda] (1-8 I M \[Omega]))-4 (m^2 (30+8 \[Lambda])-3 (\[Lambda]^2-18 I M \[Omega]+\[Lambda] (2+2 I M \[Omega]))))) c[-8+k])/(k (2 \[Lambda]+\[Lambda]^2-12 \[Omega] (-a m+I M+a^2 \[Omega])))+(I a^3 \[Omega]^6 (-192 I (54-15 k+k^2) m M^3-168 a^6 m \[Omega]^3+4 a^5 \[Omega]^2 (3459+48 k^2+93 m^2-149 \[Lambda]-732 I M \[Omega]+4 k (-204+4 \[Lambda]+21 I M \[Omega]))-2 a^4 m \[Omega] (10086+168 k^2+126 m^2-348 \[Lambda]-\[Lambda]^2-2772 I M \[Omega]+8 k (-327+4 \[Lambda]+39 I M \[Omega]))+8 I a^2 m M (3+18 m^2-37 \[Lambda]-2 \[Lambda]^2-24 I M \[Omega]+6 (-7+k)^2 (-9+2 \[Lambda]+9 I M \[Omega])-2 (-7+k) (-42+\[Lambda]+42 I M \[Omega]))+4 a M^2 (-12 (-7+k) (-15+4 m^2+4 I M \[Omega])+4 (9-12 m^2+2 \[Lambda]+\[Lambda]^2-36 I M \[Omega])+(-7+k)^2 (-108+48 m^2-2 \[Lambda]-\[Lambda]^2+60 I M \[Omega]))+a^3 (12+48 m^4+54 \[Lambda]+29 \[Lambda]^2+\[Lambda]^3-420 I M \[Omega]+236 I M \[Lambda] \[Omega]+16 I M \[Lambda]^2 \[Omega]-m^2 (156+98 \[Lambda]+\[Lambda]^2+756 I M \[Omega])+4 (-7+k)^2 (-30+36 m^2-2 \[Lambda]-\[Lambda]^2+192 I M \[Omega]-24 I M \[Lambda] \[Omega]+60 M^2 \[Omega]^2)+4 (-7+k) (\[Lambda]^2 (-1-I M \[Omega])+\[Lambda] (-2+2 I M \[Omega])+m^2 (-36+72 I M \[Omega])-6 (-9+52 I M \[Omega]+22 M^2 \[Omega]^2)))) c[-7+k])/(2 k (2 \[Lambda]+\[Lambda]^2-12 \[Omega] (-a m+I M+a^2 \[Omega])))+(a^2 \[Omega]^5 (48 I (-7+k) (-6+k) M^3+8 a m M^2 (21+21 (-6+k)+4 (-6+k)^2 (-6+\[Lambda])-12 \[Lambda])+24 a^5 m (105-14 k+\[Lambda]) \[Omega]^2+96 a^6 (-15+2 k) \[Omega]^3+4 a^3 m (-15-20 \[Lambda]-3 \[Lambda]^2+m^2 (15+2 \[Lambda])-9 I M \[Omega]-43 I M \[Lambda] \[Omega]+6 (-6+k)^2 (-2+\[Lambda]+15 I M \[Omega])+3 (-6+k) (8+I M (-39+4 \[Lambda]) \[Omega]))+4 a^4 \[Omega] (-3+11 \[Lambda]+3 \[Lambda]^2-m^2 (81+8 \[Lambda])-57 I M \[Omega]+27 I M \[Lambda] \[Omega]-6 (-6+k)^2 (-5+\[Lambda]+9 I M \[Omega])+(-6+k) (-54+36 m^2-\[Lambda]^2+135 I M \[Omega]+\[Lambda] (-2-12 I M \[Omega])))-I a^2 M (36+66 \[Lambda]+35 \[Lambda]^2+\[Lambda]^3-4 m^2 (27+16 \[Lambda])-636 I M \[Omega]+84 I M \[Lambda] \[Omega]-(-6+k) (120 m^2+10 \[Lambda]+5 \[Lambda]^2+36 (-6+7 I M \[Omega]))+2 (-6+k)^2 (-72+72 m^2-3 \[Lambda]^2+180 I M \[Omega]+\[Lambda] (-6-16 I M \[Omega])))) c[-6+k])/(k (2 \[Lambda]+\[Lambda]^2-12 \[Omega] (-a m+I M+a^2 \[Omega])))-(I a \[Omega]^4 (96 I (-1+(-5+k)^2) m M^3+216 a^6 m \[Omega]^3-8 I a^2 m M (-6+9 m^2-33 \[Lambda]-4 \[Lambda]^2+12 (-5+k)^2 (-2+\[Lambda]+3 I M \[Omega])+2 (-5+k) (15+\[Lambda]-12 I M \[Omega]))+12 a^5 \[Omega]^2 (-590-14 k^2-34 m^2+57 \[Lambda]+228 I M \[Omega]-2 k (-91+4 \[Lambda]+18 I M \[Omega]))+8 a M^2 (-10 \[Lambda]-5 \[Lambda]^2+72 I M \[Omega]+(-5+k) (-18+2 \[Lambda]+\[Lambda]^2-12 I M \[Omega])+(-5+k)^2 (18-12 m^2+2 \[Lambda]+\[Lambda]^2-24 I M \[Omega]))+6 a^4 m \[Omega] (1360+44 k^2+36 m^2-132 \[Lambda]-\[Lambda]^2-772 I M \[Omega]+4 k (-123+4 \[Lambda]+30 I M \[Omega]))+a^3 (-12-24 m^4-64 \[Lambda]-38 \[Lambda]^2-3 \[Lambda]^3+552 I M \[Omega]-132 I M \[Lambda] \[Omega]-32 I M \[Lambda]^2 \[Omega]+2 m^2 (30+50 \[Lambda]+\[Lambda]^2+288 I M \[Omega])-6 (-5+k)^2 (-10+16 m^2-\[Lambda]^2+92 I M \[Omega]+32 M^2 \[Omega]^2+\[Lambda] (-2-16 I M \[Omega]))+2 (-5+k) (48 m^2 (1-3 I M \[Omega])+\[Lambda]^2 (3+6 I M \[Omega])+\[Lambda] (6+20 I M \[Omega])+42 (-1+6 I M \[Omega]+4 M^2 \[Omega]^2)))) c[-5+k])/(2 k (2 \[Lambda]+\[Lambda]^2-12 \[Omega] (-a m+I M+a^2 \[Omega])))+(a \[Omega]^3 (8 m M^2 (k (51-14 \[Lambda])+21 (-5+\[Lambda])+2 k^2 (-3+\[Lambda]))+24 a^4 m (57-11 k+\[Lambda]) \[Omega]^2+24 a^5 (-36+7 k) \[Omega]^3+2 a^3 \[Omega] (k^2 (30-8 \[Lambda]-60 I M \[Omega])+k (-282+48 m^2+58 \[Lambda]-3 \[Lambda]^2+630 I M \[Omega]-24 I M \[Lambda] \[Omega])-2 (-321+50 \[Lambda]-9 \[Lambda]^2+7 m^2 (21+\[Lambda])+816 I M \[Omega]-72 I M \[Lambda] \[Omega]))+4 a^2 m (-6-11 \[Lambda]-3 \[Lambda]^2+m^2 (6+\[Lambda])-32 I M \[Lambda] \[Omega]+(-4+k)^2 (-3+4 \[Lambda]+42 I M \[Omega])+3 (-4+k) (2+I M (-9+4 \[Lambda]) \[Omega]))-I a M (2 ((26-16 m^2) \[Lambda]+15 \[Lambda]^2+\[Lambda]^3-168 I M \[Omega])+2 (-4+k)^2 (-18+24 m^2-3 \[Lambda]^2+84 I M \[Omega]+\[Lambda] (-6-8 I M \[Omega]))-(-4+k) (-36+24 m^2+7 \[Lambda]^2-60 I M \[Omega]+2 \[Lambda] (7+8 I M \[Omega])))) c[-4+k])/(k (2 \[Lambda]+\[Lambda]^2-12 \[Omega] (-a m+I M+a^2 \[Omega])))-(I \[Omega]^2 (120 a^5 m \[Omega]^3+4 (-4+k) k M^2 (2 \[Lambda]+\[Lambda]^2-12 I M \[Omega])-8 I a m M (-7 \[Lambda]-2 \[Lambda]^2+(-3+k) (3+2 \[Lambda])+6 I M \[Omega]+(-3+k)^2 (-3+4 \[Lambda]+6 I M \[Omega]))+4 a^4 \[Omega]^2 (-363-18 k^2-45 m^2+79 \[Lambda]+228 I M \[Omega]-2 k (-81+8 \[Lambda]+30 I M \[Omega]))+2 a^3 m \[Omega] (48 k^2+30 m^2+8 k (-39+4 \[Lambda]+21 I M \[Omega])-3 (-162+60 \[Lambda]+\[Lambda]^2+212 I M \[Omega]))+a^2 (-3 \[Lambda]^3+12 I M \[Omega] (15+11 m^2+4 I M \[Omega])+\[Lambda]^2 (-21+m^2-16 I M \[Omega])+\[Lambda] (-30+34 m^2+28 I M \[Omega])+4 (-3+k) (-3+6 I M \[Omega]+36 M^2 \[Omega]^2+\[Lambda]^2 (1+3 I M \[Omega])+6 m^2 (1-4 I M \[Omega])+2 \[Lambda] (1+5 I M \[Omega]))-4 (-3+k)^2 (-3+6 m^2-\[Lambda]^2+42 I M \[Omega]+12 M^2 \[Omega]^2+\[Lambda] (-2-8 I M \[Omega])))) c[-3+k])/(2 k (2 \[Lambda]+\[Lambda]^2-12 \[Omega] (-a m+I M+a^2 \[Omega])))-(\[Omega] (8 a^3 m (-33+12 k-\[Lambda]) \[Omega]^2-24 a^4 (-8+3 k) \[Omega]^3+M (3+5 k-2 k^2+\[Lambda]) (2 I \[Lambda]+I \[Lambda]^2+12 M \[Omega])+4 a^2 \[Omega] (-18+\[Lambda]-3 \[Lambda]^2+2 m^2 (12+\[Lambda])+75 I M \[Omega]-15 I M \[Lambda] \[Omega]+k^2 (-3+\[Lambda]+6 I M \[Omega])+k (15-6 m^2-2 \[Lambda]+\[Lambda]^2-45 I M \[Omega]+4 I M \[Lambda] \[Omega]))+4 a m (\[Lambda]^2-3 I (5-7 k+2 k^2) M \[Omega]-\[Lambda] (2+k^2-15 I M \[Omega]+k (-4+4 I M \[Omega])))) c[-2+k])/(k (2 \[Lambda]+\[Lambda]^2-12 \[Omega] (-a m+I M+a^2 \[Omega])))-(I (-\[Lambda]^3+\[Lambda]^2 (k^2-2 (2+a m \[Omega]+2 I M \[Omega])+k (-1+4 I M \[Omega]))+12 \[Omega] (2 a^3 m \[Omega]^2+a^2 \[Omega] (-6+5 k-k^2-2 m^2+4 I M \[Omega]-4 I k M \[Omega])+a m (-2-k+k^2-2 I M \[Omega]+4 I k M \[Omega])+M (2 I+I k-I k^2-4 M \[Omega]+4 k M \[Omega]))+2 \[Lambda] (-2+k^2-24 a m \[Omega]+2 I M \[Omega]+22 a^2 \[Omega]^2+k (-1+8 a m \[Omega]+4 I M \[Omega]-8 a^2 \[Omega]^2))) c[-1+k])/(2 k (2 \[Lambda]+\[Lambda]^2-12 \[Omega] (-a m+I M+a^2 \[Omega])));
  Table[c[j]=0,{j,-13,-1}];
  c[0]=1;

  (*FIXME: make this sum dependent on some accuracy threshold*)
  For[j=1, j<=20, j++,
    c[j]=recur/.k->j
  ];

  rp=M+Sqrt[M^2-a^2];
  rm=M-Sqrt[M^2-a^2];

  rs[r_]=r+(2M rp)/(rp-rm) Log[(r-rp)/(2M)]-(2M rm)/(rp-rm) Log[(r-rm)/(2M)];
  X[r_]:=Exp[I \[Omega] rs[r]]Sum[c[j](\[Omega] r)^-j,{j,0,20}];

  {X[r], X'[r], X''[r]}/.r->rout
 
 ]
 
(* Teukolsky equation*)
TeukolskyRadialEquation[s_, \[Lambda]_, m_, a_, \[Omega]_, r1_, {R_, dR_, d2R_}]:=Module[{M=1,\[CapitalDelta],r,K},

   \[CapitalDelta]=r^2-2M r+a^2;
   K=(r^2+a^2)\[Omega]-m a;
   
 ((K^2-2 I K (-M+r) s)/\[CapitalDelta]-\[Lambda]+4 I r s \[Omega]) R+\[CapitalDelta]^-s ((-2 M+2 r) \[CapitalDelta]^s (1+s) dR + \[CapitalDelta]^(1+s) d2R)

]
 
(* Sasaki-Nakamura equation for s=-2 perturbations *)
SasakiNakamuraEquation[\[Lambda]_, m_, a_,  \[Omega]_, r1_, {X_, dX_, d2X_}] := Module[{M=1,\[CapitalDelta],f,c0,c1,c2,c3,c4,\[Eta],K,V,\[Beta],\[Alpha],U1,G,F,U,r},

\[CapitalDelta]=r^2-2M r +a^2;

f=\[CapitalDelta]/(r^2+a^2);

c0 = -12I \[Omega] M+\[Lambda](\[Lambda]+2)-12 a \[Omega](a \[Omega]-m);
c1 = 8I a(3a \[Omega]-\[Lambda](a \[Omega]-m));
c2 = -24I a M(a \[Omega]-m)+12a^2 (1-2(a \[Omega]-m)^2);
c3 = 24I a^3 (a \[Omega]-m)-24 M a^2;
c4 = 12a^4;

\[Eta] = c0+c1/r+c2/r^2+c3/r^3+c4/r^4;

K=(r^2+a^2)\[Omega]-m a;
V=-((K^2+4I (r-M)K)/\[CapitalDelta])+8 I \[Omega] r+\[Lambda];
\[Beta]=2 \[CapitalDelta](-I K+r-M-2 \[CapitalDelta]/r);
\[Alpha]=-I K \[Beta]/\[CapitalDelta]^2+3 I D[K,r]+\[Lambda]+6 \[CapitalDelta]/r^2;
U1=V+\[CapitalDelta]^2/\[Beta] (D[2\[Alpha]+D[\[Beta],r]/\[CapitalDelta],r]-D[\[Eta],r]/\[Eta] (\[Alpha]+D[\[Beta],r]/\[CapitalDelta]));
G=-((2(r-M))/(r^2+a^2))+(r \[CapitalDelta])/(r^2+a^2)^2;

F=D[\[Eta],r]/\[Eta] \[CapitalDelta]/(r^2+a^2);
U=((\[CapitalDelta] U1)/(r^2+a^2)^2+G^2+(\[CapitalDelta] D[G,r])/(r^2+a^2)-F G);

f^2 d2X + f(D[f,r]-F)dX - U X/.r->r1

]
 
TeukolskyRadialUp[a_, \[Lambda]_, m_, \[Omega]_, r1_, Method -> "MST"] := Module[{},
	Print["MST method not yet implemented"]
]
 
SasakiNakamuraRadialUp[a_?NumericQ, \[Lambda]_?NumericQ, m_?NumericQ, \[Omega]_?NumericQ, r1_?NumericQ] := Module[{M=1, r, rout, outBC, XUp},

 rout = 1000;
 outBC = OuterBCs[\[Lambda], m, a, \[Omega], rout];
 XUp = X/.NDSolve[{SasakiNakamuraEquation[\[Lambda], m, a, \[Omega], r, {X[r],X'[r],X''[r]}]==0, X[rout]==outBC[[1]], X'[rout]==outBC[[2]]}, X, {r, rout, 10}][[1,1]];
  
 XUp 
]
 
TeukolskyRadialUp[s_, \[Lambda]_, m_, a_, \[Omega]_, r1_, Method->"Numerical"]/; s != -2 := Print["Only s=-2 is currently implemented for Method -> \"Numerical\""]
 
TeukolskyRadialUp[-2, \[Lambda]_?NumericQ, m_?IntegerQ, a_?NumericQ, \[Omega]_?NumericQ, r1_?NumericQ, Method -> "Numerical"]:=Module[{M=1, r, \[Chi]m, X, XUp, Rup, c0, c1, c2, c3 ,c4, K, \[Eta], \[CapitalDelta], \[Alpha], \[Beta], \[Chi]},

XUp = SasakiNakamuraRadialUp[a, \[Lambda], m, \[Omega], r1];

\[CapitalDelta]=r^2-2M r +a^2;


c0 = -12I \[Omega] M+\[Lambda](\[Lambda]+2)-12 a \[Omega](a \[Omega]-m);
c1 = 8I a(3a \[Omega]-\[Lambda](a \[Omega]-m));
c2 = -24I a M(a \[Omega]-m)+12a^2 (1-2(a \[Omega]-m)^2);
c3 = 24I a^3 (a \[Omega]-m)-24 M a^2;
c4 = 12a^4;

K=(r^2+a^2)\[Omega]-m a;

\[Eta] = c0+c1/r+c2/r^2+c3/r^3+c4/r^4; 

\[Beta] = 2 \[CapitalDelta](-I K+r-M-2 \[CapitalDelta]/r);
\[Alpha] = -I K \[Beta]/\[CapitalDelta]^2+3 I D[K,r] + \[Lambda] + 6 \[CapitalDelta]/r^2;
 
 
 \[Chi] = (X[r]\[CapitalDelta])/Sqrt[r^2+a^2];
 
 Rup = 1/\[Eta] ((\[Alpha] + D[\[Beta],r]/\[CapitalDelta])\[Chi] - \[Beta]/\[CapitalDelta] D[\[Chi],r]);

 {Rup, D[Rup,r]}/.{r->r1, X[r]->XUp[r1], X'[r] -> XUp'[r1], X''[r] -> XUp''[r1]}
 
 ]
  
 
End[];
EndPackage[];



