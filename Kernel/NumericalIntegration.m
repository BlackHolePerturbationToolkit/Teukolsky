(* ::Package:: *)

(* ::Title:: *)
(*HyperboloidalSlicing*)


(* ::Section::Closed:: *)
(*Create Package*)


(* ::Subsection::Closed:: *)
(*Begin Package*)


BeginPackage["Teukolsky`NumericalIntegration`"];


(* ::Subsection::Closed:: *)
(*Begin Private*)


Begin["`Private`"];


(* ::Section::Closed:: *)
(*Radial Solutions with HPS*)


(* ::Subsection::Closed:: *)
(*Useful Functions*)


f=1-2/r;
rs[r_]:=r+2 Log[r/2-1];
fr[r_]=1-2/r;


(* ::Subsection::Closed:: *)
(*Radial Bardeen-Press-Teukolsky Equation*)


SetAttributes[psi, {NumericFunction}];

psi[s_, \[Lambda]_, l_, m_, a_, \[Omega]_, "In", ndsolveopts___][rmax_?NumericQ] := psi[s, \[Lambda], m, a, \[Omega], "In", ndsolveopts][{Automatic, rmax}];
psi[s_, \[Lambda]_, l_, m_, a_, \[Omega]_, "Up", ndsolveopts___][rmin_?NumericQ] := psi[s, \[Lambda], m, a, \[Omega], "Up", ndsolveopts][{rmin, Automatic}];

psi[s_, \[Lambda]_, l_, m_, a_, \[Omega]_, bc_, ndsolveopts___][{rmin_, rmax_}] :=
 Module[{bcFunc, psiBC, dpsidrBC, rBC, rMin, rMax, H, amp},
    bcFunc = Lookup[<|"In" -> TeukolskyInBC, "Up" -> TeukolskyUpBC|>, bc];
    {psiBC, dpsidrBC, rBC, amp} = bcFunc[s, \[Lambda], l, m, a, \[Omega], Lookup[{ndsolveopts}, {WorkingPrecision, PrecisionGoal, AccuracyGoal}]];
    If[bc === "In" && rmin === Automatic, rMin = rBC, rMin = rmin];
    If[bc === "Up" && rmax === Automatic, rMax = rBC, rMax = rmax];
    If[bc === "In", H = -1];
    If[bc === "Up", H = +1];
    {Integrator[s, \[Lambda], m, a, \[Omega], psiBC, dpsidrBC, rBC, rMin, rMax, H, ndsolveopts], amp}
];

psi[s_, \[Lambda]_, l_, m_, a_, \[Omega]_, bc_, ndsolveopts___][All] :=
 Module[{bcFunc, psiBC, dpsidrBC, rBC, rMin, rMax, H, amp},
    bcFunc = Lookup[<|"In" -> TeukolskyInBC, "Up" -> TeukolskyUpBC|>, bc];
    {psiBC, dpsidrBC, rBC, amp} = bcFunc[s, \[Lambda], l, m, a, \[Omega], Lookup[{ndsolveopts}, {WorkingPrecision, PrecisionGoal, AccuracyGoal}]];
    If[bc === "In", H = -1];
    If[bc === "Up", H = +1];
    {AllIntegrator[s, \[Lambda], m, a, \[Omega], psiBC, dpsidrBC, rBC, H, ndsolveopts], amp}
];

psi[s_, \[Lambda]_, l_, m_, \[Omega]_, bc_, ndsolveopts___][None] := $Failed;


Integrator[s_,\[Lambda]_,m_,a_,\[Omega]_,y1BC_,y2BC_,rBC_,rmin_?NumericQ,rmax_?NumericQ,H_?NumericQ,ndsolveopts___]:=Module[{y1,y2,r},
	Quiet[NDSolveValue[
		{y1'[r]==y2[r],(((a^2-2 r+r^2) (2 a^2-2 r (1+s)-r^2 \[Lambda])-2 a (1+H) m r^2 (a^2+r^2) \[Omega]+2 I r^2 (-(1+H) (-a^2+r^2)+(1-H) r (a^2-2 r+r^2)) s \[Omega]+(1-H^2) r^2 (a^2+r^2)^2 \[Omega]^2-2 I a r (a^2-2 r+r^2) (m+a H \[Omega])) y1[r])/r^6+((2 (-a^2+r^2) (a^2-2 r+r^2))/(r^4 (a^2+r^2))-(2 (a^2-2 r+r^2) (a^2 (a^2-2 r+r^2)+(a^2+r^2) ((-1+r) r s-I r (a m+H (a^2+r^2) \[Omega]))))/(r^5 (a^2+r^2))) y2[r]+((a^2-2 r+r^2)^2 Derivative[1][y2][r])/r^4==0,y1[rBC]==y1BC,y2[rBC]==y2BC},
		y1,
		{r, rmin, rmax},
		ndsolveopts,
		Method->"StiffnessSwitching",
		MaxSteps->Infinity,
		InterpolationOrder->All
		], NDSolveValue::precw]
	];


AllIntegrator[s_,\[Lambda]_,m_,a_,\[Omega]_,y1BC_,y2BC_,rBC_,H_?NumericQ,ndsolveopts___][rval:(_?NumericQ | {_?NumericQ..})] := Module[{y1,y2,r},
	Quiet[NDSolveValue[
		{y1'[r]==y2[r],(((a^2-2 r+r^2) (2 a^2-2 r (1+s)-r^2 \[Lambda])-2 a (1+H) m r^2 (a^2+r^2) \[Omega]+2 I r^2 (-(1+H) (-a^2+r^2)+(1-H) r (a^2-2 r+r^2)) s \[Omega]+(1-H^2) r^2 (a^2+r^2)^2 \[Omega]^2-2 I a r (a^2-2 r+r^2) (m+a H \[Omega])) y1[r])/r^6+((2 (-a^2+r^2) (a^2-2 r+r^2))/(r^4 (a^2+r^2))-(2 (a^2-2 r+r^2) (a^2 (a^2-2 r+r^2)+(a^2+r^2) ((-1+r) r s-I r (a m+H (a^2+r^2) \[Omega]))))/(r^5 (a^2+r^2))) y2[r]+((a^2-2 r+r^2)^2 Derivative[1][y2][r])/r^4==0,y1[rBC]==y1BC,y2[rBC]==y2BC},
		y1[rval],
		{r, Min[rBC,rval], Max[rBC,rval]},
		ndsolveopts,
		Method->"StiffnessSwitching",
		MaxSteps->Infinity,
		InterpolationOrder->All
		], NDSolveValue::precw]
	];

Derivative[n_][AllIntegrator[s_,\[Lambda]_,m_,a_,\[Omega]_,y1BC_,y2BC_,rBC_,H_?NumericQ,ndsolveopts___]][rval:(_?NumericQ | {_?NumericQ..})] := Module[{y1,y2,r},
	Quiet[NDSolveValue[
		{y1'[r]==y2[r],(((a^2-2 r+r^2) (2 a^2-2 r (1+s)-r^2 \[Lambda])-2 a (1+H) m r^2 (a^2+r^2) \[Omega]+2 I r^2 (-(1+H) (-a^2+r^2)+(1-H) r (a^2-2 r+r^2)) s \[Omega]+(1-H^2) r^2 (a^2+r^2)^2 \[Omega]^2-2 I a r (a^2-2 r+r^2) (m+a H \[Omega])) y1[r])/r^6+((2 (-a^2+r^2) (a^2-2 r+r^2))/(r^4 (a^2+r^2))-(2 (a^2-2 r+r^2) (a^2 (a^2-2 r+r^2)+(a^2+r^2) ((-1+r) r s-I r (a m+H (a^2+r^2) \[Omega]))))/(r^5 (a^2+r^2))) y2[r]+((a^2-2 r+r^2)^2 Derivative[1][y2][r])/r^4==0,y1[rBC]==y1BC,y2[rBC]==y2BC},
		Derivative[n][y1][rval],
		{r, Min[rBC,rval], Max[rBC,rval]},
		ndsolveopts,
		Method->"StiffnessSwitching",
		MaxSteps->Infinity,
		InterpolationOrder->All
		], NDSolveValue::precw]
	];


(* ::Subsection::Closed:: *)
(*Boundary Conditions*)


TeukolskyInBC[s_Integer, \[Lambda]_, l_Integer, m_Integer, a_, \[Omega]_, {wp_, prec_, acc_}]:=
 Module[{R, amp, res, dres, r=2+10^-5, Rr, dRr},
        If[wp === MachinePrecision,
		  R = Teukolsky`TeukolskyRadial`TeukolskyRadial[s, l, m, SetPrecision[a, 32], SetPrecision[\[Omega], 32], "BoundaryConditions" -> "In", Method -> "MST", PrecisionGoal -> prec, AccuracyGoal -> acc];
  		amp = N[R["Amplitudes"]];,
		  R = Teukolsky`TeukolskyRadial`TeukolskyRadial[s, l, m, a, \[Omega], "BoundaryConditions" -> "In", Method -> "MST", WorkingPrecision -> wp, PrecisionGoal -> prec, AccuracyGoal -> acc];
  		amp = R["Amplitudes"];
		];
		Rr = R[r];
		dRr = R'[r];
		res = E^(I \[Omega] (r+((1+Sqrt[1-a^2]) Log[1/2 (-1-Sqrt[1-a^2]+r)]-(1-Sqrt[1-a^2]) Log[1/2 (-1+Sqrt[1-a^2]+r)])/Sqrt[1-a^2])) r ((-1-Sqrt[1-a^2]+r)/(-1+Sqrt[1-a^2]+r))^(-((I a m)/(2 Sqrt[1-a^2]))) (a^2-2 r+r^2)^s Rr;
		dres =E^(I \[Omega] (r+((1+Sqrt[1-a^2]) Log[1/2 (-1-Sqrt[1-a^2]+r)]-(1-Sqrt[1-a^2]) Log[1/2 (-1+Sqrt[1-a^2]+r)])/Sqrt[1-a^2])) ((-1-Sqrt[1-a^2]+r)/(-1+Sqrt[1-a^2]+r))^(-((I a m)/(2 Sqrt[1-a^2]))) (a^2-2 r+r^2)^s Rr-(I a E^(I \[Omega] (r+((1+Sqrt[1-a^2]) Log[1/2 (-1-Sqrt[1-a^2]+r)]-(1-Sqrt[1-a^2]) Log[1/2 (-1+Sqrt[1-a^2]+r)])/Sqrt[1-a^2])) m r ((-1-Sqrt[1-a^2]+r)/(-1+Sqrt[1-a^2]+r))^(-1-(I a m)/(2 Sqrt[1-a^2])) (a^2-2 r+r^2)^s (-((-1-Sqrt[1-a^2]+r)/(-1+Sqrt[1-a^2]+r)^2)+1/(-1+Sqrt[1-a^2]+r)) Rr)/(2 Sqrt[1-a^2])+E^(I \[Omega] (r+((1+Sqrt[1-a^2]) Log[1/2 (-1-Sqrt[1-a^2]+r)]-(1-Sqrt[1-a^2]) Log[1/2 (-1+Sqrt[1-a^2]+r)])/Sqrt[1-a^2])) r ((-1-Sqrt[1-a^2]+r)/(-1+Sqrt[1-a^2]+r))^(-((I a m)/(2 Sqrt[1-a^2]))) (-2+2 r) (a^2-2 r+r^2)^(-1+s) s Rr+I E^(I \[Omega] (r+((1+Sqrt[1-a^2]) Log[1/2 (-1-Sqrt[1-a^2]+r)]-(1-Sqrt[1-a^2]) Log[1/2 (-1+Sqrt[1-a^2]+r)])/Sqrt[1-a^2])) r ((-1-Sqrt[1-a^2]+r)/(-1+Sqrt[1-a^2]+r))^(-((I a m)/(2 Sqrt[1-a^2]))) (a^2+r^2) (a^2-2 r+r^2)^(-1+s) \[Omega] Rr+E^(I \[Omega] (r+((1+Sqrt[1-a^2]) Log[1/2 (-1-Sqrt[1-a^2]+r)]-(1-Sqrt[1-a^2]) Log[1/2 (-1+Sqrt[1-a^2]+r)])/Sqrt[1-a^2])) r ((-1-Sqrt[1-a^2]+r)/(-1+Sqrt[1-a^2]+r))^(-((I a m)/(2 Sqrt[1-a^2]))) (a^2-2 r+r^2)^s dRr;
		{res,dres,r,amp}
	];

TeukolskyUpBC[s_Integer, \[Lambda]_, l_Integer, m_Integer, a_, \[Omega]_, {wp_, prec_, acc_}]:=
 Module[{R, amp, res, dres, r, Rr, dRr},
		r = 1000;
        If[wp === MachinePrecision,
		  R = Teukolsky`TeukolskyRadial`TeukolskyRadial[s, l, m, SetPrecision[a, 32], SetPrecision[\[Omega], 32], "BoundaryConditions" -> "Up", Method -> "MST", PrecisionGoal -> prec, AccuracyGoal -> acc];
  		amp = N[R["Amplitudes"]];,
		  R = Teukolsky`TeukolskyRadial`TeukolskyRadial[s, l, m, a, \[Omega], "BoundaryConditions" -> "Up", Method -> "MST", WorkingPrecision -> wp, PrecisionGoal -> prec, AccuracyGoal -> acc];
  		amp = R["Amplitudes"];
		];
		Rr = R[r];
		dRr = R'[r];
		res = E^(-I \[Omega] (r+((1+Sqrt[1-a^2]) Log[1/2 (-1-Sqrt[1-a^2]+r)]-(1-Sqrt[1-a^2]) Log[1/2 (-1+Sqrt[1-a^2]+r)])/Sqrt[1-a^2])) r ((-1-Sqrt[1-a^2]+r)/(-1+Sqrt[1-a^2]+r))^(-((I a m)/(2 Sqrt[1-a^2]))) (a^2-2 r+r^2)^s Rr;
		dres = E^(-I \[Omega] (r+((1+Sqrt[1-a^2]) Log[1/2 (-1-Sqrt[1-a^2]+r)]-(1-Sqrt[1-a^2]) Log[1/2 (-1+Sqrt[1-a^2]+r)])/Sqrt[1-a^2])) ((-1-Sqrt[1-a^2]+r)/(-1+Sqrt[1-a^2]+r))^(-((I a m)/(2 Sqrt[1-a^2]))) (a^2-2 r+r^2)^s Rr-(I a E^(-I \[Omega] (r+((1+Sqrt[1-a^2]) Log[1/2 (-1-Sqrt[1-a^2]+r)]-(1-Sqrt[1-a^2]) Log[1/2 (-1+Sqrt[1-a^2]+r)])/Sqrt[1-a^2])) m r ((-1-Sqrt[1-a^2]+r)/(-1+Sqrt[1-a^2]+r))^(-1-(I a m)/(2 Sqrt[1-a^2])) (a^2-2 r+r^2)^s (-((-1-Sqrt[1-a^2]+r)/(-1+Sqrt[1-a^2]+r)^2)+1/(-1+Sqrt[1-a^2]+r)) Rr)/(2 Sqrt[1-a^2])+E^(-I \[Omega] (r+((1+Sqrt[1-a^2]) Log[1/2 (-1-Sqrt[1-a^2]+r)]-(1-Sqrt[1-a^2]) Log[1/2 (-1+Sqrt[1-a^2]+r)])/Sqrt[1-a^2])) r ((-1-Sqrt[1-a^2]+r)/(-1+Sqrt[1-a^2]+r))^(-((I a m)/(2 Sqrt[1-a^2]))) (-2+2 r) (a^2-2 r+r^2)^(-1+s) s Rr-I E^(-I \[Omega] (r+((1+Sqrt[1-a^2]) Log[1/2 (-1-Sqrt[1-a^2]+r)]-(1-Sqrt[1-a^2]) Log[1/2 (-1+Sqrt[1-a^2]+r)])/Sqrt[1-a^2])) r ((-1-Sqrt[1-a^2]+r)/(-1+Sqrt[1-a^2]+r))^(-((I a m)/(2 Sqrt[1-a^2]))) (a^2+r^2) (a^2-2 r+r^2)^(-1+s) \[Omega] Rr+E^(-I \[Omega] (r+((1+Sqrt[1-a^2]) Log[1/2 (-1-Sqrt[1-a^2]+r)]-(1-Sqrt[1-a^2]) Log[1/2 (-1+Sqrt[1-a^2]+r)])/Sqrt[1-a^2])) r ((-1-Sqrt[1-a^2]+r)/(-1+Sqrt[1-a^2]+r))^(-((I a m)/(2 Sqrt[1-a^2]))) (a^2-2 r+r^2)^s dRr;
		{res,dres,r,amp}
	];



(* ::Section::Closed:: *)
(*End Package*)


End[]
EndPackage[];
