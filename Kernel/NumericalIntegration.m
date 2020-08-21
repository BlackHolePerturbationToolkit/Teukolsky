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


(* ::Section:: *)
(*Radial Solutions with HPS*)


(* ::Subsection::Closed:: *)
(*Useful Functions*)


f=1-2/r;
rs[r_]:=r+2 Log[r/2-1];
fr[r_]=1-2/r;


(* ::Subsection:: *)
(*Radial Bardeen-Press-Teukolsky Equation*)


SetAttributes[psi, {NumericFunction}];

psi[s_, \[Lambda]_, l_, m_, a_, \[Omega]_, "In", ndsolveopts___][rmax_?NumericQ] := psi[s, \[Lambda], m, a, \[Omega], "In", ndsolveopts][{Automatic, rmax}];
psi[s_, \[Lambda]_, l_, m_, a_, \[Omega]_, "Up", ndsolveopts___][rmin_?NumericQ] := psi[s, \[Lambda], m, a, \[Omega], "Up", ndsolveopts][{rmin, Automatic}];

psi[s_, \[Lambda]_, l_, m_, a_, \[Omega]_, bc_, ndsolveopts___][{rmin_, rmax_}] :=
 Module[{bcFunc, psiBC, dpsidrBC, rBC, rMin, rMax, H, soln, r},
    bcFunc = Lookup[<|"In" -> TeukolskyInBC, "Up" -> TeukolskyUpBC|>, bc];
    {psiBC, dpsidrBC, rBC} = bcFunc[s, \[Lambda], l, m, a, \[Omega], Lookup[{ndsolveopts}, WorkingPrecision, Precision[\[Omega]]]];
    If[bc === "In" && rmin === Automatic, rMin = rBC, rMin = rmin];
    If[bc === "Up" && rmax === Automatic, rMax = rBC, rMax = rmax];
    If[bc === "In", H = -1];
    If[bc === "Up", H = +1];
    soln = Integrator[s, \[Lambda], m, a, \[Omega], psiBC, dpsidrBC, rBC, rMin, rMax, H, ndsolveopts]
];

psi[s_, \[Lambda]_, l_, m_, a_, \[Omega]_, bc_, ndsolveopts___][All] :=
 Module[{bcFunc, psiBC, dpsidrBC, rBC, rMin, rMax, H, soln, r},
    bcFunc = Lookup[<|"In" -> TeukolskyInBC, "Up" -> TeukolskyUpBC|>, bc];
    {psiBC, dpsidrBC, rBC} = bcFunc[s, \[Lambda], l, m, a, \[Omega], Lookup[{ndsolveopts}, WorkingPrecision, Precision[\[Omega]]]];
    soln = AllIntegrator[s, \[Lambda], m, a, \[Omega], psiBC, dpsidrBC, rBC, H, ndsolveopts]
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


AllIntegrator[s_,\[Lambda]_,m_,a_,\[Omega]_,y1BC_,y2BC_,rBC_,potential_,H_?NumericQ,ndsolveopts___][rval:(_?NumericQ | {_?NumericQ..})] := Module[{y1,y2,r},
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

Derivative[n_][AllIntegrator[s_,\[Lambda]_,m_,a_,\[Omega]_,y1BC_,y2BC_,rBC_,potential_,H_?NumericQ,ndsolveopts___]][rval:(_?NumericQ | {_?NumericQ..})] := Module[{y1,y2,r},
	Quiet[NDSolveValue[
		{y1'[r]==y2[r],(((a^2-2 r+r^2) (2 a^2-2 r (1+s)-r^2 \[Lambda])-2 a (1+H) m r^2 (a^2+r^2) \[Omega]+2 I r^2 (-(1+H) (-a^2+r^2)+(1-H) r (a^2-2 r+r^2)) s \[Omega]+(1-H^2) r^2 (a^2+r^2)^2 \[Omega]^2-2 I a r (a^2-2 r+r^2) (m+a H \[Omega])) y1[r])/r^6+((2 (-a^2+r^2) (a^2-2 r+r^2))/(r^4 (a^2+r^2))-(2 (a^2-2 r+r^2) (a^2 (a^2-2 r+r^2)+(a^2+r^2) ((-1+r) r s-I r (a m+H (a^2+r^2) \[Omega]))))/(r^5 (a^2+r^2))) y2[r]+((a^2-2 r+r^2)^2 Derivative[1][y2][r])/r^4,y1[rBC]==y1BC,y2[rBC]==y2BC},
		Derivative[n][y1][rval],
		{r, Min[rBC,rval], Max[rBC,rval]},
		ndsolveopts,
		Method->"StiffnessSwitching",
		MaxSteps->Infinity,
		InterpolationOrder->All
		], NDSolveValue::precw]
	];


(* ::Subsection:: *)
(*Boundary Conditions*)


TeukolskyInBC[s1_Integer, \[Lambda]1_, l1_Integer, m1_Integer, a1_, \[Omega]1_, workingprecision_]:=
	Block[{s=s1,\[Lambda]=\[Lambda]1,l=l1,m=m1,a=a1,\[Omega]=\[Omega]1,R,r,res,dres,rin=2+10^-5},
		R = Teukolsky`TeukolskyRadial`TeukolskyRadial[s, l, m, a, \[Omega]];
		res = E^(I \[Omega] (r+((1+Sqrt[1-a^2]) Log[1/2 (-1-Sqrt[1-a^2]+r)]-(1-Sqrt[1-a^2]) Log[1/2 (-1+Sqrt[1-a^2]+r)])/Sqrt[1-a^2])) r ((-1-Sqrt[1-a^2]+r)/(-1+Sqrt[1-a^2]+r))^(-((I a m)/(2 Sqrt[1-a^2]))) (a^2-2 r+r^2)^s R["In"][r];
		dres =E^(I \[Omega] (r+((1+Sqrt[1-a^2]) Log[1/2 (-1-Sqrt[1-a^2]+r)]-(1-Sqrt[1-a^2]) Log[1/2 (-1+Sqrt[1-a^2]+r)])/Sqrt[1-a^2])) ((-1-Sqrt[1-a^2]+r)/(-1+Sqrt[1-a^2]+r))^(-((I a m)/(2 Sqrt[1-a^2]))) (a^2-2 r+r^2)^s R["In"][r]-(I a E^(I \[Omega] (r+((1+Sqrt[1-a^2]) Log[1/2 (-1-Sqrt[1-a^2]+r)]-(1-Sqrt[1-a^2]) Log[1/2 (-1+Sqrt[1-a^2]+r)])/Sqrt[1-a^2])) m r ((-1-Sqrt[1-a^2]+r)/(-1+Sqrt[1-a^2]+r))^(-1-(I a m)/(2 Sqrt[1-a^2])) (a^2-2 r+r^2)^s (-((-1-Sqrt[1-a^2]+r)/(-1+Sqrt[1-a^2]+r)^2)+1/(-1+Sqrt[1-a^2]+r)) R["In"][r])/(2 Sqrt[1-a^2])+E^(I \[Omega] (r+((1+Sqrt[1-a^2]) Log[1/2 (-1-Sqrt[1-a^2]+r)]-(1-Sqrt[1-a^2]) Log[1/2 (-1+Sqrt[1-a^2]+r)])/Sqrt[1-a^2])) r ((-1-Sqrt[1-a^2]+r)/(-1+Sqrt[1-a^2]+r))^(-((I a m)/(2 Sqrt[1-a^2]))) (-2+2 r) (a^2-2 r+r^2)^(-1+s) s R["In"][r]+I E^(I \[Omega] (r+((1+Sqrt[1-a^2]) Log[1/2 (-1-Sqrt[1-a^2]+r)]-(1-Sqrt[1-a^2]) Log[1/2 (-1+Sqrt[1-a^2]+r)])/Sqrt[1-a^2])) r ((-1-Sqrt[1-a^2]+r)/(-1+Sqrt[1-a^2]+r))^(-((I a m)/(2 Sqrt[1-a^2]))) (a^2+r^2) (a^2-2 r+r^2)^(-1+s) \[Omega] R["In"][r]+E^(I \[Omega] (r+((1+Sqrt[1-a^2]) Log[1/2 (-1-Sqrt[1-a^2]+r)]-(1-Sqrt[1-a^2]) Log[1/2 (-1+Sqrt[1-a^2]+r)])/Sqrt[1-a^2])) r ((-1-Sqrt[1-a^2]+r)/(-1+Sqrt[1-a^2]+r))^(-((I a m)/(2 Sqrt[1-a^2]))) (a^2-2 r+r^2)^s Derivative[1][R["In"]][r];
		{res,dres,rin}/.r->rin
	];

TeukolskyUpBC[s1_Integer, \[Lambda]1_, l1_Integer, m1_Integer, a1_, \[Omega]1_, workingprecision_]:=
	Block[{s=s1,\[Lambda]=\[Lambda]1,l=l1,m=m1,a=a1,\[Omega]=\[Omega]1,R,r,res,dres,rout},
		rout =100\[Omega]^-1;
		R = Teukolsky`TeukolskyRadial`TeukolskyRadial[s, l, m, a, \[Omega]];
		res = E^(-I \[Omega] (r+((1+Sqrt[1-a^2]) Log[1/2 (-1-Sqrt[1-a^2]+r)]-(1-Sqrt[1-a^2]) Log[1/2 (-1+Sqrt[1-a^2]+r)])/Sqrt[1-a^2])) r ((-1-Sqrt[1-a^2]+r)/(-1+Sqrt[1-a^2]+r))^(-((I a m)/(2 Sqrt[1-a^2]))) (a^2-2 r+r^2)^s R["Up"][r];
		dres = E^(-I \[Omega] (r+((1+Sqrt[1-a^2]) Log[1/2 (-1-Sqrt[1-a^2]+r)]-(1-Sqrt[1-a^2]) Log[1/2 (-1+Sqrt[1-a^2]+r)])/Sqrt[1-a^2])) ((-1-Sqrt[1-a^2]+r)/(-1+Sqrt[1-a^2]+r))^(-((I a m)/(2 Sqrt[1-a^2]))) (a^2-2 r+r^2)^s R["Up"][r]-(I a E^(-I \[Omega] (r+((1+Sqrt[1-a^2]) Log[1/2 (-1-Sqrt[1-a^2]+r)]-(1-Sqrt[1-a^2]) Log[1/2 (-1+Sqrt[1-a^2]+r)])/Sqrt[1-a^2])) m r ((-1-Sqrt[1-a^2]+r)/(-1+Sqrt[1-a^2]+r))^(-1-(I a m)/(2 Sqrt[1-a^2])) (a^2-2 r+r^2)^s (-((-1-Sqrt[1-a^2]+r)/(-1+Sqrt[1-a^2]+r)^2)+1/(-1+Sqrt[1-a^2]+r)) R["Up"][r])/(2 Sqrt[1-a^2])+E^(-I \[Omega] (r+((1+Sqrt[1-a^2]) Log[1/2 (-1-Sqrt[1-a^2]+r)]-(1-Sqrt[1-a^2]) Log[1/2 (-1+Sqrt[1-a^2]+r)])/Sqrt[1-a^2])) r ((-1-Sqrt[1-a^2]+r)/(-1+Sqrt[1-a^2]+r))^(-((I a m)/(2 Sqrt[1-a^2]))) (-2+2 r) (a^2-2 r+r^2)^(-1+s) s R["Up"][r]-I E^(-I \[Omega] (r+((1+Sqrt[1-a^2]) Log[1/2 (-1-Sqrt[1-a^2]+r)]-(1-Sqrt[1-a^2]) Log[1/2 (-1+Sqrt[1-a^2]+r)])/Sqrt[1-a^2])) r ((-1-Sqrt[1-a^2]+r)/(-1+Sqrt[1-a^2]+r))^(-((I a m)/(2 Sqrt[1-a^2]))) (a^2+r^2) (a^2-2 r+r^2)^(-1+s) \[Omega] R["Up"][r]+E^(-I \[Omega] (r+((1+Sqrt[1-a^2]) Log[1/2 (-1-Sqrt[1-a^2]+r)]-(1-Sqrt[1-a^2]) Log[1/2 (-1+Sqrt[1-a^2]+r)])/Sqrt[1-a^2])) r ((-1-Sqrt[1-a^2]+r)/(-1+Sqrt[1-a^2]+r))^(-((I a m)/(2 Sqrt[1-a^2]))) (a^2-2 r+r^2)^s Derivative[1][R["Up"]][r];
		{res,dres,rout}/.r->rout
	];



(* ::Section::Closed:: *)
(*End Package*)


End[]
EndPackage[];
