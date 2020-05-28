(* ::Package:: *)

(* ::Title:: *)
(*HyperboloidalSlicing*)


(* ::Section:: *)
(*Create Package*)


(* ::Subsection:: *)
(*Begin Package*)


BeginPackage["Teukolsky`HyperboloidalSlicing`"];


(* ::Subsection:: *)
(*Begin Private*)


Begin["`Private`"];


(* ::Section:: *)
(*Radial Solutions with HPS*)


(* ::Subsection:: *)
(*Useful Functions*)


f=1-2/r;
rs[r_]:=r+2 Log[r/2-1];


(* ::Subsection:: *)
(*Radial Bardeen-Press-Teukolsky Equation*)


SetAttributes[psi, {NumericFunction}];

psi[s_, l_, m_, \[Omega]_, "In", ndsolveopts___][rmax_?NumericQ] := psi[s, l, \[Omega], "In", ndsolveopts][{Automatic, rmax}];
psi[s_, l_, m_, \[Omega]_, "Up", ndsolveopts___][rmin_?NumericQ] := psi[s, l, \[Omega], "Up", ndsolveopts][{rmin, Automatic}];

psi[s_, l_, m_, \[Omega]_, bc_, ndsolveopts___][{rmin_, rmax_}] :=
 Module[{bcFunc, psiBC, dpsidrBC, rBC, rMin, rMax, H, soln},
    bcFunc = Lookup[<|"In" -> TeukolskyInBC, "Up" -> TeukolskyUpBC|>, bc];
    {psiBC, dpsidrBC, rBC} = bcFunc[s, l, m, \[Omega], Lookup[{ndsolveopts}, WorkingPrecision, Precision[\[Omega]]]];
    If[bc === "In" && rmin === Automatic, rMin = rBC, rMin = rmin, H = -1];
    If[bc === "Up" && rmax === Automatic, rMax = rBC, rMax = rmax, H = +1];
    soln = r^-(2s+1) f^-s Exp[I \[Omega] H rs[r]]Integrator[s, l, \[Omega], psiBC, dpsidrBC, rBC, rMin, rMax, H, ndsolveopts]
];

psi[s_, l_, m_, \[Omega]_, bc_, ndsolveopts___][All] :=
 Module[{bcFunc, psiBC, dpsidrBC, rBC, rMin, rMax, H, soln},
    bcFunc = Lookup[<|"In" -> TeukolskyInBC, "Up" -> TeukolskyUpBC|>, bc];
    {psiBC, dpsidrBC, rBC} = bcFunc[s, l, m, \[Omega], Lookup[{ndsolveopts}, WorkingPrecision, Precision[\[Omega]]]];
    soln = r^-(2s+1) f[r]^-s Exp[I \[Omega] H rs[r]]AllIntegrator[s, l, m, \[Omega], psiBC, dpsidrBC, rBC, H, ndsolveopts]
];

psi[s_, l_, m_, \[Omega]_, bc_, ndsolveopts___][None] := $Failed;


Integrator[s_,l_,m_,\[Omega]_,y1BC_,y2BC_,rBC_,rmin_?NumericQ,rmax_?NumericQ,H_?NumericQ,ndsolveopts___]:=Module[{y1,y2,r},
	Quiet[NDSolveValue[
		{y1'[r]==y2[r],(-(((1-2/r) (l (1+l)+(2 (1+s))/r-s (1+s)))/r^2)+(2 I (-1-H+(1-H) (1-2/r) r) s \[Omega])/r^2+(1-H^2) \[Omega]^2) y1[r]+(1-2/r) (2/r^2-(2 (-1+r) s)/r^2+2 I H \[Omega]) y2[r]+(1-2/r)^2 Derivative[1][y2][r]==0,y1[rBC]==y1BC,y2[rBC]==y2BC},
		y1,
		{r, rmin, rmax},
		ndsolveopts,
		Method->"StiffnessSwitching",
		MaxSteps->Infinity,
		InterpolationOrder->All
		], NDSolveValue::precw]
	];


AllIntegrator[s_,l_,m_,\[Omega]_,y1BC_,y2BC_,rBC_,potential_,H_?NumericQ,ndsolveopts___][rval:(_?NumericQ | {_?NumericQ..})] := Module[{y1,y2,r},
	Quiet[NDSolveValue[
		{y1'[r]==y2[r],(-(((1-2/r) (l (1+l)+(2 (1+s))/r-s (1+s)))/r^2)+(2 I (-1-H+(1-H) (1-2/r) r) s \[Omega])/r^2+(1-H^2) \[Omega]^2) y1[r]+(1-2/r) (2/r^2-(2 (-1+r) s)/r^2+2 I H \[Omega]) y2[r]+(1-2/r)^2 Derivative[1][y2][r]==0,y1[rBC]==y1BC,y2[rBC]==y2BC},
		y1[rval],
		{r, Min[rBC,rval], Max[rBC,rval]},
		ndsolveopts,
		Method->"StiffnessSwitching",
		MaxSteps->Infinity,
		InterpolationOrder->All
		], NDSolveValue::precw]
	];

Derivative[n_][AllIntegrator[s_,l_,m_,\[Omega]_,y1BC_,y2BC_,rBC_,potential_,H_?NumericQ,ndsolveopts___]][rval:(_?NumericQ | {_?NumericQ..})] := Module[{y1,y2,r},
	Quiet[NDSolveValue[
		{y1'[r]==y2[r],(-(((1-2/r) (l (1+l)+(2 (1+s))/r-s (1+s)))/r^2)+(2 I (-1-H+(1-H) (1-2/r) r) s \[Omega])/r^2+(1-H^2) \[Omega]^2) y1[r]+(1-2/r) (2/r^2-(2 (-1+r) s)/r^2+2 I H \[Omega]) y2[r]+(1-2/r)^2 Derivative[1][y2][r]==0,y1[rBC]==y1BC,y2[rBC]==y2BC},
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


TeukolskyInBC[s1_Integer, l1_Integer, m1_Integer, \[Omega]1_, workingprecision_]:=
	Block[{s=s1,l=l1,m=m1,\[Omega]=\[Omega]1,R,r,res,dres,rin=2+10^-5},
		
		R = TeukolskyRadial[s, l, m, 0, \[Omega]];
		res = E^(I \[Omega] (r-Log[4]+2 Log[-2+r])) (-2+r)^s r^(1+s) R["In"][r];
		dres =E^(I \[Omega] (r-Log[4]+2 Log[-2+r])) (-2+r)^(-1+s) r^s ((r+2 r s-2 (1+s)+I r^2 \[Omega]) R["In"][r]+(-2+r) r Derivative[1][R["In"]][r]);
		{res,dres,rin}/.r->rin
	];

TeukolskyUpBC[s1_Integer, l1_Integer, m1_Integer, \[Omega]1_, workingprecision_]:=
	Block[{s=s1,l=l1,m=m1,\[Omega]=\[Omega]1,R,r,res,dres,rout},
		rout =100\[Omega]^-1;
		R=TeukolskyRadial[s, l, m, 0, \[Omega]];
		res = E^(-I \[Omega] (r+2 Log[1/2 (-2+r)])) (-2+r)^s r^(1+s) R["Up"][r];
		dres = E^(-I \[Omega] (r+2 Log[1/2 (-2+r)])) (-2+r)^(-1+s) r^s ((r+2 r s-2 (1+s)-I r^2 \[Omega]) R["Up"][r]+(-2+r) r Derivative[1][R["Up"]][r]);
		{res,dres,rout}/.r->rout
	];



(* ::Section:: *)
(*End Package*)


End[]
EndPackage[];
