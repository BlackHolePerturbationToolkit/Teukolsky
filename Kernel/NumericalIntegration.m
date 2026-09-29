(* ::Package:: *)

(* ::Title:: *)
(*HyperboloidalSlicing*)


(* ::Section::Closed:: *)
(*Create Package*)


(* ::Subsection::Closed:: *)
(*Begin Package*)


BeginPackage["Teukolsky`NumericalIntegration`", {"Teukolsky`"}];


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

psi[s_, \[Lambda]_, l_, m_, a_, \[Omega]_, "In", amps_, \[Nu]_, ndsolveopts___][rmax_?NumericQ] := psi[s, \[Lambda], m, a, \[Omega], "In", amps, \[Nu], ndsolveopts][{Automatic, rmax}];
psi[s_, \[Lambda]_, l_, m_, a_, \[Omega]_, "Up", amps_, \[Nu]_, ndsolveopts___][rmin_?NumericQ] := psi[s, \[Lambda], m, a, \[Omega], "Up", amps, \[Nu], ndsolveopts][{rmin, Automatic}];

(* Boundary data: by default the MST solution near the horizon ("In") or at large radius ("Up"), converted to
   the HPS variables (the MST series are summed with a relative goal only, since an absolute AccuracyGoal would
   truncate them too early where the solution is small); alternatively "BoundaryData" -> <|bc -> {R, R', r}|> among the options supplies the
   value and derivative of the radial function at a chosen radius (used by TeukolskyPointParticleMode
   to start the integration over the orbit's radial range from an accurate MST solution at its edge). *)
boundaryData[s_, \[Lambda]_, l_, m_, a_, \[Omega]_, bc_, amps_, \[Nu]_, {rmin_, rmax_}, ndsolveopts___] :=
 Module[{bdata = Lookup[{ndsolveopts}, "BoundaryData", None], bcFunc, Rr, dRr, rb},
    If[AssociationQ[bdata] && KeyExistsQ[bdata, bc],
      {Rr, dRr, rb} = bdata[bc];
      Return[Lookup[<|"In" -> TeukolskyInBCFromValues, "Up" -> TeukolskyUpBCFromValues|>, bc][s, m, a, \[Omega], rb, Rr, dRr]]];
    bcFunc = Lookup[<|"In" -> TeukolskyInBC, "Up" -> TeukolskyUpBC|>, bc];
    bcFunc[s, \[Lambda], l, m, a, \[Omega], amps, \[Nu], Lookup[{ndsolveopts}, {WorkingPrecision, PrecisionGoal, AccuracyGoal}], If[bc === "In", rmin, rmax]]
 ];

(* Radius at which the MST solutions are evaluated to provide boundary data for the integration: the edge
   of the requested domain, otherwise r+ + 1/5 for "In" and r = 1000 for "Up".  The MST "Up" series is
   accurate at any radius beyond about r+ + 1, and the "In" solution (evaluated with precision padding,
   and from the Coulomb-type series at large radius) at any radius, so that integrating from the domain
   edges is as accurate as integrating inwards from a large radius or outwards from near the horizon, and
   cheaper.  The "In" integration over the full domain starts at r+ + 1/5: the MST "In" series is most
   accurate near the horizon, but integrating outwards from much closer to it is unstable for positive
   spin, where the "In" solution decays like Delta^(-s) and roundoff seeds the growing solution (for s = 2
   only about 8 digits survive from r+ + 10^-3 and 2 from r+ + 10^-5), and costs a digit or two for
   s <= 0 through the 1/Delta coefficients near the regular singular point; from r+ + 1/10 outwards there
   is no loss. *)
upBoundaryRadius[a_, rmax_] := If[NumericQ[rmax] && rmax < Infinity, Max[rmax, rp[a, 1] + 1], 1000];

inBoundaryRadius[a_, rmin_] := Module[{rh = rp[a, 1] + 1/5}, If[NumericQ[rmin] && rmin > rh, rmin, rh]];

rp[a_, M_] := M + Sqrt[M^2 - a^2];

ndsolveOptions[ndsolveopts___] := Sequence @@ DeleteCases[{ndsolveopts}, Rule["BoundaryData", _]];

psi[s_, \[Lambda]_, l_, m_, a_, \[Omega]_, bc_, amps_, \[Nu]_, ndsolveopts___][{rmin_, rmax_}] :=
 Module[{psiBC, dpsidrBC, rBC, rMin, rMax, H},
    {psiBC, dpsidrBC, rBC} = boundaryData[s, \[Lambda], l, m, a, \[Omega], bc, amps, \[Nu], {rmin, rmax}, ndsolveopts];
    If[bc === "In" && rmin === Automatic, rMin = rBC, rMin = rmin];
    If[bc === "Up" && rmax === Automatic, rMax = rBC, rMax = rmax];
    If[bc === "In", H = -1];
    If[bc === "Up", H = +1];
    Integrator[s, \[Lambda], m, a, \[Omega], psiBC, dpsidrBC, rBC, rMin, rMax, H, ndsolveOptions[ndsolveopts]]
];

psi[s_, \[Lambda]_, l_, m_, a_, \[Omega]_, bc_, amps_, \[Nu]_, ndsolveopts___][All] :=
 Module[{psiBC, dpsidrBC, rBC, rMin, rMax, H, bdata = Lookup[{ndsolveopts}, "BoundaryData", None]},
    If[bc === "In", H = -1];
    If[bc === "Up", H = +1];
    (* Without explicit boundary data, the "Up" solution on the full domain is integrated inwards from the
       MST solution at the outermost radius requested in each evaluation (integrating in from a large fixed
       radius loses several digits, whereas the MST "Up" series is accurate at any radius beyond r+ + 1) *)
    If[bc === "Up" && !(AssociationQ[bdata] && KeyExistsQ[bdata, "Up"]),
      Return[AllIntegratorUp[s, \[Lambda], l, m, a, \[Omega], amps, \[Nu], H, ndsolveopts]]];
    {psiBC, dpsidrBC, rBC} = boundaryData[s, \[Lambda], l, m, a, \[Omega], bc, amps, \[Nu], {Automatic, Automatic}, ndsolveopts];
    AllIntegrator[s, \[Lambda], m, a, \[Omega], psiBC, dpsidrBC, rBC, H, ndsolveOptions[ndsolveopts]]
];

psi[s_, \[Lambda]_, l_, m_, a_, \[Omega]_, bc_, amps_, \[Nu]_, ndsolveopts___][None] := $Failed;


Integrator[s_,\[Lambda]_,m_,a_,\[Omega]_,y1BC_,y2BC_,rBC_,rmin_?NumericQ,rmax_?NumericQ,H_?NumericQ,ndsolveopts___]:=Module[{Global`y1,Global`y2,Global`r,sol},
	Quiet[NDSolveValue[
		{Global`y1'[Global`r]==Global`y2[Global`r],(((a^2-2 Global`r+Global`r^2) (2 a^2-2 Global`r (1+s)-Global`r^2 \[Lambda])-2 a (1+H) m Global`r^2 (a^2+Global`r^2) \[Omega]+2 I Global`r^2 (-(1+H) (-a^2+Global`r^2)+(1-H) Global`r (a^2-2 Global`r+Global`r^2)) s \[Omega]+(1-H^2) Global`r^2 (a^2+Global`r^2)^2 \[Omega]^2-2 I a Global`r (a^2-2 Global`r+Global`r^2) (m+a H \[Omega])) Global`y1[Global`r])/Global`r^6+((2 (-a^2+Global`r^2) (a^2-2 Global`r+Global`r^2))/(Global`r^4 (a^2+Global`r^2))-(2 (a^2-2 Global`r+Global`r^2) (a^2 (a^2-2 Global`r+Global`r^2)+(a^2+Global`r^2) ((-1+Global`r) Global`r s-I Global`r (a m+H (a^2+Global`r^2) \[Omega]))))/(Global`r^5 (a^2+Global`r^2))) Global`y2[Global`r]+((a^2-2 Global`r+Global`r^2)^2 Derivative[1][Global`y2][Global`r])/Global`r^4==0,Global`y1[rBC]==y1BC,Global`y2[rBC]==y2BC},
		Global`y1,
		{Global`r, rmin, rmax},
		ndsolveopts,
		Method->"StiffnessSwitching",
		MaxSteps->Infinity,
		InterpolationOrder->All
		], NDSolveValue::precw]
	];


AllIntegrator[s_,\[Lambda]_,m_,a_,\[Omega]_,y1BC_,y2BC_,rBC_,H_?NumericQ,ndsolveopts___][rval:(_?NumericQ | {_?NumericQ..})] := Module[{Global`y1,Global`y2,Global`r,sol},
	Quiet[NDSolveValue[
		{Global`y1'[Global`r]==Global`y2[Global`r],(((a^2-2 Global`r+Global`r^2) (2 a^2-2 Global`r (1+s)-Global`r^2 \[Lambda])-2 a (1+H) m Global`r^2 (a^2+Global`r^2) \[Omega]+2 I Global`r^2 (-(1+H) (-a^2+Global`r^2)+(1-H) Global`r (a^2-2 Global`r+Global`r^2)) s \[Omega]+(1-H^2) Global`r^2 (a^2+Global`r^2)^2 \[Omega]^2-2 I a Global`r (a^2-2 Global`r+Global`r^2) (m+a H \[Omega])) Global`y1[Global`r])/Global`r^6+((2 (-a^2+Global`r^2) (a^2-2 Global`r+Global`r^2))/(Global`r^4 (a^2+Global`r^2))-(2 (a^2-2 Global`r+Global`r^2) (a^2 (a^2-2 Global`r+Global`r^2)+(a^2+Global`r^2) ((-1+Global`r) Global`r s-I Global`r (a m+H (a^2+Global`r^2) \[Omega]))))/(Global`r^5 (a^2+Global`r^2))) Global`y2[Global`r]+((a^2-2 Global`r+Global`r^2)^2 Derivative[1][Global`y2][Global`r])/Global`r^4==0,Global`y1[rBC]==y1BC,Global`y2[rBC]==y2BC},
		Global`y1[rval],
		{Global`r, Min[rBC,rval], Max[rBC,rval]},
		ndsolveopts,
		Method->"StiffnessSwitching",
		MaxSteps->Infinity,
		InterpolationOrder->All
		], NDSolveValue::precw]
	];

(* "Up" solution on the full domain: boundary data from the MST solution one unit beyond the outermost
   requested radius, then integration inwards over the requested radii. The boundary data are cached and
   reused while the requested radii stay below the cached radius and within 20 of it, so that repeated
   evaluations and decreasing sequences of radii do not re-evaluate the MST series. *)
$upBoundaryCache = <||>;

upBoundaryFromRadii[s_, \[Lambda]_, l_, m_, a_, \[Omega]_, amps_, \[Nu]_, rval_, ndsolveopts___] :=
 Module[{key = {s, \[Lambda], l, m, a, \[Omega], amps, \[Nu], Lookup[{ndsolveopts}, {WorkingPrecision, PrecisionGoal, AccuracyGoal}]}, rmax = Max[rval], cached, bc},
  cached = Lookup[$upBoundaryCache, Key[key], None];
  If[ListQ[cached] && rmax <= cached[[3]] && cached[[3]] <= rmax + 20, Return[cached]];
  bc = TeukolskyUpBC[s, \[Lambda], l, m, a, \[Omega], amps, \[Nu], Lookup[{ndsolveopts}, {WorkingPrecision, PrecisionGoal, AccuracyGoal}], rmax + 1];
  If[Length[$upBoundaryCache] >= 50, $upBoundaryCache = <||>];
  $upBoundaryCache[key] = bc
 ];

AllIntegratorUp[s_, \[Lambda]_, l_, m_, a_, \[Omega]_, amps_, \[Nu]_, H_?NumericQ, ndsolveopts___][rval:(_?NumericQ | {_?NumericQ..})] :=
 Module[{psiBC, dpsidrBC, rBC},
  {psiBC, dpsidrBC, rBC} = upBoundaryFromRadii[s, \[Lambda], l, m, a, \[Omega], amps, \[Nu], rval, ndsolveopts];
  AllIntegrator[s, \[Lambda], m, a, \[Omega], psiBC, dpsidrBC, rBC, H, ndsolveOptions[ndsolveopts]][rval]
 ];

Derivative[n_][AllIntegratorUp[s_, \[Lambda]_, l_, m_, a_, \[Omega]_, amps_, \[Nu]_, H_?NumericQ, ndsolveopts___]][rval:(_?NumericQ | {_?NumericQ..})] :=
 Module[{psiBC, dpsidrBC, rBC},
  {psiBC, dpsidrBC, rBC} = upBoundaryFromRadii[s, \[Lambda], l, m, a, \[Omega], amps, \[Nu], rval, ndsolveopts];
  Derivative[n][AllIntegrator[s, \[Lambda], m, a, \[Omega], psiBC, dpsidrBC, rBC, H, ndsolveOptions[ndsolveopts]]][rval]
 ];

Derivative[n_][AllIntegrator[s_,\[Lambda]_,m_,a_,\[Omega]_,y1BC_,y2BC_,rBC_,H_?NumericQ,ndsolveopts___]][rval:(_?NumericQ | {_?NumericQ..})] := Module[{Global`y1,Global`y2,Global`r,sol},
	Quiet[NDSolveValue[
		{Global`y1'[Global`r]==Global`y2[Global`r],(((a^2-2 Global`r+Global`r^2) (2 a^2-2 Global`r (1+s)-Global`r^2 \[Lambda])-2 a (1+H) m Global`r^2 (a^2+Global`r^2) \[Omega]+2 I Global`r^2 (-(1+H) (-a^2+Global`r^2)+(1-H) Global`r (a^2-2 Global`r+Global`r^2)) s \[Omega]+(1-H^2) Global`r^2 (a^2+Global`r^2)^2 \[Omega]^2-2 I a Global`r (a^2-2 Global`r+Global`r^2) (m+a H \[Omega])) Global`y1[Global`r])/Global`r^6+((2 (-a^2+Global`r^2) (a^2-2 Global`r+Global`r^2))/(Global`r^4 (a^2+Global`r^2))-(2 (a^2-2 Global`r+Global`r^2) (a^2 (a^2-2 Global`r+Global`r^2)+(a^2+Global`r^2) ((-1+Global`r) Global`r s-I Global`r (a m+H (a^2+Global`r^2) \[Omega]))))/(Global`r^5 (a^2+Global`r^2))) Global`y2[Global`r]+((a^2-2 Global`r+Global`r^2)^2 Derivative[1][Global`y2][Global`r])/Global`r^4==0,Global`y1[rBC]==y1BC,Global`y2[rBC]==y2BC},
		Derivative[n][Global`y1][rval],
		{Global`r, Min[rBC,rval], Max[rBC,rval]},
		ndsolveopts,
		Method->"StiffnessSwitching",
		MaxSteps->Infinity,
		InterpolationOrder->All
		], NDSolveValue::precw]
	];


(* ::Subsection::Closed:: *)
(*Boundary Conditions*)


TeukolskyInBC[s_Integer, \[Lambda]_, l_Integer, m_Integer, a_, \[Omega]_, amps_, \[Nu]_, {wp_, prec_, acc_}, rmin_:Automatic]:=
 Module[{R, res, dres, r, Rr, dRr},
        R = Teukolsky`TeukolskyRadial[s, l, m, a, \[Omega], "BoundaryConditions" -> "In", "Amplitudes" -> amps, "Eigenvalue" -> \[Lambda], "RenormalizedAngularMomentum" -> \[Nu], Method -> "MST", WorkingPrecision -> wp, PrecisionGoal -> prec, AccuracyGoal -> Infinity];
        r = inBoundaryRadius[a, rmin];
  	  Rr = R[r];
		dRr = R'[r];
		TeukolskyInBCFromValues[s, m, a, \[Omega], r, Rr, dRr]
	];

(* HPS boundary data {psi, dpsi/dr, r} from the radial function value and derivative at r *)
TeukolskyInBCFromValues[s_Integer, m_Integer, a_, \[Omega]_, r_, Rr_, dRr_] :=
 Module[{res, dres},
		res = E^(I \[Omega] (r+((1+Sqrt[1-a^2]) Log[1/2 (-1-Sqrt[1-a^2]+r)]-(1-Sqrt[1-a^2]) Log[1/2 (-1+Sqrt[1-a^2]+r)])/Sqrt[1-a^2])) r ((-1-Sqrt[1-a^2]+r)/(-1+Sqrt[1-a^2]+r))^(-((I a m)/(2 Sqrt[1-a^2]))) (a^2-2 r+r^2)^s Rr;
		dres =E^(I \[Omega] (r+((1+Sqrt[1-a^2]) Log[1/2 (-1-Sqrt[1-a^2]+r)]-(1-Sqrt[1-a^2]) Log[1/2 (-1+Sqrt[1-a^2]+r)])/Sqrt[1-a^2])) ((-1-Sqrt[1-a^2]+r)/(-1+Sqrt[1-a^2]+r))^(-((I a m)/(2 Sqrt[1-a^2]))) (a^2-2 r+r^2)^s Rr-(I a E^(I \[Omega] (r+((1+Sqrt[1-a^2]) Log[1/2 (-1-Sqrt[1-a^2]+r)]-(1-Sqrt[1-a^2]) Log[1/2 (-1+Sqrt[1-a^2]+r)])/Sqrt[1-a^2])) m r ((-1-Sqrt[1-a^2]+r)/(-1+Sqrt[1-a^2]+r))^(-1-(I a m)/(2 Sqrt[1-a^2])) (a^2-2 r+r^2)^s (-((-1-Sqrt[1-a^2]+r)/(-1+Sqrt[1-a^2]+r)^2)+1/(-1+Sqrt[1-a^2]+r)) Rr)/(2 Sqrt[1-a^2])+E^(I \[Omega] (r+((1+Sqrt[1-a^2]) Log[1/2 (-1-Sqrt[1-a^2]+r)]-(1-Sqrt[1-a^2]) Log[1/2 (-1+Sqrt[1-a^2]+r)])/Sqrt[1-a^2])) r ((-1-Sqrt[1-a^2]+r)/(-1+Sqrt[1-a^2]+r))^(-((I a m)/(2 Sqrt[1-a^2]))) (-2+2 r) (a^2-2 r+r^2)^(-1+s) s Rr+I E^(I \[Omega] (r+((1+Sqrt[1-a^2]) Log[1/2 (-1-Sqrt[1-a^2]+r)]-(1-Sqrt[1-a^2]) Log[1/2 (-1+Sqrt[1-a^2]+r)])/Sqrt[1-a^2])) r ((-1-Sqrt[1-a^2]+r)/(-1+Sqrt[1-a^2]+r))^(-((I a m)/(2 Sqrt[1-a^2]))) (a^2+r^2) (a^2-2 r+r^2)^(-1+s) \[Omega] Rr+E^(I \[Omega] (r+((1+Sqrt[1-a^2]) Log[1/2 (-1-Sqrt[1-a^2]+r)]-(1-Sqrt[1-a^2]) Log[1/2 (-1+Sqrt[1-a^2]+r)])/Sqrt[1-a^2])) r ((-1-Sqrt[1-a^2]+r)/(-1+Sqrt[1-a^2]+r))^(-((I a m)/(2 Sqrt[1-a^2]))) (a^2-2 r+r^2)^s dRr;
				{res,dres,r}
	];

TeukolskyUpBC[s_Integer, \[Lambda]_, l_Integer, m_Integer, a_, \[Omega]_, amps_, \[Nu]_, {wp_, prec_, acc_}, rmax_:Automatic]:=
 Module[{R, res, dres, r, Rr, dRr},
		r = upBoundaryRadius[a, rmax];
        R = Teukolsky`TeukolskyRadial[s, l, m, a, \[Omega], "BoundaryConditions" -> "Up", "Amplitudes" -> amps, "Eigenvalue" -> \[Lambda], "RenormalizedAngularMomentum" -> \[Nu], Method -> "MST", WorkingPrecision -> wp, PrecisionGoal -> prec, AccuracyGoal -> Infinity];
  	  Rr = R[r];
		dRr = R'[r];
		TeukolskyUpBCFromValues[s, m, a, \[Omega], r, Rr, dRr]
	];

(* HPS boundary data {psi, dpsi/dr, r} from the radial function value and derivative at r *)
TeukolskyUpBCFromValues[s_Integer, m_Integer, a_, \[Omega]_, r_, Rr_, dRr_] :=
 Module[{res, dres},
		res = E^(-I \[Omega] (r+((1+Sqrt[1-a^2]) Log[1/2 (-1-Sqrt[1-a^2]+r)]-(1-Sqrt[1-a^2]) Log[1/2 (-1+Sqrt[1-a^2]+r)])/Sqrt[1-a^2])) r ((-1-Sqrt[1-a^2]+r)/(-1+Sqrt[1-a^2]+r))^(-((I a m)/(2 Sqrt[1-a^2]))) (a^2-2 r+r^2)^s Rr;
		dres = E^(-I \[Omega] (r+((1+Sqrt[1-a^2]) Log[1/2 (-1-Sqrt[1-a^2]+r)]-(1-Sqrt[1-a^2]) Log[1/2 (-1+Sqrt[1-a^2]+r)])/Sqrt[1-a^2])) ((-1-Sqrt[1-a^2]+r)/(-1+Sqrt[1-a^2]+r))^(-((I a m)/(2 Sqrt[1-a^2]))) (a^2-2 r+r^2)^s Rr-(I a E^(-I \[Omega] (r+((1+Sqrt[1-a^2]) Log[1/2 (-1-Sqrt[1-a^2]+r)]-(1-Sqrt[1-a^2]) Log[1/2 (-1+Sqrt[1-a^2]+r)])/Sqrt[1-a^2])) m r ((-1-Sqrt[1-a^2]+r)/(-1+Sqrt[1-a^2]+r))^(-1-(I a m)/(2 Sqrt[1-a^2])) (a^2-2 r+r^2)^s (-((-1-Sqrt[1-a^2]+r)/(-1+Sqrt[1-a^2]+r)^2)+1/(-1+Sqrt[1-a^2]+r)) Rr)/(2 Sqrt[1-a^2])+E^(-I \[Omega] (r+((1+Sqrt[1-a^2]) Log[1/2 (-1-Sqrt[1-a^2]+r)]-(1-Sqrt[1-a^2]) Log[1/2 (-1+Sqrt[1-a^2]+r)])/Sqrt[1-a^2])) r ((-1-Sqrt[1-a^2]+r)/(-1+Sqrt[1-a^2]+r))^(-((I a m)/(2 Sqrt[1-a^2]))) (-2+2 r) (a^2-2 r+r^2)^(-1+s) s Rr-I E^(-I \[Omega] (r+((1+Sqrt[1-a^2]) Log[1/2 (-1-Sqrt[1-a^2]+r)]-(1-Sqrt[1-a^2]) Log[1/2 (-1+Sqrt[1-a^2]+r)])/Sqrt[1-a^2])) r ((-1-Sqrt[1-a^2]+r)/(-1+Sqrt[1-a^2]+r))^(-((I a m)/(2 Sqrt[1-a^2]))) (a^2+r^2) (a^2-2 r+r^2)^(-1+s) \[Omega] Rr+E^(-I \[Omega] (r+((1+Sqrt[1-a^2]) Log[1/2 (-1-Sqrt[1-a^2]+r)]-(1-Sqrt[1-a^2]) Log[1/2 (-1+Sqrt[1-a^2]+r)])/Sqrt[1-a^2])) r ((-1-Sqrt[1-a^2]+r)/(-1+Sqrt[1-a^2]+r))^(-((I a m)/(2 Sqrt[1-a^2]))) (a^2-2 r+r^2)^s dRr;
				{res,dres,r}
	];



(* ::Section::Closed:: *)
(*End Package*)


End[]
EndPackage[];
