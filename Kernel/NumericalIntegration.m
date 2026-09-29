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
    bcFunc[s, \[Lambda], l, m, a, \[Omega], amps, \[Nu], Lookup[{ndsolveopts}, {WorkingPrecision, PrecisionGoal, AccuracyGoal}], If[bc === "In", rmin, rmax], Lookup[{ndsolveopts}, "BoundaryMethod", "MST"]]
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

ndsolveOptions[ndsolveopts___] := Sequence @@ DeleteCases[{ndsolveopts}, Rule["BoundaryData" | "BoundaryMethod", _]];

psi[s_, \[Lambda]_, l_, m_, a_, \[Omega]_, bc_, amps_, \[Nu]_, ndsolveopts___][{rmin_, rmax_}] :=
 Module[{psiBC, dpsidrBC, rBC, rMin, rMax, H},
    {psiBC, dpsidrBC, rBC} = boundaryData[s, \[Lambda], l, m, a, \[Omega], bc, amps, \[Nu], {rmin, rmax}, ndsolveopts];
    (* series boundary data sit at r+ + 1/5 ("In") or at the radius where the large-r series converges ("Up"),
       which may lie outside the requested domain: the integration then covers both *)
    If[bc === "In" && rmin === Automatic, rMin = rBC, rMin = Min[rmin, rBC]];
    If[bc === "Up" && rmax === Automatic, rMax = rBC, rMax = Max[rmax, rBC]];
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
      If[Lookup[{ndsolveopts}, "BoundaryMethod", "MST"] === "Series",
        With[{sol = upSeriesSolutionBuild[s, \[Lambda], m, a, \[Omega], ndsolveopts]}, If[sol =!= $Failed, Return[sol]]]];
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
 Module[{method = Lookup[{ndsolveopts}, "BoundaryMethod", "MST"], key, rmax = Max[rval], cached, bc},
  key = {s, \[Lambda], l, m, a, \[Omega], amps, \[Nu], Lookup[{ndsolveopts}, {WorkingPrecision, PrecisionGoal, AccuracyGoal}], method};
  cached = Lookup[$upBoundaryCache, Key[key], None];
  (* series data are placed where the series converges, possibly far out; the spin of the integration is then
     chosen for stability (see TeukolskyRadialNumericalIntegration), so they are reused for any radius below *)
  If[ListQ[cached] && rmax <= cached[[3]] && (method === "Series" || cached[[3]] <= rmax + 20), Return[cached]];
  bc = TeukolskyUpBC[s, \[Lambda], l, m, a, \[Omega], amps, \[Nu], Lookup[{ndsolveopts}, {WorkingPrecision, PrecisionGoal, AccuracyGoal}], rmax + 1, method];
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


TeukolskyInBC[s_Integer, \[Lambda]_, l_Integer, m_Integer, a_, \[Omega]_, amps_, \[Nu]_, {wp_, prec_, acc_}, rmin_:Automatic, method_:"MST"]:=
 Module[{R, res, dres, r, Rr, dRr},
        If[method === "Series",
          res = inBoundarySeries[s, \[Lambda], m, a, \[Omega]];
          If[res =!= $Failed, Return[res]]];
        R = Teukolsky`TeukolskyRadial[s, l, m, a, \[Omega], "BoundaryConditions" -> "In", "Amplitudes" -> amps, "Eigenvalue" -> \[Lambda], "RenormalizedAngularMomentum" -> \[Nu], Method -> "MST", WorkingPrecision -> wp, PrecisionGoal -> prec, AccuracyGoal -> Infinity];
        r = inBoundaryRadius[a, rmin];
		{Rr, dRr} = R[r, {0, 1}];
		TeukolskyInBCFromValues[s, m, a, \[Omega], r, Rr, dRr]
	];

(* HPS boundary data {psi, dpsi/dr, r} from the radial function value and derivative at r *)
TeukolskyInBCFromValues[s_Integer, m_Integer, a_, \[Omega]_, r_, Rr_, dRr_] :=
 Module[{res, dres},
		res = E^(I \[Omega] (r+((1+Sqrt[1-a^2]) Log[1/2 (-1-Sqrt[1-a^2]+r)]-(1-Sqrt[1-a^2]) Log[1/2 (-1+Sqrt[1-a^2]+r)])/Sqrt[1-a^2])) r ((-1-Sqrt[1-a^2]+r)/(-1+Sqrt[1-a^2]+r))^(-((I a m)/(2 Sqrt[1-a^2]))) (a^2-2 r+r^2)^s Rr;
		dres =E^(I \[Omega] (r+((1+Sqrt[1-a^2]) Log[1/2 (-1-Sqrt[1-a^2]+r)]-(1-Sqrt[1-a^2]) Log[1/2 (-1+Sqrt[1-a^2]+r)])/Sqrt[1-a^2])) ((-1-Sqrt[1-a^2]+r)/(-1+Sqrt[1-a^2]+r))^(-((I a m)/(2 Sqrt[1-a^2]))) (a^2-2 r+r^2)^s Rr-(I a E^(I \[Omega] (r+((1+Sqrt[1-a^2]) Log[1/2 (-1-Sqrt[1-a^2]+r)]-(1-Sqrt[1-a^2]) Log[1/2 (-1+Sqrt[1-a^2]+r)])/Sqrt[1-a^2])) m r ((-1-Sqrt[1-a^2]+r)/(-1+Sqrt[1-a^2]+r))^(-1-(I a m)/(2 Sqrt[1-a^2])) (a^2-2 r+r^2)^s (-((-1-Sqrt[1-a^2]+r)/(-1+Sqrt[1-a^2]+r)^2)+1/(-1+Sqrt[1-a^2]+r)) Rr)/(2 Sqrt[1-a^2])+E^(I \[Omega] (r+((1+Sqrt[1-a^2]) Log[1/2 (-1-Sqrt[1-a^2]+r)]-(1-Sqrt[1-a^2]) Log[1/2 (-1+Sqrt[1-a^2]+r)])/Sqrt[1-a^2])) r ((-1-Sqrt[1-a^2]+r)/(-1+Sqrt[1-a^2]+r))^(-((I a m)/(2 Sqrt[1-a^2]))) (-2+2 r) (a^2-2 r+r^2)^(-1+s) s Rr+I E^(I \[Omega] (r+((1+Sqrt[1-a^2]) Log[1/2 (-1-Sqrt[1-a^2]+r)]-(1-Sqrt[1-a^2]) Log[1/2 (-1+Sqrt[1-a^2]+r)])/Sqrt[1-a^2])) r ((-1-Sqrt[1-a^2]+r)/(-1+Sqrt[1-a^2]+r))^(-((I a m)/(2 Sqrt[1-a^2]))) (a^2+r^2) (a^2-2 r+r^2)^(-1+s) \[Omega] Rr+E^(I \[Omega] (r+((1+Sqrt[1-a^2]) Log[1/2 (-1-Sqrt[1-a^2]+r)]-(1-Sqrt[1-a^2]) Log[1/2 (-1+Sqrt[1-a^2]+r)])/Sqrt[1-a^2])) r ((-1-Sqrt[1-a^2]+r)/(-1+Sqrt[1-a^2]+r))^(-((I a m)/(2 Sqrt[1-a^2]))) (a^2-2 r+r^2)^s dRr;
				{res,dres,r}
	];

TeukolskyUpBC[s_Integer, \[Lambda]_, l_Integer, m_Integer, a_, \[Omega]_, amps_, \[Nu]_, {wp_, prec_, acc_}, rmax_:Automatic, method_:"MST"]:=
 Module[{R, res, dres, r, Rr, dRr},
		If[method === "Series",
		  res = upBoundarySeries[s, \[Lambda], m, a, \[Omega], If[NumericQ[rmax] && rmax < Infinity, rmax, rp[a, 1] + 1]];
		  If[res =!= $Failed, Return[res[[1 ;; 3]]]]];
		r = upBoundaryRadius[a, rmax];
        R = Teukolsky`TeukolskyRadial[s, l, m, a, \[Omega], "BoundaryConditions" -> "Up", "Amplitudes" -> amps, "Eigenvalue" -> \[Lambda], "RenormalizedAngularMomentum" -> \[Nu], Method -> "MST", WorkingPrecision -> wp, PrecisionGoal -> prec, AccuracyGoal -> Infinity];
		{Rr, dRr} = R[r, {0, 1}];
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
(*Boundary data from series*)


(* ::Text:: *)
(*At machine precision the boundary data come from series solutions of the integrator's equation for*)
(*psi = e^{-H i omega r*} r ((r - r+)/(r - r-))^{-i a m/(2 kappa)} Delta^s R rather than from precision-padded MST*)
(*evaluations: a convergent power series in r - r+ for "In" (radius of convergence r+ - r- = 2 kappa, evaluated at*)
(*r+ + 1/5, or at r+ + kappa for near-extremal spins), and the asymptotic series in 1/r for "Up", evaluated at the*)
(*smallest radius beyond the requested domain at which its optimal truncation reaches 1e-15 (about 18/omega for*)
(*low modes, more for high l). Both are normalised to unit transmission: psi_up -> 1 at infinity, and psi_in(r+)*)
(*is the limit of the conversion factor applied to Delta^-s e^{-i k r*}. The recurrences follow from the equation*)
(*P0 psi + P1 psi' + P2 psi'' = 0 with polynomial P_i, obtained once symbolically below.*)


(* {symbols, {P0, P1, P2}} as coefficient lists in r (degree 8), the integrator's equation times r^6 (a^2 + r^2) *)
hpsPolynomials = Block[{s = hpsS, \[Lambda] = hps\[Lambda], m = hpsM, a = hpsA, \[Omega] = hps\[Omega], H = hpsH, r = hpsR, y1 = hpsY1, y2 = hpsY2, y3 = hpsY3},   (* fixed private symbols *)
 Module[{ode, poly},
  ode = (((a^2-2 r+r^2) (2 a^2-2 r (1+s)-r^2 \[Lambda])-2 a (1+H) m r^2 (a^2+r^2) \[Omega]+2 I r^2 (-(1+H) (-a^2+r^2)+(1-H) r (a^2-2 r+r^2)) s \[Omega]+(1-H^2) r^2 (a^2+r^2)^2 \[Omega]^2-2 I a r (a^2-2 r+r^2) (m+a H \[Omega])) y1)/r^6+((2 (-a^2+r^2) (a^2-2 r+r^2))/(r^4 (a^2+r^2))-(2 (a^2-2 r+r^2) (a^2 (a^2-2 r+r^2)+(a^2+r^2) ((-1+r) r s-I r (a m+H (a^2+r^2) \[Omega]))))/(r^5 (a^2+r^2))) y2+((a^2-2 r+r^2)^2 y3)/r^4;
  poly = Cancel[Together[ode r^6 (a^2 + r^2)]];
  {{s, \[Lambda], m, a, \[Omega], H}, Table[PadRight[CoefficientList[Coefficient[poly, y], r], 9], {y, {y1, y2, y3}}]}
 ]];

hpsPolys[s_, \[Lambda]_, m_, a_, \[Omega]_, H_] := hpsPolynomials[[2]] /. Thread[hpsPolynomials[[1]] -> {s, \[Lambda], m, a, \[Omega], H}];

(* "In" boundary data {psi, psi', r} from the power series in x = r - r+, psi = Sum d_j x^j, at x = Min[1/5, kappa].
   The coefficient of x^j gives d_j from d_{j-8}..d_{j-1} (the regular singular point at r+ makes the
   equation for d_j homogeneous of degree j (j - 1 + q11/q22), where -q11/q22 is the exponent of the other
   solution); $Failed when that coefficient vanishes (the resonant case of the superradiant bound frequency
   with s >= 1, where the MST solution is used instead) or the series does not converge. *)
inBoundarySeries[s_, \[Lambda]_, m_, a_, \[Omega]_] :=
 Module[{\[Kappa], rp, x, xb, p, q, qq, d, d0, psi, dpsi, term, coef, rest, j, n},
  \[Kappa] = Sqrt[1 - a^2]; rp = 1 + \[Kappa]; xb = Min[1/5, \[Kappa]];
  p = hpsPolys[s, \[Lambda], m, a, \[Omega], -1];
  q = Table[PadRight[CoefficientList[Expand[Sum[pi[[k + 1]] (rp + x)^k, {k, 0, 8}]], x], 9], {pi, p}];
  qq[i_, k_] := If[0 <= k <= 8, q[[i + 1, k + 1]], 0];
  d0 = rp Exp[I a m/2] 2^(-I a m/(2 \[Kappa])) (2 \[Kappa])^(I a m/(2 \[Kappa])) \[Kappa]^(-I a m (1 - \[Kappa])/(2 \[Kappa] rp));
  d = {d0}; psi = d0; dpsi = 0;
  Do[
    coef = qq[0, 0] + j qq[1, 1] + j (j - 1) qq[2, 2];
    If[Abs[coef] < 10^-8 j^2 Abs[qq[2, 2]], Return[$Failed, Module]];
    rest = Sum[d[[n + 1]] (qq[0, j - n] + n qq[1, j - n + 1] + n (n - 1) qq[2, j - n + 2]), {n, Max[0, j - 8], j - 1}];
    AppendTo[d, -rest/coef];
    term = d[[-1]] xb^j; psi += term; dpsi += j d[[-1]] xb^(j - 1);
    If[j > 5 && Abs[term] < 10^-17 Abs[psi], Return[{psi, dpsi, rp + xb}, Module]],
    {j, 1, 800}];
  $Failed
 ];

(* "Up" boundary data {psi, psi', r} from the asymptotic series psi = Sum c_n r^-n (c_0 = 1), summed in the
   scaled form t_n = c_n r^-n and truncated at its smallest term, at the smallest radius r >= rmin at which that
   term is below 1e-15 relative (the radius is raised by factors of 5/4 until it is). The coefficient of
   r^(8-j) of the equation gives c_{j-1} from c_0..c_{j-2}. *)
upSeriesAt[s_, \[Lambda]_, m_, a_, \[Omega]_, r_] :=
 Module[{p, pp, t, psi, dpsi, best, coef, rest, tj, j, n},
  (* returns {psi, psi', relative size of the smallest term, the scaled coefficients t_n = c_n r^-n} *)
  p = hpsPolys[s, \[Lambda], m, a, \[Omega], 1];
  pp[i_, k_] := If[0 <= k <= 8, p[[i + 1, k + 1]], 0];
  t = {1.}; psi = 1.; dpsi = 0.; best = Infinity;
  Do[
    coef = pp[0, 7] - (j - 1) pp[1, 8];
    rest = Sum[t[[n + 1]] r^(n - (j - 1)) (pp[0, 8 - j + n] - n pp[1, 8 - j + n + 1] + n (n + 1) pp[2, 8 - j + n + 2]), {n, Max[0, j - 10], j - 2}];
    tj = -rest/coef;
    If[Abs[tj] > best, Break[]];
    AppendTo[t, tj]; best = Abs[tj];
    psi += tj; dpsi += -(j - 1) tj/r;
    If[Abs[tj] < 10^-17 Abs[psi], Break[]],
    {j, 2, 600}];
  (* a vanishing partial sum (1 - 4/r at r = 4 for a = 0, say) means the series has not converged there *)
  {psi, dpsi, If[Abs[psi] > 10^-6 Max[Abs[t]], best/Abs[psi], Infinity], t}
 ];

(* {psi, psi', r, t}: the series data at the smallest radius r >= rmin where the series converges to 1e-15 *)
upBoundarySeries[s_, \[Lambda]_, m_, a_, \[Omega]_, rmin_] :=
 Module[{r = N[Max[rmin, 2 rp[a, 1]]], res},
  While[r < 10.^7,
    res = upSeriesAt[s, \[Lambda], m, a, \[Omega], r];
    If[res[[3]] <= 10.^-15, Return[{res[[1]], res[[2]], r, res[[4]]}, Module]];
    r *= 5/4];
  $Failed
 ];

(* The "Up" solution on the full domain from the large-r series: the series itself beyond its radius rb, an
   interpolating function of a single inward integration between rc = r+ + 1 and rb, and, below rc, integration
   from rc on demand (the region next to the horizon, where the solution oscillates in r_* and the requested
   radius sets how far to go). All three pieces are in the integration variable psi. *)
upSeriesSolutionBuild[s_, \[Lambda]_, m_, a_, \[Omega]_, ndsolveopts___] :=
 Module[{rc = rp[a, 1] + 1, res, rb, t, \[Psi]b, d\[Psi]b, interp, lower, goals},
  res = upBoundarySeries[s, \[Lambda], m, a, \[Omega], rc + 1];
  If[res === $Failed, Return[$Failed]];
  {\[Psi]b, d\[Psi]b, rb, t} = res;
  goals = ndsolveOptions[ndsolveopts];
  interp = Integrator[s, \[Lambda], m, a, \[Omega], \[Psi]b, d\[Psi]b, rb, rc, rb, 1, goals];
  lower = AllIntegrator[s, \[Lambda], m, a, \[Omega], interp[rc], interp'[rc], rc, 1, goals];
  upSeriesSolution[t, rb, rc, interp, lower]
 ];

upSeriesSolution[t_, rb_, rc_, interp_, lower_][r_?NumericQ] :=
  Which[r >= rb, Sum[t[[n + 1]] (rb/r)^n, {n, 0, Length[t] - 1}], r >= rc, interp[r], True, lower[r]];
upSeriesSolution[t_, rb_, rc_, interp_, lower_][r:{__?NumericQ}] := Map[upSeriesSolution[t, rb, rc, interp, lower], r];
Derivative[k_Integer?Positive][upSeriesSolution[t_, rb_, rc_, interp_, lower_]][r_?NumericQ] :=
  Which[r >= rb, Sum[t[[n + 1]] (-1)^k Pochhammer[n, k] (rb/r)^n r^-k, {n, 0, Length[t] - 1}], r >= rc, Derivative[k][interp][r], True, Derivative[k][lower][r]];
Derivative[k_Integer?Positive][upSeriesSolution[t_, rb_, rc_, interp_, lower_]][r:{__?NumericQ}] := Map[Derivative[k][upSeriesSolution[t, rb, rc, interp, lower]], r];


(* ::Section::Closed:: *)
(*End Package*)


End[]
EndPackage[];
