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

psi[s_, \[Lambda]_, l_, m_, a_, \[Omega]_, "In", amps_, \[Nu]_, ndsolveopts___][rmax_?NumericQ] := psi[s, \[Lambda], l, m, a, \[Omega], "In", amps, \[Nu], ndsolveopts][{Automatic, rmax}];
psi[s_, \[Lambda]_, l_, m_, a_, \[Omega]_, "Up", amps_, \[Nu]_, ndsolveopts___][rmin_?NumericQ] := psi[s, \[Lambda], l, m, a, \[Omega], "Up", amps, \[Nu], ndsolveopts][{rmin, Automatic}];

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
    (* "In" of negative spin with series boundary data: the outward integration from the horizon loses digits at
       large radius when the reflection is small, so beyond the radius where the large-r series converge the
       solution is Binc R_ingoing + Bref R_up from the series and the MST amplitudes (see inSeriesSolutionBuild) *)
    If[bc === "In" && s < 0 && Lookup[{ndsolveopts}, "BoundaryMethod", "MST"] === "Series" && !(AssociationQ[bdata] && KeyExistsQ[bdata, "In"]),
      With[{sol = inSeriesSolutionBuild[s, \[Lambda], l, m, a, \[Omega], amps, \[Nu], ndsolveopts]}, If[sol =!= $Failed, Return[sol]]]];
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


(* NDSolve leaves behind temporary copies of its dependent and independent variables (niY1$123 and so on,
   in the Global` context when the variables were Global` symbols scoped by Module), one set per call; they are
   removed once the solution, which does not refer to them, has been computed. Their numbers lie between the
   values of $ModuleNumber before and after the call, so the candidates are checked with NameQ: searching with
   Names scans the whole symbol table and cost 20 ms per call, ten times the evaluation it followed. *)
SetAttributes[withoutTemporaries, HoldAll];
withoutTemporaries[expr_] := Module[{n0 = $ModuleNumber, res, ctx = Context[niY1], tmp},
  res = expr;
  (* a few hundred numbers per call; beyond 10^5 the check would cost more than the symbols, which are kept *)
  tmp = If[$ModuleNumber - n0 > 10^5, {},
    Select[Flatten[Table[ctx <> v <> "$" <> ToString[k], {v, {"niY1", "niY2", "niR"}}, {k, n0, $ModuleNumber}]], NameQ]];
  If[tmp =!= {}, Quiet[Remove @@ tmp]];
  res
 ];

Integrator[s_,\[Lambda]_,m_,a_,\[Omega]_,y1BC_,y2BC_,rBC_,rmin_?NumericQ,rmax_?NumericQ,H_?NumericQ,ndsolveopts___]:=Block[{niY1,niY2,niR,sol},
	withoutTemporaries@Quiet[NDSolveValue[
		{niY1'[niR]==niY2[niR],(((a^2-2 niR+niR^2) (2 a^2-2 niR (1+s)-niR^2 \[Lambda])-2 a (1+H) m niR^2 (a^2+niR^2) \[Omega]+2 I niR^2 (-(1+H) (-a^2+niR^2)+(1-H) niR (a^2-2 niR+niR^2)) s \[Omega]+(1-H^2) niR^2 (a^2+niR^2)^2 \[Omega]^2-2 I a niR (a^2-2 niR+niR^2) (m+a H \[Omega])) niY1[niR])/niR^6+((2 (-a^2+niR^2) (a^2-2 niR+niR^2))/(niR^4 (a^2+niR^2))-(2 (a^2-2 niR+niR^2) (a^2 (a^2-2 niR+niR^2)+(a^2+niR^2) ((-1+niR) niR s-I niR (a m+H (a^2+niR^2) \[Omega]))))/(niR^5 (a^2+niR^2))) niY2[niR]+((a^2-2 niR+niR^2)^2 Derivative[1][niY2][niR])/niR^4==0,niY1[rBC]==y1BC,niY2[rBC]==y2BC},
		niY1,
		{niR, rmin, rmax},
		ndsolveopts,
		Method->"StiffnessSwitching",
		MaxSteps->Infinity,
		InterpolationOrder->All
		], NDSolveValue::precw]
	];


AllIntegrator[s_,\[Lambda]_,m_,a_,\[Omega]_,y1BC_,y2BC_,rBC_,H_?NumericQ,ndsolveopts___][rval:(_?NumericQ | {_?NumericQ..})] := Block[{niY1,niY2,niR,sol},
	withoutTemporaries@Quiet[NDSolveValue[
		{niY1'[niR]==niY2[niR],(((a^2-2 niR+niR^2) (2 a^2-2 niR (1+s)-niR^2 \[Lambda])-2 a (1+H) m niR^2 (a^2+niR^2) \[Omega]+2 I niR^2 (-(1+H) (-a^2+niR^2)+(1-H) niR (a^2-2 niR+niR^2)) s \[Omega]+(1-H^2) niR^2 (a^2+niR^2)^2 \[Omega]^2-2 I a niR (a^2-2 niR+niR^2) (m+a H \[Omega])) niY1[niR])/niR^6+((2 (-a^2+niR^2) (a^2-2 niR+niR^2))/(niR^4 (a^2+niR^2))-(2 (a^2-2 niR+niR^2) (a^2 (a^2-2 niR+niR^2)+(a^2+niR^2) ((-1+niR) niR s-I niR (a m+H (a^2+niR^2) \[Omega]))))/(niR^5 (a^2+niR^2))) niY2[niR]+((a^2-2 niR+niR^2)^2 Derivative[1][niY2][niR])/niR^4==0,niY1[rBC]==y1BC,niY2[rBC]==y2BC},
		niY1[rval],
		{niR, Min[rBC,rval], Max[rBC,rval]},
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

Derivative[n_][AllIntegrator[s_,\[Lambda]_,m_,a_,\[Omega]_,y1BC_,y2BC_,rBC_,H_?NumericQ,ndsolveopts___]][rval:(_?NumericQ | {_?NumericQ..})] := Block[{niY1,niY2,niR,sol},
	withoutTemporaries@Quiet[NDSolveValue[
		{niY1'[niR]==niY2[niR],(((a^2-2 niR+niR^2) (2 a^2-2 niR (1+s)-niR^2 \[Lambda])-2 a (1+H) m niR^2 (a^2+niR^2) \[Omega]+2 I niR^2 (-(1+H) (-a^2+niR^2)+(1-H) niR (a^2-2 niR+niR^2)) s \[Omega]+(1-H^2) niR^2 (a^2+niR^2)^2 \[Omega]^2-2 I a niR (a^2-2 niR+niR^2) (m+a H \[Omega])) niY1[niR])/niR^6+((2 (-a^2+niR^2) (a^2-2 niR+niR^2))/(niR^4 (a^2+niR^2))-(2 (a^2-2 niR+niR^2) (a^2 (a^2-2 niR+niR^2)+(a^2+niR^2) ((-1+niR) niR s-I niR (a m+H (a^2+niR^2) \[Omega]))))/(niR^5 (a^2+niR^2))) niY2[niR]+((a^2-2 niR+niR^2)^2 Derivative[1][niY2][niR])/niR^4==0,niY1[rBC]==y1BC,niY2[rBC]==y2BC},
		Derivative[n][niY1][rval],
		{niR, Min[rBC,rval], Max[rBC,rval]},
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
        (* values of the parent solution, registered there if they were supplied *)
        R = Block[{Teukolsky`MST`MST`Private`$registerSupplied = False}, Teukolsky`TeukolskyRadial[s, l, m, a, \[Omega], "BoundaryConditions" -> "In", "Amplitudes" -> amps, "Eigenvalue" -> \[Lambda], "RenormalizedAngularMomentum" -> If[NumericQ[\[Nu]], \[Nu], Automatic], Method -> "MST", WorkingPrecision -> wp, PrecisionGoal -> prec, AccuracyGoal -> Infinity]];
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
        (* values of the parent solution, registered there if they were supplied *)
        R = Block[{Teukolsky`MST`MST`Private`$registerSupplied = False}, Teukolsky`TeukolskyRadial[s, l, m, a, \[Omega], "BoundaryConditions" -> "Up", "Amplitudes" -> amps, "Eigenvalue" -> \[Lambda], "RenormalizedAngularMomentum" -> If[NumericQ[\[Nu]], \[Nu], Automatic], Method -> "MST", WorkingPrecision -> wp, PrecisionGoal -> prec, AccuracyGoal -> Infinity]];
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
(* The precision the series are summed at: that of their inputs, with machine precision as $MachinePrecision
   (they used to be summed with machine-number literals whatever the inputs, which capped an arbitrary-precision
   integration started from them at about 15 digits) *)
seriesPrecision[\[Lambda]_, a_, \[Omega]_] := With[{p = Precision[{\[Lambda], a, \[Omega]}]}, If[p === MachinePrecision, $MachinePrecision, p]];
seriesNumber[x_, \[Lambda]_, a_, \[Omega]_] := If[Precision[{\[Lambda], a, \[Omega]}] === MachinePrecision, N[x], N[x, Precision[{\[Lambda], a, \[Omega]}]]];

inBoundarySeries[s_, \[Lambda]_, m_, a_, \[Omega]_] :=
 Module[{\[Kappa], rp, x, xb, p, q, qq, d, d0, psi, dpsi, term, coef, rest, j, n, res = $Failed, tol = 10^-(Floor[seriesPrecision[\[Lambda], a, \[Omega]]] + 2)},
  \[Kappa] = Sqrt[1 - a^2]; rp = 1 + \[Kappa]; xb = Min[1/5, \[Kappa]];
  p = hpsPolys[s, \[Lambda], m, a, \[Omega], -1];
  q = Table[PadRight[CoefficientList[Expand[Sum[pi[[k + 1]] (rp + x)^k, {k, 0, 8}]], x], 9], {pi, p}];
  qq[i_, k_] := If[0 <= k <= 8, q[[i + 1, k + 1]], 0];
  d0 = rp Exp[I a m/2] 2^(-I a m/(2 \[Kappa])) (2 \[Kappa])^(I a m/(2 \[Kappa])) \[Kappa]^(-I a m (1 - \[Kappa])/(2 \[Kappa] rp));
  d = {d0}; psi = d0; dpsi = 0;
  Do[
    coef = qq[0, 0] + j qq[1, 1] + j (j - 1) qq[2, 2];
    If[Abs[coef] < 10^-8 j^2 Abs[qq[2, 2]], Break[]];
    rest = Sum[d[[n + 1]] (qq[0, j - n] + n qq[1, j - n + 1] + n (n - 1) qq[2, j - n + 2]), {n, Max[0, j - 8], j - 1}];
    AppendTo[d, -rest/coef];
    term = d[[-1]] xb^j; psi += term; dpsi += j d[[-1]] xb^(j - 1);
    If[j > 5 && Abs[term] < tol Abs[psi], res = {psi, dpsi, rp + xb}; Break[]],
    {j, 1, 800}];
  (* a local function with definitions is not removed by Module when it is referenced from another local
     (here qq refers to q), and leaked on the first call *)
  Clear[qq];
  res
 ];

(* "Up" boundary data {psi, psi', r} from the asymptotic series psi = Sum c_n r^-n (c_0 = 1), summed in the
   scaled form t_n = c_n r^-n and truncated at its smallest term, at the smallest radius r >= rmin at which that
   term is below 1e-15 relative (the radius is raised by factors of 5/4 until it is). The coefficient of
   r^(8-j) of the equation gives c_{j-1} from c_0..c_{j-2}. *)
(* Optimal truncation of an asymptotic series given its scaled terms t_n (t_0 = 1): the sum up to the smallest
   term. The terms of these series need not decrease monotonically at first (the second term is small for m = 0
   or small a, where the first coefficient nearly or exactly vanishes), so an increase is only taken as the onset
   of the asymptotic divergence once it has persisted for three terms (seriesTerms); the list is then cut at the
   smallest term. Returns {psi, psi', relative size of the smallest term, t} for psi = r^rho Sum t_n. Stopping at
   the first increase, as before, cut the ingoing series after a small second term and pushed the join radius of
   the negative-spin "In" solution far out (1e-10 instead of 1e-14 at r = 100 for s = -2, m = 0, omega = 1). *)
seriesTerms[next_, tol_, jmax_] :=
 Module[{t = {1}, psi = 1, best = Infinity, jbest = 0, prev = Infinity, rise = 0, small = 0, tj, j},
  Do[
    tj = next[j, t]; AppendTo[t, tj];
    (* a coefficient that vanishes exactly (c_1 at a = 0, say) is neither the smallest term nor convergence *)
    If[Abs[tj] == 0, Continue[]];
    rise = If[Abs[tj] > prev, rise + 1, 0]; prev = Abs[tj];
    If[Abs[tj] < best, best = Abs[tj]; jbest = j - 1];
    If[rise >= 3, Break[]];
    psi += tj;
    small = If[Abs[tj] < tol Abs[psi], small + 1, 0];
    If[small >= 2, jbest = j - 1; Break[]],
    {j, 2, jmax}];
  {Take[t, jbest + 1], best}
 ];

seriesResult[t_, best_, \[Rho]_, r_] :=
 Module[{psi = Total[t], dpsi = Sum[(\[Rho] - n) t[[n + 1]], {n, 0, Length[t] - 1}]/r},
  (* a vanishing partial sum (1 - 4/r at r = 4 for a = 0, say) means the series has not converged there *)
  {psi r^\[Rho], dpsi r^\[Rho], If[Abs[psi] > 10^-6 Max[Abs[t]], best/Abs[psi], Infinity], t}
 ];

upSeriesAt[s_, \[Lambda]_, m_, a_, \[Omega]_, r_] :=
 Module[{p, pp, t, best, res, tol = 10^-(Floor[seriesPrecision[\[Lambda], a, \[Omega]]] + 2)},
  (* returns {psi, psi', relative size of the smallest term, the scaled coefficients t_n = c_n r^-n}; the
     arithmetic is that of the inputs (exact literals) *)
  p = hpsPolys[s, \[Lambda], m, a, \[Omega], 1];
  pp[i_, k_] := If[0 <= k <= 8, p[[i + 1, k + 1]], 0];
  {t, best} = seriesTerms[
    Function[{j, tt}, -Sum[tt[[n + 1]] r^(n - (j - 1)) (pp[0, 8 - j + n] - n pp[1, 8 - j + n + 1] + n (n + 1) pp[2, 8 - j + n + 2]), {n, Max[0, j - 10], j - 2}]/(pp[0, 7] - (j - 1) pp[1, 8])],
    tol, 600];
  res = seriesResult[t, best, 0, r];
  Clear[pp];   (* see inBoundarySeries *)
  res
 ];

(* {psi, psi', r, t}: the series data at the smallest radius r >= rmin where the series converges to the
   precision of the inputs (1e-15 at machine precision) *)
upBoundarySeries[s_, \[Lambda]_, m_, a_, \[Omega]_, rmin_] :=
 Module[{r = seriesNumber[Max[rmin, 2 rp[a, 1]], \[Lambda], a, \[Omega]], res, tol = 10^-Floor[seriesPrecision[\[Lambda], a, \[Omega]]]},
  While[r < 10.^7,
    res = upSeriesAt[s, \[Lambda], m, a, \[Omega], r];
    If[res[[3]] <= tol, Return[{res[[1]], res[[2]], r, res[[4]]}, Module]];
    r *= 5/4];
  $Failed
 ];

(* The ingoing asymptotic series in the "In" variable (H = -1): psi = r^(2s) Sum e_n r^-n with e_0 = 1, which is
   R -> e^{-i omega r*}/r, summed like the "Up" series; the coefficient of r^(2s+8-j) gives e_{j-1}, the two
   top powers vanish identically for this exponent. Returns {psi, psi', relative size of the smallest term, t}. *)
inSeriesAt[s_, \[Lambda]_, m_, a_, \[Omega]_, r_] :=
 Module[{p, pp, A, \[Rho] = 2 s, t, best, res, tol = 10^-(Floor[seriesPrecision[\[Lambda], a, \[Omega]]] + 2)},
  p = hpsPolys[s, \[Lambda], m, a, \[Omega], -1];
  pp[i_, k_] := If[0 <= k <= 8, p[[i + 1, k + 1]], 0];
  A[j_, n_] := pp[0, 8 - j + n] + (\[Rho] - n) pp[1, 8 - j + n + 1] + (\[Rho] - n) (\[Rho] - n - 1) pp[2, 8 - j + n + 2];
  {t, best} = seriesTerms[Function[{j, tt}, -Sum[tt[[n + 1]] r^(n - (j - 1)) A[j, n], {n, Max[0, j - 10], j - 2}]/A[j, j - 1]], tol, 600];
  res = seriesResult[t, best, \[Rho], r];
  Clear[pp, A];   (* see inBoundarySeries *)
  res
 ];

tortoise[r_, a_] := With[{\[Kappa] = Sqrt[1 - a^2]}, r + ((1 + \[Kappa]) Log[(r - 1 - \[Kappa])/2] - (1 - \[Kappa]) Log[(r - 1 + \[Kappa])/2])/\[Kappa]];

(* "In" data {psi, psi', r} in the H = -1 variable at the smallest radius r >= rmin where both large-r series
   converge to 1e-15: psi_in = Binc psi_ingoing + Bref e^{2 i omega r*} psi_up, with the amplitudes relative to
   unit transmission. Also returns the two scaled coefficient lists. *)
inLargeRadiusData[s_, \[Lambda]_, m_, a_, \[Omega]_, Binc_, Bref_, rmin_] :=
 Module[{r = seriesNumber[Max[rmin, 2 rp[a, 1]], \[Lambda], a, \[Omega]], ing, up, ph, dph, psi, dpsi, tol = 10^-Floor[seriesPrecision[\[Lambda], a, \[Omega]]]},
  While[r < 10.^7,
    ing = inSeriesAt[s, \[Lambda], m, a, \[Omega], r]; up = upSeriesAt[s, \[Lambda], m, a, \[Omega], r];
    If[ing[[3]] <= tol && up[[3]] <= tol,
      ph = Exp[2 I \[Omega] tortoise[r, a]]; dph = 2 I \[Omega] (r^2 + a^2)/(r^2 - 2 r + a^2) ph;
      psi = Binc ing[[1]] + Bref ph up[[1]]; dpsi = Binc ing[[2]] + Bref (dph up[[1]] + ph up[[2]]);
      Return[{psi, dpsi, r, ing[[4]], up[[4]]}, Module]];
    r *= 5/4];
  $Failed
 ];

(* The "In" solution of negative spin on the full domain: the outward integration from the horizon series below
   the radius rb where the large-r series converge (on demand, as before), and Binc R_ingoing + Bref R_up from
   the series beyond it. Both directions of integration are unstable somewhere for this solution: outwards,
   once omega r >> 1, roundoff in the outgoing solution grows like r^(-2s) relative to the incident part; inwards,
   the solution singular at the horizon grows like a high power of 1/r wherever omega r << 1, so the far-zone
   series cannot be integrated in to the horizon either, and the join is made at rb. The amplitudes are needed;
   without them ("Amplitudes" -> False) $Failed, and the outward integration is used throughout. *)
inSeriesSolutionBuild[s_, \[Lambda]_, l_, m_, a_, \[Omega]_, amps_, \[Nu]_, ndsolveopts___] :=
 Module[{Binc, Bref, res, rb, tIng, tUp, lower, psiBC, dpsidrBC, rBC},
  (* the representation is relative to unit transmission; where the transmission vanishes (a degeneracy
     2 I epsilon_+ = n >= 1 - s, at which the "In" solution is normalised to unit incidence instead) it is not
     available and the outward integration is used throughout *)
  If[!(AssociationQ[amps] && AssociationQ[amps["In"]] && AllTrue[Lookup[amps["In"], {"Incidence", "Reflection", "Transmission"}, Indeterminate], NumericQ] && amps["In"]["Transmission"] != 0), Return[$Failed]];
  Binc = amps["In"]["Incidence"]/amps["In"]["Transmission"]; Bref = amps["In"]["Reflection"]/amps["In"]["Transmission"];
  res = inLargeRadiusData[s, \[Lambda], m, a, \[Omega], Binc, Bref, rp[a, 1] + 2];
  If[res === $Failed, Return[$Failed]];
  rb = res[[3]]; {tIng, tUp} = res[[4 ;; 5]];
  {psiBC, dpsidrBC, rBC} = boundaryData[s, \[Lambda], l, m, a, \[Omega], "In", amps, \[Nu], {Automatic, Automatic}, ndsolveopts];
  lower = AllIntegrator[s, \[Lambda], m, a, \[Omega], psiBC, dpsidrBC, rBC, -1, ndsolveOptions[ndsolveopts]];
  inSeriesSolution[s, a, \[Omega], Binc, Bref, tIng, tUp, rb, lower]
 ];

inSeriesSolution[s_, a_, \[Omega]_, Binc_, Bref_, tIng_, tUp_, rb_, lower_][r_?NumericQ] :=
  If[r >= rb,
    Binc r^(2 s) Sum[tIng[[n + 1]] (rb/r)^n, {n, 0, Length[tIng] - 1}] + Bref Exp[2 I \[Omega] tortoise[r, a]] Sum[tUp[[n + 1]] (rb/r)^n, {n, 0, Length[tUp] - 1}],
    lower[r]];
(* lists: the radii below rb go to the integrator in one call (one integration), the others to the series *)
inSeriesSolution[s_, a_, \[Omega]_, Binc_, Bref_, tIng_, tUp_, rb_, lower_][r:{__?NumericQ}] :=
  splitEvaluate[inSeriesSolution[s, a, \[Omega], Binc, Bref, tIng, tUp, rb, lower], lower, rb, r];
Derivative[1][inSeriesSolution[s_, a_, \[Omega]_, Binc_, Bref_, tIng_, tUp_, rb_, lower_]][r_?NumericQ] :=
  If[r >= rb,
    With[{ing = Sum[tIng[[n + 1]] (rb/r)^n, {n, 0, Length[tIng] - 1}], ding = Sum[-n tIng[[n + 1]] (rb/r)^n/r, {n, 0, Length[tIng] - 1}],
          up = Sum[tUp[[n + 1]] (rb/r)^n, {n, 0, Length[tUp] - 1}], dup = Sum[-n tUp[[n + 1]] (rb/r)^n/r, {n, 0, Length[tUp] - 1}], ph = Exp[2 I \[Omega] tortoise[r, a]]},
      Binc (2 s r^(2 s - 1) ing + r^(2 s) ding) + Bref ph (2 I \[Omega] (r^2 + a^2)/(r^2 - 2 r + a^2) up + dup)],
    lower'[r]];
Derivative[1][inSeriesSolution[s_, a_, \[Omega]_, Binc_, Bref_, tIng_, tUp_, rb_, lower_]][r:{__?NumericQ}] :=
  splitEvaluate[Derivative[1][inSeriesSolution[s, a, \[Omega], Binc, Bref, tIng, tUp, rb, lower]], Derivative[1][lower], rb, r];

(* f on a list of radii: those below rsplit evaluated by g in a single call, the others by f one at a time *)
splitEvaluate[f_, g_, rsplit_, r_List] :=
 Module[{low = Flatten[Position[r, x_ /; x < rsplit, {1}]], res = ConstantArray[0, Length[r]]},
  If[low =!= {}, res[[low]] = g[r[[low]]]];
  res[[Complement[Range[Length[r]], low]]] = f /@ r[[Complement[Range[Length[r]], low]]];
  res
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
upSeriesSolution[t_, rb_, rc_, interp_, lower_][r:{__?NumericQ}] := splitEvaluate[upSeriesSolution[t, rb, rc, interp, lower], lower, rc, r];
Derivative[k_Integer?Positive][upSeriesSolution[t_, rb_, rc_, interp_, lower_]][r_?NumericQ] :=
  Which[r >= rb, Sum[t[[n + 1]] (-1)^k Pochhammer[n, k] (rb/r)^n r^-k, {n, 0, Length[t] - 1}], r >= rc, Derivative[k][interp][r], True, Derivative[k][lower][r]];
Derivative[k_Integer?Positive][upSeriesSolution[t_, rb_, rc_, interp_, lower_]][r:{__?NumericQ}] := splitEvaluate[Derivative[k][upSeriesSolution[t, rb, rc, interp, lower]], Derivative[k][lower], rc, r];


(* ::Section::Closed:: *)
(*End Package*)


End[]
EndPackage[];
