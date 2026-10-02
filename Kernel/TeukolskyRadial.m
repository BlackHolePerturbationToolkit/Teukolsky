(* ::Package:: *)

(* ::Title:: *)
(*TeukolskyRadial*)


(* ::Section::Closed:: *)
(*Create Package*)


(* ::Subsection::Closed:: *)
(*BeginPackage*)


BeginPackage["Teukolsky`TeukolskyRadial`",
  {
  "Teukolsky`",
  "Teukolsky`SasakiNakamura`",
   "Teukolsky`NumericalIntegration`",
   "Teukolsky`MST`RenormalizedAngularMomentum`",
   "Teukolsky`MST`MST`",
   "SpinWeightedSpheroidalHarmonics`"
  }
];


(* ::Subsection::Closed:: *)
(*Unprotect symbols*)


ClearAttributes[{TeukolskyRadial, TeukolskyRadialFunction}, {Protected, ReadProtected}];


(* ::Subsection::Closed:: *)
(*Usage messages*)


TeukolskyRadial::usage = "TeukolskyRadial[s, l, m, a, \[Omega]] computes homogeneous solutions to the radial Teukolsky equation."
TeukolskyRadialFunction::usage = "TeukolskyRadialFunction[s, l, m, a, \[Omega], assoc] is an object representing a homogeneous solution to the radial Teukolsky equation."


(* ::Subsection::Closed:: *)
(*Error Messages*)


TeukolskyRadial::precw = "The precision of `1`=`2` is less than WorkingPrecision (`3`).";
TeukolskyRadial::optx = "Unknown options in `1`";
TeukolskyRadial::params = "Invalid parameters s=`1`, l=`2`, m=`3`";
TeukolskyRadial::cmplx = "Only real values of a are allowed, but a=`1` specified.";
TeukolskyRadial::spin = "The spin a=`1` must satisfy |a| < 1.";
TeukolskyRadial::dm = "Option `1` is not valid with BoundaryConditions \[RightArrow] `2`.";
TeukolskyRadial::sopt = "Option `1` not supported for static (\[Omega]=0) modes.";
TeukolskyRadial::hc = "Method HeunC is only supported with Mathematica version 12.1 and later.";
TeukolskyRadial::hcopt = "Option `1` not supported for HeunC method.";
TeukolskyRadialFunction::dmval = "Radius `1` lies outside the computational domain.";
TeukolskyRadial::opti = "Options in set `1` are incompatible.";
TeukolskyRadial::topopt = "`1` are options of TeukolskyRadial, not of Method `2`; specify them outside Method.";
TeukolskyRadial::exact = "Exact arguments a=`1`, \[Omega]=`2` require a WorkingPrecision; specify one, or apply N or SetPrecision to the arguments.";
TeukolskyRadial::superradiant = "\[Omega] = m \[CapitalOmega]_H is the superradiant bound frequency: the horizon basis is resonant and the amplitudes `1`, which diverge there, are Indeterminate.";
TeukolskyRadial::degenerate = "2 I \[Epsilon]_+ = `1` is an integer at \[Omega] = `2`, where the horizon basis is resonant: the amplitudes `3`, which diverge there, are Indeterminate.";
TeukolskyRadial::innorm = "The transmission amplitude of the \"In\" solution vanishes at \[Omega] = `1` (2 I \[Epsilon]_+ = `2` >= 1 - s); it is normalised to unit incidence instead of unit transmission.";


(* ::Subsection::Closed:: *)
(*Begin Private section*)


Begin["`Private`"];


(* ::Section::Closed:: *)
(*Utility Functions*)


(* ::Subsection::Closed:: *)
(*Horizon Locations*)


rp[a_,M_] := M+Sqrt[M^2-a^2];

rm[a_,M_] := M-Sqrt[M^2-a^2];

(* Evaluate f[p] (a computation with its inputs set to precision p) with p padded until the result has the
   precision of the input: the MST quantities (renormalized angular momentum, asymptotic amplitudes, radial
   functions) lose digits to cancellation in their series, which is recovered by working at a higher
   precision. A non-numeric result (a precision-zero intermediate) is retried at twice the precision. For
   machine-precision input the computation is done in arbitrary precision. A result still short of the
   input precision after the retries is reported. *)
TeukolskyRadial::prec = "`1` could only be computed to a precision of `2` (`3` requested).";
TeukolskyRadial::nufail = "The renormalized angular momentum could not be computed for s=`1`, l=`2`, m=`3`, a=`4`, \[Omega]=`5`.";
TeukolskyRadial::ampfail = "The asymptotic amplitudes `2` could not be computed at \[Omega] = `1`: the MST amplitude formulae overflow or underflow at large |\[Omega]|. Radial functions from numerical integration keep the unit-transmission normalisation of their boundary data; the amplitudes, and MST radial functions whose normalising amplitude is missing, are Indeterminate.";

$lastPaddingDigits = 0;        (* digits of padding the last paddedComputation added beyond its first pass *)
$lastPaddingRetried = False;   (* whether it had to retry after a non-numeric result *)

paddedComputation[f_, wp_, name_:"The result", extra_:0] :=
 Module[{target, p, p0, res, deficit, tries = 0, minprec, grpPrec, numericQ, retried = False},
  (* The precision of the result: the smallest precision of its numbers. An inexact zero has no precision of
     its own; it counts with the digits its accuracy gives relative to the largest number of the same group
     (the amplitudes of one solution), and a group of nothing but zeros has none. The amplitudes at omega = 50
     cancel to 0``-136 at 32 digits, which was taken as an exact 0 (and reported as an overflow) although
     padding gives them to full precision. Exact zeros are exact. *)
  minprec[x_] := Module[{grp},
    grp[y_] := Which[
      AssociationQ[y] || ListQ[y], With[{v = If[AssociationQ[y], Values[y], y]},
        If[NoneTrue[v, AssociationQ[#] || ListQ[#] &], grpPrec[Select[v, NumericQ]], Min[Infinity, grp /@ v]]],
      NumericQ[y], grpPrec[{y}],
      True, Infinity];
    grp[x]];
  grpPrec[v_List] := Module[{inexact = Select[v, NumericQ[#] && Precision[#] < Infinity &], nz, scale},
    nz = Select[inexact, # != 0 &];
    scale = If[nz === {}, 0, Max[Abs[N[nz]]]];
    Min[Infinity, Precision /@ nz, If[scale > 0, Accuracy[#] + Log10[scale], 0] & /@ Select[inexact, # == 0 &]]];
  numericQ[x_] := Module[{nums},
    nums[y_] := If[AssociationQ[y], Flatten[nums /@ Values[y]], If[ListQ[y], Flatten[nums /@ y], {y}]];
    AllTrue[nums[x], NumericQ[#] || ($acceptIndeterminate && # === Indeterminate) &]];
  If[wp === MachinePrecision,
    target = $MachinePrecision; p = 2 $MachinePrecision;,
    target = wp; p = wp;
  ];
  p = p0 = Ceiling[p + extra];   (* an integer working precision; wp may be a real from Precision[...] *)
  (* a precision-zero intermediate at too low a working precision is retried, not reported; likewise a
     failure of the monodromy method at a precision the retries will raise (its convergence message is
     for a final failure, which the caller reports) *)
  res = Quiet[f[p], {Power::infy, Infinity::indet, Divide::infy, RenormalizedAngularMomentum::conv}];
  While[tries < 4 && (!numericQ[res] || (deficit = target - minprec[res]) > 1),
    If[!numericQ[res], retried = True];
    p = If[numericQ[res], p + Ceiling[deficit] + 3, 2 p];
    res = Quiet[f[p], {Power::infy, Infinity::indet, Divide::infy, RenormalizedAngularMomentum::conv}];
    tries++;
  ];
  $lastPaddingDigits = p - p0;
  $lastPaddingRetried = retried;
  If[!numericQ[res] || minprec[res] < target - 1,
    Message[TeukolskyRadial::prec, name, If[numericQ[res], minprec[res], res], target];
  ];
  (* a result that reached the target is rounded down to it; one that fell short keeps its actual precision,
     so that the shortfall just reported stays visible downstream rather than being dressed up as the target *)
  Which[
    wp === MachinePrecision, N[res],
    (* zeros keep their accuracy: SetPrecision would make them exact *)
    numericQ[res] && minprec[res] >= target - 1, With[{tg = target}, res /. z_?InexactNumberQ :> If[z == 0, z, SetPrecision[z, tg]]],
    True, res]
];


(* Relative error of the Wronskian of an "In"/"Up" pair of MST solutions against its value from the
   asymptotic amplitudes, Delta^(s+1) (R_in R_up' - R_in' R_up) = 2 I omega B^inc C^trans, at one radius. An
   independent check of the whole MST construction: for some modes (e.g. l = 36, m = 2, omega = 3 at
   a = 3/5) the recurrence for the MST coefficients silently yields the wrong solution below a certain
   working precision, with nothing in the tracked precision to show it, and this check catches it. The
   check is skipped (0 returned) when the amplitudes are not numeric, e.g. at the superradiant bound
   frequency for s >= 1. *)
mstWronskianError[R_Association, s_Integer, a_, \[Omega]_, wp_] :=
 Module[{r, W, scale, Wexact},
  r = 2 rp[a, 1];
  Wexact = 2 I \[Omega] R["In"]["Amplitudes"]["Incidence"] R["Up"]["Amplitudes"]["Transmission"];
  If[!NumericQ[Wexact] || Wexact == 0, Return[0]];
  {W, scale} = Quiet[Module[{i, di, u, du}, {i, di} = R["In"][r, {0, 1}]; {u, du} = R["Up"][r, {0, 1}]; (r^2 - 2 r + a^2)^(s + 1) {i du - di u, Abs[i du] + Abs[di u]}]];
  If[!NumericQ[W], Return[Infinity]];
  (* relative to the larger of the exact Wronskian and the size of its two terms: at a complex frequency
     the terms are of order Exp[2 |Im omega| r*] times W and the identity can only hold to that
     cancellation, which is not an error of the solutions *)
  Abs[W - Wexact]/Max[Abs[Wexact], scale]
 ];


(* The unscaled MST amplitudes at working precision p, with the eigenvalue and the renormalized angular
   momentum refined to p digits (the series amplify an error in nu by roughly the number of digits they
   lose to cancellation, so nu must carry the padded precision, not merely be set to it). The refined
   values are recorded in $refinedEigenvalue and $refinedNu. *)
(* At a degeneracy 2 I epsilon_+ = n (see epsilonPlusDegeneracy) the frequency is replaced, at the padded
   precision p, by the one at which 2 I epsilon_+ = n exactly, (m a - I n kappa)/(2 (1 + kappa)): the given
   omega satisfies the condition only to its own precision, and raised to p it would miss the MST package's
   snapping tolerance 10^(3 - p), leaving the amplitude factors at a near-pole instead of their limits
   (a machine-precision omega = m Omega_H, padded to 32 digits, is off by 1e-17). *)
degenerateFrequency[m_, a_, n_] := With[{\[Kappa] = Sqrt[1 - a^2]}, (m a - I n \[Kappa])/(2 (1 + \[Kappa]))];
$degeneracyOrder = None;

mstAmplitudes[s_, l_, m_, a_, \[Omega]_, \[Lambda]_, \[Nu]_, p_, prec_, acc_] :=
 Module[{\[Lambda]p, \[Nu]p, ap = SetPrecision[a, p], \[Omega]p},
  \[Omega]p = If[$degeneracyOrder === None, SetPrecision[\[Omega], p], degenerateFrequency[m, ap, $degeneracyOrder]];
  {\[Lambda]p, \[Nu]p} = Teukolsky`MST`MST`Private`refinedParameters[s, l, m, ap, 2 \[Omega]p, \[Lambda], \[Nu], p];
  {$refinedEigenvalue, $refinedNu} = {\[Lambda]p, \[Nu]p};
  Teukolsky`MST`MST`Private`Amplitudes[s, l, m, ap, 2 \[Omega]p, \[Nu]p, \[Lambda]p, {p, Max[prec, p - 2], acc}]
 ];


(* Degeneracies of the MST formulae: 2 I epsilon_+ = n, an integer, where epsilon_+ = (epsilon + tau)/2 =
   2 r+ k/(r+ - r-) with k = omega - m Omega_H. For real omega only n = 0 occurs, the
   superradiant bound frequency omega = m Omega_H, where the horizon basis Exp[i k r_*], Delta^-s Exp[-i k r_*]
   is resonant (coincident for s = 0, Frobenius exponents differing by the integer s otherwise); on the
   imaginary axis every n does. The MST package evaluates the amplitudes there directly, from their
   pole-free form (see Amplitudes): the "In" amplitudes are finite, with a transmission amplitude that
   vanishes for n >= 1 - s, where the regularised "In" solution is the horizon solution of larger exponent
   and is normalised to unit incidence (see normalisationKey); of the "Up" horizon coefficients the one along
   the basis function that has ceased to be independent diverges and is Indeterminate (the "Reflection",
   along Delta^-s Exp[-i k r_*], for n < 1 - s, the "Incidence", along Exp[i k r_*], for n > -1 - s, both for
   -1 - s < n < 1 - s, i.e. at the bound frequency for s = 0). Returns n, or None. Frequencies within
   10^(3 - wp) of a degeneracy in 2 I epsilon_+ are treated as being at it, the tolerance the MST package
   snaps with. *)
epsilonPlusDegeneracy[s_, m_, a_, \[Omega]_, wp_] :=
 Module[{\[Kappa] = Sqrt[1 - a^2], \[Epsilon] = 2 \[Omega], \[Tau], x, n, tol},
  \[Tau] = (\[Epsilon] - m a)/\[Kappa];
  x = I (\[Tau] + \[Epsilon]);   (* 2 I epsilon_+ *)
  n = Round[Re[x]];
  tol = 10^(3 - Floor[If[wp === MachinePrecision, $MachinePrecision, wp]]);
  If[Abs[x - n] <= tol, n, None]
 ];

(* Messages at a degeneracy: which amplitudes are Indeterminate, and the normalisation of the "In" solution *)
degeneracyMessages[n_, \[Omega]_, amps_Association] :=
 Module[{divergent},
  divergent = Flatten[Table[If[amps[bc][key] === Indeterminate, {bc, key}, Nothing], {bc, {"In", "Up"}}, {key, {"Incidence", "Transmission", "Reflection"}}], 1];
  If[n === 0, Message[TeukolskyRadial::superradiant, divergent], Message[TeukolskyRadial::degenerate, n, \[Omega], divergent]];
  If[Block[{$degenerateIn = True}, normalisationKey[amps["In"], "In"]] === "Incidence", Message[TeukolskyRadial::innorm, \[Omega], n]];
 ];

(* Non-numeric entries accepted by paddedComputation: at a degeneracy the divergent amplitudes are
   Indeterminate by design and must not trigger its retries *)
$acceptIndeterminate = False;

(* The amplitude a solution is normalised to: unit transmission, except for an "In" solution whose
   transmission amplitude vanishes at a degeneracy 2 I epsilon_+ = n >= 1 - s (see epsilonPlusDegeneracy),
   which is normalised to unit incidence. $degenerateIn, set by TeukolskyRadial while it builds the solutions
   of a mode at such a degeneracy, restricts that to the case it is meant for: a zero transmission of the "Up"
   solution, or of the "In" solution away from a degeneracy, is an overflow or underflow of the amplitude
   formulae (TeukolskyRadial::ampfail), not a normalisation convention. *)
$degenerateIn = False;
(* The eigenvalue and nu as stored in a radial function: at machine precision the values refined to the padded
   precision of the MST evaluations are kept internally, and the properties report machine numbers, as the
   working precision asks for *)
reportedValue[x_, MachinePrecision] := If[NumericQ[x], N[x], x];
reportedValue[x_, _] := x;

normalisationKey[ns_Association, bc_] := If[bc === "In" && $degenerateIn && NumericQ[ns["Transmission"]] && ns["Transmission"] == 0, "Incidence", "Transmission"];

(* The amplitudes relative to the normalising one; all Indeterminate when that one is not a nonzero number (the
   formulae overflowed or underflowed, see TeukolskyRadial::ampfail), rather than a division by zero *)
normaliseAmplitudes[ns_Association, bc_] :=
 With[{k = normalisationKey[ns, bc]}, If[NumericQ[ns[k]] && ns[k] != 0, ns/ns[k], Indeterminate & /@ ns]];



(* ::Subsection::Closed:: *)
(*Hyperboloidal Transformation Functions*)


f[r_] := 1-2/r;
rs[r_,a_]:=r+2/(rp[a,1]-rm[a,1]) (rp[a,1] Log[(r-rp[a,1])/2]-rm[a,1] Log[(r-rm[a,1])/2]);
\[CapitalDelta][r_,a_]:=r^2+a^2-2r;
\[Phi]Reg[r_,a_]:=a /(rp[a,1]-rm[a,1]) Log[(r-rp[a,1])/(r-rm[a,1])];


(* ::Section::Closed:: *)
(*TeukolskyRadial*)


(* ::Subsection::Closed:: *)
(*Numerical Integration Method*)


Options[TeukolskyRadialNumericalIntegration] = Join[
  {"Domain" -> All, "BoundaryData" -> None, "BoundaryMethod" -> Automatic},
  FilterRules[Options[NDSolve], Except[WorkingPrecision|AccuracyGoal|PrecisionGoal]]];


domainQ[domain_] := MatchQ[domain, {_?NumericQ, _?NumericQ} | (_?NumericQ) | All];


(* The Teukolsky-Starobinsky map from spin +s (s > 0) to spin -s, Delta^s (D0^dagger)^(2s) Delta^s with
   D0^dagger = d/dr + I K/Delta, K = (r^2 + a^2) omega - a m, reduced with the radial equation of spin +s and
   eigenvalue lambda (of that spin) to {f, g} with R_{-s} proportional to f R_{+s} + g R_{+s}'. For the "Up"
   solutions normalised to unit transmission the constant of proportionality is (2 I omega)^(2s). Derived
   symbolically once per spin. *)
teukolskyStarobinskyFlipSymbolic[s_Integer?Positive] := teukolskyStarobinskyFlipSymbolic[s] =
 Block[{tsR, tsr, ts\[Lambda], tsa, tsm, ts\[Omega]},   (* fixed private symbols, so that the memoised result holds no Module-generated ones *)
 Module[{K, \[CapitalDelta], Ddag, expr, rules},
  K = (tsr^2 + tsa^2) ts\[Omega] - tsa tsm; \[CapitalDelta] = tsr^2 - 2 tsr + tsa^2;
  Ddag[e_] := D[e, tsr] + I K/\[CapitalDelta] e;
  expr = \[CapitalDelta]^s Nest[Ddag, \[CapitalDelta]^s tsR[tsr], 2 s];
  rules = {Derivative[n_][tsR][tsr] :> D[(-(-ts\[Lambda] + 2 I tsr s 2 ts\[Omega] + (-2 I (-1 + tsr) s (-tsa tsm + (tsa^2 + tsr^2) ts\[Omega]) + (-tsa tsm + (tsa^2 + tsr^2) ts\[Omega])^2)/(tsa^2 - 2 tsr + tsr^2)) tsR[tsr] - (-2 + 2 tsr) (1 + s) tsR'[tsr])/(tsa^2 - 2 tsr + tsr^2), {tsr, n - 2}] /; n >= 2};
  expr = Collect[expr //. rules, {tsR[tsr], tsR'[tsr]}, Together];
  {{ts\[Lambda], tsa, tsm, ts\[Omega], tsr}, {Coefficient[expr, tsR[tsr]], Coefficient[expr, tsR'[tsr]]}}
 ]];

(* The Teukolsky-Starobinsky map (f R+ + g R+')/C of a spin s+ solution R+ (given as a pure function), with the
   derivative (f' R+ + (f + g') R+' + g R+'')/C in which R+'' comes from the radial equation of spin s+ and
   eigenvalue lambda+ *)
tsMappedFunction[f_, g_, df_, dg_, C_, Rp_, sp_, \[Lambda]p_, a_, m_, \[Omega]_][r_?NumericQ] := (f[r] Rp[r] + g[r] Rp'[r])/C;
tsMappedFunction[f_, g_, df_, dg_, C_, Rp_, sp_, \[Lambda]p_, a_, m_, \[Omega]_][r:{__?NumericQ}] := Map[tsMappedFunction[f, g, df, dg, C, Rp, sp, \[Lambda]p, a, m, \[Omega]], r];
Derivative[1][tsMappedFunction[f_, g_, df_, dg_, C_, Rp_, sp_, \[Lambda]p_, a_, m_, \[Omega]_]][r_?NumericQ] :=
 With[{R = Rp[r], dR = Rp'[r]},
  (df[r] R + (f[r] + dg[r]) dR + g[r] (-(-\[Lambda]p + 2 I r sp 2 \[Omega] + (-2 I (-1 + r) sp (-a m + (a^2 + r^2) \[Omega]) + (-a m + (a^2 + r^2) \[Omega])^2)/(a^2 - 2 r + r^2)) R - (-2 + 2 r) (1 + sp) dR)/(a^2 - 2 r + r^2))/C];
Derivative[1][tsMappedFunction[f_, g_, df_, dg_, C_, Rp_, sp_, \[Lambda]p_, a_, m_, \[Omega]_]][r:{__?NumericQ}] := Map[Derivative[1][tsMappedFunction[f, g, df, dg, C, Rp, sp, \[Lambda]p, a, m, \[Omega]]], r];

(* A radial function assembled from two pieces: fp (the Teukolsky-Starobinsky map of the flipped-spin
   integration) for r >= rc and fm (the integration at the original spin below rc) for r < rc *)
flippedFunction[fp_, fm_, rc_][r_?NumericQ] := If[r >= rc, fp[r], fm[r]];
flippedFunction[fp_, fm_, rc_][r:{__?NumericQ}] := Map[flippedFunction[fp, fm, rc], r];
Derivative[1][flippedFunction[fp_, fm_, rc_]][r_?NumericQ] := If[r >= rc, fp'[r], fm'[r]];
Derivative[1][flippedFunction[fp_, fm_, rc_]][r:{__?NumericQ}] := Map[Derivative[1][flippedFunction[fp, fm, rc]], r];

(* {f, g} for the given parameters, as expressions in the symbol r *)
teukolskyStarobinskyFlip[s_Integer?Positive, \[Lambda]_, a_, m_, \[Omega]_, r_] :=
  teukolskyStarobinskyFlipSymbolic[s][[2]] /. Thread[teukolskyStarobinskyFlipSymbolic[s][[1]] -> {\[Lambda], a, m, \[Omega], r}];


TeukolskyRadialNumericalIntegration[s_Integer, l_Integer, m_Integer, a_, \[Omega]_, \[Lambda]_, \[Nu]_, BCs_, norms_, {wp_, prec_, acc_}, opts:OptionsPattern[]] :=
 Module[{TRF, amps, ndsolveopts, psiopts, solFuncs, domains, Uptmp, Intmp, bmethod, flipUp, bdata, \[Nu]psi},
  (* Function to construct a single TeukolskyRadialFunction. For the "Up" solution of negative spin integrated at
     the flipped spin (see below), the radial function of spin s is the Teukolsky-Starobinsky map of the
     integrated one, divided by the constant (2 I omega)^(-2 s) that keeps unit transmission. *)
  TRF[bc_, ns_, sf_, domain_,  ndsolveopts___] :=
   Module[{solutionFunction, bcdir, amp, sInt, radialFunction, ft, gt, r = flipR, fPlus, fMinus, rc, Rc, dRc, \[Psi]c, d\[Psi]c, lower, goals},   (* flipR: a fixed private symbol as the variable of the pure functions built below, so that no Module-generated symbol is kept alive by them *)
    solutionFunction = sf[domain];
    bcdir = bc /. {"In" -> -1, "Up" -> +1};
    sInt = If[bc === "Up" && flipUp, -s, s];
    (*  Rescale amplitudes to give unit transmission coefficient (unit incidence where the transmission vanishes). *)
    amp = normaliseAmplitudes[ns, bc];
    radialFunction = If[sInt === s,
      Evaluate[#^-1 \[CapitalDelta][#,a]^-s Exp[bcdir I \[Omega] rs[#,a]] Exp[I m \[Phi]Reg[#,a]] solutionFunction[#]]&,
      (* Spin -s from the integrated spin +s: the Teukolsky-Starobinsky map beyond rc = r+ + 1, where the map's
         cancellation of the dominant Delta^s part of the spin +s solution costs nothing; below rc the spin s
         equation is integrated inwards from the mapped data at rc (stable over that short range), since the
         map loses about -2 s Log10[Delta] digits close to the horizon. *)
      {ft, gt} = teukolskyStarobinskyFlip[sInt, \[Lambda] + 2 s, a, m, \[Omega], r];
      rc = rp[a, 1] + 1;
      (* the map of the spin -s solution R+ (a pure function of r whose derivative involves only psi', which the
         integrators give accurately); its derivative needs R+'', which is taken from the radial equation rather
         than from the interpolating function, whose second derivative would cost two digits *)
      fPlus = With[{RpFun = Function @@ {r, r^-1 \[CapitalDelta][r,a]^-sInt Exp[bcdir I \[Omega] rs[r,a]] Exp[I m \[Phi]Reg[r,a]] solutionFunction[r]}, C = (2 I \[Omega])^(-2 s)},
        tsMappedFunction[Function @@ {r, ft}, Function @@ {r, gt}, Function @@ {r, D[ft, r]}, Function @@ {r, D[gt, r]}, C, RpFun, sInt, \[Lambda] + 2 s, a, m, \[Omega]]];
      If[domain =!= All && rc <= First[solutionFunction["Domain"]][[1]],
        fPlus,
        {Rc, dRc} = {fPlus[rc], fPlus'[rc]};
        {\[Psi]c, d\[Psi]c, rc} = Teukolsky`NumericalIntegration`Private`TeukolskyUpBCFromValues[s, m, a, \[Omega], rc, Rc, dRc];
        goals = Sequence[WorkingPrecision -> wp, PrecisionGoal -> prec, AccuracyGoal -> acc, Teukolsky`NumericalIntegration`Private`ndsolveOptions[ndsolveopts]];
        lower = If[domain === All,
          Teukolsky`NumericalIntegration`Private`AllIntegrator[s, \[Lambda], m, a, \[Omega], \[Psi]c, d\[Psi]c, rc, 1, goals],
          Teukolsky`NumericalIntegration`Private`Integrator[s, \[Lambda], m, a, \[Omega], \[Psi]c, d\[Psi]c, rc, First[solutionFunction["Domain"]][[1]], rc, 1, goals]];
        fMinus = Evaluate[#^-1 \[CapitalDelta][#,a]^-s Exp[bcdir I \[Omega] rs[#,a]] Exp[I m \[Phi]Reg[#,a]] lower[#]]&;
        flippedFunction[fPlus, fMinus, rc]]];
    TeukolskyRadialFunction[s, l, m, a, \[Omega],
     Association["s" -> s, "l" -> l, "m" -> m, "a" -> a, "\[Omega]" -> \[Omega], "Eigenvalue" -> reportedValue[\[Lambda], wp], "RenormalizedAngularMomentum" -> reportedValue[\[Nu], wp],
      "Method" -> {"NumericalIntegration", ndsolveopts},
      "BoundaryConditions" -> bc, "Amplitudes" -> amp, "UnscaledAmplitudes" -> ns,
      (* the domain requested; the integrated range can be larger (it includes the radius of the boundary data) *)
      "Domain" -> Which[domain === All, {rp[a, 1], \[Infinity]}, MatchQ[domain, {_?NumericQ, _?NumericQ}], domain, True, First[solutionFunction["Domain"]]],
      "RadialFunction" -> radialFunction
     ]
    ]
   ];

  (* Boundary data: series solutions of the integrator's equation at machine precision, precision-padded MST
     evaluations otherwise (see NumericalIntegration.m). The large-r series is evaluated where it converges, up
     to about 18/omega, and integrating the "Up" solution inwards from there is unstable for s < 0 (the ingoing
     solution grows like r^(-2s) relative to it: four digits lost per decade of radius for s = -2), so the "Up"
     solution of negative spin is integrated at the flipped spin -s, with eigenvalue lambda + 2 s, and mapped
     back with the Teukolsky-Starobinsky identity, for which that direction is stable. *)
  bmethod = OptionValue["BoundaryMethod"] /. Automatic -> If[wp === MachinePrecision, "Series", "MST"];
  If[!MatchQ[bmethod, "Series" | "MST"],
    Message[TeukolskyRadial::optx, "BoundaryMethod" -> OptionValue["BoundaryMethod"]];
    Return[$Failed];
  ];
  bdata = OptionValue["BoundaryData"];
  If[!(bdata === None || MatchQ[bdata, _Association?(AllTrue[Keys[#], MemberQ[{"In", "Up"}, #] &] && AllTrue[Values[#], MatchQ[{_?NumericQ, _?NumericQ, _?NumericQ}]] &)]),
    Message[TeukolskyRadial::optx, "BoundaryData" -> bdata];
    Return[$Failed];
  ];
  flipUp = bmethod === "Series" && s < 0 && !(AssociationQ[bdata] && KeyExistsQ[bdata, "Up"]);

  (* Domain over which the numerical solution can be evaluated *)
  domains = OptionValue["Domain"];
  If[ListQ[BCs],
    If[domains === All, domains = Thread[BCs -> All]];
    If[!MatchQ[domains, (List|Association)[Rule["In"|"Up",_?domainQ]..]],
      Message[TeukolskyRadial::dm, "Domain" -> domains, BCs];
      Return[$Failed];
    ];
    domains = Lookup[domains, BCs, None]; 
    If[!AllTrue[domains, domainQ],
      Message[TeukolskyRadial::dm, "Domain" -> OptionValue["Domain"], BCs];
      Return[$Failed];
    ];
  ,
    If[!domainQ[domains],
      Message[TeukolskyRadial::dm, "Domain" -> domains, BCs];
      Return[$Failed];
    ];
  ];
  
  (* Solution functions for the specified boundary conditions *)
  ndsolveopts = Sequence@@Join[FilterRules[{opts}, Options[NDSolve]], If[bdata =!= None, {"BoundaryData" -> bdata}, {}], If[OptionValue["BoundaryMethod"] =!= Automatic, {"BoundaryMethod" -> bmethod}, {}]];
  (* the boundary method actually used is passed to the integrators; the reported "Method" lists only the options given *)
  psiopts = Sequence[ndsolveopts, "BoundaryMethod" -> bmethod];
  (* psi has the NumericFunction attribute, so a non-numeric nu ("RenormalizedAngularMomentum" -> False gives
     Indeterminate) must be passed as a symbol or psi itself evaluates to Indeterminate; the integrators need
     nu only for their MST fall-backs, which compute it themselves when it is not numeric *)
  \[Nu]psi = If[NumericQ[\[Nu]], \[Nu], None];
  Uptmp = If[flipUp,
    (* the amplitudes and nu are those of spin s; should the series fail (a series radius beyond its cap, at
       very small omega) the MST fall-back of the boundary data must compute its own at spin -s, or the
       integrated solution would be normalised by the transmission amplitude of the other spin *)
    Teukolsky`NumericalIntegration`Private`psi[-s, \[Lambda] + 2 s, l, m, a, \[Omega], "Up", Automatic, None, WorkingPrecision -> wp, PrecisionGoal -> prec, AccuracyGoal -> acc, psiopts],
    Teukolsky`NumericalIntegration`Private`psi[s, \[Lambda], l, m, a, \[Omega], "Up", norms, \[Nu]psi, WorkingPrecision -> wp, PrecisionGoal -> prec, AccuracyGoal -> acc, psiopts]];
  Intmp = Teukolsky`NumericalIntegration`Private`psi[s, \[Lambda], l, m, a, \[Omega], "In", norms, \[Nu]psi, WorkingPrecision -> wp, PrecisionGoal -> prec, AccuracyGoal -> acc, psiopts];
  solFuncs =
   <|"Up" :> Uptmp,
     "In" :> Intmp
	 |>;
  solFuncs = Lookup[solFuncs, BCs];

  (* Select normalisation coefficients for the specified boundary conditions *)
  amps = Lookup[norms, BCs];

  If[ListQ[BCs],
    Return[Association[MapThread[#1 -> TRF[#1, #2, #3, #4, ndsolveopts]&, {BCs, amps, solFuncs, domains}]]],
    Return[TRF[BCs, amps, solFuncs, domains, ndsolveopts]]
  ];
];


(* ::Subsection::Closed:: *)
(*Sasaki-Nakamura Method*)


Options[TeukolskyRadialSasakiNakamura] = Join[
  {"Domain" -> None},
  FilterRules[Options[NDSolve], Except[WorkingPrecision|AccuracyGoal|PrecisionGoal]]];


domainQ[domain_] := MatchQ[domain, {_?NumericQ, _?NumericQ} | (_?NumericQ) | All];


TeukolskyRadialSasakiNakamura[s_Integer, l_Integer, m_Integer, a_, \[Omega]_, \[Lambda]_, \[Nu]_, BCs_, norms_, {wp_, prec_, acc_}, opts:OptionsPattern[]] :=
 Module[{TRF, amps, ndsolveopts, solFuncs, domains},
  (* Function to construct a single TeukolskyRadialFunction *)
  TRF[bc_, ns_, sf_, domain_, ndsolveopts___] :=
   Module[{solutionFunction, amp},
    If[sf === $Failed, Return[$Failed]];
    solutionFunction = sf[domain];
    (*  Rescale amplitudes to give unit transmission coefficient (unit incidence where the transmission vanishes). *)
    amp = normaliseAmplitudes[ns, bc];
    TeukolskyRadialFunction[s, l, m, a, \[Omega],
     Association["s" -> s, "l" -> l, "m" -> m, "a" -> a, "\[Omega]" -> \[Omega], "Eigenvalue" -> reportedValue[\[Lambda], wp], "RenormalizedAngularMomentum" -> reportedValue[\[Nu], wp],
      "Method" -> {"SasakiNakamura", ndsolveopts},
      "BoundaryConditions" -> bc, "Amplitudes" -> amp, "UnscaledAmplitudes" -> ns,
      (* the domain requested; the integrated range can be larger (it includes the radius of the boundary data) *)
      "Domain" -> Which[domain === All, {rp[a, 1], \[Infinity]}, MatchQ[domain, {_?NumericQ, _?NumericQ}], domain, True, First[solutionFunction["Domain"]]],
      "RadialFunction" -> solutionFunction
     ]
    ]
   ];

  (* Domain over which the numerical solution can be evaluated *)
  domains = OptionValue["Domain"];
  If[ListQ[BCs],
    If[!MatchQ[domains, (List|Association)[Rule["In"|"Up",_?domainQ]..]],
      Message[TeukolskyRadial::dm, "Domain" -> domains, BCs];
      Return[$Failed];
    ];
    domains = Lookup[domains, BCs, None]; 
    If[!AllTrue[domains, domainQ],
      Message[TeukolskyRadial::dm, "Domain" -> OptionValue["Domain"], BCs];
      Return[$Failed];
    ];
  ,
    If[!domainQ[domains],
      Message[TeukolskyRadial::dm, "Domain" -> domains, BCs];
      Return[$Failed];
    ];
  ];

  (* Solution functions for the specified boundary conditions *)
  ndsolveopts = Sequence@@FilterRules[{opts}, Options[NDSolve]];
  solFuncs =
   <|"In" :> $Failed,
     "Up" :> Teukolsky`SasakiNakamura`Private`TeukolskyRadialUp[s, \[Lambda], m, a, \[Omega](*, WorkingPrecision -> wp, PrecisionGoal -> prec, AccuracyGoal -> acc, ndsolveopts*)]
     |>;
  solFuncs = Lookup[solFuncs, BCs];

  (* Select normalisation coefficients for the specified boundary conditions *)
  amps = Lookup[norms, BCs];

  If[ListQ[BCs],
    Return[Association[MapThread[#1 -> TRF[#1, #2, #3, #4, ndsolveopts]&, {BCs, amps, solFuncs, domains}]]],
    Return[TRF[BCs, amps, solFuncs, domains, ndsolveopts]]
  ];
];


(* ::Subsection::Closed:: *)
(*MST Method*)


Options[TeukolskyRadialMST] = {};


(* Machine-precision Automatic method: the numerical-integration solutions with precision and accuracy
   goals two digits below machine precision. Their boundary data are series solutions of the integrator's
   equation, a power series about the horizon for "In" and the large-r asymptotic series for "Up" (see
   NumericalIntegration.m), so they are accurate to ~1e-12 - 1e-15, and being interpolating or integrated
   solutions they are cheap to evaluate on many radii, unlike the MST series which are summed at every
   point. The accuracy of the pair is estimated from the Wronskian, which must equal 2 i omega B^inc C^trans
   (with the amplitudes from the MST formulae) and be independent of r; a poor estimate is reported. *)
TeukolskyRadial::acc = "The estimated relative accuracy of the radial functions is only `1`; use a higher WorkingPrecision for better accuracy.";

radialAccuracyEstimate[R_Association, s_Integer, a_, \[Omega]_] :=
 Module[{rp1 = rp[a, 1], W, Wexact, w1, w2},
  (* evaluated on lists of radii, so that each function is integrated once *)
  W[rs_List] := Module[{i, di, u, du}, {i, di, u, du} = {R["In"][rs], R["In"]'[rs], R["Up"][rs], R["Up"]'[rs]}; (rs^2 - 2 rs + a^2)^(s + 1) (i du - u di)];
  Wexact = 2 I \[Omega] R["In"]["Amplitudes"]["Incidence"] R["Up"]["Amplitudes"]["Transmission"];
  {w1, w2} = Quiet[W[{2. rp1, 10. rp1}]];
  If[!(NumericQ[w1] && NumericQ[w2] && NumericQ[Wexact]) || w2 == 0 || Wexact == 0, Return[Infinity]];
  Max[Abs[w1/Wexact - 1], Abs[w2/Wexact - 1], Abs[w1/w2 - 1]]
 ];

(* The accuracy estimate needs both solutions and the amplitudes B^inc and C^trans; it is made only when
   they are all there (both solutions requested and the amplitudes computed), so that asking for one solution
   costs one integration and disabling the amplitudes does not produce a meaningless warning. *)
TeukolskyRadialAutomaticMachinePrecision[s_Integer, l_Integer, m_Integer, a_, \[Omega]_, \[Lambda]_, \[Nu]_, BCs_, norms_, {wp_, prec_, acc_}, opts:OptionsPattern[]] :=
 Module[{R, e},
  (* the user's goals; by default PrecisionGoal is already WorkingPrecision - 2, and the absolute AccuracyGoal,
     Infinity for the MST series, is set to the same for the integration *)
  R = TeukolskyRadialNumericalIntegration[s, l, m, a, \[Omega], \[Lambda], \[Nu], BCs, norms, {wp, prec, If[acc === Infinity, $MachinePrecision - 2, acc]}, opts];
  If[ListQ[BCs] && ContainsAll[BCs, {"In", "Up"}] && AssociationQ[R] && NumericQ[Lookup[norms["In"], "Incidence", None]] && NumericQ[Lookup[norms["Up"], "Transmission", None]] && norms["Up"]["Transmission"] != 0,
    e = radialAccuracyEstimate[R, s, a, \[Omega]];
    (* warned about only when the estimate falls short of what the goals asked for: with PrecisionGoal -> 6
       the expected 1e-5 was reported as a failure *)
    If[e > Max[10^-6, 10^(2 - Min[prec, If[acc === Infinity, $MachinePrecision - 2, acc]])], Message[TeukolskyRadial::acc, N[e, 2]]];
  ];
  R
 ];

Options[TeukolskyRadialAutomaticMachinePrecision] = Options[TeukolskyRadialNumericalIntegration];


TeukolskyRadialMST[s_Integer, l_Integer, m_Integer, a_, \[Omega]_, \[Lambda]_, \[Nu]_, BCs_, norms_, {wp_, prec_, acc_}, opts:OptionsPattern[]] :=
 Module[{amps, solFuncs, TRF},
  (* Function to construct a TeukolskyRadialFunction *)
  TRF[bc_, ns_, sf_] := Module[{amp},
    (*  Rescale amplitudes to give unit transmission coefficient (unit incidence where the transmission vanishes). *)
    amp = normaliseAmplitudes[ns, bc];
    TeukolskyRadialFunction[s, l, m, a, \[Omega],
     Association["s" -> s, "l" -> l, "m" -> m, "a" -> a, "\[Omega]" -> \[Omega], "Eigenvalue" -> reportedValue[\[Lambda], wp], "RenormalizedAngularMomentum" -> reportedValue[\[Nu], wp],
      "Method" -> {"MST"},
      "BoundaryConditions" -> bc, "Amplitudes" -> amp, "UnscaledAmplitudes" -> ns,
      "Domain" -> {rp[a, 1], \[Infinity]}, "RadialFunction" -> sf
     ]
    ]
  ];

  (* Solution functions for the specified boundary conditions *)
  solFuncs =
    <|"In" :> Teukolsky`MST`MST`Private`MSTRadialIn[s,l,m,a,2\[Omega],\[Nu],\[Lambda],norms["In", normalisationKey[norms["In"], "In"]], {wp, prec, acc}],
      "Up" :> Teukolsky`MST`MST`Private`MSTRadialUp[s,l,m,a,2\[Omega],\[Nu],\[Lambda],norms["Up", "Transmission"], {wp, prec, acc}]|>;
  solFuncs = Lookup[solFuncs, BCs];

  (* Select normalisation coefficients for the specified boundary conditions *)
  amps = Lookup[norms, BCs];

  If[ListQ[BCs],
    Return[Association[MapThread[#1 -> TRF[#1, #2, #3]&, {BCs, amps, solFuncs}]]],
    Return[TRF[BCs, amps, solFuncs]]
  ];
];


(* ::Subsection::Closed:: *)
(*HeunC*)


Options[TeukolskyRadialHeunC] = {};


TeukolskyRadialHeunC[s_Integer, l_Integer, m_Integer, a_, \[Omega]_, \[Lambda]_, \[Nu]_, BCs_, norms_, {wp_, prec_, acc_}, opts:OptionsPattern[]] :=
 Module[{amps, solFuncs, TRF, \[Omega]c, \[Lambda]c},
  (* The HeunC method is only supported on version 12.1 and newer *)
  If[$VersionNumber < 12.1,
    Message[TeukolskyRadial::hc];
    Return[$Failed]
  ];

  (* Function to construct a TeukolskyRadialFunction *)
  TRF[bc_, ns_, sf_] := Module[{amp},
    (*  Rescale amplitudes to give unit transmission coefficient (unit incidence where the transmission vanishes). *)
    amp = normaliseAmplitudes[ns, bc];
    If[sf === $Failed, $Failed,
      TeukolskyRadialFunction[s, l, m, a, \[Omega],
        Association["s" -> s, "l" -> l, "m" -> m, "a" -> a, "\[Omega]" -> \[Omega], "Eigenvalue" -> reportedValue[\[Lambda], wp], "RenormalizedAngularMomentum" -> reportedValue[\[Nu], wp],
          "Method" -> {"HeunC"},
          "BoundaryConditions" -> bc, "Amplitudes" -> amp, "UnscaledAmplitudes" -> ns,
          "Domain" -> {rp[a, 1], \[Infinity]}, "RadialFunction" -> sf
        ]
      ]
    ]
  ];

  (* Solution functions for the specified boundary conditions *)
  \[Omega]c = Conjugate[\[Omega]];
  \[Lambda]c = Conjugate[\[Lambda]];
  solFuncs =
    <|"In" :> ((2^((-I a m-Sqrt[1-a^2] s+2 I \[Omega])/Sqrt[1-a^2]) (1-a^2)^((I a m)/(2 (1+Sqrt[1-a^2]))-s-I \[Omega]) E^(1/2 I (a m-2 # \[Omega])) ((-1-Sqrt[1-a^2]+#)/Sqrt[1-a^2])^((I a m)/(2 Sqrt[1-a^2])-s-I (1+1/Sqrt[1-a^2]) \[Omega]) ((-1+Sqrt[1-a^2]+#)/Sqrt[1-a^2])^((I (a m+2 (-1+Sqrt[1-a^2]) \[Omega]))/(2 Sqrt[1-a^2])) HeunC[s+s^2+\[Lambda]-(a m-2 \[Omega])^2/(-1+a^2)-4 \[Omega]^2+(I (-a m-4 (-1+s) \[Omega]+2 a^2 (-1+2 s) \[Omega]))/Sqrt[1-a^2],-4 \[Omega] (-I Sqrt[1-a^2]+a m+I Sqrt[1-a^2] s-2 \[Omega]+2 Sqrt[1-a^2] \[Omega]),1-s+(I (a m-2 \[Omega]))/Sqrt[1-a^2]-2 I \[Omega],1+s+(I (a m-2 \[Omega]))/Sqrt[1-a^2]+2 I \[Omega],4 I Sqrt[1-a^2] \[Omega],(1+Sqrt[1-a^2]-#)/(2 Sqrt[1-a^2])])&),
      "Up" :> (norms["Up"]["Reflection"]/norms["Up"]["Transmission"](2^((-I a m-Sqrt[1-a^2] s+2 I \[Omega])/Sqrt[1-a^2]) (1-a^2)^((I a m)/(2 (1+Sqrt[1-a^2]))-s-I \[Omega]) E^(1/2 I (a m-2 # \[Omega])) ((-1-Sqrt[1-a^2]+#)/Sqrt[1-a^2])^((I a m)/(2 Sqrt[1-a^2])-s-I (1+1/Sqrt[1-a^2]) \[Omega]) ((-1+Sqrt[1-a^2]+#)/Sqrt[1-a^2])^((I (a m+2 (-1+Sqrt[1-a^2]) \[Omega]))/(2 Sqrt[1-a^2])) HeunC[s+s^2+\[Lambda]-(a m-2 \[Omega])^2/(-1+a^2)-4 \[Omega]^2+(I (-a m-4 (-1+s) \[Omega]+2 a^2 (-1+2 s) \[Omega]))/Sqrt[1-a^2],-4 \[Omega] (-I Sqrt[1-a^2]+a m+I Sqrt[1-a^2] s-2 \[Omega]+2 Sqrt[1-a^2] \[Omega]),1-s+(I (a m-2 \[Omega]))/Sqrt[1-a^2]-2 I \[Omega],1+s+(I (a m-2 \[Omega]))/Sqrt[1-a^2]+2 I \[Omega],4 I Sqrt[1-a^2] \[Omega],(1+Sqrt[1-a^2]-#)/(2 Sqrt[1-a^2])])+
              norms["Up"]["Incidence"]/norms["Up"]["Transmission"](#^2-2 #+a^2)^-s (2^((I a m-Sqrt[1-a^2] (-s)-2 I \[Omega]c)/Sqrt[1-a^2]) (1-a^2)^((-I a m)/(2 (1+Sqrt[1-a^2]))-(-s)+I \[Omega]c) E^(-1/2 I (a m-2 # \[Omega]c)) ((-1-Sqrt[1-a^2]+#)/Sqrt[1-a^2])^((-I a m)/(2 Sqrt[1-a^2])-(-s)+I (1+1/Sqrt[1-a^2]) \[Omega]c) ((-1+Sqrt[1-a^2]+#)/Sqrt[1-a^2])^((-I (a m+2 (-1+Sqrt[1-a^2]) \[Omega]c))/(2 Sqrt[1-a^2])) HeunC[(-s)+(-s)^2+(\[Lambda]c+2s)-(a m-2 \[Omega]c)^2/(-1+a^2)-4 \[Omega]c^2+(-I (-a m-4 (-1+(-s)) \[Omega]c+2 a^2 (-1+2 (-s)) \[Omega]c))/Sqrt[1-a^2],-4 \[Omega]c (+I Sqrt[1-a^2]+a m-I Sqrt[1-a^2] (-s)-2 \[Omega]c+2 Sqrt[1-a^2] \[Omega]c),1-(-s)+(-I (a m-2 \[Omega]c))/Sqrt[1-a^2]+2 I \[Omega]c,1+(-s)+(-I (a m-2 \[Omega]c))/Sqrt[1-a^2]-2 I \[Omega]c,-4 I Sqrt[1-a^2] \[Omega]c,(1+Sqrt[1-a^2]-#)/(2 Sqrt[1-a^2])])&) |>;
  solFuncs = Lookup[solFuncs, BCs];

  (* Select normalisation coefficients for the specified boundary conditions *)
  amps = Lookup[norms, BCs];

  If[ListQ[BCs],
    Return[Association[MapThread[#1 -> TRF[#1, #2, #3]&, {BCs, amps, solFuncs}]]],
    Return[TRF[BCs, amps, solFuncs]]
  ];
];


(* ::Subsection::Closed:: *)
(*Static modes*)


staticAmplitudes[s_, l_, m_, a_] :=
 Module[{\[Tau] = -((m a)/Sqrt[1-a^2]), \[Kappa] = Sqrt[1 - a^2], ampIn1, ampIn2, ampUp3, ampUp4, ampUp5},
  (* Return results as an Association *)
  If[\[Tau]==0,
     <|"In" -> <|
         "\[ScriptCapitalH]" -> 1,
         "\[ScriptCapitalI]" -> ((2 \[Kappa])^(-l+Abs[s]) Gamma[2l+1] Gamma[1+Abs[s]])/(Gamma[l+1] Gamma[1+l+Abs[s]])|>,
       "Up" -> <|
         "\[ScriptCapitalH]" ->
           If[s==0,
             -((2 \[Kappa])^(-1-l) Gamma[3+2l])/(2 Gamma[l+2] Gamma[1+l]),
             ((2 \[Kappa])^(-1-l+Abs[s]) Gamma[2+2 l] Gamma[Abs[s]])/(Gamma[1+l+Abs[s]] Gamma[1+l])],
         "\[ScriptCapitalI]" -> 1|>|>
  ,
     <|"In" -> <|
         "\[ScriptCapitalH]" -> 1,
         "\[ScriptCapitalI]-" -> -((2^(l-s-I \[Tau]) \[Kappa]^(1-s) (-1)^(l+s) Gamma[1+l+s] Gamma[1-s-I \[Tau]])/(Gamma[2+2 l] Gamma[-l-I \[Tau]])),
         "\[ScriptCapitalI]+" -> ((2 \[Kappa])^(-l-s-I \[Tau]) Gamma[2l+1] Gamma[1-s-I \[Tau]])/(Gamma[1+l-s] Gamma[1+l-I \[Tau]])|>,
       "Up" -> <|
         "\[ScriptCapitalH]-" -> ((2 \[Kappa])^(-1-l+s+I \[Tau]) Gamma[2+2 l] Gamma[s+I \[Tau]])/(Gamma[1+l+s] Gamma[1+l+I \[Tau]]),
         "\[ScriptCapitalH]+" -> ((2 \[Kappa])^(-1-l-s-I \[Tau]) Gamma[2+2 l] Gamma[-s-I \[Tau]])/(Gamma[1+l-s] Gamma[1+l-I \[Tau]]),
         "\[ScriptCapitalI]" -> 1|>|>
  ]
 ];


TeukolskyRadialStatic[s_Integer, l_Integer, m_Integer, a_, \[Omega]_, \[Lambda]_, \[Nu]_, BCs_, norms_] :=
 Module[{amps, solFuncs, TRF},
  (* Function to construct a TeukolskyRadialFunction *)
  TRF[bc_, amp_, sf_] :=
    TeukolskyRadialFunction[s, l, m, a, \[Omega],
     Association["s" -> s, "l" -> l, "m" -> m, "a" -> a, "\[Omega]" -> \[Omega], "Eigenvalue" -> \[Lambda], "RenormalizedAngularMomentum" -> \[Nu],
      "Method" -> {"Static"},
      "BoundaryConditions" -> bc, "Amplitudes" -> amp, "UnscaledAmplitudes" -> amp,
      "Domain" -> {rp[a, 1], \[Infinity]}, "RadialFunction" -> sf
     ]
    ];

  (* Solution functions for the specified boundary conditions *)
  With[{\[Tau] = -((m a)/Sqrt[1-a^2]), \[Kappa] = Sqrt[1 - a^2]},
    With[{normIn = If[\[Tau]==0 && s>0, Gamma[s+1]Pochhammer[l+s+1,-2s], (2\[Kappa])^(-2s-I \[Tau]) Gamma[1-s-I \[Tau]]],
          normUp = (2 \[Kappa])^(-s-l-1)},
    solFuncs =
      <|"In" :> (normIn (-(1 + \[Kappa] - #)/(2 \[Kappa]))^(-s - I \[Tau]/2) (1 - (1 + \[Kappa] - #)/(2 \[Kappa]))^(-I \[Tau]/2) Hypergeometric2F1Regularized[-l - I \[Tau], l + 1 - I \[Tau], 1 - s - I \[Tau], (1 + \[Kappa] - #)/(2 \[Kappa])]&),
        "Up" :> (normUp (-(1 + \[Kappa] - #)/(2 \[Kappa]))^(-s - (I \[Tau])/2) (1 - (1 + \[Kappa] - #)/(2 \[Kappa]))^((I \[Tau])/2 - l - 1) Hypergeometric2F1[l + 1 - I \[Tau], l + 1 - s, 2 l + 2, 1/(1 - (1 + \[Kappa] - #)/(2 \[Kappa]))]&)
       |>;
    ];
  ];
  solFuncs = Lookup[solFuncs, BCs];

  (* Select normalisation coefficients for the specified boundary conditions *)
  amps = Lookup[norms, BCs];

  If[ListQ[BCs],
    Return[Association[MapThread[#1 -> TRF[#1, #2, #3]&, {BCs, amps, solFuncs}]]],
    Return[TRF[BCs, amps, solFuncs]]
  ];
];


(* ::Subsection::Closed:: *)
(*TeukolskyRadial*)


SyntaxInformation[TeukolskyRadial] =
 {"ArgumentsPattern" -> {_, _, _, _, _, OptionsPattern[]}};


Options[TeukolskyRadial] = {
  Method -> Automatic,
  "BoundaryConditions" -> {"In", "Up"},
  "Amplitudes" -> Automatic,
  "RenormalizedAngularMomentum" -> Automatic,
  "Eigenvalue" -> Automatic,
  WorkingPrecision -> Automatic,
  PrecisionGoal -> Automatic,
  AccuracyGoal -> Automatic,
  "WronskianCheck" -> Automatic
};


TeukolskyRadial[s_?NumericQ, l_?NumericQ, m_?NumericQ, a_, \[Omega]_, OptionsPattern[]] /;
  l < Abs[s] || Abs[m] > l || !AllTrue[{2s, 2l, 2m}, IntegerQ] || !IntegerQ[l-s] || !IntegerQ[m-s] := 
 (Message[TeukolskyRadial::params, s, l, m]; $Failed);


TeukolskyRadial[s_, l_, m_, a_Complex, \[Omega]_, OptionsPattern[]] :=
 (Message[TeukolskyRadial::cmplx, a]; $Failed);

(* |a| >= 1 is not a black hole (no horizon at a > 1, degenerate horizons at a = 1, where kappa = 0): the methods
   took Sqrt[1 - a^2] and hung (Min::nord without end at a = 1.2) or hit $RecursionLimit (a = 1) *)
TeukolskyRadial[s_, l_, m_, a_?NumericQ, \[Omega]_, OptionsPattern[]] /; Im[a] == 0 && Abs[a] >= 1 :=
 (Message[TeukolskyRadial::spin, a]; $Failed);


(* ::Subsubsection::Closed:: *)
(*Static modes*)


TeukolskyRadial[s_Integer, l_Integer, m_Integer, a_, \[Omega]_, opts:OptionsPattern[]] /; AllTrue[{a, \[Omega]}, NumericQ] && \[Omega] == 0 :=
 Module[{\[Lambda], BCs, norms, wp, prec, acc},
  (* Determine which boundary conditions the homogeneous solution(s) should satisfy *)
  BCs = OptionValue["BoundaryConditions"];
  If[!MatchQ[BCs, "In"|"Up"|{("In"|"Up")..}], 
    Message[TeukolskyRadial::optx, "BoundaryConditions" -> BCs];
    Return[$Failed];
  ];

  (* Eigenvalue *)
  \[Lambda] = SpinWeightedSpheroidalEigenvalue[s, l, m, a \[Omega]];

  (* Some options are not supported for static modes; the Wronskian check is one of them (there are no MST
     solutions to check), and an invalid value of it is rejected as in the non-static case *)
  If[!MatchQ[OptionValue["WronskianCheck"], Automatic | True | False],
    Message[TeukolskyRadial::optx, "WronskianCheck" -> OptionValue["WronskianCheck"]];
    Return[$Failed];
  ];
  Do[
    If[OptionValue[opt] =!= Automatic, Message[TeukolskyRadial::sopt, opt]];,
    {opt, {"Eigenvalue", "RenormalizedAngularMomentum", Method, WorkingPrecision, PrecisionGoal, AccuracyGoal, "WronskianCheck"}}
  ];

  (* Compute the asymptotic amplitudes *)
  Which[
  OptionValue["Amplitudes"] === False,
    norms = <|"In" -> <|"Transmission" -> 1|>, "Up" -> <|"Transmission" -> 1|>|>;,
  MatchQ[OptionValue["Amplitudes"], <|"In"-><|___|>, "Up" -> <|___|>|>],
    norms = OptionValue["Amplitudes"];,
  MatchQ[OptionValue["Amplitudes"], Automatic|True],
    norms = staticAmplitudes[s, l, m, a];,
  True,
    Message[TeukolskyRadial::optx, "Amplitudes" -> OptionValue["Amplitudes"]];
    Return[$Failed];
  ];

  (* Call the chosen implementation *)
  TeukolskyRadialStatic[s, l, m, a, \[Omega], \[Lambda], \[Lambda], BCs, norms]
]


(* ::Subsubsection::Closed:: *)
(*Non-static modes*)


(* Exact arguments are evaluated at the WorkingPrecision, which must then be given: silently returning
   unevaluated let symbolic expressions propagate through calling code. *)
TeukolskyRadial[s_Integer, l_Integer, m_Integer, a_, \[Omega]_, opts:OptionsPattern[]] /; AllTrue[{a, \[Omega]}, NumericQ] && \[Omega] != 0 && !(InexactNumberQ[a] || InexactNumberQ[\[Omega]]) :=
 Module[{wp = OptionValue[WorkingPrecision]},
  If[wp === Automatic,
    Message[TeukolskyRadial::exact, a, \[Omega]];
    Return[$Failed];
  ];
  TeukolskyRadial[s, l, m, SetPrecision[a, wp], SetPrecision[\[Omega], wp], opts]
 ];


TeukolskyRadial[s_Integer, l_Integer, m_Integer, a_, \[Omega]_, opts:OptionsPattern[]] /; AllTrue[{a, \[Omega]}, NumericQ] && (InexactNumberQ[a] || InexactNumberQ[\[Omega]]) :=
 Module[{TRF, subopts, BCs, norms, \[Nu], \[Lambda], wp, prec, acc, compute, check, extra, wpn, tol, res, e, k, ampPadding, ampRetried},
  (* Extract suboptions from Method to be passed on. *)
  If[ListQ[OptionValue[Method]],
    subopts = Rest[OptionValue[Method]];,
    subopts = {};
  ];

  (* Determine which boundary conditions the homogeneous solution(s) should satisfy *)
  BCs = OptionValue["BoundaryConditions"];
  If[!MatchQ[BCs, "In"|"Up"|{("In"|"Up")..}], 
    Message[TeukolskyRadial::optx, "BoundaryConditions" -> BCs];
    Return[$Failed];
  ];
  If[!MatchQ[OptionValue["WronskianCheck"], Automatic | True | False],
    Message[TeukolskyRadial::optx, "WronskianCheck" -> OptionValue["WronskianCheck"]];
    Return[$Failed];
  ];

  (* Options associated with precision and accuracy *)
  {wp, prec, acc} = OptionValue[{WorkingPrecision, PrecisionGoal, AccuracyGoal}];
  If[wp === Automatic, wp = Precision[{a, \[Omega]}]];
  (* The MST series are truncated when the terms fall below 10^-prec relative to the sum, and the
     numerical integration uses prec as its PrecisionGoal; the former default of wp/2 limited the
     solutions to ~8 digits at machine precision (and to ~16 digits with 32-digit input), so the
     default is two digits below the working precision. *)
  If[prec === Automatic, prec = wp - 2];
  If[acc === Automatic, acc = Infinity];
  If[Precision[a] < wp, Message[TeukolskyRadial::precw, "a", a, wp]];
  If[Precision[\[Omega]] < wp, Message[TeukolskyRadial::precw, "\[Omega]", \[Omega], wp]];

  (* Decide which implementation to use *)
  Switch[OptionValue[Method],
    Automatic,
      (* At a complex frequency the numerical integration is ill-conditioned at machine precision: the two
         solutions differ by Exp[2 I omega r*], of modulus Exp[2 |Im omega| r*], so the "In" solution is
         exponentially subdominant outwards and the "Up" solution inwards, and any error in the boundary data is
         amplified by Exp[2 |Im omega| (range of tortoise coordinate)] (1e-15 to 1e-5.5 over 40 in r at Im omega = -1/5). The MST
         series do not integrate and keep 1e-13 there, so they are the default for complex frequencies. *)
      If[wp === MachinePrecision && Im[\[Omega]] == 0,
         TRF = TeukolskyRadialAutomaticMachinePrecision,
         TRF = TeukolskyRadialMST],
    "MST" | {"MST", OptionsPattern[TeukolskyRadialMST]},
      TRF = TeukolskyRadialMST,
    "NumericalIntegration" | {"NumericalIntegration", OptionsPattern[TeukolskyRadialNumericalIntegration]},
      TRF = TeukolskyRadialNumericalIntegration;,
    "HeunC",
      (* Some options are not supported for the HeunC method modes *)
      Do[
        If[OptionValue[opt] =!= Automatic, Message[TeukolskyRadial::hcopt, opt]];,
        {opt, {WorkingPrecision, PrecisionGoal, AccuracyGoal}}
      ];
      TRF = TeukolskyRadialHeunC;,
    _,
      Message[TeukolskyRadial::optx, Method -> OptionValue[Method]];
      Return[$Failed];
  ];

  (* Check only supported sub-options have been specified; options of TeukolskyRadial itself
     (e.g. "Eigenvalue", "RenormalizedAngularMomentum") are reported as being in the wrong place *)
  Module[{unknown = Complement[subopts, FilterRules[subopts, Options[TRF]]]},
    If[unknown =!= {},
      With[{misplaced = Select[unknown, MemberQ[Keys[Options[TeukolskyRadial]], First[#]] &]},
        If[misplaced =!= {}, Message[TeukolskyRadial::topopt, Keys[misplaced], First[OptionValue[Method]]]];
        (* sub-options that are neither supported nor misplaced are still reported *)
        If[Complement[unknown, misplaced] =!= {}, Message[TeukolskyRadial::optx, Method -> OptionValue[Method]]]
      ]
    ];
    subopts = FilterRules[subopts, Options[TRF]];
  ];

  (* Supplied eigenvalue and nu are used as given by the MST series, never recomputed at their padded precision *)
  Teukolsky`MST`MST`Private`registerSupplied[If[NumericQ[OptionValue["Eigenvalue"]], OptionValue["Eigenvalue"], None],
    If[NumericQ[OptionValue["RenormalizedAngularMomentum"]], OptionValue["RenormalizedAngularMomentum"], None]];

  (* Eigenvalue *)
  Which[
  OptionValue["Eigenvalue"] === False,
    \[Lambda] = Indeterminate;,
  NumericQ[OptionValue["Eigenvalue"]],
    \[Lambda] = OptionValue["Eigenvalue"];,
  True,
    \[Lambda] = SpinWeightedSpheroidalEigenvalue[s, l, m, a \[Omega]];
  ];

  (* The renormalized angular momentum, the asymptotic amplitudes and the radial functions, computed with
     extra working precision beyond the padding of the individual evaluations when the Wronskian check
     below has found that necessary for this mode (see mstWronskianError). *)
  compute[extra_, bcs_] := Module[{},
  (* Renormalized angular momentum *)
    Which[
    OptionValue["RenormalizedAngularMomentum"] === False,
      \[Nu] = Indeterminate;,
    NumericQ[OptionValue["RenormalizedAngularMomentum"]],
      \[Nu] = OptionValue["RenormalizedAngularMomentum"];,
    True,
      \[Nu] = paddedComputation[RenormalizedAngularMomentum[s, l, m, SetPrecision[a, #], SetPrecision[\[Omega], #], SetPrecision[\[Lambda], #], Method -> (OptionValue["RenormalizedAngularMomentum"] /. (Automatic|True) -> "Monodromy")] &, wp, "The renormalized angular momentum", extra];
    ];
    (* without nu neither the amplitudes nor the MST solutions can be computed: fail now, rather than let the
       padded evaluations retry at ever higher precision (over five minutes at omega = 270 before this check) *)
    If[!NumericQ[\[Nu]] && OptionValue["RenormalizedAngularMomentum"] =!= False,
      Message[TeukolskyRadial::nufail, s, l, m, a, \[Omega]];
      Return[$Failed, Module];
    ];

    (* Compute the asymptotic amplitudes *)
    Which[
    OptionValue["Amplitudes"] === False,
      norms = <|"In" -> <|"Transmission" -> 1|>, "Up" -> <|"Transmission" -> 1|>|>;,
    MatchQ[OptionValue["Amplitudes"], <|"In"-><|___|>, "Up" -> <|___|>|>],
      norms = OptionValue["Amplitudes"];,
    MatchQ[OptionValue["Amplitudes"], Automatic|True],
      If[OptionValue["RenormalizedAngularMomentum"] === False,
        Message[TeukolskyRadial::opti, {"Amplitudes" -> OptionValue["Amplitudes"], "RenormalizedAngularMomentum" -> OptionValue["RenormalizedAngularMomentum"]}];
        Return[$Failed];
      ];
      (* the eigenvalue and nu refined to the padded precision of the amplitudes are kept for the radial
         functions, which are evaluated at a similar padded precision *)
      {$refinedEigenvalue, $refinedNu} = {\[Lambda], \[Nu]};
      With[{n = epsilonPlusDegeneracy[s, m, a, \[Omega], wp]},
        norms = Block[{$acceptIndeterminate = n =!= None, $degeneracyOrder = n},
          paddedComputation[mstAmplitudes[s, l, m, a, \[Omega], \[Lambda], \[Nu], #, prec, acc] &, wp, "The asymptotic amplitudes", extra]];
        If[n =!= None, degeneracyMessages[n, \[Omega], norms]];
        (* away from a degeneracy every amplitude should be a number; the formulae overflow at large |omega|
           (omega = 50 for l = 2, say), and without this message the only symptoms were a precision of
           Indeterminate for the MST "Up" solution and an accuracy estimate of Infinity *)
        If[n === None && AssociationQ[norms],
          With[{bad = Flatten[Table[If[!NumericQ[Lookup[norms[bc], key, 0]] || (key === normalisationKey[norms[bc], bc] && norms[bc][key] == 0), {bc, key}, Nothing], {bc, {"In", "Up"}}, {key, {"Incidence", "Transmission", "Reflection"}}], 1]},
            If[bad =!= {}, Message[TeukolskyRadial::ampfail, \[Omega], bad]]]];
        (* At a degeneracy the refined values belong to the exactly degenerate frequency of mstAmplitudes, not
           to omega as given (1e-17 apart at machine precision), and the radial functions must refine their
           own: the MST series amplify an inconsistency between nu and omega by the digits they lose to
           cancellation (1e-7 in the functions at machine precision when nu was carried over). *)
        (* refined values replace the given ones only where they agree (a supplied value of another problem is kept) *)
        If[n === None,
          If[NumericQ[$refinedNu] && Precision[$refinedNu] > Precision[\[Nu]] && Teukolsky`MST`MST`Private`consistentQ[$refinedNu, \[Nu]], \[Nu] = $refinedNu];
          If[NumericQ[$refinedEigenvalue] && Precision[$refinedEigenvalue] > Precision[\[Lambda]] && Teukolsky`MST`MST`Private`consistentQ[$refinedEigenvalue, \[Lambda]], \[Lambda] = $refinedEigenvalue]];
      ];
      {ampPadding, ampRetried} = {$lastPaddingDigits, $lastPaddingRetried};,
    True,
      Message[TeukolskyRadial::optx, "Amplitudes" -> OptionValue["Amplitudes"]];
      Return[$Failed];
    ];

    Block[{$degenerateIn = epsilonPlusDegeneracy[s, m, a, \[Omega], wp] =!= None},
      TRF[s, l, m, a, \[Omega], \[Lambda], \[Nu], bcs, norms, {wp, prec, acc}, Sequence@@subopts]]
  ];

  (* Call the chosen implementation, checking the MST solutions through their Wronskian and, when the
     check fails, recomputing everything with the working precision raised by wp and then 3 wp. With
     "WronskianCheck" -> Automatic the check (two summations of the MST series, for the value and
     derivative of each solution at one radius) runs only when there is a risk
     indicator: the amplitudes needed more than 40 digits of padding, or a retry after a non-numeric
     result. The failure it guards against, the coefficient recurrence yielding its wrong solution, needs a
     deep bump in the coefficients, which shows as padding of 100 digits and more, whereas ordinary modes
     need 10-20 digits at 32 digits of working precision and none at machine precision. *)
  (* The check used to be skipped at complex frequencies, where the identity appeared to be violated by an
     amount growing like |omega|^5 on the imaginary axis; that was the "Up" solution evaluated on the wrong
     side of a branch cut (see hypergeometricU in the MST package), and the identity now holds there to
     the working precision, within the cancellation between its two terms (mstWronskianError). *)
  (* a forced check needs the amplitudes B^inc and C^trans, which it cannot have when they are supplied or
     disabled: that combination is refused rather than the check silently not running *)
  (* the check is of the MST solutions; the machine-precision default makes its own accuracy estimate *)
  If[OptionValue["WronskianCheck"] === True && !MemberQ[{TeukolskyRadialMST, TeukolskyRadialAutomaticMachinePrecision}, TRF],
    Message[TeukolskyRadial::opti, {"WronskianCheck" -> True, Method -> OptionValue[Method]}];
    Return[$Failed];
  ];
  If[OptionValue["WronskianCheck"] === True && !MatchQ[OptionValue["Amplitudes"], Automatic|True],
    Message[TeukolskyRadial::opti, {"WronskianCheck" -> True, "Amplitudes" -> OptionValue["Amplitudes"]}];
    Return[$Failed];
  ];
  check = MatchQ[OptionValue["WronskianCheck"], True|Automatic] && TRF === TeukolskyRadialMST && MatchQ[OptionValue["Amplitudes"], Automatic|True];
  extra = Teukolsky`MST`MST`Private`modePadding[s, l, m, a, 2 \[Omega]];
  If[!check, Return[compute[extra, BCs]]];
  {ampPadding, ampRetried} = {0, False};
  res = compute[extra, {"In", "Up"}];
  If[res === $Failed, Return[$Failed]];
  (* l >= 20: the high-l modes whose coefficient recurrence can converge to its wrong solution do not always
     need much padding at moderate working precision (l = 36, m = 2, omega = 3 needs 11 digits at 24 and 40
     digits and is then wrong by O(1) with full tracked precision), so they are always checked *)
  (* Near-integer nu: the downward MST coefficients pass close to a pole there and can lose all their digits
     (l = 8, omega = -0.01 I: nu = 8 - 4e-5, a_n Infinity below n = -10 at 40 digits), so that the summation
     stops early and the result is wrong by 1e-3 with a tracked precision that claims full accuracy; such modes
     are checked as well. At real frequencies nu is near l for every small omega and these modes come out right,
     so the check (two summations of the series) is not added there. *)
  If[OptionValue["WronskianCheck"] === Automatic && !(ampPadding > 40 || ampRetried || extra > 0 || l >= 20 || (Im[\[Omega]] != 0 && NumericQ[\[Nu]] && Abs[\[Nu] - Round[Re[\[Nu]]]] < 10^-2)),
    Return[If[ListQ[BCs], KeyTake[res, BCs], res[BCs]]]];
  wpn = If[wp === MachinePrecision, $MachinePrecision, wp];
  tol = 10^(4 - wpn);
  e = mstWronskianError[res, s, a, \[Omega], wp];
  k = If[epsilonPlusDegeneracy[s, m, a, \[Omega], wp] === None, 0, 2];
  While[e > tol && k < 2,
    k++;
    Teukolsky`MST`MST`Private`setModePadding[s, l, m, a, 2 \[Omega], extra + wpn (2^k - 1)];
    res = compute[extra + wpn (2^k - 1), {"In", "Up"}];
    If[res === $Failed, Return[$Failed]];
    e = mstWronskianError[res, s, a, \[Omega], wp];
  ];
  If[e > tol, Message[TeukolskyRadial::acc, N[e, 2]]];
  If[ListQ[BCs], KeyTake[res, BCs], res[BCs]]
];


(* ::Section::Closed:: *)
(*TeukolskyRadialFunction*)


(* ::Subsection::Closed:: *)
(*Output format*)


(* ::Subsubsection::Closed:: *)
(*Icons*)


icons = <|
 "In" -> Graphics[{
         Line[{{0,1/2},{1/2,1},{1,1/2},{1/2,0},{0,1/2}}],
         Line[{{3/4,1/4},{1/2,1/2}}],
         {Arrowheads[0.2],Arrow[Line[{{1/2,1/2},{1/4,3/4}}]]},
         {Arrowheads[0.2],Arrow[Line[{{1/2,1/2},{3/4,3/4}}]]}},
         Background -> White,
         ImageSize -> Dynamic[{Automatic, 3.5 CurrentValue["FontCapHeight"]/AbsoluteCurrentValue[Magnification]}]],
 "Up" -> Graphics[{
         Line[{{0,1/2},{1/2,1},{1,1/2},{1/2,0},{0,1/2}}],
         Line[{{1/4,1/4},{1/2,1/2}}],
         {Arrowheads[0.2],Arrow[Line[{{1/2,1/2},{1/4,3/4}}]]},
         {Arrowheads[0.2],Arrow[Line[{{1/2,1/2},{3/4,3/4}}]]}},
         Background -> White,
         ImageSize -> Dynamic[{Automatic, 3.5 CurrentValue["FontCapHeight"]/AbsoluteCurrentValue[Magnification]}]]
|>;


(* ::Subsubsection::Closed:: *)
(*Formatting of TeukolskyRadialFunction*)


TeukolskyRadialFunction /:
 MakeBoxes[trf:TeukolskyRadialFunction[s_, l_, m_, a_, \[Omega]_, assoc_], form:(StandardForm|TraditionalForm)] :=
 Module[{summary, extended},
  summary = {Row[{BoxForm`SummaryItem[{"s: ", s}], "  ",
                  BoxForm`SummaryItem[{"l: ", l}], "  ",
                  BoxForm`SummaryItem[{"m: ", m}], "  ",
                  BoxForm`SummaryItem[{"a: ", a}], "  ",
                  BoxForm`SummaryItem[{"\[Omega]: ", \[Omega]}]}],
             BoxForm`SummaryItem[{"Domain: ", assoc["Domain"]}],
             BoxForm`SummaryItem[{"Boundary Conditions: " , assoc["BoundaryConditions"]}]};
  If[assoc["Method"] === {"Static"},
  extended = {BoxForm`SummaryItem[{"Eigenvalue: ", assoc["Eigenvalue"]}],
              BoxForm`SummaryItem[{"Amplitudes: ", assoc["Amplitudes"]}],
              BoxForm`SummaryItem[{"Method: ", First[assoc["Method"]]}]},
  extended = {BoxForm`SummaryItem[{"Eigenvalue: ", assoc["Eigenvalue"]}],
              BoxForm`SummaryItem[{"Renormalized angular momentum: ", assoc["RenormalizedAngularMomentum"]}],
              BoxForm`SummaryItem[{"Transmission Amplitude: ", assoc["Amplitudes", "Transmission"]}],
              BoxForm`SummaryItem[{"Incidence Amplitude: ", Lookup[assoc["Amplitudes"], "Incidence", Missing]}],
              BoxForm`SummaryItem[{"Reflection Amplitude: ", Lookup[assoc["Amplitudes"], "Reflection", Missing]}],
              BoxForm`SummaryItem[{"Method: ", First[assoc["Method"]]}],
              BoxForm`SummaryItem[{"Method options: ",Column[Rest[assoc["Method"]]]}]}];
  BoxForm`ArrangeSummaryBox[
    TeukolskyRadialFunction,
    trf,
    Lookup[icons, assoc["BoundaryConditions"], None],
    summary,
    extended,
    form]
];


(* ::Subsection::Closed:: *)
(*Accessing attributes*)


TeukolskyRadialFunction[s_, l_, m_, a_, \[Omega]_, assoc_][y_String] /; !MemberQ[{"RadialFunction"}, y] :=
  assoc[y];


Keys[m_TeukolskyRadialFunction] ^:= DeleteCases[Join[Keys[m[[-1]]], {}], "RadialFunction"];


(* ::Subsection::Closed:: *)
(*Numerical evaluation*)


SetAttributes[TeukolskyRadialFunction, {NHoldAll}];


outsideDomainQ[r_, rmin_, rmax_] := Min[r]<rmin || Max[r]>rmax;


TeukolskyRadialFunction[s_, l_, m_, a_, \[Omega]_, assoc_][r:(_?NumericQ|{_?NumericQ..})] :=
 Module[{rmin, rmax},
  {rmin, rmax} = assoc["Domain"];
  If[outsideDomainQ[r, rmin, rmax],
    Message[TeukolskyRadialFunction::dmval, #]& /@ Select[Flatten[{r}], outsideDomainQ[#, rmin, rmax]&];
    Return[Indeterminate];
  ];
  Quiet[assoc["RadialFunction"][r], InterpolatingFunction::dmval]
 ];


(* R[r, n] is the n-th derivative; R[r, {0, 1}] the value and the first derivative, which for an MST
   solution come from a single summation of the series (the coefficients and the hypergeometric
   functions are shared), at about half the cost of the two separate evaluations *)
TeukolskyRadialFunction[s_, l_, m_, a_, \[Omega]_, assoc_][r:(_?NumericQ|{_?NumericQ..}), n_Integer?NonNegative] :=
  Derivative[n][TeukolskyRadialFunction[s, l, m, a, \[Omega], assoc]][r];

TeukolskyRadialFunction[s_, l_, m_, a_, \[Omega]_, assoc_][r_?NumericQ, {0, 1}] :=
 Module[{rmin, rmax, f = assoc["RadialFunction"]},
  {rmin, rmax} = assoc["Domain"];
  If[outsideDomainQ[r, rmin, rmax],
    Message[TeukolskyRadialFunction::dmval, r];
    Return[Indeterminate];
  ];
  If[MatchQ[f, _Teukolsky`MST`MST`Private`MSTRadialIn | _Teukolsky`MST`MST`Private`MSTRadialUp],
    f[r, {0, 1}],
    Quiet[{f[r], f'[r]}, InterpolatingFunction::dmval]
  ]
 ];

TeukolskyRadialFunction[s_, l_, m_, a_, \[Omega]_, assoc_][r:(_?NumericQ|{_?NumericQ..}), ns:{_Integer?NonNegative..}] :=
  TeukolskyRadialFunction[s, l, m, a, \[Omega], assoc][r, #] & /@ ns;


Derivative[n:1][TeukolskyRadialFunction[s_, l_, m_, a_, \[Omega]_, assoc_]][r:(_?NumericQ|{_?NumericQ..})] :=
 Module[{rmin, rmax},
  {rmin, rmax} = assoc["Domain"];
  If[outsideDomainQ[r, rmin, rmax],
    Message[TeukolskyRadialFunction::dmval, #]& /@ Select[Flatten[{r}], outsideDomainQ[#, rmin, rmax]&];
    Return[Indeterminate];
  ];
  Quiet[Derivative[n][assoc["RadialFunction"]][r], InterpolatingFunction::dmval]
 ];


Derivative[n_Integer/;n>1][trf:(TeukolskyRadialFunction[s_, l_, m_, a_, \[Omega]_, assoc_])][r0:(_?NumericQ|{_?NumericQ..})] :=
 Module[{Rderivs, R, r, i, res},
  Rderivs = D[R[r_], {r_, i_}] :> D[(-(-trf["Eigenvalue"] + 2 I r s 2 \[Omega] + (-2 I (-1 + r) s (-a m + (a^2 + r^2) \[Omega]) + (-a m + (a^2 + r^2) \[Omega])^2)/(a^2 - 2 r + r^2)) R[r] - (-2 + 2 r) (1 + s) R'[r])/(a^2 - 2 r + r^2), {r, i - 2}] /; i >= 2;
  Do[Derivative[i][R][r] = Collect[D[Derivative[i - 1][R][r], r] /. Rderivs, {R'[r], R[r]}, Simplify];, {i, 2, n}];
  res = Derivative[n][R][r] /. {
    R'[r] -> trf'[r0],
    R[r] -> trf[r0], r -> r0};
  Clear[Rderivs, i];
  Remove[R, r];
  res
];


(* ::Section::Closed:: *)
(*End Package*)


(* ::Subsection::Closed:: *)
(*Protect symbols*)


SetAttributes[{TeukolskyRadial, TeukolskyRadialFunction}, {Protected, ReadProtected}];


(* ::Subsection::Closed:: *)
(*End*)


End[];
EndPackage[];
