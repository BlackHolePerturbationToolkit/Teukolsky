(* ::Package:: *)

(* ::Input:: *)
(*SetOptions[EvaluationNotebook[],StyleDefinitions->$UserBaseDirectory<>"/SystemFiles/FrontEnd/StyleSheets/maTHEMEatica.nb"]*)


(* ::Section:: *)
(*Beginning*)


(* ::Subsection:: *)
(*Beginning*)


BeginPackage["Teukolsky`PN`Tools`",{"Teukolsky`","Teukolsky`PN`"}]


(* ::Subsection:: *)
(*Unprotecting*)


ClearAttributes[{SeriesTake,SeriesMap,SeriesCoefficientList,SeriesMinOrder,SeriesMaxOrder,SeriesLength,SeriesCollect,SeriesExpand,SeriesTerms,IgnoreSeriesParameter,ChangeSeriesParameter,PowerCounting,StraightenSeries,SeriesPlusSimplify,DropZeroSeries,InactiveSeriesPrefactor}, {Protected, ReadProtected}];


ClearAttributes[{Scalings, RemovePN}, {Protected, ReadProtected}];


ClearAttributes[{ExpandLog, ExpandGamma,ExpandPolyGamma,PochhammerToGamma,GammaToPochhammer,ExpandDiracDelta,CollectDerivatives}, {Protected, ReadProtected}];


ClearAttributes[{TeukolskyPointParticleSource}, {Protected, ReadProtected}];


ClearAttributes[{Paint,CowboyConjugate,ChangeContext,ChooseSide}, {Protected, ReadProtected}];


ClearAttributes[{AngularTeukolskyEquation,RadialTeukolskyEquation,RadialTeukolskyEquationPN}, {Protected, ReadProtected}];


(* ::Section:: *)
(*Public*)


(* ::Subsection:: *)
(*General Tools for Series*)


SeriesTake::usage="SeriesTake[series, n] takes the first n terms of series"
SeriesMap::usage="SeriesMap[function,series] maps f onto the coefficients of series"
SeriesMinOrder::usage="SeriesMinOrder[series] gives the leading order of series"
SeriesMaxOrder::usage="SeriesMaxOrder[series] gives the first suppressed order of series"
SeriesLength::usage="SeriesLenght[series] gives the number of terms in series"
SeriesCollect::usage="SeriesCollect[expr, var, func] works like Collect but applied to each order individually. Crucially, unlike Collect it keeps the SeriesData structure."
SeriesExpand::usage="SeriesExpand[expr] works like Expand but applied to each order individually."
SeriesTerms::usage="SeriesTerms[series, {x, x0, n}] works exactly like Series, with the difference that n gives the desired number of terms instead of a maximum order"
IgnoreSeriesParameter::usage="IgnoreSeriesParameter[series,x] sets all occurences of the expansion parameter in the series coefficients to x. If no value is entered x defaults to 1."
ChangeSeriesParameter::usage="ChangeSeriesParameter[series,expr] changes the expansion parameter in series to be expr." 
PowerCounting::usage="PowerCounting[series,symbol] replaces the expansion parameter in series with symbol. Unlike ChangeSeriesParameter it keeps the original expansion parameter as a constant in each coefficient."
StraightenSeries::usage="StraightenSeries[expr] straightens out SeriesData objects with redundant denominator arguments, i.e., it removes counting in powers of roots, if all respective coefficients are zero."
DropZeroSeries::usage="DropZeroSeries sets the 0 Series O[x\!\(\*SuperscriptBox[\(]\), \(n\)]\) to 0 without Normaling the entire expressions."
SeriesPlusSimplify::usage="SeriesPlusSimplify[expr,assumptions] simplifies sums of SeriesData objects under assumptions."
SeriesCoefficientList::usage="SeriesCoefficientList[series] returns the Series coefficients as a list."
InactiveSeriesPrefactor::usage="InactiveSeriesPrefactor[series] pulls out the leading order of a series."


(* ::Subsection:: *)
(*Tools for PN Scalings*)


Scalings::usage="Scalings[expr,params,var] applies the given scalings params with power counting parameter var to expr."
(*PNScalings::usage="Same as Scalings but with different input. Just here to not break my older code but you should use Scalings instead"*)
RemovePN::usage="PNScalings[expr,var] takes the Normal[] and sets var to 1"
(*Zero::usage="Zero[expr,vars] sets all vars in expr to 0"
One::usage="One[expr,vars] sets all vars in expr to 1"*)


(* ::Subsection:: *)
(*Tools for Logs, Gammas, and PolyGammas*)


ExpandLog::usage="ExpandLog[expr,assumptions] expands all logaritms in expr under assumptions. Crucially unlike PowerExpand it does not make unprompted assumptions."
ExpandGamma::usage="ExpandGamma[expr] factors out all Integer facors out of the Gammas in expr. E.g. Gamma[x+1]->x Gamma[x]"
ExpandPolyGamma::usage="ExpandPolyGamma[expr] factors out all Integer facors out of the PolyGammas in expr. E.g. PolyGamma[x+1]->\!\(\*FractionBox[\(1\), \(x\)]\) PolyGamma[x]"
PochhammerToGamma::usage="PochhammerToGamma[expr] replaces all Pochhammer in expr with the respecive Gamma."
GammaToPochhammer::usage="GammaToPochhammer[expr,n] replaces all Gamma in expr that contain n with the respective Pochhammer[__,n]"



(* ::Subsection:: *)
(*Tools for DiracDelta *)


ExpandDiracDelta::usage="ExpandDiracDelta[expr,r] applies identities for Dirac deltas and it's derivatives to expr."


(* ::Subsection:: *)
(*Tools  for SpinWeightedSpheroidalHarmonics *)


(*ExpandSpheroidals::usage="ExpandSpheroidal[expr,{param,order}] returns a all SpinWeightedSpheroidalHarmonicS in expr have been Series expanded around param->0 to order."*)


(* ::Subsection:: *)
(*Teukolsky Equation*)


AngularTeukolskyEquation::usage="AngularTeukolskyEquation[s,l,m,\[Gamma],\[Theta],\[Phi]] returns the angular Teukolsky equation equation. It is solved by the spin weighted spheroidal harmonics"


RadialTeukolskyEquation::usage="RadialTeukolskyEquation[\[ScriptS],\[ScriptL],\[ScriptM],a,\[Omega],R[r]] gives the radial Teukolsky equation."
RadialTeukolskyEquationPN::usage="RadialTeukolskyEquationPN[\[ScriptS],\[ScriptL],\[ScriptM],a,\[Omega],R[r],{\[Eta],order}] gives the PN expanded radial Teukolsky equation."


(* ::Subsection:: *)
(*Misc*)


CollectDerivatives::usage="CollectDerivatives[expr,f] works exactly like Collect[] but also collects for derivatives of f."
Paint::usage="Paint[expr,var] paints all occurences of var in expr Red."
CowboyConjugate::usage="\"Shoot first, ask questions later\". CowboyConjugate[expr] performs the complex conjugate by assuming everything but \[ImaginaryI] is real."
ChangeContext::usage="ChangeContext[expr,context1,context2] is a debugging tool that allows to change the context of all symbols in expr"
ChooseSide::usage="ChooseSide[expr,assumptions] Simplifies all HeavisideTheta and DiracDelta according to assumptions."


(* ::Subsection:: *)
(*Source*)


TeukolskyPointParticleSource::usage="TeukolskyPointParticleSource[\[ScriptS],\[ScriptL],\[ScriptM],orbit][r] gives an analytical expression for the Teukolsky point particle source for a given {\[ScriptS],\[ScriptL],\[ScriptM]} mode. orbit needs to be a KerrGeoOrbit object "


(* ::Section:: *)
(*Private*)


Begin["Private`"]


(* ::Subsection:: *)
(*MST Coefficients*)


(*MSTCoefficients=Teukolsky`PN`Private`MSTCoefficients*)


(* ::Subsection:: *)
(*General  Tools  for  Series*)


SeriesTake=Teukolsky`PN`Private`SeriesTake
SeriesMinOrder=Teukolsky`PN`Private`SeriesMinOrder
SeriesMaxOrder=Teukolsky`PN`Private`SeriesMaxOrder
SeriesLength=Teukolsky`PN`Private`SeriesLength
SeriesCollect=Teukolsky`PN`Private`SeriesCollect
SeriesExpand=Teukolsky`PN`Private`SeriesExpand
SeriesTerms=Teukolsky`PN`Private`SeriesTerms
IgnoreSeriesParameter=Teukolsky`PN`Private`IgnoreExpansionParameter
ChangeSeriesParameter=Teukolsky`PN`Private`ChangeSeriesParameter
PowerCounting=Teukolsky`PN`Private`PowerCounting
StraightenSeries=Teukolsky`PN`Private`StraightenSeries
SeriesMap=Teukolsky`PN`Private`SeriesMap
DropZeroSeries=Teukolsky`PN`Private`DropZeroSeries
SeriesPlusSimplify=Teukolsky`PN`Private`SeriesPlusSimplify
SeriesCoefficientList=Teukolsky`PN`Private`SeriesCoefficientList
InactiveSeriesPrefactor=Teukolsky`PN`Private`InactiveSeriesPrefactor
InactiveSeriesPrefactor


(* ::Subsection::Closed:: *)
(*Tools for PN Scalings*)


(*PNScalings=Teukolsky`PN`Private`PNScalings*)
Scalings=Teukolsky`PN`Private`Scalings
RemovePN=Teukolsky`PN`Private`RemovePN
(*Zero=Teukolsky`PN`Private`Zero
One=Teukolsky`PN`Private`One*)


(* ::Subsection::Closed:: *)
(*Tools for Logs, Gammas, and PolyGammas*)


ExpandLog=Teukolsky`PN`Private`ExpandLog
ExpandGamma=Teukolsky`PN`Private`ExpandGamma
ExpandPolyGamma=Teukolsky`PN`Private`ExpandPolyGamma
PochhammerToGamma=Teukolsky`PN`Private`PochhammerToGamma
GammaToPochhammer=Teukolsky`PN`Private`GammaToPochhammer


(* ::Subsection::Closed:: *)
(*Tools for DiracDelta *)


ExpandDiracDelta=Teukolsky`PN`Private`ExpandDiracDelta


(* ::Subsection::Closed:: *)
(*Tools  for SpinWeightedSpheroidalHarmonics *)


(*ExpandSpheroidals=Teukolsky`PN`Private`ExpandSpheroidals*)


(* ::Subsection::Closed:: *)
(*Misc*)


CollectDerivatives=Teukolsky`PN`Private`CollectDerivatives
Paint=Teukolsky`PN`Private`Paint
CowboyConjugate=Teukolsky`PN`Private`CowboyConjugate
ChangeContext=Teukolsky`PN`Private`ChangeContext
ChooseSide=Teukolsky`PN`Private`ChooseSide


(* ::Subsection:: *)
(*Source*)


TeukolskyPointParticleSource=Teukolsky`PN`Private`TeukolskyPointParticleSource


(* ::Subsection:: *)
(*Teukolsky Equation*)


AngularTeukolskyEquation=Teukolsky`PN`Private`AngularTeukolskyEquation


RadialTeukolskyEquation=Teukolsky`PN`Private`RadialTeukolskyEquation
RadialTeukolskyEquationPN=Teukolsky`PN`Private`RadialTeukolskyEquationPN


(* ::Section:: *)
(*Ending Package*)


(* ::Subsection:: *)
(*Protecting*)


SetAttributes[{SeriesTake,SeriesMap,SeriesCoefficientList,SeriesMinOrder,SeriesMaxOrder,SeriesLength,SeriesCollect,SeriesExpand,SeriesTerms,IgnoreSeriesParameter,ChangeSeriesParameter,PowerCounting,StraightenSeries,SeriesPlusSimplify,DropZeroSeries,InactiveSeriesPrefactor}, {Protected, ReadProtected}];


SetAttributes[{Scalings, RemovePN}, {Protected, ReadProtected}];


SetAttributes[{ExpandLog, ExpandGamma,ExpandPolyGamma,PochhammerToGamma,GammaToPochhammer,ExpandDiracDelta,CollectDerivatives}, {Protected, ReadProtected}];


SetAttributes[{TeukolskyPointParticleSource}, {Protected, ReadProtected}];


SetAttributes[{Paint,CowboyConjugate,ChangeContext}, {Protected, ReadProtected}];


SetAttributes[{AngularTeukolskyEquation,RadialTeukolskyEquation,RadialTeukolskyEquationPN}, {Protected, ReadProtected}];


(* ::Subsection:: *)
(*Ending*)


End[]
EndPackage[]
