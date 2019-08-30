{% include head.html %}

# Teukolsky

A Mathematica package for computing solutions to the Teukolsky equation. Note this package depends upon the [SpinWeightedSpheroidalHarmonics](https://bhptoolkit.org/SpinWeightedSpheroidalHarmonics/) and the [KerrGeodesics](https://bhptoolkit.org/KerrGeodesics/) package to run.

Explicitly the package computes solutions to:

$\Delta^{-s} \dfrac{d}{dr} \bigg[\Delta^{s+1}\dfrac{d R}{dr}\bigg] + \bigg[\frac{K^2 - 2 i s (r-M)K}{\Delta} + 4 i s \omega r - \lambda \bigg]R = \text{source}$
 
where

$K=$
$s$ is the spin-weight of the perturbing field.
$\lambda$ is the spin-weighted spheroidal eigenvalue.
$\omega$ is the mode frequency
