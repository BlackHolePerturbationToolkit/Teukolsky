{% include head.html %}

# Teukolsky

A Mathematica package for computing solutions to the Teukolsky equation. Note this package depends upon the [SpinWeightedSpheroidalHarmonics](https://bhptoolkit.org/SpinWeightedSpheroidalHarmonics/) and the [KerrGeodesics](https://bhptoolkit.org/KerrGeodesics/) package to run.

Explicitly the package computes solutions to:

$\Delta^{-s} \dfrac{d}{dr} \bigg[\Delta^{s+1}\dfrac{d R}{dr}\bigg] + \bigg[\frac{K^2 - 2 i s (r-M)K}{\Delta} + 4 i s \omega r - \lambda \bigg]R = \mathcal{T}$
 
where

$\Delta = r^2 - 2Mr + a^2$  
$K=(r^2 + a^2)\omega - a m$  
$s$ is the spin-weight of the perturbing field 
$\lambda$ is the spin-weighted spheroidal eigenvalue 
$\omega$ is the mode frequency  
$\mathcal{T}$ is the source

Currently the source has been implemented for a point particle moving along a circular orbit in Kerr spacetime. As an example, the flux in this case for the $l=2,m=2$ mode is easily computed using:  
```
a = 0.9`32;
r0 = 10.`32;

orbit = KerrGeoOrbit[a, r0, 0, 1];

s = -2; l = 2; m = 2; n = 0; k = 0;
mode = TeukolskyPointParticleMode[s, l, m, n, k, orbit];

mode["Fluxes"]
```  
This returns an association with the results:  
```
<|"FluxInf" -> 0.000044546001102360994, "FluxHor" -> -1.1967358426766726*10^-7, "FluxTotal" -> 0.000044426327518093327|>
```
Note the high precision of the input values for $a$ and $r_0$. Currently this is often a requirement to get an accurate result.


