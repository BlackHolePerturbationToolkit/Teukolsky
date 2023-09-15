{% include head.html %}

<p>
 <h1 style="display:inline">Teukolsky</h1> <span style="float:right;"><a href="https://bhptoolkit.org/mathematica-install.html" class = "code_btn">Install this package!</a></span>
</p>

A Mathematica package for computing solutions to the Teukolsky equation. Note this package depends upon the [SpinWeightedSpheroidalHarmonics](https://bhptoolkit.org/SpinWeightedSpheroidalHarmonics/) and the [KerrGeodesics](https://bhptoolkit.org/KerrGeodesics/) package to run.

Explicitly the package computes solutions to:

$$\Delta^{-s} \dfrac{d}{dr} \bigg[\Delta^{s+1}\dfrac{d R}{dr}\bigg] + \bigg[\frac{K^2 - 2 i s (r-M)K}{\Delta} + 4 i s \omega r - \lambda \bigg]R = \mathcal{T} \nonumber $$
 
where

$\Delta = r^2 - 2Mr + a^2$  
$K=(r^2 + a^2)\omega - a m$  
$s$ is the spin-weight of the perturbing field  
$\lambda$ is the spin-weighted spheroidal eigenvalue  
$\omega$ is the mode frequency  
$\mathcal{T}$ is the source

Currently the source has been implemented for a point particle moving along a generic bound orbit in Kerr spacetime for perturbations with spin-weight $s={0,\pm 1, \pm 2}$. As an example, the gravitational wave flux for the $l=2,m=2$ mode for a circular, equatorial orbit is easily computed using:  
```Mathematica
With[{a = 0.9, p = 10.0, e=0, x=1, s = -2, l = 2, m = 2},
  orbit = KerrGeoOrbit[a, p, e, x];
  ψ4 = TeukolskyPointParticleMode[s, l, m, orbit];
  ψ4["Fluxes"]
]
```  
This returns an association with the results:  
```Mathematica
<|"Energy" -> <|"ℐ" -> 0.000022273, "ℋ" -> -5.9836*10^-8|>,
  "AngularMomentum" -> <|"ℐ" -> 0.00072438, "ℋ" -> -1.94603*10^-6|>|>
```

## Homogeneous solutions

The homogeneous solutions are also easily computed. They can be extracted from the `ψ4` object above using `R = ψ4["RadialFunctions"]`. This returns a pair `TeukolskyRadialFunction` objects which can be evaluated at a given radius, i.e., `R["In"][10.]`. The homogeneous solutions can also be computed directly via the `TeukolskyRadial[s, l, m, a, ω]` function.

## Renormalized angular momentum

Under the hood the Teukolsky package uses the MST method for computing homogeneous solutions. A key part of the MST method is the calculation of the renormalized angular momentum, $\nu$. This can be computed directly via
```
 ν = RenormalizedAngularMomentum[s, l, m, a, ω]
```

### Further examples

See the Mathematica Documentation Centre for a tutorial and documentation on individual functions.

## Citing

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7037850.svg)](https://doi.org/10.5281/zenodo.7037850)

In addition to acknowledging the Black Hole Perturbation Toolkit as suggested on the [front page](https://bhptoolkit.org) we also recommend citing the specific package version you use via the citation information on the package’s Zenodo page linked from the above DOI.

## Authors and contributors

Barry Wardell, Niels Warburton, Marc Casals, Adrian Ottewill, Chris Kavanagh, Leanne Durkan, Ben Leather, Theo Torres
