# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed
 - TeukolskyRadial with exact arguments (e.g. a = 6/10, omega = 1/2) evaluates them at the given WorkingPrecision instead of returning unevaluated, and fails with the new TeukolskyRadial::exact message when no WorkingPrecision is given. N[RenormalizedAngularMomentum[s, l, m, a, omega], p] now works with exact arguments.
 - At the superradiant bound frequency omega = m Omega_H the radial functions were silently mis-normalised (the unscaled transmission amplitude was taken to be 1). The asymptotic amplitudes there are now obtained as the limit from neighbouring frequencies; the ones that genuinely diverge (the "Up" horizon coefficients, and for s >= 1 the "In" amplitudes) are Indeterminate and the new TeukolskyRadial::superradiant message lists them. Frequency grids commensurate with Omega_H hit this point (a = 3/5 gives Omega_H = 1/6).
 - Public symbols (TeukolskyRadial, TeukolskyRadialFunction, TeukolskyMode, TeukolskyPointParticleMode, RenormalizedAngularMomentum) now live in the Teukolsky` context rather than in sub-contexts, so that packages depending on Teukolsky` (BeginPackage["X`", {"Teukolsky`"}]) see them instead of silently creating private symbols of the same name.

### Changed
 - Improved accuracy of the radial functions and of the source integrals:
   - The default PrecisionGoal is now WorkingPrecision - 2 for all methods (previously WorkingPrecision / 2, which limited machine-precision solutions to ~8 digits and 32-digit solutions to ~16 digits).
   - With Method -> Automatic at machine precision, TeukolskyRadial now returns numerically integrated solutions with goals of WorkingPrecision - 2 whose boundary data are precision-padded MST solutions (near the horizon for "In", one unit beyond the outermost requested radius for "Up", cached between evaluations), accurate to ~1e-11 - 1e-15 and cheap to evaluate on many radii; the accuracy of the pair is checked through the Wronskian, which must equal 2 I omega B^inc C^trans, and the new TeukolskyRadial::acc message reports a poor estimate.
   - The MST series (renormalized angular momentum, asymptotic amplitudes, radial functions) are now evaluated with the working precision padded by the number of digits lost to cancellation in the series, measured from the tracked precision of a first evaluation, so that the results carry the precision of the input; at machine precision the series are summed in arbitrary precision. Previously the MST radial functions lost roughly 0.8 omega (r - r+) digits ("In") and an omega-dependent number of digits ("Up"), which at machine precision made them unusable for omega M >~ 0.5.
   - At large omega (r - r+) the MST "In" radial function is now evaluated as a combination of the outgoing ("Up") and ingoing Coulomb-type series, with connection coefficients given analytically by the coefficients K_nu and A_+/-^nu of Sasaki & Tagoshi (Eqs. (157), (158), (165), (166)), so that its accuracy no longer degrades with radius. The hypergeometric series in 1/(1-x) of Sasaki & Tagoshi Eq. (138) is available as an alternative large-radius representation (Teukolsky`MST`MST`Private`$MSTInLargeRadiusRepresentation = "Hypergeometric").
   - For a finite "Domain", the NumericalIntegration method now takes the "Up" boundary data from the MST solution at the outer edge of the domain rather than at r = 1000 (the MST "Up" series is accurate at any radius beyond about r+ + 1), and the "In" boundary data from the MST solution at the inner edge rather than near the horizon (the precision-padded "In" solution, evaluated from the Coulomb-type series at large radius, is accurate at any radius).
   - For non-circular orbits, TeukolskyPointParticleMode now seeds the numerical integration over the radial libration region with the values of the global solutions (same choice of radii), with goals of WorkingPrecision - 2, rather than with MST data at large radius and goals of WorkingPrecision / 2 (new "BoundaryData" option of the NumericalIntegration method).


## [1.1.1] - 2025-06-26

### Fixed
- Corrected PN calculation of inhomogeneous mode amplitudes
- Improved performance of PN calculation

## [1.1.0] - 2025-01-29

### Added
 - Support for complex frequencies.
 - Post-Newtonian series expansions.
 - Support for computing homogeneous solutions using confluent Heun functions.
 - Tutorial on self-force calculations.

### Fixed
 - Various numerical corner cases addressed


## [1.0.0] - 2022-09-17

### Added
 - Support for additional point particle orbits:
   - Generic orbits in Kerr spacetime for s=-2, 0, +2.
   - Circular orbits in Kerr spacetime for s=+1.
 - Improvements to TeukolskyRadial and TeukolskyRadialFunction:
   - Support for computing amplitudes is now available for all Methods.
   - Options have been added to not compute amplitudes, renormalized angular momentum, and the eigenvalue.
   - Options have been added to pass in amplitudes, renormalized angular momentum, and the eigenvalue.
   - Performance improvements by avoiding repeatedly computing amplitudes, nu, and the eigenvalue.
   - Second and higher derivatives are now more accurately and quickly computed by using the field equation.
   - Unscaled amplitudes with non-unit transmission coefficient are now availalbe in a TeukolskyRadialFunction.
   - Keys now works with a TeukolskyRadialFunction.
   - Improvements to NumericalIntegration method:
     - Asymptotic amplitudes are now available.
     - Domain -> All now works and is the default. 
     - Significant improvement in speed when creating a TeukolskyRadialFunction.
     - Private symbols are no longer exposed through RadialFunction.
     - WorkingPrecision, PrecisionGoal and AccuracyGoal options are respected in boundary conditions.
     - High precision is used for setting MST boundary conditions when working at machine precision.
 - Improvements to TeukolskyPointParticleMode and TeukolskyMode:
   - Support in TeukolskyMode for numerical evaluation of the inhomogeneous solution outside the source region.
   - Support in TeukolskyMode for numerical evaluation with extended homogeneous solutions.
   - Keys now works with a TeukolskyMode.
   - TeukolskyPointParticleMode now has a "Domain" option.
 

### Changed
 - Default Method for TeukolskyRadial changed to NumericalIntegration when working at machine precision.
 - Default for AccuracyGoal in TeukolskyRadial is now Infinity.
 - With NumericalIntegration the "up" boundary condition is now applied at r=1000.
 - Implementation of sources using TeukolskySource has been simplified.
 - The n and k arguments to TeukolskyPointParticleMode are now optional for special (circular, spherical, eccentric) orbits.

### Fixed
 - Resolved some memory leaks.
 - Fixed problem with evaluating derivatives using numerical integration.
 - Fixed Flux calculation when omega = 0.
 - Fixed problem with inclination <0.


## [0.3.0] - 2020-09-14

### Added
 - Support for solving for a point particle with electric charge on a circular orbit in Kerr.
 - Support for computing solutions using numerical integration on a hyperbolical slice.
 - Support for computing "In" solutions using Mathematica's HeunC function (available since version 12.1).
 - All asymptotic amplitudes (incidence, transmission and reflection) are now computed and available in a TeukolskyRadialFunction.

### Fixed
 - Fixed several problems with static modes:
   - For s = 0, m=0 the static "Up" solutions incorrectly evaluated to 0.
   - For arbitrary s the static "Up" solutions for m != 0 but a = 0 returned "Indeterminate."
   - The static "Up" solutions when m != 0 and a !=0 experienced large cancellations.
   - The "In" solutions were not consistent with their values for small non-zero omega.
 - Fixed several memory leaks.
 - Fixed a problem where loading the package with Needs would generate an error message.


## [0.2.0] - 2020-05-24

### Added
 - Working precision now tracks the precision of both a and omega.
  
## [0.1.0] - 2020-05-23
 - Initial release.
