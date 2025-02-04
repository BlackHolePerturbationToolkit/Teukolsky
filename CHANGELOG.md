# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
