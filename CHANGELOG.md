# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed
 - At complex frequencies RenormalizedAngularMomentum now returns the representative of nu continuous with nu = l at omega = 0 (l - ArcCos[Cos[2 Pi nu]]/(2 Pi)), as for real frequencies, instead of the principal value near 0, on which the amplitude formulae hit the exact pole of Pochhammer[2 nu + 2, n] (Indeterminate amplitudes at small |omega| on the imaginary axis). The Regge-Wheeler "Up" reflection amplitude, obtained from a symmetry valid for real frequencies only, is Indeterminate at complex frequencies (new message) instead of wrong.
 - TeukolskyRadial now checks the MST solutions through their Wronskian, which must equal 2 I omega B^inc C^trans, at every working precision, and when the check fails recomputes everything with the working precision raised by wp and then 3 wp before reporting the estimated accuracy (TeukolskyRadial::acc). This catches modes such as l = 36, m = 2, omega = 3 at a = 3/5, for which the recurrence for the MST coefficients silently yields the wrong solution below about 300 digits with a tracked precision that claims the working precision; at 80 digits that mode was wrong by 60 orders of magnitude and is now right. The check costs two summations of the MST series (value and derivative of the "In" and of the "Up" solution at one radius), so by default ("WronskianCheck" -> Automatic) it runs only when the amplitudes needed more than 40 digits of padding or a retry, which the modes at risk do (over 100 digits) and ordinary modes do not (10-20 digits at 32 digits of working precision, none at machine precision); "WronskianCheck" -> True forces it and False disables it. The error is measured relative to the larger of the exact Wronskian and the size of its two terms, so that at a complex frequency, where the terms are of order Exp[2 |Im omega| r*] times the Wronskian, the cancellation between them is not mistaken for an error of the solutions.
 - RenormalizedAngularMomentum (Method "Monodromy") no longer loops until the kernel dies at the degeneracies 2 I epsilon = n (omega = -I n/4 M for any a and s), where mu1 - mu2 = 2 I epsilon - 2 s is an integer and Gamma[mu1 - mu2] or Gamma[mu2 - mu1] in the monodromy formula has a pole: the products Gamma[mu1 - mu2] Pochhammer[mu1 - mu2, k] are now combined into Gamma[mu1 - mu2 + k], which is finite there (the Pochhammer symbol carries the cancelling zero), so nu is evaluated directly and to the same precision as at any other frequency (an intermediate version extrapolated it from neighbouring frequencies, with an O(h^2) error of 1e-24 at 32 digits, and reported that with a RenormalizedAngularMomentum::degenerate message, now gone); the memoised recurrences are cleared on every path, evaluated iteratively (no $RecursionLimit at large nmax), and bounded (RenormalizedAngularMomentum::conv). On the imaginary axis the degeneracies 2 I epsilon_+ = n of the MST amplitudes are handled like the superradiant bound frequency (TeukolskyRadial::degenerate).
 - High-l modes (e.g. l = m = 40 for a circular orbit at a = 3/5) came back wrong, or Indeterminate, with a precision estimate that claimed many correct digits: the Monodromy renormalized angular momentum carried far fewer digits than the working precision and the MST series, which lose about 100 digits to cancellation for such a mode, amplify an error in nu by tens of digits, which SetPrecision hid. The eigenvalue and nu are now recomputed at the padded working precision of every MST evaluation, a non-numeric result is retried at twice the precision, and a result still short of the requested precision is reported by the new TeukolskyRadial::prec and TeukolskyRadialFunction::prec messages. For a mode deep under its potential barrier the Coulomb-type representation of the "In" solution, which cancels the barrier suppression, is abandoned for the hypergeometric series.
 - TeukolskyRadial with exact arguments (e.g. a = 6/10, omega = 1/2) evaluates them at the given WorkingPrecision instead of returning unevaluated, and fails with the new TeukolskyRadial::exact message when no WorkingPrecision is given. N[RenormalizedAngularMomentum[s, l, m, a, omega], p] now works with exact arguments.
 - At the superradiant bound frequency omega = m Omega_H the radial functions were silently mis-normalised (the unscaled transmission amplitude was taken to be 1). The asymptotic amplitudes there are now evaluated directly from a pole-free form of the MST formulae: the singular part is two explicit factors, 1/(Gamma[1 - s - x] Sin[Pi x]) in the "Up" reflection and 1/(Gamma[1 + s + x] Sin[Pi x]) in the "Up" incidence with x = 2 I epsilon_+, whose limits at integer x follow from Gamma[x] Gamma[1 - x] = Pi/Sin[Pi x] (an intermediate version extrapolated all amplitudes from four neighbouring frequencies, which cost 1.5 s per mode). The "Up" coefficient along the horizon basis function that has ceased to be independent (the larger-exponent one) genuinely diverges and is Indeterminate, and the new TeukolskyRadial::superradiant message lists the Indeterminate amplitudes. Frequency grids commensurate with Omega_H hit this point (a = 3/5 gives Omega_H = 1/6). For s >= 1 the "In" solution has no unit-transmission limit there (it used to hang evaluating the series); see the "In" normalisation change below.
 - Public symbols (TeukolskyRadial, TeukolskyRadialFunction, TeukolskyMode, TeukolskyPointParticleMode, RenormalizedAngularMomentum) now live in the Teukolsky` context rather than in sub-contexts, so that packages depending on Teukolsky` (BeginPackage["X`", {"Teukolsky`"}]) see them instead of silently creating private symbols of the same name.

 - The MST "Up" solution at a frequency with Re omega < 0 was wrong by O(1) (it violated the symmetry R[s, l, m, a, omega] = Conjugate[R[s, l, -m, a, -Conjugate[omega]]] of the radial equation by up to a factor 8, and the Coulomb-type "In" representation at large radius, the "Up" reflection amplitude and the Wronskian identity 2 I omega B^inc C^trans with it), and on the negative imaginary axis by a smaller amount growing with |omega| (1e-6 at |omega| = 0.3, O(1) at 1.7: the "violation of the Wronskian identity growing like |omega|^5" that had been attributed to the amplitude formulae, and the failure of the Coulomb-type "In" representation at complex frequencies, which had been disabled there). The Coulomb-type series of Sasaki and Tagoshi are derived for Re epsilon > 0: for Re omega < 0 the powers of -2 I zhat and of epsilon in them and in the amplitude formulae cross the cuts of their principal branches, and on the negative imaginary axis the argument -2 I zhat of the Tricomi functions is negative real, on their cut, where Mathematica's principal value is the limit from the wrong side. The MST radial functions and amplitudes at Re omega < 0 are now the complex conjugates of those at (-m, -Conjugate[omega]), and on the imaginary axis the Tricomi functions and the prefactor of the incoming Coulomb series are evaluated on the side of the cut continuous with Re omega > 0 (the Tricomi function just below its cut, with the precision of its arguments raised so that the shift is representable: the connection formula between the two sides, DLMF 13.2.12, loses up to fifteen digits to cancellation). The Coulomb-type representation of the "In" solution is used again at large radius for complex frequencies, the Wronskian check of the MST solutions is applied there too, and the symmetry, the Wronskian identity and the continuity of "Up" across Re omega = 0 are verified to 1e-28 at 32 digits by the new ComplexFrequency tests.

 - Tests/AllTests.wls loads only the repository it lives in (and the dependency checkouts inside it) instead of the repository's parent directory, which on a machine with other paclets there loaded them all and let a Teukolsky checkout of a higher version shadow the repository, so that its code and its three test files were tested instead of the repository's; it now also refuses to run when the Teukolsky paclet that would be loaded is not the repository.

 - At machine precision TeukolskyRadial with "Amplitudes" -> False issued a spurious TeukolskyRadial::acc warning ("accuracy Infinity"), since its accuracy estimate needs B^inc and C^trans; with "RenormalizedAngularMomentum" -> False the radial functions were Indeterminate, although the integration from series boundary data does not need nu (the integrator has the NumericFunction attribute and evaluated to Indeterminate on the Indeterminate nu); and "BoundaryConditions" -> "In" or "Up" cost the same as both, because both were built for the estimate. The estimate is now made only when both solutions and the amplitudes are there, a non-numeric nu is passed on as a symbol (the MST fall-backs compute their own), and only the requested solution is built. Construction plus evaluation at six radii then takes a median 24 ms without the amplitudes and nu (both solutions), 9 ms for the "In" solution alone and 17 ms for the "Up" solution alone, against 112 ms by default (36 modes, machine precision).

 - TeukolskyRadial at a frequency far outside the range of the methods (omega = 270, say) ran for many minutes before giving up: the monodromy recurrences overflowed to an Indeterminate that the precision checks did not recognise, RenormalizedAngularMomentum returned an unevaluated expression instead of $Failed, and TeukolskyRadial went on to the amplitudes and their padded retries with that nu. It now fails within a second: the monodromy method returns $Failed (RenormalizedAngularMomentum::conv) on a non-numeric result, TeukolskyRadial stops with the new TeukolskyRadial::nufail message when nu could not be computed, and the sums of the amplitude formulae are bounded. At frequencies where nu is fine but the amplitude formulae overflow or underflow (omega = 50 for l = 2), the only symptoms used to be division-by-zero messages, an accuracy estimate of Infinity and a precision of Indeterminate for the MST "Up" solution; the new TeukolskyRadial::ampfail message now names the amplitudes that could not be computed, all amplitudes of that solution are Indeterminate, and the numerically integrated functions, which keep the unit-transmission normalisation of their boundary data, are returned.

 - Review fixes: the padding that the Wronskian check of TeukolskyRadial sets for a mode is now keyed by the parameters the MST series are evaluated with, so that at Re omega < 0, where the evaluations map to the conjugate partner (-m, -Conjugate[omega]), a failed check raises the precision of the radial series as well as of the amplitudes (before, only the amplitudes were recomputed and the retry could not repair the failure). RenormalizedAngularMomentum returns the representative l - ArcCos[Cos[2 Pi nu]]/(2 Pi) on every branch: for a real frequency with the monodromy cosine outside [-1, 1] it used to return 1/2 + i y or i y, equivalent values that differ from that representative by the integer l - 1 or l, so nu jumped where it turned complex. The series boundary data of the numerical integration are summed at the precision of their inputs, so "BoundaryMethod" -> "Series" at an arbitrary WorkingPrecision gives boundary data at that precision instead of machine numbers (which silently capped the solution at about 15 digits). Invalid values of "WronskianCheck" are rejected with TeukolskyRadial::optx instead of silently disabling the check.

### Changed
 - The MST "Up" solution below r+ + 1 is evaluated in the horizon basis of hypergeometric series, R_up = c1 R_in + c2 R_out, with R_out the series built on the second Kummer solution of the hypergeometric equation (normalised by Gamma[a-c+1] Gamma[b-c+1]/(Gamma[a] Gamma[b]) so that it obeys the same contiguous relations, and hence the same recurrences and coefficients, as the "In" series) and c1 = C^ref/B^trans, c2 = C^inc Exp[I (epsilon + tau) kappa (1/2 + Log[kappa]/(1 + kappa))]/Sum[a_n g_n] from the horizon amplitudes. The Coulomb-type series, an expansion about infinity, converges slowly near the horizon: an evaluation at r+ + 1/100 took 71 s at machine precision and 7 minutes at 40 digits, and now takes 0.1 s at any precision, at full accuracy (below 1e-40 against the Coulomb-type series at 32 digits, 1e-14 at machine precision, Wronskian identity satisfied down to r+ + 1/100). The Coulomb-type series is kept beyond r+ + 1, where the two cost the same, and at the degeneracies 2 I epsilon_+ = n, where the second Kummer solution has a logarithm.
 - MST evaluations at arbitrary precision are two to six times faster. Every evaluation used to cost two passes of the series, the first at the working precision to measure the loss to cancellation and the second with the padding that loss requires; the mode-dependent part of the loss is now cached per mode and series (for the hypergeometric "In" series after removing its radius-dependent part, about 0.8 |epsilon| (r - r+)/2 digits), and the next evaluation starts at the predicted padding, so that one pass suffices. The Coulomb-type "In" representation is now used from omega (r - r+) = 20 instead of 2: one of its passes costs two to four times one pass of the hypergeometric series at the same precision, independently of the radius, so it only pays once the series' loss exceeds about 20 digits (measured for s = -2, 0, 2 and omega = 0.1 to 2 at 32 digits, where the break-even lies between omega (r - r+) = 10 and 35). The accuracy is unchanged.

 - At machine precision the default method for a complex frequency is now MST rather than numerical integration. The two solutions of the radial equation differ by Exp[2 I omega r*], of modulus Exp[2 |Im omega| r*], so with Im omega < 0 the "In" solution is exponentially subdominant outwards and the "Up" solution inwards; the integration in either direction amplifies any error in the boundary data by Exp[2 |Im omega| (range in r*)], 1e-15 to 1e-6 over 40 in r at Im omega = -1/5, which over a limited range looks like a power of |omega| (this was the "violation of the Wronskian identity growing like |omega|^5" noted earlier: the identity holds, the integrated functions were wrong). At 32 digits the loss is absorbed and both methods agree to 1e-25; the MST solutions keep 1e-13 at machine precision. Method -> "NumericalIntegration" is still available there.
 - At machine precision the numerical integration (the default method there) now starts from series boundary data instead of precision-padded MST evaluations: a power series in r - r+ for the "In" solution at r+ + 1/5 and the large-r asymptotic series for the "Up" solution, evaluated at the radius (about 18/omega, chosen per mode) where its optimal truncation reaches 1e-15 and used directly beyond it. Integrating the "Up" solution inwards from there is unstable for negative spin (four digits lost per decade of radius for s = -2), so it is integrated at the flipped spin and mapped back with the Teukolsky-Starobinsky identity down to r+ + 1, and at its own spin from there to the horizon. The "In" solution of negative spin has the mirror-image problem beyond r+ + 1 whenever the reflection is small (at omega = 2 it lost four digits per decade of radius, 1e-8 at r = 50): beyond the radius where the large-r series converge (about 20 at omega = 2) it is now Binc R_ingoing + Bref R_up from the asymptotic series and the MST amplitudes (1e-14 at r = 50 and 300 instead of 1e-8 and worse); inside that radius the outward integration stays, since the inward integration is unstable wherever omega r << 1. Construction takes 0.07-0.3 s instead of 0.1-2 s (the remaining cost is the MST amplitudes), evaluation 2-9 ms, at the same 1e-12 to 1e-15 accuracy, now also next to the horizon and beyond the integration range (the derivative of the mapped "Up" solution, and the hand-over to the integration at the original spin, use the radial equation for R'' rather than the second derivative of the interpolating function, which cost two digits). The new "BoundaryMethod" sub-option of the NumericalIntegration method ("Series" or "MST", Automatic selects by working precision) restores the MST boundary data.
 - At the degeneracies 2 I epsilon_+ = n >= 1 - s of the MST method (for real frequencies the superradiant bound frequency omega = m Omega_H with s >= 1, on the imaginary axis also other n) the "In" solution is now returned normalised to unit incidence, with a vanishing transmission amplitude (TeukolskyRadial::innorm message), instead of failing: it is the horizon solution of larger exponent, the limit of the unit-incidence "In" solution from neighbouring frequencies, evaluated directly at the degenerate frequency. To make this possible the MST "In" series uses regularised hypergeometric functions, and the "UnscaledAmplitudes" of an "In" solution are those of Sasaki & Tagoshi divided by Gamma(1 - s - 2 I epsilon_+) at every frequency; the normalised amplitudes and radial functions are unchanged.
 - R[r, n] on a TeukolskyRadialFunction gives the n-th derivative (the same as Derivative[n][R][r]), and R[r, {0, 1}] the value and the first derivative, which for an MST solution come from a single summation of the series at about half the cost of the two separate evaluations. The boundary data of the numerical integration, the hand-over to the orbit-range integration and the Wronskian check use it.
 - Options of TeukolskyRadial given inside Method (e.g. Method -> {"MST", "RenormalizedAngularMomentum" -> nu}) are reported by the new TeukolskyRadial::topopt message as being in the wrong place.
 - Improved accuracy of the radial functions and of the source integrals:
   - The default PrecisionGoal is now WorkingPrecision - 2 for all methods (previously WorkingPrecision / 2, which limited machine-precision solutions to ~8 digits and 32-digit solutions to ~16 digits).
   - With Method -> Automatic at machine precision, TeukolskyRadial now returns numerically integrated solutions with goals of WorkingPrecision - 2, started from series boundary data (see the entry above), accurate to ~1e-12 - 1e-15 and cheap to evaluate on many radii; the accuracy of the pair is checked through the Wronskian, which must equal 2 I omega B^inc C^trans, and the new TeukolskyRadial::acc message reports a poor estimate.
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
