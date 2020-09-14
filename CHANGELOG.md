# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]


## [0.3.0] - 2020-09-14

### Added
 - Support for solving for a point particle with electric charge on a circular orbit in Kerr.
 - Support for computing solutions using numerical integration on a hyperbolical slice.
 - Support for computin "In" solutions using Mathematica's HeunC function (available since version 12.1).
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
