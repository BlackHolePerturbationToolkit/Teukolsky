VerificationTest[
  orbit = KerrGeoOrbit[0.3`32, 10.`32, 0, 1];
  TeukolskyPointParticleMode[-2, 2, 2, 0, 0, orbit];
  Names["Teukolsky`*`Private`*$" ~~ DigitCharacter ..],
  {},
  TestID -> "Memory leaks"
]
