(* ::Package:: *)

Paclet[
  Name -> "Teukolsky",
  Version -> "0.1.0",
  MathematicaVersion -> "8+",
  Creator -> "Barry Wardell",
  Description -> "A set of functions for computing solutions to the Teukolsky equation.",
  Extensions ->
  {
    { "Kernel",
	    "Context" -> {
        "Teukolsky`",
        "Teukolsky`RenormalizedAngularMomentum`",
        "Teukolsky`TeukolskyRadialUp`",
        "Teukolsky`MST`"
      }
	  },

    { "Documentation",
      Language -> "English", 
      MainPage -> "Guides/Teukolsky",
      Resources -> 
     	{
        "Guides/Teukolsky"
      }
    }
  }
]
