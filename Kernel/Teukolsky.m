(* ::Package:: *)

BeginPackage["Teukolsky`"];

$TeukolskyInformation::usage = "$TeukolskyInformation is a list of rules that gives information about the version of the Teukolsky package you are running.";
$TeukolskyInstallationDirectory::usage = "$TeukolskyInstallationDirectory gives the top-level directory in which the Teukolsky package is installed.";

$TeukolskyVersionNumber::usage = "$TeukolskyVersionNumber is a real number which gives the current version number for the Teukolsky package.";
$TeukolskyReleaseNumber::usage = "$TeukolskyReleaseNumber is an integer which gives the current release number for the Teukolsky package.";
$TeukolskyVersion::usage = "$TeukolskyVersionNumber is a string that gives the version of the Teukolsky package you are running.";

Begin["`Private`"];

$TeukolskyInstallationDirectory = FileNameDrop[FindFile["Teukolsky`"], -2];

$TeukolskyVersionNumber        = 1.0;
$TeukolskyReleaseNumber        = 0;

$TeukolskyVersion :=
 Module[{path, version, release, buildid, gitrev, gitdir},
  path = $TeukolskyInstallationDirectory;
  version = ToString[NumberForm[$TeukolskyVersionNumber, {Infinity, 1}]];
  release = ToString[$TeukolskyReleaseNumber];

  buildid = Quiet@ReadList[FileNameJoin[{path, "BUILD_ID"}], "String"];
  If[SameQ[buildid, $Failed],
    buildid = "";
  ,
    buildid = " (" <> First[buildid] <> ")";
  ];

  (* First, check for a GIT_REVISION file. If it exists, use its contents as the revision. *)
  gitrev = Quiet@ReadList[FileNameJoin[{path, "GIT_REVISION"}],"String"];

  (* Otherwise, try to determine the git revision directly *)
  If[SameQ[gitrev, $Failed],
    gitdir = FileNameJoin[{path, ".git"}];
    If[FileType[gitdir] === Directory,
      gitrev = Quiet@ReadList["!git --git-dir "<>gitdir<>" rev-parse HEAD", String];
      If[gitrev === {}, gitrev = $Failed];
    ];
  ];

  (* If it worked, ReadList returns a list but we just want the first element (line) *)
  If[Head[gitrev] === List, gitrev = First[gitrev]];

  (* Check we have a git revision and otherwise give up trying *)
  If[Head[gitrev] === String && StringMatchQ[gitrev, RegularExpression["[0-9a-f]{5,40}"]], gitrev = " (" <> gitrev <> ")", gitrev = ""];

  version <> "." <> release <> buildid <> gitrev
]

$TeukolskyInformation :=
  {"InstallationDirectory" -> $TeukolskyInstallationDirectory,
   "Version" -> $TeukolskyVersion,
   "VersionNumber" -> $TeukolskyVersionNumber,
   "ReleaseNumber" -> $TeukolskyReleaseNumber}

End[];
EndPackage[];
