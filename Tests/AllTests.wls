#!/usr/bin/env wolframscript

Print["Mathematica Version: ", $Version];

withinRoundoff[a_, b_] :=
  (* Ignore the last N machine digits *)
  Block[{Internal`$SameQTolerance = 4},
    SameQ[a,b]];

swshdir = FileNameJoin[{FileNameDrop[ExpandFileName[First[$ScriptCommandLine]],-2], "SpinWeightedSpheroidalHarmonics"}];
kerrgeodir = FileNameJoin[{FileNameDrop[ExpandFileName[First[$ScriptCommandLine]],-2], "KerrGeodesics"}];

If[$VersionNumber >= 12.1,
  PacletDirectoryLoad[FileNameDrop[ExpandFileName[First[$ScriptCommandLine]],-3]];
  If[DirectoryQ[swshdir], PacletDirectoryLoad[swshdir]];
  If[DirectoryQ[kerrgeodir], PacletDirectoryLoad[kerrgeodir]];
  testFiles = FileNames["*.wlt",
    FileNameJoin[{PacletObject["Teukolsky"]["Location"], "Tests"}]];
  ,
  PacletDirectoryAdd[FileNameDrop[ExpandFileName[First[$ScriptCommandLine]],-3]];
  If[DirectoryQ[swshdir], PacletDirectoryAdd[swshdir]];
  If[DirectoryQ[kerrgeodir], PacletDirectoryAdd[kerrgeodir]];
  testFiles = FileNames["*.wlt",
    FileNameJoin[{Lookup[PacletInformation["Teukolsky"], "Location"], "Tests"}]];
];

<< Teukolsky`

reports = TestReport /@ testFiles;

numTestsSucceeded = Total[#["TestsSucceededCount"] & /@ reports];
numTestsFailed    = Total[#["TestsFailedCount"]& /@ reports];
numTests          = numTestsSucceeded + numTestsFailed;

time = QuantityMagnitude[Total[(#["TimeElapsed"] & /@ reports)], "Seconds"];

(* Generate XML for a single test result *)
testResultXML[result_] :=
  XMLElement["testcase",
    {"name" -> ToString[result["TestID"]], 
     "time" -> ToString[QuantityMagnitude[result["AbsoluteTimeUsed"], "Seconds"]]}, 
    Switch[result["Outcome"],
      "Failure",
        {XMLElement["failure", {}, {XMLObject["CDATASection"][
          "Input: " <> ToString[result["Input"], InputForm] <> "\n" <>
          "Expected output: " <> ToString[result["ExpectedOutput"], InputForm] <> "\n" <>
          "Actual output: " <> ToString[result["ActualOutput"], InputForm]]}]},
      "MessagesFailure",
        {XMLElement["failure", {}, {XMLObject["CDATASection"][
          "Expected Messages: " <> ToString[result["ExpectedMessages"], InputForm] <> "\n" <>
          "Actual Messages: " <> ToString[result["ActualMessages"], InputForm]]}]},
      "Error",
          {XMLElement["failure", {}, {XMLObject["CDATASection"]["Error"]}]},
      "Success",
          {}
    ]
];

testsuiteXML[report_] :=
 Module[{xml},
  If[MatchQ[report["TestResults"], <| {"Event", "Fatal"} -> __ |>],
    xml = XMLElement["testsuite",
      {"name" -> report["Title"],
       "tests" -> "1",
       "failures" -> "1",
       "time" -> ToString[QuantityMagnitude[report["TimeElapsed"], "Seconds"]]},
      {XMLElement["testcase", {"name" -> "Loading test file"},
                  {XMLElement["failure", {}, {XMLObject["CDATASection"]["Error reading file"]}]}]}
      ];,
    xml = XMLElement["testsuite",
      {"name" -> report["Title"],
       "tests" -> ToString[report["TestsSucceededCount"] + report["TestsFailedCount"]],
       "failures" -> ToString[report["TestsFailedCount"]],
       "time" -> ToString[QuantityMagnitude[report["TimeElapsed"], "Seconds"]]},
        Map[testResultXML, Values[report["TestResults"]]]];
  ];
  xml
];

xml = XMLObject["Document"][{XMLObject["Declaration"]["Version" -> "1.0",  "Encoding" -> "UTF-8"]}, 
  XMLElement[
    "testsuites",
    {"id" -> DateString["ISODateTime"],
     "name" -> "Teukolsky Tests",
     "time" -> ToString[time],
     "tests" -> ToString[numTests], 
     "failures" -> ToString[numTestsFailed]}, 
     testsuiteXML /@ reports
     ], {}];

Export["TestReport.xml", xml, "XML"]
