<< Teukolsky`;
<< ApplicationTools`;

packages =
{ 
  "MST`RenormalizedAngularMomentum",
  "TeukolskyMode",
  "TeukolskyRadial",
  "TeukolskySource"
};

packageSymbols = Map[# -> DocumentedSymbols["Teukolsky", #] &, packages];

undocumentedSymbols = Map[# -> UndocumentedSymbols["Teukolsky", #] &, packages] /. (_ -> {}) -> Sequence[];
Map[Print["Undocumented symbols for package "<>#[[1]]<>" skipped:\n", #[[2]]]&, undocumentedSymbols];

Print["Building symbol reference pages"];
docPackage[package_ -> symbols_] :=
  Map[(Print[#]; BuildSymbolReference["Teukolsky", #, "Source"]) &, symbols];
Scan[docPackage, packageSymbols];

Print["Building guides"];
sourceGuides = FileNames["*.md", FileNameJoin[{"Source", "Documentation", "English", "Guides"}], Infinity];
destGuides =
  FileNameJoin[{Directory[], FileNameDrop[DirectoryName[#], 1],
      FileBaseName[#] <> ".nb"}] & /@ sourceGuides;
MapThread[BuildGuide, {sourceGuides, destGuides}];

Print["Building tutorials"];
tutorialSources = FileNames["*.md", FileNameJoin[{"Source", "Documentation", "English", "Tutorials"}], Infinity];
Map[(Print[#]; BuildTutorial[FileNameJoin[{Directory[], #}]])&, tutorialSources];

Print["Indexing Documentation"];
BuildIndex["Teukolsky"];

Print["Done"];
