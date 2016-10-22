(* ::Package:: *)

basedir = Directory[];


(* ::Section:: *)
(*Meap creation functions*)


(* ::Text:: *)
(*meapCreator is a panel of buttons that generate a meap file from a ctf file*)


meapCreator:=Module[{},
	ctfFiles = FileNames["*.ctf", "*", Infinity];
	Panel[Column[{Style["Generate new meap file", Large], 
				Row[Button[FileNameSplit[#][[-3]], 
							Print[#];
							ctf2meap[#];
		]&/@ctfFiles]}]]
]


(* ::Text:: *)
(*ctf2meap converts ctf files to .meap files, which are used in this mathematica ebsd analysis package to quickly process information*)
(*input: .ctf file*)
(*output: exports .meap file, returns array*)
(*If you want to select a subset, create a .bmp image file called lens.bmp, where pixels to be included in the subset are white. the file must be the same dimensions of the sample*)


ctf2meap[file_]:=Module[{currentFile, title, fileBaseName, eulerInfo, width, lensFile, aragPhase, 
						height, stepSize, numPhases, phases, eulerList, eulerArray, misL, misU, vecArray,
						vecList, vecArrayRGB, imsize, zImage, yImage, outputList, outputArray,
						xImage, bcImage, bandImage, errorImage, madImage, bsImage, eulerImage, phaseImage, directoryLevel},
Catch[
		currentFile = file;
		
		(*-------CUSTOMIZE THIS-------*) 
		directoryLevel = -3;  (*if your ctf file is in the main directory for this sample, set this to -2. If it's in a subfolder of your main directory, set to -3*)

		fileBaseName = FileNameSplit[currentFile][[1;;directoryLevel]]; (*fileBaseName is the name of your main directory for this sample*)
		If[!DirectoryQ[FileNameJoin[Append[fileBaseName, #]]],
				CreateDirectory[FileNameJoin[Append[fileBaseName, #]]];
		]&/@{"Maps", "hkl", "Pole Figures"}; (*make directories for maps, raw data, and pole figures if necessary*)

		(*-------CUSTOMIZE THIS-------*) 
		lensFile = FileNameJoin[Join[fileBaseName, {"Maps", "lens.bmp"}]]; 	
		(*this is the name of your main subset file. if you name your subset files something else or put them somewhere else, change this. points in this image file should be pure white iff they are part of your main subset*)


		eulerInfo = Import[currentFile, "tsv"][[1;;20]]; (*euler info is a table of stats that help with processing*)
		width = eulerInfo[[5,2]];
		height = eulerInfo[[6,2]];
		stepSize = eulerInfo[[7,2]]; (*this assumes that the vertical and horizontal step sizes are the same*)
		numPhases = eulerInfo[[13,2]];
		phases = eulerInfo[[13+#, 3]]&/@Range[numPhases]; (*phases is a list of phases*)
		eulerList = Import[currentFile, "tsv"][[15+numPhases;;]]; (*eulerList is the actual table of data*)
		outputList = ConstantArray[0, {Length[eulerList]+2, 30}]; (*outputList creates a framework to export new information*)

		(*-------CUSTOMIZE THIS-------*)
		aragPhase = Position[phases, "aragonite-Pmcn"][[1,1]]; (*this program is specialized for aragonite. change this to your main phase. if you are mainly looking at multiple phases, play with this*)
		outputList[[1, 1;;8]] = {"width", width, "height", height, "step", stepSize, "aragPhase", aragPhase}; (*line 1 is basic info about the scan. if you want more phase numbers stored, add them to this line*)


		outputList[[2]] = {"Phase","X","Y","Bands","Error","Euler1","Euler2","Euler3",
							"MAD","BC","BS", "<100>", "<010>", "<001>", "Left angle", "Left axis", "Up angle", "Up axis", 
							"band im", "error im", "euler im", "MAD im", "BC im", "BS im", "Grain # (1 deg)", "Grain # (5 deg)", "Lens", "Subset 2", "X index", "Y index"};
		outputList[[3;;, 1;;11]] = eulerList; (*this copies data from the ctf file to the meap file. the remaining columns are calculated here and saved in the meap file to save time later*)

		(*-------CUSTOMIZE THIS-------*)
		outputList[[3;;]] = If[#[[1]]==aragPhase, #, Flatten[{{0, #[[2]], #[[3]], 0, 3, 0,0,0,0}, #[[10;;]]}]]&/@outputList[[3;;]]; (*this clears data that is not in your main phase. delete this if you want info about multiple phases*)

		vecList = eulerList;
		vecList[[;;, 6;;8]] = euler2vecs/@(vecList[[;;, 6;;8]]); (*veclist is a list of x,y,z vectors*)
		vecArray = Partition[vecList[[;;, 6;;8]], width];
		outputList[[3;;, 12;;14]] = vecList[[;;, 6;;8]];
	
		outputList[[3;;, 29;;30]] = Flatten[Array[{#2,#1}&, {height, width}],1]; (*this finds integer indices for x and y locations*)
		{misL, misU} = misMap2[vecArray]; (*misL and misU are arrays representing the misorientation left and up angle and vectors for each pixel*)
		outputList[[3;;, 15]] = misL[[;;, 3]];
		outputList[[3;;, 16]] = misL[[;;, 4;;6]];
		outputList[[3;;, 17]] = misU[[;;, 3]];
		outputList[[3;;, 18]] = misU[[;;, 4;;6]];

		outputList[[3;;, 19]] = SetPrecision[eulerList[[;;,4]]/Max[eulerList[[;;,4]]], 6]; (*this ensures that the data is stored as floats, not fractions*)
		outputList[[3;;, 20]] = SetPrecision[eulerList[[;;,5]]/Max[eulerList[[;;,5]]], 6];
		outputList[[3;;, 21]] = SetPrecision[eulerList[[;;,6;;8]]/180, 6];
		outputList[[3;;, 22]] = SetPrecision[eulerList[[;;,9]]/Max[eulerList[[;;,9]]], 6];
		outputList[[3;;, 23]] = SetPrecision[eulerList[[;;,10]]/Max[eulerList[[;;,10]]], 6];
		outputList[[3;;, 24]] = SetPrecision[eulerList[[;;,11]]/Max[eulerList[[;;,11]]], 6];

		If[FileExistsQ[lensFile], 
			outputList[[3;;, 27]] = If[#=={1,1,1}, 1, 0]&/@ Flatten[ImageData[ImageResize[Import[lensFile, "bmp"], width]],1]; (*if you defined a file on line 7, it imports it and stores it as a binary*)
			,
			outputList[[3;;, 27]] = 1;
			Print["No lens file found for ", fileBaseName]; (*if you did not define a file, everything is considered to be in your main subset*)
		];

		Export[FileNameJoin[Flatten[{fileBaseName, {"hkl"}, {StringJoin[fileBaseName[[-1]], ".meap"]}}]], outputList, "tsv"]; (*all work done so far is saved as a table, labeled as a .meap file*)
		ccl[{FileNameJoin[Flatten[{fileBaseName, {"hkl"}, {StringJoin[fileBaseName[[-1]], ".meap"]}}]]}]; (*this function uses connected component labeling to label grains and saves its work. if you do not want grains labeled, remove this line*)

]
]


(* ::Text:: *)
(*importLensSubset adds a lens subset to a meap file*)


importLensSubset[meapFile_]:=Module[{fileBaseName, lensFile, outputList, width, directoryLevel},
Catch[
		(*-------CUSTOMIZE THIS-------*) 
		directoryLevel = -3;  (*if your meap file is in the main directory for this sample, set this to -2. If it's in a subfolder of your main directory, set to -3*)

		fileBaseName = FileNameSplit[meapFile][[1;;directoryLevel]]; (*fileBaseName is the name of your main directory for this sample*)

		(*-------CUSTOMIZE THIS-------*) 
		lensFile = FileNameJoin[Join[fileBaseName, {"Maps", "lens.bmp"}]]; (*find corresponding subset file. if your main subset file is not named lens.bmp or in the folder Maps, change this*)

		outputList = Import[meapFile, "tsv"];
		width = outputList[[1, 2]];
		If[FileExistsQ[lensFile], 
			outputList[[3;;, 27]] = If[#=={1,1,1}, 1, 0]&/@ Flatten[ImageData[ImageResize[Import[lensFile, "bmp"], width]],1]; (*if the file exists, it is imported and incorporated into the table*)
			,
			outputList[[3;;, 27]] = 1;
			Print["No lens file found for ", fileBaseName];
		];
		Export[meapFile, outputList, "tsv"]; (*exports the file*)
]]


(* ::Text:: *)
(*ccl labels grains using connected component labeling*)
(*Input: .meap file*)
(*Output: .meap file with grains labeled using 1 degree and 5 degree critical misorientations*)


ccl[data_]:=Module[{nextGrain, array, width, array2, eqList, arr, j, k, left, up,  set, mcrit, ind, height, meapFile},
Catch[
If[Length[data]<2, (*the input value, data, can come in two formats: a meap file (string), or {meapfile, array containing all data in meapfile}*)
	array = Import[data, "tsv"];
	,
	meapFile = data[[1]];
	array = data[[2]];
];
width = array[[1,2]];
height = array[[1,4]];
array2 = array[[3;;]]; (*2d array*)
arr = Flatten[Array[{#1,#2}&, {height, width}],1]; (*arr is an array of indices*)

Do[
ind = 25+l-1; (*l is for 1 and 5, the critical misorientations used here. ind is the column number that these misorientations will go into*)

(*----------CUSTOMIZE THIS---------*)
mcrit = 1 + 4*(l-1); (*1 on first iteration, 5 on second iteration. if you want different critical misorientations to be permanently saved, change these values*)


nextGrain = 1; (*nextGrain indicates the number of the next new grain to be labeled*)
eqList = Transpose[{Range[Length[array]-2], Range[Length[array]-2]}]; (*eqlist is where you will store lists of equivalent grains which are mislabeled because of the sequential nature of connected component labeling*)
array2 = Partition[array2, width]; (*array2 is now 2d*)
(
j = #[[1]];
k = #[[2]];
If[array2[[j,k,1]]>0, (*only label grains that have value*)
	If[array2[[j,k,15]]>mcrit, (*check the edge for the critical misorientation*)
		If[array2[[j,k,17]]>mcrit,
			(*new grain*)
			array2[[j,k,ind]] = nextGrain;
			nextGrain = nextGrain+1;
			,
			(*same as up*)
			array2[[j,k,ind]] = array2[[array2[[j,k,30]]-1, array2[[j,k,29]], ind]];
		];
	,
		If[array2[[j,k,17]]>mcrit,
			(*same as left*)
			array2[[j,k,ind]] = array2[[array2[[j,k,30]], array2[[j,k,29]]-1, ind]];
			,
			(*same as left and up*)
			left = array2[[array2[[j,k,30]], array2[[j,k,29]]-1, ind]];
			up = array2[[array2[[j,k,30]]-1, array2[[j,k,29]], ind]];
			If[left!=up,
				array2[[j,k,ind]] = Min[eqList[[left,2]], eqList[[up,2]]];
				eqList[[Max[left, up], 2]] = Min[eqList[[left, 2]], eqList[[up, 2]]];
				,
				array2[[j,k,ind]] = eqList[[left, 2]];
			];
		];
	];
])&/@arr;

array2 = Flatten[array2, 1];

eqList = Select[eqList, #[[1]]>#[[2]]&];
eqList = Sort/@DeleteDuplicates/@Flatten/@Gather[eqList, #1[[2]]==#2[[2]]&]; (*sorts the equivalency list into arrays of equivalent numbers*)
	Do[
	(array2[[;;,25+l-1]]= array2[[;;, 25+l-1]]/.#[[i]]->#[[1]];)
	, {i, 2, Length[#]}]&/@eqList; (*relabels grains according to equivalency list*)
, {l, 2}];
array[[3;;]] = array2;
Export[meapFile, array, "tsv"];
(*Return[{array, eqList}];*)
]
]


(* ::Text:: *)
(*euler2vecs turns a single set of euler angles and turns it into a, b, and c axes*)
(*Input : euler angles e1, e2, and e3 in degrees*)
(*Output : Rotation matrix representing {a, b, c}*)


euler2vecs[eList__]:=Module[{phi1, phi, phi2, c1, c, c2, s1, s, s2, e1, e2, e3},
Catch[
	{e1, e2, e3} = eList;
	If[{e1, e2, e3} == {0,0,0}, Return[{{0,0,0},{0,0,0},{0,0,0}}];];
	phi1 = e1*Degree;
	phi = e2*Degree;
	phi2 = e3*Degree;
	c1 = Cos[phi1];
	c = Cos[phi];
	c2 = Cos[phi2];
	s1 = Sin[phi1];
	s = Sin[phi];
	s2 = Sin[phi2];
	Return[{{c1*c2-s1*s2*c, s1*c2 + c1*s2*c, s2*s}, {-c1*s2-s1*c2*c, -s1*s2+c1*c2*c, c2*s}, {s1*s, -c1*s, c}}];
];
]


(* ::Text:: *)
(*fixPrecision fixes a previous error in the code that stored values as fractions instead of floats*)


fixPrecision[meapFile_]:=Module[{x1},
Catch[
	x1 = Import[meapFile, "tsv"];
	x1[[3;;, {19, 20, 23, 24}]] = SetPrecision[N[ToExpression[x1[[3;;,{19, 20, 23, 24}]]]],6];
	Export[N[meapFile], x1, "tsv"];
]];


(* ::Section:: *)
(*Point selection functions*)


(* ::Text:: *)
(*pointClean cleans single-pixel grains out of data sets from meap files*)
(*Input: .meap file, list of indices to return, critical grain size (px)*)
(*Output: list of points*)


pointClean[meapFile_, iList_, critGSize_]:=Module[{data, grainHash, outputList},
Catch[
If[!NumberQ[critGSize], Return["List collection failed in pointClean"];];
	data = Import[meapFile, "tsv"][[3;;]];
	grainHash = SparseArray[(#[[1]]->#[[2]])&/@Tally[Select[data[[;;, 25]], #>0&]]];
	outputList = ToExpression[Select[data[[;;, Append[iList, 25]]], grainHash[[#[[Length[iList]+1]]]]>critGSize&][[;;, Range[Length[iList]]]]];
	Return[outputList];
]];


(* ::Text:: *)
(*selectLargestGrain selects the largest grain from a sample*)
(*Input: .meap file, list of indices to return*)
(*Output: list of points *)


selectLargestGrain[meapFile_, iList_]:=Module[{data, grainNum, outputList},
Catch[
	data = Import[meapFile, "tsv"][[3;;]];
	grainNum = MaximalBy[Tally[Select[data[[;;, 25]], #>0&]], Last][[1,1]];
	outputList = ToExpression[Select[data[[;;, Prepend[iList, 25]]], #[[1]]==grainNum&][[;;, Range[2,Length[iList]+1]]]];
	Return[outputList];
]];


(* ::Text:: *)
(*pointsFrommeap gets points from a meapfile, given certain parameters*)
(*Input: .meap file, bool only aragonite, number lens (0=all, 1=lens, 2=shell), bool clean, bool critical grain size (px), bool grain*)
(*Output: list of a,b,c vectors*)


pointsFrommeap[meapFile_, onlyArag_, onlyLens_, clean_, critSize_, grain_]:=Module[{pointList, aragNum},
Catch[
	If[grain,
		pointList = selectLargestGrain[meapFile, {1,12,13,14,27}];
		,
		If[clean,
			pointList = pointClean[meapFile, {1,12,13,14,27}, critSize];
			,
			pointList = ToExpression[Import[meapFile, "tsv"][[3;;, {1,12,13,14,27}]]];
		];
	];
	If[onlyArag,
		aragNum = Import[meapFile, "tsv"][[1,8]];
		pointList = Select[pointList, #[[1]]==aragNum&];
	];
	Switch[onlyLens
		,1,
			pointList = Select[pointList, #[[5]]==1&];
		,2,
			pointList = Select[pointList, #[[5]]==0&];
	];
	Return[pointList[[;;, 2;;4]]];
]]


(* ::Text:: *)
(*mapFrommeap collects points given certain parameters and returns an array*)
(*Input: .meap file, bool only aragonite, bool clean, number lens (0 = all, 1 = lens, 2 = shell), bool largest grain, bool critical grain size*)
(*Output: array of points containing all possible maps*)


mapFrommeap[meapFile_, onlyArag_, clean_,  lens_Integer, grain_, critSize_]:=Module[{pointList, aragNum, grainHash, width, grainNum} ,
Catch[
	pointList = fullListFromMeap[meapFile, onlyArag, clean, lens, grain, critSize, {1,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27}];
	width = getWidth[meapFile];
	Return[Partition[pointList, width]];
]]


(* ::Text:: *)
(*fullListFromMeap collects points given certain parameters and returns a list*)
(*Input: .meap file, bool only aragonite, bool clean, number lens (0 = all, 1 = lens, 2 = shell), bool largest grain, bool critical grain size*)
(*Output: list of points containing all data*)


fullListFromMeap[meapFile_, onlyArag_, clean_,  lens_Integer, grain_, critSize_, outList_]:=Module[{pointList, aragNum, grainHash, width, grainNum, fullList, default, x, y, z, w,
																	indexLength, indexList, grainIndex, lensIndex, index2, index3, index29, index30, index1} ,
Catch[
	fullList = Import[meapFile, "tsv"];
	indexLength = Length[outList];
	indexList = Join[outList, {25, 27, 2, 3, 29, 30, 1}];
	default = {0,x, y, 0,0,0,0,0,0,0,0,{0,0,0},{0,0,0},{0,0,0},361,{0,0,0},361,{0,0,0},0,0,{0,0,0},0,0,0,0,0,0,0,z, w}[[indexList]];
	grainIndex = indexLength + 1;
	lensIndex = indexLength +2;
	index2 = indexLength +3;
	index3 = indexLength+4;
	index29 = indexLength+5;
	index30 = indexLength+6;
	index1 = indexLength+7;
	pointList = fullList[[3;;, indexList]];
	If[grain,
		grainNum = Commonest[Select[pointList[[;;, grainIndex]], #>0&]][[1]];
		pointList= 
			ParallelMap[If[#[[grainIndex]]==grainNum, 
				#, 
				default/.{x->#[[index2]], y->#[[index3]], z->#[[index29]], w->#[[index30]]}
			]&,pointList];
		,
		If[clean,
			grainHash = Counts[pointList[[;;, grainIndex]]];
			pointList=
				ParallelMap[If[#[[grainIndex]]>0 && grainHash[#[[grainIndex]]]>critSize
					,
					#
					,
				default/.{x->#[[index2]], y->#[[index3]], z->#[[index29]], w->#[[index30]]}
				]&,pointList];
		];
	];
	If[onlyArag,
		aragNum = ToExpression[fullList[[1,8]]];
		pointList =
			ParallelMap[If[#[[index1]]==aragNum
				,
				#
				,
				default/.{x->#[[index2]], y->#[[index3]], z->#[[index29]], w->#[[index30]]}
			]&,pointList];
	];

	Switch[lens,
		1,
		pointList = 
			ParallelMap[If[#[[lensIndex]]==1
				,#,
				default/.{x->#[[index2]], y->#[[index3]], z->#[[index29]], w->#[[index30]]}
			]&,pointList];
		,2,
		pointList = 
			ParallelMap[If[#[[lensIndex]]==0
				,#,
				default/.{x->#[[index2]], y->#[[index3]], z->#[[index29]], w->#[[index30]]}
			]&,pointList];
	];
	Return[pointList[[;;, Range[indexLength]]]];
]]


(* ::Section:: *)
(*Interfaces*)


(* ::Subsection:: *)
(*EBSD*)


(* ::Text:: *)
(*EBSDBulkCreator is an interface that modifies global variables for use with kentBulk*)
(*kentBulk calculates kent distributions for several samples at once*)
(*Input: list of .meap files, bool only aragonite, bool only lens, bool clean, direction {xval, yval, zval}, critical grain size (px), bool grain*)
(*Output: popup window with table of results and export button*)


EBSDBulkCreator:=Module[{},
xVal = 0; yVal = 0; zVal = 1;
critSize = 1;
onlyLens = 1;
grain = False;
onlyArag = True; 
clean = True;
Panel[
Column[{
	Style["Bulk Kent Distribution Module for EBSD", "Section"],
	Row[Button[#, xVal = #[[1]]; yVal = #[[2]]; zVal = #[[3]];]&/@{{1,0,0}, {0,1,0}, {0,0,1}, {1,1,0}, {1,1,1}}],
				Row[{InputField[Dynamic[xVal], FieldSize->10], InputField[Dynamic[yVal], FieldSize->10], InputField[Dynamic[zVal], FieldSize->10]}],
	Row[{"Critical grain size (px) ", InputField[Dynamic[critSize], FieldSize->10]}], 
	SetterBar[Dynamic[onlyLens], {0->" All points ",1->" Only lens ",2->" Only shell "}],
	Row[{Checkbox[Dynamic[grain]], " Only largest grain\t", Checkbox[Dynamic[onlyArag]], " Only aragonite\t", Checkbox[Dynamic[clean]], " Clean"}]
}]]

]


kentBulk[meapFiles_, onlyArag_, onlyLens_, clean_, xVal_, yVal_, zVal_, critSize_, grain_]:=DynamicModule[{pointList, aragNum, kentResults},
Catch[
kentResults = ConstantArray[0, {Length[meapFiles]+1, 15}];
kentResults[[1]] = {"Sample", "N", "mean direction", "major axis", "minor axis", 
					"mean resultant length","q", "kappa", "beta", "kappa/beta", "Kstat",
					"StDev (theta) (Deg)", "StDev (phi) (Deg)", "mode", "dist"};
CreateDocument[
	Dynamic[Grid[kentResults]], WindowSize->{Scaled[1], Scaled[1]}, 
WindowTitle->StringJoin["Kent Distributions", If[onlyArag, "_arag_", "_"], Switch[onlyLens, 0, "_all_", 1, "_lens_", 2, "_shell_"],
						If[clean, StringJoin["_clean_", ToString[critSize]], "_"],  
						If[grain, "_grain_", "_all_"], ToString[xVal], ToString[yVal], ToString[zVal]]];
Do[
	kentResults[[i+1, 1]] = FileNameSplit[meapFiles[[i]]][[-3]];
	pointList = pointsFrommeap[meapFiles[[i]], onlyArag, onlyLens, clean, critSize, grain];
	If[Length[pointList]>0,
		pointList = {xVal, yVal, zVal}.#&/@(pointList);
		kentResults[[i+1, 2;;]] = kentSp[pointList][[;;,2]];
	];
,{i, Length[meapFiles]}]
Export[FileNameJoin[{basedir, 
				StringJoin["Kent Distributions", If[onlyArag, "_arag_", "_"], Switch[onlyLens, 0, "_all_", 1, "_lens_", 2, "_shell_"],
						If[clean, StringJoin["_clean_", ToString[critSize]], "_"],  
						If[grain, "_grain_", "_all_"], ToString[xVal], ToString[yVal], ToString[zVal], ".xlsx"]}], kentResults];
]];


(* ::Subsection:: *)
(*Pole figures*)


pfCreator[meapFiles_]:=DynamicModule[{pl,meapFile, xVal, yVal, zVal, critSize, kr,onlyLens, grain, onlyArag, clean, pointList, contour, ipf, pf1, pf2, pf3, kpf, buttColor},
Catch[
xVal = 0; yVal = 0; zVal = 1;
critSize = 1;
onlyLens = 1;
grain = False;
onlyArag = True; 
clean = True;
contour = False;
ipf = False;
buttColor=True;
ControlActive[onlyLens, buttColor=True;]

Panel[Column[{
	Style["Pole Figure Creator", "Section"]
	,
	Panel[Column[{
		SetterBar[Dynamic[meapFile, (buttColor=True; meapFile=#)&], (#->FileNameSplit[#][[-3]])&/@meapFiles, Appearance->"Row"]
	}]]
	,
	Panel[
		Column[{
			Row[{"Critical grain size (px) ", InputField[Dynamic[critSize, (buttColor=True; critSize=#)&], FieldSize->10]}]
			,
			SetterBar[Dynamic[onlyLens, (buttColor=True; onlyLens=#)&], {0->" All points ",1->" Only lens ",2->" Only shell "}]
			,
			Row[{Checkbox[Dynamic[grain, (buttColor=True; grain=#)&]], " Only largest grain\t"}]
			,
			Row[{Checkbox[Dynamic[onlyArag, (buttColor=True; onlyArag=#)&]], " Only aragonite\t"}] 
			,
			Row[{Checkbox[Dynamic[clean, (buttColor=True; clean=#)&]], " Clean"}]
			,
			Button["Select subset", 
					pointList = pointsFrommeap[meapFile, onlyArag, onlyLens, clean, critSize, grain]; buttColor=False;
					,Enabled->Dynamic[buttColor]
					,Background->Dynamic[If[buttColor, Pink, Automatic]]
				]
			,
			Row[Button[#, xVal = #[[1]]; yVal = #[[2]]; zVal = #[[3]];]&/@{{1,0,0}, {0,1,0}, {0,0,1}, {1,1,0}, {1,1,1}}]
			,
			Row[{InputField[Dynamic[xVal], FieldSize->10], InputField[Dynamic[yVal], FieldSize->10], InputField[Dynamic[zVal], FieldSize->10]}]
			,
			Row[{Checkbox[Dynamic[ipf]], " Inverse Pole Figure"}]
			,
			Row[{Checkbox[Dynamic[contour]], " Contour"}]
			,
			Button["Print point list", 
					Print[{xVal, yVal, zVal}.#&/@pointList];
				]
		}]]
	,
	Panel[Row[{
				Column[{Button["Generate pole figure", 
							If[ListQ[pointList], pf1 = pf[{xVal, yVal, zVal}.#&/@pointList, False, "", 
									StringJoin[FileNameSplit[meapFile][[-3]], " <", ToString[xVal], ToString[yVal], ToString[zVal], ">"], 
									False, contour, ipf]];],
				Dynamic[pf1]}]
				,
				Column[{Button["Generate pole figure", 
							If[ListQ[pointList],pf2 = pf[{xVal, yVal, zVal}.#&/@pointList, False, "", 
									StringJoin[FileNameSplit[meapFile][[-3]], " <", ToString[xVal], ToString[yVal], ToString[zVal], ">"], 
									False, contour, ipf]];],
				Dynamic[pf2]}]
				,
				Column[{Button["Generate pole figure", 
							If[ListQ[pointList],pf3 = pf[{xVal, yVal, zVal}.#&/@pointList, False, "", 
									StringJoin[FileNameSplit[meapFile][[-3]], " <", ToString[xVal], ToString[yVal], ToString[zVal], ">"], 
									False, contour, ipf]];],
				Dynamic[pf3]}]
	}]]
	,
	Panel[Column[{
			Button["Calculate Kent Distribution", 
						If[ListQ[pointList],pl = {xVal, yVal, zVal}.#&/@pointList;
						pl =  If[#[[3]]<0, -#, #]&/@pl;
						kr = kentSp[pl];
						kpf = Show[pf[pl, False, "", 
									StringJoin[FileNameSplit[meapFile][[-3]], " <", ToString[xVal], ToString[yVal], ToString[zVal], ">"], 
									False, False, False], ellSpPara[kr[[2,2]], kr[[12,2]], kr[[11,2]], kr[[4,2]], kr[[3,2]]]];
						]]
		,
			Row[{
			Dynamic[Grid[SetPrecision[kr,4]]],
			Dynamic[kpf], 
			Button["Export", 
						Export[FileNameJoin[Join[FileNameSplit[meapFile][[1;;-3]], {"pole figures", StringJoin[ "Kent stats ", ToString[xVal], ToString[yVal], ToString[zVal], ".tiff"]}]],  
								Row[{Grid[kr],kpf}]]]}]
	}]]
}]]]]


(* ::Text:: *)
(*pf creates pole figures*)
(*Input: list of vectors, bool export, string dir, string namestr, bool printf, bool contour, bool ipf*)
(*Output: if export, exports a pole figure. if print, prints a pole figure. returns a pole figure.*)


pf[directions_, export_, dir_, namestr_, printf_, contour_, ipf_]:=Module[{euler, Direction2D, directions2, ps, binned},

Catch[
	If[ipf, directions2 = Abs[directions], directions2 = directions];
	directions2 = Select[directions2, #[[3]]>0&];
	If[!contour, 
		directions2 = DeleteDuplicates[directions2];
		ps = 0.02/Log[Length[directions2]];
		directions2 = Round[directions2, ps/10];
		directions2 = DeleteDuplicates[directions2];
		Direction2D = ParallelMap[Graphics[{PointSize[ps], Point[{#[[1]], If[ipf, -1, 1]*#[[2]]}/(Norm[#]+#[[3]])]}]&, directions2];
	, 
		ps = 0.015;
		directions2 = Round[directions2, ps];
		binned = Reverse[Sort[{#[[1,1]]/Norm[#[[1]]+#[[1,3]]], #[[1,2]]/Norm[#[[1]]+#[[1,3]]], Length[#]}&/@Gather[directions2], #2[[3]]<#1[[3]]&]];
		binned[[;;,3]] = Log[binned[[;;,3]]]/Max[Log[binned[[;;,3]]]];
		Direction2D = ParallelMap[Graphics[{PointSize[ps], ColorData["TemperatureMap"][#[[3]]], Point[{#[[1]], If[ipf, -1, 1]*#[[2]]}]}]&, binned];
	];
	If[!ipf, 
		Direction2D = Show[Graphics[{EdgeForm[Thin], If[contour, ColorData["TemperatureMap"][0], White], Disk[]}], 
				Graphics[{Line[{{-1,0},{1,0}}], Line[{{0, -1}, {0,1}}],  Text[Style[namestr, Large], {0, 1.2}],
				Text[Style["x", Large], {1.2, 0}], Text[Style["y", Large], {0, -1.2}]}],  
				Direction2D, ImageSize->Medium];
		,
		Direction2D = Show[Graphics[{EdgeForm[Thin], If[contour, ColorData["TemperatureMap"][0], White],  Disk[{0,0}, 1, {3*Pi/2, 2*Pi}]}], 
				Graphics[{Line[{{0,0},{1,0}}], Line[{{0, -1}, {0,0}}],  Text[Style[namestr, Large], {0.5, 0.1}],
				Text[Style["x", Large], {1.1, -0.1}], Text[Style["y", Large], {0.1, -1.1}]}],  
				Direction2D, ImageSize->Medium];
	];

	If[printf, Print[Direction2D]];
	If[export,
	Export[FileNameJoin[{dir, StringJoin[{namestr,"_", "pf",".tiff"}]}], Show[Direction2D, Background->White, ImageSize->Full]];
	];
	Return[Direction2D];
]
]


(* ::Subsection:: *)
(*Plot deviation*)


plotDevCreator[meapFiles_]:=DynamicModule[{meapFile, xVal, yVal, zVal, critSize, kr,onlyLens, grain, pt,onlyArag, clean, pointList, contour, ipf, pf1, pf2, pf3, kpf, centroid, plotIndex, map1, kent, cf, misO},
Catch[
xVal = 0; yVal = 0; zVal = 1;
critSize = 1;
onlyLens = 1;
grain = False;
onlyArag = False; 
clean = True;
contour = True;
misO = False;
ipf = False;
plotIndex = 21;
cf = "TemperatureMap";
pt=5;

Panel[Column[{
	Style["Axis Deviation Map Creator", "Section"]
	,
	Panel[Column[{
		SetterBar[Dynamic[meapFile], (#->FileNameSplit[#][[-3]])&/@meapFiles, Appearance->"Row"]
	}]]
	,
	Panel[
		Column[{
			Row[{"Critical grain size (px) ", InputField[Dynamic[critSize], FieldSize->10]}]
			,
			Row[{"Max deviation (degrees)", InputField[Dynamic[pt], FieldSize->10]}]
			,
			Row[{"Average taken from: ", SetterBar[Dynamic[onlyLens], {0->" All points ",1->" Only lens ",2->" Only shell "}]}]	
			,
			Row[{Checkbox[Dynamic[grain]], " Only largest grain\t"}]
			,
			Row[{Checkbox[Dynamic[onlyArag]], " Only aragonite\t"}] 
			,
			Row[{Checkbox[Dynamic[clean]], " Clean"}]
			,
			Row[{Checkbox[Dynamic[contour]], " Contour (as opposed to 3D)"}]
			,
			Row[{Checkbox[Dynamic[misO]], " Misorientation overlay"}]
			,
			Row[Button[#, xVal = #[[1]]; yVal = #[[2]]; zVal = #[[3]];]&/@{{1,0,0}, {0,1,0}, {0,0,1}, {1,1,0}, {1,1,1}}]
			,
			Row[{InputField[Dynamic[xVal], FieldSize->10], InputField[Dynamic[yVal], FieldSize->10], InputField[Dynamic[zVal], FieldSize->10]}]
			,
			Dynamic[If[contour,
					Row[{"Plot color scheme ", SetterBar[Dynamic[cf], {"TemperatureMap", "AlpineColors", "GrayYellowTones", "GrayTones"}]}]
					, 
					 Row[{"Plot texture ", SetterBar[Dynamic[plotIndex], {21->" Euler ",19->" Band ",20->" Error ", 23->" BC ", 24->" BS "}]}]]]
			,
			Button["Generate map", 
							map1 = plotDevMeap[meapFile,{xVal, yVal, zVal}, plotIndex, onlyLens, grain, onlyArag, clean, critSize, contour, cf, misO, pt];
				]
		}]
	]
	,
	Panel[Dynamic[map1]]
	}]]
]]


getWidth[meapFile_String]:=ToExpression[StringSplit[Import[meapFile,{"Lines",1}], "\t"][[2]]];


plotDevMeap[meapFile_, vector_, ipfzFile_, lens_, grain_, arag_, clean_, critSize_,contour_, cf_, misO_, pt_, format___]:=
	Module[{stat, rat, ipfz,  plot, width, grainHash, c2, grainNum, centroid2, stat2, indexList, axes, str},

Catch[
	If[contour || NumberQ[ipfzFile], width = getWidth[meapFile];];
	If[NumberQ[ipfzFile],
		indexList = {2,3,12,13,14, ipfzFile};
		stat = fullListFromMeap[meapFile, arag, clean,  lens, grain, critSize, indexList];		
		If[!contour, 
			ipfz = Image[Partition[ToExpression[stat[[;;, 6]]], width]];
			If[misO, ipfz = Show[ipfz, misMap2Images[meapFile, 10, 31]];];];
		,
		indexList = {2,3,12,13,14};
		stat = fullListFromMeap[meapFile, arag, clean,  lens, grain, critSize, indexList];
		If[!contour, ipfz = Import[ipfzFile, format]];
	];

	stat2 = ToExpression[{#[[1]], #[[2]], vector.#[[3;;5]]}&/@(stat)];
	axes = Select[stat2[[;;, 3]], #!={0,0,0}&];
	If[Length[axes]==0, Return["No points selected"];];
	centroid2 = kentSp[axes][[2,2]];
	
	If[contour && lens>0,
		stat = fullListFromMeap[meapFile, arag, clean,  0, grain, critSize, indexList];
		stat2 = ToExpression[{#[[1]], #[[2]], vector.#[[3;;5]]}&/@(stat)];	
	];	
	stat2 = If[#[[3]]=={0,0,0}, {#[[1]], #[[2]], centroid2}, #]&/@stat2;
	stat2[[;;,3]] = VectorAngle[centroid2, #[[3]]]&/@stat2/Degree;

	If[!contour, stat2[[;;, 2]] = Max[stat2[[;;, 2]]] - stat2[[;;,2]]];


	If[contour,
	plot = ArrayPlot[Partition[stat2[[;;, 3]], width], PlotLegends->Automatic, PlotRange->{0,pt}, ColorFunction->cf, Frame->False, ImageSize->Large];
	,
	plot = ListPlot3D[stat2, 
			Mesh->None, 
			PlotStyle->{Directive[Texture[ipfz]]}, 
			AxesLabel->{"microns", "microns", "Divergence (\[Degree])"}, 
			PlotRange->{0, pt}, 
			ClippingStyle->None, 
			ImageSize->Large,
			ViewPoint->{0, -1, 2},
			ViewVertical->{0,1,0},
			PlotLabel->meapFile];
	];
	If[misO,	plot = Show[plot, misMap2Images[meapFile, 10, 31]];];
	If[NumberQ[ipfzFile],
		str = StringJoin["PlotDev", 
					If[arag, "_arag", ""], 
					Switch[lens, 0, "_all", 1, "_lens", 2, "_shell"],
					If[clean, StringJoin["_clean", "_",ToString[critSize]]],  
					If[grain, "_grain_", "_all_"], 
					If[contour, "contour_", ToString[ipfzFile]],"_",
					If[misO, StringJoin["mis_", cf, "_"], ""],
					ToString[vector], ".tiff"];
		Export[getFileName[meapFile, "Maps", str], plot]; 
	];

	Return[plot];
Return[stat2];
];
]


getFileName[meapFile_, subdir_, str_]:=Module[{dir},Catch[
	dir = FileNameSplit[meapFile][[1;;-3]];
	Return[FileNameJoin[Join[dir, {subdir}, {str}]]]
]]


(* ::Subsection:: *)
(*Maps*)


mapCreator[meapFiles_]:=DynamicModule[{meapFile, xVal, yVal, zVal, critSize, lens,
										 grain, onlyArag, clean, pointList, contour, ipf, 
										pf1, pf2, pf3, kpf, phaseIm, xIm, yIm, zIm, bandIm, errorIm, 
										eulerIm, MADIm, bcIm, bsIm, grain1Im, grain5Im, lensIm},
Catch[
xVal = 0; yVal = 0; zVal = 1;
critSize = 1;
onlyArag = True; 
clean = False;
lens=0;
grain=False;
phaseIm= xIm=yIm=zIm= bandIm= errorIm=eulerIm= MADIm= bcIm= bsIm= grain1Im= grain5Im= lensIm="";

SetOptions[Image, ImageSize->Medium];
SetOptions[Colorize, ImageSize->Medium];

Panel[Column[{
	Style["Map Creator", "Section"]
	,
	Panel[Column[{
		SetterBar[Dynamic[meapFile, (buttColor=True; meapFile=#)&], (#->FileNameSplit[#][[-3]])&/@meapFiles, Appearance->"Row"]
	}]]
	,
	Panel[
		Column[{
			Row[{"Critical grain size (px) ", InputField[Dynamic[critSize, (buttColor=True; critSize=#)&], FieldSize->10]}]
			,
			SetterBar[Dynamic[lens, (buttColor=True; lens=#)&], {0->" All points ",1->" Only lens ",2->" Only shell "}]
			,
			Row[{Checkbox[Dynamic[onlyArag, (buttColor=True; onlyArag=#)&]], " Only aragonite\t"}] 
			,
			Row[{Checkbox[Dynamic[clean, (buttColor=True; clean=#)&]], " Clean"}]

			,
			Row[{Checkbox[Dynamic[grain, (buttColor=True; grain=#)&]], " Only largest grain"}]
			,
			Button["Select subset", 
					pointList = mapFrommeap[meapFile, onlyArag, clean,  lens, grain, critSize]; buttColor=False;
					,Enabled->Dynamic[buttColor]
					,Background->Dynamic[If[buttColor, Pink, Automatic]]
				]
		}]]
	, 
	Panel[Grid[{
				{Column[{Button["x image", 
							xIm = Image[ToExpression[pointList[[;;,;;,2]]]];]
							,Dynamic[xIm]}],
				Column[{Button["y image", 
							yIm = Image[ToExpression[pointList[[;;,;;,3]]]];]
							,Dynamic[yIm]}],
				Column[{Button["z image", 
							zIm = Image[ToExpression[pointList[[;;,;;,4]]]];]
							,Dynamic[zIm]}]},
				{Column[{Button["Euler image", 
							eulerIm = Image[ToExpression[pointList[[;;,;;,11]]]];]
							,Dynamic[eulerIm]}],
				Column[{Button["Grain (1 deg) image", 
							grain1Im = Colorize[pointList[[;;,;;,15]]];]
							,Dynamic[grain1Im]}],
				Column[{Button["Grain (5 deg) image", 
							grain5Im = Colorize[pointList[[;;,;;,16]]];]
							,Dynamic[grain5Im]}]},
				{Column[{Button["Phase image", 
							phaseIm = Image[pointList[[;;,;;,1]]];]
							,Dynamic[phaseIm]}],
				Column[{Button["Lens image", 
							lensIm = Image[pointList[[;;,;;,17]]];]
							,Dynamic[lensIm]}],
				Column[{Button["MAD image", 
							MADIm = Image[pointList[[;;,;;,12]]];]
							,Dynamic[MADIm]}]},
				{Column[{Button["Band image", 
							bandIm = Image[pointList[[;;,;;,9]]];]
							,Dynamic[bandIm]}],
				Column[{Button["Band contrast image", 
							bcIm = Image[pointList[[;;,;;,13]]];]
							,Dynamic[bcIm]}],
				Column[{Button["Band slope image", 
							bsIm = Image[pointList[[;;,;;,14]]];]
							,Dynamic[bsIm]}]}
	}]]
}]]]]


(* ::Subsection:: *)
(*Misorientation overlays*)


overlayCreator[meapFiles_]:=DynamicModule[{critM, imIn, bcoIm},
critM = 10; imIn = 23;
Panel[Column[{Style["Generate misorientation angle overlay", Large],
Row[{"Critical misorientation angle (degrees) ", InputField[Dynamic[critM], Number], 
			"\t Image index ", InputField[Dynamic[imIn], Number]}],
Row[{Grid[Transpose[{Range[30], {"Phase","X","Y","Bands","Error","Euler1","Euler2","Euler3",
							"MAD","BC","BS", "<100>", "<010>", "<001>", 
						"Left angle", "Left axis", "Up angle", "Up axis", 
							"band im", "error im", "euler im", "MAD im", "BC im", "BS im", 
							"Grain # (1 deg)", "Grain # (5 deg)", "Lens", "Subset 2", "X index", 
							"Y index"}}]], 
			Column[{Grid[Partition[Button[FileNameSplit[#][[-3]], 
								bcoIm = misMap2Images[#, critM, imIn]]&/@meapFiles, 3]], 
						SwatchLegend[Hue[#/120]&/@{12, 28, 54, 64, 72, 90, 104}, {12, 28, 54, 64, 72, 90, 104}, LegendMarkerSize->30, LegendLayout->"Row"]}]
			,"\t" ,
								Dynamic[Show[bcoIm, ImageSize->Large]]}]}]]
]


(* ::Subsection:: *)
(*Data*)


dataPlay[meapFiles_]:=DynamicModule[{critSize=1, lens=0, grain=False, onlyArag=False, clean=False, meapFile},
Catch[
critSize = 1;
lens = 0;
grain=onlyArag=clean=False;
Panel[Row[{
	Grid[Transpose[{Range[30], 
			{"Phase","X","Y","Bands","Error","Euler1","Euler2","Euler3",
			"MAD","BC","BS", "<100>", "<010>", "<001>", "Left angle", "Left axis", "Up angle", "Up axis", 
			"band im", "error im", "euler im", "MAD im", "BC im", "BS im", "Grain # (1 deg)", 
			"Grain # (5 deg)", "Lens", "Subset 2", "X index", "Y index"}}]]
	,
	"\t"
	,
	Column[{
	SetterBar[Dynamic[meapFile], (#->FileNameSplit[#][[-3]])&/@meapFiles, Appearance->"Vertical"->{Automatic, 3}]
	,
	""
	,
	Panel[Column[{
			Row[{"Critical grain size (px) ", InputField[Dynamic[critSize], FieldSize->10]}]
			,
			SetterBar[Dynamic[lens], {0->" All points ",1->" Only lens ",2->" Only shell "}]
			,
			Row[{Checkbox[Dynamic[grain]], " Only largest grain\t"}]
			,
			Row[{Checkbox[Dynamic[onlyArag]], " Only aragonite\t"}] 
			,
			Row[{Checkbox[Dynamic[clean]], " Clean"}]
			,
			Button["set data", 
				data =fullListFromMeap[meapFile, onlyArag, clean,  lens, grain, critSize, Range[30]];]
		}]]
	}]
}]]
]]
