(* ::Package:: *)

(* ::Text:: *)
(*These functions are tailored to orthorhombic systems!*)


(* ::Section:: *)
(*Misorientation*)


getMisAngle[mat1_, mat2_]:=Module[{MM, MMatrix, MAngle, M1, M2, MM2, tempAngle},

MM = mat1.Inverse[mat2];
			MMatrix = MM;
			MAngle = ArcCos[(MMatrix[[1,1]] + MMatrix[[2,2]] + MMatrix[[3,3]] -1)/2];
			Do[
				Do[
					If[l!=k,
						M1 = orthoSym[k];
						M2 = orthoSym[l];
						MM2 = M1.MM.Inverse[M2];
						tempAngle = ArcCos[(Tr[MM2] -1)/2];
						If[Im[tempAngle]==0, 
							If[tempAngle<MAngle, MAngle=tempAngle; MMatrix = MM2;];
						];
					];
				,{l,k}];
			,{k, 8}];

Return[{MAngle, MMatrix}];
]


(* ::Subsection:: *)
(*Misorientation Map Processor*)


(* ::Text:: *)
(*misMap outputs a file that can be used by misMapImages to make maps and pole figures*)
(*Input: OIM output file ctfFile, critical misorientation in degrees mcrit*)
(*Output: Exports a _band.txt and a _List.txt file*)


misMap[ctfFile_, mcrit_]:=Module[{str, dir, euler, eulers, phi1, phi, phi2, c1, c, c2, s1, s, s2, minx, miny,maxx, maxy, output,
									misorientationList,criticalMisorientation, upNeighbor, leftNeighbor,matrix, MMatrix, MAngle,
									Lines, Lines2, Directions, Direction2D, AxisAnglePlot, MMatrix2, MAngle2, AngleHistogram, bandImage, band, minmin, M1, M2, MM, tempAngle
									, MM2},
str = StringSplit[FileNameSplit[ctfFile][[-3]], " "][[1]];
dir = FileNameJoin[{ParentDirectory[DirectoryName[ctfFile]], "Maps"}];
If[!DirectoryQ[dir], CreateDirectory[dir]];


(*-------------------Point method---------------*)
(*Input data from file and turn into an array of position and direction vectors*)
euler = Import[ctfFile, "tsv"][[17;;,{2,3,4,6,7,8}]];  
euler = Select[euler, #[[3]]>0&];

(*find rotation of axes: xvecs is the rotated <100>, etc., ref. Euler angle wiki*)
(*eulers = ConstantArray[0, {Length[euler], 2}];
Do[
	eulers[[i,1]] = euler[[i, 1;;3]];
	eulers[[i,2]] = euler2vecs[euler[[i,4]], euler[[i,5]], euler[[i,6]]];
,{i, Length[euler]}];*)
eulers = {#[[1;;3]], euler2vecs[euler[[i,4;;6]]]}&/@eulers;

(*convert x and y indices to integers*)
minx = Min[Select[eulers[[;;,1,1]], #>0&]];
miny = Min[Select[eulers[[;;,1,2]], #>0&]];
minmin = Min[minx, miny];
eulers[[;;, 1, 1]] = Round[eulers[[;;,1,1]]/minmin]+2;
eulers[[;;,1,2]] = Round[eulers[[;;,1,2]]/minmin]+2;
Clear[minx, miny];
(*rearrange eulers into an array that reflects locations*)
maxx = Max[eulers[[;;,1,1]]];
maxy = Max[eulers[[;;,1,2]]];
output = ConstantArray[{{},{},{}}, {maxx+2, maxy+2}];
band = ConstantArray[Max[eulers[[;;,1,3]]]+1, {maxx+2, maxy+2}];
Do[
	output[[eulers[[i, 1, 1]], eulers[[i, 1, 2]]]] = eulers[[i,2]];
	band[[eulers[[i, 1, 1]], eulers[[i, 1, 2]]]] = eulers[[i,1,3]];
,{i, Length[eulers]}];
Export[FileNameJoin[{dir, StringJoin[{str, "_band.txt"}]}], band, "tsv"];

PrintTemporary["Starting misorientation calculations. There are " , Length[eulers], " points."];

(*misorientation list = {point1x point1y point2x point2y angle vector1 vector2 vector3}*)
misorientationList = ConstantArray["", {Length[eulers]*2, 8}];
misoCounter = 1;
PrintTemporary["Collected ", Dynamic[misoCounter], " orientations."];

Do[
	(*find neighbors*)
	upNeighbor = output[[eulers[[i,1,1]], eulers[[i,1,2]]-1]];
	leftNeighbor = output[[eulers[[i,1,1]]-1, eulers[[i,1,2]]]];
	(*find misorientation angle and axis*)
	Do[
		If[j==1, matrix = upNeighbor;, matrix=leftNeighbor;];
		If[matrix!={{},{},{}}, 
			(*find disorientation (minimum angle) using symmetry operators*)

			(*MAngle = getMisAngle[eulers[[i,2]], matrix];*)

			MM = eulers[[i,2]].Inverse[matrix];
			MMatrix = MM;
			MAngle = ArcCos[(MMatrix[[1,1]] + MMatrix[[2,2]] + MMatrix[[3,3]] -1)/2];
			Do[
				Do[
					If[l!=k,
						M1 = orthoSym[k];
						M2 = orthoSym[l];
						MM2 = M1.MM.Inverse[M2];
						tempAngle = ArcCos[(Tr[MM2] -1)/2];
						If[Im[tempAngle]==0, 
							If[tempAngle<MAngle, MAngle=tempAngle; MMatrix = MM2;];
						];
					];
				,{l,k}];
			,{k, 8}];			

			If[MAngle>(mcrit*Degree),
				If[j>1, (*left neighbor*)
					misorientationList[[misoCounter, 1]] = eulers[[i, 1,1]]-1;
					misorientationList[[misoCounter, 2]] = maxy-(eulers[[i, 1,2]])+2;
					,
					misorientationList[[misoCounter, 1]] = eulers[[i, 1,1]];
					misorientationList[[misoCounter, 2]] = maxy-(eulers[[i, 1,2]])+3;
				];
				misorientationList[[misoCounter, 3]] = eulers[[i, 1,1]]-1;
				misorientationList[[misoCounter, 4]] = maxy-(eulers[[i, 1,2]])+3;
				misorientationList[[misoCounter, 5]] = MAngle/Degree;
				misorientationList[[misoCounter, 6]] = (MMatrix[[2,3]] - MMatrix[[3,2]])/2*Csc[MAngle];
				misorientationList[[misoCounter, 7]] = (MMatrix[[3,1]] - MMatrix[[1,3]])/2*Csc[MAngle];			
				misorientationList[[misoCounter, 8]] = (MMatrix[[1,2]] - MMatrix[[2,1]])/2*Csc[MAngle];
				If[misorientationList[[misoCounter,8]]<0, misorientationList[[misoCounter,6;;8]] = -misorientationList[[misoCounter, 6;;8]];];
				misoCounter = misoCounter+1;
			];
		];
	,{j, 2}];
, {i, Length[eulers]}];
misorientationList = misorientationList[[1;;misoCounter-1]];
Export[FileNameJoin[{dir, StringJoin[{str, "_List.txt"}]}], misorientationList, "tsv"];

];


(* ::Text:: *)
(*misMap2 finds misorientations between pixels in an array*)
(*Input: array vecArray (produced by ctf2meap)*)
(*Output: {array of left misorientations, array of up misorientations}*)


misMap2[vecArray_]:=Module[{vecArray2, arr, misArray, 
											h, w, misArrayLeft, misArrayUp, AngleLinesLeft, AngleLinesUp, bco},
Catch[
vecArray2 = ConstantArray[0, Dimensions[vecArray]+{2,2,0,0}];
vecArray2[[2;;-2, 2;;-2]] = vecArray;
arr = Array[{#1,#2}&, Dimensions[vecArray2][[1;;2]]-{2,2}];
arr = arr + ConstantArray[1, Dimensions[arr]];
arr = Flatten[arr, 1];
misArray = ParallelMap[findMis[vecArray2, #[[1]], #[[2]]]&,arr];

misArrayLeft = Flatten[{arr[[#]], misArray[[#, 1]]}]&/@Range[Length[arr]];
misArrayUp = Flatten[{arr[[#]], misArray[[#, 2]]}]&/@Range[Length[arr]];

misArray = {misArrayLeft, misArrayUp};

Return[misArray];

];
];


(* ::Text:: *)
(*findMis is a lower-level function that finds misorientation angles and axes between neighboring pixels*)
(*input: an array of 3x3 matrices and integer indices pointing to a place in the array to probe.*)
(*Output: a 2x4 matrix that contains angles (in degrees) and axes (in direction cosines) of border between pixel and neighbors above and left *)


findMis[vecArray2_, i_, j_]:=Module[{upNeighbor, leftNeighbor, matrix, MAngle, MMatrix, retval},
Catch[
retval = {{361,0,0,0}, {361,0,0,0}}; (*this allows higher level functions to detect borders with non-indexed pixels - it marks the misorientation angle as 361 degrees*)
If[vecArray2[[i,j]]!={{0,0,0},{0,0,0},{0,0,0}}, 
		(*find neighbors*)
		upNeighbor = vecArray2[[i, j-1]];
		leftNeighbor = vecArray2[[i-1, j]];
		(*find misorientation angle and axis*)
		Do[
			If[k==1, matrix = upNeighbor;, matrix=leftNeighbor;];
			If[matrix!={{0,0,0},{0,0,0},{0,0,0}}, 
			(*find disorientation (minimum angle) using symmetry operators*)
				{MAngle, MMatrix} = getMisAngle[vecArray2[[i,j]], matrix];
				retval[[k, 1]] = MAngle/Degree;
				retval[[k, 2]] = (MMatrix[[2,3]] - MMatrix[[3,2]])/2*Csc[MAngle];
				retval[[k, 3]] = (MMatrix[[3,1]] - MMatrix[[1,3]])/2*Csc[MAngle];
				retval[[k, 4]] = (MMatrix[[1,2]] - MMatrix[[2,1]])/2*Csc[MAngle];
				If[retval[[k,4]]<0, retval[[k,2;;4]] = -retval[[k,2;;4]];];
			];
		,{k, 2}];
		];
Return[retval];
];
]


(* ::Subsection:: *)
(*Misorientation map/pole figure creator*)


misMap2Images[meapFile_, mcrit_, index_]:=Module[{array, misArrayLeft, misArrayUp, AngleLinesLeft, AngleLinesUp, bco, info, width, height, bcImage},
Catch[
	info = Import[meapFile, "tsv"];
	width = info[[1, 2]];
	height = info[[1, 4]];
	array = info[[3;;]];
	If[index<31,bcImage = Image[Partition[ToExpression[array[[;;, index]]], width]];];
	misArrayLeft = Select[ToExpression[array[[;;, {30, 29, 15}]]], #[[3]]>mcrit && #[[3]]<361&];
	misArrayUp = Select[ToExpression[array[[;;, {30, 29, 17}]]], #[[3]]>mcrit && #[[3]]<361&];
	AngleLinesLeft = Graphics[{Thick, Hue[(#[[3]])/120], Line[{{#[[2]]-1, height-#[[1]]}, {#[[2]]-1, height-#[[1]]+1}}] }]&/@misArrayLeft;
	AngleLinesUp = Graphics[{Thick, Hue[(#[[3]])/120], Line[{{#[[2]], height-#[[1]]+1}, {#[[2]]-1, height-#[[1]]+1}}] }]&/@misArrayUp;
	If[index<31, 
		bco = Show[{ImageAdjust@bcImage, AngleLinesLeft, AngleLinesUp}, ImageSize->Large];
		,
		bco = Show[{AngleLinesLeft, AngleLinesUp}, ImageSize->Large];
	];
	Return[bco];
]]


(* ::Text:: *)
(*misMapImages collects data from _List.txt and _band.txt files and creates images*)
(*Input: OIM ctf output file ctfFile, bool export, bool to show images showbulk, size imageSize (small, medium, large)*)
(*Output: Exports files if export=True, shows images if showbulk=true *)


misMapImages[ctfFile_, export_, showbulk_, imageSize_]:=Module[{str, dir, Lines, Lines2, Directions, Direction2D, bandImage, AngleHistogram, 
																	AxisAnglePlot, misorientationList, band},
		str = StringSplit[FileNameSplit[file][[-3]], " "][[1]];
		dir = FileNameJoin[{ParentDirectory[DirectoryName[ctfFile]], "Maps"}];
		If[!DirectoryQ[dir], CreateDirectory[dir]];

misorientationList = Import[FileNameJoin[{dir, StringJoin[{str, "_List.txt"}]}],  "tsv"];
band = Import[FileNameJoin[{dir, StringJoin[{str, "_band.txt"}]}], "tsv"];

If[export||showbulk, 
	Lines = Graphics[{Thick, Hue[#[[5]]/120], Line[{{#[[1]], #[[2]]}, {#[[3]], #[[4]]}}]}]&/@misorientationList;
	Lines2 = Graphics[{Thick, RGBColor @@ Abs[{#[[8]], #[[7]], #[[6]]}], Line[{{#[[1]], #[[2]]}, {#[[3]], #[[4]]}}]}]&/@misorientationList;
	Direction2D = Graphics[{Hue[#[[5]]/120], Point[{#[[6]], #[[7]]}/(Norm[#[[6;;8]]]+#[[8]])]}]&/@misorientationList;
	Direction2D = Show[Graphics[{EdgeForm[Thin], White, Disk[], Black, Line[{{-1,0},{1,0}}], Line[{{0, -1}, {0,1}}], 
						Text[Style["Misorientation Axis", Large], {0, 1.2}], Text[Style["x", Large], {1.2, 0}], Text[Style["y", Large], {0, -1.2}]}],  Direction2D];
	bandImage = Image[Transpose[(Max[band] - band)/(Max[band])]];
];

AngleHistogram = Histogram[misorientationList[[;;,5]],Ceiling[Max[misorientationList[[;;,5]]]]];
AxisAnglePlot = ListPlot[Transpose[{misorientationList[[;;,5]],ArcCos[#[[8]]/Norm[#[[6;;8]]]]/Degree&/@misorientationList}]];

If[export,
	Export[FileNameJoin[{dir, StringJoin[{str, "_angles.tiff"}]}], Show[{bandImage, Lines}, Background->White, ImageSize->Full]];
	Export[FileNameJoin[{dir, StringJoin[{str, "_axes.tiff"}]}], Show[{bandImage, Lines2}, Background->White, ImageSize->Full]];
	Export[FileNameJoin[{dir, StringJoin[{str, "_axes_angles.tiff"}]}], Show[Direction2D, Background->White, ImageSize->Full]];
	Export[FileNameJoin[{dir, StringJoin[{str, "_axes_angle_plot.tiff"}]}], Show[AxisAnglePlot, Background->White, ImageSize->Full]];
	Export[FileNameJoin[{dir, StringJoin[{str, "_angle_hist.tiff"}]}], Show[AngleHistogram, Background->White, ImageSize->Full]];
];
If[showbulk, 
	Print[Show[{bandImage, Lines}, Background->White, ImageSize->imageSize], 
		Show[{bandImage, Lines2}, Background->White, ImageSize->imageSize], 
		Show[Direction2D, Background->White, ImageSize->imageSize]
];
];
Print[Show[AxisAnglePlot, Background->White, ImageSize->imageSize], Show[AngleHistogram, Background->White, ImageSize->imageSize]];


]


(* ::Subsection:: *)
(*Agglomerated misorientation data*)


(* ::Text:: *)
(*misMapAgglom finds aggregated misorientation data for all samples*)
(*Input: list of .ctf files ctfFileList, bool hist, bool plot, bool pole, bool dpole (for density pole figure)*)
(*Output: prints requested images*)
(**)
(*misMapAgKent finds Kent distribution of misorientation axes for all samples*)
(*Input: list of ctf files ctfFileList*)
(*Output: table of Kent distribution statistics*)


misMapAgglom[ctfFileList_, hist_, plot_, pole_, dpole_]:=Module[{Direction2D, dim, Orientations, ctfFile, misorientationList, AngleHistogram, AxisAnglePlot, x,y, DensityPoleFigure},

Orientations = ConstantArray[0, Length[ctfFileList]];
Do[	
	ctfFile = ctfFileList[[i]];
	misorientationList = Import[ctfFile,  "tsv"];	
	Orientations[[i]] = misorientationList;
, {i, Length[ctfFileList]}];

Orientations = Flatten[Orientations, 1];
If[hist, 
	AngleHistogram = Histogram[Orientations[[;;,5]], Ceiling[Max[Orientations[[;;,5]]]], ChartStyle->{EdgeForm[Thick], White}];
	Print[Show[AngleHistogram, ImageSize->Large]];
];


If[plot, 
	AxisAnglePlot = ListPlot[Transpose[{Orientations[[;;,5]],ArcCos[#[[8]]/Norm[#[[6;;8]]]]/Degree&/@Orientations}]];
	Print[Show[AxisAnglePlot, ImageSize->Medium]];
];
If[pole, 
	Direction2D = Graphics[{Hue[#[[5]]/120], Point[{#[[6]], #[[7]]}/(Norm[#[[6;;8]]]+#[[8]])]}]&/@Orientations;
	Direction2D = Show[Graphics[{EdgeForm[Thin], White, Disk[], Black, Line[{{-1,0},{1,0}}], Line[{{0, -1}, {0,1}}], 
						Text[Style["Misorientation Axis", Large], {0, 1.2}], Text[Style["x", Large], {1.2, 0}], Text[Style["y", Large], {0, -1.2}]}], Direction2D];
	Print[Show[Direction2D, ImageSize->Medium]];
];
(*Density pole figure*)
If[dpole, 
	dim = 300;
	DensityPoleFigure = ConstantArray[0, {dim+1,dim+1}];
	Do[
		{x,y} = Round/@(dim/2*{Orientations[[i,6]], Orientations[[i,7]]}/(Norm[Orientations[[i,6;;8]]]+Orientations[[i,8]])) + {dim/2+1,dim/2+1};
		DensityPoleFigure[[x, y]] = DensityPoleFigure[[x,y]] + 1;
	,{i, Length[Orientations]}];
	DensityPoleFigure = DensityPoleFigure/Max[DensityPoleFigure];
	Print[Image[DensityPoleFigure]];
];
]


misMapAgKent[ctfFileList_]:=Module[{Direction2D, dim, Orientations, ctfFile, misorientationList, 
									mins, AngleHistogram, AxisAnglePlot, x,y, DensityPoleFigure},

Orientations = ConstantArray[0, Length[ctfFileList]];
Do[	
	ctfFile = ctfFileList[[i]];
	misorientationList = Import[ctfFile,  "tsv"];	
	Orientations[[i]] = misorientationList[[;;, 5;;8]];
, {i, Length[ctfFileList]}];
Orientations = Flatten[Orientations,1];
mins = {11, 27, 53, 63, 71, 89, 103};
kentresults = ConstantArray[0, {Length[mins]+1, 15}];
kentresults[[1,1]] = "Angle";
kentresults[[1, 2;;]] = kentSp[{{1,1,1},{1,0,1}}][[;;,1]];
Do[
	kentresults[[i+1, 1]] = mins[[i]]+1;
	kentresults[[i+1, 2;;]] = kentSp[Select[Orientations, #[[1]]>mins[[i]] && #[[1]]<mins[[i]]+2&][[;;, 2;;4]]][[;;,2]];
,{i, Length[mins]}];

kentresults
]


(* ::Subsection:: *)
(*Symmetry Operator Functions*)


(* ::Text:: *)
(*aragSym outputs symmetry operator matrix*)
(*Input: operator number num*)
(*Output: 4x4 symmetry operator matrix*)
(**)
(*orthoSym outputs symmetry operator matrix*)
(*Input: operator number num*)
(*Output: 3x3 symmetry operator matrix*)


aragSym[num_]:=Module[{h, O, b},
	(*applicable for #62 Pnma systems*)
h = 0.5;

Switch[num
,1,
	O = {{1,0,0,0},{0,1,0,0},{0,0,1,0}, {0,0,0,1}};
,2,
	O = {{1,0,0,0}, {h,-1,0,0},{0,0,-1,0},{h,0,0,1}};
,3,
	O = {{1,0,0,0},{0,-1,0,0},{h,0,1,0},{0,0,0,-1}};
,4,
	O = {{1,0,0,0},{h,1,0,0},{h,0,-1,0},{h,0,0,-1}};
,5,
	O = {{1,0,0,0},{0,-1,0,0},{0,0,-1,0},{0,0,0,-1}};
,6,
	O = {{1,0,0,0},{h,1,0,0},{0,0,1,0},{h,0,0,-1}};
,7,
	O = {{1,0,0,0},{0,1,0,0},{h,0,-1,0},{0,0,0,1}};
,8,
	O = {{1,0,0,0},{h,-1,0,0},{h,0,1,0},{h,0,0,1}};
];

MOut = O
];


orthoSym[num_]:=
	(*applicable for #47 Pmmm systems*)
Switch[num
,1,
	{{1,0,0},{0,1,0},{0,0,1}}
,2,
	{{-1,0,0},{0,-1,0},{0,0,1}}
,3,
	{{-1,0,0},{0,1,0},{0,0,-1}}
,4,
	{{1,0,0},{0,-1,0},{0,0,-1}}
,5,
	{{-1,0,0},{0,-1,0},{0,0,-1}}
,6,
	{{1,0,0},{0,1,0},{0,0,-1}}
,7,
	{{1,0,0},{0,-1,0},{0,0,1}}
,8,
	{{-1,0,0},{0,1,0},{0,0,1}}
]

