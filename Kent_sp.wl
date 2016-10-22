(* ::Package:: *)

(* ::Title:: *)
(*Kent Distribution Functions*)


(* ::Subsection:: *)
(*Last  updated 3/5/15 by Leanne Friedrich*)


(* ::Section:: *)
(*kentSp*)


(* ::Text:: *)
(*kentSp takes in a set of points in direction cosines and outputs statistics about their correlation*)
(*Kent's paper is J.R. Statist. Soc 1982, 44, No 1, pp. 71-80 *)
(*input: nx3 array of points in direction cosines*)
(*output: 2x14 array of statistics*)
(**)
(*adapted from SPAK*)


kentSp[v_]:=Module[{d,n,avp,a, b ,s, r, meandir ,rprime, 
			T, S, H, B, grph,x,output,d2,roated,di,mo,kb,
				eigenvals,psi, K, G, V, q, kappa, beta, 
				Kstat, rd, uv, sigv, ellz, ell, PP, CONFIDENCE, theta, phi},

Catch[
If[Length[v]<2,
	Throw[{{"number of points", 0},{"mean direction", 0},{"major axis", 0}, {"minor axis", 0}, 
 {"mean resultant length",0},  {"q", 0}, {"kappa", 0}, {"beta", 0},{"kappa/beta", 0}, {"Kstat", 0}, {"StDev (theta) (Deg)", 0}, 
{"StDev (phi) (Deg)", 0}, {"mode", 0}, {"dist", 0}}];
];
CONFIDENCE = 0.05;	(*% do our tests to 95% confidence*)

d = v;
n = Dimensions[d][[1]]; (*n is number of data points*)
meandir = Mean[v]; (*meandir is mean direction of v in direction cosines*)
avp = dc2tp[{meandir}]; (*avp is mean direction of v in radians*) 
theta = avp[[1,1]]; 
phi = avp[[1,2]];

(*% calc mean resultant length*)
s = Total[d];
r = Sqrt[Total[s * s]];
rprime = r / n  ;                      (*% mean resultant length*)

(*% compute Kent parameters*)
(*G columns: 1 is mean direction of distribution, 2 is major axis, 3 is minor axis*)
T = Table[Total[d[[;;, i]]*d[[;;,j]]], {i, 1, 3}, {j, 1, 3}];
T = T / n;
G = Transpose[Eigenvectors[T]];
If[(meandir.G)[[1]]<0, G = -G];
eigenvals = Eigenvalues[T];
(*eigenvalue 1 is an indication of the concentration of vectors*)

q = eigenvals[[2]]-eigenvals[[3]];(*q is the difference between the two lowest eigenvalues*)
(*kappa describes degree of concentration, beta is ovalness *)
kappa = (2 - 2 * rprime - q)^(-1) + (2 - 2 * rprime + q)^(-1);
beta = 0.5 * ((2 - 2 * rprime - q)^(-1) - (2 - 2 * rprime + q)^(-1));

(*% Report on a few things*)
kb = kappa/beta;
If [kb >= 2, 
	mo = "Unimodal";, 
	mo = "Multimodal";];


(*% is it Fisherian or Kent?*)

Kstat = n * ((0.5 * kappa) ^ 2) * BesselI[0.5, kappa]/BesselI[2.5, kappa] * q * q;
If [Kstat > (-2 * Log[CONFIDENCE]), 
	di = "Kent";,
	di = "Fisher"];
d2 = d.G; (*rotates system so long/short axes are aligned with lab axes*)
Do[
	If[d2[[i,1]]<0,
		d2[[i]] = -d2[[i]];
	];
,{i, Length[d2]}];
(*Print[Graphics3D[{Sphere[], Point[#]&/@d2}]];*)

sigv = StandardDeviation[dc2tp[d2]]; (*standard deviation in radians*)

output = {{"number of points", n},{"mean direction", meandir},
{"major axis", G[[;;, 3]]}, {"minor axis", G[[;;, 2]]}, 
 {"mean resultant length",rprime},  {"q", q}, {"kappa", kappa}, {"beta", beta},
{"kappa/beta", kb}, {"Kstat", Kstat}, {"StDev (theta) (Deg)", sigv[[1]]/Degree}, 
{"StDev (phi) (Deg)", sigv[[2]]/Degree}, {"mode", mo}, {"dist", di}};
Return[output];
]
];



(* ::Text:: *)
(*ellSpPara creates a graphics object representing the Kent distribution*)
(*Input: 3d vector mean, scalar majr in Degrees, scalar minr in Degrees, vector maj, vector min*)
(*output: graphics object ellipse on sphere projected onto plane*)


ellSpPara[mean_, majr_, minr_, maj_, min_]:=Module[{height, ang, rotmat, rotmat2, eq, retval, t, major, minor, retval2}, 
Catch[
	major = majr*Degree;
	minor = minr*Degree;
	ang = VectorAngle[{1,0,0}, maj];
	rotmat = RotationMatrix[{{0,0,1}, mean}];
	rotmat2 = RotationMatrix[ang/2, mean];
	eq = {major*Cos[t], minor*Sin[t], Sqrt[1-(major*Cos[t])^2-(minor*Sin[t])^2]};
	eq = rotmat2.rotmat.eq;
(*	retval = {ParametricPlot3D[eq, {t,0,360}, Boxed\[Rule]False](*, 
			Graphics3D[{Line[{{0,0,0}, mean}], Line[{mean, mean+maj*0.5}], Line[{mean, mean+min*0.2}]}]*)};*)
	retval2 = ParametricPlot[{eq[[1]], eq[[2]]}/(Norm[eq]+eq[[3]]), {t, 0, 360}];
	Return[retval2];
]];


(*function p = dc2tp(d)
% converts between direction cosines and theta phi coords in radians*)
dc2tp[d_]:=Module[{p, x}, 
Catch[
	p = {If[#[[1]]==0&&#[[2]]==0, 0, ArcTan[#[[1]], #[[2]]]], ArcCos[#[[3]]]}&/@d;
	p[[;;, 2]] = p[[;;, 2]]/.x_/;x<0->x+2*Pi;
	Return[p];
]]


(* ::Section:: *)
(*Cluster selection*)


(* ::Text:: *)
(*gather takes all data within r degrees of mean {theta, phi}*)
(*Input: angle theta in degrees, angle phi in degrees, angle r in degrees, array dat*)
(*Output: array of data*)


(*takes all data within phi degrees of mean {theta, phi}*)
gather[theta_,phi_,r_, dat_]:=Module[{k},(
k=dat;
Select[k, (VectorAngle[#, hp2cart[tp2hp[{{theta,phi}}*Degree]][[1]]])/Degree<=r&]
)]


(* ::Text:: *)
(*itgather takes all data within r degrees of the centroid of the dominant cluster. It finds the centroid of the set, then takes all data within r degrees of the centroid, then recalculates the centroid and takes all data within r degrees of the centroid. It iterates until the number of elements converges.*)
(*Input: angle r in degrees, array dat, which should be a set of cartesian coordinates*)
(*Output: array of cartesian vectors*)


itgather[r_, dat_]:=Module[{vectors, numelements, numelements2, centroid, r2, vectors2},
	r2 = r*Degree;
	vectors = dat;
	vectors2 = vectors;
	numelements = 2;
	numelements2 = 1;
	While[numelements2!=numelements,
		numelements = Length[vectors2];
		centroid = Mean[vectors2];
		vectors2 = Select[vectors, VectorAngle[centroid, #]<r2&];
		numelements2 = Length[vectors2];
	];
	vectors2
]


(* ::Text:: *)
(*angSel creates a module that selects data within angle x from center m*)
(*Input: n x 9 array of locations and angles dat*)
(*Output: interface that allows the user to manually select a region on a sphere*)


angSel[dat_]:=Module[{vecx, vecplot, t, p, ra, retval},
	vecx = dat;
	vecplot = Graphics3D[{White, Point[#]&/@vecx}];
	pointsx = Show[Graphics3D[{Blue, Sphere[{0,0,0}]}], vecplot];
{t, p, ra} = DialogInput[Manipulate[
		cone = Graphics3D[Cone[{hp2cart[tp2hp[{{theta,phi}}*Degree]][[1]], {0,0,0}}, Sin[r*Degree]]];
		Column[{
			Show[pointsx, 
				(*Graphics3D[{PointSize[.05], Red, Point[hp2cart[tp2hp[{{theta,phi}}*Degree]][[1]]]}], 
				Graphics3D[Cone[{{0,0,Cos[theta*Degree]}, {0,0,0}}, Abs[Sin[theta*Degree]]+.01]],
				Graphics3D[Cone[{hp2cart[tp2hp[{{theta,phi}}*Degree]][[1]], {0,0,0}}, Sin[r*Degree]]]*)
				cone, 
				ViewVertical->{0,0,1}, ViewPoint->{1,1,1},
				Boxed->False, PlotRange->{{0,1}, {0,1}, {0,1}}, ImageSize->Medium
			]
			, 
			Button["Save selection", DialogReturn[{theta, phi, r}];],
			Button["Skip sample", DialogReturn[{300,300,300}];],
			Button["Exit", DialogReturn[{500,500,500}];]
		}]
	, {theta, 0, 90,5, Appearance->"Open", ControlPlacement->Bottom}
	, {phi, 0, 90,5,Appearance->"Open", ControlPlacement->Bottom}
	, {r,5,45,5, Appearance->"Open", ControlPlacement->Bottom}
	, LocalizeVariables->False
	]
];
	If[{t,p,ra} == {500,500,500}, Return[5]];
	If[{t,p,ra} == {300,300,300}, Return[{0,0}]];
	retval = gather[t,p,ra, vecx];
	Return[retval];

]


(* ::Text:: *)
(*chooseClusters creates a module that automatically selects clusters for a given number of clusters*)
(*Input: n x 9 array of locations and angles dat*)
(*Output: interface that allows the user to choose how many clusters to use*)


chooseClusters[dat_]:=DynamicModule[{vecplot, retval, vecx},
Catch[
	vecx = dat;
	vecplot = Graphics3D[{White, Point[#]&/@vecx}];
(*	sphere = Show[Graphics3D[{Blue, Opacity[.1], Sphere[{0,0,0}]}], Boxed-> False];*)
	retval = DialogInput[
		Manipulate[
			vecx2 = FindClusters[vecx,num, Method-> "Agglomerate"]; 
			allvecs = plotvecs[Hue[#/num],vecx2[[#]]]&/@ Table[n, {n,1,num}]; 
			Column[{Show[(*sphere, *)allvecs, ViewVertical->{0,0,1}, ViewPoint->{1,1,1}, ImageSize->Medium, Boxed->False, PlotRange->{{0,1}, {0,1}, {0,1}}], 
					Button["Save selection", DialogReturn[{vecx2[[1]], num}];],
					Button["Skip sample", DialogReturn[{0,0}];],
					Button["Exit and lose selections", DialogReturn[5];]}]
		, {{num, 8}, 1, Min[30, Length[vecx]], 1, Appearance->"Open", ControlPlacement->Bottom}

		]
	];
	Return[retval];
]
]


(* ::Section:: *)
(*Kent Distribution Subfunctions*)


mismax[files_]:=Module[{hist, mat},(
mat = ConstantArray[0, {Dimensions[files][[1]],2}];
Do[
	hist = Import[files[[i]], "tsv"][[2;;, ;;]];
	mat[[i,1]]=DirectoryName[files[[i]]];
	mat[[i,2]] = hist[[Ordering[hist[[;;,2]],-1],1]][[1]];
,{i, Dimensions[mat][[1]]}];
mat
)];


Keval[n_, v_, theta1_, theta2_, theta3_, phi1_, phi2_, phi3_, r1_, r2_, r3_, dir_]:=
Module[{red1, red2, red3, pics, nums}, (
	red1=Prepend[kentSp[gather[theta1, phi1, r1, v]], "peak 1"];
	red2=Prepend[kentSp[gather[theta2, phi2, r2, v]], "peak 2"];
	red3=Prepend[kentSp[gather[theta3, phi3, r3, v]], "peak 3"];
	nums = Join[{n}, red1[[1;;15]], red2[[1;;15]], red3[[1;;15]]])
];


autoKeval[n_, v_, numel_]:=Module[{v2, nums},
	nums = {n}; 
	Do[
		If[Dimensions[v[[i]]]>.1*numel, Join[nums, kentSp[v[[i]]]];]
	, {i, Dimensions[v][[1]]}]
	nums
]


(* ::Subsection:: *)
(*Kent Interface.nb*)


(* ::Text:: *)
(*initializeKentFile creates a starter file for Kent distributions*)
(*Input: Array of grain files (just need the length) grainfiles, directory string nb, method string str*)
(*Output: Exports a file if it doesn't exist already*)


initializeKentFile[grainfiles_, nb_, str_]:=Module[{kentresults},
	kentresults = ConstantArray[0, {Dimensions[grainfiles][[1]]+1, 15}]; 
	kentresults[[1]]={"File", "Number of points", "mean direction",
					"major axis","minor axis","mean resultant length",
					"q","kappa","beta","kappa/beta","Kstat",
					"StDev (long)(Deg)","StDev (short) (Deg)","mode","dist"};
	If[Length[FileNames["*Kent Distributions4.txt"]]==0, 
		Export[FileNameJoin[{nb, 
						StringJoin["/Summaries/summary4/Kent Distributions_", str, ".txt"]}], 
							kentresults, "tsv"],
		Print["File already exists"]];
]


(* ::Text:: *)
(*grainKent finds kent distributions of <001> axes*)
(*Input: array of files from getFiles, notebook directory name, bool to export files, string mode, bool lens (true if you're just looking at the lens, false if you're not looking at the lens)*)
(*Output: table of statistics*)


findKentsC[files_, nb_, export_, mode_, lens_]:=DynamicModule[{dir, bitmap, caxes, caxes2, modestr, kentresults},
Catch[
	If[!(mode == "grain" || mode == "all" || mode =="5deg" || mode=="cluster" || mode=="manual"),
		Print["Mode must be one of {'grain', 'all', '5deg', 'cluster', 'manual'}"];
		Throw[0];
	];
	If[lens, modestr = StringJoin[mode, "_lens"];, modestr=StringJoin[mode, "_shell"];];
	If[Length[FileNames[StringJoin["*Kent Distributions_", modestr, ".txt"], "*", Infinity]]==0, 
		initializeKentFile[files[[;;,4]], nb, modestr];
	];
	kentresults = Import[FileNames[StringJoin["*Kent Distributions_", modestr, ".txt"], "*", Infinity][[1]], "tsv"];
	CreateDocument[{
		Button["Save Results", Export[FileNameJoin[{nb, 
			StringJoin["Summaries/summary5/Kent Distributions_", modestr, ".txt"]}], 
			kentresults, "tsv"];
		Export[FileNameJoin[{nb, 
			StringJoin["Summaries/summary5/Kent Distributions_", modestr, ".xlsx"]}], 
			kentresults];], 
		Dynamic[Grid[kentresults]]}, 
		
	WindowTitle->StringJoin["Kent results for ", modestr], WindowSize -> {Scaled[1], Scaled[1]}];

	If[mode=="cluster", 
		kentresults = Transpose[Append[Transpose[kentresults], ConstantArray[0, Length[kentresults]]]];
		kentresults[[1,-1]] = "# of clusters"];

(*-------------collect data---------------------*)
(*universally, caxes = {0,0} is a failure condition and indicates that the program should not proceed with the kent distribution*)
	Do[
		kentresults[[i+1]] = findKentsCSingle[mode, lens, files[[i]], False];
	, {i, Length[files]}];


(*--------------export results----------------*)
	If[export,
		Export[FileNameJoin[{nb, 
			StringJoin["Summaries/summary5/Kent Distributions_", modestr, ".txt"]}], 
			kentresults, "tsv"];
		Export[FileNameJoin[{nb, 
			StringJoin["Summaries/summary5/Kent Distributions_", modestr, ".xlsx"]}], 
			kentresults];
	];
	Return[kentresults];
];
]


(* ::Text:: *)
(*findKentsCSingle finds kent distributions for a single sample*)
(*Input: string mode, bool lens, 1x4 vector of strings file, bool printf*)
(*Output: if printf, an array of kent stats and a panel of diagnostics. else, an array of kent stats*)


findKentsCSingle[mode_, lens_, file_, printf_]:=Module[{bitmap, caxes, dir, row, x, poleF, poleFContour},
Catch[
		dir = FileNameSplit[file[[4]]][[-3]]; 
		row = ConstantArray[0, 15 + If[mode=="cluster", 1, 0]];
		Switch[mode,
		"grain",  (*--------------single grain ---------------*)
			bitmap = largeGrainSelect[{file[[1]], file[[2]]}, 1*Degree, False, False, False, lens];
			If[bitmap!={{0}},
				(*grain selection succeeded - continue*)		
				bitmap = Transpose[ConstantArray[bitmap, 3], {3,1,2}];	
				caxes = rearrangeCTF[file[[2]]][[1, 3;;-2, 3;;-2]];
				caxes = caxes[[;;, ;;, ;;, 3]]*bitmap;
				caxes = Flatten[caxes,1];
				caxes = Select[caxes, #!={0,0,0}&];
				,
				caxes = {0,0};
			];

		, "all", (*--------------all points in lens ---------------*)
			caxes = adjustSizes[file[[4]], lens];
			If[caxes!={0,0},
				caxes = caxes[[1, ;;, -4;;-2]];
				caxes = euler2c/@caxes;
				caxes = Abs[caxes];
			];

		, "5deg", (*--------------only points within 5 degrees of the centroid---------------*)
			caxes = adjustSizes[file[[4]], lens];
			If[caxes!={0,0},
				caxes = caxes[[1, ;;, -4;;-2]];
				caxes = euler2c/@caxes;
				caxes = Abs[caxes];
				caxes = itgather[5, caxes];
			];

		, "cluster", (*--------------automatically choose cluster of points---------------*)
			caxes = adjustSizes[file[[4]], lens];
			If[caxes!={0,0},
				caxes = caxes[[1, ;;, -4;;-2]];
				caxes = euler2c/@caxes;
				caxes = Abs[caxes];
				caxes = chooseClusters[caxes];
				If[caxes==5, Print["Function aborted by user"]; Throw[0]];
				row[[-1]] = caxes[[2]]; (*extract number of clusters*)
				caxes = caxes[[1]];
				,
				Print["No points detected for ", dir];
			];	

		, "manual", (*--------------manually choose points within a region---------------*)
			caxes = adjustSizes[file[[4]], lens];
			If[caxes!={0,0},
				caxes = caxes[[1, ;;, -4;;-2]];
				caxes = euler2c/@caxes;
				caxes = Abs[caxes];
				caxes = angSel[caxes];
				If[caxes==5, Print["Function aborted by user"]; Throw[0]];
				,
				Print["No points detected for ", dir];
			];	
		];

		If[caxes!={0,0}, 
			caxes = Abs[caxes]; (*only applies for orthorhombic systems*)
			x = kentSp[caxes];
			If[mode!="cluster", 
				row[[2;;]] = x[[;;,2]]; (*take second row of kentSp output b/c it outputs headers*)
				,
				row[[2;;-2]] = x[[;;,2]]; (*take second row of kentSp output b/c it outputs headers*)
			];
			If[printf, 
				poleF = Show[{pf[caxes, False, dir, "<001>", False, False, True], 
								ellSpPara[x[[2,2]], x[[11,2]], x[[12,2]], x[[3,2]], x[[4,2]]]}];
				poleFContour = pf[caxes, False, dir, "<001>", False, True, True];
				Return[{row, Panel[Row[{Column[{Style[dir, Large], Grid[x]}], poleF, poleFContour}]]}];
				,
				Return[row];
			];
		];
		If[printf, 
			Return[{row, Panel["No points selected"]}];
			,
			Return[row];
		];
]]


(* ::Text:: *)
(*rotateXSc rotates c axes of cross-sections, so the lab axes of cross-sections align with the lab axes of plan sections*)
(*Input: array of c axis vectors cdirections*)
(*Output: array of angles between c axis vectors and <001>*)


rotateXSc[cdirections_]:=Module[{dirs, dirsrots, rotatedxs},
	dirs = PadRight[(Import[FileNameJoin[
								{nb, "Summaries/summary/corneal directions.txt"}
								], "tsv"]), {7,3}];
	dirsrots = (RotationMatrix[{#,{0,0,1}}])&/@dirs;
	rotatedxs = ConstantArray[{0,0,0}, 7];
	Do[
		rotatedxs[[i]] = Min[VectorAngle[{0,0,1},dirsrots[[i]].(-cdirections[[18+i]])],VectorAngle[{0,0,1},dirsrots[[i]].(cdirections[[18+i]])]]
	,{i, 7}];
	rotatedxs = rotatedxs/Degree
]


(* ::Text:: *)
(*angleGridPre prints out a summary table for angles between the c axis and the lab <001>*)
(*Input: array of scalars cangles*)
(*Output: Prints a grid*)


angleGridPre[cangles_]:=Module[{},
	Print[Grid[{{"", "Mean", "STDEV", "Min", "Max", "Median"},
	{"Plan", Mean[cangles[[1;;17]]],
	StandardDeviation[cangles[[1;;17]]],
	Min[cangles[[1;;17]]],
	Max[cangles[[1;;17]]],
	Median[cangles[[1;;17]]]
	}
	,
	{"XS", Mean[cangles[[19;;25]]],
	StandardDeviation[cangles[[19;;25]]],
	Min[cangles[[19;;25]]],
	Max[cangles[[19;;25]]],
	Median[cangles[[19;;25]]]}
	,
	{"All", Mean[cangles],
	StandardDeviation[cangles],
	Min[cangles],
	Max[cangles],
	Median[cangles]}}]]
]


(* ::Text:: *)
(*angleGridPost prints out a summary table for angles between the c axis and the lab <001>*)
(*Input: array of scalars cangles and array of scalars rotatedxs*)
(*Output: Prints a grid*)


angleGridPost[cangles_, rotatedxs_]:=Module[{}, 
Print[Grid[Transpose[Join[{{"", "Plan", "XS", "All"}}, 
	Transpose[Join[{{ "Mean", "STDEV", "Min", "Max", "Median"}},
					{Mean[#], StandardDeviation[#], Min[#], Max[#], Median[#]}
					&/@{cangles[[1;;17]], rotatedxs, 
					Join[rotatedxs, cangles[[1;;17]]]}]]]]] ]
]


(* ::Subsection:: *)
(*Selecting largest grain*)


(* ::Text:: *)
(*adjustSizes fixes an overflow error in grain sizes*)
(*Input: "lens.txt" file, bool lens (true if only lens, false if only not lens)*)
(*Output: {table of euler angles, step size}*)


adjustSizes[gfile_, lens_]:=Module[{angles, pix, gstep,  ratio, rat, exponent, selection},
	Catch[
		If[lens, selection = "True", selection = "False"];
		angles = Import[gfile, "tsv"][[2;;, 2;;]];
		angles = Select[angles, ToString[#[[1]]]=="aragonite-Pmcn"&]; (*only select aragonite*)
		If[Length[angles]<1, Return[{0, 0}];];
		angles = Select[angles, #[[13]] == selection&]; (*select only files in lens*)
		If[Length[angles]<1, Return[{0, 0}];];
		angles = angles[[;;, 2;;]]; 
		angles = ToExpression/@angles;

		(*an error in hkl's software causes it to export any number above 6 chars 
		(inc. decimal) in sci note, but it doesn't include the exponent because it 
		only exports 6 chars *)
		pix = Commonest[angles[[;;,1]]][[1]]; (*under sci note failure, this would be the pixel size*)
		gstep = Min[Select[angles[[;;, 1]], #>0&]]; (*under normal conditions, this would be the pixel size*)
		If[pix<10, pix = gstep;];  (*the sci note error happens with very large scan areas. If this is the case, gstep will fail, and pix will be larger than gstep, which must be less than 10*)
	
		ratio = angles[[2,1]]/angles[[2,2]]^2;
		Do[
			If[angles[[i,1]]<pix,
				rat = angles[[i,1]]/angles[[i,2]]^2;
				exponent = Round[Log[10, ratio/rat]];
				angles[[i,1]] = angles[[i,1]]*10^exponent;
			]
		,{i,Length[angles]}];
	Return[{angles, pix}];

	];
];


(* ::Text:: *)
(*findLargestGrain determines the center and size of the largest grain found by hkl software*)
(*Input: grain sizes .txt file given by Tango, bool lens (true if only lens, false if only shell/cornea)*)
(*Output: scalars {size, x center, y center, step size}*)


findLargestGrain[gfile_, lens_]:=Module[{angles, ratio, rat, exponent, 
								gstep, pix, xcg, ycg, eSize},
Catch[
(*-------------Part 1: Determine Largest Grain--------------*)
	{angles, pix} = adjustSizes[gfile, lens];
	If[ToExpression[pix]==0
		,
		(*grain selection failed*)
		Return[{0,0,0,0}];
		,
		angles = Sort[angles, #1[[1]]>#2[[1]]&];
		xcg = angles[[1,4]];
		ycg = angles[[1,5]];
		eSize = angles[[1,1]]/pix;
		Return[{eSize, xcg, ycg, Sqrt[pix]}];
	];
];
]


(* ::Text:: *)
(*rearrangeCTF creates an array that represents a map of the scan*)
(*Input: .ctf file given by Aztec*)
(*Output: {array of 3x3 arrays, scalar step size}*)


rearrangeCTF[ctfFile_]:=Module[{euler, eulers, minmin, maxx, maxy, output, sRow, pNum},
Catch[	
	euler = Import[ctfFile, "tsv"][[14,3]];
	If[euler=="aragonite-Pmcn",
		sRow = 16;
		pNum = 1;
		,
		sRow = 17;
		pNum = 2;
	];

	euler = Import[ctfFile, "tsv"][[sRow;;,{1,2,3,4,6,7,8}]];  (*phase, x, y, bands, euler1, euler2, euler3*)
	euler = Select[ToExpression/@euler, #[[1]]==pNum&][[;;, 2;;]]; (*only select aragonite*)
	(*find rotation of axes: xvecs is the rotated <100>, etc., ref. Euler angle wiki*)
	If[Length[euler]==0, Print["ctf file failed"]; Throw[{}];];
	(*eulers is an array of cartesian vectors*)
	eulers = ConstantArray[0, {Length[euler], 2}];
	Do[
		eulers[[i,1]] = euler[[i, 1;;3]];
		eulers[[i,2]] = euler2vecs[euler[[i,4]], euler[[i,5]], euler[[i,6]]];
	,{i, Length[euler]}];

	(*convert x and y indices to integers*)
	minmin = Import[ctfFile, "tsv"][[7,2]];

	eulers[[;;, 1, 1]] = Round[eulers[[;;,1,1]]/minmin]+2;
	eulers[[;;,1,2]] = Round[eulers[[;;,1,2]]/minmin]+2;
	
	(*rearrange eulers into an array that reflects locations*)
	maxx = Max[eulers[[;;,1,1]]];
	maxy = Max[eulers[[;;,1,2]]];

	output = ConstantArray[{{0,0,0},{0,0,0},{0,0,0}}, {maxx+2, maxy+2}];
	Do[
		If[eulers[[i, 1, 3]]>0,
			output[[eulers[[i, 1, 1]]+1, eulers[[i, 1, 2]]+1]] = eulers[[i,2]];		
		]
	,{i, Length[eulers]}];
	Return[{output, minmin}];
]

]


(* ::Text:: *)
(*largeGrainSelect finds a bitmap of the largest grain in the scan*)
(*Input: {grain size file given by Tango, ctf file given by Aztec, bmp image ipfz}, critical angle in radians, bool corners (true if corners are neighbors), bool watch (true if you want to watch the grain be selected), bool printf (true if you want the function to be able to print information), bool lens (true if you only want to select from the lens) *)
(*Output: bitmap*)


largeGrainSelect[files_, angCrit_, corners_, watch_, printf_, lens_]:=Module[
		{eSize, xcg, ycg, gstep, output, minmin, adds, tempctr, tempctr2, checkQueue, im, checked,
		loopCounter, nextPix, nextPixM, thisPix, thisPixM, MAngle, gfile, ctfFile, tempctrLimit, temp},

Catch[

If[Length[files]>2,
	ipfz = files[[3]];
	im = Import[ipfz];
,
	im = Image[ConstantArray[0, {100, 100}]];
];
gfile = files[[1]];
ctfFile = files[[2]];



If[printf, Print["Selecting largest grain for file ", gfile]];

{eSize, xcg, ycg, gstep} = findLargestGrain[gfile, lens];
If[eSize==0, 
	If[printf, Print["ERROR: Largest grain failed."];]; 
	Return[{{0}}];
];

(*-----------Part 2: Determine points in largest grain------------*)

{output, minmin} = rearrangeCTF[ctfFile];
If[(minmin/gstep>1.5 || gstep<minmin/1.5),
	If[printf, Print["ERROR: Step size mismatch."];];
	Throw[ConstantArray[0, Dimensions[output][[1;;2]]]];
];
xcg = Round[xcg*minmin/gstep];
ycg = Round[ycg*minmin/gstep];

(*----------------Accumulate points in grain--------------------*)
If[!corners,
	adds = {{1,0}, {0,1}, {-1, 0}, {0, -1}};
	,
	adds = {{1,0}, {0,1}, {-1, 0}, {0, -1}, {1,1}, {1,-1}, {-1,1}, {-1,-1}};
];




(*make sure the middle point is valid*)
tempctr = 1;
tempctr2 = 1;
tempctrLimit = Min[Dimensions[output][[1]]-xcg, xcg, Dimensions[output][[2]]-ycg, ycg];

While[output[[xcg, ycg]]=={{0,0,0},{0,0,0},{0,0,0}} && tempctr2<tempctrLimit,
	While[output[[xcg, ycg]]=={{0,0,0},{0,0,0},{0,0,0}} && tempctr<Length[adds],
		{xcg, ycg} = {xcg, ycg} + tempctr2*adds[[tempctr]];
		tempctr = tempctr+1;
	];
	tempctr2 = tempctr2+1;
];
Clear[tempctr, tempctr2];
If[tempctr2==tempctrLimit,
	Print["Center is nowhere near grain."];
	Return[{{0}}];
];

(*make lists*)
checkQueue = {{xcg, ycg}};  (*check is the list of points on the edge still to be checked*)
checked = ConstantArray[0, Dimensions[output][[1;;2]]]; (*checked keeps track of checked pixels*)
checked[[xcg, ycg]] = 1; (*set checked to 1 when the pixel has been checked*)
inGrain = checked; (*inGrain is the bitmap of points in the grain*)


If[watch && printf,
imDim = ImageDimensions[im];
While[imDim[[1]]>400, imDim = imDim/2];
	temp=PrintTemporary[Row[{Show[im, ImageSize->imDim], Dynamic[Image[Transpose[inGrain[[3;;-2, 3;;-2]]], ImageSize->imDim]], "   Size: ", Dynamic[Total[inGrain,2]], 
		". Expected: ", eSize}]];
];

loopCounter = 0; (*prevents infinite loops*)


While[Length[checkQueue]>0 && loopCounter<1000,
	(*for each pixel in stack*)
	loopCounter = loopCounter+1;
	Do[
		(*pop the top pixel off of the stack*)
		thisPix = checkQueue[[1]];
		thisPixM = output[[thisPix[[1]], thisPix[[2]]]];
		(*for each surrounding pixel*)
		Do[
			nextPix = thisPix + adds[[j]];
			If[checked[[nextPix[[1]], nextPix[[2]]]]==0, 
				nextPixM = output[[nextPix[[1]], nextPix[[2]]]];
			
			(*if next pixel is not already checked, and it is indexed*)
			If[nextPixM!={{0,0,0},{0,0,0},{0,0,0}}, 
			(*find disorientation (minimum angle) using symmetry operators*)
				MAngle = getMisAngle[thisPixM, nextPixM];

				(*Print["Checking point ", nextPix, " around ", thisPix, " angle ", MAngle/Degree];*)
				(*if the angle is small enough, this point is in the grain. add it to the queue*)
				If[MAngle<angCrit,
					checkQueue = Append[checkQueue, nextPix];
					inGrain[[nextPix[[1]], nextPix[[2]]]] =1 ;
				];	
			];
			(*you've checked this pixel, so don't check it again*)
			checked[[nextPix[[1]], nextPix[[2]]]] = 1;
			];
		,{j, Length[adds]}];
		If[Length[checkQueue]>1, 
			checkQueue = checkQueue[[2;;, ;;]];
			,
			checkQueue = {};
		]
	,{i, Length[checkQueue]}]
]
If[printf,
	NotebookDelete[temp];
	Print[Row[{Show[im, ImageSize->imDim], Image[Transpose[inGrain[[3;;-2, 3;;-2]]], ImageSize->imDim], "   Size: ", Total[inGrain,2], 
		". Expected: ", eSize}]];
];
If[loopCounter==1000, Print["Infinite loop reached in grain selection."];];
Return[inGrain[[3;;-2, 3;;-2]]];
];

]
