(* ::Package:: *)

(* ::Subtitle:: *)
(*Default Settings and Utilities*)


(* ::Section::GrayLevel[0]:: *)
(*1	General setup*)


(* ::Subsection::GrayLevel[0]:: *)
(*1.1	Math utilities: rescale, unflatten, FCflatten, FCunflatten, findContour*)


$HistoryLength=3;
rescale[x_,{x1_,x2_},{y1_,y2_}]:=y1+((y2-y1) (x-x1))/(x2-x1);


(* ::Text:: *)
(*unflatten is the inverse of the standard Mathematica  function Flatten.*)
(*The following functions: FC = "Fortran convention (first index varies fastest)".*)


unflatten[e_,{d__?((IntegerQ[#]&&Positive[#])&)}]:= Fold[Partition,e,Take[{d},{-1,2,-1}]] /;(Length[e]===Times[d]);
FCflatten[tensor_]:=tensor//Transpose[#,Range[ArrayDepth@tensor,1,-1]]&//Flatten;
FCunflatten[array_,cmaxd_]:=array//unflatten[#,cmaxd//Reverse]&//Transpose[#,Range[Length@cmaxd,1,-1]]&;


(* ::Text:: *)
(*The following function extracts a contour from a contour plot and returns it as a depth-three array, that is, a List of Lists of Lists, that is, a SET of LINES that consist of POINTS which are defined as ORDERED PAIRS.  Using "// Line // Graphics" or "ListPlot" restores the graphical representation.  If the contour is a single curve, the first dimension of the array will be 1.*)


findContour[xyzData_List,contourLevel_]:=Module[{myGraphicsComplex,myBoundary},myGraphicsComplex=ListContourPlot[xyzData,Contours->{contourLevel},ContourShading->False,Frame->False,Axes->False];myBoundary=Reverse[Cases[myGraphicsComplex[[1,2]],_Line,\[Infinity]]]/. {i_Integer:>myGraphicsComplex[[1,1,i]]}/. Line[stuff_]->stuff;
Return@If[myBoundary=={},{{{0,0}}},myBoundary];  (*don't feed empty list to ListPlot *)
];


(* ::Text:: *)
(*It would be more efficient to find contours for many values of U at once, but that is too dangerous, because ListContourPlot outputs the contours in an unpredictable order, which makes it impossible to parse the result mess automatically.  The current depth-3 array is bad enough.*)


logspace[{xmin_,xmax_},nxmax_]:=Exp@Rescale[Range[nxmax],{1,nxmax},{Log@xmin,Log@xmax}];
Attributes[myDensityPlot]=HoldAll;


(* ::Subsection::GrayLevel[0]::Closed:: *)
(*1.2	String utilities: dump, toString, substVars*)


dump[vals_]:=StringJoin[Riffle[
Function[expr,ToString@expr<>"="<>ToString@Evaluate@ToExpression@expr]/@vals,", "]];

toString[x_]:=x//Map[ToString@NumberForm[#,8]&,#,{-1}]&//List//Flatten//List//ExportString[#,"Table","FieldSeparators"->" "]&;
toString::usage="toString[x] returns a string representation of x.  x can be a number, vector, or matrix.
Examples:
0.3 // ToExpression // toString
(0.1*3) // ToExpression // toString
{0.1, 0.2, 0.1*3} // ToExpression // toString
{{0.1, 0.2}, {0.1*3, 0.2*4}} // ToExpression // toString
";

substVars[s_String]:=StringReplace[s,{
RegularExpression@"\\$(\\w+)":>toString@ToExpression@"$1",
RegularExpression@"\\$\\((\\w+)\\)":>toString@ToExpression@"$1"
}];
substVars2[s_String]:=s//substVars//StringReplace[#," "->"_"]&;
substVars::usage="substVars[s] substitutes every occurrence of $var or $(var) in a string 
by the current value of the variable var in the current kernel.
Lists of numbers are output as space-delimited lists.

Warning: works with local variables in a Block, but not with local variables in a Module.";
substVars2::usage="substVars2[s] substitutes every occurrence of $var or $(var) in a string 
by the current value of the variable var in the current kernel.
Lists of numbers are output as underscore-delimited lists.

Warning: works with local variables in a Block, but not with local variables in a Module.";


(* ::Subsection::GrayLevel[0]::Closed:: *)
(*1.3	File utilities: trackInTerminal, myStringWrite, myBinaryWrite, myAsciiWrite*)


(* ::Text::GrayLevel[0]:: *)
(*trackInTerminal["file"] launches a terminal window that tracks the contents of "file" using the UNIX command "tail -f file".*)


trackInTerminal[fileName_]:=
Run["xterm -geometry 132x60 -e '
echo Press Ctrl-C to close this window.
tail --lines=+0  --retry -f "<>fileName<>" ' &"];


myStringWrite[file_,string_]:=Module[{stream},stream=OpenWrite[file];WriteString[stream,string<>"\n"];
Close[stream];];
myBinaryWrite[file_,data_,type_]:=Module[{stream=OpenWrite[file,BinaryFormat->True]},
BinaryWrite[stream,data,type];Close[stream];
];
myAsciiWrite[file_,list_List]:=myStringWrite[file,#]&@ExportString[list,"Table"];
Attributes[clobberSave]={HoldRest};
clobberSave[file_,stuff__]:=(If[FileExistsQ@file,DeleteFile@file];Save[file,stuff]);


(* ::Subsection::GrayLevel[0]::Closed:: *)
(*1.4	Graphics utilities: Solid, Thickness1, DarkGreen, color schemes, arrayPlot, colorBar*)


Solid=Dashing@{};
Thickness1=AbsoluteThickness@1;
DarkGreen=Darker[Green,0.67];
resistorColorList={Black,Brown,Red,Orange,Darker@Yellow,Darker@Green,Blue,Purple,Gray,Cyan};
resistorStyle[n_Integer]:={Thick,resistorColorList[[Mod[n,Length@resistorColorList,1]]]};
alternateColorList={Red,Darker@Green,Blue,Darker@Cyan,Magenta,Darker@Yellow,Gray,Black};
alternateStyles=({Thick,#}&/@alternateColorList);
blackStyles={{Thick,Black},{Dotted,Gray}};
myStyles=({Thick,#}&/@resistorColorList);
stemStyle={PlotMarkers->Automatic,PlotRange->All,PlotRangePadding->Scaled@0.05,Joined->False,Filling->Axis,FillingStyle->AbsoluteThickness@1.5,Axes->True};
purpleWhite[v_]:=Hue[0.75,0.25,Abs[v]^0.5];
blueWhiteRed[v_]:=Hue[If[v>0.5,0.0,0.5],(2Abs[v-0.5]),0.95];
colorFromComplex[z_,OptionsPattern[{Brightness->1.}]]:=Hue[Rescale[Arg@z,{0,2 \[Pi]}],0.8,Tanh@Abs[z*OptionValue[Brightness]]];
manualColor[zmin_,zmax_]:=ColorFunction->(ColorData["FuchsiaTones"][Rescale[#,{zmin,zmax}]]&);



(* ::Text:: *)
(*Set default options for plotting routines:*)


myIS=ImageSize->{384,256};
myISsq=ImageSize->{256,256};
myGS={
RotateLabel->False,Axes->False,Frame->True,GridLinesStyle->LightGray,  (*Plot elements*)
BaseStyle->{14,FontFamily->"Helvetica"},LabelStyle->{14,FontFamily->"Helvetica"},FrameStyle->Thick,myIS,(*LabelStyle->18*)
AspectRatio->Full
};
myPS={
PlotStyle->myStyles  (*Styling of curves -- not all graphics routines have this option! *)
};
SetOptions[Graphics,BaseStyle->{14,FontFamily->"Helvetica"},LabelStyle->{14,FontFamily->"Helvetica"}];
SetOptions[Plot,myIS,myGS,myPS,ExclusionsStyle->Automatic];
SetOptions[ParametricPlot,myIS,myGS,myPS,ExclusionsStyle->Automatic];
SetOptions[ListPlot,myIS,myGS,myPS,Joined->True];
SetOptions[ListPlot,myIS,myGS,myPS,Joined->True];
SetOptions[ListLogPlot,myIS,myGS,myPS,Joined->True];
SetOptions[ListLogLogPlot,myIS,myGS,myPS,Joined->True];
SetOptions[Histogram,myIS,myGS];
SetOptions[DensityPlot,myISsq,myGS,ColorFunctionScaling->False,AspectRatio->Full,PlotRange->{{0,1},{0,1},All}];
SetOptions[ListDensityPlot,myISsq,myGS,ColorFunctionScaling->False,InterpolationOrder->0];
SetOptions[ContourPlot,myISsq,myGS];
SetOptions[ListContourPlot,myISsq,myGS];
SetOptions[ArrayPlot,myISsq,myGS,ColorFunctionScaling->False,FrameTicks->{{All,None},{All,None}}];
SetOptions[PolarPlot,myPS,myGS];
SetOptions[ListContourPlot,PlotRange->Full,PlotRangePadding->0,PlotRangeClipping->False,ColorFunctionScaling->True];
SetOptions[Histogram,myGS,PlotRange->All,PlotRangePadding->{{0,0},{0,Scaled@0.05}}];
SetOptions[Histogram,Frame->True,PlotRange->All,PlotRangePadding->{{0,0},{0,Scaled@0.05}}];
SetOptions[ArrayPlot,ImageMargins->0,ColorFunctionScaling->False,FrameTicks->None,(*FrameTicks->Automatic,*)
AspectRatio->Full,PlotRangePadding->Scaled[0.]];
SetOptions[DensityPlot,ImageMargins->0,ColorFunctionScaling->False,FrameTicks->{{None,All},{None,None}},PlotPoints->{2,8}];
Clear[myIS,myISsq,myGS,myPS];



(* ::Subsection::GrayLevel[0]:: *)
(*1.5	Graphics utilities: enlargeFont,  makeLineLegend, arrayPlot, colorBar, myBarChart, mySmoothedBarChart*)


(*enlargeFont=Module[{s},(#/.s_String->Style[s,FontFamily->"Times",20]&)];*)
enlargeFont=#1/.s_String->Style[s,24]&;
makeLineLegend[styles_,labels_,OptionsPattern[{LegendPosition->{0.75,0.75},LegendImageSize->64}]]:=
Module[{trimmedStyles,trimmedLabels,numKeys},
(*----- Idiot-proofing -----*)
numKeys=Min[Length@styles,Length@labels];
trimmedStyles=styles[[1;;numKeys]];
trimmedLabels=labels[[1;;numKeys]];
(*----- Generate the graphics -----*)
Graphics@Inset[#,Scaled@OptionValue@LegendPosition]&@Grid[MapThread[{Graphics[{Directive@#1,Line@{{-1,0},{1,0}}},ImageSize->OptionValue@LegendImageSize,AspectRatio->0.01],#2}&,{trimmedStyles,trimmedLabels}],Alignment->Left]
];
makeLineLegend2[styles_,labels_]:=
Module[{trimmedStyles,trimmedLabels,numKeys},
(*----- Idiot-proofing -----*)
numKeys=Min[Length@styles,Length@labels];
trimmedStyles=styles[[1;;numKeys]];
trimmedLabels=labels[[1;;numKeys]];
(*----- Generate the graphics -----*)
Graphics@Text@Grid[
MapThread[{Graphics[{Directive@#1,Line@{{-1,0},{1,0}}},AspectRatio->0.01,ImageSize->24],#2}&,{trimmedStyles,trimmedLabels}],Alignment->Left,Spacings->0]
];
makeLineLegend2::usage="makeLineLegend2[{style1, style2, ...},{label1, label2, ...}]
 returns a legend consisting of styled lines and labels.  This legend is a Graphics object.";



Options[arrayPlot]=Join[{},Options@Plot];
arrayPlot$ImagePadding={{36,8},{44,8}};
arrayPlot[gxy_,opts:OptionsPattern[]]/;ArrayDepth@gxy==2:=
Module[{l,r,b,t,ll=0,rr=64,w=192,h,xmax,ymax,cf},
(*---- Extract options ----*)
{{l,r},{b,t}}=arrayPlot$ImagePadding;
{gmin,gmax}=arrayPlot$ZRange;
{xmax,ymax}=Dimensions[gxy];
neato`DensityPlotHeight=h=Round[w*ymax/xmax];(* choose image height acc/to aspect ratio of data, and save in GLOBAL variable! *)
ArrayPlot[gxy//Transpose//Reverse,
ColorFunction->(arrayPlot$ColorFunction@Rescale[#,arrayPlot$ZRange,{0,1}]&),ColorFunctionScaling->False,
FilterRules[{opts},Options@ArrayPlot],
FrameLabel->Reverse@{"\!\(\*
StyleBox[\"x\",\nFontSlant\[Rule]\"Italic\"]\)","\!\(\*
StyleBox[\"y\",\nFontSlant\[Rule]\"Italic\"]\)"},RotateLabel->False,
(*FrameTicks->Reverse@{{1,Floor[xmax/2],xmax},{1,Floor[ymax/2],ymax}},*)
DataRange->{{1,xmax},{1,ymax}},
ImageSize->{w+l+r,h+b+t},ImagePadding->{{l,r},{b,t}}
]
];
Options[colorBar]=Join[{},Options@Plot];
colorBar[opts:OptionsPattern[]]:=
Module[{l,r,b,t,ll=0,rr=64,ww=24,w=192,h,xmax,ymax,ngmax=50},
{{l,r},{b,t}}=arrayPlot$ImagePadding;
h=neato`DensityPlotHeight;(* read global variable  *)
ArrayPlot[Table[Rescale[ng,{1,ngmax},arrayPlot$ZRange],{ng,1,ngmax},{ndummy,2}]//Reverse,
ColorFunction->(arrayPlot$ColorFunction@Rescale[#,arrayPlot$ZRange,{0,1}]&),ColorFunctionScaling->False,
FilterRules[{opts},Options@ArrayPlot],
DataRange->{All,arrayPlot$ZRange},
FrameTicks->{None,None,Union[arrayPlot$ZRange,FindDivisions[arrayPlot$ZRange,4]]//N,None},
ImageSize->{ww+ll+rr,h+b+t},ImagePadding->{{ll,rr},{b,t}}]
];


arrayPlot[gxyz_,opts:OptionsPattern[]]/;ArrayDepth@gxyz==3:=
Module[{l,r,b,t,ll=0,rr=64,w=192,h,xmax,ymax,zmax,cf},
(*---- Extract options ----*)
{{l,r},{b,t}}=arrayPlot$ImagePadding;
{gmin,gmax}=arrayPlot$ZRange;
{xmax,ymax,zmax}=Dimensions[gxyz];
neato`DensityPlotHeight=h=Round[w*ymax/xmax];(* choose image height acc/to aspect ratio of data, and save in GLOBAL variable! *)

GraphicsRow[Table[
ArrayPlot[gxyz[[All,All,z]]//Transpose//Reverse,
ColorFunction->(arrayPlot$ColorFunction@Rescale[#,arrayPlot$ZRange,{0,1}]&),ColorFunctionScaling->False,
FilterRules[{opts},Options@ArrayPlot],
DataRange->{{1,xmax},{1,ymax}},
ImageSize->{w+0+1,h+0+1},ImagePadding->{{0,1},{1,1}}
],{z,zmax}],Spacings->0]
];


Options[myBarChart]={}~Join~Options@Plot;
myBarChart[\[Xi]Categories_,A\[Xi]_,opts:OptionsPattern[]]:=Module[{points,n\[Xi]max},
n\[Xi]max=Length@A\[Xi];
points=Join[
{{\[Xi]Categories[[1]],0}},
Join@@Table[{{\[Xi]Categories[[n\[Xi]]],A\[Xi][[n\[Xi]]]/(\[Xi]Categories[[n\[Xi]+1]]-\[Xi]Categories[[n\[Xi]]])},{\[Xi]Categories[[n\[Xi]+1]],A\[Xi][[n\[Xi]]]/(\[Xi]Categories[[n\[Xi]+1]]-\[Xi]Categories[[n\[Xi]]])}},{n\[Xi],n\[Xi]max}],
{{\[Xi]Categories[[n\[Xi]max+1]],0}}
];
ListPlot[points,
FilterRules[{opts}, Options@Plot],
PlotRange->All,Filling->Axis]
];
(*-------- To smear a PDF (prob. dist. func.) without spoiling sum rule, we should smear the CDF. --------*)
Options[mySmoothedBarChart]={}~Join~Options@Plot;
mySmoothedBarChart[\[Xi]Categories_,A\[Xi]_,opts:OptionsPattern[]]:=Module[{points,n\[Xi]max},
cdf=FoldList[Plus,0,A\[Xi]];  (*cumulative distribution*)
cdf\[LongDash]interpolated=Interpolation[Transpose@{\[Xi]Categories,cdf},InterpolationOrder->5];  (*find cubic interpolant to CDF*)
pdf\[LongDash]smoothed[x_]=\!\(
\*SubscriptBox[\(\[PartialD]\), \(x\)]\(cdf\[LongDash]interpolated[x]\)\);
Plot[pdf\[LongDash]smoothed[x],{x,Min[\[Xi]Categories],Max[\[Xi]Categories]},
FilterRules[{opts}, Options@Plot],PlotRange->All]
];


stemStyle={PlotMarkers->Automatic,PlotRange->All,PlotRangePadding->Scaled@0.05,Joined->False,Filling->Axis,FillingStyle->AbsoluteThickness@1.5,Axes->True};


colorFromComplex[z_,OptionsPattern[{Brightness->1.}]]:=Hue[Rescale[Arg@z,{0,2 \[Pi]}],0.8,Tanh@Abs[z*OptionValue[Brightness]]]


myDensityPlot[func_,{x_,xmin_,xmax_,xpoints_},{y_,ymin_,ymax_,ypoints_},opts:OptionsPattern[]]:=
arrayPlot[
Table[func[x,y],
{x,Rescale[Range[xpoints],{1,xpoints},{xmin,xmax}]},
{y,Rescale[Range[ypoints],{1,ypoints},{ymin,ymax}]}
],
ImageSize->{256+64+8,256+64+32},ImagePadding->{{64,8},{64,32}},opts,
PlotRange->All,PlotRangePadding->None,DataRange->{{xmin,xmax},{ymin,ymax}},
(*FrameTicks->{All,All,None,None}*)
FrameTicks->{Automatic,Automatic,None,None}
];


(* ::Section::GrayLevel[0]:: *)
(*2	Lattice utilities*)


(* ::Subsection:: *)
(*2.1	packIndex, unpackIndex, setupHypercubicLattice*)


(*\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash] INDEX PACKING/UNPACKING \[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]*)
packIndex::usage=unpackIndex::usage="
These index packing and unpacking routines are suitable for any dimension.
The following examples are the 3-dimensional versions.
Note that the indices are packed in Fortran order, as opposed to C order.

packIndex[{x,y,z}, L, 3] returns x+L*(y+L*z).
packIndex[{x,y,z}, {xmax,ymax,zmax}] returns x+xmax*(y+ymax*z).

unpackIndex[i, L, 3] performs the inverse transformation and returns {x,y,z}.
unpackIndex[i, {xmax,ymax,zmax}] returns {x,y,z}.
";
packIndex[c_,cmaxd_Integer]:=FromDigits[Reverse@c,cmaxd];
unpackIndex[c_,cmaxd_Integer,dmax_Integer]:=Reverse@IntegerDigits[c,cmaxd,dmax];
packIndex[cd_,cmaxd_List]:=Module[{dmax=Length@cmaxd,c=cd[[-1]]},
Do[c=cd[[d]]+cmaxd[[d]]*c,{d,dmax-1,1,-1}];
Return[c];
];
unpackIndex[c_,cmaxd_List]:=Module[{q=c,r,dmax=Length@cmaxd},
Table[{q,r}=QuotientRemainder[q,cmaxd[[d]]];r,{d,dmax}]
];

(*----  c\[LongDash]from\[LongDash]cd is included for backwards compatibility  ----*)
c\[LongDash]from\[LongDash]cd[cd_,cmaxd_]:=Module[{dmax=Length@cmaxd,c=cd[[-1]]},
Do[c=cd[[d]]+cmaxd[[d]]*c,{d,dmax-1,1,-1}];
Return[c];
];
cd\[LongDash]from\[LongDash]c[c_,cmaxd_]:=Module[{q=c,r,dmax=Length@cmaxd},
Table[{q,r}=QuotientRemainder[q,cmaxd[[d]]];r,{d,dmax}]
];
(*cd\[LongDash]from\[LongDash]c[c_,cmaxd_]:=Last/@Rest@FoldList[QuotientRemainder[#1[[1]],#2]&,{c,0},cmaxd];*)

(*\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash] HYPERCUBIC LATTICE \[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]*)
setupHypercubicLattice[]:=Module[{a,b,c,d,cd,ctargetd,A,B},
kmax=Times@@kmaxd;  (*Do this as a courtesy*)

dmax=Length@cmaxd;

amax=1;
bmax=dmax;
cbd=IdentityMatrix[dmax];
ab0=ab1=Table[1,{bmax}];
tb=Table[1.,{bmax}];
cad=Table[0.,{amax},{dmax}];

cmax=Times@@cmaxd;
AMAX=amax*cmax;
BMAX=bmax*cmax;

BMAXAMAX=2*bmax;
AAB=Table[0,{AMAX},{BMAXAMAX}];
ABA=Table[0,{BMAX},{2}];
BAB=Table[0,{AMAX},{BMAXAMAX}];
CAD=Table[0,{AMAX},{dmax}];
CBD=Table[0,{BMAX},{dmax}];
RDD=DiagonalMatrix[cmaxd]//N;

Do[
Do[
cd=cd\[LongDash]from\[LongDash]c[c,cmaxd];
ctargetd=Table[Mod[cd[[d]]+KroneckerDelta[b,d],cmaxd[[d]]],{d,dmax}];
ctarget=c\[LongDash]from\[LongDash]cd[ctargetd,cmaxd];

B=b+bmax*c;
CBD[[B]]=Table[Quotient[cd[[d]]+KroneckerDelta[b,d],cmaxd[[d]]],{d,dmax}];

A=1+c;
ATARGET=1+ctarget;
AAB[[A,1+2(b-1)]]=ATARGET;
AAB[[ATARGET,2+2(b-1)]]=A;

ABA[[B,1]]=A;
ABA[[B,2]]=ATARGET;

BAB[[A,1+2(b-1)]]=B;
BAB[[ATARGET,2+2(b-1)]]=BitNot@B;

CAD[[A]]=cd/cmaxd//N;

,{b,bmax}];
,{c,0,cmax-1}];

getCBD[B_]:=If[B>=0,CBD[[B]],-CBD[[BitNot@B]]];
(*getABA[B_,A_]:=If[B>=0,ABA[[B,A]],ABA[[BitNot@B,3-A]]];*)
startAtom[B_]:=If[B>=0,ABA[[B,1]],ABA[[BitNot@B,2]]];
endAtom[B_]:=If[B>=0,ABA[[B,2]],ABA[[BitNot@B,1]]];
neighbor[A_,b_]:=AAB[[A,b]];
outgoingBond[A_,b_]:=BAB[[A,b]];
neighbors[A_]:=BMAXAMAX;
]




(* ::Subsubsection:: *)
(*Deprecated/unimplemented routines*)


(*\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash] CUBIC LATTICE (UGLY OLD CODE) \[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]\[LongDash]*)
(*-------- Reads cmaxd AND kmaxd --------*)
setupCubicLattice[]:=Module[{a,b,c,d,cd,ctargetd,A,B},
dmax=3;    (* chain, square or cubic -- all will be embedded in 3D space *)
amax=1;
Switch[Count[cmaxd*kmaxd,x_/;x>1],   (* number of dimensions that are greater than 1... uGLY HACK *)
1,    bmax=1;cbd=({
 {1, 0, 0}
});ab0=ab1={1};tb={1.};,
2,    bmax=2;cbd=({
 {1, 0, 0},
 {0, 1, 0}
});ab0=ab1={1,1};tb={1.,1.};,
3,    bmax=3;cbd=({
 {1, 0, 0},
 {0, 1, 0},
 {0, 0, 1}
});ab0=ab1={1,1,1};tb={1.,1.,1.};
];
cad={0.,0.,0.};
{cmax1,cmax2,cmax3}=cmaxd;cmax=Times@@cmaxd;
AMAX=amax*cmax;
BMAX=bmax*cmax;

BMAXAMAX=2*bmax;
AAB=Table[0,{AMAX},{BMAXAMAX}];
ABA=Table[0,{BMAX},{2}];
BAB=Table[0,{AMAX},{BMAXAMAX}];
CAD=Table[0,{AMAX},{dmax}];
CBD=Table[0,{BMAX},{dmax}];
RDD=DiagonalMatrix[cmaxd]//N;


Do[
Do[
cd=cd\[LongDash]from\[LongDash]c[c,cmaxd];
ctargetd=Table[Mod[cd[[d]]+KroneckerDelta[b,d],cmaxd[[d]]],{d,dmax}];
ctarget=c\[LongDash]from\[LongDash]cd[ctargetd,cmaxd];

B=b+bmax*c;
CBD[[B]]=Table[Quotient[cd[[d]]+KroneckerDelta[b,d],cmaxd[[d]]],{d,dmax}];

A=1+c;
ATARGET=1+ctarget;
AAB[[A,1+2(b-1)]]=ATARGET;
AAB[[ATARGET,2+2(b-1)]]=A;

ABA[[B,1]]=A;
ABA[[B,2]]=ATARGET;

BAB[[A,1+2(b-1)]]=B;
BAB[[ATARGET,2+2(b-1)]]=BitNot@B;

CAD[[A]]=cd/cmaxd//N;

,{b,bmax}];
,{c,0,cmax-1}];

getCBD[B_]:=If[B>=0,CBD[[B]],-CBD[[BitNot@B]]];
(*getABA[B_,A_]:=If[B>=0,ABA[[B,A]],ABA[[BitNot@B,3-A]]];*)
startAtom[B_]:=If[B>=0,ABA[[B,1]],ABA[[BitNot@B,2]]];
endAtom[B_]:=If[B>=0,ABA[[B,2]],ABA[[BitNot@B,1]]];
neighbor[A_,b_]:=AAB[[A,b]];
outgoingBond[A_,b_]:=BAB[[A,b]];
neighbors[A_]:=BMAXAMAX;

];
setupBravaisLattice[]:=Print["Not implemented!"];
setupGeneralLattice[]:=Print["Not implemented!"];


(* ::Section::GrayLevel[0]:: *)
(*3	Obsolete versions*)


(*toString[x_]:=ExportString[x//List//Flatten//List,"Table","FieldSeparators"->" "];*)
(*substVars[s_String]:=StringReplace[s,RegularExpression@"\\$(\\w+)":>toString@ToExpression@"$1"];*)

