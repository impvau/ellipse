Remove[createAndSavePlotSR]
BeginPackage["Plots`"];

createAndSavePlotSR::usage = 
  "createAndSavePlotSR[filename, selectedKeys, results, x, scalingX, scalingY, plotTitle] creates a plot based on keys excluding 'x', applies specified scaling, and exports the plot with 'plotTitle' to 'filename'.";

Begin["`Private`"];

createAndSavePlotSR[filename_String, selectedKeys_List, results_, x_, scalingX_: "Linear", scalingY_: "Linear", plotTitle_: ""] := Module[
  {
    errorTypes, colors, plotData, legendLabels, plot, validData
  },
  
  errorTypes = DeleteCases[selectedKeys, x];

  colors = Table[If[i == 1, Black, 
    ColorData["DarkRainbow"][(i - 2)/(Length[errorTypes] - 2)]], {i, Length[errorTypes]}];

  plotData = Table[Tooltip[
    Table[{results[[idx, x]], results[[idx, errorType]]}, {idx, Length[results]}], 
    Style[Subscript[StringTake[errorType, 1], 
      StringTake[errorType, {2, First[Flatten[{StringPosition[errorType, " "], Length[errorType] + 1}]] - 1}]], 16]
  ], {errorType, errorTypes}];

  legendLabels = Style[Subscript[ToLowerCase[StringTake[#, 1]], 
    StringTake[#, {2, First[Flatten[{StringPosition[#, " "], Length[#] + 1}]] - 1}]], FontSize -> 30] & /@ errorTypes;

  plot = ListPlot[
    plotData,
    ScalingFunctions -> {scalingX, scalingY},
    Axes -> False,
    Frame -> True,
    FrameLabel -> {{Style[ToString[Subscript["Log","10"],StandardForm] <> "(Relative Absolute Error)", 30, FontFamily -> "Helvetica"], None}, {Style[ToString[Subscript["Log","10"],StandardForm] <> "(" <> x <> ")", 30, FontFamily -> "Helvetica"], None}},
    (*FrameLabel -> {{Style["Relative Absolute Error", 30, FontFamily -> "Helvetica"], None}, {Style[x, 30, FontFamily -> "Helvetica"], None}},*)
    (*PlotLabel -> If[plotTitle != "", Style[plotTitle, 30, FontFamily -> "Helvetica"], None],  (* New line for title *)*)
    PlotLabel -> If[plotTitle =!= Null, Style[plotTitle, 30, FontFamily -> "Helvetica"], None],
    LabelStyle -> {FontSize -> 30, FontFamily -> "Helvetica"},
    Joined -> True,
    (*PlotRange -> {Automatic, {0.750005,0.749995}},*)
    (*PlotRange -> {{0,80}, All},*)
    PlotRange -> {Automatic, Automatic},
    Ticks -> {Automatic, Automatic},
    ImageSize -> {969, 603},
    GridLines -> None,
    PlotStyle -> (Directive[Thick, #] & /@ colors),
    PlotLegends -> Placed[LineLegend[
      (Directive[Thick, AbsoluteThickness[4], #] & /@ colors),
      (Style[#, 30] & /@ legendLabels), 
      LegendLayout -> "Row"], 
      (*{Right, Top},*)
      {Scaled[{0.95, 0.05}], {Right, Bottom}},
      LegendFunction -> (Framed[#, FrameMargins -> 10, Background -> White] &)]
  ];

  Export[filename, plot];
  
  plot

];

End[]; 

EndPackage[];
