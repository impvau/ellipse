
Remove[DataTrain];
Remove[DataAll];

BeginPackage["EllipseData`"];

(* Export all symbols that will be defined in subpackages *)
DataTrain::usage = "DataTrain[] generates data used to train on.";
DataAll::usage = "DataAll[minZ, maxZ, n] generates data used to plot and evaluate over a larger range, with uniform z values from minZ to maxZ and corresponding a values.";

Begin["`Private`"];

(* Function to generate and return sorted combinations *)
DataTrain[] := Module[
  {additionalCombinations, aRangeCombinations, combinations},
  
  (* Define additional specific combinations *)
  additionalCombinations = {
    (*{1,1},*){105/100,1},{115/100,1},{125/100,1},{135/100,1},{145/100,1},{155/100,1},{165/100,1},{175/100,1},{185/100,1},{195/100,1},{40,1},{50,1},{60,1},{70,1},{80,1},{90,1},{100,1},{500,1},{1000,1},{10000,1},{1000000,1},{1000000000,1}
  };
  
  (* Generate a range of combinations for a = 2 to 30 *)
  aRangeCombinations = Table[{a, 1}, {a, 2, 30}];
  
  (* Join and sort combinations by increasing a *)
  combinations = Join[additionalCombinations, aRangeCombinations];
  SortBy[combinations, First]
];

(* Function to generate combinations based on a transformation formula *)
DataAll[minZ_, maxZ_, n_Integer] := Module[
  {zValues, aValues, combinations},
  zValues = Subdivide[minZ, maxZ, n - 1]; (* Uniform z values from minZ to maxZ *)
  aValues = #/(1 - #) & /@ zValues;       (* Calculate corresponding a values using the transformation formula *)
  combinations = Table[{aValues[[i]], zValues[[i]]}, {i, 1, n}]; (* Pair each a value with its corresponding z value *)
  combinations
];

End[];
EndPackage[];
