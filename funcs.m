Remove[GetEllipseFuncs];

BeginPackage["EllipseFuncs`"];

GetEllipseFuncs::usage = 
  "GetEllipseFuncs[data_List] calculates various ellipse approximations \
expressions based on a list of inputs {a, b}.";

Begin["`Private`"];

(* Set precision for all calculations *)
prec = 200;

(* Redefine GetEllipseFuncs to accept a list of {a, b} pairs *)
GetEllipseFuncs[data_List] := Module[
   {results, pi = SetPrecision[Pi, prec]},
   results = Table[
     Module[
       {a = SetPrecision[pair[[1]], prec], b = SetPrecision[pair[[2]], prec], f = <||>},

       (* Direct calculations with consistent precision *)
       f["a"] = a;
       f["z"] = SetPrecision[a/(a+1),prec];
       f["b"] = b;
       f["y"] = SetPrecision[4*EllipticE[1 - b^2/a^2]*a, prec];
       f["ynp"] = SetPrecision[f["y"]/(pi*(a+b)), prec];
       f["h"] = SetPrecision[((a - b)/(a + b))^2, prec];
       f["k"] = SetPrecision[Sqrt[1 - Min[a, b]^2/Max[a, b]^2], prec];

		

       (* Calculations depending on previously computed values, all with explicit precision setting *)
       f["MR"] = SetPrecision[pi*(a + b)*(1 + (3*f["h"])/(10 + Sqrt[4 - 3*f["h"]])), prec];
       f["Mu"] = SetPrecision[pi*(a + b)*(1 + (44/pi - 11)*f["h"]/(10 + Sqrt[4 - 3*f["h"]])), prec];
       f["h1"] = SetPrecision[(((f["MR"] - f["Mu"])^2)/((f["MR"] + f["Mu"])^2)), prec];
       f["MRu"] = SetPrecision[f["MR"] + 9938*f["Mu"]*f["h"]^7*f["h1"], prec];
       f["Mp"] = SetPrecision[pi*(a + b)*(135168 - 85760*f["h"] - 5568*f["h"]^2 + 3867*f["h"]^3)/(135168 - 119552*f["h"] + 22208*f["h"]^2 - 345*f["h"]^3), prec];
       f["h2"] = SetPrecision[(((f["Mp"] - f["Mu"])^2)/((f["Mp"] + f["Mu"])^2)), prec];
       f["M1"] = SetPrecision[f["Mu"] - ((((-125494663)/(a + 59)) + 2392263)*f["h1"]), prec];
       f["M2"] = SetPrecision[f["Mp"]*(9686.95*(f["h1"] - f["h2"])) + f["Mp"], prec];
       f["M3"] = SetPrecision[f["Mp"]*(-581184*f["h1"]/(39 + a)) + f["Mu"], prec];
       f["M4"] = SetPrecision[f["Mp"]*Power[f["Mu"]/f["Mp"],Power[(f["h1"]/615)*f["h1"],f["h1"]/(a*f["h2"])]], prec];
       f["M5"] = SetPrecision[f["Mp"]*Power[f["Mu"]/f["Mp"],Power[22,((-13)/(f["h2"]*a))*f["h1"]]], prec];
       (* Old one with less effective coefficient
       f["M6"] = SetPrecision[f["Mp"]*Power[f["Mu"]/f["Mp"],Power[f["h1"]*f["h1"]/435,(f["h1"]/f["h2"])/a]], prec];
       *)
       f["M6"] = SetPrecision[f["Mp"]*Power[f["Mu"]/f["Mp"],Power[f["h1"]*f["h1"]/615,(f["h1"]/f["h2"])/a]], prec];
       
       f["hd"] = SetPrecision[ 1/f["h"], prec];
		f["h1d"] = SetPrecision[ 1/f["h1"], prec];
		f["h2d"] = SetPrecision[ 1/f["h2"], prec];
		
       
       (* Equations from Sykora *)
       (* Keplerian *)
       f["K1"] = SetPrecision[2*pi*Sqrt[a*b],prec];
       f["K2"] = SetPrecision[2*pi*Power[(a+b),2]/Power[(Sqrt[a]+Sqrt[b]),2],prec];
       f["K3"] = SetPrecision[pi*(a+b),prec];
       f["K4"] = SetPrecision[pi*(a+b)*((3-Sqrt[1-f["h"]])/2),prec];
       f["K5"] = SetPrecision[pi*Sqrt[2*(Power[a,2]+Power[b,2])],prec];
       f["K6"] = SetPrecision[2*Pi*((2*(a + b)^2 - (Sqrt[a] - Sqrt[b])^4) / ((Sqrt[a] + Sqrt[b])^2 + 2*Sqrt[2*(a + b)]*Power[a*b,1/4])), prec];
       f["K7"] = SetPrecision[(pi/2)*Sqrt[6*(Power[a,2]+Power[b,2])+4*a*b],prec];
       f["K8"] = SetPrecision[(2*pi)*Power[(Power[a,3/2]+Power[b,3/2])/2,2/3],prec];
       f["K9"] = SetPrecision[pi*(a+b)*Power[1+f["h"]/8,2],prec];
       f["K10"] = SetPrecision[pi*(3*(a+b)-Sqrt[(a+3*b)*(3*a+b)]),prec];
       f["K11"] = SetPrecision[(pi/4)*(6+(1/2)*(Power[a-b,2])/(Power[a+b,2])),prec];
       f["K12"] = f["MR"];
       (* Pade *)
       f["Ka"] = SetPrecision[ pi*(a+b)*((16+3*f["h"])/(16-f["h"])), prec];
       f["Kb"] = SetPrecision[ pi*(a+b)*((64+16*f["h"])/(64-Power[f["h"],2])), prec];
       f["Kc"] = SetPrecision[ pi*(a+b)*((64-3*Power[f["h"],2])/(64-16*f["h"])), prec];
       f["Kd"] = SetPrecision[ pi*(a+b)*((256-48*f["h"]-21*Power[f["h"],2])/(256-112*f["h"]+3*Power[f["h"],2])), prec];
       f["Ke"] = SetPrecision[ pi*(a+b)*((3072-1280*f["h"]-252*Power[f["h"],2]+33*Power[f["h"],3])/(3072-2048*f["h"]+212*Power[f["h"],2])), prec];
       (*Optimized Peano*)
       f["O2"] = SetPrecision[ pi*Sqrt[2*(Power[a,2]+Power[b,2])-Power[a-b,2]/2.458338], prec];
       (*Extremes*)
       f["E1"] = SetPrecision[pi*(a-b)/ArcTan[(a-b)/(a+b)],prec];
       f["E2"] = SetPrecision[4*(a+b)-((8-2*pi)*a*b)/(0.410117*(a+b)+(1-2*0.410117)*(Sqrt[(a+74*b)*(74*a+b)]/(1+74))),prec];
       f["E3"] = SetPrecision[2*Sqrt[Power[pi,2]*a*b+4*Power[a-b,2]],prec];
       f["E4"] = SetPrecision[4*((Power[b,2]/a)*ArcTan[a/b]+(Power[a,2]/b)*ArcTan[b/a]),prec];
       f["E5"] = SetPrecision[(4*pi*a*b+Power[a-b,2])/(a+b), prec];
       f["s"] = SetPrecision[Log[2]/Log[pi/2],prec];
       f["E6"] = SetPrecision[4*Power[Power[a,f["s"]]+Power[b,f["s"]],1/f["s"]], prec];
       (*Extreme Pades*)
       f["Ead1"] = SetPrecision[(pi/4)*(81/64)-1,prec];
       f["Ead2"] = SetPrecision[(pi/4)*(19/15)-1,prec];
       f["Eap"]=f["Ead1"]/(f["Ead1"]-f["Ead2"]);
       f["Ea"] = SetPrecision[pi*(a+b)*(f["Eap"]*(16+3*f["h"])/(16-f["h"])+(1-f["Eap"])*Power[1+f["h"]/8,2]), prec];
       
       f["Ebd1"] = SetPrecision[(pi/4)*(80/63)-1,prec];
       f["Ebd2"] = SetPrecision[(pi/4)*(61/48)-1,prec];
       f["Ebp"]=f["Ebd1"]/(f["Ebd1"]-f["Ebd2"]);
       f["Eb"] = SetPrecision[pi*(a+b)*(f["Ebp"]*(64-3*Power[f["h"],2])/(64-16*f["h"])+(1-f["Ebp"])*((64+16*f["h"])/(64-f["h"]^2))), prec];
       
       f["Ecd1"] = SetPrecision[(pi/4)*(61/48)-1,prec];
       f["Ecd2"] = SetPrecision[(pi/4)*(187/147)-1,prec];
       f["Ecp"]=f["Ecd1"]/(f["Ecd1"]-f["Ecd2"]);
       f["Ec"] = SetPrecision[pi*(a+b)*(f["Ecp"]*(256-48*f["h"]-21*Power[f["h"],2])/(256-112*f["h"]+3*Power[f["h"],2])+(1-f["Ecp"])*(64-3*Power[f["h"],2])/(64-16*f["h"])), prec];
       
       f["Edd1"] = SetPrecision[(pi/4)*(187/147)-1,prec];
       f["Edd2"] = SetPrecision[(pi/4)*(1573/1236)-1,prec];
       f["Edp"]=f["Edd1"]/(f["Edd1"]-f["Edd2"]);
       f["Ed"] = SetPrecision[pi*(a+b)*(f["Edp"]*(3072-1280*f["h"]-252*Power[f["h"],2]+33*Power[f["h"],3])/(3072-2048*f["h"]+212*Power[f["h"],2])+(1-f["Edp"])*(256-48*f["h"]-21*Power[f["h"],2])/(256-112*f["h"]+3*Power[f["h"],2])), prec];
       
       f["Eed1"] = SetPrecision[(pi/4)*(1573/1236)-1,prec];
       f["Eed2"] = SetPrecision[(pi/4)*(47707/37479)-1,prec];
       f["Eep"]=f["Eed1"]/(f["Eed1"]-f["Eed2"]);
       f["Ee"] = SetPrecision[pi*(a + b)*((f["Eep"])*(135168 - 85760*f["h"] - 5568*f["h"]^2 + 3867*f["h"]^3)/(135168 - 119552*f["h"] + 22208*f["h"]^2 - 345*f["h"]^3)+(1-f["Eep"])((3072-1280*f["h"]-252*Power[f["h"],2]+33*Power[f["h"],3])/(3072-2048*f["h"]+Power[212*f["h"],2]))), prec];
       
       f["C1t"] = SetPrecision[(pi/4)*((a-b)/b),prec];
       f["C1"] = SetPrecision[pi*Sqrt[2*(Power[a,2]+Power[b,2])]*(Sin[f["C1t"]]/f["C1t"]),prec];
       
       f["C2"] = SetPrecision[4*a+2*(pi-2)*a*Power[b/a,1.456],prec];
       f["C3"] = SetPrecision[4*((pi*a*b+Power[a-b,2])/(a+b))-(89/146)*Power[(b*Sqrt[a]-a*Sqrt[b])/(a+b),2],prec];
       
       f["C4s"]=0.825056176207;
       f["C4"] = SetPrecision[4*(a+b)-(2*(4-pi)*a*b)/Power[(Power[a,f["C4s"]]+Power[b,f["C4s"]])/2,1/f["C4s"]],prec];
       
       f["C5"] = SetPrecision[4*((pi*a*b+Power[a-b,2])/(a+b))-(1/2)*(a*b/(a+b))*(Power[a-b,2]/(pi*a*b+Power[a+b,2])),prec];
       f["C6"] = SetPrecision[pi*(a + b)*(1 + (3*f["h"])/(10 + Sqrt[4 - 3*f["h"]])+(4/pi-14/11)*Power[f["h"],12]),prec];
       
       f["C7p"] = 3.982901;
       f["C7q"] = 66.71674;
       f["C7r"] = 56.2007;
       f["C7s"] = 18.31287;
       f["C7t"] = 23.39728;
       f["C7"] = SetPrecision[4*(a+b)-((a*b)/(a+b))*((f["C7p"]*Power[a+b,2]+f["C7q"]*a*b+f["C7r"]*Power[a*b/(a+b),2])/(Power[a+b,2]+f["C7s"]*a*b+f["C7t"]*Power[a*b/(a+b),2])),prec];
       
       f["A1"] = SetPrecision[4*(Power[a,2]+Power[b,2]+(pi-2)*a*b)/(a+b),prec];

	   f["A2t1"] = SetPrecision[2.49808365277126, prec];
	   f["A2s1"] = SetPrecision[1.22694921875000, prec];
       f["A2"] = SetPrecision[(4*(a^3 + b^3 + f["A2t1"]*a*b*(a+b))/(a^2 + b^2 + 2*f["A2s1"]*a*b)), prec];

		f["A3t1"] = SetPrecision[6.16881239339582, prec];
		f["A3t2"] = SetPrecision[4.06617730084445, prec];
		f["A3s1"] = SetPrecision[6.15241658169936, prec];
		f["A3"] = SetPrecision[(4*(a^4 + b^4 + f["A3t1"]*a*b*(a^2 + b^2)+2*f["A3t2"]*(a^2*b^2))/(a^3 + b^3 + f["A3s1"]*a*b*(a + b))), prec];

		f["A4t1"] = SetPrecision[13.02487942169925, prec];
		f["A4t2"] = SetPrecision[28.56997512074272, prec];
		f["A4s1"] = SetPrecision[13.01750519704827, prec];
		f["A4s2"] = SetPrecision[13.09922140579137, prec];
		f["A4"] = SetPrecision[(4*(a^5 + b^5 + f["A4t1"]*a*b*(a^3 + b^3) + f["A4t2"]*a^2*b^2*(a+b))/(a^4 + b^4 + f["A4s1"]*a*b*(a^2 + b^2)+2* f["A4s2"]*a^2*b^2)), prec];
	
		f["A5t1"] = SetPrecision[27.301243680755, prec];
		f["A5t2"] = SetPrecision[113.302483206429, prec];
		f["A5t3"] = SetPrecision[76.50091476282086, prec];
		f["A5s1"] = SetPrecision[27.297854333670, prec];
		f["A5s2"] = SetPrecision[110.551872985869, prec];
		f["A5"] = SetPrecision[(4*(a^6 + b^6 + f["A5t1"]*a*b*(a^4 + b^4) + f["A5t2"]*a^2*b^2(a^2 + b^2) + 2*f["A5t3"]*a^3*b^3))/(a^5 + b^5 + f["A5s1"]*a*b*(a^3 + b^3) + f["A5s2"]*a^2*b^2(a + b)), prec];

		f["A6t1"] = SetPrecision[51.447782789130, prec];
		f["A6t2"] = SetPrecision[385.327854851892, prec];
		f["A6t3"] = SetPrecision[790.4535392309255, prec];
		f["A6s1"] = SetPrecision[51.445996674310, prec];
		f["A6s2"] = SetPrecision[382.256974433855, prec];
		f["A6s3"] = SetPrecision[347.212007887276, prec];
		f["A6"] = SetPrecision[(4*(a^7 + b^7 + f["A6t1"]*a*b*(a^5 + b^5) + f["A6t2"]*a^2*b^2*(a^3 + b^3) + f["A6t3"]*a^3*b^3*(a + b)))/
                        (a^6 + b^6 + f["A6s1"]*a*b*(a^4 + b^4) + f["A6s2"]*a^2*b^2*(a^2 + b^2) + 2*f["A6s3"]*a^3*b^3), prec];
                        
         f["A7t1"] = SetPrecision[93.49235523473, prec];
		f["A7t2"] = SetPrecision[1262.73239571330, prec];
		f["A7t3"] = SetPrecision[4296.45229646421, prec];
		f["A7t4"] = SetPrecision[2903.735611540449, prec];
		f["A7s1"] = SetPrecision[93.49135794687, prec];
		f["A7s2"] = SetPrecision[1259.36473022183, prec];
		f["A7s3"] = SetPrecision[4093.96201082922, prec];
		f["A7"] = SetPrecision[(4*(a^8 + b^8 + f["A7t1"]*a*b*(a^6 + b^6) + f["A7t2"]*a^2*b^2*(a^4 + b^4) + f["A7t3"]*a^3*b^3*(a^2 + b^2) + 2*f["A7t4"]*a^4*b^4))/
                        (a^7 + b^7 + f["A7s1"]*a*b*(a^5 + b^5) + f["A7s2"]*a^2*b^2*(a^3 + b^3) + f["A7s3"]*a^3*b^3*(a+b)), prec];
            
        f["S0q"] = SetPrecision[4.16102118885517, prec];
		f["S0u"] = SetPrecision[(4 - pi)/(4 - 2*Sqrt[2+f["S0q"]]),prec];
		f["S0p"] = SetPrecision[1 - f["S0u"], prec];
		f["S0"] = SetPrecision[f["S0p"]*(a + b) + f["S0u"]*Sqrt[a^2 + b^2 + f["S0q"]*a*b], prec];    
		
		f["S1q"] = SetPrecision[92.28480788617108, prec];
		f["S1v"] = SetPrecision[0.04522028227769, prec];
		f["S1p2"] = SetPrecision[0.99983439391729, prec];
		f["S1u"] = SetPrecision[1 + f["S1v"] - f["S1p2"], prec];
		f["S1p1"] = SetPrecision[(pi/2)*(2 + f["S1v"]*Sqrt[2+ f["S1q"]]) - (2*f["S1p2"] + 2*f["S1u"]*Sqrt[2+f["S1q"]]), prec];
		
		(* Function S1 *)
		f["S1"] = SetPrecision[(f["S1p2"]*(a^2 + b^2) + f["S1p1"]*a*b + f["S1u"]*(a + b)*Sqrt[a^2 + b^2+ f["S1q"]*a*b])/
		                        ((a + b) + f["S1v"]*Sqrt[a^2 + b^2 + f["S1q"]*a*b]), prec];
	
		(* Constants for S2 *)
		f["S2q1"] = SetPrecision[13.66022044408346, prec];
		f["S2q2"] = SetPrecision[37.30921886231118, prec];
		f["S2v"] = SetPrecision[-1.03788930003090, prec];
		f["S2t"] = SetPrecision[5.51954143485218, prec];
		f["S2p2"] = SetPrecision[1.24957869093182, prec];
		f["S2u"] = SetPrecision[1 + f["S2v"] - f["S2p2"], prec];
		f["S2p1"] = SetPrecision[(pi/4)*(2 + f["S2t"] + f["S2v"]*Sqrt[2 + 2*f["S2q1"] + f["S2q2"]]) - 
		                          (f["S2p2"] + f["S2u"]*Sqrt[2 + 2*f["S2q1"] + f["S2q2"]]), prec];
		f["S2R"] = SetPrecision[Sqrt[(a^4 + b^4) + f["S2q1"]*a*b*(a^2 + b^2) + f["S2q2"]*a^2*b^2], prec];
		f["S2"] = SetPrecision[(f["S2p2"]*(a^3 + b^3) + f["S2p1"]*a*b*(a + b) + f["S2u"]*(a + b)*f["S2R"])/
		                        ((a^2 + b^2) + f["S2t"]*a*b + f["S2v"]*f["S2R"]), prec];

		f["S1rq"] = SetPrecision[74.01745125408363, prec];
		f["S1rv"] = SetPrecision[0.05027328304233, prec];
		f["S1ru"] = f["S1rv"];
		f["S1rp1"] = SetPrecision[(pi/2)*(2 + f["S1rv"]*Sqrt[2 + f["S1rq"]]) - (2 + 2*f["S1ru"]*Sqrt[2 + f["S1rq"]]), prec];
		f["S1r"] = SetPrecision[((a^2 + b^2) + f["S1rp1"]*a*b + f["S1ru"]*(a + b)*Sqrt[a^2 + b^2 + f["S1rq"]*a*b])/
		                         ((a + b) + f["S1rv"]*Sqrt[a^2 + b^2 + f["S1rq"]*a*b]), prec];
		
		f["Sac1"] = SetPrecision[Pi - 3, prec];
		f["Sac2"] = SetPrecision[Pi, prec];
		f["Sac3"] = SetPrecision[0.5, prec];
		f["Sac4"] = SetPrecision[(1 + pi)/2, prec];
		f["Sac5"] = SetPrecision[4, prec];
		f["Sak"] = SetPrecision[
		  1 - (f["Sac1"]*a*b)/(
		    (a^2 + b^2) + f["Sac2"]*Sqrt[f["Sac3"]*(a*b)^2 + a*b*Sqrt[a*b*(f["Sac4"]*(a^2 + b^2) + f["Sac5"]*a*b)]]
		  ), 
		  prec
		];

		(* Define the function S_a using constants from the association *)
		f["Sa"] = SetPrecision[
		  4 * ((pi*a*b + f["Sak"]*(a - b)^2)/(a + b)),
		  prec
		];     
		
		f["Saoc1"] = SetPrecision[0.14220038049945, prec];
		f["Saoc2"] = SetPrecision[3.30596250119242, prec];
		f["Saoc3"] = SetPrecision[0.00135657637724, prec];
		f["Saoc4"] = SetPrecision[2.00637978782056, prec];
		f["Saoc5"] = SetPrecision[5.3933761426286, prec];
		
		f["Saok"] = SetPrecision[
		  1 - (f["Saoc1"]*a*b)/(
		    (a^2 + b^2) + f["Saoc2"]*Sqrt[f["Saoc3"]*(a*b)^2 + a*b*Sqrt[a*b*(f["Saoc4"]*(a^2 + b^2) + f["Saoc5"]*a*b)]]
		  ), 
		  prec
		];

		(* Define the function S_a using constants from the association *)
		f["Sao"] = SetPrecision[
		  4 * ((pi*a*b + f["Saok"]*(a - b)^2)/(a + b)),
		  prec
		];                        	                                                                    	                                      	                                      
                           	                                      	                                                                    	                                      	                                                 	                                      	                                                                    	                                      	                                                	                                      	                                                                    	                                      	                                      
		f["Sard1"] = 0.14220096377128;
		f["Sard2"] = 3.93490847789660;
		f["Sard3"] = 2.691437204515743;
		f["Sark"] = SetPrecision[1-(f["Sard1"]*a*b)/((Power[a,2]+Power[b,2])+f["Sard2"]*Sqrt[Sqrt[Power[a*b,3]*(Power[a,2]+Power[b,2]+f["Sard3"]*a*b)]]),prec];
		f["Sar"] = SetPrecision[4*(pi*a*b+f["Sark"]*Power[a-b,2])/(a+b),prec];
		
		f["Mn"] = SetPrecision[ Power[f["MR"]/f["Mp"],148/a]*f["Mu"]*Power[f["Mp"]/f["Mu"],16/f["Mp"]],prec];

		
		f["Mcfr"] = SetPrecision[(-0.000541034776116510819*f["h"]+1.00000000884005358*f["M6"])+(-0.096787856655402682*f["h"]+0.000697814328424526341)/((-34.1723518423179868*f["M6"]-2.76394591853239291)+(1.34595617011879276*f["h"]-16.3886566776292*f["M6"]+16.5756304424531677)/((0.726375130874769348*a-0.204371274323853713*f["M6"]+0.378209596686388871)+(-51.0191665578208955*a+66.201015508602751*f["h"]+14.7882894581736259*f["M6"]-121.181029470213161)/((113.02054659946748)))),prec];
		f["Mn1"] = SetPrecision[ Power[a,f["h1"]/Sqrt[a]]*f["Sar"], prec];
		f["Mn2"] = SetPrecision[ Sqrt[f["M6"]*f["Sa"]], prec];
		f["Mn3"] = SetPrecision[ f["Sa"]/Power[f["Sa"]/f["M6"],0.360886], prec];
		f["Mn4"] = SetPrecision[ f["M6"]*Power[f["Mu"],Power[f["h2"],f["h2"]/f["h1"]]*(Sqrt[41/a]/a)], prec];
		f["Mn6"] = SetPrecision[ f["M6"]*(1+(Power[f["h"],120]/(a*a))), prec];
		f["Mn7"] = SetPrecision[ f["M6"]*(1+(Power[f["h"],91]/(a*Sqrt[a*39*Sqrt[a]]))),prec];
		f["Mn8"] = SetPrecision[ f["M6"]*(1+(Power[f["h"],20]/(a*a*(1+(f["MR"]/(-367401066*f["h2"])))))),prec];
		f["Mn9"] = SetPrecision[ f["M6"]*(1+(Power[f["h"],21]/(a*a*(1+(f["Mu"]/(-10)))))),prec];
		f["Mn10"] = SetPrecision[ f["M6"]*(1+(Power[f["h"],14]/(a*a*(1+(-32))))), prec];
		
		(* Probably the best candidate *)
		f["MMC"] = SetPrecision[ f["M6"]*(1+(Power[f["h"],18]/(a*a*(1+(-1*a))))), prec];
		f["err"] = SetPrecision[(f["y"]-f["MMC"])*10^15,prec];
		f["Mn12"] = SetPrecision[ f["M6"]*(1+((f["h"]^2/(a^2+267246))*(1+(19*(-1+f["h"]))))), prec];
		f["MA1"] = SetPrecision[ 6.44122277800000 + (f["a"] - 1.05000000000000)/(0.262839503508289 + (f["a"] - 9.)/(-883.186059069206 + (f["a"]- 10.)/(-0.0126877977870027 + (f["a"] - 25.)/(-106841.053594433 - 6591.70925503389*f["a"])))), prec];
		f["MA2"] = SetPrecision[ 6.44122277811400412 + (f["a"] - 1.05000000000000004)/(
    0.29666450294061225061345548085397649689403700831403554006642836042624810610605578554176575752449825035203152625408954906070724926578140660623 + (f["a"] - 
       1.75000000000000000)/(-127.26349228410171894690500849968774172785832429797525523802705789246006656499861591753638677782243197816614441082503221979809856955362167718062 + (f["a"] - 
          5.)/(-0.04663700142212516854516472648158800557985283876386537611117668894428022781390298955212243314541179083813612496310037646649173908043662722975976727434336584113400898285152156727178732402060616254145040 + (f["a"] - 
             500.)/(-1.135055152118334015386182340378333538639784621549089979720354304604491956876598789375630287931833705652647780796558172296470974045806160613982375439414246443662138525593883928583404839404781916862*10^6 - 
             36361.628558908774755148948001494216629448130019937942743268302509857392244377180674358345606012856740643717989049570071802301424199530678455637786326898109716971730621534825895891347073532932915581626*f["a"]))))
, prec];
		f["MA3"] = SetPrecision[6.44122277811400412 + (f["a"] - 1.05000000000000004)/(0.308074227724945333990283788085790343765364425466513046338239 + (f["a"] - 1.25000000000000000)/(-38.3725001538332332513872473700243520738530250464336326868800 + (f["a"] - 1.44999999999999996)/(-0.054977489868208621269882879106620863090766222686015460800803 + (f["a"] - 1.64999999999999991)/(-7442.34873328859863350948819810490389849089651352467406759070 + (f["a"] - 19.)/(-0.00302457649123594380521216484343666943551487813312798964101292 - 7.21613608673293174579455263322258134436491578183359859552094*10^(-14)*f["a"])))))
,prec];

f["MA4"] = SetPrecision[6.44122277811400412 + (f["MMC"] - 6.44122277811400412)/(0.999999999999639395820547030342709088509838558180644730177540 + (f["MMC"] - 8.80079054932078364)/(4.86976638890347583775776965568442648892709003712287725873807*10^7 + (f["MMC"] - 72.4201271945980238)/(-2.31790149399261070983378688359340617268243921604440091023237*10^(-6) + 5.79475464632943812565648563386999295558300150755221119223119*10^(-16)*f["MMC"])))

,prec];

(*
f["MA5"] = SetPrecision[ pi*(a+b)*(1.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 + f["h"]/(3.9735487424059950893908295838190842474135325224481680713362810500367712386642485513570478199726282695410143306514411616051755988224241363505016209385779491304585134157539091557793342694802373775922632 + (f["h"] - 0.10370583165756966388968687158862395863257684573398448721631715024418270611893134156851479459925308819304797471990807239299051996552714737144498707268026429187015225509910945130709566216604424016087331)/(-3.4901312695381510935834432334750764867213683579662899489053415817919963449480442300370064192524063878734016370578351827075025292293496181785991239486845484374864380542238644676362938244832129849315805 + (f["h"] - 0.51020408163265306122448979591836734693877551020408163265306122448979591836734693877551020408163265306122448979591836734693877551020408163265306122448979591836734693877551020408163265306122448979591837)/(1.4963585175575033595071335352493257395517438269871716747026516730112653080855049099922196029921918085895730751274918174428441950014784283356327841451294745198344736377869084248736357607870047512569569 - 0.72316929671977993446090144996072895429383051161316437497378119954691145734159840161666942695941155370408525743233442164607737484696394879069722263362061018501418582653061723632043525941389563974251443*f["h"])))
),prec];
*)
f["MA5"] = SetPrecision[pi*(a+b)*(1.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 + f["h"]/(3.9944312459433201854792656834599057564004983328748771079150632744336086766996818692158099764614839949717197902516091910689193749126701403621279686788372735298978068946150963472134922757343647534910083 + (f["h"] - 0.02218198279764599366229062924400181077410593028519692168401991851516523313716613852421910366681756450882752376641014033499320959710276143051154368492530556813037573562698053417836124943413309189678588)/(-3.9154234838942453340526494249558392703306364101655345327640015442610617721915949086605770465428035011241978479958883180251535102524888076986662077523104338228185207811456630107240208197026993386755257 + (f["h"] - 0.08895044629116651277316097260695598645737149892274546014158202523853493382579255155432440751000307787011388119421360418590335487842413050169282856263465681748230224684518313327177593105570944906124962)/(1.1034253078148054035243055754447147552300333340431623591657412513462498270152039000073412906275812898907831365927337393384544798106321145959586821493384116510646312811300106296342546329919722063628320 + (f["h"] - 0.66942148760330578512396694214876033057851239669421487603305785123966942148760330578512396694214876033057851239669421487603305785123966942148760330578512396694214876033057851239669421487603305785123967)/(-2.2407956071413060541473429359757432152055371780303295818845790821384717805364813675416515042559195027061380693489067368586847120763285445376022061107230691130685033510188921089322940356700269556861436 + (f["h"] - 0.81859410430839002267573696145124716553287981859410430839002267573696145124716553287981859410430839002267573696145124716553287981859410430839002267573696145124716553287981859410430839002267573696145125)/(1.3428431979228680490787340129823410645493447825699093922769655763888791045633089734152534922924455962081471484243830008671935653915446353458892536550286804134852033458191518846479242962895469448005111 - 0.95413882658148728619619098818265865477477577185576282372817886316912803798331899276677950034350733944929514486366994400606284025388106235030972480613968239065162570226727070275701665893851042797372742*f["h"]))))))
,prec];


(*f["err"] = SetPrecision[(f["err"] - Min[f["err"]]) / (Max[f["err"]] - Min[f["err"]]),prec];
f["yd"] = SetPrecision[f["err"],prec];*)
(*f["yd"] = SetPrecision[f["y"] / f["maxY"], prec];*)
(*f["yd"] = SetPrecision[ (f["y"]/f["MMC"]), prec];*)

       keys = {
       "MR", "Mu", "MRu", "Mp", 
       "M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9",
       "K1", "K2", "K3", "K4", "K5", "K6", "K7", "K8", "K9", "K10", "K11", "K12",
       "Ka","Kb","Kc","Kd","Ke",
       "E1","E2","E3","E4","E5","E6", "O2",
       "Ea","Eb","Ec","Ed","Ee",
       "C1","C2","C3","C4","C5","C6","C7",
       "A1","A2","A3","A4","A5","A6","A7",
       "S0","S1","S2","S1r","Sa","Sao","Sar",
       "Mn", "Mcfr", "Mn1", "Mn2", "Mn3", "Mn4", "Mn6", "Mn7", "Mn8", "Mn9", "Mn10", "Mn11", "Mn12", "MMC", "MA1", "MA2", "MA3", "MA4", "MA5"
       };
		
       calculateErrors = Function[key,
         f[key <> " Err"] = SetPrecision[Abs[f["y"] - f[key]]/f["y"], prec]
       ];
       Map[calculateErrors, keys];
       
       f
     ],
     {pair, data}
   ];
   results
];

End[];
EndPackage[];