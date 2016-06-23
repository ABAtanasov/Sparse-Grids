(* ::Package:: *)

switch[x_]:=If[x>0,x,0];


legendre[n_,x_]:=LegendreP[n-1,2 x-1] Sqrt[2 (n-1)+1]
h[1,1,x_]:=If[Abs[x]>1,0,If[x>=0,1/Sqrt[2],-(1/Sqrt[2])]]
h[2,1,x_]:=If[Abs[x]>1,0,If[x>=0,Sqrt[3/2] (-1+2 x),h[2,1,-x]]]
h[2,2,x_]:=If[Abs[x]>1,0,If[x>=0,Sqrt[1/2] (-2+3 x),-h[2,2,-x]]]
h[3,1,x_]:=If[Abs[x]>1,0,If[x>=0,1/3 Sqrt[1/2] (1-24 x+30 x^2),-h[3,1,-x]]]
h[3,2,x_]:=If[Abs[x]>1,0,If[x>=0,1/2 Sqrt[3/2] (3-16 x+15 x^2),h[3,2,-x]]]
h[3,3,x_]:=If[Abs[x]>1,0,If[x>=0,1/3 Sqrt[5/2] (4-15 x+12 x^2),-h[3,3,-x]]]
h[4,1,x_]:=If[Abs[x]>1,0,If[x>=0, Sqrt[15/34] (1+4x-30x^2+28 x^3),h[4,1,-x]]]
h[4,2,x_]:=If[Abs[x]>1,0,If[x>=0, Sqrt[1/42] (-4+105x-300x^2+210 x^3),-h[4,2,-x]]]
h[4,3,x_]:=If[Abs[x]>1,0,If[x>=0,1/2 Sqrt[35/34] (-5+48x-105x^2+64 x^3),h[4,3,-x]]]
h[4,4,x_]:=If[Abs[x]>1,0,If[x>=0, 1/2 Sqrt[5/42] (-16+105x-192x^2+105 x^3),-h[4,4,-x]]]
h[5,1,x_]:=If[Abs[x]>1,0,If[x>=0,Sqrt[1/186] (1+30x+210x^2-840 x^3+630x^4),-h[5,1,-x]]]
h[5,2,x_]:=If[Abs[x]>1,0,If[x>=0,1/2 Sqrt[1/38] (-5-144x+1155x^2-2240 x^3+1260x^4),h[5,2,-x]]]
h[5,3,x_]:=If[Abs[x]>1,0,If[x>=0,Sqrt[35/14694] (22-735x+3504x^2-5460 x^3+2700x^4),-h[5,3,-x]]]
h[5,4,x_]:=If[Abs[x]>1,0,If[x>=0,1/8 Sqrt[21/38] (35-512x+1890x^2-2560 x^3+1155x^4),h[5,4,-x]]]
h[5,5,x_]:=If[Abs[x]>1,0,If[x>=0,1/2 Sqrt[7/158] (32-315x+960x^2-1155 x^3+480x^4),-h[5,5,-x]]]


fk[k_,fnumber_]:=h[k,fnumber,#]&


vk[k_,fnumber_,level_,place_]:=If[level==0,legendre[fnumber,#1]&,2^(level/2) fk[k,fnumber][2^level #1-(2 place-1)]&]


fk[k_,fnumber_]:=h[k,fnumber,#]&


vk[k_,fnumber_,level_,place_]:=If[level==0,legendre[fnumber,#1]&,2^(level/2) fk[k,fnumber][2^level #1-(2 place-1)]&]


innerProduct[f_,g_,level_,place_]:=NIntegrate[f[x] g[x],{x,(place-1) 2^-switch[level-1],place 2^-switch[level-1]}]

iteratek[k_,n_]:=Module[{fnumber,level,place,hashMap,count},
hashMap=Table[0,{i,1,k 2^n}];count=1;
For[level=0,level<=n,level++,
For[place=1,place<=2^switch[level-1],place++,
For[fnumber=1,fnumber<=k,fnumber++,
hashMap[[count++]]={fnumber,level+1,place};
]]];

Return[hashMap];]
(*working*)


getCoefficientsk[k_,f_,fnumber_,level_,place_]:=innerProduct[f,vk[k,fnumber,level,place],level,place]

fullCoefficientsk[k_,f_,n_]:=Module[{iters,i,hash,count,coeff},
iters=iteratek[k,n];
hash=Table[0,{i,1,Length[iters]}];
For[i=1,i<=Length[iters],i++,
coeff=getCoefficientsk[k,f,iters[[i]][[1]],iters[[i]][[2]]-1,iters[[i]][[3]]];hash[[i]]=iters[[i]]->coeff;];Return[SparseArray[hash]]]

reconstructk[k_,coeffs_,x_,n_]:=Module[{fnumber,level,place,value},
value=0;
For[level=0,level<=n,level++,
For[place=1,place<=2^switch[level-1],place++,
For[fnumber=1,fnumber<=k,fnumber++,
value+=coeffs[[fnumber,level+1,place]] vk[k,fnumber,level,place][x];
]]];Return[value]]
