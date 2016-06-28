(* ::Package:: *)

(* ::Subsection:: *)
(*Standard Definitions*)


\[Phi][x_]:=If[Abs[x]>1,0,1-Abs[x]]
\[Phi][l_Integer,i_Integer]:=If[l==-1,1,If[l==0,#,\[Phi][2^l (#-2^-l (2i-1))]]]&
\[Phi][lx_,ly_,i_,j_]:=\[Phi][lx,i][#1]\[Phi][ly,j][#2]&
\[Phi][l_List,i_List]:=Product[\[Phi][l[[d]],i[[d]]][#[[d]]],{d,1,Length[l]}]&
\[Psi][l_,i_]:=\[Phi][2^l (#1-2^-l i)]&
\[Psi][lx_,ly_,i_,j_]:=\[Phi][2^lx (#1-2^-lx i)]\[Phi][2^ly (#2-2^-ly j)]&
\[Psi][l_List,i_List]:=Product[\[Phi][2^l[[d]] (#[[d]]-2^-l[[d]] i[[d]])],{d,1,Length[l]}]&

proj[f_,x_]:=f[Join[#,{x}]]&

ceil[x_,l_Integer]:=If[x>=1,2^(l-1),1+Floor[2^(l-1) x]]
ceil[x_List,l_List]:=Table[If[x[[i]]>=1,2^(l[[i]]-1),1+Floor[2^(l[[i]]-1) x[[i]]]],{i,1,Length[x]}]

atPosition[x_List,l_List,i_List]:=Module[{j,prod,xpos},xpos=ceil[x,l];
											prod=Product[If[xpos[[j]]==i[[j]] || l[[j]]<1,1,0],{j,1,Length[x]}];
											If[prod==1,Return[True];,Return[False];];]
(*Perhaps this could be made more efficient*)

switch[x_]:=If[x>=1,x,1]
switch2[x_]:=If[x==-1,0,If[x==0,1,2^-x]]
switch3[x_,i_]:=If[x<=0,1,(-1)^(i+1) 2^-x]

area[kx_]:=Which[kx==-1,1,kx==0,1/2,kx>0,2^-kx]

column2[diag_,row_]:=If[diag<1,diag, If[row<1,diag,diag-row+1]]
column1[diag_,row_]:=If[diag==row,-1,column2[diag,row]]

getrow[perp_,diag_]:=Which[diag==-1,-1,diag==0,Min[perp,0],perp<=diag,perp,True,diag]
getcol[perp_,diag_]:=Which[diag==-1,-1,diag==0,Min[diag-perp,0],perp> 0,diag-perp,True,diag]

endperp[n_]:=Which[n==-1,-1,n==0,1,True,n+1]

d\[Phi][x_]:=Which[Abs[x]>1,0,
			x==1 ,-1/2,
			x==-1,1/2,
			x==0,0,
			x>0,-1,
			True,1](* Don't think I need this*)


(* ::Subsection:: *)
(*Standard Lagrange Basis*)


standardCoefficients[f_,l_List]:=If[Length[l]==1,Table[f[{j 2^-l[[1]]}],{j,0,2^l[[1]]}],Table[standardCoefficients[proj[f,j 2^-l[[-1]]],l[[;;-2]]],{j,0,2^l[[-1]]}]]
(*Working*)

standardReconstruct[coefficients_,l_List,x_List]:=If[Length[l]==1,Sum[
coefficients[[i+1]]  \[Psi][l[[1]],i] [x[[1]]],{i,0,2^l[[1]]}],
Sum[standardReconstruct[coefficients[[i+1]],l[[;;-2]],x[[;;-2]]] \[Psi][l[[-1]],i] [x[[-1]]],{i,0,2^l[[1]]}]
]
(*The hierarchical and sparse basis will be more efficient*)


standardProject[f_,l_List,x_]:=Module[{},coeffs=standardCoefficients[f,l];Return[standardReconstruct[coeffs,l,x]];]


(* ::Subsection:: *)
(*Hierarchical Iteration in n-D*)


hierIterate[l_List,ki_List]:=Module[{k,i,hashEntry,hashMap,count},
hashMap=Table[0,{i,0,2^l[[1]]}];
count=1;
If[Length[l]==1,
For[k=-1,k<=l[[1]],k++,
For[i=1,i<=2^(switch[k]-1),i++,hashEntry=Append[Append[ki,k+2],i];
hashMap[[count++]]=hashEntry;
]];
Return[Partition[Flatten[hashMap],2Length[l]]];
,
For[k=-1,k<=l[[1]],k++,
For[i=1,i<=2^(switch[k]-1),i++,
hashEntry=hierIterate[l[[2;;]],Append[Append[ki,k+2],i]];hashMap[[count++]]=hashEntry]];
Return[Partition[Flatten[hashMap],2Length[l]]];
]]
(*This is working*)

hierIterate[l_List]:=hierIterate[l,{}]
(*This is working*)


(* ::Subsection:: *)
(*Getting Hierarchical Coefficients*)


getCoefficient[f_,x_List,l_List]:=If[Length[l]==1,Which[l[[1]]==-1,f[{0}],l[[1]]==0,f[{1}]-f[{0}],True,f[x]-1/2 f[x+2^-l[[1]]]-1/2 f[x-2^-l[[1]]]],
				Which[l[[-1]]==-1,getCoefficient[proj[f,0],x[[;;-2]],l[[;;-2]]],
					l[[-1]]==0,getCoefficient[proj[f,1],x[[;;-2]],l[[;;-2]]]-getCoefficient[proj[f,0],x[[;;-2]],l[[;;-2]]],
					True,getCoefficient[proj[f,x[[-1]]],x[[;;-2]],l[[;;-2]]]-1/2 getCoefficient[proj[f,x[[-1]]+2^-l[[-1]]],x[[;;-2]],l[[;;-2]]]-1/2getCoefficient[proj[f,x[[-1]]-2^-l[[-1]]],x[[;;-2]],l[[;;-2]]]
				]
]

pointIndex[x_,l_Integer]:=Module[{},Table[If[k<=1,1,ceil[x,k]],{k,-1,l}]]
pointIndex[x_List,l_List]:=Table[pointIndex[x[[i]],l[[i]]],{i,1,Length[l]}]
spointIndex[x_List,n_Integer]:=Table[pointIndex[x[[i]],n],{i,1,Length[x]}]



(* ::Subsection:: *)
(*Full Basis*)


fullCoefficients[f_,l_List]:=Module[{li,x,i,k,diff,hash,count},li=hierIterate[l];hash=Table[0,{i,1,Length[li]}];
count=1;
For[i=1,i<=Length[li],i++,
k=li[[i]][[1;;;;2]]-2;
x=2^-k (2li[[i]][[2;;;;2]]-1);
diff=getCoefficient[f,x,k];
hash[[count++]]=li[[i]]-> diff;];
Return[SparseArray[hash]];
]
(*This is working*)

fullReconstruct[coeffs_,x_,l_]:=Module[{a,k,value,i,d,indices,ki},
value=0;
d=Length[x]; 
indices=pointIndex[x,l];
Do[
ki=Flatten[Table[{k[a]+2,indices[[a,k[a]+2]]},{a,1,d}]];
value+=(coeffs[[##]]&@@Sequence[ki])\[Phi][ki[[1;;;;2]]-2,ki[[2;;;;2]]][x];
,##]&
@@Sequence[Table[{k[i],-1,l[[i]]},{i,1,d}]];
Return[value];
]
(*This is working,but you can make it slightly more efficient by*)

(*
fullReconstruct[li_,x_]:=Module[{a,k,value,i},value=0;
For[a=1,a\[LessEqual] Length[li],a++,
k=Keys[li[[a]]][[1;;;;2]];
i=Keys[li[[a]]][[2;;;;2]];
If[atPosition[x,k-2,1/2 (i+1)],
value+=Values[li[[a]]] \[Phi][k-2,1/2 (i+1)][x]]
];Return[value];]

(*This is working, but you can make it much more effici*)ent*)


(* ::Subsection:: *)
(*Sparse Iteration in n-D*)


sparseIterate[n_,d_,ki_List]:=Module[{k,i,hashEntry,hashMap,count},
hashMap=Table[0,{i,0,2^n}];
count=1;
If[d==1,
For[k=-1,k<=n,k++,
For[i=1,i<=2^(switch[k]-1),i++,
hashEntry=Append[Append[ki,k+2],i];
hashMap[[count++]]=hashEntry;
]];
Return[hashMap];
,
For[k=-1,k<=n,k++,
For[i=1,i<=2^(switch[k]-1),i++,
hashEntry=sparseIterate[n-k-1,d-1,Append[Append[ki,k+2],i]];hashMap[[count++]]=hashEntry;]];Return[hashMap],2d(*All the work with Flatten/Partition is like.. no cost*)
]]
(*This is working*)

sparseIterate[n_,d_]:=Partition[Flatten[sparseIterate[n,d,{}]],2d]
(*This is working*)

sparsekIterate[n_,d_,l_List]:=Module[{k,hashEntry,hashMap,count},
If[d==1,
count=1;
hashMap=Table[0,{i,1,n+2}];
For[k=-1,k<=n,k++,
hashEntry=Append[l,k+2];
hashMap[[count++]]=hashEntry;
];Return[hashMap]
,
count=1;
hashMap=Table[0,{i,1,(n+2)}];
For[k=-1,k<=n,k++,
hashEntry=sparsekIterate[n-k-1,d-1,Append[l,k+2]];hashMap[[count++]]=hashEntry;];
Return[hashMap]
]]

sparsekIterate[n_,d_]:=Partition[sparsekIterate[n,d,{}]//Flatten,d]


(* ::Subsection:: *)
(*Sparse Basis*)


sparseCoefficients[f_,n_,d_]:=Module[{li,hash,k,i,diff,x,count},li=sparseIterate[n,d];
hash=Table[0,{i,1,Length[li]}];
count=1;
For[i=1,i<=Length[li],i++,
k=li[[i]][[1;;;;2]]-2;
x=2^-k (2li[[i]][[2;;;;2]]-1);
diff=getCoefficient[f,x,k];
hash[[count++]]=li[[i]]-> diff;];
Return[SparseArray[hash]];
]
(*Working efficiently*)



sparseReconstruct[coefficients_,x_,n_]:=Module[{a,k,value,i,d,indices,ki,iters},value=0;d=Length[x]; indices=spointIndex[x,n];
iters=sparsekIterate[n,d];

For[i=1,i<=Length[iters],i++,
ki=Flatten[Table[{(iters[[i]])[[a]],indices[[a,(iters[[i]])[[a]]]]},{a,1,d}]];
value+=(coefficients[[##]]&@@Sequence[ki])\[Phi][ki[[1;;;;2]]-2,ki[[2;;;;2]]][x];
];

Return[value];
]
(*Finally Working*)

sparseReconstruct[coefficients_,x_]:=fullReconstruct[coefficients,x]
(*Finally Working*)

reverseFlattenSparse[coefficients_,n_]:=Module[{k},k=0;Return[Table[Table[Table[coefficients[[++k]] ,{i,1,2^switch[getrow[perp,diag]]-1,2},{j,1,2^switch[getcol[perp,diag]]-1,2}],{perp,-1,endperp[diag]}],{diag,-1,n}]];]



(* ::Subsection:: *)
(*Differentiation*)


diffx[f_,l_Integer,x_List]:=Module[{ans,temp},temp=x;Which[x[[1]]< 2^-l,temp[[1]]+=2^-l; ans=f[temp]-f[x]; ans*=2^l;,
x[[1]]>1-2^-l ,temp[[1]]-=2^-l; ans=f[x]-f[temp]; ans*=2^l;,
True,temp[[1]]+=2^-l;ans=f[temp];temp[[1]]-=2 2^-l; ans-=f[temp];ans*=(2^(l-1));];Return[ans];]
(*This is working.. I think it's efficient. We could do a general one: *)

diff[f_,l_Integer,x_List,i_:1]:=Module[{ans,temp},temp=x;Which[x[[i]]< 2^-l,temp[[i]]+=2^-l; ans=f[temp]-f[x]; ans*=2^l;,
x[[i]]>1-2^-l ,temp[[i]]-=2^-l; ans=f[x]-f[temp]; ans*=2^l;,
True,temp[[i]]+=2^-l;ans=f[temp];temp[[i]]-=2 2^-l; ans-=f[temp];ans*=(2^(l-1));];Return[ans];]

(*And just in case, a gradient:*)
grad[f_,l_List,x_List]:=Table[diff[f,l[[i]],x,i],{i,1,Length[x]}]



(* ::Input:: *)
(*(*Past this point idk... *)*)


(* ::Subsection:: *)
(*Showing the Grids*)


makeFull[lx_,ly_]:=Table[Table[{switch2[kx] i,switch2[ky] j},{i,1,2^switch[kx]-1,2},{j,1,2^switch[ky]-1,2}],{kx,-1,lx},{ky,-1,ly}]


makeSparse[n_]:=Table[Table[Table[{switch2[getrow[perp,diag]] i,switch2[getcol[perp,diag]] j},{i,1,2^switch[getrow[perp,diag]]-1,2},{j,1,2^switch[getcol[perp,diag]]-1,2}],{perp,-1,endperp[diag]}],{diag,-1,n}]


(* ::Subsection:: *)
(*Integration on the Grids*)


posW[x_]:=If[(x==0||x==1),1/2,1]


fullIntegrate[f_,kx_,ky_]:=Sum[2^-kx 2^-ky  posW[i]posW[j] f[i,j],{i,0,1,2^-kx},{j,0,1,2^-ky}]

fullIntegrate[f_,kx_]:=Sum[2^-kx  posW[i]f[i],{i,0,1,2^-kx}]

sparseIntegrate[f_,n_]:=Sum[Sum[Sum[posW[switch2[getrow[perp,diag]] i]posW[switch2[getcol[perp,diag]] j]switch2[getcol[perp,diag]]switch2[n]f[switch2[getrow[perp,diag]] i,switch2[getcol[perp,diag]] j],{i,1,2^switch[getrow[perp,diag]],2},{j,1,2^switch[getcol[perp,diag]]-1,2}],{perp,-1,endperp[diag]}],{diag,-1,n}]

sparseIntegrateByCoefficients[f_,n_]:=Sum[Sum[Sum[getCoefficient2D[f,switch2[getrow[perp,diag]] i,switch2[getcol[perp,diag]] j,getrow[perp,diag],getcol[perp,diag]] area[getrow[perp,diag]] area[getcol[perp,diag]],{i,1,2^switch[getrow[perp,diag]],2},{j,1,2^switch[getcol[perp,diag]],2}]
,{perp,-1,endperp[diag]}],{diag,-1,n}]


sparseBoundaryXIntegrate[f_,n_]:=0 (*do later, if you find a better definition*)

(*Sum[Sum[posW[i]posW[j] 1/ GCD[2^n i,2^(n)]1/2^n f[i, j ],{j,0,1,1/GCD[2^n i,2^(n)]}],{i,0,1,2^(-n)}]*)


(*The area weights, although they sum to 1, are not totally right for sparse grids*)


(* ::Subsection:: *)
(*Weights for the Sparse Grid*)


weights1D[i_,kx_,lx_]:=If[Mod[i,2^(lx-kx)]==0,switch3[kx,2^(kx-lx) i],0]

showFullWeights[lx_,ly_]:=Sum[Table[weights1D[i,kx,lx] weights1D[j,ky,ly],{i,0,2^lx},{j,0,2^ly}],{kx,0,lx},{ky,0,ly}]

showSparseWeights[n_]:=Table[weights1D[i,0,n] weights1D[j,0,n],{i,0,2^n},{j,0,2^n}]+Sum[Sum[Table[weights1D[i,Max[0,getrow[perp,diag]],n]weights1D[j,Max[0,getcol[perp,diag]],n],{i,0,2^n},{j,0,2^n}],{perp,-1,endperp[diag]}],{diag,1,n}]


(* ::Subsection:: *)
(*Redefining differentiation (Sparse Differentiation)*)


sparseDiffx[f_,n_,kx_,ky_,i_,j_]:=Which[#1< 2^-n,(f[#1+2^-l,#2]-f[#1,#2])/2^-l,#1>1-2^-n ,(f[#1,#2]-f[#1-2^-l,#2])/2^-l,True,(f[#1+2^-l,#2]-f[#1-2^-l,#2])/(2 2^-l) ]&


sparseDiffy[f_,l_]:=
Which[#2< 2^-l,(f[#1,#2+2^-l]-f[#1,#2])/2^-l,#2> 1-2^-l,(f[#1,#2]-f[#1,#2-2^-l])/2^-l,True,(f[#1,#2+2^-l]-f[#1,#2-2^-l])/(2 2^-l)]&
