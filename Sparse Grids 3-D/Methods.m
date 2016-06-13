(* ::Package:: *)

\[Phi][x_]:=If[Abs[x]>1,0,1-Abs[x]]
\[Phi][l_,i_]:=If[l==-1,1,If[l==0,#,\[Phi][2^l (#-2^-l (2i-1))]]]&
\[Phi][lx_,ly_,i_,j_]:=\[Phi][lx,i][#1]\[Phi][ly,j][#2]&
\[Psi][lx_,ly_,i_,j_]:=\[Phi][2^lx (#1-2^-lx i)]\[Phi][2^ly (#2-2^-ly j)]&
ceil[x_,l_]:=If[x>=1,2^(l-1),1+Floor[2^(l-1) x]]
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



standardCoefficients2D[f_,lx_,ly_]:=Table[f[2^-lx i,2^-ly j],{i,0,2^lx},{j,0,2^ly}]
standardReconstruct2D[coefficients_,lx_,ly_]:=Sum[coefficients[[i+1,j+1]] \[Psi][lx,ly,i,j][#1,#2],{i,0,2^lx},{j,0,2^ly}]&

standardProject2D[f_,lx_,ly_]:=standardReconstruct2D[standardCoefficients2D[f,lx,ly],lx,ly]


getCoefficient[f_,x_,l_]:=If[l==-1,f[0],If[l==0,f[1]-f[0],f[x]-1/2. f[x+2^-l]-1/2. f[x-2^-l]]]

getCoefficient2D[f_,x_,y_,lx_,ly_]:=getCoefficient[Function[y1,getCoefficient[Function[x1,f[x1,y1]],x,lx]],y,ly]
(* Tensor product construction
*)

fullCoefficients2D[f_,lx_,ly_]:=Table[Table[getCoefficient2D[f,2^-kx i,2^-ky j,kx,ky],{i,1,2^switch[kx]-1,2},{j,1,2^switch[ky]-1,2}],{kx,-1,lx},{ky,-1,ly}]


Reconstruct2D[coefficients_]:= Sum[coefficients[[kx+2,ky+2,ceil[#1,switch[kx]],ceil[#2,switch[ky]]]] 
\[Phi][kx,ky,ceil[#1,switch[kx]],ceil[#2,switch[ky]]][#1,#2]
,{kx,-1,Length[coefficients]-2},{ky,-1,Length[coefficients]-2}]&

reverseFlatten[coefficients_,kx_,ky_]:=Module[{k},k=0;Return[Table[Table[coefficients[[++k]] ,{i,1,2^(switch[lx]-1)},{j,1,2^(switch[ly]-1)}],{lx,-1,kx},{ly,-1,ky}]];]



project2D[f_,lx_,ly_]:=Reconstruct2D[fullCoefficients2D[f,lx,ly]]


sparseCoefficients2D[f_,n_]:=Table[Table[Table[getCoefficient2D[f,switch2[getrow[perp,diag]] i,switch2[getcol[perp,diag]] j,getrow[perp,diag],getcol[perp,diag]],{i,1,2^switch[getrow[perp,diag]]-1,2},{j,1,2^switch[getcol[perp,diag]]-1,2}],{perp,-1,endperp[diag]}],{diag,-1,n}]

sparseReconstruct[coefficients_]:=Sum[Sum[coefficients[[diag+2,perp+2,ceil[#1,switch[getrow[perp,diag]]],ceil[#2,switch[getcol[perp,diag]]]]] 
\[Phi][getrow[perp,diag],getcol[perp,diag],ceil[#1,switch[getrow[perp,diag]]],ceil[#2,switch[getcol[perp,diag]]]][#1,#2]
,{perp,-1,endperp[diag]}],{diag,-1,Length[coefficients]-2}]&

sparseProject[f_,n_]:=sparseReconstruct[sparseCoefficients2D[f,n]]


reverseFlattenSparse[coefficients_,n_]:=Module[{k},k=0;Return[Table[Table[Table[coefficients[[++k]] ,{i,1,2^switch[getrow[perp,diag]]-1,2},{j,1,2^switch[getcol[perp,diag]]-1,2}],{perp,-1,endperp[diag]}],{diag,-1,n}]];]



diffx[f_,l_]:=Which[#1< 2^-l,(f[#1+2^-l,#2]-f[#1,#2])/2^-l,#1>1-2^-l ,(f[#1,#2]-f[#1-2^-l,#2])/2^-l,True,(f[#1+2^-l,#2]-f[#1-2^-l,#2])/(2 2^-l) ]&


diffy[f_,l_]:=
Which[#2< 2^-l,(f[#1,#2+2^-l]-f[#1,#2])/2^-l,#2> 1-2^-l,(f[#1,#2]-f[#1,#2-2^-l])/2^-l,True,(f[#1,#2+2^-l]-f[#1,#2-2^-l])/(2 2^-l)]&


dxMatrix[lx_,ly_]:=Flatten[Table[Table[fullCoefficients2D[diffx[\[Phi][kx,ky,i,j],lx],lx,ly]//Flatten,{i,1,2^(switch[kx]-1)},{j,1,2^(switch[ky]-1)}],{kx,-1,lx},{ky,-1,ly}],3]//Transpose
dyMatrix[lx_,ly_]:=Flatten[Table[Table[fullCoefficients2D[diffy[\[Phi][kx,ky,i,j],lx],lx,ly]//Flatten,{i,1,2^(switch[kx]-1)},{j,1,2^(switch[ky]-1)}],{kx,-1,lx},{ky,-1,ly}],3]//Transpose


dxMatrixStandard[lx_,ly_]:=Flatten[Table[standardCoefficients2D[diffx[\[Psi][lx,ly,i,j],lx],lx,ly]//Flatten,{i,0,2^lx},{j,0,2^ly}],1]
dyMatrixStandard[lx_,ly_]:=Flatten[Table[standardCoefficients2D[diffy[\[Psi][lx,ly,i,j],lx],lx,ly]//Flatten,{i,0,2^lx},{j,0,2^ly}],1]


dxMatrixSparse[n_]:=Flatten[Table[Table[Table[
sparseCoefficients2D[diffx[\[Phi][getrow[perp,diag],getcol[perp,diag],i,j],n],n]//Flatten,{i,1,2^(switch[getrow[perp,diag]]-1)},{j,1,2^(switch[getcol[perp,diag]]-1)}],{perp,-1,endperp[diag]}],{diag,-1,n}],3]//Transpose

dyMatrixSparse[n_]:=Flatten[Table[Table[Table[
sparseCoefficients2D[diffy[\[Phi][getrow[perp,diag],getcol[perp,diag],i,j],n],n]//Flatten,{i,1,2^(switch[getrow[perp,diag]]-1)},{j,1,2^(switch[getcol[perp,diag]]-1)}],{perp,-1,endperp[diag]}],{diag,-1,n}],3]//Transpose



makeFull[lx_,ly_]:=Table[Table[{switch2[kx] i,switch2[ky] j},{i,1,2^switch[kx]-1,2},{j,1,2^switch[ky]-1,2}],{kx,-1,lx},{ky,-1,ly}]


makeSparse[n_]:=Table[Table[Table[{switch2[getrow[perp,diag]] i,switch2[getcol[perp,diag]] j},{i,1,2^switch[getrow[perp,diag]]-1,2},{j,1,2^switch[getcol[perp,diag]]-1,2}],{perp,-1,endperp[diag]}],{diag,-1,n}]


posW[x_]:=If[(x==0||x==1),1/2,1]


fullIntegrate[f_,kx_,ky_]:=Sum[2^-kx 2^-ky  posW[i]posW[j] f[i,j],{i,0,1,2^-kx},{j,0,1,2^-ky}]

fullIntegrate[f_,kx_]:=Sum[2^-kx  posW[i]f[i],{i,0,1,2^-kx}]

sparseIntegrate[f_,n_]:=Sum[Sum[Sum[posW[switch2[getrow[perp,diag]] i]posW[switch2[getcol[perp,diag]] j]switch2[getcol[perp,diag]]switch2[n]f[switch2[getrow[perp,diag]] i,switch2[getcol[perp,diag]] j],{i,1,2^switch[getrow[perp,diag]],2},{j,1,2^switch[getcol[perp,diag]]-1,2}],{perp,-1,endperp[diag]}],{diag,-1,n}]

sparseIntegrateByCoefficients[f_,n_]:=Sum[Sum[Sum[getCoefficient2D[f,switch2[getrow[perp,diag]] i,switch2[getcol[perp,diag]] j,getrow[perp,diag],getcol[perp,diag]] area[getrow[perp,diag]] area[getcol[perp,diag]],{i,1,2^switch[getrow[perp,diag]],2},{j,1,2^switch[getcol[perp,diag]],2}]
,{perp,-1,endperp[diag]}],{diag,-1,n}]


sparseBoundaryXIntegrate[f_,n_]:=0 (*do later, if you find a better definition*)

(*Sum[Sum[posW[i]posW[j] 1/ GCD[2^n i,2^(n)]1/2^n f[i, j ],{j,0,1,1/GCD[2^n i,2^(n)]}],{i,0,1,2^(-n)}]*)


(*The area weights, although they sum to 1, are not totally right for sparse grids*)


weights1D[i_,kx_,lx_]:=If[Mod[i,2^(lx-kx)]==0,switch3[kx,2^(kx-lx) i],0]

showFullWeights[lx_,ly_]:=Sum[Table[weights1D[i,kx,lx] weights1D[j,ky,ly],{i,0,2^lx},{j,0,2^ly}],{kx,0,lx},{ky,0,ly}]

showSparseWeights[n_]:=Table[weights1D[i,0,n] weights1D[j,0,n],{i,0,2^n},{j,0,2^n}]+Sum[Sum[Table[weights1D[i,Max[0,getrow[perp,diag]],n]weights1D[j,Max[0,getcol[perp,diag]],n],{i,0,2^n},{j,0,2^n}],{perp,-1,endperp[diag]}],{diag,1,n}]


lagrangeW[lx_,ly_]:=Flatten[Table[Table[fullIntegrate[\[Psi][lx,ly,i1,j1][#1,#2] \[Psi][lx,ly,i2,j2][#1,#2]&,lx,ly],{i2,0,2^lx},{j2,0,2^ly}]//Flatten,{i1,0,2^lx},{j1,0,2^ly}],1]

lagrangedW[lx_,ly_]:=Flatten[Table[Table[fullIntegrate[\[Psi][lx,ly,i1,j1][1,#1] \[Psi][lx,ly,i2,j2][1,#1]-\[Psi][lx,ly,i1,j1][0,#1] \[Psi][lx,ly,i2,j2][0,#1]&,ly],{i2,0,2^lx},{j2,0,2^ly}]//Flatten,{i1,0,2^lx},{j1,0,2^ly}],1]

fullW[lx_,ly_]:=Flatten[Table[Table[Table[Table[fullIntegrate[\[Phi][kx1,ky1,i1,j1][#1,#2] \[Phi][kx2,ky2,i2,j2][#1,#2]&,lx,ly],{i2,1,2^(switch[kx2]-1)},{j2,1,2^(switch[ky2]-1)}],{kx2,-1,lx},{ky2,-1,ly}]//Flatten,{i1,1,2^(switch[kx1]-1)},{j1,1,2^(switch[ky1]-1)}],{kx1,-1,lx},{ky1,-1,ly}],3]

fulldW[lx_,ly_]:=Flatten[Table[Table[Table[Table[fullIntegrate[\[Phi][kx1,ky1,i1,j1][1,#] \[Phi][kx2,ky2,i2,j2][1,#]-\[Phi][kx1,ky1,i1,j1][0,#] \[Phi][kx2,ky2,i2,j2][0,#]&,ly],{i2,1,2^(switch[kx2]-1)},{j2,1,2^(switch[ky2]-1)}],{kx2,-1,lx},{ky2,-1,ly}]//Flatten,{i1,1,2^(switch[kx1]-1)},{j1,1,2^(switch[ky1]-1)}],{kx1,-1,lx},{ky1,-1,ly}],3]

sparseW[n_]:=Flatten[Table[Table[Table[Table[Table[Table[
sparseIntegrateByCoefficients[\[Phi][getrow[perp1,diag1],getcol[perp1,diag1],i1,j1][#1,#2]\[Phi][getrow[perp2,diag2],getcol[perp2,diag2],i2,j2][#1,#2]&,n],{i2,1,2^(switch[getrow[perp2,diag2]]-1)},{j2,1,2^(switch[getcol[perp2,diag2]]-1)}],{perp2,-1,endperp[diag2]}],{diag2,-1,n}]//Flatten,{i1,1,2^(switch[getrow[perp1,diag1]]-1)},{j1,1,2^(switch[getcol[perp1,diag1]]-1)}],{perp1,-1,endperp[diag1]}],{diag1,-1,n}],3]

sparsedW[n_]:=Flatten[Table[Table[Table[Table[Table[Table[
sparseBoundaryXIntegrate[\[Phi][getrow[perp1,diag1],getcol[perp1,diag1],i1,j1][#1,#2]\[Phi][getrow[perp2,diag2],getcol[perp2,diag2],i2,j2][#1,#2]&,n],{i2,1,2^(switch[getrow[perp2,diag2]]-1)},{j2,1,2^(switch[getcol[perp2,diag2]]-1)}],{perp2,-1,endperp[diag2]}],{diag2,-1,n}]//Flatten,{i1,1,2^(switch[getrow[perp1,diag1]]-1)},{j1,1,2^(switch[getcol[perp1,diag1]]-1)}],{perp1,-1,endperp[diag1]}],{diag1,-1,n}],3]


sparseDiffx[f_,n_,kx_,ky_,i_,j_]:=Which[#1< 2^-n,(f[#1+2^-l,#2]-f[#1,#2])/2^-l,#1>1-2^-n ,(f[#1,#2]-f[#1-2^-l,#2])/2^-l,True,(f[#1+2^-l,#2]-f[#1-2^-l,#2])/(2 2^-l) ]&


sparseDiffy[f_,l_]:=
Which[#2< 2^-l,(f[#1,#2+2^-l]-f[#1,#2])/2^-l,#2> 1-2^-l,(f[#1,#2]-f[#1,#2-2^-l])/2^-l,True,(f[#1,#2+2^-l]-f[#1,#2-2^-l])/(2 2^-l)]&
