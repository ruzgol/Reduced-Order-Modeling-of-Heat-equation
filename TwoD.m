function [Grid,Modes,ModeGrads] = TwoD(File,modInfo,Grid, Modes, ModeGrads)
Grid=[];
load Running_space\points.mat;
[Grid.dim,Grid.nnode]=size(p);
Grid.coord(:,1)=p(1,:);
Grid.coord(:,2)=p(2,:);
load Running_space\Triangles.mat
[~,Grid.nelem]=size(t);
Grid.LEM=t(1:3,1:Grid.nelem);
[Grid] = Area_tri(Grid)
Modes.U.Mean=0;
Modes.V.Mean=0
dt=modInfo.Phys.dt;
%t=modInfo.Phys.t;
k=modInfo.Phys.nu;
x=Grid.coord(:,1);
y=Grid.coord(:,2);
[~,nt]=size(File.Step);
NModes=modInfo.POD.Config.NMODES;
%%%%%%%%%%%%%%%%%%%% Time Avaraging of the modes %%%%%%%%%%%%%%%%%%%%%%%%%%
load Running_space\solution.mat
T_exact=T;
for j=1:nt
[T_dx_exact(:,j,1), T_dx_exact(:,j,2), zxx, zyy, zxy] = trigradient2(x,y,T_exact(:,j),(Grid.LEM)',1);    
end                  
%%%%%%%%%%%%%%%%%%%%  Calculate POD modes of solution and its deruvative%%%% 
FOMSnap=[];
FOMSnap_dx=[];
FOMSnap_dy=[];
for i=1:nt
     FOMSnap=[FOMSnap,T_exact(:,i)-mean(T_exact,2)];
     FOMSnap_dx=[FOMSnap_dx,T_dx_exact(:,i,1)-mean(T_dx_exact(:,:,1),2)];
     FOMSnap_dy=[FOMSnap_dy,T_dx_exact(:,i,1)-mean(T_dx_exact(:,:,2),2)];
end
% FOMSnap=FOMSnap';
% FOMSnap_dx=FOMSnap_dx';
[Utemp,Stemp,Vtemp]=svd(FOMSnap,'econ');
Modes.T.Modes=Utemp(:,1:NModes);
Modes.T.Sigma=Stemp(1:NModes,1:NModes);
Modes.T.coeffs=Vtemp(:,1:NModes);
[Udxtemp,st,vt]=svd(FOMSnap_dx,'econ');
ModeGrads.T.Modes(:,1:NModes,1)=Udxtemp(:,1:NModes);
[Vdxtemp,st,vt]=svd(FOMSnap_dy,'econ');
ModeGrads.T.Modes(:,1:NModes,2)=Vdxtemp(:,1:NModes);
%%%%%%%%%%%%%%%%%%%%  Calculate POD modes of solution and its deruvative%%%% 
ModeGrads.T.Mean(:,1)=mean(ModeGrads.T.Modes(:,:,1),2);
ModeGrads.T.Mean(:,2)=mean(ModeGrads.T.Modes(:,:,2),2);
Modes.T.Mean=mean(Modes.T.Modes,2);

tempxxx(Grid.nnode,1)=0;
Modes.U.Mean=tempxxx;
Modes.V.Mean=tempxxx;
%Modes.T.a0=(Modes.T.Modes)'*T_exact(:,1);
Modes.T.a0=(Modes.T.Modes)'*T_exact;

load Running_space\edges.mat
Grid.Boundry=e(1,:);
end



% st=Modes.T.Sigma;
% vt=Modes.T.coeffs;
% vt=vt';
% temp=st*vt;
% T=Modes.T.Modes*ai';
% plot(x,-T_exact(:,10))