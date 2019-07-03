function [Grid,Modes,ModeGrads] = OneD(File,modInfo,Grid, Modes, ModeGrads)
Grid=[];
Grid.dim=1;
Grid.nnode=41;
Grid.dim=1;
L=1;
dx=L/(Grid.nnode-1);
Grid.coord(:,Grid.dim)=(0:dx:L);
Grid.CellVolumes=ones(Grid.nnode,1)*dx;
Modes.U.Mean=0;

%%%%%%%%%%%%%%%%%%%% Exact Solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    I.C.  T(x,0)=sin(2*pi*x/L)
%    B.C.L  T(0,t)=0
%    B.C.R  T(L,t)=0
%    T=sin(2*pi*x/L)*exp(-k*(2*pi/L)^2)*t
%%%%%%%%%%%%%%%%%%%% Exact Solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt=modInfo.Phys.dt;
t=modInfo.Phys.t;
T_exact(Grid.nnode,File.Step)=0;
k=modInfo.Phys.nu;
x=Grid.coord(:,Grid.dim);
[~,nt]=size(File.Step);
NModes=modInfo.POD.Config.NMODES;
%%%%%%%%%%%%%%%%%%%% Time Avaraging of the modes %%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:nt
    T_exact(:,i)=sin(2*pi*x/L).*exp(-k*(2*pi/L)^2*t(i));
    T_dx_exact(:,i)=2*pi/L*cos(2*pi*x/L).*exp(-k*(2*pi/L)^2*t(i));
end
%%%%%%%%%%%%%%%%%%%%  Calculate POD modes of solution and its deruvative%%%% 
FOMSnap=[];
FOMSnap_dx=[];
for i=1:nt
     FOMSnap=[FOMSnap,T_exact(:,i)-mean(T_exact,2)];
     FOMSnap_dx=[FOMSnap_dx,T_dx_exact(:,i)-mean(T_dx_exact,2)];
end
FOMSnap=FOMSnap';
FOMSnap_dx=FOMSnap_dx';
[Utemp,Stemp,Vtemp]=svd(FOMSnap,'econ');
Modes.T.Modes=Utemp(1:NModes,:);
Modes.T.Sigma=Stemp(1:NModes,1:NModes);
Modes.T.coeffs=Vtemp(1:NModes,:);
[Udxtemp,st,vt]=svd(FOMSnap_dx,'econ');
ModeGrads.T.Modes=Udxtemp(1:NModes,:);

% 
% [Utemp,Stemp,Vtemp]=svd(FOMSnap,'econ');
% Modes.T.Modes=Utemp(:,1:NModes);
% Modes.T.Sigma=Stemp(1:NModes,1:NModes);
% Modes.T.coeffs=Vtemp(:,1:NModes);
% [Udxtemp,st,vt]=svd(FOMSnap_dx,'econ');
% ModeGrads.T.Modes=Udxtemp(:,1:NModes);
%%%%%%%%%%%%%%%%%%%%  Calculate POD modes of solution and its deruvative%%%% 
ModeGrads.T.Mean=mean(ModeGrads.T.Modes,2);
Modes.T.Mean=mean(Modes.T.Modes,2);
[nmodes,~]=size(Modes.T.Mean);
tempxxx(nmodes,1)=0;
Modes.U.Mean=tempxxx;
Modes.T.a0=(Modes.T.Modes)*T_exact;
end



% st=Modes.T.Sigma;
% vt=Modes.T.coeffs;
% vt=vt';
% temp=st*vt;
% T=Modes.T.Modes*ai';
% plot(x,-T_exact(:,10))