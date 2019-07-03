%#########################################################################%
% Galerkin POD solution for heat tranfer problems 
% Roozbeh Golshan @ University of South Carolina                                  %
% Email your questions to GolshanR@cec.sc.edu
%#########################################################################%

% Main codes for reading data and calculation
clc;
close all;
clear all;
format short;

% Add paths
addpath (genpath('./'));

%#########################################################################%
File.WorkingPath =  './Running_space/'; % Any user-specified path
File.Name    = 'lbstep';
File.Type    =  'ACE';
File.Suffix  = '.dat';
%File.Folder  = './Running_space/CFD_Data/Backstep/';
File.Step    = 1:1:1001;


% ROM Configuration
modInfo.POD.Config.EnergyRatio = 0.9999;
modInfo.POD.Config.Order = 0;
modInfo.POD.Config.Method = 'EVD'; % or 'EVD'
modInfo.POD.Config.NMODES=10;
% Simulation Parameters
modInfo.Phys.nu = 1.0e-2;
modInfo.Phys.dt = 0.01;
modInfo.Phys.t=(modInfo.Phys.dt)*(File.Step);
modInfo.Simu.tf = modInfo.Phys.t(end);
modInfo.Simu.nt=1000;
% Initialize the POD analysis
[ Grid, Modes, ModeGrads, ODEs ] = Initialize(File);

% Compute the grid information if "Grid" is not loaded
%Grid        = GridEngine(Grid, File);

% POD Modes
%Modes       = PODEngine(Grid, Modes, File, modInfo);

% POD Mode gradient
%ModeGrads    = ModeGradEngine(Grid, Modes, ModeGrads, File);

% Calculate ROM coefficients
[ODEs,Modes]        = ODEsEngine(Grid, Modes, ModeGrads, ODEs, File, modInfo); 

% ROM Computation
[tt ai] = ROMSolve(Modes, ODEs, File, modInfo);

Solution=(Modes.T.Modes)*ai'; 
Grid.dim=1
if (Grid.dim == 1)
    L=1;
    dx=1/40;
    x=(0:dx:L);
    k=modInfo.Phys.nu;
    for i=1:length(modInfo.Phys.t)
        T_exact(:,i)=sin(2*pi*x/L).*exp(-k*(2*pi/L)^2*tt(i));
    end
    nt=modInfo.Simu.nt;
    h1=plot(x,Solution(:,nt),'-.')
    hold
    h2=plot(x,T_exact(:,nt),':')
    [h1,h2]=legend('G-ROM','Exact')
    ylabel('Temperature')
    xlabel('$x_1$','interpreter','latex')
    title('t=1.5 sec k=0.01 number of modes 2')
elseif (Grid.dim == 2) 
x=Grid.coord(:,1);
y=Grid.coord(:,2);
nt=modInfo.Simu.nt;
[X,Y]=meshgrid(x,y,Solution(:,nt));

end    