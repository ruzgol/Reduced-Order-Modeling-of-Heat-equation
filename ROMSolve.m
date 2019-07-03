function [t, ai] = ROMSolve( Modes, ODEs, File, modInfo )
%ROMSolve Compute the POD coefficients as a function of time or parameters
%   Modes:      POD Modes
%   ODEs:       ODE coefficients
%   File:       Information for the simulation
    
    % Extract the initial condition of the POD coefficients
%     st=Modes.T.Sigma;
%     vt=Modes.T.coeffs;
%     vt=vt';
%     temp=st*vt;
    
    ai0 = Modes.T.a0(:,modInfo.Simu.nt);
    ai0 = Modes.T.a0(:,1);
    options = [];
    %options = odeset('Jacobian',-1e9,'Stats','on');
    [t,ai]              = ode45(@dai_dt_LamFlowROM,[0:(modInfo.Phys.dt):modInfo.Simu.tf], ...
                                 ai0, options, ODEs);

                             
    WriteOutData(File, ODEs, 'ODEs', 'mat'); 
    
end