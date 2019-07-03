function [ODEs,Modes] = ODEsEngine( Grid, Modes, ModeGrads, ODEs, File, modInfo )
%ODEsEngine Compute the coefficients of ROM to assemble ODE for computation
%   ModeGrad:   Structure to store the gradient of modes. 
%   Modes:      POD Modes
%   File:       Information for the simulation
%     if (~isempty(ODEs))
%         return;
%     end
     Grid.dim =1;
    if(Grid.dim == 1)
      [Grid,Modes,ModeGrads] = OneD(File,modInfo, Grid, Modes, ModeGrads);
    else(Grid.dim == 2)
      [Grid,Modes,ModeGrads] = TwoD(File,modInfo,Grid, Modes, ModeGrads);
    end  
    %Term1 = ComputeTerm1(File,Grid, Modes, ModeGrads, modInfo);
    Term2 = ComputeTerm2(File,Grid, Modes, ModeGrads, modInfo);  
    %Term3 = ComputeTerm3(File,Grid, Modes, ModeGrads, modInfo);
    %Term4 = ComputeTerm4(File,Grid, Modes, ModeGrads, modInfo);
    %Term5 = ComputeTerm5(File,Grid, Modes, ModeGrads, modInfo);
    Term6 = ComputeTerm6(File,Grid, Modes, ModeGrads, modInfo);
    
    % Save the GradMatrix information into the model
%     ODEs.Ai     = Term1 + Term2;
%     ODEs.Aij    = Term3 + Term4+ Term5+ Term6;
    
    ODEs.Ai     =  Term2;
    ODEs.Aij    =  Term6;
    WriteOutData(File, ODEs, 'ODEs', 'mat'); 
    
end