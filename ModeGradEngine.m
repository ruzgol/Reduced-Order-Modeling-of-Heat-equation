function ModeGrads = ModeGradEngine( Grid, Modes, ModeGrads, File )
%ModeGradEngine Compute the gradient of modes
%   ModeGrads:   Structure to store the gradient of modes. 
%   Modes:      POD Modes
%   File:       Information for the simulation
    if (~isempty(ModeGrads))
        return;
    end
    
    % Calculate the gradients of the modes using the Least Squares Method
%     ModeGrads = ModeGradEngineLeastSquares( Grid, Modes, ModeGrads, File );
    
    % Calculate the gradients of the modes using the dimensional gradient 
    % data output by the CFD solvers
    ModeGrads = ModeGradEngineFromGradData( Grid, Modes, ModeGrads, File );
    
    % Save the GradMatrix information into the model
    WriteOutData(File, ModeGrads, 'ModeGrads', 'mat'); 
    
end

