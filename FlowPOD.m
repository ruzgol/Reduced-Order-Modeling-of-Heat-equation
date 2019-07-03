function Modes = FlowPOD(Tray, Grid, Config)
%FlowPOD Compute POD modes of the velocity and pressure using snapshot data
%and regular POD (i.e., batch POD) via SVD or EVD approach
%   SnapshotData: Data from which to extract POD modes
%   Config: Configuration parameters
%   Modes: The data structure contains the POD modes. Modes have several
%   subfield, "U", "V", "W", "P", "rho", and T
%   Each subfield will have subfield, including: Mean, Modes; Coefficients,
%   and Lambda (eigenvalues)
%   Tray: The snapshotdata
    Modes.U = [];
    Modes.V = [];
    Modes.W = [];
    
    % Temprature
    Modes.T             = PODScalar(Grid.CellVolumes, Config, Tray.T);
    if (Grid.dim == 2)
        [Modes.U, Modes.V]  = PODVector(Grid.CellVolumes, Config, Tray.U, Tray.V);
    else
        [Modes.U, Modes.V, Modes.W] = PODVecotor(Grid.CellVolumes, Config, ...
                                        Tray.U, Tray.V, Tray.W);  
    end

end

