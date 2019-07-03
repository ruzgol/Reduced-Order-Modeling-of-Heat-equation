function Modes = PODEngine(Grid, Modes, File, modInfo)

% Purpose:  Compute the modes of pressure, density, velocity, etc.
% Input:    Grid Information to tell 2D/3D simulation
%           Modes to check if Modes is already available
%           File to indicate simulation information
% Output: modes and coefficients
    DataFileLoaded = 0;
    
    if (~isempty(Modes))
        return;
    else
        PathToFile  = [File.WorkingPath File.Name '_Data.mat'];  
        if(exist(PathToFile,'file'))
            load(PathToFile, '-mat');
            DataFileLoaded = 1;
        else
            % settings of IPOD
            nnode = Grid.nnode;

            % Read the simulation primitive data and perform incremental POD
            nstep   = 0;

            % Initialize the DataTray
            Tray.U = zeros(nnode, length(File.Step));
            Tray.V = zeros(nnode, length(File.Step));
            if (Grid.dim == 3)
                Tray.W = zeros(nnode, length(File.Step));
            end
            Tray.T = zeros(nnode, length(File.Step));
            

            % Reserve for future compressible flow uses
            % Tray.rho = zeros(nnode, File.Step);
            % Tray.T = zeros(nnode, File.Step);
            %Tray.P = zeros(nnode, length(File.Step));
            for istep = File.Step

                % load data file
                fileName = [File.Folder File.Name '_SimData_' num2str(istep) File.Suffix];
                [rawData, variables, ios] = ReadData(Grid, fileName, File);

                sprintf('%s read done ... \n', fileName)

                if (ios == 0)
                    nvariables = size(variables, 1);        
                    nstep = nstep + 1;

                    for ivar = 1 : nvariables
                        switch strtrim(variables(ivar, :))
                            case 'U'
                                Tray.U(:, nstep) = rawData(:, ivar);
                            case 'V'
                                Tray.V(:, nstep) = rawData(:, ivar);
                            case 'W'
                                Tray.W(:, nstep) = rawData(:, ivar);
                            case 'Pressure'
                                Tray.P(:, nstep) = rawData(:, ivar);
                            case 'Density'
                                Tray.P(:, nstep) = rawData(:, ivar);
                            case 'Temperature'
                                Tray.T(:, nstep) = rawData(:, ivar);
                            otherwise
                                error('Variables for simulation data is not wrong');
                        end
                    end
                end
            end
        end

        % Calculate POD Modes
        Modes = FlowPOD(Tray, Grid, modInfo.POD.Config);

        % Save the Solution Data into the model
        if(~DataFileLoaded)
            WriteOutData(File, Tray, 'SolData', 'mat');
        end
        
        % Save the Mode information into the model
        WriteOutData(File, Modes, 'Modes', 'mat');
     end
end
    


