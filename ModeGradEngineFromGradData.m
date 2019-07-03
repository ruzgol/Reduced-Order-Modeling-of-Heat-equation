function ModeGrads = ModeGradEngineFromGradData( Grid, Modes, ModeGrads, File )
%ModeGradEngineFromGradData Calculate the gradient of the modes using the
%gradients of the dimensional data, e.g., du/dx, du/dy, dv/dx, dv/dy, etc. 
%   Detailed explanation goes here
    DataFileLoaded = 0;
    PathToFile  = [File.WorkingPath File.Name '_GradData.mat'];  
    if(exist(PathToFile,'file'))
        load(PathToFile, '-mat');
        DataFileLoaded = 1;
    else
        % settings of IPOD
        nnode = Grid.nnode;

        % Read the simulation primitive data and perform incremental POD
        nstep   = 0;

        % Initialize the DataTray
        if (Grid.dim == 3)
            GradTray.U = zeros(nnode, length(File.Step), 3);
            GradTray.V = zeros(nnode, length(File.Step), 3);
            GradTray.W = zeros(nnode, length(File.Step), 3);
            GradTray.P = zeros(nnode, length(File.Step), 3);
        else
            GradTray.U = zeros(nnode, length(File.Step), 2);
            GradTray.V = zeros(nnode, length(File.Step), 2);
            GradTray.P = zeros(nnode, length(File.Step), 2);
        end

        
        % Reserve for future compressible flow uses
        % GradTray.rho = zeros(nnode, File.Step, 3);
        % GradTray.T = zeros(nnode, File.Step, 3);
        
        for istep = File.Step

            % load data file
            fileName = [File.Folder File.Name '_GRAD_SimData_' num2str(istep) File.Suffix];
            [rawData, variables, ios] = ReadGradData(Grid, fileName, File);

            sprintf('%s read done ... \n', fileName)

            if (ios == 0)
                nvariables = size(variables, 1);        
                nstep = nstep + 1;

                for ivar = 1 : nvariables
                    colrange = (ivar-1)*Grid.dim+1:ivar*Grid.dim;
                    switch strtrim(variables(ivar, :))    
                        case 'U'
                            GradTray.U(:, nstep, :) = rawData(:, colrange);
                        case 'V'
                            GradTray.V(:, nstep, :) = rawData(:, colrange);
                        case 'W'
                            GradTray.W(:, nstep, :) = rawData(:, colrange);
                        case 'Pressure'
                            GradTray.P(:, nstep, :) = rawData(:, colrange);
                        case 'Density'
                            GradTray.P(:, nstep, :) = rawData(:, colrange);
                        case 'Temperature'
                            GradTray.T(:, nstep, :) = rawData(:, colrange);
                        otherwise
                            error('Variables for simulation data is not wrong');
                    end
                end
            end
        end
    end

  
    % Compute the gradients of modes using the gradients of the dimensional
    % data
    ModeGrads.P = ModeGradsFromGradData( Grid, Modes.P, GradTray.P );
    
    % Compute the gradients for the mean value and modes of U velocity;
    ModeGrads.U = ModeGradsFromGradData( Grid, Modes.U, GradTray.U );
    
    % Compute the gradients for the mean value and modes of V velocity;
    ModeGrads.V = ModeGradsFromGradData( Grid, Modes.V, GradTray.V );
    
    if (Grid.dim == 3)
        % Compute the gradients for the mean value and modes of W velocity;
        ModeGrads.V = ModeGradsFromGradData( Grid, Modes.W, GradTray.W );
    end
    
    % Save the GradTray Data into the model
    if(~DataFileLoaded)
        WriteOutData(File, GradTray, 'GradData', 'mat');
    end
        
end

