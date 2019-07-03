function WriteOutData( File, Data, DataType, fmt )
%WriteOutData Write out the data according to specified DataType and format
%   modInfo: Model information structure
%   Data: The data structure
%   DataType: {Grid, Modes, ModeGrads, GradMatrix, PODModes, ODEs}
%   fmt: {'mat', 'bin'}
    
    PathToFile  = [File.WorkingPath File.Name '_ROM.mat'];
    
    if (strcmp(fmt, 'mat'))
        switch DataType
            case 'Grid'
                Grid = Data;
                save(PathToFile, 'Grid', '-mat');
            case 'ModeGrads'
                ModeGrads = Data;
                save(PathToFile, 'ModeGrads', '-mat', '-append');
            case 'Modes'
                Modes = Data;
                save(PathToFile, 'Modes', '-mat', '-append');
            case 'GradMatrix'
                GradMatrix = Data;
                save(PathToFile, 'GradMatrix', '-mat', '-append');
            case 'ODEs'
                ODEs = Data;
                save(PathToFile, 'ODEs', '-mat', '-append');
            case 'SolData'
                Tray = Data;
                PathToFile = [File.WorkingPath File.Name '_Data.mat'];
                save(PathToFile, 'Tray', '-mat');
            case 'GradData'
                GradTray = Data;
                PathToFile = [File.WorkingPath File.Name '_GradData.mat'];
                save(PathToFile, 'GradTray', '-mat','-v7.3');
            otherwise
                warning('Unexpected output datatype. No data is written out.')
        end
        
    else
        disp('output format .bin will be implemented in the future')
    end

end

