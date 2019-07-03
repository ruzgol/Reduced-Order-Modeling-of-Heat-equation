function ODEs = ReadODEs(File)

% Purpose: Read the ODE equations (previous generated)
% Input: File type, name and folder
% Output: ODE equations information
    ODEs = [];
    PathToFile  = [File.WorkingPath File.Name '_ROM.mat'];
    if(exist(PathToFile, 'file'))
        try
            load(PathToFile, '-mat', 'ODEs');
        catch
            ODEs = [];
        end
    end 
end
