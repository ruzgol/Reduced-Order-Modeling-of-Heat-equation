function GradMatrix = ReadGradMatrix(File)

% Purpose: Read the modes and gradient of modes from .mat data
% Input: File type, name and folder
% Output: Modes and Mode Gradients information
    GradMatrix = [];
    
    PathToFile  = [File.WorkingPath File.Name '_ROM.mat'];
    if(exist(PathToFile, 'file'))
        try
            load(PathToFile, '-mat', 'GradMatrix');
        catch
            GradMatrix = [];
        end
    end 
end