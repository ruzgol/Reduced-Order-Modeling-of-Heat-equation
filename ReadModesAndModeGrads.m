function [Modes, ModeGrads] = ReadModesAndModeGrads(File)

% Purpose: Read the modes and gradient of modes from .mat data
% Input: File type, name and folder
% Output: Modes and Mode Gradients information
    Modes = [];
    ModeGrads = [];
    
    PathToFile  = [File.WorkingPath File.Name '_ROM.mat'];
    if(exist(PathToFile, 'file'))
        try
            load(PathToFile, '-mat', 'Modes', 'ModeGrads');
        catch
            Modes = [];
            ModeGrads = [];
        end
    end 
end