function Grid = ReadGrid(File)
% Purpose: Read the grid information from .mat data
% Input: File type, name and folder
% Output: Grid information

    Grid = [];
    PathToFile  = [File.WorkingPath File.Name '_ROM.mat'];
    if(exist(PathToFile, 'file'))
        try
            load(PathToFile, '-mat', 'Grid');
        catch
            Grid = [];
        end
    end 
end