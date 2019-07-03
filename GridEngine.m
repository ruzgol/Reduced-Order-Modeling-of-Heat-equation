function Grid = GridEngine( Grid, File )
%GridEngine Compute the grid information of the model, e.g., neighbors,
%volumes, connectivities, weights for gradient calculation, etc.
%   Grid:       Structure to store the grid information. 
%   File:       File system information for the model

    if (~isempty(Grid))
        return;
    else
        switch File.Type
            case 'Plot3D'
                fileName  = [File.Folder File.Name '.xyz'];
                Grid = ReadPlot3Dxyz(fileName, File);
                Grid.Type = 'Plot3D';
            case 'ACE'
                PathToFile  = [File.Folder File.Name];
                Grid = ReadACEModel(PathToFile);
                Grid.Type = 'ACE';
            otherwise
                error('Unkown file type!');
        end
    end
    
    % Save the Grid information into the model
    WriteOutData(File, Grid, 'Grid', 'mat');

end

