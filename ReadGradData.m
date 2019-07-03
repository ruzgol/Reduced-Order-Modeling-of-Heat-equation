function  [rawData, variables, ios] = ReadGradData(Grid, fileName, File)

% Purpose: load gradient data file from simulation
% Input:  fileName
%         grid  -- Information about grid
% Output: rawData -- content in data file
%         variables -- the varialbes in the data fil

switch Grid.Type   
    case 'ACE'
        if(strcmp(File.Suffix, '.bin'))
            [rawData, variables, ios]  =  ReadACEGradient(Grid, fileName);
        elseif(strcmp(File.Suffix, '.dat'))
            [rawData, variables, ios]  =  ReadACEGradientDat(Grid, fileName);
        else
            error('ACE Solution File.Suffix is wrong !');
        end
    otherwise
        error('Unknow file type');
end