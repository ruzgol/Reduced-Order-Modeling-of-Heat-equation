function  [rawData, variables, ios] = ReadData(Grid, fileName, File)

% Purpose: load data file from simulation
% Input:  fileName
%         grid  -- Information about grid
% Output: rawData -- content in data file
%         variables -- the varialbes in the data fil

    if(strcmp(File.Suffix, '.bin'))
        [rawData, variables, ios]  =  ReadACESolution(Grid, fileName);
    elseif(strcmp(File.Suffix, '.dat'))
        [rawData, variables, ios]  =  ReadACESolutionDat(Grid, fileName);
    else
        error('ACE Solution File.Suffix is wrong !');
    end
end