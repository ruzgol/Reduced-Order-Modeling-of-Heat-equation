function [rawData, variables, ios] = ReadACESolutionDat(Grid, fileName)

    % Read ACE solution binary file

    ios = 0;
    variables = [];

    if exist(fileName, 'file')
        rawData     = importdata(fileName);
        variables   = char(Grid.DataVar);
    else
        ios = 1;
        fprintf('%s is not existing. \n', fileName);
    end
end
