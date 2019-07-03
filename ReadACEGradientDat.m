function [rawData, variables, ios] = ReadACEGradientDat(Grid, fileName)

    % Read ACE solution binary file

    ios = 0;
    variables = [];

    if exist(fileName, 'file')
        rawData     = importdata(fileName);
        variables   = char(Grid.GradDataVar);
    else
        ios = 1;
        fprintf('%s is not existing. \n', fileName);
    end
end
