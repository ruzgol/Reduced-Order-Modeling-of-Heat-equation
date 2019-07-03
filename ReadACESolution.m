function [rawData, variables, ios] = ReadACESolution(Grid, fileName)

% Read ACE solution binary file

ios = 0;
variables = [];

if exist(fileName, 'file')
    fid     = fopen(fileName,'r');
    hr1     = fread(fid,1,'int32');
    rawData = fread(fid,[Grid.nnode,Grid.nDataVar],'double');
    hr1     = fread(fid,1,'int32');
    fclose(fid);
    variables = char(Grid.DataVar);
else
    ios = 1;
    fprintf('%s is not existing. \n', fileName);
end
    