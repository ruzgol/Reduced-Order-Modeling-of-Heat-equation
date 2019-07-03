function [rawData, variables, ios] = ReadACEGradient(Grid, fileName)

% Read ACE solution binary file

ios = 0;
variables = [];

if exist(fileName, 'file')
    fid     = fopen(fileName,'r');
    hr1     = fread(fid,1,'int32');
    if(Grid.dim == 3)
        rawData = fread(fid,[Grid.nnode,Grid.nGradDataVar*3],'double');
    else
        rawData = fread(fid,[Grid.nnode,Grid.nGradDataVar*2],'double');
    end
    hr1     = fread(fid,1,'int32');
    fclose(fid);
    variables = char(Grid.DataVar);
else
    ios = 1;
    fprintf('%s is not existing. \n', fileName);
end
    