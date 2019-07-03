function Grid = ReadACEModel(PathToFile)

% Read Grid
fileName = [PathToFile '_SimInfo.dat'];
if exist(fileName, 'file')
    raw              = importdata(fileName,' ',11);
    Grid.nnode      = str2num(raw.textdata{2});
    Grid.CellVolumes     = raw.data;
    geomIdx  = str2num(raw.textdata{8});
    if (geomIdx == 1)
        Grid.dim = 3;
    else
        Grid.dim = 2;
    end

    % fprintf('%s has been loaded. \n', fileName);
else
    error('%s is not existing. \n', fileName);
end

% Read model information
fileName = [PathToFile '_SimInfoSupplement.dat'];
if exist(fileName, 'file')
    fid       = fopen(fileName);
    
    if (Grid.dim == 3)
        dataArray = textscan (fid,'%d64 %f %f %f %d64 %d64 %d64 %d64 %d64 %d64','HeaderLines',2);
        fclose(fid);
        Grid.coord(:,1)    = dataArray{2};
        Grid.coord(:,2)    = dataArray{3};
        Grid.coord(:,3)    = dataArray{4};
        Grid.neighbor(:,1) = dataArray{5};
        Grid.neighbor(:,2) = dataArray{6};
        Grid.neighbor(:,3) = dataArray{7};
        Grid.neighbor(:,4) = dataArray{8};
        Grid.neighbor(:,5) = dataArray{9};
        Grid.neighbor(:,6) = dataArray{10};
    else      
        dataArray = textscan (fid,'%d64 %f %f %d64 %d64 %d64 %d64','HeaderLines',2);
        fclose(fid);
        Grid.coord(:,1)    = dataArray{2};
        Grid.coord(:,2)    = dataArray{3};
        Grid.neighbor(:,1) = dataArray{4};
        Grid.neighbor(:,2) = dataArray{5};
        Grid.neighbor(:,3) = dataArray{6};
        Grid.neighbor(:,4) = dataArray{7};      
    end

else
    fprintf('%s is not existing. \n', fileName);
    stop;
end


fileName = [PathToFile '_variable_info.dat'];
if exist(fileName, 'file')
    fid     = fopen(fileName, 'r');
    txtline = textscan(fid, '%d', 1);
    Grid.nDataVar = txtline{1};
    txtline = textscan(fid, '%s', Grid.nDataVar);
    Grid.DataVar = txtline{1};
    % Check if Grad Data needs to be read
    if(~feof(fid))
        txtline = textscan(fid, '%d', 1);
        Grid.nGradDataVar = txtline{1};
        txtline = textscan(fid, '%s', Grid.nGradDataVar);
        Grid.GradDataVar = txtline{1};
    end
else
    fprintf('%s is not existing. \n', fileName);
    error('code stopped');
end