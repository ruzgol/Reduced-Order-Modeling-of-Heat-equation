function Cut = PODOrder( Lambda, Config )
%PODOrder Determine the POD order (number of basis) based on POD
%configuration or energy ratio
%   Config:     Configuration parameters for POD
%   Lambda:     Eigenvalues to determine the number of POD modes
    
    % Individual Energy Level
    cumLambda = cumsum(Lambda)/sum(Lambda);

    if(isfield(Config, 'Order') && ~isempty(Config.Order) && Config.Order > 0)
        Cut = Config.Order;
    elseif (isfield(Config, 'EnergyRatio') && ~isempty(Config.EnergyRatio) ...
            && Config.EnergyRatio > 0 && Config.EnergyRatio < 1)
        % Find the order to trucate using energy ratio
        Cut = find(cumLambda>Config.EnergyRatio, 1, 'first');
    else
        error('Method to determine number of POD modes is not set correctly...');
    end


end

