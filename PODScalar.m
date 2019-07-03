function Modes = PODScalar( Weight, Config, Data )
%PODScale Calculate the mean, POD modes, and coefficient of a scalar
%variable
%   Data: The snapshot data for the scalar variable
%   Weight: The weight for POD
%   Config: The configuration for POD
    Modes.Mean = mean(Data, 2);
    DeMeanData = bsxfun(@minus, Data, Modes.Mean);

    
    switch Config.Method
        case 'SVD'    
            [Modes.Modes, Modes.Coeffs, Modes.Sigma, Modes.Z] = SVD(DeMeanData, Weight, Config);
        case 'EVD'        
            [Modes.Modes, Modes.Coeffs, Modes.Sigma, Modes.Z] = EVD(DeMeanData, Weight, Config);
        otherwise
            error('Unexpected method for POD');
    end
    
end

