function varargout = PODVector( Weight, Config, varargin )
%PODScale Calculate the mean, POD modes, and coefficient of a scalar
%variable
%   Data: The snapshot data for the scalar variable
%   Weight: The weight for POD
%   Config: The configuration for POD
    nVar = size(varargin, 2);
    [nnode, ~] = size(varargin{1,1});
    
    if (nVar == 2)
        Data = [varargin{1,1}; varargin{1,2}];
        Weight = repmat(Weight, 2, 1);
    else
        Data = [varargin{1,1}; varargin{1,2}; varargin{1,3}];
        Weight = repmat(Weight, 3, 1);
    end
    Mean = mean(Data, 2);
    
    DeMeanData = bsxfun(@minus, Data, Mean);

    
    switch Config.Method
        case 'SVD'    
            [Modes, Coeffs, Sigma, Z] = SVD(DeMeanData, Weight, Config);
        case 'EVD'        
            [Modes, Coeffs, Sigma, Z] = EVD(DeMeanData, Weight, Config);
        otherwise
            error('Unexpected method for POD');
    end
    
    
    for k = 1:nVar
        id = (k-1)*nnode+1 : k*nnode;
        varargout{k}.Mean       = Mean(id,:);
        varargout{k}.Modes      = Modes(id,:);
        % All vector components share common POD coefficients
        varargout{k}.Coeffs     = Coeffs;
        varargout{k}.Sigma      = Sigma;
        varargout{k}.Z          = Z;
    end
    
end

