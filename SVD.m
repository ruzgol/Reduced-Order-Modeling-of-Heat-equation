function [phi, A, sigma, V] = SVD(DeMeanData, Weight, Config)
%SVD Compute the POD modes, Coefficients, and Lambda for the retained modes
%using Singular Value Decomposition
%   DeMeanData:     Data with mean removed 
%   Weight:         Weights to weight the demeaned data
%   Config:         POD Configuration parameters.
    [~, Ns] = size(DeMeanData);
    WeightedData = bsxfun(@times, DeMeanData, sqrt(Weight));
    sqrt_Nt    = sqrt(Ns);
    sqrt_weight     = sqrt(Weight);
    inv_sqrt_weight = 1./sqrt_weight;
    
    % Singular value decomposition of the weighted snapshot data
    [phi,sigma,V] = svd(WeightedData/sqrt_Nt,'econ');

    % extract the singular values from signular matrix
    S               = diag(sigma);

    % Obtian the eigenvalues from singular values
    Lambda          = S.*S;

    % Determine the POD Order
    Cut             = PODOrder(Lambda, Config);  
    
    % Truncate the singular values at the dominant order
    sigma           = S(1:Cut); 

    % Obtain the weighted and trunctaed POD modes
    phi             = bsxfun(@times, phi(:,1:Cut), inv_sqrt_weight);

    % Obtain the truncated V matrix
    V               = V(:,1:Cut);
    
    % Temporal POD coefficients
    A               = sqrt_Nt*bsxfun(@times, V', sigma);
    
    fprintf('*** POD using SVD is completed ***\n\n')
end
