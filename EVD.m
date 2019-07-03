function [phi, A, sqrt_eValc, V] = EVD(DeMeanData, Weight, Config)
%EVD Compute the POD modes, Coefficients, and Lambda for the retained modes
%using eigen-value decomposition
%   DeMeanData:     Data with mean removed 
%   Weight:         Weights to weight the demeaned data
%   Config:         POD Configuration parameters. 
    
    [~, Ns ]   = size(DeMeanData);
    WeightedData = bsxfun(@times, DeMeanData, sqrt(Weight));
    sqrt_Nt    = sqrt(Ns);

    % construct the empirical correlation matrix
    C = (WeightedData'*WeightedData)/Ns;

    % Calculate the right eight vector
    [eVectC,Lambda] = eig(C, 'vector');
    [Lambda,I]  = sort(Lambda, 'descend');
    eVectC      = eVectC(:, I);
    
    % Determine the POD Order
    Cut     = PODOrder(Lambda, Config); 

    % Truncate the eigen values at the dominant order
    Lambda  = Lambda(1:Cut); 

    % Snapshot data projected onto the eigenvector with truncated eigenvector
    % Non-truncated version: S       = SnapshotData * eVectC
    V       = eVectC(:,1:Cut);
    S       = (DeMeanData * V)/sqrt_Nt;

    sqrt_eValc = sqrt(Lambda);
    inv_sqrt_eValc = 1./sqrt(Lambda);

    % Orthonormal POD modes with respect to volume (weights)
    phi = bsxfun(@times, S, inv_sqrt_eValc');

    % Temporal POD coefficients
%     A   = sqrt_Nt*sqrt_eValc*eVectC(:,1:TOrd)';
    A  = sqrt_Nt*bsxfun(@times, V', sqrt_eValc);

fprintf('*** POD Engine using EVD is completed *** \n')