function ModeGrad = ModeGradsFromGradData( Grid, ModesData, GradData )
%ModeGradsFromGradData Calculate the gradients of the mode using the
%gradients of the dimensional data for each individual variable
%   Grid: Structure to store the grid information
%   ModesData: Structure to store the mode data
%   GradDat: Structure to store the gradients of the dimensional data
    Ns      = size(GradData, 2);
    Nt      = 1/sqrt(Ns);
    Z       = ModesData.Z;
    Sigma   = ModesData.Sigma;
    
    for i = 1:Grid.dim
        Mean = mean(GradData(:,:,i), 2);
        DeMeanData = bsxfun(@minus, GradData(:,:,i), Mean);
        ModeGrad.Mean(:,i) = Mean;
        ModeGrad.Modes(:,:,i) = bsxfun(@times, Nt*DeMeanData*Z, (1./Sigma)');
    end

end

