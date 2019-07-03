function aidot = dai_dt_LamFlowROM(t, aj, ODEs)
%dai_dt_LamFlow: The function to calculate the dai_dt given aj and ODEs
%   t:          Time
%   aj:         The POD coefficients
%   ODEs:       The vector, matrix, and tensor for calculating dai_dt

    % Unpack the matrix
    Ai      = ODEs.Ai;
    Aij     = ODEs.Aij;
    
    aidot = Ai + Aij*aj;

end