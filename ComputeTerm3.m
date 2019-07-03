function Term3 = ComputeTerm3(File, Grid, Modes, ModeGrads, modInfo )
%ComputeTerm3: Compute the second term in Yi Wang's derivation: u0(dot)del(u) 
%   Grid: Structure to store grid information
%   Modes: Structure to store mode information
%   ModeGrads: Structure to store mode gradient information

%   Note that in the current version we didn't consider the boundary
%   integral, which will be inclued in the future version
%   This is the $u_0.\nabla T_$ term
    NModes=modInfo.POD.Config.NMODES;
    Term3   = zeros(NModes, NModes);

    u0x = Modes.U.Mean;
    u0y = Modes.V.Mean;
    % Add z-compnent as needed
    if (Grid.dim == 3)
      u0z = Modes.W.Mean;
    end

    % Iterate to calculate the constant for Term1
    for i = 1:NModes
        for j = 1:NModes
            if(Grid.dim == 3)
               Bi = Modes.T.Modes(:,i);
               dBj_dx = ModeGrads.T.Modes(:,j,1);
               dBj_dy = ModeGrads.T.Modes(:,j,2);
               dBj_dz = ModeGrads.T.Modes(:,j,3);
               tmp = Bix.*(u0x.*dBjx_dx+u0y.*dBjx_dy+u0z.*dBjx_dz);
               Term3(i,j) = -sum(tmp.*Grid.CellVolumes);
            elseif(Grid.dim == 2)
               Bi = Modes.T.Modes(:,i);
               dBj_dx = ModeGrads.T.Modes(:,j,1);
               dBj_dy = ModeGrads.T.Modes(:,j,2);
               tmp = Bi.*(u0x.*dBj_dx+u0y.*dBj_dy);
               [Term3(i,j)] = Integrate(Grid,tmp);
               Term3(i,j) = -Term3(i,j);
            else
               Bi = Modes.T.Modes(:,i);
               dBjx_dx=ModeGrads.T.Modes(:,j);
               tmp = Bi.*(u0x.*dBjx_dx);
               Term3(i,j) = -sum(tmp.*Grid.CellVolumes);   
            end
        end
    end

end

