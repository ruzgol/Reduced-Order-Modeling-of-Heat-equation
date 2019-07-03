function Term4 = ComputeTerm4( File,Grid, Modes, ModeGrads, modInfo )
%ComputeTerm4: Compute the forth term in Yi Wang's derivation: -u(dot)del(u0) 
%   Grid: Structure to store grid information
%   Modes: Structure to store mode information
%   ModeGrads: Structure to store mode gradient information
%   This is the $T_0.\nabla u_0$ term
    NModes=modInfo.POD.Config.NMODES;
    Term4   = zeros(NModes, NModes);

    dT0_dx = ModeGrads.T.Mean(:,1);
    % Add z-compnent as needed
    if (Grid.dim == 3)
       dT0_dz = ModeGrads.T.Mean(:,3);
       dT0_dy = ModeGrads.T.Mean(:,2);
    elseif(Grid.dim == 2)  
       dT0_dy = ModeGrads.T.Mean(:,2);
    end

    % Iterate to calculate the constant for Term1
    for i = 1:NModes
        for j = 1:NModes
            if(Grid.dim == 3)
               Bi = Modes.T.Modes(:,i);
               dBj_dx = ModeGrads.T.Modes(:,j,1);
               dBj_dy = ModeGrads.T.Modes(:,j,2);
               dBj_dz = ModeGrads.T.Modes(:,j,3);
               tmp = Bi.*(Bjx.*dT0x_dx+Bjy.*dT0x_dy+Bjz.*dT0x_dz);
               Term4(i,j) = -sum(tmp.*Grid.CellVolumes);
            elseif(Grid.dim == 2) 
               Bi = Modes.T.Modes(:,i);
               dBj_dx = ModeGrads.T.Modes(:,j,1);
               dBj_dy = ModeGrads.T.Modes(:,j,2);
               tmp = Bi.*(dBj_dx.*dT0_dx+dBj_dy.*dT0_dy);
               [Term4(i,j)] = Integrate(Grid,tmp);
               Term4(i,j) = -Term4(i,j);
            else
               Bi = Modes.T.Modes(:,i);
               dBjx_dx=ModeGrads.T.Modes(:,j);
               tmp = Bi.*(dBjx_dx.*dT0_dx);
               Term4(i,j) = -sum(tmp.*Grid.CellVolumes); 
            end
        end
    end

end

