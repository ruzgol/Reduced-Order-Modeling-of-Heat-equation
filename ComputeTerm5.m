function Term5 = ComputeTerm5( File,Grid, Modes, ModeGrads, modInfo )
%ComputeTerm5: Compute the fifth term in Yi Wang's derivation: u(dot)del(u) 
%   Grid: Structure to store grid information
%   Modes: Structure to store mode information
%   ModeGrads: Structure to store mode gradient information
%   This is the $u.\nabla T$ term
    NModes=modInfo.POD.Config.NMODES;
    Term5   = zeros(NModes, NModes);
    T0 = Modes.T.Mean;

    % Iterate to calculate the constant for Term1
    for i = 1:NModes
        Bi = Modes.T.Modes(:,i);
        for j = 1:NModes
             Bj = Modes.T.Modes(:,j);
          
             if (Grid.dim == 3)
                dBj_dz = ModeGrads.T.Modes(:,j,3);
                dBj_dx = ModeGrads.T.Modes(:,j,1);
                dBj_dy = ModeGrads.T.Modes(:,j,2);     
             elseif (Grid.dim == 2)   
                dBj_dx = ModeGrads.T.Modes(:,j,1);
                dBj_dy = ModeGrads.T.Modes(:,j,2);
             else
                dBjx_dx=ModeGrads.T.Modes(:,j);
             end
           
            if(Grid.dim == 3)
               tmp = Bi.*(Bj.*dBj_dx + Bj.*dBj_dy + Bj.*dBj_dz);
               Term5(i,j) = -sum(tmp.*Grid.CellVolumes);
            elseif(Grid.dim == 2)
               tmp = Bi.*(Bj.*dBj_dx + Bj.*dBj_dy);
               [Term5(i,j)] = Integrate(Grid,tmp);
               Term5(i,j) = -Term5(i,j);
            else
               tmp = Bi.*(Bj.*dBjx_dx);
               Term5(i,j) = -sum(tmp.*Grid.CellVolumes);  
            end
           
        end
    end

end

