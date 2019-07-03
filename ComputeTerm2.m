function Term2 = ComputeTerm2( File,Grid, Modes, ModeGrads, modInfo )
%ComputeTerm2: Compute the second term in Yi Wang's derivation: nu*Laplacian(u0) 
%   Grid: Structure to store grid information
%   Modes: Structure to store mode information
%   ModeGrads: Structure to store mode gradient information

%   Note that in the current version we didn't consider the boundary
%   integral, which will be inclued in the future version
%   This is the $ \delta (D.\delta T_0)$ term
    nu      = modInfo.Phys.nu;
    NModes=modInfo.POD.Config.NMODES;
    Term2   = zeros(NModes, 1);



    % Add z-compnent as needed
    if (Grid.dim == 3)
       dT0_dx = ModeGrads.T.Mean(:,1);
       dT0_dy = ModeGrads.T.Mean(:,2);
       dT0_dz = ModeGrads.T.Mean(:,3);
    elseif (Grid.dim == 2)  
       dT0_dx = ModeGrads.T.Mean(:,1);
       dT0_dy = ModeGrads.T.Mean(:,2); 
    else
       dT0_dx = ModeGrads.T.Mean; 
    end

    % Iterate to calculate the constant for Term1
    for i = 1:NModes
        if(Grid.dim == 3)
            dBi_dx = ModeGrads.T.Modes(:,i,1);
            dBi_dy = ModeGrads.T.Modes(:,i,2);
            dBi_dz = ModeGrads.T.Modes(:,i,3);
            tmp = dBix_dx.*dT0x_dx+dBix_dy.*dT0x_dy+dBix_dz.*dT0x_dz;
            Term2(i) = -nu*sum(tmp.*Grid.CellVolumes);
        elseif(Grid.dim == 2) 
            Bi = Modes.T.Modes(:,i);
            dBi_dx = ModeGrads.T.Modes(:,i,1);
            dBi_dy = ModeGrads.T.Modes(:,i,2);
            tmp = dBi_dx.*dT0_dx+dBi_dy.*dT0_dy;
            BoundSet=Grid.Boundry;
            tmp(BoundSet)=tmp(BoundSet)-Bi(BoundSet).*dT0_dx(BoundSet)...
                          -Bi(BoundSet).*dT0_dy(BoundSet);  
            [Term2(i)] = Integrate(Grid,tmp);
            Term2(i) = (modInfo.Phys.nu)*Term2(2);
        else
            dBi_dx=ModeGrads.T.Modes(:,i);
            %dBi_dx=dBi_dx./dBi_dx;
            %dBi_dx=ModeGrads.T.Modes(i,:);
            Bi = Modes.T.Modes(:,i);
            %Bi = Modes.T.Modes(i,:);
            %tmp = dBi_dx.*dT0_dx;
            tmp = dBi_dx.*dBi_dx;
            BoundSet=[1,Grid.nnode];
            %tmp(BoundSet)=tmp(BoundSet)-Bi(BoundSet).*dT0_dx(BoundSet);
            %tmp(BoundSet)=tmp(BoundSet)-Bi(BoundSet).*dBi_dx(BoundSet);
            Term2(i) = (modInfo.Phys.nu)*sum(tmp.*Grid.CellVolumes);   

        end
    end

end

