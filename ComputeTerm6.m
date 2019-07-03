function Term6 = ComputeTerm6( File,Grid, Modes, ModeGrads, modInfo )
%ComputeTerm2: Compute the second term in Yi Wang's derivation: nu*Laplacian(u0) 
%   Grid: Structure to store grid information
%   Modes: Structure to store mode information
%   ModeGrads: Structure to store mode gradient information

%   Note that in the current version we didn't consider the boundary
%   integral, which will be inclued in the future version
%   This is the $D \delta (D.\delta T)$ term
    nu      = modInfo.Phys.nu;
    NModes=modInfo.POD.Config.NMODES;
    Term6   = zeros(NModes, NModes);




    % Iterate to calculate the constant for Term1
    for i = 1:NModes
        dBix_dx=ModeGrads.T.Modes(:,i,1);
        %dBix_dx=dBix_dx./dBix_dx;
        for j = 1:NModes
            if(Grid.dim == 3)
                dBi_dx = ModeGrads.T.Modes(:,i,1);
                dBi_dy = ModeGrads.T.Modes(:,i,2);
                dBi_dz = ModeGrads.T.Modes(:,i,3);
                tmp = dBix_dx.*dBi_dx +dBix_dy.*dBix_dy+dBix_dz.*dBix_dz;
                Term6(i) = -nu*sum(tmp.*Grid.CellVolumes);
            elseif(Grid.dim == 2)
                dBi_dx=ModeGrads.T.Modes(:,i,1);
                dBi_dy=ModeGrads.T.Modes(:,i,2);
                dBj_dx = ModeGrads.T.Modes(:,j,1);
                dBj_dy = ModeGrads.T.Modes(:,j,2);
                Bj = Modes.T.Modes(:,j);
                tmp = dBi_dx.*dBj_dx+dBi_dy.*dBj_dy;
                BoundSet=Grid.Boundry;
                tmp(BoundSet)=tmp(BoundSet)-...
                    Bj(BoundSet).*dBj_dx(BoundSet,:)+Bj(BoundSet).*dBj_dx(BoundSet,:); 
                [Term6(i)] = Integrate(Grid,tmp);
                Term6(i) = (modInfo.Phys.nu)*Term6(i);
            else
                dBjx_dx=ModeGrads.T.Modes(:,j);
                %dBjx_dx=dBjx_dx./dBjx_dx;
                Bj = Modes.T.Modes(:,j);
                %Bj = Modes.T.Modes(j,:);
                tmp = dBix_dx.*dBjx_dx;
                BoundSet=[1,Grid.nnode];
                tmp(BoundSet)=tmp(BoundSet)-Bj(BoundSet).*dBjx_dx(BoundSet,:);    
                Term6(i,j) = (modInfo.Phys.nu)*sum(tmp.*Grid.CellVolumes);  
            end
        end 
    end

end

