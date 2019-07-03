function Term1 = ComputeTerm1( File,Grid, Modes, ModeGrads, modInfo)
%ComputeTerm1: Compute the first term in Yi Wang's derivation: -u0(dot)del(u0) 
%   Grid: Structure to store grid information
%   Modes: Structure to store mode information
%   ModeGrads: Structure to store mode gradient information

%   This is the $u_0.\nabla T_0$ term
    NModes=modInfo.POD.Config.NMODES;
    Term1   = zeros(NModes, 1);



    % Add z-compnent as needed
    if (Grid.dim == 3)
       T0 = Modes.T.Mean;
       u0x = Modes.U.Mean;
       u0y = Modes.V.Mean;
       u0z = Modes.W.Mean;
       dT0_dx = ModeGrads.T.Mean(:,1);
       dT0_dy = ModeGrads.T.Mean(:,2);
       dT0_dz = ModeGrads.T.Mean(:,3);
    elseif (Grid.dim == 2)
        T0 = Modes.T.Mean;
        u0x = Modes.U.Mean;
        u0y = Modes.V.Mean;
        dT0_dx = ModeGrads.T.Mean(:,1);
        dT0_dy = ModeGrads.T.Mean(:,2);
    else
        T0 = Modes.T.Mean;
        u0x = Modes.U.Mean;
        dT0_dx = ModeGrads.T.Mean;
    end

    % Iterate to calculate the constant for Term1
    for i = 1:NModes
        if(Grid.dim == 3)
            Bi = Modes.T.Modes(:,i);
            tmp = Bi.*(u0x.*dT0_dx+u0y.*dT0_dy+u0z.*dT0_dz);
            Term1(i) = -sum(tmp.*Grid.CellVolumes);
        elseif (Grid.dim == 2)
            Bi = Modes.T.Modes(:,i);
            tmp = Bi.*(u0x.*dT0_dx+u0y.*dT0_dy);
            [Term1(i)] = Integrate(Grid,tmp);
             Term1(i)=-Term1(i);
        else 
            Bi = Modes.T.Modes(:,i);
            tmp = Bi.*(u0x.*dT0_dx);
            Term1(i) = -sum(tmp.*Grid.CellVolumes);
        end
    end

end

