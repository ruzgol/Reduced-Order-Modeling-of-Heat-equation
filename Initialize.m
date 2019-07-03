function [ Grid, Modes, ModeGrads, ODEs ] = Initialize( File )
%Initialize Initialize the entire simulation
%   File: The file used to initilize the POD-ROM for CFD settings and
%   simulation DATA
%   Grid: The structure to store the grid information
%   Modes:  POD modes of the simulation including pressure, velocity, etc.
%   Grid:   Grid information
%   ModesGrad:  Gradient of modes
%   ODEs: Structure to store the ODE Set
    
    % Read ODEs: if found exit the program immidiately
    ODEs = ReadODEs(File);
    
    % Read Grid information
    Grid = ReadGrid(File);
    
    % Read POD modes and POD gradients
    [Modes, ModeGrads] = ReadModesAndModeGrads(File);


end

