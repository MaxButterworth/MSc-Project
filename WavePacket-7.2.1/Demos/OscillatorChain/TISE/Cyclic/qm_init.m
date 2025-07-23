% Copyright (C) 2023-.... Burkhard Schmidt's group
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init (n)
global hamilt plots space time

prt.disp ( '*****************************************************' )
prt.disp ( 'Cyclic chain of harmonic oscillators in one dimension' )
prt.disp ( ['Number of sites : ' int2str(n)])
prt.disp ( '*****************************************************' )

% Spatial discretization
for j=1:n
    space.dof{j}       = dof.hermite;    % Gauss-Hermite DVR/FBR
    space.dof{j}.n_pts = 8;              % Number of grid points 
    space.dof{j}.omega = 1e-3 * sqrt(2); % Harmonic frequency
    space.dof{j}.mass  =  1;             % Mass for the kinetic energy
end

% Initial wave functions
for j=1:n
    time.dof{j}       = init.gauss;      % Gaussian-shaped wavepacket
    time.dof{j}.width =  2.0;            % Width
    time.dof{j}.pos_0 =  0.0;            % Center in position representation
    time.dof{j}.mom_0 =  0.0;            % Center in momentum representation
end

% Temporal discretization
time.steps.m_start  = 000;               % Index of initial time step
time.steps.m_stop   = 020;               % Index of final time step
time.steps.m_delta  = 200;               % Size of time steps 

hamilt.pot{1,1}     = pot.chain;         % Chain of harmonic oscillators
hamilt.pot{1,1}.oNN = 1e-3 * sqrt(2);    % Harmonic frequency: nearest neighbors
hamilt.pot{1,1}.oPR = 1e-3;              % Harmonic frequency: position restraints
hamilt.pot{1,1}.pbc = true;              % Toggle periodic boundary conditions

% Plots of densities
plots.density         = vis.reduced_2d;  % Reduced density plots

% Plots of expectation values
plots.expect        = vis.expect;   
plots.expect.e_max  = 0.02;              % Tune the display of the energy values
