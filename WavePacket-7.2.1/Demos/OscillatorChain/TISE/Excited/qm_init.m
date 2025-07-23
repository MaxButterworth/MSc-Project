% Copyright (C) 2023-.... Burkhard Schmidt's group
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init (n)
global hamilt plots space

prt.disp ( '*****************************************************' )
prt.disp ( 'Cyclic chain of harmonic oscillators in one dimension' )
prt.disp ( ['Number of sites : ' int2str(n)])
prt.disp ( '*****************************************************' )

% Spatial discretization
for d=1:n
    space.dof{d}       = dof.fft;        % Using FFT methods
    space.dof{d}.mass  = 1.0;            % Effective mass
    space.dof{d}.x_min = -120;           % Lower bound
    space.dof{d}.x_max = +120;           % Upper bound
    space.dof{d}.n_pts = 024;            % Number of grid points / basis functionsend
end

% Potential enregy function
hamilt.pot{1,1}     = pot.chain;         % Chain of harmonic oscillators
hamilt.pot{1,1}.oNN = 1e-3 * sqrt(2);    % Harmonic frequency: nearest neighbors
hamilt.pot{1,1}.oPR = 1e-3;              % Harmonic frequency: position restraints
hamilt.pot{1,1}.pbc = true;              % Toggle periodic boundary conditions

% Eigenfunctions
hamilt.eigen.start = 000;
hamilt.eigen.stop  = 008;                % Toroidal shape for cyclic trimer
hamilt.eigen.storage = 's';              % sparse storage scheme

% Density plots
plots.density          = vis.surface;
plots.density.srf_view = [115 20];        % Viewing angle [az el]

% Plots of expectation values
plots.expect        = vis.expect;   
plots.expect.e_max  = 0.007;              % Tune the display of the energy values
