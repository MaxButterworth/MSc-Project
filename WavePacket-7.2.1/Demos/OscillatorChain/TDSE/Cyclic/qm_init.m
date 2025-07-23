% Copyright (C) 2023-.... Burkhard Schmidt's group
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init (N)
global hamilt plots space time

prt.disp ( '*****************************************************' )
prt.disp ( 'Cyclic chain of harmonic oscillators in one dimension' )
prt.disp ( ['Number of sites : ' int2str(N)])
prt.disp ( '*****************************************************' )

% Physical parameters
mass  = + 1.0000;                        % Particle masses
nu    = + 0.0010;                        % On-site harmonic vibration
omega = sqrt(2)*nu;                      % NN-pair harmonic vibration
nu_E = sqrt(nu^2+omega^2);               % Effective frequencies

% Spatial discretization
for d=1:N
    space.dof{d}       = dof.fft;        % FFT-based DVR/FBR
    space.dof{d}.n_pts = 12;             % Number of grid points 
    space.dof{d}.x_min = -100;           % Lower bound
    space.dof{d}.x_max = +100;           % Upper bound
    space.dof{d}.mass  = mass;           % Mass for the kinetic energy
end

% Initial wave functions
for d=1:N
    time.dof{d}       = init.gauss;      % Gaussian-shaped wavepacket
    if d == round((N+1)/2)
        time.dof{d}.pos_0   = 20;        % Set initial displacement
    else
        time.dof{d}.pos_0   = 00;        % No initial displacement
    end
    time.dof{d}.mom_0 = 0.0;             % Center in momentum representation
    time.dof{d}.width = (2*mass*nu_E)^(-1/2);
end

% Temporal discretization
time.steps.m_start  = 000;               % Index of initial time step
time.steps.m_stop   = 050;               % Index of final time step
time.steps.m_delta  = 200;               % Size of time steps
time.steps.s_number = 050;               % Number of sub steps per time step

hamilt.pot{1,1}     = pot.chain;         % Chain of harmonic oscillators
hamilt.pot{1,1}.oNN = omega;             % Harmonic frequency: nearest neighbors
hamilt.pot{1,1}.oPR = nu;                % Harmonic frequency: position restraints
hamilt.pot{1,1}.pbc = true;              % Toggle periodic boundary conditions

% Plots of densities
plots.density = vis.reduced_1d;          % Reduced density plots
plots.density.represent = 'wig';         % (R,P) Wigner representation
plots.density.expect = true;             % 

% Plots of expectation values
plots.expect        = vis.expect;   
plots.expect.e_max  = 0.006;             % Tune the display of the energy values
