%--------------------------------------------------------------------------
%
% Quantum dynamics of coupled excitons and phonons on a homogeneous 
% cyclic chain, subject to the following simplifying assumptions:
% - Only two-state models: electronic ground and excited states
% - Harmonic approximation for phonons: on-site and NN-pair
% - Linear (sigma) coupling of the two subsystems
%
% see our work at DOI:10.1063/5.0074948 (TISE)
%
% N : Number of interacting sites, each with electron and vibrational dof's
%
%--------------------------------------------------------------------------
function qm_init (N)
global hamilt plots space state time

% Parameters estimate, fallen from heaven
alpha = + 0.1000;   % on-site energy for excitons
beta  = - 0.0100;   % NN-pair coupling for excitons
eta   = + 0.0000;   % constant energy offset
mass  = + 1.0000;   % particle mass
nu    = + 0.0010;   % on-site harmonic vibration
omega = sqrt(2)*nu; % NN-pair harmonic vibration
chi   = + 0.0000;   % exciton-phonon tuning, localized
rho   = + 0.0000;   % exciton-phonon tuning, non-symmetric
sig   = + 0.0002;   % exciton-phonon tuning, symmetrized
tau   = + 0.0000;   % exciton-phonon coupling

% Log file output
prt.disp ( '***************************************************************' )
prt.disp ( 'Fröhlich-Holstein-Peierls model for excitons coupled to phonons' )
prt.disp ( 'using direct (local mode) representation of excitons/phonons   ' )
prt.disp ( '***************************************************************' )
prt.disp ( ' ' )
prt.disp (['Electronic site energy         : ' num2str(alpha )]) 
prt.disp (['Electronic coupling strength   : ' num2str(beta  )]) 
prt.disp (['Constant energy offset         : ' num2str(eta   )]) 
prt.disp (['Particle mass                  : ' num2str(mass  )])
prt.disp (['Site frequency (for restraint) : ' num2str(nu    )])
prt.disp (['Pair frequency (for phonons)   : ' num2str(omega )]) 
prt.disp (['Electron-phonon tuning   (chi) : ' num2str(chi   )]) 
prt.disp (['Electron-phonon tuning   (rho) : ' num2str(rho   )]) 
prt.disp (['Electron-phonon tuning   (sig) : ' num2str(sig   )]) 
prt.disp (['Electron-phonon coupling (tau) : ' num2str(tau   )]) 
prt.disp (' ')

% Toggle saving wavefunctions to files
state.save_export=false; 

% Number of sites; electronic coupling scheme
hamilt.coupling.n_eqs      = N;          % Number of (coupled) Schrödinger equations
hamilt.coupling.represent  = 'dia';      % Diabatic or adiabatic representation

% Equal spatial discretization for cyclic chains
for d = 1:N
    space.dof{d}       = dof.fft;        % Using FFT methods
    space.dof{d}.mass  = mass;           % Effective mass
    space.dof{d}.x_min = -120;           % Lower bound
    space.dof{d}.x_max = +120;           % Upper bound
    space.dof{d}.n_pts = 08;             % Number of grid points / basis functions
%     space.dof{d}       = dof.hermite;    % Using Gauss-Hermite DVR/FBR
%     space.dof{d}.mass  = mass;           % Effective mass
%     space.dof{d}.omega = 2*nu;           % Harmonic frequency
%     space.dof{d}.n_pts = 10;             % Number of grid points / basis functions
end


% Set up diabatic potential energy matrix for coupled electrons vs. phonons
for m = 1:N
    for n = m:N
        if n==m % Diagonal
            hamilt.pot{m,n} = pot.chain; % Chain of harmonic oscillators
            hamilt.pot{m,n}.oPR = nu;    % harmonic frequency: position restraints
            hamilt.pot{m,n}.oNN = omega; % harmonic frequency: nearest neighbors
            hamilt.pot{m,n}.oCM = 0.0;   % harmonic frequency: c-of-m restraint
            hamilt.pot{m,n}.r_e = 0.0;   % equilibrium distance
            hamilt.pot{m,n}.pbc = true;  % Periodic boundary conditions
            hamilt.pot{m,n}.off = alpha+eta; % Vertical energy offset
            hamilt.pot{m,n}.chi = chi;   % Electron-phonon tuning, localized
            hamilt.pot{m,n}.rho = rho;   % Electron-phonon tuning, non-symmetric
            hamilt.pot{m,n}.sig = sig;   % Electron-phonon tuning, symmetrized
        elseif n==m+1 % Off-diagonal: NN coupling
            if N==2 % including periodic boundary conditions into super-diagonal (1,2)
                hamilt.pot{m,n} = pot.chain;  % Chain of harmonic oscillators
                hamilt.pot{m,n}.off = 2*beta; % Vertical energy offset
                hamilt.pot{m,n}.tau = 2*tau;  % Electron-phonon coupling
                hamilt.pot{m,n}.pbc = true;   % Periodic boundary conditions
             else 
                hamilt.pot{m,n} = pot.chain;
                hamilt.pot{m,n}.off = beta;   % Vertical energy offset
                hamilt.pot{m,n}.tau = tau;    % Electron-phonon coupling
                hamilt.pot{m,n}.pbc = true;   % Periodic boundary conditions
             end
        end
    end
end

% Periodic boundary conditions
if N>2
    hamilt.pot{1,N} = pot.chain;         % Chain of harmonic oscillators
    hamilt.pot{1,N}.off = beta;          % Vertical energy offset
    hamilt.pot{1,N}.tau = tau;           % Electron-phonon coupling
    hamilt.pot{1,N}.pbc = true;          % Periodic boundary conditions
end

% QM_PROPA only: Discretization of imaginary(!) time
time.steps.m_start  = 000;               % Index of initial time step
time.steps.m_stop   = 025;               % Index of final time step
time.steps.m_delta  = 100;               % Size of time steps

% QM_PROPA only: Initial wave function for relaxation
hamilt.coupling.ini_rep    = 'dia';      % Initial WF in dia or adi representation
hamilt.coupling.ini_coeffs = ones(1,N);  % Initially equal on ALL diabatic states
for d = 1:N
    time.dof{d}       = init.gauss;      % Gaussian-shaped wavepacket
    time.dof{d}.pos_0 =  0.0;            % Center in position representation
    time.dof{d}.mom_0 =  0.0;            % Center in momentum representation
    time.dof{d}.width =  20.0;           % Widths of wavepacket
end
    
% Visualize time evolution of densities
plots.density = vis.reduced_1d;  % Reduced densities in 1D
plots.density.represent = 'dvr'; % R representation

% Plot expectation values versus time
plots.expect = vis.expect;
plots.expect.e_min = 0.080;              % Minimum of energy in expct. plot
plots.expect.e_max = 0.085;              % Maximum of energy in expct. plot
% plots.expect.p_min = -0.1;
% plots.expect.p_max = +1.1;
