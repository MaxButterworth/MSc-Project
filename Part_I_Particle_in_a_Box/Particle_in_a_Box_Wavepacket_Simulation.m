%%%%%%%%%% Preamble %%%%%%%%%%
% Part I - Particle in a Box Wavepacket Simulation
% Author: Max L Butterworth
% MSc in Theoretical and Computational Chemistry Project
% University of Oxford

%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Define constants %%%%%%%%%%

L = 1e-10; % Length of the 1D box in m
m = 9.11e-31; % Mass of electron in kg
h = 6.626e-34; % Planck's constant in Js
hbar = h/(2*pi); % Definition of h bar
N_steps = 1000; % Number of discretisation points
N_superposition = 100; % The number of basis functions in the wavepacket superposition

%%%%%%%%%% Discretise the spatial domain, x, and time domain, t %%%%%%%%%%

x = linspace(0, L, N_steps); % Define the domain of the infinite potential well
dx = x(2) - x(1); % Calculate the spatial step size

dt = 1e-18; % Define the time step size
N_t = 1000; % Define the number of time steps to simulatie

%%%%%%%%%% Solve the Schrödinger equation using the finite difference method %%%%%%%%%%

% Construct the Hamiltonian inside the infinite potential well
laplacian = (1/dx^2) * spdiags([1, -2, 1], -1:1, N, N); % Define the Laplacian operator
H = -((hbar^2)/(2*m)) * laplacian; % Define the Hamiltonian operator
H = H(2:N-1, 2:N-1); % Impose boundary conditions: psi(0) = psi(L) = 0

% Find the eigenvalues and eigenvectors of the Hamiltonian matrix
[psi, E] = eigs(H, N_eigenvectors, 'smallestabs');

% Normalise eigenvectors
x_internal = x(2:N-1); % Exclude the boundaries from the normalisation process

for j = 1:N_eigenvectors
    psi_norm(:, j) = psi(:, j)/sqrt(trapz(x_internal, abs(psi(:, j)).^2));
end

%%%%%%%%%% Generate an initial wavepacket %%%%%%%%%%

x0 = L/2; % Start evolving the wavepacket from the centre of the box at t = 0
sigma = L/10; % Set the initial width of the wavepacket
psi0 = exp(-(x - x0).^2/(2 * sigma^2)); % Define the initial Gaussian wavepacket
psi0_norm = psi0/sqrt(trapz(x, abs(psi0).^2)); % Normalise the initial Gaussian wavepacket

%%%%%%%%%% Plot the first four eigenfunctions %%%%%%%%%%

figure % Generate a figure

subplot(2, 2, 1) % Top left plot
scatter(x_internal, psi_norm(:,1)) % Plot the ground state wavefunction
title('Ground State Eigenfunction')

subplot(2, 2, 2) % Top right plot
scatter(x_internal, psi_norm(:,2)) % Plot the first excited state wavefunction
title('First Excited State Eigenfunction')

subplot(2, 2, 3) % Bottom left plot
scatter(x_internal, psi_norm(:,3)) % Plot the second excited state wavefunction
title('Second Excited State Eigenfunction')

subplot(2, 2, 4) % Bottom right plot
scatter(x_internal, psi_norm(:,4)) % Plot the third excited state wavefunction
title('Third Excited State Eigenfunction')