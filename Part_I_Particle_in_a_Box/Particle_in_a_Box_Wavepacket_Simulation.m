%%%%% Preamble %%%%%
% Part I - Particle in a Box Wavepacket Simulation
% Author: Max L Butterworth
% MSc in Theoretical and Computational Chemistry Project
% University of Oxford

%%%%%%%%%%

%%%%% Define constants %%%%%

L = 1e-10; % Length of the 1D box in m
m = 9.11e-31; % Mass of electron in kg
h = 6.626e-34; % Planck's constant in Js
hbar = h/(2*pi); % Definition of h bar
N_steps = 1000; % Number of discretisation points
N_superposition = 100; % The number of basis functions in the wavepacket superposition

%%%%% Discretise the spatial domain, x, and time domain, t %%%%%

x = linspace(0, L, N_steps); % Define the domain of the infinite potential well
dx = x(2) - x(1); % Calculate the spatial step size

dt = 1e-18; % Define the time step size
N_t = 1000; % Define the number of time steps to simulatie

%%%%% Solve the Schrödinger equation using the finite difference method %%%%%

% Construct the Hamiltonian inside the infinite potential well
laplacian = (1/dx^2) * spdiags([1, -2, 1], -1:1, N, N); % Define the Laplacian operator
H = -((hbar^2)/(2*m)) * laplacian; % Define the Hamiltonian operator

% Find the eigenvalues and eigenvectors of the Hamiltonian matrix
[psi, E] = eigs(H, N_steps, 'smallestabs');

% Normalise eigenvectors
psi_norm = psi/sqrt(trapz(x, abs(psi).^2));

%%%%% Generate a superposition of eigenstates %%%%%

% Initialise the wavepacket
wavepacket = zeros(N_steps, N_steps);

% Determine random coefficients for the basis functions
coeff = rand(N_superposition, 1);
coeff_norm = coeff/norm(coeff);

for i = 1:N_superposition
    wavepacket(:, i) = coeff_norm(i) * psi(:, i);
end

%%%%% Plot the first four eigenfunctions %%%%%
figure % Generate a figure

subplot(2, 2, 1) % Top left plot
scatter(x, psi(:,1)) % Plot the ground state wavefunction
title('Ground State Eigenfunction')

subplot(2, 2, 2) % Top right plot
scatter(x, psi(:,2)) % Plot the first excited state wavefunction
title('First Excited State Eigenfunction')

subplot(2, 2, 3) % Bottom left plot
scatter(x, psi(:,3)) % Plot the second excited state wavefunction
title('Second Excited State Eigenfunction')

subplot(2, 2, 4) % Bottom right plot
scatter(x, psi(:,4)) % Plot the third excited state wavefunction
title('Third Excited State Eigenfunction')