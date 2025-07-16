%%%%%%%%%% Preamble %%%%%%%%%%
% Part I - Particle in a Box Wave Packet Simulation
% Superposition of particle in a box eigenstates modulated by a Gaussian
% Time Propagation conducted with the split operator method

% Author: Max L Butterworth
% MSc in Theoretical and Computational Chemistry Project
% University of Oxford

%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Define constants and variables %%%%%%%%%%

% Natural units have been used throughout
L = 1;% Length of the 1D box
m = 1; % Mass
h = 1; % Planck's constant in J s
hbar = 1; % Definition of h bar

N_steps = 1000; % Number of discretisation points on the x-axis

basis_funcs_indices = [1, 2, 3]; % Create an array of the indices of PIB_eigenstates_norm that form the superposition
basis_funcs_coeffs = rand(1, length(basis_funcs_indices)); % Weightings of PIB eigenstates in the superposition
N_PIB_eigenfuncs = max(basis_funcs_indices); % The number of basis functions in the wave packet superposition

%%%%%%%%%% Discretise the spatial domain, x, and time domain, t %%%%%%%%%%

x = linspace(0, L, N_steps); % Define the domain of the infinite potential well
dx = x(2) - x(1); % Calculate the spatial step size

dt = 1e-3; % Define the time step size
N_t = 1000; % Define the number of time steps to simulate

%%%%%%%%% Solve the Schrödinger equation using the finite difference method %%%%%%%%%%

% Construct the Hamiltonian inside the infinite potential well
laplacian = (1/dx^2) * spdiags([1, -2, 1], -1:1, N_steps, N_steps); % Define the Laplacian operator
H = -((hbar^2)/(2*m)) * laplacian; % Define the Hamiltonian operator

% Find the eigenvalues and eigenvectors of the Hamiltonian matrix
[PIB_eigenstates, E_PIB] = eigs(H, N_PIB_eigenfuncs, 'smallestabs');

PIB_eigenstates_norm = zeros(N_steps, N_PIB_eigenfuncs); % Set up an array to store normalised PIB eigenfunctions

% Normalise eigenvectors
for m = 1:N_PIB_eigenfuncs
    PIB_eigenstates_norm(:, m) = PIB_eigenstates(:, m)/sqrt(trapz(x, abs(PIB_eigenstates(:, m)).^2));
end

%%%%%%%%%% Generate an initial wave packet composed of a superposition of PIB eigenfunctions modulated by a Gaussian %%%%%%%%%%

x0 = L/2; % Start evolving the wave packet from the centre of the box at t = 0
sigma = L/20; % Set the initial width of the wave packet

psi0 = zeros(N_steps, 1); % Initialise an empty array to store the initial wave packet

% Generate the superposition of PIB basis functions
for l = 1:length(basis_funcs_indices)
    psi0 = psi0 + (basis_funcs_coeffs(l) * PIB_eigenstates_norm(:, basis_funcs_indices(l)));
end

psi0 = exp(-(x - x0).^2/(2 * sigma^2)) .* psi0; % Modulate the superposition by a Gaussian

psi0_norm = psi0/sqrt(trapz(x, abs(psi0).^2)); % Normalise the initial Gaussian wave packet

%%%%%%%%%% Impose boundary conditions %%%%%%%%%%

% Set the wavefunction to zero at the boundaries
psi0_norm(1) = 0;
psi0_norm(N_steps) = 0;

x_internal = x(2:N_steps - 1); % Truncate the x array to account for boundary conditions

%%%%%%%%%% Propagate the wave packet through time using the split operator method %%%%%%%%%%

% Define kinetic and potential operators
p = (2 * pi * (0:N_steps-1)/L); % Define the momentum space grid
V = zeros(N_steps, 1); % Define the potential energy array (zero for all 0 < x < L for particle in a box)

T_op = exp(-(1i * (p.^2) * dt)/(2 * m * hbar)).'; % Kinetic energy operator (full time step)
V_op = exp(-(1i * V * dt)/(2 * hbar)); % Potential energy operator (half time step)

% Initialise arrays to store fluxes and a first derivativ operator to calculate the fluxes
J = zeros(N_steps, N_t); % Initialise an array to store probability currents
first_deriv = spdiags([-1, 1], 0:1, N_steps, N_steps)/dx; % Define a first derivative operator using the finite difference method

% Set up the wave packet on the internal coordinates for propagation
psi = psi0_norm(:); % Set the initial value of the wavefunction
psi_t = zeros(N_steps, N_t); % Initialise an array to store the wavefunction as it evolves in time
psi_t(:, 1) = psi; % Store the initial wavefunction in the time evolution array

% Propagate the wave packet defined on the internal coordinates
for t = 2:N_t
    psi = V_op .* psi; % Operate a half time step in real space
    psi_k = fft(psi); % Fourier transform the wavefunction into k-space
    psi_k = T_op .* psi_k; % Operate a full time step in k-space
    psi = ifft(psi_k); % Inverse Fourier transform into real space
    psi = V_op .* psi; % Operate a half time step in real space

    psi(1, 1) = 0; % Impose boundary condition at x = 0
    psi(N_steps, 1) = 0; % Impose boundary condition at x = L
    psi = psi/sqrt(trapz(x, abs(psi).^2)); % Normalise the time-evolved wavefunction

    psi_t(:, t) = psi; % Store the time-evolved wavefunction in the time evolution array
    J(:, t) = -((1i * hbar)/(2 * m)) * ((conj(psi) .* (first_deriv * psi)) - ((first_deriv * conj(psi)) .* psi)); % Calculate probability current at each point along x
end

%%%%%%%%%% Plot the time evolution of the wave packet, probability density, and flux %%%%%%%%%%

figure; % Generate a figure

subplot(2, 2, 1) % Top Left subfigure
real_wavefunction = plot(x, real(psi_t(:, 1))); % Plot the real wavefunction
xlabel('$x$', 'Interpreter', 'latex'); % Label the x-axis
ylabel('$\mathrm{Re}(\psi(x, t))$', 'Interpreter', 'latex'); % Label the y-axis
ylim([min(real(psi_t(:))) max(real(psi_t(:)))]); % Set the y-limits for convenience
title('Real Component of the Wavefunction') % Add a title
grid on; % Add a grid to the plot

subplot(2, 2, 2) % Top right subfigure
imag_wavefunction = plot(x, imag(psi_t(:, 1))); % Plot the imaginary wavefunction
xlabel('$x$', 'Interpreter', 'latex'); % Label the x-axis
ylabel('$\mathrm{Im}(\psi(x, t))$', 'Interpreter', 'latex'); % Label the y-axis
ylim([min(imag(psi_t(:))) max(imag(psi_t(:)))]); % Set the y-limits for convenience
title('Imaginary Component of the Wavefunction') % Add a title
grid on; % Add a grid to the plot

subplot(2, 2, 3) % Bottom Right subfigure
prob_density = plot(x, abs(psi_t(:, 1)).^2); % Plot the initial probability density
xlabel('$x$', 'Interpreter', 'latex'); % Label the x-axis
ylabel('$|\psi(x, t)|^2$', 'Interpreter', 'latex'); % Label the y-axis
ylim([min(abs(psi_t(:)).^2) max(real(abs(psi_t(:)).^2))]); % Set the y-limits for convenience
title('Probability Density') % Add a title
grid on; % Add a grid to the plot

subplot(2, 2, 4) % Bottom left subfigure
flux_plot = plot(x, J(:, 1)); % Plot the initial probability density
xlabel('$x$', 'Interpreter', 'latex'); % Label the x-axis
ylabel('$J(x, t)$', 'Interpreter', 'latex'); % Label the y-axis
ylim([min(J(:)) max(J(:))]); % Set the y-limits for convenience
title('Probability Current') % Add a title
grid on; % Add a grid to the plot

% Animate the figures
for n = 1:N_t % Loop over all timesteps
    set(real_wavefunction, 'YData', real(psi_t(:, n))) % Update the real part of the wavefunction
    set(imag_wavefunction, 'YData', imag(psi_t(:, n))) % Update the imaginary part of the wavefunction
    set(prob_density, 'YData', abs(psi_t(:, n)).^2); % Update the probability density
    set(flux_plot, 'YData', J(:, n)); % Update the flux plot
    pause(0.05); % Pause to create an animation
    drawnow; % Update the figures and display immediately
end