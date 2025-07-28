% ======================================================================================================================================
%%%%%%%%%% Preamble %%%%%%%%%%
% ======================================================================================================================================

% Part Ic - Free Particle Gaussian Wave Packet Simulation
% A single free particle Gaussian wave packet

% Author: Max L Butterworth
% MSc in Theoretical and Computational Chemistry Project
% University of Oxford

% ======================================================================================================================================
%%%%%%%%%% Define constants and variables %%%%%%%%%%
% ======================================================================================================================================

% Natural units adopted throughout
L = 50; % Length of the 1D box in m
m = 1; % Mass of electron in kg
h = 1; % Planck's constant in Js
hbar = 1; % Definition of h bar
N_steps = 1000; % Number of discretisation points

x0 = L/2; % Start evolving the wave packet from the centre of the box at t = 0
sigma = L/50; % Set the initial width of the wave packet

wp_energy = 20; % Set the wave packet energy
k = sqrt((wp_energy * 2 * m)/(hbar^2)); % Calculate the wavenumber from wp_energy; k = 0 gives a stationary Gaussian wave packet

% ======================================================================================================================================
%%%%%%%%%% Discretise the spatial domain, x, and time domain, t %%%%%%%%%%
% ======================================================================================================================================

x = linspace(0, L, N_steps); % Define the domain of the infinite potential well
dx = x(2) - x(1); % Calculate the spatial step size

dt = 1e-2; % Define the time step size
N_t = 1000; % Define the number of time steps to simulate

% ======================================================================================================================================
%%%%%%%%%% Construct the Hamiltonian using the finite difference method %%%%%%%%%%
% ======================================================================================================================================

% Construct the Hamiltonian inside the infinite potential well
laplacian = (1/dx^2) * spdiags([1, -2, 1], -1:1, N_steps, N_steps); % Define the Laplacian operator
H = -((hbar^2)/(2*m)) * laplacian; % Define the Hamiltonian operator

% ======================================================================================================================================
%%%%%%%%%% Generate an initial wave packet composed of one plane wave modulated by a Gaussian %%%%%%%%%%
% ======================================================================================================================================

psi0 = exp(-(x - x0).^2/(2 * sigma^2)) .* exp(1i * k * x); % Define the initial Gaussian wave packet
psi0_norm = psi0/sqrt(trapz(x, abs(psi0).^2)); % Normalise the initial Gaussian wave packet

% ======================================================================================================================================
%%%%%%%%%% Implement the Crank-Nicolson method to evolve the wavefunction and calculate probability current %%%%%%%%%%
% ======================================================================================================================================

J = zeros(N_steps, N_t); % Initialise an array to store probability currents
first_deriv = spdiags([-1, 1], 0:1, N_steps, N_steps);

psi = psi0_norm(1:N_steps).'; % Set the initial value of the wavefunction
psi_t = zeros(N_steps, N_t); % Initialise an array to store the wavefunction as it evolves in time
psi_t(:, 1) = psi; % Store the initial wavefunction in the time evolution array

% Pre-compute matrices required for the Crank-Nicolson method
A = eye(N_steps) + (((1i * dt)/(2 * hbar)) * H);
B = eye(N_steps) - (((1i * dt)/(2 * hbar)) * H);

for t = 2:N_t % Loop over all time steps
    psi = A \ (B * psi); % Evolve the wavefunction over time
    psi = psi/sqrt(trapz(x, abs(psi).^2)); % Normalise the time-evolved wavefunction
    psi_t(:, t) = psi; % Store the time-evolved wavefunction in the time evolution array
    J(:, t) = -((1i * hbar)/(2 * m)) * ((conj(psi) .* (first_deriv * psi)) - ((first_deriv * conj(psi)) .* psi)); % Calculate probability current at each point along x
end

% ======================================================================================================================================
%%%%%%%%%% Plot the time evolution of the wave packet probability density %%%%%%%%%%
% ======================================================================================================================================

figure; % Generate a figure

subplot(2, 2, 1) % Top left subfigure
real_wavefunction = plot(x, real(psi_t(:, 1))); % Plot the real wavefunction
xlabel('$x$', 'Interpreter','latex'); % Label the x-axis
ylabel('$\mathrm{Re}(\psi(x, t))$', 'Interpreter','latex'); % Label the y-axis
xlim([min(x) max(x)]) % Set the y-limits for convenience
ylim([min(real(psi_t(:))) max(real(psi_t(:)))]); % Set the y-limits foe convenience
title('Real Component of the Wavefunction') % Add a title
grid on; % Add a grid to the plot

subplot(2, 2, 2) % Top right subfigure
imag_wavefunction = plot(x, imag(psi_t(:, 1))); % Plot the imaginary wavefunction
xlabel('$x$', 'Interpreter','latex'); % Label the x-axis
ylabel('$\mathrm{Im}(\psi(x, t))$', 'Interpreter','latex'); % Label the y-axis
xlim([min(x) max(x)]) % Set the y-limits for convenience
ylim([min(imag(psi_t(:))) max(imag(psi_t(:)))]); % Set the y-limits for convenience
title('Imaginary Component of the Wavefunction') % Add a title
grid on; % Add a grid to the plot

subplot(2, 2, 3) % Bottom left subfigure
prob_density = plot(x, abs(psi_t(:, 1)).^2); % Plot the initial probability density
xlabel('$x$', 'Interpreter','latex'); % Label the x-axis
ylabel('$|\psi(x, t)|^2$', 'Interpreter','latex'); % Label the y-axis
xlim([min(x) max(x)]) % Set the y-limits for convenience
ylim([min(abs(psi_t(:)).^2) max(real(abs(psi_t(:)).^2))]); % Set the y-limits foe convenience
title('Probability Density') % Add a title
grid on; % Add a grid to the plot

subplot(2, 2, 4) % Bottom right subfigure
flux_plot = plot(x, J(:, 1)); % Plot the initial probability current
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
    set(flux_plot, 'YData', J(:, n)); % Update the probability density
    pause(0.1); % Pause to create an animation effect
    drawnow; % Update the relevant figures
end