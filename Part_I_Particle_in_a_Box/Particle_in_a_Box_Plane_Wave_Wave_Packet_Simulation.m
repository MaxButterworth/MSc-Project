%%%%%%%%%% Preamble %%%%%%%%%%
% Part I - Particle in a Box Wave Packet Simulation
% Superposition of a free particle eigenstates modulated by a Gaussian
% Author: Max L Butterworth
% MSc in Theoretical and Computational Chemistry Project
% University of Oxford

%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Define constants %%%%%%%%%%

L = 1e-10; % Length of the 1D box in m
m = 9.110e-31; % Mass of electron in kg
h = 6.626e-34; % Planck's constant in Js
hbar = h/(2*pi); % Definition of h bar
N_steps = 1000; % Number of discretisation points
N_superposition = 100; % The number of basis functions in the wave packet superposition

%%%%%%%%%% Discretise the spatial domain, x, and time domain, t %%%%%%%%%%

x = linspace(0, L, N_steps); % Define the domain of the infinite potential well
dx = x(2) - x(1); % Calculate the spatial step size

dt = 1e-20; % Define the time step size
N_t = 1000; % Define the number of time steps to simulate

%%%%%%%%%% Construct the Hamiltonian using the finite difference method %%%%%%%%%%

% Construct the Hamiltonian inside the infinite potential well
laplacian = (1/dx^2) * spdiags([1, -2, 1], -1:1, N_steps, N_steps); % Define the Laplacian operator
H = -((hbar^2)/(2*m)) * laplacian; % Define the Hamiltonian operator

%%%%%%%%%% Generate an initial wave packet composed of one plane wave modulated by a Gaussian %%%%%%%%%%

k = (50 * pi)/L; % Set the wavenumber
x0 = L/2; % Start evolving the wave packet from the centre of the box at t = 0
sigma = L/20; % Set the initial width of the wave packet
psi0 = exp(-(x - x0).^2/(2 * sigma^2)) .* exp(1i * k * x); % Define the initial Gaussian wave packet
psi0_norm = psi0/sqrt(trapz(x, abs(psi0).^2)); % Normalise the initial Gaussian wave packet

%%%%%%%%%% Impose boundary conditions %%%%%%%%%%

% Set the wavefunction to zero at the boundaries
psi0_norm(1) = 0;
psi0_norm(N_steps) = 0;
psi0_norm = sparse(psi0_norm); % Define the initial wavefunction as a sparse matrix to speed up the calculation

H = H(2:N_steps-1, 2:N_steps-1); % Impose boundary conditions: psi(0) = psi(L) = 0
H = sparse(H); % Define the Hamiltonian as a sparse matrix to speed up the calcualtion

x_internal = x(2:N_steps - 1); % Truncate the x array to account for boundary conditions

%%%%%%%%%% Implement the Crank-Nicolson method to evolve the wavefunction and calculate probability current %%%%%%%%%%

J = zeros(1, N_t); % Initialise an array to store probability currents
first_deriv = spdiags([-1, 1], 0:1, N_steps-2, N_steps-2);

psi = psi0_norm(2:N_steps-1).'; % Set the initial value of the wavefunction
psi_t = zeros(N_steps, N_t); % Initialise an array to store the wavefunction as it evolves in time
psi_t(2:N_steps-1, 1) = psi; % Store the initial wavefunction in the time evolution array

% Pre-compute matrices required for the Crank-Nicolson method
A = eye(N_steps-2) + (((1i * dt)/(2 * hbar)) * H);
B = eye(N_steps-2) - (((1i * dt)/(2 * hbar)) * H);

for t = 2:N_t % Loop over all time steps
    psi = A \ (B * psi); % Evolve the wavefunction over time
    psi = psi/sqrt(trapz(x_internal, abs(psi).^2)); % Normalise the time-evolved wavefunction
    psi_t(2:N_steps-1, t) = psi; % Store the time-evolved wavefunction in the time evolution array
    J(t) = -((1i * hbar)/(2 * m)) * ((psi' * (first_deriv * psi)) -  ((first_deriv * conj(psi))' * psi)); % Calculate probability current
end

%%%%%%%%%% Plot the time evolution of the wave packet probability density %%%%%%%%%%

x_ang = x * 1e10; % Generate an array of x-values in angstroms

figure; % Generate a figure

subplot(1, 3, 1) % Left subfigure
real_wavefunction = plot(x_ang, real(psi_t(:, 1))); % Plot the real wavefunction
xlabel('$x\ (\AA)$', 'Interpreter','latex'); % Label the x-axis
ylabel('$\mathrm{Re}(\psi(x, t))$', 'Interpreter','latex'); % Label the y-axis
title('Real Component of the Wavefunction') % Add a title
grid on; % Add a grid to the plot

subplot(1, 3, 2) % Middle subfigure
imag_wavefunction = plot(x_ang, imag(psi_t(:, 1))); % Plot the imaginary wavefunction
xlabel('$x\ (\AA)$', 'Interpreter','latex'); % Label the x-axis
ylabel('$\mathrm{Im}(\psi(x, t))$', 'Interpreter','latex'); % Label the y-axis
title('Imaginary Component of the Wavefunction') % Add a title
grid on; % Add a grid to the plot

subplot(1, 3, 3) % Right subfigure
prob_density = plot(x_ang, abs(psi_t(:, 1)).^2); % Plot the initial probability density
xlabel('$x\ (\AA)$', 'Interpreter','latex'); % Label the x-axis
ylabel('$|\psi(x, t)|^2$', 'Interpreter','latex'); % Label the y-axis
title('Probability Density') % Add a title
grid on; % Add a grid to the plot

% Animate the figures

for n = 1:N_t % Loop over all timesteps
    set(real_wavefunction, 'YData', real(psi_t(:, n))) % Update the real part of the wavefunction
    set(imag_wavefunction, 'YData', imag(psi_t(:, n))) % Update the imaginary part of the wavefunction
    set(prob_density, 'YData', abs(psi_t(:, n)).^2); % Update the probability density
    pause(0.1); % Pause to create an animation effect
    drawnow; % Update the relevant figures
end