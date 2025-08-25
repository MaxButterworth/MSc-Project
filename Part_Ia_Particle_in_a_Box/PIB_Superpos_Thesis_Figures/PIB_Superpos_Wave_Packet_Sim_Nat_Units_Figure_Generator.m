% ======================================================================================================================================
%%%%%%%%%% Preamble %%%%%%%%%%
% ======================================================================================================================================

% Part Ia - Particle in a Box Wave Packet Simulation
% Superposition of particle in a box eigenstates to generate a wave packet
% Propagation of the wavefunction is performed using the Crank-Nicolson Method
% This code generates the figures for the thesis

% Author: Max L Butterworth
% MSc in Theoretical and Computational Chemistry Project
% University of Oxford

% ======================================================================================================================================
%%%%%%%%%% Define constants and variables %%%%%%%%%%
% ======================================================================================================================================

% Natural units have been adopted throughout
L = 1; % Length of the 1D box
m = 1; % Mass
h = 1; % Planck's constant in J s
hbar = 1; % Definition of h bar

N_steps = 7500; % Number of discretisation points

basis_funcs_indices = [1, 2]; % Create an array of the indices of PIB_eigenstates_norm that form the superposition
basis_funcs_coeffs = [sqrt(1/2), sqrt(1/2)]; % Equal Weightings of PIB eigenstates in the superposition
N_PIB_eigenfuncs = max(basis_funcs_indices); % The number of basis functions in the wave packet superposition

save_figures = true; % Variable to determine whether figures are saved or not

% ======================================================================================================================================
%%%%%%%%%% Discretise the spatial domain, x, and time domain, t %%%%%%%%%%
% ======================================================================================================================================

x = linspace(0, L, N_steps); % Define the domain of the infinite potential well
dx = x(2) - x(1); % Calculate the spatial step size

N_t = 1001; % Define the number of time steps to simulate

dt = (4 * m * L^2)/(3 * pi * hbar * N_t); % Define the time step size

% ======================================================================================================================================
%%%%%%%%%% Solve the Schrödinger equation using the finite difference method %%%%%%%%%%
% ======================================================================================================================================

% Construct the Hamiltonian inside the infinite potential well
laplacian = (1/dx^2) * spdiags([1, -2, 1], -1:1, N_steps, N_steps); % Define the Laplacian operator
H = -((hbar^2)/(2*m)) * laplacian; % Define the Hamiltonian operator

% Find the eigenvalues and eigenvectors of the Hamiltonian matrix
[PIB_eigenstates, E_PIB] = eigs(H, N_PIB_eigenfuncs, 'smallestabs');

PIB_eigenstates_norm = zeros(N_steps, N_PIB_eigenfuncs); % Set up an array to store normalised PIB eigenfunctions

% Normalise eigenvectors
for i = 1:N_PIB_eigenfuncs
    PIB_eigenstates_norm(:, i) = PIB_eigenstates(:, i)/sqrt(trapz(x, abs(PIB_eigenstates(:, i)).^2));
end

% ======================================================================================================================================
%%%%%%%%%% Generate an initial wave packet composed of a superposition of PIB eigenfunctions modulated by a Gaussian %%%%%%%%%%
% ======================================================================================================================================

x0 = L/2; % Start evolving the wave packet from the centre of the box at t = 0
sigma = L/20; % Set the initial width of the wave packet

psi0 = zeros(N_steps, 1); % Initialise an empty array to store the initial wave packet

% Generate the superposition of PIB basis functions
for l = 1:length(basis_funcs_indices)
    psi0 = psi0 + (basis_funcs_coeffs(l) * PIB_eigenstates_norm(:, basis_funcs_indices(l)));
end

psi0_norm = psi0/sqrt(trapz(x, abs(psi0).^2)); % Normalise the initial Gaussian wave packet

% ======================================================================================================================================
%%%%%%%%%% Impose boundary conditions %%%%%%%%%%
% ======================================================================================================================================

% Set the wavefunction to zero at the boundaries
psi0_norm(1) = 0;
psi0_norm(N_steps) = 0;
psi0_norm = sparse(psi0_norm); % Define the initial wavefunction as a sparse matrix to speed up the calculation

H = H(2:N_steps-1, 2:N_steps-1); % Impose boundary conditions: psi(0) = psi(L) = 0
H = sparse(H); % Define the Hamiltonian as a sparse matrix to speed up the calcualtion

x_internal = x(2:N_steps - 1); % Truncate the x array to account for boundary conditions

% ======================================================================================================================================
%%%%%%%%%% Implement the Crank-Nicolson method to evolve the wavefunction and calculate probability current %%%%%%%%%%
% ======================================================================================================================================

J = zeros(N_steps, N_t); % Initialise an array to store probability currents
first_deriv = spdiags([-1, 1], 0:1, N_steps-2, N_steps-2)/dx; % Define a first derivative operator using the finite difference method

psi = psi0_norm(2:N_steps-1); % Set the initial value of the wavefunction defined on the internal coordinates
psi_t = zeros(N_steps, N_t); % Initialise an array to store the wavefunction as it evolves in time
psi_t(2:N_steps-1, 1) = psi; % Store the initial wavefunction in the time evolution array

% Pre-compute matrices required for the Crank-Nicolson method
A = speye(N_steps-2) + (((1i * dt)/(2 * hbar)) * H);
B = speye(N_steps-2) - (((1i * dt)/(2 * hbar)) * H);

for t = 2:N_t % Loop over all time steps
    psi = A \ (B * psi); % Evolve the wavefunction defined on the internal coordinates over time
    psi = psi/sqrt(trapz(x_internal, abs(psi).^2)); % Normalise the time-evolved wavefunction
    psi_t(2:N_steps-1, t) = psi; % Store the time-evolved wavefunction in the time evolution array
    J(2:N_steps - 1, t) = -((1i * hbar)/(2 * m)) * ((conj(psi) .* (first_deriv * psi)) - ((first_deriv * conj(psi)) .* psi)); % Calculate probability current at each point along x
end

% ======================================================================================================================================
%%%%%%%%%% Plot the time evolution of the wave packet, probability density, and flux %%%%%%%%%%
% ======================================================================================================================================

figure; % Generate a figure

t_array = dt * (0:N_t - 1); % Create a time array

subplot(2, 2, 1) % Top Left subfigure

yyaxis('left')
real_wavefunction = plot(x, real(psi_t(:, 1)), 'LineWidth', 3); % Plot the real wavefunction
hold on
imag_wavefunction = plot(x, imag(psi_t(:, 1)), 'LineWidth', 3); % Plot the imaginary wavefunction
ylabel('$\psi(x, t)$', 'Interpreter','latex'); % Label the wavefunction y-axis
ylim([min(min(real(psi_t(:))), min(imag(psi_t(:)))) max(max(real(psi_t(:))), max(imag(psi_t(:))))]); % Set the y-limits for wavefunction plot

yyaxis('right')
prob_density = plot(x, abs(psi_t(:, 1)).^2, 'LineWidth', 3); % Plot the initial probability density
ylabel('$|\psi(x, t)|^2$', 'Interpreter','latex'); % Label the prob density y-axis
ylim([min(abs(psi_t(:)).^2) max(abs(psi_t(:)).^2)]); % Set the y-limits of prob density axis
hold off

xlabel('$x$', 'Interpreter','latex'); % Label the x-axis
grid on; % Add a grid to the plot

set(groot, 'DefaultAxesFontSize', 24); % Set the font size for axes
set(groot, 'DefaultTextFontSize', 24); % Set the font size for other text

% Animate the figures

for n = 1:N_t % Loop over all timesteps
    set(real_wavefunction, 'YData', real(psi_t(:, n))) % Update the real part of the wavefunction
    set(imag_wavefunction, 'YData', imag(psi_t(:, n))) % Update the imaginary part of the wavefunction
    set(prob_density, 'YData', abs(psi_t(:, n)).^2); % Update the probability density
    
    pause(0.005); % Pause to create an animation
    drawnow; % Update the figures and display immediately
    
    if save_figures == true
        if ismember(n, [1, 250, 500, 750])
            time = t_array(1, n); % Assign the current time to a variable
            filename = sprintf('PIB_2_State_WP_%.2f.png', time); % Create the file name for the figure
            exportgraphics(gcf, filename, 'ContentType', 'image', 'Resolution', 300); % Save the figure
    
        end
    end
end