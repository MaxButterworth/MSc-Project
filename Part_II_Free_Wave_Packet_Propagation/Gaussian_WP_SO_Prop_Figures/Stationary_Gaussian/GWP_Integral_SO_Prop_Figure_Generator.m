% ======================================================================================================================================
%%%%%%%%%% Preamble %%%%%%%%%%
% ======================================================================================================================================

% Part II - Free Particle Gaussian Wave Packet Simulations
% A single free particle Gaussian wave packet constructed using an integral over free particle eigenfunctions
% Time propagation is performed using the split operator method
% Figure generator for thesis

% Author: Max L Butterworth
% MSc in Theoretical and Computational Chemistry Project
% University of Oxford

% ======================================================================================================================================
%%%%%%%%%% Define constants and variables %%%%%%%%%%
% ======================================================================================================================================

% Natural units adopted throughout
L = 50; % Length of the 1D box
m = 1; % Mass of electron
h = 1; % Planck's constant
hbar = 1; % Definition of h bar
N_steps = 1001; % Number of discretisation points

x0 = L/2; % Set the starting position of wave packet on the x-axis
k0 = 0; % Set the expectation value for k for the wave packet
sigma = L/50; % Set the initial width of the wave packet

include_elapsed_time = false; % Define a variable to show elapsed time on figure or not

save_figures = false; % Define a variable to save figures at various points in the simulation or not

% ======================================================================================================================================
%%%%%%%%%% Discretise the spatial domain, x; time domain, t; and k-space domain, k %%%%%%%%%%
% ======================================================================================================================================

x = linspace(0, L, N_steps); % Define the x-domain
dx = x(2) - x(1); % Calculate the spatial step size

dk = (2 * pi)/L; % Define spacing in k-space

if rem(N_steps, 2) == 0 % Define k-space grid if N_steps is even
    k = dk * (-(N_steps/2):((N_steps/2)-1)).';

else % Define k-space grid if N_steps is odd
    k = dk * (-((N_steps-1)/2):((N_steps-1)/2)).';

end

dt = 1e-2; % Define the time step size
N_t = 1000; % Define the number of time steps to simulate

% ======================================================================================================================================
%%%%%%%%%% Construct the Hamiltonian using the finite difference method %%%%%%%%%%
% ======================================================================================================================================

laplacian = (1/dx^2) * spdiags([1, -2, 1], -1:1, N_steps, N_steps); % Define the Laplacian operator for infinitely high boundaries

V_vector = zeros(N_steps, 1); % Define the potential energy vector for a free wave packet
V_matrix = zeros(N_steps, N_steps); % Define the potential energy matrix for a free wave packet

H = (-((hbar^2)/(2*m)) * laplacian) + V_matrix; % Define the Hamiltonian operator

% ======================================================================================================================================
%%%%%%%%%% Generate an initial wave packet by applying the inverse Fourier transform to a Gaussian %%%%%%%%%%
% ======================================================================================================================================

a_0 = 1; % Prefactor for the Gaussian distribution
phase = 0; % Define the phase term in the Gaussian distribution
phase_coeff = 0; % Define the coefficients of the components of the phase term

% Construct the Gaussian phase term in k-space
for index = 1:length(phase_coeff)
    phase = phase + (phase_coeff(index) * (k - k0).^(index - 1));
end

a_k = a_0 * exp((-(1/(2 * sigma^2)) * (k - k0).^2) + (1i * phase)); % Construct the whole Gaussian distribution in k-space
psi0 = ifft(ifftshift(a_k .* exp(-1i * k * x0))); % Initial Gaussian wave packet in real space

psi0_norm = psi0/sqrt(trapz(x, abs(psi0).^2)); % Normalise the initial Gaussian wave packet

% ======================================================================================================================================
%%%%%%%%%% Implement the split operator method to evolve the wave packet and calculate the probability current %%%%%%%%%%
% ======================================================================================================================================

% Define kinetic and potential operators
p = hbar * k; % Calcualte the momentum at each point in k-space

T_op = exp(-(1i * (p.^2) * dt)/(2 * m * hbar)); % Kinetic energy operator (full time step)
V_op = exp(-(1i * V_vector * dt)/(2 * hbar)); % Potential energy operator (half time step)

% Initialise arrays to store fluxes and a first derivative operator to calculate the fluxes
J = zeros(N_steps, N_t); % Array to store probability currents
first_deriv = spdiags([-1, 1], 0:1, N_steps, N_steps)/dx; % A matrix to calculte the first derivative using the finite difference method

% Set up the wave packet for propagation
psi = psi0_norm(1:N_steps); % Set the initial value of the wave packet
psi_t = zeros(N_steps, N_t); % Initialise an array to store the wave packet as it evolves in time
psi_t(:, 1) = psi; % Store the initial wave packet in the time evolution array

% Propagate the wave packet
for t = 2:N_t
    psi = V_op .* psi; % Operate a half time step in real space
    psi_k = fftshift(fft(psi)); % Fourier transform the wave packet into k-space
    psi_k = T_op .* psi_k; % Operate a full time step in k-space
    psi = ifft(ifftshift(psi_k)); % Inverse Fourier transform into real space
    psi = V_op .* psi; % Operate a half time step in real space

    psi = psi/sqrt(trapz(x, abs(psi).^2)); % Normalise the time-evolved wave packet

    psi_t(:, t) = psi; % Store the time-evolved wave packet in the time evolution array
    J(:, t) = -((1i * hbar)/(2 * m)) * ((conj(psi) .* (first_deriv * psi)) - ((first_deriv * conj(psi)) .* psi)); % Calculate probability current at each point along x
end

% ======================================================================================================================================
%%%%%%%%%% Plot the time evolution of the wave packet probability density %%%%%%%%%%
% ======================================================================================================================================

time_indices_plot = [1, 141, 741]; % Define the time indices which data are to be obtained for
t_array = dt * (0:N_t - 1); % Create a time array for the simulation

% Initialise arrays for plotting
WP_re_part = zeros(size(x, 2), size(time_indices_plot, 2)); % Initialise an array to store real parts of wave packet
WP_im_part = zeros(size(x, 2), size(time_indices_plot, 2)); % Initialise an array to store imaginary parts of wave packet
prob_density_plot = zeros(size(x, 2), size(time_indices_plot, 2)); % Initialise an array to store the wave packet probability density
J_plot = zeros(size(x, 2), size(time_indices_plot, 2)); % Initialise an array to store the flux
time_plot = zeros(1, size(time_indices_plot, 2)); % Initialise an array to store imaginary parts of wave packet

counter = 1; % Initialise a counter to append data to arrays

% Search for the data required based on the time index of the simulation
for n = time_indices_plot
    time_plot(1, counter) = t_array(1, n); % Append time to the time array
    WP_re_part(:, counter) = real(psi_t(:, n)); % Append real part of wave packet to the relevant array
    WP_im_part(:, counter) = imag(psi_t(:, n)); % Append imaginary part of wave packet to the relevant array
    prob_density_plot(:, counter) = abs(psi_t(:, n)).^2; % Append probability density to the relevant array
    J_plot(:, counter) = J(:, n); % Append the flux to the relevant array

    counter = counter + 1; % Increment the counter by one
end

C = orderedcolors('gem'); % Set the colour of the plots

t = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact'); % Generate a figure

set(groot, 'DefaultAxesFontSize', 14); % Set the font size for axes
set(groot, 'DefaultTextFontSize', 14); % Set the font size for other text

ax1 = nexttile([1 2]); % Top Subfigure
% yyaxis('left')

hold on;
for index_re_plot = 1:size(WP_re_part, 2)
    plot(ax1, x, WP_re_part(:, index_re_plot), 'LineWidth', 2, 'LineStyle', '-', 'Color', C(index_re_plot, :)); % Plot the real wavefunction
end
hold off;

% ax1.YColor = 'k'; % Set the colour of the first y-axis
xlim(ax1, [min(x(:)) max(x(:))]); % Set the x-limits for convenience
ylabel(ax1, '$\mathrm{Re}\psi(x, t)$', 'Interpreter','latex'); % Label the wavefunction y-axis
ylim(ax1, [min(WP_re_part(:)) max(WP_re_part(:))]); % Set the y-limits for wavefunction plot

% Uncomment all the below when the imaginary part of the wave packet needs to be plotted
% yyaxis('right')
% 
% hold on;
% for index_im_plot = 1:size(WP_im_part, 2)
%    plot(ax1, x, WP_im_part(:, index_im_plot), 'LineWidth', 2, 'LineStyle', '-', 'Color', C(index_im_plot, :)); % Plot the real wavefunction
% end
% hold off;
% 
% ax1.YColor = 'k'; % Set the colour of the second y-axis
% ylabel(ax1, '$\mathrm{Im}\psi(x, t)$', 'Interpreter','latex'); % Label the wavefunction y-axis
% ylim(ax1, [min(WP_im_part(:)) max(WP_im_part(:))]); % Set the y-limits for wavefunction plot

xlabel(ax1, '$x$', 'Interpreter','latex'); % Label the x-axis
grid on; % Add a grid to the plot

set(groot, 'DefaultAxesFontSize', 16); % Set the font size for axes
set(groot, 'DefaultTextFontSize', 16); % Set the font size for other text

ax2 = nexttile; % Bottom Left Subfigure

hold on;
for index_prob_density_plot = 1:size(prob_density_plot, 2)
    plot(ax2, x, prob_density_plot(:, index_prob_density_plot), 'LineWidth', 2, 'LineStyle', '-', 'Color', C(index_prob_density_plot, :)); % Plot the initial probability density
end
hold off;

xlabel(ax2, '$x$', 'Interpreter','latex'); % Label the x-axis
ylabel(ax2, '$|\psi(x, t)|^2$', 'Interpreter','latex'); % Label the y-axis
xlim(ax2, [min(x(:)) max(x(:))]); % Set the x-limits for convenience
ylim(ax2, [min(prob_density_plot(:)) max(prob_density_plot(:))]); % Set the y-limits for convenience
grid on; % Add a grid to the plot

ax3 = nexttile; % Bottom Left Subfigure

hold on;
for index_J_plot = 1:size(J_plot, 2)
    plot(ax3, x, J_plot(:, index_J_plot), 'LineWidth', 2, 'LineStyle', '-', 'Color', C(index_J_plot, :)); % Plot the initial probability current
end
hold off;

xlabel(ax3, '$x$', 'Interpreter', 'latex'); % Label the x-axis
ylabel(ax3, '$J(x, t)$', 'Interpreter', 'latex'); % Label the y-axis
xlim(ax3, [min(x(:)) max(x(:))]); % Set the x-limits for convenience
ylim(ax3, [min(J_plot(:)) max(J_plot(:))]); % Set the y-limits for convenience
grid on; % Add a grid to the plot