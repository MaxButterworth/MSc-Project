% ======================================================================================================================================
%%%%%%%%%% Preamble %%%%%%%%%%
% ======================================================================================================================================

% Part Ic - Free Particle Gaussian Wave Packet Incident on a Potential Barrier Simulation
% Wave packet propagation is performed using the split operator method

% Author: Max L Butterworth
% MSc in Theoretical and Computational Chemistry Project
% University of Oxford

% ======================================================================================================================================
%%%%%%%%%% Define constants and variables %%%%%%%%%%
% ======================================================================================================================================

% Natural units adopted throughout
L = 100; % Length of the 1D box
m = 1; % Mass of electron
h = 1; % Planck's constant
hbar = 1; % Definition of h bar
N_steps = 10001; % Number of discretisation points

wp_energy = 25; % Set the wave packet energy
barrier_energy = 0.1*wp_energy; % Set the magnitude of the potential barrier height
barrier_width = 100; % Set the barrier width in units of dx

x0 = 20; % Set the starting position of wave packet on the x-axis
k0 = sqrt((wp_energy * 2 * m)/(hbar^2)); % Calculate the wavenumber from wp_energy; k = 0 gives a stationary Gaussian wave packet
sigma = L/100; % Set the initial width of the wave packet

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

% Construct the Hamiltonian inside the infinite potential well
laplacian = (1/dx^2) * spdiags([1, -2, 1], -1:1, N_steps, N_steps); % Define the Laplacian operator

if rem(barrier_width, 2) == 0 % If the barrier width is even

    lower_index = (N_steps - barrier_width + 1)/2; % Lower x-index for which the potential barrier starts
    upper_index = (N_steps - barrier_width - 1)/2; % Upper x-index for which the potential barrier ends

    % Define the potential energy accross the x-domain; potential step at halfway accross domain
    V_vector = [zeros(lower_index, 1); repmat(barrier_energy, barrier_width, 1); zeros(upper_index, 1)];

else

    lower_index = (N_steps - barrier_width)/2; % Lower x-index for which the potential barrier starts
    upper_index = (N_steps - barrier_width)/2; % Upper x-index for which the potential barrier ends

    % Define the potential energy accross the x-domain; potential step at halfway accross domain
    V_vector = [zeros(lower_index, 1); repmat(barrier_energy, barrier_width, 1); zeros(upper_index, 1)];

end

V_matrix = diag(V_vector); % Define a diagonal potential matrix

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
%%%%%%%%%% Implement the Crank-Nicolson method to evolve the wavefunction and calculate probability current %%%%%%%%%%
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

figure; % Generate a figure

subplot(2, 2, 1) % Top left subfigure
real_wavefunction = plot(x, real(psi_t(:, 1))); % Plot the real wavefunction
rectangle('Position', [x(lower_index), min(real(psi_t(:))), (barrier_width * dx), ...
         (abs(min(real(psi_t(:)))) + max(real(psi_t(:))))], EdgeColor='red', LineWidth=2) % Illustrate the barrier on the plot

xlabel('$x$', 'Interpreter','latex'); % Label the x-axis
ylabel('$\mathrm{Re}(\psi(x, t))$', 'Interpreter','latex'); % Label the y-axis
xlim([min(x) max(x)]) % Set the x-limits for convenience
ylim([min(real(psi_t(:))) max(real(psi_t(:)))]); % Set the y-limits foe convenience
title('Real Component of the Wavefunction') % Add a title
grid on; % Add a grid to the plot

subplot(2, 2, 2) % Top right subfigure
imag_wavefunction = plot(x, imag(psi_t(:, 1))); % Plot the imaginary wavefunction
rectangle('Position', [x(lower_index), min(imag(psi_t(:))), (barrier_width * dx), ...
         (abs(min(imag(psi_t(:)))) + max(imag(psi_t(:))))], EdgeColor='red', LineWidth=2) % Illustrate the barrier on the plot

xlabel('$x$', 'Interpreter','latex'); % Label the x-axis
ylabel('$\mathrm{Im}(\psi(x, t))$', 'Interpreter','latex'); % Label the y-axis
xlim([min(x) max(x)]) % Set the x-limits for convenience
ylim([min(imag(psi_t(:))) max(imag(psi_t(:)))]); % Set the y-limits for convenience
title('Imaginary Component of the Wavefunction') % Add a title
grid on; % Add a grid to the plot

subplot(2, 2, 3) % Bottom left subfigure
prob_density = plot(x, abs(psi_t(:, 1)).^2); % Plot the initial probability density
rectangle('Position', [x(lower_index), min(abs(psi_t(:)).^2), (barrier_width * dx), ...
         (abs(min(abs(psi_t(:)).^2)) + max(abs(psi_t(:)).^2))], EdgeColor='red', LineWidth=2) % Illustrate the barrier on the plot

xlabel('$x$', 'Interpreter','latex'); % Label the x-axis
ylabel('$|\psi(x, t)|^2$', 'Interpreter','latex'); % Label the y-axis
xlim([min(x) max(x)]) % Set the x-limits for convenience
ylim([min(abs(psi_t(:)).^2) max(abs(psi_t(:)).^2)]); % Set the y-limits for convenience
title('Probability Density') % Add a title
grid on; % Add a grid to the plot

subplot(2, 2, 4) % Bottom right subfigure
flux_plot = plot(x, J(:, 1)); % Plot the initial probability current
rectangle('Position', [x(lower_index), min(J(:)), (barrier_width * dx), ...
         (abs(min(J(:))) + max(J(:)))], EdgeColor='red', LineWidth=2) % Illustrate the barrier on the plot

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
    pause(0.05); % Pause to create an animation effect
    drawnow; % Update the relevant figures
end