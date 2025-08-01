% ======================================================================================================================================
%%%%%%%%%% Preamble %%%%%%%%%%
% ======================================================================================================================================

% Part II - Free Particle Gaussian Wave Packet Simulations
% A single free particle Gaussian wave packet constructed using an integral over free particle eigenfunctions
% Time propagation is performed using the split operator method

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

x0 = L/4; % Set the starting position of wave packet on the x-axis
k0 = 10; % Set the expectation value for k for the wave packet
sigma = L/50; % Set the initial width of the wave packet

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
N_t = 300; % Define the number of time steps to simulate

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

figure; % Generate a figure

t_array = dt * (0:N_t - 1); % Create a time array

subplot(2, 2, 1) % Top left subfigure
real_wavefunction = plot(x, real(psi_t(:, 1))); % Plot the real wavefunction
xlabel('$x$', 'Interpreter','latex'); % Label the x-axis
ylabel('$\mathrm{Re}(\psi(x, t))$', 'Interpreter','latex'); % Label the y-axis
xlim([min(x) max(x)]) % Set the y-limits for convenience
ylim([min(real(psi_t(:))) max(real(psi_t(:)))]); % Set the y-limits for convenience
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
ylim([min(abs(psi_t(:)).^2) max(real(abs(psi_t(:)).^2))]); % Set the y-limits for convenience
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

    sgtitle(sprintf('Time Elapsed: %.3f seconds', t_array(n))); % Update time elpased in the overall title for the figure
    
    pause(0.1); % Pause to create an animation effect
    drawnow; % Update the relevant figures
end