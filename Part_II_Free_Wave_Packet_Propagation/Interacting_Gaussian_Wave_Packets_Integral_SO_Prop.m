% ======================================================================================================================================
%%%%%%%%%% Preamble %%%%%%%%%%
% ======================================================================================================================================

% Part II - Free Particle Gaussian Wave Packet Simulations
% The interaction between two free particle Gaussian wave packets
% Wave packets are constructed using an integral over free particle eigenfunctions
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
N_steps = 1000; % Number of discretisation points

% Define variables for wave packet A
x0_A = (3*L)/4; % Set the starting position of wave packet A on the x-axis
k0_A = -10; % Set the expectation value for k for wave packet A
sigma_A = L/50; % Set the initial width of wave packet A

% Define variables for wave packet B
x0_B = L/4; % Set the starting position of wave packet B on the x-axis
k0_B = 10; % Set the expectation value for k for wave packet B
sigma_B = L/50; % Set the initial width of wave packet B
t_delay = 0; % Set the time delay (in units of dt) to specify when wave packet B should be introduced into the system

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
%%%%%%%%%% Generate initial wave packet A by applying the forward Fourier transform to a Gaussian %%%%%%%%%%
% ======================================================================================================================================

a_0_A = 1; % Prefactor for the Gaussian distribution
phase_A = 0; % Define the phase term in the Gaussian distribution
phase_coeff_A = 0; % Define the coefficients of the components of the phase term

% Construct the Gaussian phase term in k-space
for index_A = 1:length(phase_coeff_A)
    phase_A = phase_A + (phase_coeff_A(index_A) * (k - k0_A).^(index_A - 1));
end

a_k_A = a_0_A * exp((-(1/(2 * sigma_A^2)) * (k - k0_A).^2) + (1i * phase_A)); % Construct the whole Gaussian distribution in k-space
psi0_A = ifft(ifftshift(a_k_A .* exp(-1i * k * x0_A))); % Initial Gaussian wave packet in real space

psi0_A_norm = psi0_A/sqrt(trapz(x, abs(psi0_A).^2)); % Normalise the initial Gaussian wave packet

% ======================================================================================================================================
%%%%%%%%%% Generate initial wave packet B by applying the forward Fourier transform to a Gaussian %%%%%%%%%%
% ======================================================================================================================================

a_0_B = 1; % Prefactor for the Gaussian distribution
phase_B = 0; % Define the phase term in the Gaussian distribution
phase_coeff_B = 0; % Define the coefficients of the components of the phase term

% Construct the Gaussian phase term in k-space
for index_B = 1:length(phase_coeff_B)
    phase_B = phase_B + (phase_coeff_B(index_B) * (k - k0_B).^(index_B - 1));
end

a_k_B = a_0_B * exp((-(1/(2 * sigma_B^2)) * (k - k0_B).^2) + (1i * phase_B)); % Construct the whole Gaussian distribution in k-space
psi0_B = ifft(ifftshift(a_k_B .* exp(-1i * k * x0_B))); % Initial Gaussian wave packet in real space

psi0_B_norm = psi0_B/sqrt(trapz(x, abs(psi0_B).^2)); % Normalise the initial Gaussian wave packet

% ======================================================================================================================================
%%%%%%%%%% Implement the split operator method to evolve the total wave packet and calculate its probability current %%%%%%%%%%
% ======================================================================================================================================

% Initialise arrays to store fluxes and a first derivative operator to calculate the fluxes
J = zeros(N_steps, N_t); % Array to store probability currents
first_deriv = spdiags([-1, 1], 0:1, N_steps, N_steps)/dx; % A matrix to calculte the first derivative using the finite difference method

% Define kinetic and potential operators
p = hbar * k; % Calcualte the momentum at each point in k-space

T_op = exp(-(1i * (p.^2) * dt)/(2 * m * hbar)); % Kinetic energy operator (full time step)
V_op = exp(-(1i * V_vector * dt)/(2 * hbar)); % Potential energy operator (half time step)

% Set up the wave packet for propagation
if t_delay == 0
    psi = psi0_A_norm(1:N_steps) + psi0_B_norm(1:N_steps); % Set the initial value of the wave packet

else
    psi = psi0_A_norm(1:N_steps); % Set the initial value of the wave packet

end

psi_t = zeros(N_steps, N_t); % Initialise an array to store the wave packet as it evolves in time
psi_t(:, 1) = psi; % Store the initial wave packet in the time evolution array

% Propagate the wave packet
for t = 2:N_t
    if t == (t_delay + 1)
        psi = psi + psi0_B_norm;
    end

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
%%%%%%%%%% Plot the time evolution of the wave packet, and its probability density and flux %%%%%%%%%%
% ======================================================================================================================================

figure; % Generate a figure

t_array = dt * (0:N_t - 1); % Create a time array

subplot(2, 2, 1) % Top left subfigure
real_wavefunction = plot(x, real(psi_t(:, 1))); % Plot the real wave packet
xlabel('$x$', 'Interpreter','latex'); % Label the x-axis
ylabel('$\mathrm{Re}(\psi(x, t))$', 'Interpreter','latex'); % Label the y-axis
xlim([min(x) max(x)]) % Set the x-limits for convenience
ylim([min(real(psi_t(:))) max(real(psi_t(:)))]); % Set the y-limits for convenience
title('Real Component of the Wave Packet') % Add a title
grid on; % Add a grid to the plot

subplot(2, 2, 2) % Top right subfigure
imag_wavefunction = plot(x, imag(psi_t(:, 1))); % Plot the imaginary wave packet
xlabel('$x$', 'Interpreter','latex'); % Label the x-axis
ylabel('$\mathrm{Im}(\psi(x, t))$', 'Interpreter','latex'); % Label the y-axis
xlim([min(x) max(x)]) % Set the x-limits for convenience
ylim([min(imag(psi_t(:))) max(imag(psi_t(:)))]); % Set the y-limits for convenience
title('Imaginary Component of the Wave Packet') % Add a title
grid on; % Add a grid to the plot

subplot(2, 2, 3) % Bottom left subfigure
prob_density = plot(x, abs(psi_t(:, 1)).^2); % Plot the initial probability density
xlabel('$x$', 'Interpreter','latex'); % Label the x-axis
ylabel('$|\psi(x, t)|^2$', 'Interpreter','latex'); % Label the y-axis
xlim([min(x) max(x)]) % Set the x-limits for convenience
ylim('auto') % Set the y-limits for convenience
title('Probability Density') % Add a title
grid on; % Add a grid to the plot

subplot(2, 2, 4) % Bottom right subfigure
flux_plot = plot(x, J(:, 1)); % Plot the initial probability current
xlabel('$x$', 'Interpreter', 'latex'); % Label the x-axis
ylabel('$J(x, t)$', 'Interpreter', 'latex'); % Label the y-axis
ylim('auto') % Set the y-limits for convenience
title('Probability Current') % Add a title
grid on; % Add a grid to the plot

% Animate the figures

for n = 1:N_t % Loop over all timesteps
    set(real_wavefunction, 'YData', real(psi_t(:, n))) % Update the real part of the wave packet
    set(imag_wavefunction, 'YData', imag(psi_t(:, n))) % Update the imaginary part of the wave packet
    set(prob_density, 'YData', abs(psi_t(:, n)).^2); % Update the probability density
    set(flux_plot, 'YData', J(:, n)); % Update the probability density

    sgtitle(sprintf('Time Elapsed: %.3f seconds', t_array(n))); % Update time elapsed in the overall title for the figure

    pause(0.1); % Pause to create an animation effect
    drawnow; % Update the relevant figures
end