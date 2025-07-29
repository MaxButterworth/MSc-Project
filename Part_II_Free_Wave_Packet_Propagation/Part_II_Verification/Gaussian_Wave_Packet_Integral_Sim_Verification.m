% ======================================================================================================================================
%%%%%%%%%% Preamble %%%%%%%%%%
% ======================================================================================================================================

% Part II Verification - Free Particle Gaussian Wave Packet Simulations
% A single free particle Gaussian wave packet
% The wave packet is constructed using an integral over free particle eigenfunctions
% Time propagation is performed using the Crank-Nicolson method
% This script aims to verify numerical time-propagation using the analytical form of the time-propagated wave packet

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

x0 = L/4; % Set the starting position of wave packet on the x-axis
k0 = 10; % Set the expectation value for k for the wave packet
sigma = L/50; % Set the initial width of the wave packet

set_PBC = false; % Determine whether periodic boundary conditions are activated or not

% ======================================================================================================================================
%%%%%%%%%% Discretise the spatial domain, x; time domain, t; and k-space domain, k %%%%%%%%%%
% ======================================================================================================================================

x = linspace(0, L, N_steps); % Define the domain of the infinite potential well
dx = x(2) - x(1); % Calculate the spatial step size

dk = (2 * pi)/L; % Define spacing in k-space

if rem(N_steps, 2) == 0 % Define k-space grid if N_steps is even
    k = dk * (-(N_steps/2):((N_steps/2)-1)).';
    k = ifftshift(k); % Shift the position of zero to operate in k-space and make the calculation compatible with the inverse FFT

else % Define k-space grid if N_steps is odd
    k = dk * (-((N_steps-1)/2):((N_steps-1)/2)).';
    k = ifftshift(k); % Shift the position of zero to operate in k-space and make the calculation compatible with the inverse FFT

end

dt = 1e-3; % Define the time step size
N_t = 1000; % Define the number of time steps to simulate

% ======================================================================================================================================
%%%%%%%%%% Construct the Hamiltonian using the finite difference method %%%%%%%%%%
% ======================================================================================================================================

% Construct the Hamiltonian inside the infinite potential well

if set_PBC == false
    laplacian = (1/dx^2) * spdiags([1, -2, 1], -1:1, N_steps, N_steps); % Define the Laplacian operator for infinitely high boundaries

else
    laplacian = spdiags([1, -2, 1], -1:1, N_steps, N_steps); % Define the Laplacian operator for infinitely high boundaries
    
    % Impose the periodic boundary conditions
    laplacian(1, N_steps) = 1;
    laplacian(N_steps, 1) = 1;
    
    laplacian = (1/dx^2) * laplacian; % Divide by dx^2; define the Laplacian for periodic boundaries

end

H = -((hbar^2)/(2*m)) * laplacian; % Define the Hamiltonian operator

% ======================================================================================================================================
%%%%%%%%%% Generate an initial wave packet by applying the forward Fourier transform to a Gaussian %%%%%%%%%%
% ======================================================================================================================================

a_0 = 1; % Prefactor for the Gaussian distribution
phase = 0; % Define the phase term in the Gaussian distribution
phase_coeff = 0; % Define the coefficients of the components of the phase term

% Construct the Gaussian phase term in k-space
for index = 1:length(phase_coeff)
    phase = phase + (phase_coeff(index) * (k - k0).^(index - 1));
end

a_k = a_0 * exp((-(1/(2 * sigma^2)) * (k - k0).^2) + (1i * phase)); % Construct the whole Gaussian distribution in k-space
psi0 = ifft(a_k .* exp(-1i * k * x0)); % Initial Gaussian wave packet in real space

psi0_norm = psi0/sqrt(trapz(x, abs(psi0).^2)); % Normalise the initial Gaussian wave packet

% ======================================================================================================================================
%%%%%%%%%% Implement the Crank-Nicolson method to evolve the wavefunction and calculate probability current %%%%%%%%%%
% ======================================================================================================================================

J = zeros(N_steps, N_t); % Initialise an array to store probability currents
first_deriv = spdiags([-1, 1], 0:1, N_steps, N_steps)/dx; % Define a matrix to calculte the first derivative using the finite difference method

psi = psi0_norm(1:N_steps); % Set the initial value of the wavefunction
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
%%%%%%%%%% Generate the analytical time-propagated Gaussian wave packet %%%%%%%%%%
% ======================================================================================================================================

psi_analytical_t = zeros(N_steps, 1);
J_analytical = zeros(N_steps, N_t); % Initialise an array to store probability currents

for t_analytical = 0:(N_t - 1)

    % Construct the whole Gaussian distribution in k-space
    a_k_t_analytical = a_0 * exp((-(1/(2 * sigma^2)) * (k - k0).^2) + (1i * phase) - (1i * t_analytical * dt * ((hbar * k.^2)/(2 * m))));

    psi_analytical = ifft(a_k_t_analytical .* exp(-1i * k * x0)); % Generate the Gaussian wave packet in real space
    psi_analytical_norm = psi_analytical/sqrt(trapz(x, abs(psi_analytical).^2)); % Normalise the wave packet
    psi_analytical_t(:, t_analytical + 1) = psi_analytical_norm; % Store the analytical time-evolved wavefunction in the time evolution array

    % Calculate probability current at each point along x
    J_analytical(:, t_analytical + 1) = -((1i * hbar)/(2 * m)) * ((conj(psi_analytical) .* (first_deriv * psi_analytical)) ...
                                        - ((first_deriv * conj(psi_analytical)) .* psi_analytical));
end

% ======================================================================================================================================
%%%%%%%%%% Plot the time evolution of the wave packet probability density %%%%%%%%%%
% ======================================================================================================================================

figure; % Generate a figure

subplot(2, 2, 1) % Top left subfigure
real_wavefunction = plot(x, real(psi_t(:, 1))); % Plot the real wavefunction
hold on
real_wavefunction_analytical = plot(x, real(psi_analytical_t(:, 1)));
hold off
xlabel('$x$', 'Interpreter','latex'); % Label the x-axis
ylabel('$\mathrm{Re}(\psi(x, t))$', 'Interpreter','latex'); % Label the y-axis
xlim([min(x) max(x)]) % Set the y-limits for convenience
ylim([min(real(psi_t(:))) max(real(psi_t(:)))]); % Set the y-limits for convenience
title('Real Component of the Wavefunction') % Add a title
grid on; % Add a grid to the plot

subplot(2, 2, 2) % Top right subfigure
imag_wavefunction = plot(x, imag(psi_t(:, 1))); % Plot the imaginary wavefunction
hold on
imag_wavefunction_analytical = plot(x, imag(psi_analytical_t(:, 1)));
hold off
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
    set(real_wavefunction_analytical, 'YData', real(psi_analytical_t(:, n)))

    set(imag_wavefunction, 'YData', imag(psi_t(:, n))) % Update the imaginary part of the wavefunction
    set(imag_wavefunction_analytical, 'YData', imag(psi_analytical_t(:, n)))

    set(prob_density, 'YData', abs(psi_t(:, n)).^2); % Update the probability density
    set(flux_plot, 'YData', J(:, n)); % Update the probability density
    pause(0.1); % Pause to create an animation effect
    drawnow; % Update the relevant figures
end