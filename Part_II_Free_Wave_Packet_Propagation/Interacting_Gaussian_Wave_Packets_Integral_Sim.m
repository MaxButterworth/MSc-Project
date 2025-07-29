% ======================================================================================================================================
%%%%%%%%%% Preamble %%%%%%%%%%
% ======================================================================================================================================

% Part II - Free Particle Gaussian Wave Packet Simulations
% The interaction between two free particle Gaussian wave packets
% Wave packets are constructed using an integral over free particle eigenfunctions
% Time propagation is performed using the Crank-Nicolson method

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
k0_A = 10; % Set the expectation value for k for wave packet A
sigma_A = L/50; % Set the initial width of wave packet A

% Define variables for wave packet B
x0_B = L/4; % Set the starting position of wave packet B on the x-axis
k0_B = 5; % Set the expectation value for k for wave packet B
sigma_B = L/50; % Set the initial width of wave packet B
t_delay = 0; % Set the time delay from when wave packet B should be introduced into the system

% Determine whether periodic boundary conditions are activated or not
set_PBC = false;

% ======================================================================================================================================
%%%%%%%%%% Discretise the spatial domain, x, and time domain, t; determine the k-space domain %%%%%%%%%%
% ======================================================================================================================================

x = linspace(0, L, N_steps); % Define the domain of the infinite potential well
dx = x(2) - x(1); % Calculate the spatial step size

dk = (2 * pi)/L; % Define spacing in k-space

if rem(N_steps, 2) == 0 % Define k-space grid if N_steps is even
    k = dk * (-(N_steps/2):((N_steps/2)-1)).';
    k = ifftshift(k); % Shift the position of zero to operate in k-space and make the calculation compatible with the FFT

else % Define k-space grid if N_steps is odd
    k = dk * (-((N_steps-1)/2):((N_steps-1)/2)).';
    k = ifftshift(k); % Shift the position of zero to operate in k-space and make the calculation compatible with the FFT

end

dt = 1e-2; % Define the time step size
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
psi0_A = ifft(a_k_A .* exp(-1i * k * x0_A)); % Initial Gaussian wave packet in real space

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
psi0_B = ifft(a_k_B .* exp(-1i * k * x0_B)); % Initial Gaussian wave packet in real space

psi0_B_norm = psi0_B/sqrt(trapz(x, abs(psi0_B).^2)); % Normalise the initial Gaussian wave packet

% ======================================================================================================================================
%%%%%%%%%% Implement the Crank-Nicolson method to evolve wave packet A and calculate its probability current %%%%%%%%%%
% ======================================================================================================================================

J = zeros(N_steps, N_t); % Initialise an array to store probability currents
first_deriv = spdiags([-1, 1], 0:1, N_steps, N_steps);

if t_delay == 0
    psi = psi0_A_norm(1:N_steps) + psi0_B_norm(1:N_steps); % Set the initial value of the wavefunction

else
    psi = psi0_A_norm(1:N_steps); % Set the initial value of the wavefunction

end

psi_t = zeros(N_steps, N_t); % Initialise an array to store the wavefunction as it evolves in time
psi_t(:, 1) = psi; % Store the initial wavefunction in the time evolution array

% Pre-compute matrices required for the Crank-Nicolson method
A = eye(N_steps) + (((1i * dt)/(2 * hbar)) * H);
B = eye(N_steps) - (((1i * dt)/(2 * hbar)) * H);

for t = 2:N_t % Loop over all time steps
    if t == t_delay
        psi = psi + psi0_B_norm;
    end

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
    pause(0.1); % Pause to create an animation effect
    drawnow; % Update the relevant figures
end