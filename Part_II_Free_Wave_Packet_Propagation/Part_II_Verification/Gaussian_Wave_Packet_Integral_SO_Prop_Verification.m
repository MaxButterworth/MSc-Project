% ======================================================================================================================================
%%%%%%%%%% Preamble %%%%%%%%%%
% ======================================================================================================================================

% Part II Verification - Free Particle Gaussian Wave Packet Simulations
% A single free particle Gaussian wave packet
% The wave packet is constructed using an integral over free particle eigenfunctions
% Time propagation is performed using the split operator method
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

% ======================================================================================================================================
%%%%%%%%%% Discretise the spatial domain, x; time domain, t; and k-space domain, k %%%%%%%%%%
% ======================================================================================================================================

x = linspace(0, L, N_steps); % Define the domain of the infinite potential well
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
%%%%%%%%%% Generate the analytical time-propagated Gaussian wave packet %%%%%%%%%%
% ======================================================================================================================================

psi_analytical_t = zeros(N_steps, 1);
J_analytical = zeros(N_steps, N_t); % Initialise an array to store probability currents

for t_analytical = 0:(N_t - 1)

    % Construct the whole Gaussian distribution in k-space
    a_k_t_analytical = a_0 * exp((-(1/(2 * sigma^2)) * (k - k0).^2) + (1i * phase) - (1i * t_analytical * dt * ((hbar * k.^2)/(2 * m))));

    psi_analytical = ifft(ifftshift(a_k_t_analytical .* exp(-1i * k * x0))); % Generate the Gaussian wave packet in real space
    psi_analytical_norm = psi_analytical/sqrt(trapz(x, abs(psi_analytical).^2)); % Normalise the wave packet
    psi_analytical_t(:, t_analytical + 1) = psi_analytical_norm; % Store the analytical time-evolved wave packet in the time evolution array

    % Calculate probability current at each point along x
    J_analytical(:, t_analytical + 1) = -((1i * hbar)/(2 * m)) * ((conj(psi_analytical) .* (first_deriv * psi_analytical)) ...
                                        - ((first_deriv * conj(psi_analytical)) .* psi_analytical));
end

% ======================================================================================================================================
%%%%%%%%%% Calculate the measures of error on the split operator propagation method %%%%%%%%%%
% ======================================================================================================================================

% Arrays to store errors over time
norm_squared_error_t = zeros(N_steps, N_t); % Array to store the norm sqaured error over time
x_avg_error_t = zeros(1, N_t); % Array to store the error on the average position over time
overlap_squared_t = zeros(1, N_t); % Array to store the value of the overlap of the wave packet
group_velocity_error = zeros(1, N_t); % Array to store difference in group velocity over time
dispersion_error = zeros(1, N_t); % Array to store difference in dispersion over time

for error_index = (1:N_t)
    % Norm squared of the error on the values of the real and imaginary parts of the C-N propagated wave packet
    norm_squared_error_t(:, N_steps) = abs(psi_t(:,error_index) - psi_analytical_t(:,error_index)).^2;
    
    % Error on the average position
    x_avg_num = trapz(x, (x.' .* abs(psi_t(:, error_index)).^2)); % Average position for the numerical wave packet
    x_avg_anal = trapz(x, (x.' .* abs(psi_analytical_t(:, error_index)).^2)); % Average position for the numerical wave packet

    x_avg_error_t(1, error_index) = x_avg_num - x_avg_anal; % Error on the average position

    % Calculate the overlap between the two wave packets
    overlap_squared = abs(trapz(x, (conj(psi_t(:,error_index)) .* psi_analytical_t(:, error_index))))^2; % Calculate the squared overlap
    overlap_squared_t(1, error_index) = overlap_squared; % Assign the squared overlap to the overlap_squared array

    % Calculate the difference in group velocity over time
    k_avg_num = trapz(k, (k .* abs(a_k).^2)); % Average k-value for the numerical wave packet
    k_avg_anal = trapz(k, (k .* abs(a_k_t_analytical).^2)); % Average k-value for the analytical wave packet
    group_velocity_error(1, error_index) = (hbar/m) * (k_avg_num - k_avg_anal); % Assign difference in group velocity to the group_velocity_error array

    % Calculate the difference in dispersion over time
    x_sq_avg_num = trapz(x, ((x.^2).' .* abs(psi_t(:, error_index)).^2)); % Average squared position for the numerical wave packet
    x_sq_avg_anal = trapz(x, ((x.^2).' .* abs(psi_analytical_t(:, error_index)).^2)); % Average squared position for the analytical wave packet
    
    % Assign the difference in dispersion to the dispersion_error arrray
    dispersion_error(1, error_index) = sqrt((x_sq_avg_num - (x_avg_num.^2))) - sqrt((x_sq_avg_anal - (x_avg_anal.^2)));

end

% ======================================================================================================================================
%%%%%%%%%% Plot the time evolution of the wave packet probability density %%%%%%%%%%
% ======================================================================================================================================

figure; % Generate a figure

t_array = dt * (0:N_t - 1); % Create a time array

subplot(3, 2, 1) % Top left subfigure
real_wavefunction = plot(x, real(psi_t(:, 1))); % Plot the real component of the wave packet
hold on
real_wavefunction_analytical = plot(x, real(psi_analytical_t(:, 1)));
hold off
xlabel('$x$', 'Interpreter','latex'); % Label the x-axis
ylabel('$\mathrm{Re}(\psi(x, t))$', 'Interpreter','latex'); % Label the y-axis
xlim([min(x) max(x)]) % Set the y-limits for convenience
ylim([min(real(psi_t(:))) max(real(psi_t(:)))]); % Set the y-limits for convenience
title('Real Component of the Wave Packet', 'Interpreter', 'latex') % Add a title
grid on; % Add a grid to the plot
legend('Numerical Wave Packet', 'Analytical Wave Packet')

subplot(3, 2, 2) % Top middle subfigure
imag_wavefunction = plot(x, imag(psi_t(:, 1))); % Plot the imaginary component of the wave packet
hold on
imag_wavefunction_analytical = plot(x, imag(psi_analytical_t(:, 1)));
hold off
xlabel('$x$', 'Interpreter','latex'); % Label the x-axis
ylabel('$\mathrm{Im}(\psi(x, t))$', 'Interpreter','latex'); % Label the y-axis
xlim([min(x) max(x)]) % Set the y-limits for convenience
ylim([min(imag(psi_t(:))) max(imag(psi_t(:)))]); % Set the y-limits for convenience
title('Imaginary Component of the Wave Packet', 'Interpreter', 'latex') % Add a title
grid on; % Add a grid to the plot
legend('Numerical Wave Packet', 'Analytical Wave Packet')

subplot(3, 2, 3) % Top right subfigure
norm_squared_error_plot = plot(t_array, x_avg_error_t.'); % Plot the error on the real component of the wave packet
xlabel('$t$', 'Interpreter','latex'); % Label the x-axis
ylabel('$\Delta\langle x\rangle$', 'Interpreter','latex'); % Label the y-axis
ylim([min(x_avg_error_t(:)) max(x_avg_error_t(:))]); % Set the y-limits for convenience
title('Error on the Average Position', 'Interpreter', 'latex') % Add a title
grid on; % Add a grid to the plot

subplot(3, 2, 4) % Bottom left subfigure
overlap_squared_plot = plot(t_array, overlap_squared_t.'); % Plot the error on the imaginary component of the wave packet
xlabel('$t$', 'Interpreter', 'latex'); % Label the x-axis
ylabel('$\langle\psi_\mathrm{num}|\psi_\mathrm{anal}\rangle$', 'Interpreter', 'latex'); % Label the y-axis
title('Square Overlap of the Wave Packets', 'Interpreter', 'latex') % Add a title
grid on; % Add a grid to the plot

subplot(3, 2, 5) % Bottom middle subfigure
grou_velocity_error_plot = plot(t_array, group_velocity_error.'); % Plot the error on the imaginary component of the wave packet
xlabel('$x$', 'Interpreter', 'latex'); % Label the x-axis
ylabel('$\Delta v_g$', 'Interpreter', 'latex'); % Label the y-axis
title('Error on the Group Velocity', 'Interpreter', 'latex') % Add a title
grid on; % Add a grid to the plot

subplot(3, 2, 6) % Bottom right subfigure
error_imag_plot = plot(t_array, dispersion_error.'); % Plot the error on the imaginary component of the wave packet
xlabel('$t$', 'Interpreter', 'latex'); % Label the x-axis
ylabel('$\Delta(\Delta x)$', 'Interpreter', 'latex'); % Label the y-axis
title('Error on the Dispersion, $\Delta x$', 'Interpreter', 'latex') % Add a title
grid on; % Add a grid to the plot

% Animate the figures

for n = 1:N_t % Loop over all timesteps
    set(real_wavefunction, 'YData', real(psi_t(:, n))) % Update the real part of the numerical wave packet
    set(real_wavefunction_analytical, 'YData', real(psi_analytical_t(:, n))) % Update the real part of the analytical wave packet

    set(imag_wavefunction, 'YData', imag(psi_t(:, n))) % Update the imaginary part of the numerical wave packet
    set(imag_wavefunction_analytical, 'YData', imag(psi_analytical_t(:, n))) % Update the real part of the analytical wave packet

    sgtitle(sprintf('Time Elapsed: %.3f seconds', t_array(n))); % Update time in the overall title for the figure

    pause(0.01); % Pause to create an animation effect
    drawnow; % Update the relevant figures

end