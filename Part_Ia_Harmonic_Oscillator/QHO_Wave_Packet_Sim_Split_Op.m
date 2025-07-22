% ======================================================================================================================================
%%%%%%%%%% Preamble %%%%%%%%%%
% ======================================================================================================================================

% Part Ia - Quantum Harmonic Oscillator Wave Packet Simulation
% Superposition of quantum harmonic oscillator eigenstates
% Time Propagation conducted with the split operator method

% Author: Max L Butterworth
% MSc in Theoretical and Computational Chemistry Project
% University of Oxford

% ======================================================================================================================================
%%%%%%%%%% Define constants and variables %%%%%%%%%%
% ======================================================================================================================================

% Natural units have been adopted throughout
L = 5; % Length of the 1D box
m = 1; % Mass
h = 1; % Planck's constant in J s
hbar = 1; % Definition of h bar

N_steps = 10000; % Number of discretisation points on the x-axis

basis_funcs_indices = [1, 2]; % Create an array of the indices of PIB_eigenstates_norm that form the superposition
basis_funcs_coeffs = rand(1, length(basis_funcs_indices)); % Weightings of PIB eigenstates in the superposition
N_PIB_eigenfuncs = max(basis_funcs_indices); % The number of basis functions in the wave packet superposition

% ======================================================================================================================================
%%%%%%%%%% Discretise the spatial domain, x, and time domain, t %%%%%%%%%%
% ======================================================================================================================================

x = linspace(-L, L, N_steps); % Define the domain of the infinite potential well
dx = abs(x(2) - x(1)); % Calculate the spatial step size

dt = 5e-2; % Define the time step size
N_t = 1000; % Define the number of time steps to simulate

% ======================================================================================================================================
%%%%%%%%% Solve the Schrödinger equation using the finite difference method %%%%%%%%%%
% ======================================================================================================================================

% Construct the Hamiltonian inside the infinite potential well
laplacian = (1/dx^2) * spdiags([1, -2, 1], -1:1, N_steps, N_steps); % Define the Laplacian operator

force_const = 1; % Define the force constant for the QHO
V_vector = 0.5 * force_const * (x.').^2; % Define the potential energy vector for a QHO
V_matrix = diag(V_vector); % Define the potential energy matrix for a QHO

H = (-((hbar^2)/(2*m)) * laplacian) + V_matrix; % Define the Hamiltonian operator

% Find the eigenvalues and eigenvectors of the Hamiltonian matrix
[PIB_eigenstates, E_PIB] = eigs(H, N_PIB_eigenfuncs, 'smallestabs');
 
PIB_eigenstates_norm = zeros(N_steps, N_PIB_eigenfuncs); % Set up an array to store normalised PIB eigenfunctions

% Normalise eigenvectors
for i = 1:N_PIB_eigenfuncs
    PIB_eigenstates_norm(:, i) = PIB_eigenstates(:, i)/sqrt(trapz(x, abs(PIB_eigenstates(:, i)).^2));
end

% ======================================================================================================================================
%%%%%%%%%% Generate an initial superposition of QHO eigenfunctions %%%%%%%%%%
% ======================================================================================================================================

psi0 = zeros(N_steps, 1); % Initialise an empty array to store the initial wave packet

% Generate the superposition of PIB basis functions
for l = 1:length(basis_funcs_indices)
    psi0 = psi0 + (basis_funcs_coeffs(l) * PIB_eigenstates_norm(:, basis_funcs_indices(l)));
end

psi0_norm = psi0/sqrt(trapz(x, abs(psi0).^2)); % Normalise the initial Gaussian wave packet

% ======================================================================================================================================
%%%%%%%%%% Propagate the wave packet through time using the split operator method %%%%%%%%%%
% ======================================================================================================================================

% Define kinetic and potential operators
dk = pi/L; % Define spacing in k-space

if rem(N_steps, 2) == 0 % Define k-space grid if N_steps is even
    k = dk * (-(N_steps/2):((N_steps/2)-1)).';

else % Define k-space grid if N_steps is odd
    k = dk * (-((N_steps-1)/2):((N_steps-1)/2)).';

end

p = hbar * k; % Calcualte the momentum at each point in k-space

T_op = exp(-(1i * (p.^2) * dt)/(2 * m * hbar)); % Kinetic energy operator (full time step)
V_op = exp(-(1i * V_vector * dt)/(2 * hbar)); % Potential energy operator (half time step)

% Initialise arrays to store fluxes and a first derivative operator to calculate the fluxes
J = zeros(N_steps, N_t); % Initialise an array to store probability currents
first_deriv = spdiags([-1, 1], 0:1, N_steps, N_steps)/dx; % Define a first derivative operator using the finite difference method

% Set up the wave packet for propagation
psi = psi0_norm; % Set the initial value of the wavefunction
psi_t = zeros(N_steps, N_t); % Initialise an array to store the wavefunction as it evolves in time
psi_t(:, 1) = psi; % Store the initial wavefunction in the time evolution array

% Propagate the wave packet
for t = 2:N_t
    psi = V_op .* psi; % Operate a half time step in real space
    psi_k = fftshift(fft(psi)); % Fourier transform the wavefunction into k-space
    psi_k = T_op .* psi_k; % Operate a full time step in k-space
    psi = ifft(ifftshift(psi_k)); % Inverse Fourier transform into real space
    psi = V_op .* psi; % Operate a half time step in real space

    psi = psi/sqrt(trapz(x, abs(psi).^2)); % Normalise the time-evolved wavefunction

    psi_t(:, t) = psi; % Store the time-evolved wavefunction in the time evolution array
    J(:, t) = -((1i * hbar)/(2 * m)) * ((conj(psi) .* (first_deriv * psi)) - ((first_deriv * conj(psi)) .* psi)); % Calculate probability current at each point along x
end

% ======================================================================================================================================
%%%%%%%%%% Plot the time evolution of the wave packet, probability density, and flux %%%%%%%%%%
% ======================================================================================================================================

figure; % Generate a figure

subplot(2, 2, 1) % Top Left subfigure
real_wavefunction = plot(x, real(psi_t(:, 1))); % Plot the real wavefunction
xlabel('$x$', 'Interpreter', 'latex'); % Label the x-axis
ylabel('$\mathrm{Re}(\psi(x, t))$', 'Interpreter', 'latex'); % Label the y-axis
ylim([min(real(psi_t(:))) max(real(psi_t(:)))]); % Set the y-limits for convenience
title('Real Component of the Wavefunction') % Add a title
grid on; % Add a grid to the plot

subplot(2, 2, 2) % Top right subfigure
imag_wavefunction = plot(x, imag(psi_t(:, 1))); % Plot the imaginary wavefunction
xlabel('$x$', 'Interpreter', 'latex'); % Label the x-axis
ylabel('$\mathrm{Im}(\psi(x, t))$', 'Interpreter', 'latex'); % Label the y-axis
ylim([min(imag(psi_t(:))) max(imag(psi_t(:)))]); % Set the y-limits for convenience
title('Imaginary Component of the Wavefunction') % Add a title
grid on; % Add a grid to the plot

subplot(2, 2, 3) % Bottom left subfigure
prob_density = plot(x, abs(psi_t(:, 1)).^2); % Plot the initial probability density
xlabel('$x$', 'Interpreter', 'latex'); % Label the x-axis
ylabel('$|\psi(x, t)|^2$', 'Interpreter', 'latex'); % Label the y-axis
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
    set(flux_plot, 'YData', J(:, n)); % Update the flux plot
    pause(0.05); % Pause to create an animation
    drawnow; % Update the figures and display immediately
end