% ======================================================================================================================================
%%%%%%%%%% Preamble %%%%%%%%%%
% ======================================================================================================================================

% Part Ia - Particle in a Box Wave Packet Simulation
% Verification of the numerical calculation of particle in a box eigenstates

% Author: Max L Butterworth
% MSc in Theoretical and Computational Chemistry Project
% University of Oxford

% ======================================================================================================================================
%%%%%%%%%% Define constants and variables %%%%%%%%%%
% ======================================================================================================================================

% Natural units have been adopted throughout
L = 1;% Length of the 1D box
m = 1; % Mass
h = 1; % Planck's constant in J s
hbar = 1; % Definition of h bar

N_steps = 1000; % Number of discretisation points

basis_funcs_indices = [1, 2, 3, 4]; % Create an array of the indices of PIB_eigenstates_norm that form the superposition
N_PIB_eigenfuncs = max(basis_funcs_indices); % The number of basis functions in the wave packet superposition

% ======================================================================================================================================
%%%%%%%%%% Discretise the spatial domain, x, and time domain, t %%%%%%%%%%
% ======================================================================================================================================

x = linspace(0, L, N_steps); % Define the domain of the infinite potential well
dx = x(2) - x(1); % Calculate the spatial step size

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
    
    % Set the wavefunction to zero at the boundaries
    PIB_eigenstates_norm(1, i) = 0;
    PIB_eigenstates_norm(N_steps, i) = 0;
end

% ======================================================================================================================================
%%%%%%%%%% Generate the analytical solutions to the Schrödinger equation for the PIB %%%%%%%%%%
% ======================================================================================================================================

PIB_eigenstates_analytical = zeros(length(x), max(basis_funcs_indices)); % Define an array to store the analytical PIB eigenstates

% Create an array of analytical eigenfunctions of the PIB Hamiltonian
for n = 1:max(basis_funcs_indices)
    PIB_eigenstates_analytical(:, n) = sqrt(2/L) * sin((n * pi * x)/L).';
end

% ======================================================================================================================================
%%%%%%%%%% Calculate the absolute error of the numerical solutions %%%%%%%%%%
% ======================================================================================================================================

absolute_error = zeros(length(x), max(basis_funcs_indices)); % Initialise an array to store the errors

for j = 1:length(basis_funcs_indices)
    absolute_error(:, j) = abs(PIB_eigenstates_norm(:, j)) - abs(PIB_eigenstates_analytical(:, j));
end

% ======================================================================================================================================
%%%%%%%%%% Plot the analytical and numerical eigenstates for comparison %%%%%%%%%%
% ======================================================================================================================================

figure; % Generate a figure

% Ground state comparison
subplot(2, 2, 1) % Top Left subfigure
plot(x, -real(PIB_eigenstates_norm(:, 1))); % Plot the numerically determined real wavefunction
hold on
plot(x, PIB_eigenstates_analytical(:, 1)); % Plot the analytical wavefunction
plot(x, absolute_error(:, 1)); % Plot the absolute error
hold off
xlabel('$x$', 'Interpreter', 'latex'); % Label the x-axis
ylabel('$\mathrm{Re}(\psi(x, t))$', 'Interpreter', 'latex'); % Label the y-axis
title('Ground State') % Add a title
legend('Numerical Solution', 'Analytical Solution', 'Absolute Error') % Add a legend and label the data
grid on; % Add a grid to the plot

% First excited state comparison
subplot(2, 2, 2) % Top Left subfigure
plot(x, real(PIB_eigenstates_norm(:, 2))); % Plot the numerically determined real wavefunction
hold on
plot(x, PIB_eigenstates_analytical(:, 2)); % Plot the analytical wavefunction
plot(x, absolute_error(:, 2)); % Plot the absolute error
hold off
xlabel('$x$', 'Interpreter', 'latex'); % Label the x-axis
ylabel('$\mathrm{Re}(\psi(x, t))$', 'Interpreter', 'latex'); % Label the y-axis
title('First Excited State') % Add a title
legend('Numerical Solution', 'Analytical Solution', 'Absolute Error') % Add a legend and label the data
grid on; % Add a grid to the plot

% Second excited state comparison
subplot(2, 2, 3) % Bottom Right subfigure
plot(x, real(PIB_eigenstates_norm(:, 3))); % Plot the numerically determined real wavefunction
hold on
plot(x, PIB_eigenstates_analytical(:, 3)); % Plot the analytical wavefunction
plot(x, absolute_error(:, 3)); % Plot the absolute error
hold off
xlabel('$x$', 'Interpreter', 'latex'); % Label the x-axis
ylabel('$\mathrm{Re}(\psi(x, t))$', 'Interpreter', 'latex'); % Label the y-axis
title('Second Excited State') % Add a title
legend('Numerical Solution', 'Analytical Solution', 'Absolute Error') % Add a legend and label the data
grid on; % Add a grid to the plot

% Third excited state comparison
subplot(2, 2, 4) % Bottom left subfigure
plot(x, real(PIB_eigenstates_norm(:, 4))); % Plot the numerically determined real wavefunction
hold on
plot(x, PIB_eigenstates_analytical(:, 4)); % Plot the analytical wavefunction
plot(x, absolute_error(:, 4)); % Plot the absolute error
hold off
xlabel('$x$', 'Interpreter', 'latex'); % Label the x-axis
ylabel('$\mathrm{Re}(\psi(x, t))$', 'Interpreter', 'latex'); % Label the y-axis
title('Third Excited State') % Add a title
legend('Numerical Solution', 'Analytical Solution', 'Absolute Error') % Add a legend and label the data
grid on; % Add a grid to the plot