m = 1;
omega = 1;
hbar = 1;
n = 3;
x = linspace(-5, 5, 1000);

y = sqrt(1/((2^n) * factorial(n))) * ((m * omega)/(pi * hbar))^0.25 * exp(-(m * omega * x.^2)/(2 * hbar)) .* hermiteH(n, (sqrt((m * omega)/hbar)) * x);
% y = 0.5 * m * omega^2 * x.^2;

figure;
plot(x, y, 'k', 'LineWidth', 2); % 'k' = black line
axis equal;
axis off;

% Save as EPS (vector format)
print('-depsc2', 'QHO_Eigenfunction_n_3.eps');
% print('-depsc2', 'Quadratic_Potential.eps');