x = linspace(0, 2*pi, 1000); % Assume L = 2pi for the PIB
y = sin(3*x/2);

figure;
plot(x, y, 'k', 'LineWidth', 2); % 'k' = black line
axis equal;
axis off;

% Save as EPS (vector format)
print('-depsc2', 'sine_wave_PIB_n_3.eps');