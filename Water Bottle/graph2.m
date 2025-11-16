% G(f) as a function of M and f
G = @(f, M) (M^2 .* f.^4 + 4*M .* f.^3 - 6*M .* f.^2 + 4*M .* f + 1) ./ (1 + M .* f).^2;

% Parameter and sampling
M        = 20;
f_values = linspace(0, 1, 2000);

% Values for plotting
G_values = G(f_values, M);

% Continuous minimum
f_min   = fminbnd(@(f) G(f, M), 0, 1);
G_min   = G(f_min, M);

% Plot
figure;
plot(f_values, G_values, 'b-', 'LineWidth', 2); hold on;
plot(f_min, G_min, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');

xlabel('f (filling fraction)');
ylabel('G(f)');
title(sprintf('G(f) for M = %d', M));
legend('G(f)', 'min G(f)', 'Location', 'best');
grid on; box on;
set(gca, 'FontSize', 12);

% Print 
fprintf('Minimum G(f): %.12f\n', G_min);
fprintf('f at minimum G(f): %.12f\n', f_min);
