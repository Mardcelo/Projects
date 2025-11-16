% hCM/H as a function of M and uhh f
hCM_over_H = @(M, f) 0.5 * ((1 + M .* f.^2) ./ (1 + M .* f));

% Parameter 
M        = 20;
f_values = linspace(0, 1, 2000);

% Values
h_values = hCM_over_H(M, f_values);

%  f* = (sqrt(1+M) - 1)/M
f_min   = (sqrt(1 + M) - 1) / M;
h_min   = hCM_over_H(M, f_min);

% Plot
figure;
plot(f_values, h_values, 'r-', 'LineWidth', 2); hold on;
plot(f_min, h_min, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');

xlabel('f (filling fraction)');
ylabel('h_{CM}/H', 'Color', 'r');
title(sprintf('h_{CM}/H when M = %d', M));
ylim([0.1, 0.5]);
grid on; box on;
set(gca, 'FontSize', 12);

% Print 
fprintf('Minimum h_{CM}/H: %.12f\n', h_min);
fprintf('Corresponding f value: %.12f\n', f_min);
