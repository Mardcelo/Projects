% Define the functions for hCM/H and G(f)
hCM_over_H = @(M, f) 0.5 * ((1 + M .* f.^2) ./ (1 + M .* f));
G          = @(f, M) (M^2 .* f.^4 + 4*M .* f.^3 - 6*M .* f.^2 + 4*M .* f + 1) ./ (1 + M .* f).^2;

% Parameter and sampling
M        = 20;
f_values = linspace(0, 1, 2000);  

hCM_over_H_values = hCM_over_H(M, f_values);
G_values          = G(f_values, M);

% Numerical minima
[min_hCM_over_H, min_f_index_hCM] = min(hCM_over_H_values);
f_min_hCM = f_values(min_f_index_hCM);

f_min_G   = fminbnd(@(f) G(f, M), 0, 1);  % refine G minimum
G_min_val = G(f_min_G, M);
figure;

% G(f) on left axis
yyaxis left;
hG = plot(f_values, G_values, 'b-', 'LineWidth', 2); hold on;
ylabel('G(f)', 'Color', 'b');
ylim([0.3, 1.2]);

% Mark minimum of G(f)
hGmin = plot(f_min_G, G_min_val, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');

% hCM/H on right axis
yyaxis right;
hh = plot(f_values, hCM_over_H_values, 'r-', 'LineWidth', 2);
ylabel('h_{CM}/H', 'Color', 'r');
ylim([0.1, 0.5]);

% Mark minimum of hCM/H
hHmin = plot(f_min_hCM, min_hCM_over_H, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');

% X-axis label
xlabel('f (filling fraction)');

% Vertical lines at specific f values
f1 = 0.1791;
f2 = 0.4;
xline(f1, 'r--', 'LineWidth', 1.5);
xline(f2, 'b--', 'LineWidth', 1.5);

% Shaded region between f1 and f2 
y_limits = ylim;
hPatch = patch([f1, f2, f2, f1], [y_limits(1), y_limits(1), y_limits(2), y_limits(2)], ...
               [0.8 0.8 0.8], 'FaceAlpha', 0.15, 'EdgeColor', 'none');

% Title and legend
title(sprintf('Comparison of G(f) and h_{CM}/H when M = %d', M));
legend([hG, hGmin, hh, hHmin, hPatch], ...
       {'G(f)', 'min G(f)', 'h_{CM}/H', 'min h_{CM}/H', sprintf('%.4f \x2264 f \x2264 %.4f', f1, f2)}, ...
       'Location', 'best');

grid on;
box on;
set(gca, 'FontSize', 12);

fprintf('Minimum h_{CM}/H: %.10f\n', min_hCM_over_H);
fprintf('Corresponding f value (hCM/H): %.10f\n', f_min_hCM);
fprintf('Minimum G(f): %.10f\n', G_min_val);
fprintf('Corresponding f value (G(f)): %.10f\n', f_min_G);
