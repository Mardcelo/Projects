% Define the functions for hCM/H and G(f)
hCM_over_H = @(M, f) (1/2) * ((1 + M * f.^2) ./ (1 + M * f));
G = @(f, M) (M^2 * f.^4 + 4 * M * f.^3 - 6 * M * f.^2 + 4 * M * f + 1) ./ (1 + M * f).^2;

% Set the value of M and a fine range of f values for higher accuracy
M = 20;
f_values = linspace(0, 1, 1000);  % Increase number of points to 1000 for better resolution

% Calculate hCM/H and G(f) for the given M and f values
hCM_over_H_values = hCM_over_H(M, f_values);
G_values = G(f_values, M);

% Create the figure and hold it for multiple plots
figure;

% Plot the function G(f) versus f (left y-axis)
yyaxis left;
plot(f_values, G_values, 'b-', 'LineWidth', 2);
ylabel('G(f)', 'Color', 'b');
ylim([0.3, 1.2]);  % Adjust limits based on the image
hold on;

% Mark the specific minimum point for G(f) at f = 0.41041 and G(f) = 0.359696
f_min_G = 0.41041;
G_min_value = 0.359696;
plot(f_min_G, G_min_value, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');

% Mark the minimum point for hCM/H
[min_hCM_over_H, min_f_index_hCM] = min(hCM_over_H_values);
plot(f_values(min_f_index_hCM), min_hCM_over_H, 'ro', 'MarkerSize', 10);

% Plot the function hCM/H versus f (right y-axis)
yyaxis right;
plot(f_values, hCM_over_H_values, 'r-', 'LineWidth', 2);
ylabel('h_{CM}/H', 'Color', 'r');
ylim([0.1, 0.5]);  % Adjust limits based on the image

% Label the x-axis
xlabel('f (filling fraction)');

% Add vertical lines at specific f values (dashed lines)
xline(0.1791, 'r--', 'LineWidth', 1.5);  % Dashed line at f = 0.1791
xline(0.4, 'b--', 'LineWidth', 1.5);     % Dashed line at f = 0.4

% Add shaded region between f = 0.1791 and f = 0.4
patch([0.1791, 0.4, 0.4, 0.1791], [0.1, 0.1, 0.5, 0.5], [0.8 0.8 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Add a title
title('Comparison of G(f) and h_{CM}/H when M = 20');

% Add a legend (optional, for more clarity)
legend('G(f)', 'Minimum G(f)', 'h_{CM}/H', 'Minimum h_{CM}/H', 'Location', 'Best');

% Display the minimum values with high precision
disp(['Minimum h_{CM}/H: ', num2str(min_hCM_over_H, 10)]);
disp(['Corresponding f value (hCM/H): ', num2str(f_values(min_f_index_hCM), 10)]);
disp(['Minimum G(f): ', num2str(G_min_value)]);
disp(['Corresponding f value (G(f)): ', num2str(f_min_G)]);

% Set the figure properties to maintain the aspect ratio and limits
axis tight; % Adjusts the axes to fit the data
grid on;   % Adds grid for better readability
