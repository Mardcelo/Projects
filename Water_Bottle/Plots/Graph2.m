% Define the function
G = @(f, M) (M^2 * f.^4 + 4 * M * f.^3 - 6 * M * f.^2 + 4 * M * f + 1) ./ (1 + M * f).^2;
G_values = G(f_values, M);

% Find the minimum value of G(f)
[min_G, min_f] = min(G_values);

% Plot the function
figure;
plot(f_values, G_values, 'b-', 'LineWidth', 2);

% Mark the minimum point
hold on;
plot(f_values(min_f), min_G, 'ro', 'MarkerSize', 10);

% Label the axes and set the title
xlabel('f (filling fraction)');
ylabel('G(f)');
title('G(f) for M = 20');

% Add a legend
legend('G(f)', 'Minimum G(f)');

% Display the minimum value
disp(['Minimum G(f): ', num2str(min_G)]);
disp(['f at minimum G(f): ', num2str(f_values(min_f))]);
