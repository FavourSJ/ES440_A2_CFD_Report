% Definition of the function and its exact second derivative
f = @(x) sin(x);
exact_second_derivative = @(x) -sin(x);

% Grid spacings for this example
grid_spacings = [0.1, 0.05, 0.01];
x_i = pi / 4; % Point of interest

% Arrays for errors
errors_chosen = zeros(size(grid_spacings));
errors_second_order = zeros(size(grid_spacings));
errors_fourth_order = zeros(size(grid_spacings));

% For Loop over each grid spacing that is being used (0.1, 0.05 and 0.01)
for idx = 1:length(grid_spacings)
    dx = grid_spacings(idx);
    
    % Neighboring points (±0.1, ±0.05 and ±0.01)
    x_i_minus_1 = x_i - dx;
    x_i_minus_2 = x_i - 2*dx;
    x_i_minus_3 = x_i - 3*dx;
    x_i_plus_1 = x_i + dx;
    x_i_plus_2 = x_i + 2*dx;
    
    % Chosen scheme
    numerical_chosen = (2*f(x_i) - 5*f(x_i_minus_1) + 4*f(x_i_minus_2) - f(x_i_minus_3)) / dx^2;
    errors_chosen(idx) = abs(numerical_chosen - exact_second_derivative(x_i));
    
    % Second-order scheme
    numerical_second_order = (f(x_i_plus_1) - 2*f(x_i) + f(x_i_minus_1)) / dx^2;
    errors_second_order(idx) = abs(numerical_second_order - exact_second_derivative(x_i));
    
    % Fourth-order scheme
    numerical_fourth_order = (-f(x_i_plus_2) + 16*f(x_i_plus_1) - 30*f(x_i) + 16*f(x_i_minus_1) - f(x_i_minus_2)) / (12*dx^2);
    errors_fourth_order(idx) = abs(numerical_fourth_order - exact_second_derivative(x_i));
end

% Results Table
fprintf('Grid Spacing (dx) | Chosen Error | 2nd-Order Error | 4th-Order Error\n');
fprintf('---------------------------------------------------------------\n');
for idx = 1:length(grid_spacings)
    fprintf('%15.5f | %12.6e | %15.6e | %15.6e\n', grid_spacings(idx), errors_chosen(idx), errors_second_order(idx), errors_fourth_order(idx));
end

% Plot errors against grid spacing on a log-log scale
figure(1);
loglog(grid_spacings, errors_chosen, '-o', 'LineWidth', 1.5, 'DisplayName', 'Chosen Scheme');
hold on;
loglog(grid_spacings, errors_second_order, '-s', 'LineWidth', 1.5, 'DisplayName', '2nd-Order Scheme');
loglog(grid_spacings, errors_fourth_order, '-d', 'LineWidth', 1.5, 'DisplayName', '4th-Order Scheme');
grid on;
xlabel('Grid Spacing (\Delta x)', 'FontSize', 12);
ylabel('Leading Error', 'FontSize', 12);
title('Error Comparison of Finite Difference Schemes', 'FontSize', 14);
legend('Location', 'best');
