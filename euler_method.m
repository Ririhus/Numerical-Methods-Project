% MATLAB Code: Euler's Method with Separate Figures for Each Step Size
clc;
clear;

% Given differential equation: dy/dx = -2y + exp(-4*x), y(0) = 1
f = @(x, y) -2*y + exp(-4*x); % Function for dy/dx

% Initial conditions
x0 = 0;
y0 = 1;
x_end = 2;

% Analytical solution
syms y(x)
Dy = diff(y, x);
eqn = Dy == -2*y + exp(-4*x);
cond = y(0) == 1;
y_analytical = dsolve(eqn, cond);
fprintf('Analytical Solution: y(x) = %s\n', y_analytical);

% Step sizes
h_values = [0.1, 0.05, 0.01];

% Convert analytical solution to MATLAB function
y_exact = matlabFunction(y_analytical);

% Loop through each step size
for j = 1:length(h_values)
    h = h_values(j);
    n = (x_end - x0) / h; % Number of steps
    x = x0:h:x_end; % x values
    y = zeros(1, length(x));
    y(1) = y0; % Initial condition
    
    % Apply Euler's method
    for i = 1:n
        y(i+1) = y(i) + h * f(x(i), y(i));
    end
    
    % Analytical solution at grid points
    y_exact_at_x = y_exact(x);
    
    % Plot the numerical and analytical solutions
    figure;
    plot(x, y_exact_at_x, 'k-', 'LineWidth', 2, 'DisplayName', 'Analytical Solution');
    hold on;
    plot(x, y, 'r--o', 'LineWidth', 1.5, 'DisplayName', sprintf("Euler's Method (h = %.2f)", h));
    xlabel('x');
    ylabel('y');
    title(sprintf("Comparison of Analytical and Numerical Solution (h = %.2f)", h));
    legend('Location', 'Best');
    grid on;
    hold off;
    
    % Calculate error
    error = abs(y - y_exact_at_x);
    
    % Calculate Mean Absolute Error (MAE)
    MAE = mean(error);
    
    % Display error for this step size
    fprintf('\nStep Size h = %.2f\n', h);
    fprintf('x\t\tNumerical y\tExact y\t\tError\n');
    fprintf('%.2f\t%.6f\t%.6f\t%.6f\n', [x; y; y_exact_at_x; error]);
    fprintf('Mean Absolute Error (MAE) = %.6f\n', MAE);
end