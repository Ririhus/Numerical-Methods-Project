% 3rd order taylor's series

% Define the differential equation
f = @(x, y) -2*y + exp(-4*x);

% Analytical solution
analytical_solution = @(x) -0.5*exp(-4*x) + 1.5*exp(-2*x);

% Initial conditions
x0 = 0;
y0 = 1;

% Step sizes to use
step_sizes = [0.1, 0.05, 0.01];
x_end = 2;

% Loop over different step sizes
for h = step_sizes
    % Interval and step size
    x = x0:h:x_end;
    
    % Initialize arrays to store solutions
    y = zeros(size(x));
    y(1) = y0;
    
    % Compute the third-order Taylor series solution
    for i = 1:(length(x) - 1)
        % First derivative
        f1 = f(x(i), y(i));
        
        % Second derivative
        f2 = -2*f1 - 4*exp(-4*x(i));
        
        % Third derivative
        f3 = -2*f2 + 16*exp(-4*x(i));
        
        % Third-order Taylor series approximation
        y(i + 1) = y(i) + h * f1 + (h^2 / 2) * f2 + (h^3 / 6) * f3;
    end
    
    % Compute the analytical solution
    y_analytical = analytical_solution(x);
    
    % Calculate the error at each step
    error = abs(y - y_analytical);

    % Calculate the mean absolute error
    mae = mean(error);

    % Display the error and mean absolute error
    disp(['Step size h = ', num2str(h)]);
    disp(' x         Numerical y      Analytical y     Error');
    disp([x', y', y_analytical', error']);
    disp(['Mean Absolute Error (MAE) for h = ', num2str(h), ': ', num2str(mae)]);
    
    % Plot the numerical solution vs analytical solution graph
    figure;
    hold on;
    plot(x, y, 'o', 'LineWidth', 1.5, 'DisplayName', ['Numerical Solution (h = ', num2str(h), ')']); % Numerical solution
    plot(x, y_analytical, '-r', 'LineWidth', 1.5, 'DisplayName', 'Analytical Solution'); % Analytical solution
    title(['Third-Order Taylor Series vs Analytical Solution (h = ', num2str(h), ')']);
    xlabel('x');
    ylabel('y');
    grid on;
    legend('Numerical Solution (3rd Order Taylor Series)','Analytical Solution')
    hold off;
end

