%% 4th Order Runge Kutta Method
f = @(x, y) -2*y + exp(-4*x); 

syms y(x)
ode = diff(y, x) == -2*y + exp(-4*x);
cond = y(0) == 1;
ySol(x) = dsolve(ode,cond);

% Display the analytical solution
disp('Analytical Solution:');
disp(ySol);

x0 = 0;
x_end = 2;
y0 = 1;
h_value = [0.1, 0.05, 0.01]; 

for h = h_value
    x = x0:h:x_end;
    N = length(x);

    y_RK4 = zeros(1, N);
    y_RK4(1) = y0;

    for i = 1:N-1
        k1 = h * f(x(i), y_RK4(i));
        k2 = h * f(x(i) + h/2, y_RK4(i) + k1/2);
        k3 = h * f(x(i) + h/2, y_RK4(i) + k2/2);
        k4 = h * f(x(i) + h, y_RK4(i) + k3);
        y_RK4(i+1) = y_RK4(i) + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
    end

    % Plot Results
    figure;
    g1 = fplot(ySol, [x0, x_end], 'LineWidth', 1.5);
    hold on;
    g2 = plot(x, y_RK4, '-o', 'LineWidth', 1.2);
    legend([g1, g2], {'Analytical Solution', sprintf('RK4 (h=%.2f)', h)}, 'Location', 'best');

    % Title and labels for each plot
    title(sprintf('4th Order Runge-Kutta Method vs Analytical Solution (h=%.2f)', h));
    xlabel('x');
    ylabel('y');
    grid on;
    hold off;

    fprintf('Result for step size h = %.2f:\n', h);
    fprintf('x\t\t y\n');
    for i = 1:N
        fprintf('%.2f\t %.6f\n', x(i), y_RK4(i));
    end
    
    % Calculate Error
    y_analytical = double(ySol(x)); % Compute analytical solution at each x
    error = abs(y_analytical - y_RK4); % Calculate the error

    % Calculate the Mean Error
    mean_error = mean(error); % Calculate the mean error

    % Display mean error for each step size
    fprintf('Mean error for h = %.2f: %.6f\n', h, mean_error);
    fprintf('\n');
end
