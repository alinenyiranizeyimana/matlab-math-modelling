% Parameters (from Table 3)
b = 0.0784;
mu = 0.0545;
beta = 0.09091;
gamma = 0.125;
alpha = 0.14286;
sigma = 0.30;

% Time parameters
tmax = 10;
dt = 0.1;   % Time step size (larger than Euler to reduce iterations)
time = 0:dt:tmax;
n = length(time);

% Dimension of spatial discretization or system size
N = 50;      % Example: discretization points in D matrix

% Define matrix D (e.g. discrete differentiation operator or zero matrix if none)
D = zeros(N); % Replace with your actual matrix D

% Identity matrix
I_mat = eye(N);

% Initial populations - vectors of length N
S_r = ones(N,1) * 800/N;
E_r = ones(N,1) * 180/N;
I_r = ones(N,1) * 20/N;
R_r = zeros(N,1);

% Storage for solutions over time
S_all = zeros(N,n);
E_all = zeros(N,n);
I_all = zeros(N,n);
R_all = zeros(N,n);

S_all(:,1) = S_r;
E_all(:,1) = E_r;
I_all(:,1) = I_r;
R_all(:,1) = R_r;

% Relaxation parameters
omega = 1.2;  % Over-relaxation factor
tol = 1e-6;   % Convergence tolerance
max_iter = 100;

for k = 2:n
    % Compute right hand sides for linear systems
    R1 = b - beta .* (S_r .* I_r);
    R2 = beta .* (S_r .* I_r);
    R3 = alpha .* E_r;
    R4 = gamma .* I_r + sigma .* E_r;
    
    % Define coefficient matrices
    A1 = D + mu * I_mat;
    A2 = D + (mu + alpha + sigma) * I_mat;
    A3 = D + (mu + gamma) * I_mat;
    A4 = D + mu * I_mat;
    
    % Solve for S_{r+1} by relaxation
    S_next = relaxation_solver(A1, R1, S_r, omega, tol, max_iter);
    % Solve for E_{r+1} by relaxation
    E_next = relaxation_solver(A2, R2, E_r, omega, tol, max_iter);
    % Solve for I_{r+1} by relaxation
    I_next = relaxation_solver(A3, R3, I_r, omega, tol, max_iter);
    % Solve for R_{r+1} by relaxation
    R_next = relaxation_solver(A4, R4, R_r, omega, tol, max_iter);
    
    % Update for next iteration
    S_r = S_next;
    E_r = E_next;
    I_r = I_next;
    R_r = R_next;
    
    % Store results
    S_all(:,k) = S_r;
    E_all(:,k) = E_r;
    I_all(:,k) = I_r;
    R_all(:,k) = R_r;
end

% Define relaxation solver function
function x = relaxation_solver(A, b, x0, omega, tol, max_iter)
    n = length(b);
    x = x0;
    for iter = 1:max_iter
        x_old = x;
        for i = 1:n
            sigma = A(i,1:i-1)*x(1:i-1) + A(i,i+1:n)*x_old(i+1:n);
            x(i) = (1 - omega)*x_old(i) + omega*(b(i) - sigma)/A(i,i);
        end
        if norm(x - x_old, inf) < tol
            break;
        end
    end
end

% Plot example: sum over spatial domain to get total population curves
figure;
plot(time, sum(S_all), 'b', 'LineWidth', 2); hold on;
plot(time, sum(E_all), 'r', 'LineWidth', 2);
plot(time, sum(I_all), 'g', 'LineWidth', 2);
plot(time, sum(R_all), 'm', 'LineWidth', 2);
xlabel('Time');
ylabel('Population');
title('SEIR Model Dynamics via Relaxation Linear Systems');
legend('Susceptible','Exposed','Infective','Recovered');
grid on;
