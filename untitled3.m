%% =================== SEIR Model via Chebyshev Spectral Relaxation with Jacobian ===================
clc; clear; close all;

%% ---------------- Parameters ----------------
b     = 0.0784; 
mu    = 0.0545; 
beta  = 0.09091; 
gamma = 0.125; 
alpha = 0.14286;
sigma = 0.3;          % single sigma

%% ---------------- Chebyshev setup ----------------
N = 30;                % number of spectral points
iteration = 20;        % relaxation iterations

[D1, x] = cheb(N);

% Map Chebyshev nodes to [0,10] days
t0 = 0; tf = 10;
t = ((tf-t0)/2)*x + ((tf+t0)/2);
D = 2/(tf-t0) * D1;
Ide = eye(N+1);

%% ---------------- Initial conditions ----------------
S0 = 800; E0 = 200; I0 = 10; R0 = 0;

S = S0 * ones(N+1,1);
E = E0 * ones(N+1,1);
I = I0 * ones(N+1,1);
R = R0 * ones(N+1,1);

%% ---------------- Relaxation iterations with Jacobian ----------------
for k = 1:iteration
    %% Susceptible
    % Jacobian of RHS wrt S
    JS = mu * Ide; 
    A1 = D - JS; 
    R1 = b + mu * S - beta*S.*I;
    A1(N+1,:) = 0; A1(N+1,N+1) = 1; R1(N+1) = S0;
    S = A1 \ R1;

    %% Exposed
    JE = -(mu + alpha + sigma) * Ide;   % Jacobian for E
    A2 = D - JE;
    R2 = beta*S.*I - (mu + alpha + sigma)*E;
    A2(N+1,:) = 0; A2(N+1,N+1) = 1; R2(N+1) = E0;
    E = A2 \ R2;

    %% Infective
    JI = -(mu + gamma) * Ide;
    A3 = D - JI;
    R3 = alpha*E - (mu + gamma)*I;
    A3(N+1,:) = 0; A3(N+1,N+1) = 1; R3(N+1) = I0;
    I = A3 \ R3;

    %% Recovered
    JR = -mu * Ide;
    A4 = D - JR;
    R4 = gamma*I + sigma*E - mu*R;
    A4(N+1,:) = 0; A4(N+1,N+1) = 1; R4(N+1) = R0;
    R = A4 \ R4;
end

%% ---------------- Plot ----------------
figure;
plot(t, S, 'b', 'LineWidth',2); hold on;
plot(t, E, 'r', 'LineWidth',2);
plot(t, I, 'g', 'LineWidth',2);
plot(t, R, 'm', 'LineWidth',2);
xlim([0 10]); ylim([0 1000]);
xlabel('Time (days)'); ylabel('Population');
title('SEIR Model Dynamics (\sigma = 0.3) with Jacobian stabilization');
legend('Susceptible','Exposed','Infective','Recovered','Location','best');
grid on;

%% ---------------- Chebyshev Differentiation Matrix ----------------
function [D, x] = cheb(N)
    if N == 0
        D = 0; x = 1; return;
    end
    x = cos(pi*(0:N)/N)'; 
    c = [2; ones(N-1,1); 2] .* (-1).^(0:N)'; 
    X = repmat(x,1,N+1); 
    dX = X - X'; 
    D = (c * (1./c)') ./ (dX + eye(N+1)); 
    D = D - diag(sum(D'));
end
