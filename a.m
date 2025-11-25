% Parameters (from Table 3 of the article)
b = 0.0784;
mu = 0.0545;
beta = 0.09091;
gamma = 0.125;
alpha = 0.14286;

% Time parameters
tmax = 10;
dt = 0.01;  % smaller dt for better accuracy
time = 0:dt:tmax;
n = length(time);

% Initial conditions
S0 = 800;
E0 = 180;
I0 = 20;
R0 = 0;

% Function to simulate SEIR with Euler for given sigma
function Y = simulate_SEIR_euler(b, mu, beta, gamma, alpha, sigma, S0, E0, I0, R0, dt, n)
    Y = zeros(n,4);
    Y(1,:) = [S0, E0, I0, R0];
    for i = 2:n
        S = Y(i-1,1);
        E = Y(i-1,2);
        I = Y(i-1,3);
        R = Y(i-1,4);
        dS = b - beta * S * I - mu * S;
        dE = beta * S * I - (mu + alpha + sigma) * E;
        dI = alpha * E - (mu + gamma) * I;
        dR = gamma * I + sigma * E - mu * R;
        Y(i,1) = S + dt * dS;
        Y(i,2) = E + dt * dE;
        Y(i,3) = I + dt * dI;
        Y(i,4) = R + dt * dR;
    end
end

% Simulate for sigma = 0.00 (Figure 2)
sigma1 = 0.00;
Y1 = simulate_SEIR_euler(b, mu, beta, gamma, alpha, sigma1, S0, E0, I0, R0, dt, n);

% Simulate for sigma = 0.30 (Figure 3)
sigma2 = 0.30;
Y2 = simulate_SEIR_euler(b, mu, beta, gamma, alpha, sigma2, S0, E0, I0, R0, dt, n);

% Plot results for sigma = 0.00
figure(1);
plot(time, Y1(:,1), 'b', 'LineWidth', 2); hold on;
plot(time, Y1(:,2), 'r', 'LineWidth', 2);
plot(time, Y1(:,3), 'g', 'LineWidth', 2);
plot(time, Y1(:,4), 'm', 'LineWidth', 2);
xlim([0,10]);
ylim([0,800]);
yticks(0:200:800);
xticks(0:2:10);
title('SEIR Model Dynamics for \sigma = 0.00');
xlabel('Time');
ylabel('Population');
legend('Susceptible','Exposed','Infective','Recovered','Location','best');
grid on;

% Plot results for sigma = 0.30
figure(2);
plot(time, Y2(:,1), 'b--', 'LineWidth', 2); hold on;
plot(time, Y2(:,2), 'r--', 'LineWidth', 2);
plot(time, Y2(:,3), 'g--', 'LineWidth', 2);
plot(time, Y2(:,4), 'm--', 'LineWidth', 2);
xlim([0,10]);
ylim([0,1000]);
yticks(0:200:1000);
xticks(0:2:10);
title('SEIR Model Dynamics for \sigma = 0.30');
xlabel('Time');
ylabel('Population');
legend('Susceptible','Exposed','Infective','Recovered','Location','best');
grid on;
