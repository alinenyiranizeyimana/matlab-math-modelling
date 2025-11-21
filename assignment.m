% SEIR Model Simulation to Reproduce Figures 2 and 3
% Parameters (from Table 3 of the article)
b = 0.0784;
mu = 0.0545;
beta = 0.09091;
gamma = 0.125;
alpha = 0.14286;

% Define the SEIR ODE system
seir_ode = @(t, y, sigma) [
    b - beta * y(1) * y(3) - mu * y(1);                      % dS/dt
    beta * y(1) * y(3) - (mu + alpha + sigma) * y(2);        % dE/dt
    alpha * y(2) - (mu + gamma) * y(3);                      % dI/dt
    gamma * y(3) + sigma * y(2) - mu * y(4)                  % dR/dt
];

% Initial Conditions (Assumed)
S0 = 800;
E0 = 180;
I0 = 20;
R0 = 0;
y0 = [S0; E0; I0; R0];

% Time span (in years, adjusted to show dynamics)
tspan = [0 10];

% Case 1: sigma = 0.00 (Figure 2)
sigma = 0.30;
[t1, Y1] = ode45(@(t, y) seir_ode(t, y, sigma), tspan, y0);

% Case 2: sigma = 0.30 (Figure 3)
sigma = 0.6;
[t2, Y2] = ode45(@(t, y) seir_ode(t, y, sigma), tspan, y0);



figure(1);
plot(t1, Y1(:,1), 'b', 'LineWidth', 2); hold on;   % Susceptible
plot(t1, Y1(:,2), 'r', 'LineWidth', 2);            % Exposed
plot(t1, Y1(:,3), 'g', 'LineWidth', 2);            % Infective
plot(t1, Y1(:,4), 'm', 'LineWidth', 2);    
% Recovered
xlim([0,10])
ylim([0,800])
yticks(0:200:800);
xticks(0:2:10);
title('SEIR Model Dynamics for \sigma = 0.00');
xlabel('Time');
ylabel('Population');

legend('Susceptible','Exposed','Infective','Recovered','Location','best');
grid on;


figure(2);
plot(t2, Y2(:,1), 'b--', 'LineWidth', 2); hold on;   % Susceptible
plot(t2, Y2(:,2), 'r--', 'LineWidth', 2);            % Exposed
plot(t2, Y2(:,3), 'g--', 'LineWidth', 2);            % Infective
plot(t2, Y2(:,4), 'm--', 'LineWidth', 2);            % Recovered
xlim([0,10])
yticks(0:200:1000);
xticks(0:2:10);
ylim([0,1000])

title('SEIR Model Dynamics for \sigma = 0.30');
xlabel('Time');
ylabel('Population');
legend('Susceptible','Exposed','Infective','Recovered','Location','best');
grid on;

length(Y2)
 