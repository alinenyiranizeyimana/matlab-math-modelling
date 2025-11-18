% ************* mu varying ***********************

y0 = [9900, 100, 0];
tspan = [0, 250];
beta = 0.5;
gamma = 0.2;
N = 10000;
mu_values = [0.00, 0.02, 0.06, 0.08];

%*********** no semi colon on last R**********

SIR_mu =@(t,y,my) [
    -beta* y(1)*y(2) + mu*(N - y(1));
    beta * y(1)*y(2) - (gamma + mu)*y(2);
    gamma * y(2) - mu * y(3)
]

for i = 1:length(mu_values)
    ano_mu_ode = @(t,y) SIR_mu(t,y,mu_values(i));
    [t,y] = ode45(ano_mu_ode, tspan, y0);

    figure(i)
    plot(t,y(:,1),'b',t,y(:,2),'r',t,y(:,3),'g')
    legend('Susceptible','Infected','Recovered');
end
