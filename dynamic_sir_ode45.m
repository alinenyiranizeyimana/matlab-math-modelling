%% *********** QUESTION  *****************
figure(1)
image =imread("work_home.png")
imshow(image)

%% ************ ODE45_Solution*****************


%% ****************** Beta varying*************
y0 = [9900, 100, 0];
tspan = [0, 250];
betas = [0.5, 0.45, 0.4, 0.35];
gamma = 0.2;
N = 10000;
mu = 0.02;

SIR = @(t,y,beta) [
    -beta*y(1)*y(2) + mu*(N - y(1));
     beta*y(1)*y(2) - (gamma + mu)*y(2);
     gamma*y(2) - mu*y(3)
];

for i = 1:length(betas)
    
    anoy_ode = @(t,y) SIR(t,y,betas(i));
    [t,Y] = ode45(anoy_ode, tspan, y0);

    figure(i+1)
    plot(t,Y(:,1),'b',t,Y(:,2),'r',t,Y(:,3),'g')
    legend('Susceptible','Infected','Recovered');

    title(['SIR MODEL, \beta = ', num2str(betas(i))], 'Interpreter','latex')
    xlabel('Time');
    ylabel('Population');
end

