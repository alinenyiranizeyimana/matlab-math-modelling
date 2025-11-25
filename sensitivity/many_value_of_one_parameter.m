%% ************************ Anonymous function ***************************

r_parameter = [0.5 0.45 0.4 0.35 0.3 0.25 0.21 0.201 0.2005];
K= 1000;
beta = 0.1;
delta = 0.2;
gamma = 0.02;
mu = 0.1;
 
y_0= [599 1 0];
t_span=[0,100];



seir_sens_anonymous = @(t,y,r)[...
    r*(y(1) + y(2)+ y(3))*(1-((y(1) + y(2)+ y(3))/K))-beta* y(1) * y(2) - delta*y(1);
    beta* y(1)*y(2) - gamma *y(2) - mu * y(2) -delta*y(2);
    gamma*y(2) - delta*y(3)
];

for i= 1:length(r_parameter)
    sir_ode= @(t,y) seir_sens_anonymous(t,y,r_parameter(i));
    [t,Y] = ode45(sir_ode, t_span, y_0);
end




