%% ************************ Anonymous function ***************************
% ******* To write a function you should put space between - and other ****



seir_sens_anonymous = @(t,y,r,beta)[
    r*(y(1) + y(2) + y(3))*(1 - (y(1) + y(2)+ y(3))/K) - beta* y(1) * y(2) - delta*y(1);
    beta* y(1)*y(2) - gamma *y(2) - mu * y(2) - delta*y(2);
    gamma*y(2) - delta*y(3)
];


%********** Prints 45 graphs for each beta compares with 9r*************
beta_values = [0.1, 0.2, 0.4, 0.6, 0.8]

for j=1:length(beta_values)
    beta=beta_values(j)
    r_parameter = [0.5 0.45 0.4 0.35 0.3 0.25 0.21 0.201 0.2005];
    K= 1000;
    delta = 0.2;
    gamma = 0.02;
    mu = 0.1;
    t_span=[0,100];


    for i= 1:length(r_parameter)
        sir_ode= @(t,y) seir_sens_anonymous(t,y,r_parameter(i),beta);
 
        [t,y] = ode45(sir_ode, t_span, y_0);
   
        figure(i)
   
        hold on
    
        plot(t, y(:,1),"mo-",MarkerSize= 5, MarkerFaceColor="r", ...
    MarkerEdgeColor='k')
   
        hold on
    
        plot(t, y(:,2),"k*--",MarkerSize= 3,MarkerEdgeColor='k')
   
        plot(t, y(:,3),"b*-",MarkerSize= 3,MarkerEdgeColor='g')
    
        hold off
    
        title("SIR MODEL USING ODE45 SOLVER")
   
        xlabel("time points")
    
        ylabel(" Number of pop in each comp. ")
    
        legend(["Susceptible","Infected","Recovered"])
    
        grid off 
    end

end





