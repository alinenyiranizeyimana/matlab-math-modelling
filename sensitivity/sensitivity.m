%***********To make all variable not displayed put a semicolon 
%***********on the end

r_parameter = [0.5 0.45 0.4 0.35 0.3 0.25 0.21 0.201 0.2005];

K= 1000;
beta = 0.1;
delta = 0.2;
gamma = 0.02;
mu = 0.1;

 R_0=zeros(1,length(r_parameter));
 sens_index =zeros(1,length(r_parameter));
 

for i=1:length(r_parameter)
    r= r_parameter(i);
    R_0(i) = beta*K* (r -delta) / (r*(gamma+mu+delta));
    sens_index(i) = delta/(r-delta);
    disp(R_0);
    disp(sens_index);
end


%% code for SIR Model

function dydt = sir_sens(t,y)
%N= S+I+R

r=0.5;
K= 1000;
beta = 0.1;
delta = 0.2;
gamma = 0.02;
mu = 0.1;
s_0 = 599;
I_0 = 1;
R_0 = 0; 

N= y(1) + y(2)+ y(3)
dydt =[r*N*(1-(N/K))- beta* y(1) * y(2) - delta*y(1);beta* y(1)*y(2) - gamma *y(2) - mu * y(2) - delta*y(2);gamma*y(2) - delta*y(3)]
end


y_0= [599 1 0];
t_span=[0,100]

[t,y] = ode45(@sir_sens,t_span,y_0)

figure()
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






