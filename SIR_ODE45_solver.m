%% Anonymous function
%**** take any input of x
y= @(x)(x^2)
y(2)

y_two = @(x,y) (x+3*y)

y_two(1,0)


%% Equation to solve

% ***** Figure numbering to be able to see all at same times**********
figure(1)
image =imread("sir.png")
imshow(image)


%% Solution 

function return_dydt = SIR(t,y)

gamma = 0.1;
beta = 0.0004;

%********** Susceptible, Infected and Recovered functions ***********
%********** y(1) substitute S previous and I previous. no need ****** 
%********** we have to mention only their corresponding places only
%           y(1) = S, y(2) = I , y(3) = R

return_dydt = [-beta *y(1) *y(2) ; beta*y(1)*y(2) - gamma*y(2); gamma*y(2)]   

end

tspan = [0,250]

%*********** Same order as (equation in return_dydt list) not anony.
y0 = [990;10;0]

[t,y] = ode45(@SIR ,tspan, y0)


%*********** Solution are form of col1(S), col2(I) and Col3(R)

%            Susceptible length 145 
%            because it choose its time points to meet error tolerance



figure(2)
hold on
length(y(:,1))
plot(t, y(:,1),"mo-",MarkerSize= 3, MarkerFaceColor="r", ...
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

%*********** no png on filename **************
%print("SIR_ODE45",'-dpng','-r500')
%gcf get current figure

saveas(gcf,"gcf_SIR_ODE45.png")

