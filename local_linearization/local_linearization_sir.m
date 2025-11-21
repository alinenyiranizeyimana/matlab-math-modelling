%%**************Local linearization approach*************************


%% cheb function

function [D ,x] = cheb (N)
if N==0
    D= 0;
    x =1; 
    return
end

x = cos(pi * (0:N)/N)' ;
c = [2;ones(N-1,1);2].* (-1).^(0:N)';
x =repmat(x, 1, N+1)
dx = x - x' 

% ***scaling dx to remove zero******
D = (c* (1./c)') ./ (dx +(eye(N+1)))
D = D - diag(sum(D'));
end


%% SIR MODEL USING RELAXATION APPROACH

% dS/dt =  - Beta I S
% dI/dt = Beta I S - \gamma I
% dR/dt = \gamma I 
% S' = D S

%**********************************************************
N = 100
[D1, x] = cheb(N)
t0 = 0; tf = 250;
t = ((tf -t0)/2) *x + ((tf + t0 )/2); 
D = 2/(tf - t0)* D1
Ide= eye(N+1)



%************** Initial Conditions **************

S0 = 990; I0 = 10 ; R0 = 0;
gamma = 0.1; beta = 0.0004


%********* making them as a column

S = S0* ones(N+1,1);
I = I0*ones(N+1,1);
R = R0*ones(N+1,1);

iteration = 50;

for i = 1:iteration
    %************** SUSCEPTIBLE***********
   

    A1 = D + diag(beta*I)
    R1 = zeros(N+1,1)
    A1(N+1,:)=0
    A1(N+1,N+1)=1
    R1(N+1) = S0
    %***** inv(A1)
    S = A1 \ R1;

    %************** INFECTED ***********

    A2 = D - diag(beta*S) + gamma*Ide
    R2 = zeros(N+1,1)
    A2(N+1, :)=0
    A2(N+1,N+1)=1 
    R2(N+1) = I0

    I = A2\R2;

     %************** RECOVERED *********

     A3 = D; 
     R3 = gamma * I
     A3(N+1,:) = 0;
     A3(N+1,N+1) =1;
     R3(N+1) = R0

     R = A3\R3;
end


hold on
plot(t,S,'b')
plot(t, I,'r')
plot(t,R,'g')
xlabel("Time in days")
ylabel("Number")
legend("Susceptible","Infected","Recovered")
hold off



 %******** Steps:

 % where there is previously computed value make the S E I REACHED
 % linear