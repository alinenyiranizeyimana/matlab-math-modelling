%************ Local_linearization_example 2*********************

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



N = 100
[D1, x] = cheb(N)
t0 = 0; tf = 250;
t = ((tf -t0)/2) *x + ((tf + t0 )/2); 
D = 2/(tf - t0)* D1
Ide= eye(N+1)


S0 = 0.5;
E0= 0.2;
I0 = 0.1;
Q0 = 0.1;
R0 = 0.1;

delta = 0.5
beta1 = 1.05
beta2 = 0.005,
phi = 0.5,
q1 = 0.001;
epsilon = 0.00398; 
v1 = 0.0047876;
v2 = 0.000001231;
gamma = 0.085432;
kappa = 0.09871; 
chi = 0.1243;



S = S0* ones(N+1,1);
E = E0*ones(N+1,1);
I = I0*ones(N+1,1);
Q = Q0*ones(N+1,1);
R = R0*ones(N+1,1);


iteration = 50;

for i = 1:iteration
    %************** SUSCEPTIBLE***********
   
    A1 = D + diag(beta1*I) + diag(beta2*E) + phi*Ide
    R1 = delta *Ide
    A1(N+1,:)=0
    A1(N+1,N+1)=1
    R1(N+1) = S0
    S = A1 \ R1;

     %************** EXPOSED ***********

    A2 = D - diag(beta2*S) + (q1+epsilon+gamma+phi)
    R2 = beta1* S.*I
    A2(N+1, :)=0
    A2(N+1,N+1)=1 
    R2(N+1) = E0

    E = A2\R2;


    %************** INFECTED ***********

    A3 = D + diag(kappa + phi + v1) 
    R3 = gamma*E
    A3(N+1, :)=0
    A3(N+1,N+1)=1 
    R3(N+1) = I0

    I = A3\R3;


     %************** QUARANTINE *********

     A4 = D + diag(chi + phi + v2); 
     R4 = q1 * E
     A4(N+1,:) = 0;
     A4(N+1,N+1) =1;
     R4(N+1) = Q0

     Q = A4\R4;



     %************** RECOVERED *********

     A5 = D+ diag(phi); 
     R5= epsilon *E + kappa * I + chi*Q
     A5(N+1,:) = 0;
     A5(N+1,N+1) =1;
     R5(N+1) = R0

     R = A5\R5;
end


figure(1)
plot(t,S,'b')
hold on
plot(t, E,'r')
hold on
plot(t,I,'g')
hold on
plot(t,Q,'k')
hold on
plot(t, R,'y')
hold off
xlabel("Time in days")
%ylim([0,300])
%xlim([0,20])
ylabel("Number")
legend(["S","E","I","Q","R"])








