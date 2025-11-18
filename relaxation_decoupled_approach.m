% =================== powering with 2 by the values of the vector
x = 2.^(0:3)

% =================== gives the same value puted rows and columns
one = ones(2,1)

N = 0

% if on online   .................................................
if N== 0, disp("D= 0"); end


% ******* for loop
%for index = startvalue:step:endvalue

%end

for i= 1:4
    disp(i)
end




%% cheb function

function [D ,x] = cheb (N)
% ======= taking values of x as symmetric distribution==================
% ======= Chebyvs gives accurate answer than euler(same steps)==========
% ======= Having eqtn f(t) = e^t, taking N=4 we can copute t and then after
% =======use lagrange to show the exact graph
% ======= cos(i*pi / N)

if N==0
    D= 0;
    x =1; 
    return
end

x = cos(pi * (0:N)/N)' ;


% ****** to find L'(t) we use bary centre interpolation=================
%        D_ij = (ci/cj)* (1/dx), dx = xi - xj = x - x'

% ****** end boundary 2 other 1 as columns =============================
% ******element by elemt multiplication(.*) ===========================
%         c: weight 

c = [2;ones(N-1,1);2].* (-1).^(0:N)';

% ******* Making a matrix *********************
x =repmat(x, 1, N+1)

% *******  D_ij = (ci/cj)* (1/dx), dx = xi - xj = x - x'

dx = x - x' 

% ********scaling dx to remove zero ************
D = (c* (1./c)') ./ (dx +(eye(N+1)))

% ==========   D_ii = - (D0,1+ D0,2) in notes ...............

D = D - diag(sum(D'));
end





%% SIR MODEL USING RELAXATION APPROACH

% dS/dt =  - Beta I S
% dI/dt = Beta I S - \gamma I
% dR/dt = \gamma I 

% S' = D S

% ================RELAXATION=======================

N = 100
[D1, x] = cheb(N)

% ===============Linear transformation (scaling)===

t0 = 0; tf = 250;

% ========= since (t0 -tf)/2)  scalar no need for ./
t = ((tf -t0)/2) *x + ((tf + t0 )/2); 

% ========= Since D [-1,1] SCale chain rule=========
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
    %***** D S = R1 , S = D^-1 R1
    A1 = D
    R1 = - beta * S.* I
    A1(N+1,:)=0
    A1(N+1,N+1)=1
    R1(N+1) = S0
    %***** inv(A1)
    S = A1 \ R1;

    %************** INFECTED ***********
    %***** D I = R2 - gamma I(current) since linear

    A2 = D + gamma*Ide
    R2 = beta * S.* I
    A2(N+1, :)=0
    A2(N+1,N+1)=1 
    R2(N+1)= I0;

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



%% Another equations

%% ***** SIMULATION OF COVID-19 USING MULTIDOMAIN SPECTRAL RELAXATION TECHNIQUE *********

image = imread("equations.png")
imshow(image)


%% **** simulation codes and graph of COVID-19
% Parameters
alpha = 0.5; 
beta1 = 1.05;
beta2 = 0.005; 
gamma = 0.5; 
q1 = 0.001;
delta = 0.00398;
lambda1 = 0.0047876; 
lambda2 = 0.000001231;
omega = 0.085432;
K = 0.09871; 
mu = 0.1243;

N = 6; % number of spectral points
iteration = 50;

% Initial conditions
S0 = 0.5; E0 = 0.2; I0 = 0.1; Q0 = 0.1; R0 = 0.1;

% Compute differentiation matrix and grid points
[D, x] = cheb(N);

% Identity matrix
Ide = eye(N+1);

% Initialize solution vectors
S = S0 * ones(N+1,1);
E = E0 * ones(N+1,1);
I = I0 * ones(N+1,1);
Q = Q0 * ones(N+1,1);
R = R0 * ones(N+1,1);

for i = 1:iteration
    % Susceptible update
    A1 = D;
    R1 = - beta1 * S .* I - beta2 * S .* E - delta * S + omega; % replace with model's S eqn
    A1(N+1,:) = 0; A1(N+1,N+1) = 1; R1(N+1) = S0;
    S = A1 \ R1;

    % Exposed update
    A2 = D + (q1 + lambda1 + gamma + delta) * Ide;
    R2 = beta1 * S .* I + beta2 * S .* E;
    A2(N+1,:) = 0; A2(N+1,N+1) = 1; R2(N+1) = E0;
    E = A2 \ R2;

    % Infected update
    A3 = D + (K + delta + lambda2) * Ide;
    R3 = gamma * E;
    A3(N+1,:) = 0; A3(N+1,N+1) = 1; R3(N+1) = I0;
    I = A3 \ R3;

    % Quarantined update
    A4 = D + (mu + delta) * Ide;
    R4 = q1 * E;
    A4(N+1,:) = 0; A4(N+1,N+1) = 1; R4(N+1) = Q0;
    Q = A4 \ R4;

    % Recovered update
    A5 = D + delta * Ide;
    R5 = lambda1 * E + K * I + mu * Q;
    A5(N+1,:) = 0; A5(N+1,N+1) = 1; R5(N+1) = R0;
    R = A5 \ R5;
end

% Linear transform for time range [0,75]
t0 = 0; tf = 75;
t = ((tf - t0)/2)*x + (tf + t0)/2;

% Plot results
figure
hold on
plot(t, S, 'b', 'DisplayName', 'Susceptible');
plot(t, E, 'k', 'DisplayName', 'Exposed');
plot(t, I, 'r', 'DisplayName', 'Infected');
plot(t, Q, 'm', 'DisplayName', 'Quarantined');
plot(t, R, 'g', 'DisplayName', 'Recovered');
xlabel('Time (days)');
ylabel('Population fraction');
legend;
grid on;
hold off;
