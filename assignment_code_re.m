%% ************** QUESTION 1.C*************************


%% ***************Cheb function************************

function [D, x] = chebfunction(N)
    if N == 0
        D = 0; x = 1; return;
    end
    x = cos(pi*(0:N)/N)'; 
    c = [2; ones(N-1,1); 2] .* (-1).^(0:N)'; 
    X = repmat(x,1,N+1); 
    dX = X - X'; 
    D = (c * (1./c)') ./ (dX + eye(N+1)); 
    D = D - diag(sum(D'));
end


%% ************** Solving ******************

b     = 0.0784; 
mu    = 0.0545; 
beta  = 0.09091; 
gamma = 0.125; 
alpha = 0.14286;
sigma_values = [0.00,0.30]; 

S0 = 800; 
E0 = 200; 
I0 = 10;
R0 = 0;

N = 200;              
iteration = 20;      

[D1, x] = chebfunction(N);
t0 = 0; tf = 10;
t = ((tf-t0)/2)*x + ((tf+t0)/2);
D = 2/(tf-t0) * D1;
Ide = eye(N+1);



for value = 1:length(sigma_values)
    sigma = sigma_values(value); 

    S = S0 * ones(N+1,1);
    E = E0 * ones(N+1,1);
    I = I0 * ones(N+1,1);
    R = R0 * ones(N+1,1);

    for iter = 1:iteration
     
        A1 = D + mu * Ide + diag(beta * I); 
        R1 = b * ones(N+1,1);
        A1(N+1,:) = 0; A1(N+1,N+1) = 1; R1(N+1) = S0; 
        S = A1 \ R1;

        A2 = D + (mu + alpha + sigma) * Ide; 
        R2 = beta * S .* I;
        A2(N+1,:) = 0; A2(N+1,N+1) = 1; R2(N+1) = E0;
        E = A2 \ R2;

        A3 = D + (mu + gamma) * Ide;
        R3 = alpha * E;
        A3(N+1,:) = 0; A3(N+1,N+1) = 1; R3(N+1) = I0;
        I = A3 \ R3;

        A4 = D + mu * Ide;
        R4 = gamma * I + sigma * E;
        A4(N+1,:) = 0; A4(N+1,N+1) = 1; R4(N+1) = R0;
        R = A4 \ R4;
    end

    
    figure(value+1);
    plot(t, S, 'b', 'LineWidth', 2); hold on;
    plot(t, E, 'r', 'LineWidth', 2);
    plot(t, I, 'g', 'LineWidth', 2);
    plot(t, R, 'm', 'LineWidth', 2);
    xlim([0 10]); ylim([0 1000]);
    xlabel('Time (days)');
    ylabel('Population');
    title(['SEIR Model Dynamics, \sigma = ', num2str(sigma)]);
    legend('Susceptible','Exposed','Infective','Recovered','Location','best');
    grid on;

end

