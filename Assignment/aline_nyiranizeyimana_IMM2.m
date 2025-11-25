%% ************** QUESTION 1.C*************************


%***************Cheb function************************

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


% ************** Solving and Ploting ****************

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

    
    figure(value);
    set(gcf, 'Position', [100 100 800 600]);
    plot(t, S, 'b',  'LineWidth', 1); hold on;  
    plot(t, E, 'color','red',  'LineWidth', 1);           
    plot(t, I, 'color','yellow',  'LineWidth', 1);            
    plot(t, R, 'color','magenta',  'LineWidth', 1);          
    
    xlim([0 10]);
    xticks(0:2:10);
     if sigma == 0.3
        ylim([0 800]);
        yticks(0:200:800);
    else
        ylim([0 1000]);
        yticks(0:200:1000);
     end

    title(['Measles Dynamics when \sigma = ', num2str(sigma)]);
    xlabel('Time(years)');
    ylabel('Population');
    
    legend('S(t)','E(t)','I(t)','R(t)');
    grid on;
    set(gca, 'GridColor','k', 'GridAlpha',1, 'LineWidth', 1.5);
    saveas(gcf, ['Figure ', num2str(value), '.png'])

end

%% ********QUESTION 2  ode45 to reproduce  Figure 3 and Figure 4 **********


b = 0.0784;
mu = 0.0545;
beta = 0.09091;
gamma = 0.125;
alpha = 0.14286;


%******* SEIR anonymous function to be able to pass a third parameter
seir_ode_anonymous = @(t, y, sigma) [
    b - beta * y(1) * y(3) - mu * y(1);                     
    beta * y(1) * y(3) - (mu + alpha + sigma) * y(2);      
    alpha * y(2) - (mu + gamma) * y(3);                      
    gamma * y(3) + sigma * y(2) - mu * y(4)                  
];


tspan = [0 10];
y0 = [800 200 10 0];   
sigma_values = [0.3, 0.6];   

for value = 1:length(sigma_values)
    sigma = sigma_values(value);
    [t, Y] = ode45(@(t, y) seir_ode_anonymous(t, y, sigma), tspan, y0);

    figure(value + 2);  
    set(gcf, 'Position', [100 100 800 600]);
    plot(t, Y(:,1), 'b',  'LineWidth', 1); hold on;  
    plot(t, Y(:,2), 'color','red',  'LineWidth', 1);           
    plot(t, Y(:,3), 'color','yellow',  'LineWidth', 1);            
    plot(t, Y(:,4), 'color','magenta',  'LineWidth', 1);          
    
    xlim([0 10]);
    xticks(0:2:10);
    ylim([0 800]);
    yticks(0:200:800);
    
    title(['Measles Dynamics when \sigma = ', num2str(sigma)]);
    xlabel('Time(years)');
    ylabel('Population');
    
    legend('S(t)','E(t)','I(t)','R(t)');
    grid on;
    set(gca, 'GridColor','k', 'GridAlpha',1, 'LineWidth', 1.5);
    saveas(gcf, ['Figure ', num2str(value + 2), '.png'])
end






%% ******************* QUESTION 3. e) Eigen values**********************
L = [
    0.0  0.8  1.2  0.6  0.2;
    0.6  0.0  0.0  0.0  0.0;
    0.0  0.5  0.0  0.0  0.0;
    0.0  0.0  0.4  0.0  0.0;
    0.0  0.0  0.0  0.2  0.0
];
[V,D] = eig(L)
absolute_D = abs(D)
maximum_real_eigen_values = max(diag(absolute_D))

