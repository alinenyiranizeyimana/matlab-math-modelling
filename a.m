function [D, x] = cheb(N)

    if N==0
        D = 0; 
        x = 1;
        return;
    end

    x = cos(pi*(0:N)/N)';      % Chebyshev-Gauss-Lobatto nodes
    c = [2; ones(N-1,1); 2].*(-1).^(0:N)'; 

    X = repmat(x,1,N+1);
    dX = X - X';

    D = (c*(1./c)')./(dX + eye(N+1));
    D = D - diag(sum(D'));
end




%% ---------------- Parameters ----------------
b     = 0.0784;
mu    = 0.0545;
beta  = 0.09091;
gamma = 0.125;
alpha = 0.14286;

%% ---------------- Initial Conditions ----------------
S0 = 800; 
E0 = 180; 
I0 = 20; 
R0 = 0;

%% ---------------- Chebyshev Points & Differentiation ----------------
N = 240;                      % Number of spectral points                         % identity matrix for linear terms

%% ---------------- Sigma Cases ----------------
sigma_values = [0.00, 0.30];

S = S0 * ones(N+1,1);
E = E0 * ones(N+1,1);
I = I0 * ones(N+1,1);
R = R0 * ones(N+1,1);

for s = 1:length(sigma_values)
    
    sigma = sigma_values(s);

    % Initialize solution vectors
    S = S0 * ones(N+1,1);
    E = E0 * ones(N+1,1);
    I = I0 * ones(N+1,1);
    R = R0 * ones(N+1,1);

    iteration = 100;   % relaxation iterations

    for k = 1:iteration

        % --- Susceptible ---
       
        R1 = b - beta * S .* I;
        A1 = D  - mu *Ide
        A1(N+1,:) = 0; A1(N+1,N+1) = 1; R1(N+1) = S0;
        S = A1 \ R1;

        % --- Exposed ---
        A2 = D + (mu + alpha + sigma) * Ide;
        R2 = beta * S .* I;
        A2(N+1,:) = 0; A2(N+1,N+1) = 1; R2(N+1) = E0;
        E = A2 \ R2;

        % --- Infective ---
        A3 = D + (mu + gamma) * Ide;
        R3 = alpha * E;
        A3(N+1,:) = 0; A3(N+1,N+1) = 1; R3(N+1) = I0;
        I = A3 \ R3;

        % --- Recovered ---
        A4 = D + mu * Ide;
        R4 = gamma * I + sigma * E;
        A4(N+1,:) = 0; A4(N+1,N+1) = 1; R4(N+1) = R0;
        R = A4 \ R4;

    end

   
    % --- Plot results ---
    figure(s)
    hold on
    plot(t, S, 'b', 'LineWidth', 2)
    plot(t, E, 'r', 'LineWidth', 2)
    plot(t, I, 'g', 'LineWidth', 2)
    plot(t, R, 'm', 'LineWidth', 2)

    xlim([0,10])
    ylim([0,800])
    xticks(0:2:10)
    yticks(0:200:800)
    title(['SEIR Model Dynamics for \sigma = ', num2str(sigma)]);
    xlabel('Time (days)');
    ylabel('Population');
    legend('Susceptible','Exposed','Infective','Recovered','Location','best');
    grid on
    hold off

end


