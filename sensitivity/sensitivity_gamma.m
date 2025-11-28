r = 0.5;
K= 1000;
beta = 0.1;
delta = 0.2;
gamma_values = [0.014 0.016 0.018 0.02 0.022 0.024 0.026];
mu = 0.1;


 R_0=zeros(1,length(gamma_values));
 sens_index =zeros(1,length(gamma_values));
 

for i=1:length(gamma_values)
    gamma= gamma_values(i);
    R_0(i) = beta*K* (r -delta) / (r*(gamma+mu+delta));
    sens_index(i) = -gamma/(gamma+mu+delta);
    
end

R_0
sens_index
