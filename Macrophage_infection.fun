function  var = Macrophage_infection(~,y,k1,k2,k_m1,k_m2,delta_M,g_M,...
                        s, q1, q_m1, q_m2, q2, gamma, g_T,Tmax, beta, delta_I,...
                        p, p_M, delta_V,...
                        alpha1, gamma_E, I50, n_E, tau_E, phi_E, delta_E,...
                        alpha2, gamma_B, V50, n_B, tau_B, phi_P, delta_P,...
                        mu_AS, delta_AS, mu_AL, delta_AL,...
                        kappa_N, kappa_E, kappa_AS, kappa_AL)

var = zeros(8+n_E+2+n_B+2,1);


var(1) = g_M - delta_M * y(1) - k1 * y(1) - q1 * y(6) * y(1) +  k_m1 * y(2) - ...
         + q_m1 * (y(3) + gamma * y(8+n_E)) * y(2) -...
         k2 * y(1) - q2 * (y(3) + gamma * y(8+n_E))* y(1) + k_m2 * y(3) + q_m2 * y(6) * y(3); % MR


var(2) = k1 * y(1) + q1 * y(6) * y(1) + s * y(2) - k_m1 * y(2) - q_m1 * (y(3) + gamma * y(8+n_E)) * y(2) - delta_M * y(2); % M1
 
var(3) = k2 * y(1) + q2 * (y(3) + gamma * y(8+n_E))* y(1) - k_m2 * y(3) - delta_M * y(3) - q_m2 * y(6) * y(3); % M2


% =================================== % 
var(4) = g_T * (y(4) + y(5)) * (1 - (y(4) + y(5))/ Tmax) - beta * y(4) * y(6); %T

var(5) = beta * y(4) * y(6) - delta_I * y(5) - kappa_N * y(2) * y(5) - kappa_E * y(8+n_E) * y(5); %I

var(6) = p * y(5) + p_M * y(2) - delta_V * y(6) - kappa_AS * y(8+n_E+2+n_B+1) * y(6) - kappa_AL * y(8+n_E+2+n_B+2) * y(6); %V


% =================================== % 
var(7) = alpha1 * y(2) - gamma_E * y(5) / (y(5) + I50) * y(7); % E0

var(8) = gamma_E * y(5) / (y(5) + I50) * y(7) - n_E/tau_E * y(8); % E1 

for i = 8+1 : 8+(n_E -1) % E2...E_{n_E}
    
    var(i) = n_E/tau_E * (y(i-1) - y(i));
    
end 

var(8+n_E) = phi_E * n_E/tau_E * y(8+n_E-1) - delta_E * y(8+n_E); % E

var(8+n_E+1) = alpha2 * y(3) - gamma_B * y(6) / (y(6) + V50) * y(8+n_E+1); % B0

var(8+n_E+2) = gamma_B * y(6) / (y(6) + V50) * y(8+n_E+1) - n_B/tau_B * y(8+n_E+2); % B1

for j = 8+n_E+2+1 : 8+n_E+2+(n_B - 1) % B2...B_{n_B}
    
    var(j) = n_B/tau_B * (y(j-1) - y(j));
    
end 

var(8+n_E+2+n_B) = phi_P * n_B/tau_B *  y(8+n_E+2+n_B-1) - delta_P * y(8+n_E+2+n_B); % P
 
var(8+n_E+2+n_B+1) = mu_AS * y(8+n_E+2+n_B) - delta_AS * y(8+n_E+2+n_B+1); % A_S

var(8+n_E+2+n_B+2) = mu_AL * y(8+n_E+2+n_B) - delta_AL * y(8+n_E+2+n_B+2); % A_L


end 
