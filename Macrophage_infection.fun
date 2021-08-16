function  var = Macrophage_infection2(~,y,k1,k2,k_m1,k_m2,delta_MA,s1,s2,M_max,g_M,delta_MR,...
                       q_MRV,q_MRI,rho,alpha_M,f_D,...
                       p_CM1,p_CM2,D_50,delta_C,...
                       g_T,T_max,beta,beta_prime,q_prime,phi,xi_R,delta_I,kappa_N,kappa_E,...
                       p_I,p_M,delta_V,kappa_MV,kappa_AS,kappa_AL,...
                       q_FI,q_FM,delta_F,kappa_D,delta_D,...
                       gamma_E,V_50E,n_E,tau_E,phi_E,delta_E,...
                       gamma_B,V_50B,n_B,tau_B,phi_p,delta_p,mu_AS,delta_AS,mu_AL,delta_AL)
            
n_total = 13 + n_E + n_B + 3;
var = zeros(n_total,1);


var(1) = alpha_M*(y(6)+f_D*y(10)) + g_M * (1 - (y(1) + y(2) + y(3))/M_max) * y(1) - k1 /(1 + s1*(y(3)/M_max))*y(1) - k2 * (1 + s2* (y(2)/M_max))*y(1) + k_m1 * y(2) + k_m2 * y(3) - q_MRV*y(1)*y(7) - q_MRI*y(1)*y(6) - rho*y(4)*y(1) - delta_MR*y(1);
var(2) = k1/(1 + s1*(y(3)/M_max))*y(1) - k_m1 * y(2) - delta_MA * y(2) + q_MRV*y(1)*y(7) + q_MRI*y(1)*y(6);
var(3) = k2 * (1 + s2* (y(2)/M_max))*y(1) - k_m2 * y(3) - delta_MA * y(3) + rho*y(4)*y(1);
var(4) = p_CM1*y(2)*y(10)/(y(10)+D_50) + p_CM2*y(3) - delta_C*y(4);

var(5) = g_T*(y(5)+y(9))*(1-(y(5)+y(6)+y(9))/T_max) - beta*y(5)*y(7) - phi*y(8)*y(5) + xi_R*y(9); % T
var(6) = beta*y(5)*y(7) - delta_I*y(6)  - kappa_E*y(6)*y(13+n_E-1) - kappa_N*y(6)*y(8); % I
var(7) = p_I*y(6) + p_M*y(3) - delta_V*y(7) - beta_prime*y(7)*y(5)- q_prime*y(1)*y(7) - kappa_MV*y(2)*y(7) - kappa_AS*y(7)*y(13+n_E+n_B+2) - kappa_AL*y(7)*y(13+n_E+n_B+3); % V

var(8) = q_FI*y(6) + q_FM*y(2) - delta_F*y(8); % F
var(9) = phi*y(8)*y(5) - xi_R*y(9); % R
var(10) = delta_I*y(6)  + kappa_E*y(13+n_E-1)*y(6) + kappa_N*y(6)*y(8)- kappa_D*y(2)*y(10) - delta_D*y(10); % D


var(11) = -gamma_E * y(7) / (y(7) + V_50E) * y(11); % E0

var(12) = gamma_E * y(7) / (y(7) + V_50E) * y(11) -  n_E/tau_E * y(12); % E1 

for i = 13 : 13+n_E-2
    
    var(i) = n_E/tau_E * (y(i-1) - y(i));
    
end 
var(13+n_E-1) = phi_E * n_E/tau_E * y(13+n_E-2) - delta_E * y(13+n_E-1); % E



var(13+n_E) = -gamma_B * y(7) / (y(7) + V_50B) * y(13+n_E); % B0

var(13+n_E+1) = gamma_B * y(7) / (y(7) + V_50B) * y(13+n_E) -  n_B/tau_B * y(13+n_E+1); % B1

for j = 13+n_E+2 : 13+n_E+2+n_B-2 % B2...B_{n_B}
    
    var(j) = n_B/tau_B * (y(j-1) - y(j));
    
end 

var(13+n_E+n_B+1) = phi_p * n_B/tau_B *  y(13+n_E+2+n_B-2) - delta_p * y(13+n_E+n_B+1); % P
 
var(13+n_E+n_B+2) = mu_AS * y(13+n_E+n_B+1) - delta_AS * y(13+n_E+n_B+2); % A_S
var(13+n_E+n_B+3) = mu_AL * y(13+n_E+n_B+1) - delta_AL * y(13+n_E+n_B+3); % A_L
end 




