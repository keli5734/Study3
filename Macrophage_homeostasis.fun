function  var = Macrophages_Homeostasis2(~,y,k1,k2,k_m1,k_m2,delta_MA,s1,s2,M_max,g_M,delta_MR)

var = zeros(3,1);

var(1) = g_M * (1 - (y(1) + y(2) + y(3))/M_max) * y(1) - k1 /(1 + s1*(y(3)/M_max))*y(1) - k2 * (1 + s2* (y(2)/M_max))*y(1) + k_m1 * y(2) +  k_m2 * y(3) - delta_MR*y(1);
var(2) = k1/(1 + s1*(y(3)/M_max))*y(1) - k_m1 * y(2) - delta_MA * y(2);
var(3) = k2 * (1 + s2* (y(2)/M_max))*y(1) - k_m2 * y(3) - delta_MA * y(3);


end 
