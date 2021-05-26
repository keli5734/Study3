function  var = Macrophage_homeostasis(~,y,k1,k2,k_m1,k_m2,delta_M,g_M)

var = zeros(3,1);

var(1) = g_M - delta_M * y(1) - k1 * y(1) - k2 * y(1) + k_m1 * y(2) + k_m2 * y(3);

var(2) = k1 * y(1) - k_m1 * y(2) - delta_M * y(2);

var(3) =  k2 * y(1) - k_m2 * y(3) - delta_M * y(3);


end 
