%{ 
  This code numerically integrates the differential equations for the SynExSlow strain (simple model).
%}

function dz=integrator_SynExSlow(~,z)

global alpha_k alpha_m alpha_s beta_k beta_m beta_s k_k k_m k_s gamma_k gamma_s p h delta lambda
%dz = zeros(#vectors,1);
dz=zeros(3,1);

%a=Z(:,1);

dz(1)= alpha_k + beta_k*(z(1)^h)/((k_k^h)+(z(1)^h)) - (delta*z(1)*z(2)/(1+z(1)/gamma_k+z(3)/gamma_s)) - lambda*z(1);
dz(2)= alpha_m + beta_m*(z(1)^p)/((k_m^p)+(z(1)^p)) - lambda*z(2);
dz(3)= alpha_s + beta_s*(z(1)^h)/((k_s^h)+(z(1)^h)) - (delta*z(3)*z(2)/(1+z(1)/gamma_k+z(3)/gamma_s)) - lambda*z(3);

end


