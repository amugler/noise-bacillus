%{ 
  This code numerically integrates the differential equations for the SynEx strain (simple model).
%}

function dz=integrator_SynEx_simple(~,z)

global alpha_k alpha_m beta_k beta_m k_k k_m p h delta lambda
%dz = zeros(#vectors,1);
dz=zeros(2,1);

%a=Z(:,1);


dz(1)= alpha_k + beta_k*(z(1)^h)/((k_k^h)+(z(1)^h)) - (delta*z(2) + lambda)*z(1);
dz(2)= alpha_m + beta_m*(z(1)^p)/((k_m^p)+(z(1)^p)) - lambda*z(2);

end


