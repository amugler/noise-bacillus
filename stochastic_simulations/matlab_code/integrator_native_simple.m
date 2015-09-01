%{ This code numerically integrates the differential equations for the native strain (simple model).
%}

function dz=integrator_native_simple(~,z)

global alpha_k alpha_s beta_k beta_s k_k k_s p h delta lambda_k lambda_s Gamma_k Gamma_s
%dz = zeros(#vectors,1);
dz=zeros(2,1);

%a=Z(:,1);



dz(1)=alpha_k + (beta_k*z(1)^h)/(k_k^h+z(1)^h) - (delta/(1+z(1)/Gamma_k+z(2)/Gamma_s)+lambda_k)*z(1);
dz(2)=alpha_s + (beta_s)/(1+(z(1)/k_s)^p) - (delta/(1+z(1)/Gamma_k+z(2)/Gamma_s)+lambda_s)*z(2);


end
