%{ This code determines the eigenvalues at the fixed points of the SynExSlow strain (simple model).
%}

a_k = .3;

alpha_k = a_k*.04;
alpha_m = 0.075;
alpha_s = 0.5;
beta_k = 7.5;
beta_m = 2.5;
beta_s = 0.5;
k_k = 5000;
k_m = 2500;
k_s = 500;
gamma_k = 25000;
gamma_s = 20;
delta = 2e-6;
lambda = 1e-4;
h = 2;
p = 2;

% eqn1 = (1+Kvals/gamma_k + Svals/gamma_s)/(delta*Kvals)*(alpha_k + (beta_k*Kvals^h)/(k_k^h+Kvals^h) - lambda*Kvals)- Mvals;
% eqn2 = (alpha_m + (beta_m*Kvals^p)/(k_m^p+Kvals^p))/lambda -Mvals;
% eqn3 = (1+Kvals/gamma_k + Svals/gamma_s)/(delta*Svals)*(alpha_s + (beta_s*Kvals^h)/(k_s^h+Kvals^h) - lambda*Svals)- Mvals;
               
syms Kvals Svals Mvals
[sol_Kvals, sol_Mvals, sol_Svals] = vpasolve([(1+Kvals/gamma_k + Svals/gamma_s)/(delta*Kvals)*(alpha_k + (beta_k*Kvals^h)/(k_k^h+Kvals^h) - lambda*Kvals)- Mvals ==0, ...
                                 Mvals == (alpha_m + (beta_m*Kvals^p)/(k_m^p+Kvals^p))/lambda, ...
                                 Mvals == (1+Kvals/gamma_k + Svals/gamma_s)/(delta*Svals)*(alpha_s + (beta_s*Kvals^h)/(k_s^h+Kvals^h) - lambda*Svals)], [Kvals Mvals Svals]);

% syms Kvals Svals Mvals
% [sol_Kvals, sol_Mvals, sol_Svals] = vpasolve([alpha_k + beta_k*Kvals^h/(k_k^h+Kvals^h)-delta*Kvals*Mvals/(1+Kvals/gamma_k+Svals/gamma_s)-lambda*Kvals ==0, ...
%                                               alpha_m + beta_m*Kvals^p/(k_m^p+Kvals^p)-lambda*Mvals ==0, ...
%                                               alpha_s + beta_s*Kvals^h/(k_s^h+Kvals^h)-delta*Svals*Mvals/(1+Kvals/gamma_k+Svals/gamma_s)-lambda*Svals ==0], [Kvals Svals Mvals]);
                                                      
                             
loopcountreal = 1;
clear vec_intersects
for j = 1:numel(sol_Kvals)
    clear isrealchecker
    isrealchecker = isreal(sol_Kvals(j));
    if isrealchecker == 1 && sol_Kvals(j)>0 && sol_Svals(j)>0 && sol_Mvals(j)>0
        vec_intersects(loopcountreal,1)=sol_Kvals(j);
        vec_intersects(loopcountreal,2)=sol_Mvals(j);
        vec_intersects(loopcountreal,3)=sol_Svals(j);
        loopcountreal = loopcountreal +1;
    end
end
vec_intersects
numintersects = numel(vec_intersects(:,1));


syms x1 x2 x3
F_SynExSlow_high = alpha_k + beta_k*x1^h/(k_k^h + x1^h) - delta*x1*x2/(1+x1/gamma_k + x3/gamma_s)-lambda*x1;
G_SynExSlow_high = alpha_m + beta_m*x1^p/(k_m^p + x1^p) - lambda*x2 + 0*x3;
H_SynExSlow_high = alpha_s + beta_s*x1^h/(k_s^h + x1^h) - delta*x2*x3/(1+x1/gamma_k + x3/gamma_s)-lambda*x3;

jakF_SynExSlow_high = jacobian(F_SynExSlow_high);
jakG_SynExSlow_high = jacobian(G_SynExSlow_high);
jakH_SynExSlow_high = jacobian(H_SynExSlow_high);
jak_native_simple = [jakF_SynExSlow_high(1,1) jakF_SynExSlow_high(1,2) jakF_SynExSlow_high(1,3)
                     jakG_SynExSlow_high(1,1) jakG_SynExSlow_high(1,2) 0
                     jakH_SynExSlow_high(1,1) jakH_SynExSlow_high(1,2) jakH_SynExSlow_high(1,3) ];

clear xxx
for i = 1:numintersects              
    x1 = vec_intersects(i,1);
    x2 = vec_intersects(i,2);
    x3 = vec_intersects(i,3);

evaljak_native_simple = eval(jak_native_simple);
xxx(1:3,i) = eig(evaljak_native_simple);

end


xxx


