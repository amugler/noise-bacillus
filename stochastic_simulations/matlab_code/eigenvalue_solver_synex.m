%{ 
    This code determines the eigenvalues at the fixed points of the SynEx strain (simple model).
%}

global a_k a_m b_k b_m rho mu h p   

a_k = 5;
a_m = 0.3;
b_k = 15;
b_m = 10;
rho = 0.5;
mu = 1;
h = 2;
p = 2;

%SynEx
%y1 = a_k + b_k*x(1)^h/(1+x(1)^h) - x(1)*(1+mu*x(2));
%y2 = a_m + b_m*x(1)^p/(rho^p+x(1)^p) - x(2);

grid=0:.001:4;
S1=zeros(1,numel(grid));
S2=zeros(1,numel(grid));
for j=1:numel(grid)
    
   
  S1(j)=1/mu * ((a_k+b_k*grid(j)^h/(1+grid(j)^h))/grid(j)-1);
  S2(j)= a_m + b_m*grid(j)^p/(rho^p+grid(j)^p);
end
    
loglog(grid,S1,grid,S2);
axis([10e-4 10e0 10e-2 10e1]);
%plot((grid),S1,(grid),S2,a,b);

[X0,Y0] = intersections(grid,S1,grid,S2,0);
numintersects = numel(X0);

syms x1 x2
F_SynEx_simple = a_k + b_k*x1^h/(1+x1^h) - x1*(1+mu*x2);
G_SynEx_simple = a_m + b_m*x1^p/(rho^p+x1^p) - x2;

jakF_SynEx_simple = jacobian(F_SynEx_simple);
jakG_SynEx_simple = jacobian(G_SynEx_simple);
jak_SynEx_simple = [jakF_SynEx_simple(1,1) jakF_SynEx_simple(1,2)
                    jakG_SynEx_simple(1,1) jakG_SynEx_simple(1,2)];
   
                
clear xxx
for i = 1:numel(X0)                 
    x1 = X0(i);
    x2 = Y0(i);

evaljak_SynEx_simple = eval(jak_SynEx_simple);
xxx(1:2,i) = eig(evaljak_SynEx_simple);

end


xxx


