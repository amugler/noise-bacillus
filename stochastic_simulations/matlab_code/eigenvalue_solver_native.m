%{ 
    This code determines the eigenvalues at the fixed points of the native strain (simple model).
%}

global a_k a_s b_k b_s k_0 k_1 Delta_k Delta_s h p   

a_k = .0002;
a_s = 0;
b_k = .3;
b_s = 3;
k_0 = .2;
k_1 = 1/30;
Delta_k = .1;
Delta_s = .1;
h = 2;
p = 5;

%Native
%y1 = a_k + (b_k)*K^h/(k_0^h+K^h) - K*(Delta_k + 1/(1+K+S));
%y2 = a_s + (b_s)/(1+(K/k_1)^p) - S*(Delta_s + 1/(1+K+S));

grid=0:.0002:1;
S1=zeros(1,numel(grid));
S2=zeros(1,numel(grid));
for j=1:numel(grid)
    
   
  S1(j)=((a_k + (b_k)*grid(j)^h/(k_0^h+grid(j)^h))/grid(j)-Delta_k)^(-1) - grid(j)-1;
  S2(j)= 1/(2*Delta_s)*( ((1+(1+grid(j))*Delta_s-(a_s + (b_s)/(1+(grid(j)/k_1)^p)))^2+4*Delta_s*((1+grid(j))*(a_s + (b_s)/(1+(grid(j)/k_1)^p))))^(1/2)- (1+(1+grid(j))*Delta_s- (a_s + (b_s)/(1+(grid(j)/k_1)^p))));
end
    
loglog(grid,S1,grid,S2);
% axis([.2*10e-3 2*10e-1 10e-5 3*10e1]);

[X0,Y0] = intersections(grid,S1,grid,S2,0);
numintersects = numel(X0);

syms x1 x2
F_native_simple = a_k + (b_k)*x1^h/(k_0^h+x1^h) - x1*Delta_k - x1/(1+x1+x2);
G_native_simple =  a_s + (b_s)/(1+(x1/k_1)^p) - x2*(Delta_s + 1/(1+x1+x2));

jakF_native_simple = jacobian(F_native_simple);
jakG_native_simple = jacobian(G_native_simple);
jak_native_simple = [jakF_native_simple(1,1) jakF_native_simple(1,2)
                     jakG_native_simple(1,1) jakG_native_simple(1,2)];

clear xxx
for i = 1:numel(X0)                 
    x1 = X0(i);
    x2 = Y0(i);

evaljak_native_simple = eval(jak_native_simple);
xxx(1:2,i) = eig(evaljak_native_simple);

end


xxx


