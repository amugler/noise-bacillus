%{ 
    This code creates stochastic and deterministic time series, nullclines, and power spectra plots for the native strain
    (simple model, low molecule number).
%}

global xvals yvals alpha_k alpha_s beta_k beta_s k_k k_s p h delta lambda_k lambda_s Gamma_k Gamma_s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_k = .00013;
ak = 00013;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%loop data
step_size = .00005;
loop_length = 1;


%general data for deterministic plotter
a_s = 0;
b_k = 0.3;
b_s = 3;
k_0 = 0.2;
k_1 = 1;
k_1 = k_1/30;
Delta_k = 0.1;
Delta_s = 0.1;

Gamma_k = 100;
Gamma_s = 1;
delta = .001;
h = 2;
p = 5;
alpha_k = a_k*Gamma_k*delta;
alpha_s = a_s*Gamma_s*delta;
beta_k = b_k*Gamma_k*delta;
beta_s = b_s*Gamma_s*delta;
k_k = k_0*Gamma_k;
k_s = k_1*Gamma_k;
lambda_k = Delta_k*delta;
lambda_s = Delta_s*delta;


% vec_ak = zeros(1,loop_length);
% vec_A = zeros(1,loop_length);
% vec_mu = zeros(1,loop_length);
% vec_sigma = zeros(1,loop_length);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:loop_length;
    
%Reads the proper data file and creates variables
ak_str=num2str(ak, '%.05i');
ak_filename = sprintf('gns_low_TKS_%s.txt',ak_str);
gns_var_TKS = dlmread(ak_filename,',',0,0);
ak_varname_time_00013 = gns_var_TKS(:,1);
ak_varname_Kvals_00013 = gns_var_TKS(:,2);
ak_varname_Svals_00013 = gns_var_TKS(:,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Fourier transform and power spectrum calculations
adjust_ak_varname_Kvals_00013 = zeros(1,numel(ak_varname_Kvals_00013));
mean_ak_varname_Kvals_00013 = mean(ak_varname_Kvals_00013);
for j = 1:numel(ak_varname_time_00013)
   adjust_ak_varname_Kvals_00013(j) = ak_varname_Kvals_00013(j)-mean_ak_varname_Kvals_00013; 
end
adjust_ak_varname_time_00013 = ak_varname_time_00013;

tops = max(adjust_ak_varname_time_00013);
fs = .01;   % sample frequency (Hz)
fourier_grid=0:1/fs:tops-1/fs;  % length of sample
vq=interp1(adjust_ak_varname_time_00013,adjust_ak_varname_Kvals_00013,fourier_grid); % interpolating data on fixed separation grid
m = numel(vq);  % window length
n = pow2(nextpow2(m));  % transform length
y = fft(vq,n);  % DFT
f = (0:n-1)*(fs/n); % frequency range - x axis
power = y.*conj(y)/n;   % power of the DFT

% Adjusts data for proper length and to ignore low frequencies
f_cutoff = 1e-4;
adjust_power_00013=power(f<f_cutoff);
adjust_power_00013(1:15)=[];
adjust_f_00013=f(f<f_cutoff);
adjust_f_00013(1:15)=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calculates the Gaussian fit. Input the guesses for A, mu, and sigma2
xvals = adjust_f_00013;
yvals = adjust_power_00013;
options = optimset('MaxIter',50000,'MaxFunEvals',50000);
Amusigma2 = fminsearch(@gaussian_fitter,[1,1e-5,1e-10],options);

%stores A,mu, and sigma in vector
% vec_ak(i) = ak/100000;
% vec_A(i) = Amusigma2(1);
% vec_mu(i) = Amusigma2(2);
% vec_sigma(i) = sqrt(Amusigma2(3));

%creates Gaussian curve
yyy_00013=zeros(1, numel(xvals));
for j=1:numel(xvals)
yyy_00013(j) = Amusigma2(1)/sqrt(2*pi*Amusigma2(3))*exp(-(xvals(j)-Amusigma2(2))^2/2/Amusigma2(3));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% = starttime:stepsize:maxtime for the deterministic
time_length_hours = 500;
timesint=0:10:time_length_hours*3600;

z0=[3,1];

options = odeset('RelTol',2.22045e-14,'AbsTol',1e-16,'MaxStep',00013);
[deterministic_time,Z]=ode113(@integrator_native_simple,timesint,z0,options);

%defines variable for the deterministic y coordinate values
gns_deterministic_Kvals_00013 = Z(:,1);
gns_deterministic_Svals_00013 = Z(:,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%nullcline equations
% grid=0:.0001:2;
% gns_deterministic_null1_00013=zeros(1,numel(grid));
% gns_deterministic_null2_00013=zeros(1,numel(grid));
% 
% for j=1:numel(grid)
%   gns_deterministic_null1_00013(j)=((a_k + (b_k)*grid(j)^h/(k_0^h+grid(j)^h))/grid(j)-Delta_k)^(-1) - grid(j)-1;
%   gns_deterministic_null2_00013(j)= 1/(2*Delta_s)*( ((1+(1+grid(j))*Delta_s-(a_s + (b_s)/(1+(grid(j)/k_1)^p)))^2+4*Delta_s*((1+grid(j))*(a_s + (b_s)/(1+(grid(j)/k_1)^p))))^(1/2)- (1+(1+grid(j))*Delta_s- (a_s + (b_s)/(1+(grid(j)/k_1)^p))));
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%adjust viewing window for deterministic and stochastic figures
% ak_varname_Kvals_00013(1:5861)=[];
% ak_varname_Svals_00013(1:5861)=[];
% ak_varname_time_00013(numel(ak_varname_time_00013)-5861+1:numel(ak_varname_time_00013))=[];
% deterministic_time(1:5861)=[];
% gns_deterministic_Kvals_00013(numel(gns_deterministic_Kvals_00013)-5861+1:numel(gns_deterministic_Kvals_00013))=[];
% gns_deterministic_Svals_00013(numel(gns_deterministic_Svals_00013)-5861+1:numel(gns_deterministic_Svals_00013))=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Plots nullclines
% figure(1);
% loglog((grid),gns_deterministic_null1_00013,(grid),gns_deterministic_null2_00013,gns_deterministic_Kvals_00013/100,gns_deterministic_Svals_00013,'k');
% axis([.2*10e-3 2*10e-1 10e-5 2*10e1]);
% title('a_k = 0.00013');
% ax=gca;
% xlabel('Kbar','FontSize',8);
% ylabel('Sbar','FontSize',8);
% ax.XTick = [10e-3,10e-2,10e-1,10e0];
% ax.YTick = [10e-3,10e-1, 10e1];


%plots deterministic time series for ComK
figure(2);
plot(deterministic_time/3600,gns_deterministic_Kvals_00013,'k')
axis([0 100 0 115])
title('Deterministic','FontSize',15);
ax=gca;
xlabel('Time, \tau','FontSize',9);
ylabel('Molecules of ComK','FontSize',9);
ax.XTick = [0, 50, 100];
ax.YTick = [0, 50, 100];


%plots deterministic time series for ComS
figure(3);
plot(deterministic_time/3600,gns_deterministic_Svals_00013,'k')
axis([0 100 0 50])
title('Deterministic','FontSize',15);
ax=gca;
ylabel('Molecules of ComS','FontSize',9);
ax.XTick = [0, 50, 100];
ax.YTick = [0, 50, 100];


%Plots stochastic time series values for ComK
figure(4);
plot(ak_varname_time_00013/3600,ak_varname_Kvals_00013,'b')
axis([0 200 0 160])
title('Stochastic','FontSize',15);
ax=gca;
ylabel('Molecules of ComK','FontSize',9);
% ax.XTick = [0, 25, 50];
% ax.YTick = [0, 50, 100];


%Plots stochastic time series values for ComS
figure(5);
plot(ak_varname_time_00013/3600,ak_varname_Svals_00013,'r')
axis([0 200 0 110])
ax=gca;
xlabel('Time, \tau','FontSize',8);
ylabel('Molecules of ComS','FontSize',9);
% ax.XTick = [0, 25, 50];
% ax.YTick = [0, 50, 100];


%Plots Gaussian fit to periodogram
figure(6);
plot(adjust_f_00013,adjust_power_00013,'b',adjust_f_00013,yyy_00013,'r');
title('Periodogram');
ax=gca;
xlabel('frequency (Hz), \tau','FontSize',9);
ylabel('Power','FontSize',9);
legend('Power spectrum','Gaussian curve fit','Location','northeast');



%loop propagation
i
a_k = a_k + step_size;
end

