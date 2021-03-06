%{ 
    This code creates stochastic and deterministic time series, nullclines, and power spectra plots for the SynExSlow strain
    (simple model, high molecule number).
%}

global xvals yvals alpha_k alpha_m alpha_s beta_k beta_m beta_s k_k k_m k_s gamma_k gamma_s p h delta lambda
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_k = 150.0;
ak = 1500000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%general data for deterministic plotter

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


% vec_ak = zeros(1,loop_length);
% vec_A = zeros(1,loop_length);
% vec_mu = zeros(1,loop_length);
% vec_sigma = zeros(1,loop_length);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Reads the proper data file and creates variables
ak_str=num2str(ak, '%.05i');
ak_filename = sprintf('gssslow_high_TKMS_%s.txt',ak_str);
gssslow_high_var_TKMS = dlmread(ak_filename,',',0,0);
ak_varname_time_1500000 = gssslow_high_var_TKMS(:,1);
ak_varname_Kvals_1500000 = gssslow_high_var_TKMS(:,2);
ak_varname_Mvals_1500000 = gssslow_high_var_TKMS(:,3);
ak_varname_Svals_1500000 = gssslow_high_var_TKMS(:,4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Fourier transform and power spectrum calculations
adjust_ak_varname_Kvals_1500000 = zeros(1,numel(ak_varname_Kvals_1500000));
mean_ak_varname_Kvals_1500000 = mean(ak_varname_Kvals_1500000);
for j = 1:numel(ak_varname_time_1500000)
   adjust_ak_varname_Kvals_1500000(j) = ak_varname_Kvals_1500000(j)-mean_ak_varname_Kvals_1500000; 
end
adjust_ak_varname_time_1500000 = ak_varname_time_1500000;

tops = max(adjust_ak_varname_time_1500000);
fs = .01;   % sample frequency (Hz)
fourier_grid=0:1/fs:tops-1/fs;  % length of sample
vq=interp1(adjust_ak_varname_time_1500000,adjust_ak_varname_Kvals_1500000,fourier_grid); % interpolating data on fixed separation grid
m = numel(vq);  % window length
n = pow2(nextpow2(m));  % transform length
y = fft(vq,n);  % DFT
f = (0:n-1)*(fs/n); % frequency range - x axis
power = y.*conj(y)/n;   % power of the DFT

%Adjusts data for proper length and to ignore low frequencies
f_cutoff = 1e-4;
adjust_power_1500000=power(f<f_cutoff);
adjust_power_1500000(1:15)=[];
adjust_f_1500000=f(f<f_cutoff);
adjust_f_1500000(1:15)=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calculates the Gaussian fit. Input the guesses for A, mu, and sigma2
xvals = adjust_f_1500000;
yvals = adjust_power_1500000;
options = optimset('MaxIter',50000,'MaxFunEvals',50000);
Amusigma2 = fminsearch(@gaussian_fitter,[1e6,1.2e-5,6e-12],options);

%creates Gaussian curve
yyy_1500000=zeros(1, numel(xvals));
for j=1:numel(xvals)
yyy_1500000(j) = Amusigma2(1)/sqrt(2*pi*Amusigma2(3))*exp(-(xvals(j)-Amusigma2(2))^2/2/Amusigma2(3));
end

% stores A,mu, and sigma in vector
% vec_ak(i) = ak/10000;
% vec_A(i) = Amusigma2(1);
% vec_mu(i) = Amusigma2(2);
% vec_sigma(i) = sqrt(Amusigma2(3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % = starttime:stepsize:maxtime for the deterministic
time_length_hours = 500;
timesint=0*3600:100:time_length_hours*3600;

z0=[1000, 1000, 1000];

options = odeset('RelTol',2.22045e-14,'AbsTol',1e-16,'MaxStep',1000);
[deterministic_time,Z]=ode113(@integrator_SynExSlow,timesint,z0,options);

%defines variable for the deterministic y coordinate values
gssslow_deterministic_Kvals_1500000 = Z(:,1);
gssslow_deterministic_Mvals_1500000 = Z(:,2);
gssslow_deterministic_Svals_1500000 = Z(:,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%nullcline equations
% grid=0:.001:10;
% gss_deterministic_null1_1500000=zeros(1,numel(grid));
% gss_deterministic_null2_1500000=zeros(1,numel(grid));
% 
% for j=1:numel(grid)
%   gss_deterministic_null1_1500000(j)=1/mu * ((a_k+b_k*grid(j)^h/(1+grid(j)^h))/grid(j)-1);
%   gss_deterministic_null2_1500000(j)= a_m + b_m*grid(j)^p/(rho^p+grid(j)^p);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%adjust viewing window for deterministic and stochastic figures
% ak_varname_Kvals_1500000(1:120)=[];
% ak_varname_Mvals_1500000(1:120)=[];
% ak_varname_time_1500000(1001:1120)=[];
% deterministic_time(1001:1100)=[];
% gss_deterministic_Kvals_1500000(1:100)=[];
% gss_deterministic_Svals_1500000(1:100)=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Plots nullclines
% figure(1);
% loglog((grid),Gillespie_SynEx_simple_deterministic_null1_01,(grid),Gillespie_SynEx_simple_deterministic_null2_01,Gillespie_SynEx_simple_deterministic_Kvals_01/k_k,Gillespie_SynEx_simple_deterministic_Mvals_01/k_m,'k');
% axis([.2*10e-3 10e0 10e-3 5*10e1]);
% title('a_k = 0.01');
% ax=gca;
% xlabel('Kbar','FontSize',8);
% ylabel('Mbar','FontSize',8);
% ax.XTick = [10e-3,10e-2,10e-1,10e0];
% ax.YTick = [10e-3,10e-1, 10e1];

%plots deterministic time series for ComK
figure(2);
plot(deterministic_time/3600,gssslow_deterministic_Kvals_1500000,'k','LineWidth',1)
% axis([0 100 0 110 ]);
title('Deterministic ComK','FontSize',15);
ax=gca;
xlabel('Time, \tau','FontSize',9);
ylabel('Molecules of ComK','FontSize',9);
% ax.XTick = [0, 50, 100];
% ax.YTick = [0, 50, 100];

%plots deterministic time series for MecA
% figure(3);
% plot(deterministic_time/3600,gss_deterministic_Mvals_1500000,'k','LineWidth',1)
% axis([0 100 0 110 ]);
% title('Deterministic MecA','FontSize',15);
% ax=gca;
% xlabel('time \tau','FontSize',9);
% ylabel('Molecules of MecA','FontSize',9);
% ax.XTick = [0, 50, 100];
% ax.YTick = [0, 50, 100];

%plots deterministic time series for ComS
% figure(4);
% plot(deterministic_time/3600,gss_deterministic_Svals_1500000,'k','LineWidth',1)
% axis([0 100 0 110 ]);
% title('Deterministic ComS','FontSize',15);
% ax=gca;
% xlabel('time \tau','FontSize',9);
% ylabel('Molecules of ComS','FontSize',9);
% ax.XTick = [0, 50, 100];
% ax.YTick = [0, 50, 100];

%Plots stochastic time series values for ComK
figure(5);
plot(ak_varname_time_1500000/3600,ak_varname_Kvals_1500000,'b')
axis([0 200 0 3e4])
ax=gca;
ylabel('Molecules of ComK','FontSize',9);
% ax.XTick = [0, 25, 50];
% ax.YTick = [0, 50, 100];

%Plots stochastic time series values for MecA
% figure(6);
% plot(ak_varname_time_1500000/3600,ak_varname_Mvals_1500000,'r')
% axis([0 200 0 3e4])
% ax=gca;
% xlabel('Time, \tau','FontSize',8);
% ylabel('Molecules of MecA','FontSize',9);
% % ax.XTick = [0, 25, 50];
% % ax.YTick = [0, 50, 100];

%Plots stochastic time series values for ComS
% figure(7);
% plot(ak_varname_time_1500000/3600,ak_varname_Svals_1500000,'k')
% axis([0 200 0 3e4])
% ax=gca;
% xlabel('Time, \tau','FontSize',8);
% ylabel('Molecules of ComS','FontSize',9);
% % ax.XTick = [0, 25, 50];
% % ax.YTick = [0, 50, 100];

%Plots Gaussian fit to periodogram
figure(8);
plot(adjust_f_1500000,adjust_power_1500000,'b',adjust_f_1500000,yyy_1500000,'r');
title('Periodogram');
ax=gca;
xlabel('frequency (Hz), \tau','FontSize',9);
ylabel('Power','FontSize',9);
legend('Power spectrum','Gaussian curve fit','Location','northeast');


