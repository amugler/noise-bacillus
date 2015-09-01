%{ 
    This code creates stochastic and deterministic time series, nullclines, and power spectra plots for the SynEx strain
    (simple model, low molecule number).
%}

global xvals yvals alpha_k alpha_m beta_k beta_m k_k k_m p h delta lambda
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_k = 0.03;
ak = 00300;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%loop data
step_size = .00010;
loop_length = 1;


%general data for deterministic plotter
a_m = 0.3;
b_k = 15;
b_m = 10;
rho = 0.5;
mu = 1;

h = 2;
p = 2;
k_k = 10;
k_m = rho*k_k;
lambda = .0001;
alpha_k = a_k*k_k*lambda;
alpha_m = a_m*k_m*lambda;
beta_k = b_k*k_k*lambda;
beta_m = b_m*k_m*lambda;
delta = mu*lambda/k_m;

% vec_ak = zeros(1,loop_length);
% vec_A = zeros(1,loop_length);
% vec_mu = zeros(1,loop_length);
% vec_sigma = zeros(1,loop_length);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:loop_length

%Reads the proper data file and creates variables
ak_str=num2str(ak, '%.05i');
ak_filename = sprintf('gss_low_TKM_%s.txt',ak_str);
gss_var_TKM = dlmread(ak_filename,',',0,0);
ak_varname_time_00300 = gss_var_TKM(:,1);
ak_varname_Kvals_00300 = gss_var_TKM(:,2);
ak_varname_Mvals_00300 = gss_var_TKM(:,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Fourier transform and power spectrum calculations
adjust_ak_varname_Kvals_00300 = zeros(1,numel(ak_varname_Kvals_00300));
mean_ak_varname_Kvals_00300 = mean(ak_varname_Kvals_00300);
for j = 1:numel(ak_varname_time_00300)
   adjust_ak_varname_Kvals_00300(j) = ak_varname_Kvals_00300(j)-mean_ak_varname_Kvals_00300; 
end
adjust_ak_varname_time_00300 = ak_varname_time_00300;

tops = max(ak_varname_time_00300);
fs = .01;   % sample frequency (Hz)
fourier_grid=0:1/fs:tops-1/fs;  % length of sample
vq=interp1(ak_varname_time_00300,ak_varname_Kvals_00300,fourier_grid); % interpolating data on fixed separation grid
m = numel(vq);  % window length
n = pow2(nextpow2(m));  % transform length
y = fft(vq,n);  % DFT
f = (0:n-1)*(fs/n); % frequency range - x axis
power = y.*conj(y)/n;   % power of the DFT

% Adjusts data for proper length and to ignore low frequencies
f_cutoff = 1e-4;
adjust_power_00300=power(f<f_cutoff);
adjust_power_00300(1:15)=[];
adjust_f_00300=f(f<f_cutoff);
adjust_f_00300(1:15)=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calculates the Gaussian fit. Input the guesses for A, mu, and sigma2
xvals = adjust_f_00300;
yvals = adjust_power_00300;
options = optimset('MaxIter',50000,'MaxFunEvals',50000);
Amusigma2 = fminsearch(@gaussian_fitter,[1,3e-5,9e-9],options);

%stores A,mu, and sigma in vector
% vec_ak(i) = ak/10000;
% vec_A(i) = Amusigma2(1);
% vec_mu(i) = Amusigma2(2);
% vec_sigma(i) = sqrt(Amusigma2(3));

%creates Gaussian curve
yyy_00300=zeros(1, numel(xvals));
for j=1:numel(xvals)
yyy_00300(j) = Amusigma2(1)/sqrt(2*pi*Amusigma2(3))*exp(-(xvals(j)-Amusigma2(2))^2/2/Amusigma2(3));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% = starttime:stepsize:maxtime for the deterministic
time_length_hours = 500;
timesint=0:10:time_length_hours*3600;

z0=[.08*k_k,0.01*k_m];

options = odeset('RelTol',2.22045e-14,'AbsTol',1e-16,'MaxStep',1000);
[deterministic_time,Z]=ode113(@integrator_SynEx_simple,timesint,z0,options);

%defines variable for the deterministic y coordinate values
gss_deterministic_Kvals_00300 = Z(:,1);
gss_deterministic_Mvals_00300 = Z(:,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%nullcline equations
% grid=0:.001:10;
% gss_deterministic_null1_00300=zeros(1,numel(grid));
% gss_deterministic_null2_00300=zeros(1,numel(grid));
% 
% for j=1:numel(grid)
%   gss_deterministic_null1_00300(j)=1/mu * ((a_k+b_k*grid(j)^h/(1+grid(j)^h))/grid(j)-1);
%   gss_deterministic_null2_00300(j)= a_m + b_m*grid(j)^p/(rho^p+grid(j)^p);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%adjust viewing window for deterministic and stochastic figures
% ak_varname_Kvals_00300(1:120)=[];
% ak_varname_Mvals_00300(1:120)=[];
% ak_varname_time_00300(1001:1120)=[];
% deterministic_time(1001:1100)=[];
% gss_deterministic_Kvals_00300(1:100)=[];
% gss_deterministic_Svals_00300(1:100)=[];

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
plot(deterministic_time/3600,gss_deterministic_Kvals_00300,'k','LineWidth',1)
axis([0 200 0 110 ]);
title('Deterministic','FontSize',15);
ax=gca;
xlabel('Time, \tau','FontSize',9);
ylabel('Molecules of ComK','FontSize',9);
% ax.XTick = [0, 50, 100];
% ax.YTick = [0, 50, 100];

%plots deterministic time series for MecA
% figure(3);
% plot(deterministic_time/3600,gss_deterministic_Mvals_00300,'k','LineWidth',1)
% axis([0 200 0 110 ]);
% title('Deterministic (a_k = 0.00300)','FontSize',15);
% ax=gca;
% xlabel('time \tau','FontSize',9);
% ylabel('Molecules of MecA','FontSize',9);
% ax.XTick = [0, 50, 100];
% ax.YTick = [0, 50, 100];

%Plots stochastic time series values for ComK
figure(4);
plot(ak_varname_time_00300/3600,ak_varname_Kvals_00300,'b')
axis([0 200 0 110])
title('a_k = 0.0250','FontSize',15);
ax=gca;
ylabel('Molecules of ComK','FontSize',9);
% ax.XTick = [0, 25, 50];
% ax.YTick = [0, 50, 100];

%Plots stochastic time series values for MecA
figure(5);
plot(ak_varname_time_00300/3600,ak_varname_Mvals_00300,'r')
axis([0 200 0 110])
ax=gca;
xlabel('Time, \tau','FontSize',8);
ylabel('Molecules of MecA','FontSize',9);
% ax.XTick = [0, 25, 50];
% ax.YTick = [0, 50, 100];

%Plots Gaussian fit to periodogram
figure(6);
plot(adjust_f_00300,adjust_power_00300,'b',adjust_f_00300,yyy_00300,'r');
title('Periodogram');
ax=gca;
xlabel('frequency (Hz), \tau','FontSize',9);
ylabel('Power','FontSize',9);
legend('Power spectrum','Gaussian curve fit','Location','northeast');


%loop propagation
i
if a_k >= .05000
    step_size = .000100;
end
a_k = a_k + step_size;
end

