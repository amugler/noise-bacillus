%{ 
    This code creates stochastic and deterministic time series, nullclines, and power spectra plots for the SynEx strain
    (full model, high molecule number).
%}

global xvals yvals alpha_k alpha_m beta_k beta_m k_k k_m p h delta lambda
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k1 = 0.0036;
k1_st = 00036;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%loop data
step_size = .00010;
loop_length = 1;


%general data for deterministic plotter
h = 2;
p = 2;
k_k = 5000;
k_m = 2500;
lambda = .0001;
a_k = k1/.0125;
alpha_k = a_k*k_k*lambda;
alpha_m = .075;
beta_k = 7.5;
beta_m = 2.5;
delta = .00000004;

% vec_ak = zeros(1,loop_length);
% vec_A = zeros(1,loop_length);
% vec_mu = zeros(1,loop_length);
% vec_sigma = zeros(1,loop_length);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:loop_length

%Reads the proper data file and creates variables
k1_str=num2str(k1_st, '%.05i');
ak_filename = sprintf('gsf_TKM_%s.txt',k1_str);
gsf_var_TKM = dlmread(ak_filename,',',0,0);
ak_varname_time_00036 = gsf_var_TKM(:,1);
ak_varname_Kvals_00036 = gsf_var_TKM(:,2);
ak_varname_Mvals_00036 = gsf_var_TKM(:,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Fourier transform and power spectrum calculations
adjust_ak_varname_Kvals_00036 = zeros(1,numel(ak_varname_Kvals_00036));
mean_ak_varname_Kvals_00036 = mean(ak_varname_Kvals_00036);
for j = 1:numel(ak_varname_time_00036)
   adjust_ak_varname_Kvals_00036(j) = ak_varname_Kvals_00036(j)-mean_ak_varname_Kvals_00036; 
end
adjust_ak_varname_time_00036 = ak_varname_time_00036;

tops = max(adjust_ak_varname_time_00036);
fs = .01;   % sample frequency (Hz)
fourier_grid=0:1/fs:tops-1/fs;  % length of sample
vq=interp1(adjust_ak_varname_time_00036,adjust_ak_varname_Kvals_00036,fourier_grid); % interpolating data on fixed separation grid
m = numel(vq);  % window length
n = pow2(nextpow2(m));  % transform length
y = fft(vq,n);  % DFT
f = (0:n-1)*(fs/n); % frequency range - x axis
power = y.*conj(y)/n;   % power of the DFT

%Adjusts data for proper length and to ignore low frequencies
f_cutoff = 1e-4;
adjust_power_00036=power(f<f_cutoff);
adjust_power_00036(1:15)=[];
adjust_f_00036=f(f<f_cutoff);
adjust_f_00036(1:15)=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calculates the Gaussian fit. Input the guesses for A, mu, and sigma2
xvals = adjust_f_00036;
yvals = adjust_power_00036;
options = optimset('MaxIter',50000,'MaxFunEvals',50000);
Amusigma2 = fminsearch(@gaussian_fitter,[2e8,-1e-6,1e-10],options);

%creates Gaussian curve
yyy_00036=zeros(1, numel(xvals));
for j=1:numel(xvals)
yyy_00036(j) = Amusigma2(1)/sqrt(2*pi*Amusigma2(3))*exp(-(xvals(j)-Amusigma2(2))^2/2/Amusigma2(3));
end

%stores A,mu, and sigma in vector
% vec_ak(i) = ak/10000;
% vec_A(i) = Amusigma2(1);
% vec_mu(i) = Amusigma2(2);
% vec_sigma(i) = sqrt(Amusigma2(3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% = starttime:stepsize:maxtime for the deterministic
time_length_hours = 500;
timesint=0:1000:time_length_hours*3600;

z0=[500,150];

options = odeset('RelTol',2.22045e-14,'AbsTol',1e-16,'MaxStep',1000);
[deterministic_time,Z]=ode113(@integrator_SynEx_simple,timesint,z0,options);

%defines variable for the deterministic y coordinate values
gsf_deterministic_Kvals_00036 = Z(:,1);
gsf_deterministic_Mvals_00036 = Z(:,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %nullcline equations
% grid=1:.1:20000;
% gsf_deterministic_null1_00036=zeros(1,numel(grid));
% gsf_deterministic_null2_00036=zeros(1,numel(grid));
% 
% for j=1:numel(grid)
%   gsf_deterministic_null1_00036(j)= (alpha_k+beta_k*grid(j)^h/(k_k^h+grid(j)^h)-lambda *grid(j))/(delta*grid(j));
%   gsf_deterministic_null2_00036(j)= (alpha_m + beta_m*grid(j)^p/(k_m^p+grid(j)^p))/lambda;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%adjust viewing window for deterministic and stochastic figures
% ak_varname_Kvals_00036(1:120)=[];
% ak_varname_Mvals_00036(1:120)=[];
% ak_varname_time_00036(1001:1120)=[];
% deterministic_time(1001:1100)=[];
% gsf_deterministic_Kvals_00036(1:100)=[];
% gsf_deterministic_Svals_00036(1:100)=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %Plots nullclines
% figure(1);
% loglog((grid),gsf_deterministic_null1_00036,'r',(grid),gsf_deterministic_null2_00036,'g',gsf_deterministic_Kvals_00036,gsf_deterministic_Mvals_00036,'k',ak_varname_Kvals_00036,ak_varname_Mvals_00036,'b');
% % axis([.2*10e-3 10e0 10e-3 5*10e1]);
% % title('a_k = 0.01');
% ax=gca;
% xlabel('Kbar','FontSize',8);
% ylabel('Mbar','FontSize',8);
% % ax.XTick = [10e-3,10e-2,10e-1,10e0];
% % ax.YTick = [10e-3,10e-1, 10e1];

%plots deterministic time series for ComK
% figure(2);
% plot(deterministic_time/3600,gsf_deterministic_Kvals_00036,'k','LineWidth',1)
% axis([0 200 0 2e4 ]);
% title('Deterministic ComK (a_k = 0.8000)','FontSize',15);
% ax=gca;
% xlabel('Time, \tau','FontSize',9);
% ylabel('Molecules of ComK','FontSize',9);
% ax.XTick = [0, 50, 100];
% ax.YTick = [0, 50, 100];

%plots deterministic time series for MecA
% figure(3);
% plot(deterministic_time/3600,gsf_deterministic_Mvals_00036,'k','LineWidth',1)
% axis([0 200 0 5e4 ]);
% title('Deterministic MecA (a_k = 0.8000)','FontSize',15);
% ax=gca;
% xlabel('time \tau','FontSize',9);
% ylabel('Molecules of MecA','FontSize',9);
% ax.XTick = [0, 50, 100];
% ax.YTick = [0, 50, 100];

%Plots stochastic time series values for ComK
figure(4);
plot(ak_varname_time_00036/3600,ak_varname_Kvals_00036,'b')
axis([0 200 0 2e4])
% title('a_k = 0.8000','FontSize',15);
ax=gca;
ylabel('Molecules of ComK','FontSize',9);
% ax.XTick = [0, 25, 50];
% ax.YTick = [0, 50, 100];

%Plots stochastic time series values for MecA
figure(5);
plot(ak_varname_time_00036/3600,ak_varname_Mvals_00036,'r')
axis([0 200 0 5e4])
ax=gca;
xlabel('Time, \tau','FontSize',8);
ylabel('Molecules of MecA','FontSize',9);
% ax.XTick = [0, 25, 50];
% ax.YTick = [0, 50, 100];

%Plots Gaussian fit to periodogram
figure(6);
plot(adjust_f_00036,adjust_power_00036,'b',adjust_f_00036,yyy_00036,'r');
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

