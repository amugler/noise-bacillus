%{ 
    This code creates stochastic and deterministic time series, nullclines, and power spectra plots for the native strain
    (full model, high molecule number).
%}

global xvals yvals alpha_k alpha_s beta_k beta_s k_k k_s p h delta lambda_k lambda_s Gamma_k Gamma_s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k1 = .00060;
k1_st = 00060;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%general data for deterministic plotter
Gamma_k = 25000;
Gamma_s = 20;
delta = .001;
h = 2;
p = 5;
a_k = k1/.625;
alpha_k = a_k*Gamma_k*delta;
alpha_s = .0004;
beta_k = 7.5;
beta_s = .06;
k_k = 5000;
k_s = 833;
lambda_k = .0001;
lambda_s = .0001;


% vec_ak = zeros(1,loop_length);
% vec_A = zeros(1,loop_length);
% vec_mu = zeros(1,loop_length);
% vec_sigma = zeros(1,loop_length);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Reads the proper data file and creates variables
k1_str=num2str(k1_st, '%.05i');
ak_filename = sprintf('gnf_TKS_%s.txt',k1_str);
gnf_var_TKS = dlmread(ak_filename,',',0,0);
ak_varname_time_00060 = gnf_var_TKS(:,1);
ak_varname_Kvals_00060 = gnf_var_TKS(:,2);
ak_varname_Svals_00060 = gnf_var_TKS(:,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Fourier transform and power spectrum calculations
adjust_ak_varname_Kvals_00060 = zeros(1,numel(ak_varname_Kvals_00060));
mean_ak_varname_Kvals_00060 = mean(ak_varname_Kvals_00060);
for j = 1:numel(ak_varname_time_00060)
   adjust_ak_varname_Kvals_00060(j) = ak_varname_Kvals_00060(j)-mean_ak_varname_Kvals_00060; 
end
adjust_ak_varname_time_00060 = ak_varname_time_00060;

tops = max(adjust_ak_varname_time_00060);
fs = .01;   % sample frequency (Hz)
fourier_grid=0:1/fs:tops-1/fs;  % length of sample
vq=interp1(adjust_ak_varname_time_00060,adjust_ak_varname_Kvals_00060,fourier_grid); % interpolating data on fixed separation grid
m = numel(vq);  % window length
n = pow2(nextpow2(m));  % transform length
y = fft(vq,n);  % DFT
f = (0:n-1)*(fs/n); % frequency range - x axis
power = y.*conj(y)/n;   % power of the DFT

%Adjusts data for proper length and to ignore low frequencies
f_cutoff = 1e-4;
adjust_power_00060=power(f<f_cutoff);
adjust_power_00060(1:15)=[];
adjust_f_00060=f(f<f_cutoff);
adjust_f_00060(1:15)=[];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calculates the Gaussian fit. Input the guesses for A, mu, and sigma2
xvals = adjust_f_00060;
yvals = adjust_power_00060;
options = optimset('MaxIter',50000,'MaxFunEvals',50000);
Amusigma2 = fminsearch(@gaussian_fitter,[1e7,1e-5,1e-9],options);

%stores A,mu, and sigma in vector
% vec_ak(i) = ak/100000;
% vec_A(i) = Amusigma2(1);
% vec_mu(i) = Amusigma2(2);
% vec_sigma(i) = sqrt(Amusigma2(3));

%creates Gaussian curve
gaussian_curve_fit_00060=zeros(1, numel(xvals));
for j=1:numel(xvals)
gaussian_curve_fit_00060(j) = Amusigma2(1)/sqrt(2*pi*Amusigma2(3))*exp(-(xvals(j)-Amusigma2(2))^2/2/Amusigma2(3));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% = starttime:stepsize:maxtime for the deterministic
% time_length_hours = 500;
% timesint=0:1000:time_length_hours*3600;
% 
% z0=[350,7000];
% 
% options = odeset('RelTol',2.22045e-14,'AbsTol',1e-16,'MaxStep',00060);
% [deterministic_time,Z]=ode113(@integrator_native_simple,timesint,z0,options);
% 
% %defines variable for the deterministic y coordinate values
% gnf_deterministic_Kvals_00060 = Z(:,1);
% gnf_deterministic_Svals_00060 = Z(:,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%nullcline equations
% grid=0:.0001:2;
% gnf_deterministic_null1_00060=zeros(1,numel(grid));
% gnf_deterministic_null2_00060=zeros(1,numel(grid));
% 
% for j=1:numel(grid)
%   gnf_deterministic_null1_00060(j)=((a_k + (b_k)*grid(j)^h/(k_0^h+grid(j)^h))/grid(j)-Delta_k)^(-1) - grid(j)-1;
%   gnf_deterministic_null2_00060(j)= 1/(2*Delta_s)*( ((1+(1+grid(j))*Delta_s-(a_s + (b_s)/(1+(grid(j)/k_1)^p)))^2+4*Delta_s*((1+grid(j))*(a_s + (b_s)/(1+(grid(j)/k_1)^p))))^(1/2)- (1+(1+grid(j))*Delta_s- (a_s + (b_s)/(1+(grid(j)/k_1)^p))));
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%adjust viewing window for deterministic and stochastic figures
% ak_varname_Kvals_00060(1:4000)=[];
% ak_varname_Svals_00060(1:4000)=[];
% ak_varname_time_00060(000601:104000)=[];
% deterministic_time(1001:1100)=[];
% gnf_deterministic_Kvals_00060(1:100)=[];
% gnf_deterministic_Svals_00060(1:100)=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Plots nullclines
% figure(1);
% loglog((grid),gnf_deterministic_null1_00060,(grid),gnf_deterministic_null2_00060,gnf_deterministic_Kvals_00060/100,gnf_deterministic_Svals_00060,'k');
% axis([.2*10e-3 2*10e-1 10e-5 2*10e1]);
% title('a_k = 0.00060');
% ax=gca;
% xlabel('Kbar','FontSize',8);
% ylabel('Sbar','FontSize',8);
% ax.XTick = [10e-3,10e-2,10e-1,10e0];
% ax.YTick = [10e-3,10e-1, 10e1];


%plots deterministic time series for ComK
% figure(2);
% plot(deterministic_time/3600,gnf_deterministic_Kvals_00060,'k')
% axis([0 200 0 110])
% title('Deterministic','FontSize',15);
% ax=gca;
% ylabel('Molecules of ComK','FontSize',9);
% ax.XTick = [0, 50, 100];
% ax.YTick = [0, 50, 100];


%plots deterministic time series for ComS
% figure(3);
% plot(deterministic_time/3600,gnf_deterministic_Svals_00060,'k')
% axis([0 200 0 110])
% title('Deterministic','FontSize',15);
% ax=gca;
% ylabel('Molecules of ComS','FontSize',9);
% ax.XTick = [0, 50, 100];
% ax.YTick = [0, 50, 100];


%Plots stochastic time series values for ComK
figure(4);
plot(ak_varname_time_00060/3600,ak_varname_Kvals_00060,'b')
axis([0 200 0 1.05*max(ak_varname_Kvals_00060)])
title('Stochastic','FontSize',15);
ax=gca;
ylabel('Molecules of ComK','FontSize',9);
% ax.XTick = [0, 25, 50];
% ax.YTick = [0, 50, 100];


%Plots stochastic time series values for ComS
figure(5);
plot(ak_varname_time_00060/3600,ak_varname_Svals_00060,'r')
axis([0 200 0 1.05*max(ak_varname_Svals_00060)])
ax=gca;
xlabel('Time, \tau','FontSize',8);
ylabel('Molecules of ComS','FontSize',9);
% ax.XTick = [0, 25, 50];
% ax.YTick = [0, 50, 100];


%Plots Gaussian fit to periodogram
figure(6);
plot(adjust_f_00060,adjust_power_00060,'b',adjust_f_00060,gaussian_curve_fit_00060,'r');
title('Periodogram');
ax=gca;
xlabel('frequency (Hz), \tau','FontSize',9);
ylabel('Power','FontSize',9);
legend('Power spectrum','Gaussian curve fit','Location','northeast')

