%{ 
    This code generates a figure of the time series data and power spectra for the native full model (high molecule number).
%}

%plots
figure(1);

%ComK values vs Time
subplot(3,5,1);
plot(ak_varname_time_00008/3600,ak_varname_Kvals_00008/10^4,'b','LineWidth',1.5)
axis([0 200 0 4.6])
ax = gca;
title('a_k = ');
ylabel('ComK (count/10^4)','FontSize',8);
ax.XTick = [];
ax.YTick = [0, 2, 4];

subplot(3,5,2);
plot(ak_varname_time_00060/3600,ak_varname_Kvals_00060/10^4,'b','LineWidth',1.5)
axis([0 200 0 4.6])
ax = gca;
title('a_k = ');
ylabel('ComK (count/10^4)','FontSize',8);
ax.XTick = [];
ax.YTick = [0, 2, 4];

subplot(3,5,3);
plot(ak_varname_time_00096/3600,ak_varname_Kvals_00096/10^4,'b','LineWidth',1.5)
axis([0 200 0 4.6])
ax = gca;
title('a_k = ');
ylabel('ComK (count/10^4)','FontSize',8);
ax.XTick = [];
ax.YTick = [0, 2, 4];

subplot(3,5,4);
plot(ak_varname_time_01095/3600,ak_varname_Kvals_01095/10^4,'b','LineWidth',1.5)
axis([0 200 0 4.6])
ax = gca;
title('a_k = ');
ylabel('ComK (count/10^4)','FontSize',8);
ax.XTick = [];
ax.YTick = [0, 2, 4];

subplot(3,5,5);
plot(ak_varname_time_10950/3600,ak_varname_Kvals_10950/10^4,'b','LineWidth',1.5)
axis([0 200 0 4.6])
ax = gca;
title('a_k = ');
ylabel('ComK (count/10^4)','FontSize',8);
ax.XTick = [];
ax.YTick = [0, 2, 4];


%ComS values vs Time
subplot(3,5,6);
plot(ak_varname_time_00008/3600,ak_varname_Svals_00008/100,'g','LineWidth',1.5)
axis([0 200 0 11.5])
ax = gca;
xlabel('Time, \tau','FontSize',8);
ylabel('ComS (count/10^2)','FontSize',8);
ax.XTick = [0, 100, 200];
ax.YTick = [0, 5, 10];

subplot(3,5,7);
plot(ak_varname_time_00060/3600,ak_varname_Svals_00060/100,'g','LineWidth',1.5)
axis([0 200 0 11.5])
ax = gca;
xlabel('Time, \tau','FontSize',8);
ylabel('ComS (count/10^2)','FontSize',8);
ax.XTick = [0, 100, 200];
ax.YTick = [0, 5, 10];

subplot(3,5,8);
plot(ak_varname_time_00096/3600,ak_varname_Svals_00096/100,'g','LineWidth',1.5)
axis([0 200 0 11.5])
ax = gca;
xlabel('Time, \tau','FontSize',8);
ylabel('ComS (count/10^2)','FontSize',8);
ax.XTick = [0, 100, 200];
ax.YTick = [0, 5, 10];

subplot(3,5,9);
plot(ak_varname_time_01095/3600,ak_varname_Svals_01095/100,'g','LineWidth',1.5)
axis([0 200 0 4.6])
ax = gca;
xlabel('Time, \tau','FontSize',8);
ylabel('ComS 0(count/10^2)','FontSize',8);
ax.XTick = [0, 100, 200];
ax.YTick = [0, 2, 4];

subplot(3,5,10);
plot(ak_varname_time_10950/3600,ak_varname_Svals_10950/100,'g','LineWidth',1.5)
axis([0 200 0 2.3])
ax = gca;
xlabel('Time, \tau','FontSize',8);
ylabel('ComS (count/10^2)','FontSize',8);
ax.XTick = [0, 100, 200];
ax.YTick = [0, 1, 2];


%Periodograms
subplot(3,5,11);
plot(adjust_f_00008,adjust_power_00008,'b',adjust_f_00008,gaussian_curve_fit_00008,'r');
% axis([0 1e-4 0 1.1*max(adjust_power_00008)]);
xlabel('frequency (Hz)','FontSize',9);
ylabel('Power','FontSize',9);
legend('Power spectrum','Gaussian curve fit','Location','northeast')

subplot(3,5,12);
plot(adjust_f_00060,adjust_power_00060,'b',adjust_f_00060,gaussian_curve_fit_00060,'r');
% axis([0 1e-4 0 1.1*max(adjust_power_00008)]);
xlabel('frequency (Hz)','FontSize',9);
ylabel('Power','FontSize',9);

subplot(3,5,13);
plot(adjust_f_00096,adjust_power_00096,'b',adjust_f_00096,gaussian_curve_fit_00096,'r');
% axis([0 1e-4 0 1.1*max(adjust_power_00096)]);
xlabel('frequency (Hz)','FontSize',9);
ylabel('Power','FontSize',9);

subplot(3,5,14);
plot(adjust_f_01095,adjust_power_01095,'b',adjust_f_01095,gaussian_curve_fit_01095,'r');
% axis([0 3e-6 0 1.1*max(adjust_power_01095)]);
xlabel('frequency (Hz)','FontSize',9);
ylabel('Power','FontSize',9);

subplot(3,5,15);
plot(adjust_f_10950,adjust_power_10950,'b',adjust_f_10950,gaussian_curve_fit_10950,'r');
% axis([0 3e-6 0 1.1*max(adjust_power_10950)]);
xlabel('frequency (Hz)','FontSize',9);
ylabel('Power','FontSize',9);



