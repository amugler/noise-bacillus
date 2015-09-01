%{ This code generates a figure, plotting time series and periodogram data for the native strain (full model) for various
    data points of a_k.
%}

%plots
figure(1);

%ComK values vs Time
subplot(3,7,1);
plot(ak_varname_time_00008/3600,ak_varname_Kvals_00008,'b')
axis([0 200 0 5e4])
ax = gca;
title('a_k = 0.00008');
% xlabel('time \tau','FontSize',8);
ylabel('Molecules of ComK','FontSize',9);
ax.XTick = [];
% ax.YTick = [0, 50, 100];

subplot(3,7,2);
plot(ak_varname_time_00038/3600,ak_varname_Kvals_00038,'b')
axis([0 200 0 5e4])
ax = gca;
title('a_k = 0.00038');
% xlabel('time \tau','FontSize',8);
ylabel('Molecules of ComK','FontSize',9);
ax.XTick = [];
% ax.YTick = [0, 50, 100];

subplot(3,7,3);
plot(ak_varname_time_00096/3600,ak_varname_Kvals_00096,'b')
axis([0 200 0 5e4])
ax = gca;
title('a_k = 0.00096');
% xlabel('time \tau','FontSize',8);
ylabel('Molecules of ComK','FontSize',9);
ax.XTick = [];
% ax.YTick = [0, 50, 100];

subplot(3,7,4);
plot(ak_varname_time_00241/3600,ak_varname_Kvals_00241,'b')
axis([0 200 0 5e4])
ax = gca;
title('a_k = 0.00241');
% xlabel('time \tau','FontSize',8);
ylabel('Molecules of ComK','FontSize',9);
ax.XTick = [];
% ax.YTick = [0, 50, 100];

subplot(3,7,5);
plot(ak_varname_time_01095/3600,ak_varname_Kvals_01095,'b')
axis([0 200 0 5e4])
ax = gca;
title('a_k = 0.01095');
% xlabel('time \tau','FontSize',8);
ylabel('Molecules of ComK','FontSize',9);
ax.XTick = [];
% ax.YTick = [0, 50, 100];

subplot(3,7,6);
plot(ak_varname_time_02190/3600,ak_varname_Kvals_02190,'b')
axis([0 200 0 5e4])
ax = gca;
title('a_k = 0.01095');
% xlabel('time \tau','FontSize',8);
ylabel('Molecules of ComK','FontSize',9);
ax.XTick = [];
% ax.YTick = [0, 50, 100];

subplot(3,7,7);
plot(ak_varname_time_10950/3600,ak_varname_Kvals_10950,'b')
axis([0 200 0 5e4])
ax = gca;
title('a_k = 0.01095');
% xlabel('time \tau','FontSize',8);
ylabel('Molecules of ComK','FontSize',9);
ax.XTick = [];
% ax.YTick = [0, 50, 100];



%ComS values vs Time
subplot(3,7,8);
plot(ak_varname_time_00008/3600,ak_varname_Svals_00008,'r')
axis([0 200 0 1.1*max(ak_varname_Svals_00008)])
ax = gca;
xlabel('Time, \tau','FontSize',8);
ylabel('Molecules of ComS','FontSize',9);
ax.XTick = [0, 50, 100];
% ax.YTick = [0, 500, 1000];

subplot(3,7,9);
plot(ak_varname_time_00038/3600,ak_varname_Svals_00038,'r')
axis([0 200 0 1.1*max(ak_varname_Svals_00038)])
ax = gca;
xlabel('Time, \tau','FontSize',8);
ylabel('Molecules of ComS','FontSize',9);
ax.XTick = [0, 50, 100];
% ax.YTick = [0, 500, 1000];

subplot(3,7,10);
plot(ak_varname_time_00096/3600,ak_varname_Svals_00096,'r')
axis([0 200 0 1.1*max(ak_varname_Svals_00096)])
ax = gca;
xlabel('Time, \tau','FontSize',8);
ylabel('Molecules of ComS','FontSize',9);
ax.XTick = [0, 50, 100];
% ax.YTick = [0, 50, 100];

subplot(3,7,11);
plot(ak_varname_time_00241/3600,ak_varname_Svals_00241,'r')
axis([0 200 0 1.1*max(ak_varname_Svals_00241)])
ax = gca;
xlabel('Time, \tau','FontSize',8);
ylabel('Molecules of ComS','FontSize',9);
ax.XTick = [0, 50, 100];
% ax.YTick = [0, 500, 1000];

subplot(3,7,12);
plot(ak_varname_time_01095/3600,ak_varname_Svals_01095,'r')
axis([0 200 0 1.1*max(ak_varname_Svals_01095)])
ax = gca;
xlabel('Time, \tau','FontSize',8);
ylabel('Molecules of ComS','FontSize',9);
ax.XTick = [0, 50, 100];
% ax.YTick = [0, 500, 1000];

subplot(3,7,13);
plot(ak_varname_time_02190/3600,ak_varname_Svals_02190,'r')
axis([0 200 0 1.1*max(ak_varname_Svals_01095)])
ax = gca;
xlabel('Time, \tau','FontSize',8);
ylabel('Molecules of ComS','FontSize',9);
ax.XTick = [0, 50, 100];
% ax.YTick = [0, 500, 1000];

subplot(3,7,14);
plot(ak_varname_time_10950/3600,ak_varname_Svals_10950,'r')
axis([0 200 0 1.1*max(ak_varname_Svals_01095)])
ax = gca;
xlabel('Time, \tau','FontSize',8);
ylabel('Molecules of ComS','FontSize',9);
ax.XTick = [0, 50, 100];
% ax.YTick = [0, 500, 1000];

%Periodograms
subplot(3,7,15);
plot(adjust_f_00008,adjust_power_00008,'b',adjust_f_00008,gaussian_curve_fit_00008,'r');
% axis([0 1e-4 0 1.1*max(adjust_power_00008)]);
title('Periodogram');
xlabel('frequency (Hz)','FontSize',9);
ylabel('Power','FontSize',9);
legend('Power spectrum','Gaussian curve fit','Location','northeast')

subplot(3,7,16);
plot(adjust_f_00038,adjust_power_00038,'b',adjust_f_00038,gaussian_curve_fit_00038,'r');
% axis([0 1e-4 0 1.1*max(adjust_power_00038)]);
title('Periodogram');
xlabel('frequency (Hz)','FontSize',9);
ylabel('Power','FontSize',9);
% legend('Power spectrum','Gaussian curve fit','Location','northeast')

subplot(3,7,17);
plot(adjust_f_00096,adjust_power_00096,'b',adjust_f_00096,gaussian_curve_fit_00096,'r');
% axis([0 1e-4 0 1.1*max(adjust_power_00096)]);
title('Periodogram');
xlabel('frequency (Hz)','FontSize',9);
ylabel('Power','FontSize',9);
% legend('Power spectrum','Gaussian curve fit','Location','northeast')

subplot(3,7,18);
plot(adjust_f_00241,adjust_power_00241,'b',adjust_f_00241,gaussian_curve_fit_00241,'r');
axis([0 1e-4 0 1.1*max(adjust_power_00241)]);
title('Periodogram');
xlabel('frequency (Hz)','FontSize',9);
ylabel('Power','FontSize',9);
% legend('Power spectrum','Gaussian curve fit','Location','northeast')

subplot(3,7,19);
plot(adjust_f_01095,adjust_power_01095,'b',adjust_f_01095,gaussian_curve_fit_01095,'r');
% axis([0 3e-6 0 1.1*max(adjust_power_01095)]);
title('Periodogram');
xlabel('frequency (Hz)','FontSize',9);
ylabel('Power','FontSize',9);
% legend('Power spectrum','Gaussian curve fit','Location','northeast')

subplot(3,7,20);
plot(adjust_f_02190,adjust_power_02190,'b',adjust_f_02190,gaussian_curve_fit_02190,'r');
% axis([0 3e-6 0 1.1*max(adjust_power_02190)]);
title('Periodogram');
xlabel('frequency (Hz)','FontSize',9);
ylabel('Power','FontSize',9);
% legend('Power spectrum','Gaussian curve fit','Location','northeast')

subplot(3,7,21);
plot(adjust_f_10950,adjust_power_10950,'b',adjust_f_10950,gaussian_curve_fit_10950,'r');
% axis([0 3e-6 0 1.1*max(adjust_power_10950)]);
title('Periodogram');
xlabel('frequency (Hz)','FontSize',9);
ylabel('Power','FontSize',9);
% legend('Power spectrum','Gaussian curve fit','Location','northeast')



