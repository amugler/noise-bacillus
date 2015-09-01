%{ 
    This code generates a figure of the time series data for the SynEx full model (high molecule number).
%}

%plots
figure(1);

%ComK values vs Time
subplot(3,5,1);
plot(ak_varname_time_00008/3600,ak_varname_Kvals_00008/10^4,'b','LineWidth',1.5)
axis([0 200 0 2.3])
ax = gca;
title('a_k = ');
ylabel('ComK (count/10^4)','FontSize',8);
ax.XTick = [];
ax.YTick = [0, 1, 2];

subplot(3,5,2);
plot(ak_varname_time_00036/3600,ak_varname_Kvals_00036/10^4,'b','LineWidth',1.5)
axis([0 200 0 2.3])
ax = gca;
title('a_k = ');
ylabel('ComK (count/10^4)','FontSize',8);
ax.XTick = [];
ax.YTick = [0, 1, 2];

subplot(3,5,3);
plot(ak_varname_time_00197/3600,ak_varname_Kvals_00197/10^4,'b','LineWidth',1.5)
axis([0 200 0 2.3])
ax = gca;
title('a_k = ');
ylabel('ComK (count/10^4)','FontSize',8);
ax.XTick = [];
ax.YTick = [0, 1, 2];

subplot(3,5,4);
plot(ak_varname_time_01462/3600,ak_varname_Kvals_01462/10^4,'b','LineWidth',1.5)
axis([0 200 0 2.3])
ax = gca;
title('a_k = ');
ylabel('ComK (count/10^4)','FontSize',8);
ax.XTick = [];
ax.YTick = [0, 1, 2];

subplot(3,5,5);
plot(ak_varname_time_04880/3600,ak_varname_Kvals_04880/10^4,'b','LineWidth',1.5)
axis([0 200 0 2.3])
ax = gca;
title('a_k = ');
ylabel('ComK (count/10^4)','FontSize',8);
ax.XTick = [];
ax.YTick = [0, 1, 2];


%MecA values vs Time
subplot(3,5,6);
plot(ak_varname_time_00008/3600,ak_varname_Mvals_00008/10^4,'g','LineWidth',1.5)
axis([0 200 0 4.6])
ax = gca;
xlabel('Time, \tau','FontSize',8);
ylabel('MecA (count/10^4)','FontSize',8);
ax.XTick = [0, 100, 200];
ax.YTick = [0, 2, 4];

subplot(3,5,7);
plot(ak_varname_time_00036/3600,ak_varname_Mvals_00036/10^4,'g','LineWidth',1.5)
axis([0 200 0 4.6])
ax = gca;
xlabel('Time, \tau','FontSize',8);
ylabel('MecA (count/10^4)','FontSize',8);
ax.XTick = [0, 100, 200];
ax.YTick = [0, 2, 4];

subplot(3,5,8);
plot(ak_varname_time_00197/3600,ak_varname_Mvals_00197/10^4,'g','LineWidth',1.5)
axis([0 200 0 4.6])
ax = gca;
xlabel('Time, \tau','FontSize',8);
ylabel('MecA (count/10^4)','FontSize',8);
ax.XTick = [0, 100, 200];
ax.YTick = [0, 2, 4];

subplot(3,5,9);
plot(ak_varname_time_01462/3600,ak_varname_Mvals_01462/10^4,'g','LineWidth',1.5)
axis([0 200 0 4.6])
ax = gca;
xlabel('Time, \tau','FontSize',8);
ylabel('MecA (count/10^4)','FontSize',8);
ax.XTick = [0, 100, 200];
ax.YTick = [0, 2, 4];

subplot(3,5,10);
plot(ak_varname_time_04880/3600,ak_varname_Mvals_04880/10^4,'g','LineWidth',1.5)
axis([0 200 0 4.6])
ax = gca;
xlabel('Time, \tau','FontSize',8);
ylabel('MecA (count/10^4)','FontSize',8);
ax.XTick = [0, 100, 200];
ax.YTick = [0, 2, 4];


%Periodograms
subplot(3,5,11);
plot(adjust_f_00008,adjust_power_00008,'b',adjust_f_00008,yyy_00008,'r');
% axis([0 1e-4 0 1.1*max(adjust_power_00008)]);
xlabel('frequency (Hz)','FontSize',9);
ylabel('Power','FontSize',9);
legend('Power spectrum','Gaussian curve fit','Location','northeast')

subplot(3,5,12);
plot(adjust_f_00036,adjust_power_00036,'b',adjust_f_00036,yyy_00036,'r');
xlabel('frequency (Hz)','FontSize',9);
ylabel('Power','FontSize',9);

subplot(3,5,13);
plot(adjust_f_00197,adjust_power_00197,'b',adjust_f_00197,yyy_00197,'r');
xlabel('frequency (Hz)','FontSize',9);
ylabel('Power','FontSize',9);

subplot(3,5,14);
plot(adjust_f_01462,adjust_power_01462,'b',adjust_f_01462,yyy_01462,'r');
xlabel('frequency (Hz)','FontSize',9);
ylabel('Power','FontSize',9);

subplot(3,5,15);
plot(adjust_f_04880,adjust_power_04880,'b',adjust_f_04880,yyy_04880,'r');
xlabel('frequency (Hz)','FontSize',9);
ylabel('Power','FontSize',9);



