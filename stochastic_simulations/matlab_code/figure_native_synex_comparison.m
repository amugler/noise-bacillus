%{ 
    This code generates a figure comparing the time series data for the native and synex strains
    (simple models, low molecule numbers).
%}


figure(1);

subplot(3,4,1);
plot(ak_varname_time_00065/3600,ak_varname_Kvals_00065,'b','LineWidth',1.5);
hold on;
plot(deterministic_time/3600,gns_deterministic_Kvals_00065,'k','LineWidth',2);
axis([0 100 0 115])
ax = gca;
title('a_k = ');
% xlabel('time \tau','FontSize',8);
ylabel('Molecules of ComK','FontSize',9);
ax.XTick = [];
ax.YTick = [0, 50, 100];
legend('Stochastic','Deterministic','Location','northeast')


subplot(3,4,2);
plot(ak_varname_time_00300/3600,ak_varname_Kvals_00300,'r','LineWidth',1.5);
hold on;
plot(deterministic_time/3600,gss_deterministic_Kvals_00300,'k','LineWidth',2);
axis([0 100 0 115])
ax = gca;
title('\alpha_k = ');
% xlabel('time \tau','FontSize',8);
ylabel('Molecules of ComK','FontSize',9);
ax.XTick = [];
ax.YTick = [0, 50, 100];
legend('Stochastic','Deterministic','Location','northeast')

subplot(3,4,3);
plot(ak_varname_time_01000/3600,ak_varname_Kvals_01000,'b','LineWidth',1.5);
hold on;
plot(deterministic_time/3600,gns_deterministic_Kvals_01000,'k','LineWidth',2);
axis([0 100 0 115])
ax = gca;
title('a_k = ');
% xlabel('time \tau','FontSize',8);
ylabel('Molecules of ComK','FontSize',9);
ax.XTick = [];
ax.YTick = [0, 50, 100];

subplot(3,4,4);
plot(ak_varname_time_11000/3600,ak_varname_Kvals_11000,'r','LineWidth',1.5);
hold on;
plot(deterministic_time/3600,gss_deterministic_Kvals_11000,'k','LineWidth',2);
axis([0 100 0 57.5])
ax = gca;
title('a_k = ');
% xlabel('time \tau','FontSize',8);
ylabel('Molecules of ComK','FontSize',9);
ax.XTick = [];
ax.YTick = [0, 25, 50];




subplot(3,4,5);
plot(ak_varname_time_00065/3600,ak_varname_Svals_00065,'g','LineWidth',1.5);
hold on;
plot(deterministic_time/3600,gns_deterministic_Svals_00065,'k','LineWidth',2);
axis([0 100 0 115])
ax = gca;
xlabel('Time (hours)','FontSize',8);
ylabel('Molecules of ComS','FontSize',9);
ax.XTick = [0, 50, 100];
ax.YTick = [0, 50, 100];
legend('Stochastic','Deterministic','Location','northeast')

subplot(3,4,6);
plot(ak_varname_time_00300/3600,ak_varname_Mvals_00300,'r','LineWidth',1.5);
hold on;
plot(deterministic_time/3600,gss_deterministic_Mvals_00300,'k','LineWidth',2);
axis([0 100 0 115])
ax = gca;
xlabel('Time (hours)','FontSize',8);
ylabel('Molecules of MecA','FontSize',9);
ax.XTick = [0, 50, 100];
ax.YTick = [0, 50, 100];
legend('Stochastic','Deterministic','Location','northeast')

subplot(3,4,7);
plot(ak_varname_time_01000/3600,ak_varname_Svals_01000,'g','LineWidth',1.5);
hold on;
plot(deterministic_time/3600,gns_deterministic_Svals_01000,'k','LineWidth',2);
axis([0 100 0 23])
ax = gca;
xlabel('Time (hours)','FontSize',8);
ylabel('Molecules of ComS','FontSize',9);
ax.XTick = [0, 50, 100];
ax.YTick = [0, 10, 20];

subplot(3,4,8);
plot(ak_varname_time_11000/3600,ak_varname_Mvals_11000,'r','LineWidth',1.5);
hold on;
plot(deterministic_time/3600,gss_deterministic_Mvals_11000,'k','LineWidth',2);
axis([0 100 0 57.5])
ax = gca;
xlabel('Time (hours)','FontSize',8);
ylabel('Molecules of MecA','FontSize',9);
ax.XTick = [0, 50, 100];
ax.YTick = [0, 25, 50];





subplot(3,4,9);
plot(adjust_f_00065,adjust_power_00065,'b',adjust_f_00065,yyy_00065,'r');
% axis([0 1e-4 0 1.1*max(adjust_power_00065)]);
xlabel('frequency (Hz)','FontSize',9);
ylabel('Power','FontSize',9);
legend('Power spectrum','Gaussian fit','Location','northeast')

subplot(3,4,10);
plot(adjust_f_00300,adjust_power_00300,'b',adjust_f_00300,yyy_00300,'r');
% axis([0 1e-4 0 1.1*max(adjust_power_00300)]);
xlabel('frequency (Hz)','FontSize',9);
ylabel('Power','FontSize',9);

subplot(3,4,11);
plot(adjust_f_01000,adjust_power_01000,'b',adjust_f_01000,yyy_01000,'r');
% axis([0 3e-6 0 1.1*max(adjust_power_01000)]);
xlabel('frequency (Hz)','FontSize',9);
ylabel('Power','FontSize',9);

subplot(3,4,12);
plot(adjust_f_11000,adjust_power_11000,'b',adjust_f_11000,yyy_11000,'r');
% axis([0 3e-6 0 1.1*max(adjust_power_11000)]);
xlabel('frequency (Hz)','FontSize',9);
ylabel('Power','FontSize',9);
