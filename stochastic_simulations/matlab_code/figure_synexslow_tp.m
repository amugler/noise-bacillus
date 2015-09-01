%{ 
    This code generates a figure of the time series and power spectra data for the SynExSlow simple model
    (high molecule number).
%}

%plots
figure(1);

%Deterministic ComK values vs Time
subplot(5,5,1);
plot(deterministic_time/3600,gssslow_deterministic_Kvals_00476/1e4,'k')
axis([0 200 0 4.6])
ax = gca;
title('a_k = .2 x a_k^{(1)} ');
ax = gca;
xlabel('Time (hours)','FontSize',8);
ylabel('ComK (count/10^4)','FontSize',8);
ax.XTick = [0, 100, 200];
ax.YTick = [0, 2, 4];

subplot(5,5,2);
plot(deterministic_time/3600,gssslow_deterministic_Kvals_02142/1e4,'k')
axis([0 200 0 4.6])
ax = gca;
title('a_k = .2 x a_k^{(1)} ');
ax = gca;
xlabel('Time (hours)','FontSize',8);
ylabel('ComK (count/10^4)','FontSize',8);
ax.XTick = [0, 100, 200];
ax.YTick = [0, 2, 4];

subplot(5,5,3);
plot(deterministic_time/3600,gssslow_deterministic_Kvals_08450/1e4,'k')
axis([0 200 0 4.6])
ax = gca;
title('a_k = ');
xlabel('Time (hours)','FontSize',8);
ylabel('ComK (count/10^4)','FontSize',8);
ax.XTick = [0, 100, 200];
ax.YTick = [0, 2, 4];

subplot(5,5,4);
plot(deterministic_time/3600,gssslow_deterministic_Kvals_60000/1e4,'k')
axis([0 200 0 4.6])
ax = gca;
title('a_k = 2 x a_k^{(2)} ');
xlabel('Time (hours)','FontSize',8);
ylabel('ComK (count/10^4)','FontSize',8);
ax.XTick = [0, 100, 200];
ax.YTick = [0, 2, 4];

subplot(5,5,5);
plot(deterministic_time/3600,gssslow_deterministic_Kvals_150000/1e4,'k')
axis([0 200 0 4.6])
ax = gca;
title('a_k = 50 x a_k^{(2)}');
xlabel('Time (hours)','FontSize',8);
ylabel('ComK (count/10^4)','FontSize',8);
ax.XTick = [0, 100, 200];
ax.YTick = [0, 2, 4];


%Stochastic ComK values vs Time
subplot(5,5,6);
plot(ak_varname_time_00476/3600,ak_varname_Kvals_00476/1e4,'b')
axis([0 200 0 4.6])
ax = gca;
ylabel('ComK (count/10^4)','FontSize',8);
ax.XTick = [];
ax.YTick = [0, 2, 4];

subplot(5,5,7);
plot(ak_varname_time_02142/3600,ak_varname_Kvals_02142/1e4,'b')
axis([0 200 0 4.6])
ax = gca;
ylabel('ComK (count/10^4)','FontSize',8);
ax.XTick = [];
ax.YTick = [0, 2, 4];

subplot(5,5,8);
plot(ak_varname_time_08450/3600,ak_varname_Kvals_08450/1e4,'b')
axis([0 200 0 4.6])
ax = gca;
ylabel('ComK (count/10^4)','FontSize',8);
ax.XTick = [];
ax.YTick = [0, 2, 4];

subplot(5,5,9);
plot(ak_varname_time_60000/3600,ak_varname_Kvals_60000/1e4,'b')
axis([0 200 0 4.6])
ax = gca;
ylabel('ComK (count/10^4)','FontSize',8);
ax.XTick = [];
ax.YTick = [0, 2, 4];

subplot(5,5,10);
plot(ak_varname_time_150000/3600,ak_varname_Kvals_150000/1e4,'b')
axis([0 200 0 4.6])
ax = gca;
ylabel('ComK (count/10^4)','FontSize',8);
ax.XTick = [];
ax.YTick = [0, 2, 4];


%MecA values vs Time
subplot(5,5,11);
plot(ak_varname_time_00476/3600,ak_varname_Mvals_00476/1e4,'r')
axis([0 200 0 4.6])
ax = gca;
ylabel('MecA (count/10^4)','FontSize',8);
ax.XTick = [];
ax.YTick = [0, 2, 4];

subplot(5,5,12);
plot(ak_varname_time_02142/3600,ak_varname_Mvals_02142/1e4,'r')
axis([0 200 0 4.6])
ax = gca;
ylabel('MecA (count/10^4)','FontSize',8);
ax.XTick = [];
ax.YTick = [0, 2, 4];

subplot(5,5,13);
plot(ak_varname_time_08450/3600,ak_varname_Mvals_08450/1e4,'r')
axis([0 200 0 4.6])
ax = gca;
ylabel('MecA (count/10^4)','FontSize',8);
ax.XTick = [];
ax.YTick = [0, 2, 4];

subplot(5,5,14);
plot(ak_varname_time_60000/3600,ak_varname_Mvals_60000/1e4,'r')
axis([0 200 0 4.6])
ax = gca;
ylabel('MecA (count/10^4)','FontSize',8);
ax.XTick = [];
ax.YTick = [0, 2, 4];

subplot(5,5,15);
plot(ak_varname_time_150000/3600,ak_varname_Mvals_150000/1e4,'r')
axis([0 200 0 4.6])
ax = gca;
ylabel('MecA (count/10^4)','FontSize',8);
ax.XTick = [];
ax.YTick = [0, 2, 4];


%ComS values vs Time
subplot(5,5,16);
plot(ak_varname_time_00476/3600,ak_varname_Svals_00476/1e4,'k')
axis([0 200 0 2.3])
ax = gca;
xlabel('Time (hours)','FontSize',8);
ylabel('ComS (count/10^4)','FontSize',8);
ax.XTick = [0, 100, 200];
ax.YTick = [0, 1, 2];

subplot(5,5,17);
plot(ak_varname_time_02142/3600,ak_varname_Svals_02142/1e4,'k')
axis([0 200 0 2.3])
ax = gca;
xlabel('Time (hours)','FontSize',8);
ylabel('ComS (count/10^4)','FontSize',8);
ax.XTick = [0, 100, 200];
ax.YTick = [0, 1, 2];

subplot(5,5,18);
plot(ak_varname_time_08450/3600,ak_varname_Svals_08450/1e4,'k')
axis([0 200 0 2.3])
ax = gca;
xlabel('Time (hours)','FontSize',8);
ylabel('ComS (count/10^4)','FontSize',8);
ax.XTick = [0, 100, 200];
ax.YTick = [0, 1, 2];

subplot(5,5,19);
plot(ak_varname_time_60000/3600,ak_varname_Svals_60000/1e4,'k')
axis([0 200 0 2.3])
ax = gca;
xlabel('Time (hours)','FontSize',8);
ylabel('ComS (count/10^4)','FontSize',8);
ax.XTick = [0, 100, 200];
ax.YTick = [0, 1, 2];

subplot(5,5,20);
plot(ak_varname_time_150000/3600,ak_varname_Svals_150000/1e4,'k')
axis([0 200 0 2.3])
ax = gca;
xlabel('Time (hours)','FontSize',8);
ylabel('ComS (count/10^4)','FontSize',8);
ax.XTick = [0, 100, 200];
ax.YTick = [0, 1, 2];

%Periodograms
subplot(5,5,21);
plot(adjust_f_00476/1e3,adjust_power_00476,'b',adjust_f_00476/1e3,yyy_00476,'r');
xlabel('frequency (mHz)','FontSize',9);
ylabel('Power','FontSize',9);
legend('Power spectrum','Gaussian fit','Location','northeast')

subplot(5,5,22);
plot(adjust_f_02142/1e3,adjust_power_02142,'b',adjust_f_02142/1e3,yyy_02142,'r');
xlabel('frequency (mHz)','FontSize',9);
ylabel('Power','FontSize',9);

subplot(5,5,23);
plot(adjust_f_08450/1e3,adjust_power_08450,'b',adjust_f_08450/1e3,yyy_08450,'r');
xlabel('frequency (mHz)','FontSize',9);
ylabel('Power','FontSize',9);

subplot(5,5,24);
plot(adjust_f_60000/1e3,adjust_power_60000,'b',adjust_f_60000/1e3,yyy_60000,'r');
xlabel('frequency (mHz)','FontSize',9);
ylabel('Power','FontSize',9);

subplot(5,5,25);
plot(adjust_f_150000/1e3,adjust_power_150000,'b',adjust_f_150000/1e3,yyy_150000,'r');
xlabel('frequency (mHz)','FontSize',9);
ylabel('Power','FontSize',9);

