%% Plot Single Particle Model w/ Electrolyte (SPMe) Results
%   Published November 5, 2016 by Professor Scott Moura
%   Energy, Controls, and Applications Lab (eCAL)
%   University of California, Berkeley
%   http://ecal.berkeley.edu/

close all;
fs = 16;

%% Plot Current, SOC, Voltage
%%% Static Plant Data
figure(1); clf;
set(gcf,'Position',[99,29,623,669]);

% Current
subplot(311);
plot(t,I/p.OneC,'LineWidth',2);
ylabel('Current [C-Rate]','FontSize',fs);
legend({'$$I(t)$$'},'interpreter','latex','Fontsize',fs)
set(gca,'FontSize',fs)
xlim([0, t(end)])

% Surface Concentrations
subplot(312);
plot(t,c_ss_n/p.c_s_n_max,'b-',t,c_ss_p/p.c_s_p_max,'r-','LineWidth',2);
hold on;
plot(t,SOC_n,'b--',t,SOC_p,'r--','LineWidth',2);
leg_css = {'$$\theta^-(t)$$';'$$\theta^+(t)$$';'$$\overline{\theta}^-(t)$$';'$$\overline{\theta}^+(t)$$'};
legend(leg_css,'FontSize',fs,'Interpreter','latex','Fontsize',fs)
ylabel('Surface & Bulk SOC. [-]','FontSize',fs);
set(gca,'FontSize',fs)
xlim([0, t(end)])

% Voltage
subplot(313);
plot(t,V,'LineWidth',2);
ylabel('Voltage [V]','FontSize',fs);
xlabel('Time [sec]','FontSize',fs);
legend({'$$V(t)$$'},'interpreter','latex','Fontsize',fs)
set(gca,'FontSize',fs)
xlim([0, t(end)])