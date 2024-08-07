clear all;close all;clc

global Ntotal N_end Vmin

Vmin = 2.5;
%% nominal parameters, 3A under constant temperature, test1005
% Cb = 10037; Cs = 973; Rb = 0.019; Rs = 0;
% Cb = 10500; Cs = 927.7; Rb = 0.005; Rs = 0;
% Cb = 1684; Cs = 7372; Rb = 0.144; Rs = 0;
Cb = 5777; Cs = 3273.9; Rb = 0.0294; Rs = 0; %2.5V min
% Cb = 8963; Cs = 87.08; Rb = 0.017; Rs = 0; %2.5V min
x_real = [ -3.0338*1e4    9.6609*1e4    0.1148*1e4   -1.1439*1e4    2.6867*1e4 ...
              0.0120    0.0163    0.0163    0.0535   14.6658    0.0122    1.2343    0.0084    0.0524]; %2.5V min
% x_real = [ -3.0338*1e4    9.6609*1e4    0.1148*1e4   -1.1439*1e4    2.6867*1e4 ...
%               0.0171    0.6665   -1.7683   -0.6136    0.9164    2.0283   -0.1580    0.4708    0.0000]; %2.5V min         
R1 = x_real(13); C1 = 1/x_real(14)/R1;
%% variant current data
% clear all; clc
% % load('data/UDDSx2_batt_ObsData');
% % load('data/US06x3_batt_ObsData');
% % current = -Load.I*1.2;
% % I(1768:1774) = I(1768:1774)/7;
% % load('data/LA92x2_batt_ObsData');
% load('data/SC04x4_batt_ObsData');
% % I = load('data/(sam)LA92_2.csv');
% % I = [I;I];
% % plot(-I*1.02/2.5)
% I = -I*1.205;
% 
% % I = [I(1:2731);I(1:1381);I(2732:3482)];
% % current = readmatrix('data/WLTC_Ning.txt');
% current = I;%1.12
% % current = current*3.7;
% for k = 1:length(I)
%     if current(k)>4
%         current(k) = 4;
%     end
%     if current(k)<-2.5*7.8
%         current(k) = -2.5*7.8;
%     end
% end
% % current = current - min(current);
% total = sum(current)/3600 % total discharge capacity
% min = min(current)/2.5 % max dischrage C-rate
% max = max(current)/2.5 % max charge C-rate
% current = round(current,3);
% % time = 0:1:length(current)-1;
% figure
% plot(current/2.5)
% writematrix(current,'(sam)SC04_life_8C.csv')
%% current-voltage data
% data = readmatrix("Samsung cell\(sam)0.08A_discharge.csv"); %1/30C %end
% data = readmatrix("Samsung cell\(sam)2.5A_discharge.csv"); %1C %end
% data = readmatrix("Samsung cell\(sam)5A_discharge.csv"); %2C %1682
% data = readmatrix("Samsung cell\(sam)7.5A_discharge.csv"); %3C %1080 %1134
% data = readmatrix("Samsung cell\(sam)10A_discharge_off.csv"); %4C
% %779 %823
% data = readmatrix("Samsung cell\(sam)10A_discharge_off.csv"); %4C %779
% data = readmatrix("Samsung cell\(sam)12.5A_discharge.csv"); %5C %598 %637
% data = readmatrix("Samsung cell\(sam)15A_discharge.csv"); %6C %477
% data = readmatrix("Samsung cell\(sam)17.5A_discharge.csv"); %7C %389 %424
% data = readmatrix("Samsung cell\(sam)20A_discharge.csv"); %8C %322 %358
data = readmatrix("Samsung cell\(sam)UDDS2_discharge.csv"); %UDDS_1
% data = readmatrix("Samsung cell\(sam)US062_discharge.csv"); %US06_1
% data = readmatrix("Samsung cell\(sam)LA922_discharge.csv"); %LA92_1
% data = readmatrix("Samsung cell\(sam)SC042_discharge.csv"); %SC04_1
% data = readmatrix("Samsung cell\(sam)WLTC2_discharge_off.csv"); %WLTC_1

% for ct = 1:length(data(:,5))
%     if data(ct,5)<3.1
%         break
%     end
% end

data = data(2:358,:);
voltage = data(:,5);current = data(:,6);  Vh = voltage(1); % column vector
temperature = data(:,8); temperature(end) = temperature(end-1);
time = 0:1:length(current)-1;
% Init_SoC = 0.5;
% Ntotal = length(time);
% % N_end = fun_breaking_point(current);
%% state-space equations and initial state setting
A = [-1/Cb/(Rb+Rs) 1/Cb/(Rb+Rs); 1/Cs/(Rb+Rs) -1/Cs/(Rb+Rs)];
B = [Rs/Cb/(Rb+Rs); Rb/Cs/(Rb+Rs)];
N = length(current); dt = 1;
[F,G] = c2d(A,B,dt);
xh = zeros(2,N+1);
Vsh = 1;
xh(:,1) = [Vsh; Vsh];
Vt = zeros(1,N);
Vs = zeros(1,N);
Vb = zeros(1,N);
fVs = zeros(1,N);
% SoC(1) = 1;
V0 = zeros(1,N);
%%
A1 = -1/R1/C1; B1 = -1/C1;
[F1,G1] = c2d(A1,B1,dt);
V1 = zeros(1,N+1);
%% compute terminal voltage Vt and SoC
% SoC = 0;
for i = 1:N
    Vs(i) = [0 1]*xh(:,i);
    Vb(i) = [1 0]*xh(:,i);
    SoC(i) = (Vb(i)*Cb+Vs(i)*Cs)/(Cb+Cs);
    Vt(i) = current(i)*(x_real(8)+x_real(9)*exp(-x_real(10)*SoC(i))+x_real(11)*exp(-x_real(12)*(1-SoC(i))))...
            -V1(i)+...
            (x_real(1)*Vs(i)^2 + x_real(2)*Vs(i) + Vmin*x_real(3))/(Vs(i)^3 + x_real(4)*Vs(i)^2 + x_real(5)*Vs(i) + x_real(3));
%             [Vmin x_real(1:5)]*[1;Vs(i);Vs(i)^2;Vs(i)^3;Vs(i)^4;Vs(i)^5];
    R0(i) = x_real(8)+x_real(9)*exp(-x_real(10)*SoC(i))+x_real(11)*exp(-x_real(12)*(1-SoC(i)));
    V0(i) = -current(i)*R0(i);
%     fVs(i) = [Vmin x_real(1:5)]*[1;Vs(i);Vs(i)^2;Vs(i)^3;Vs(i)^4;Vs(i)^5];
    fVs(i) = (x_real(1)*Vs(i)^2 + x_real(2)*Vs(i) + Vmin*x_real(3))/(Vs(i)^3 + x_real(4)*Vs(i)^2 + x_real(5)*Vs(i) + x_real(3));
    xh(:,i+1) = F*xh(:,i)+G*current(i);
    V1(1,i+1) = F1*V1(1,i)+G1*current(i);
end
%% result presentation
%% figure1 for test1001
% M_end = 9677; %test1001
% figure;plot(time,voltage,'b-','LineWidth',1.5);
% hold on;plot(time,Vt,'r:','LineWidth',3);
% set(gca,'Units','normalized','Position',[.15 .2 .75 .6],...
% 'FontUnits','points','FontWeight','normal','FontSize',22,'FontName','Times');
% xlabel('Time (s)');ylabel('Voltage (V)');
% % this is to generate the square in the figure
% m1=3100;m2=3300;n1=3.5;n2=4; % four points for the square
% m = [m1, m2, m2, m1, m1];n = [n1, n1, n2, n2, n1];
% hold on; plot(m, n, 'k-', 'LineWidth', 0.5); 
% 
% % this is to generate the arrow in the figure
% x4 = [0.6 0.58];y4 = [0.75 0.75];
% annotation('textarrow',x4,y4,'String','','FontSize',16,'FontName','Times',...
%              'FontUnits','points','FontWeight','normal')
%          
% legend('Measured terminal voltage','Predicted terminal voltage');axis([0 M_end 3.1 4.2])
% 
% axes('position',[0.25 0.25 0.1 0.1]);
% box on; your_index = 3100<time & time<3300; % [yy_axis,LLine1,LLine2] = plotyy(x,y1,x,y2);
% plot(time(your_index),voltage(your_index),'b-','LineWidth',1.5);
% hold on;plot(time(your_index),Vt(your_index),'r:','LineWidth',3);
% set(gca,'FontSize',18,'FontName','Times') 

%% figure1 for test1003 
% M_end = 4000; %test1003
% % DoD = 1 - SoC;
% % M_end = DoD(end);
% % figure
% % plot(DoD,voltage,'k-','LineWidth',1.3);
% % hold on;plot(DoD,Vt,'b:','LineWidth',1.5);
% % set(gca,'Units','normalized','Position',[.15 .2 .75 .6],...
% % 'FontUnits','points','FontWeight','normal','FontSize',22,'FontName','Times');
% % xlabel('DoD (%)');ylabel('Voltage (V)');
% % legend('True voltage','NDC voltage');
% % this is to generate the square in the figure
% m1=1200;m2=1300;n1=3.5;n2=4; % four points for the square
% m = [m1, m2, m2, m1, m1];n = [n1, n1, n2, n2, n1];
% hold on; plot(m, n, 'k-', 'LineWidth', 0.5); 
% 
% % this is to generate the arrow in the figure
% x4 = [0.6 0.58];y4 = [0.75 0.75];
% annotation('textarrow',x4,y4,'String','','FontSize',16,'FontName','Times',...
%              'FontUnits','points','FontWeight','normal')
%% figure1 for test1003 
% M_end = 4000; %test1003
M_end = 1 - SoC(end);
figure
% plot(1 - SoC,data2(2:779,5),'r-','LineWidth',1.3);
% hold on
plot(1 - SoC,voltage,'k-','LineWidth',1.3);
hold on;
plot(1 - SoC,Vt,'b:','LineWidth',1.5);
set(gca,'Units','normalized','Position',[.15 .2 .75 .6],...
'FontUnits','points','FontWeight','normal','FontSize',22,'FontName','Times');
xlabel('Time (s)');ylabel('Voltage (V)');
% legend('True voltage','NDC voltage');
%%
% Vhb = readmatrix("predict_Exp_UDDS2_Train_sam_off.csv");
% hold on; plot(time,Vhb,'r-.','LineWidth',1.3)
% legend('True voltage','NDC voltage','Hybrid voltage');axis([0 M_end 3 4.5])
%%
% legend('Measured terminal voltage','Predicted terminal voltage');axis([0 M_end 2.5 4.3])

% axes('position',[0.25 0.25 0.1 0.1]);
% box on; your_index = 1200<time & time<1300; % [yy_axis,LLine1,LLine2] = plotyy(x,y1,x,y2);
% plot(time(your_index),voltage(your_index),'b-','LineWidth',1.5);
% hold on;plot(time(your_index),Vt(your_index),'r:','LineWidth',3);
% set(gca,'FontSize',18,'FontName','Times') 
% %% figure2 
% figure;plot(time,current,'b-','LineWidth',1.5);
% set(gca,'Units','normalized','Position',[.15 .2 .75 .6],...
% 'FontUnits','points','FontWeight','normal','FontSize',22,'FontName','Times');
% xlabel('Time (s)');ylabel('Current (A)');axis([0 M_end -3 0.25])
% 
% fitdata(1,:) = time(1:N_end-1)';
% fitdata(2,:) = voltage(1:N_end-1)';
% fitdata(3,:) = Vt(1:N_end-1)';

%% Output data for hybrid modeling
% Output = [Vt' Vb' Vs' V1(1:end-1)' SoC' R0' current temperature];
% writematrix(Output,'NDC_8C_sam_etd.csv')
% writematrix(voltage,'Real_8C_sam_etd.csv')