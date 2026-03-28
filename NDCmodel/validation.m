clear all;close all;clc

global Ntotal N_end Vmin

Vmin = 2.5;
%% Model parameters

Cb = 5777; Cs = 3273.9; Rb = 0.0294; Rs = 0; %2.5V min

% x_real: 1st row: SoC-OCV polynominal; 2nd row: NDC parameters
x_real = [ -3.0338*1e4    9.6609*1e4    0.1148*1e4   -1.1439*1e4    2.6867*1e4 ...
              0.0120    0.0163    0.0163    0.0535   14.6658    0.0122    1.2343    0.0084    0.0524]; %2.5V min
        
R1 = x_real(13); C1 = 1/x_real(14)/R1;

%% current-voltage data
% data = readmatrix("..\Dataset\NDC\(sam)0.08A_discharge.csv"); %1/30C
% data = readmatrix("..\Dataset\NDC\(sam)2.5A_discharge.csv"); %1C
% data = readmatrix("..\Dataset\NDC\(sam)5A_discharge.csv"); %2C
% data = readmatrix("..\Dataset\NDC\(sam)7.5A_discharge.csv"); %3C
% data = readmatrix("..\Dataset\NDC\(sam)10A_discharge_off.csv"); %4C
% data = readmatrix("..\Dataset\NDC\(sam)12.5A_discharge.csv"); %5C
% data = readmatrix("..\Dataset\NDC\(sam)15A_discharge.csv"); %6C
% data = readmatrix("..\Dataset\NDC\(sam)17.5A_discharge.csv"); %7C
data = readmatrix("..\Dataset\NDC\(sam)20A_discharge.csv"); %8C
% data = readmatrix("..\Dataset\NDC\(sam)UDDS2_discharge.csv"); %UDDS_1
% data = readmatrix("..\Dataset\NDC\(sam)US062_discharge.csv"); %US06_1
% data = readmatrix("..\Dataset\NDC\(sam)LA922_discharge.csv"); %LA92_1
% data = readmatrix("..\Dataset\NDC\(sam)SC042_discharge.csv"); %SC04_1
% data = readmatrix("..\Dataset\NDC\(sam)WLTC2_discharge_off.csv"); %WLTC_1

data = data(2:358,:);
voltage = data(:,5);current = data(:,6);  Vh = voltage(1); % column vector
temperature = data(:,8); temperature(end) = temperature(end-1);
time = 0:1:length(current)-1;

%% state-space equations and initial state setting
A = [-1/Cb/(Rb+Rs) 1/Cb/(Rb+Rs); 1/Cs/(Rb+Rs) -1/Cs/(Rb+Rs)];
B = [Rs/Cb/(Rb+Rs); Rb/Cs/(Rb+Rs)];
N = length(current); dt = 1;
[F,G] = c2d(A,B,dt);
xh = zeros(2,N+1);
Vsh = 1; % SoC(1) = 1;
xh(:,1) = [Vsh; Vsh];
Vt = zeros(1,N);
Vs = zeros(1,N);
Vb = zeros(1,N);
fVs = zeros(1,N);
V0 = zeros(1,N);
%%
A1 = -1/R1/C1; B1 = -1/C1;
[F1,G1] = c2d(A1,B1,dt);
V1 = zeros(1,N+1);
%% compute terminal voltage Vt and SoC
for i = 1:N
    Vs(i) = [0 1]*xh(:,i);
    Vb(i) = [1 0]*xh(:,i);
    SoC(i) = (Vb(i)*Cb+Vs(i)*Cs)/(Cb+Cs);
    Vt(i) = current(i)*(x_real(8)+x_real(9)*exp(-x_real(10)*SoC(i))+x_real(11)*exp(-x_real(12)*(1-SoC(i))))...
            -V1(i)+...
            (x_real(1)*Vs(i)^2 + x_real(2)*Vs(i) + Vmin*x_real(3))/(Vs(i)^3 + x_real(4)*Vs(i)^2 + x_real(5)*Vs(i) + x_real(3));
    R0(i) = x_real(8)+x_real(9)*exp(-x_real(10)*SoC(i))+x_real(11)*exp(-x_real(12)*(1-SoC(i)));
    V0(i) = -current(i)*R0(i);
    fVs(i) = (x_real(1)*Vs(i)^2 + x_real(2)*Vs(i) + Vmin*x_real(3))/(Vs(i)^3 + x_real(4)*Vs(i)^2 + x_real(5)*Vs(i) + x_real(3));
    xh(:,i+1) = F*xh(:,i)+G*current(i);
    V1(1,i+1) = F1*V1(1,i)+G1*current(i);
end

%% figure
figure
plot(time,voltage,'k-','LineWidth',1.3);
hold on;
plot(time,Vt,'b:','LineWidth',1.5);
set(gca,'Units','normalized','Position',[.15 .2 .75 .6],...
'FontUnits','points','FontWeight','normal','FontSize',22,'FontName','Times');
xlabel('Time (s)');ylabel('Voltage (V)');
legend('True voltage','NDC voltage');

%% Output data for hybrid modeling
Output = [Vt' Vb' Vs' V1(1:end-1)' SoC' R0' current temperature];
writematrix(Output,'NDC_8C_sam_etd.csv')
writematrix(voltage,'Real_8C_sam_etd.csv')