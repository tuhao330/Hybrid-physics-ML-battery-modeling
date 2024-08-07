% this code is to identify SoC-OCV from a 0.1 A discharging current,  
% i.e., cell20_test983_discharge_0.1A
close all;clear all; clc
global SoC V
% data = load('C:\Users\n259t568.HOME\Dropbox\Dataset\txtdata\cell20\cell20_test983_discharge_0.1A.txt');
% data = readmatrix("cell 30\0.1A_discharge.csv");
data = readmatrix("Samsung cell\(sam)0.08A_discharge.csv");

global Vmin
Vmin = 2.5;
for ct = 1:1e6
    if data(ct,5)<Vmin
        break
    end
end

data = data(2:ct-1,:);
% time = data(:,1);voltage = data(:,2);current = data(:,3); Vh = voltage(1); % column vector
voltage = data(:,5);current = data(:,6); Vh = voltage(1);

% Ntotal = length(time); % total points in time/voltage/current
% % N_start = 3; 
% for i = 4:Ntotal
%     if current(i) < current(5)+0.05 && current(i+1) > current(5)+0.05
%         N_end = i;
%     end
% end
N_end = length(voltage);

V = voltage(1:N_end)';
% fun = @(x,udata)2.0+x(1)*udata+x(2)*udata.^2+x(3)*udata.^3+x(4)*udata.^4+x(5)*udata.^5+x(6)*udata.^6+...
%     x(7)*udata.^7+x(8)*udata.^8+x(9)*udata.^9+x(10)*udata.^10+x(11)*udata.^11+x(12)*udata.^12+x(13)*udata.^13+...
%     x(14)*udata.^14+x(15)*udata.^15+x(16)*udata.^16+x(17)*udata.^17+x(18)*udata.^18+x(19)*udata.^19+...
%     x(20)*udata.^20+x(21)*udata.^21+x(22)*udata.^22+x(23)*udata.^23+x(24)*udata.^24+x(25)*udata.^25+...
%                 (Vh-x(1)-x(2)-x(3)-x(4)-x(5)-x(6)-x(7)-x(8)-x(9)-x(10)-x(11)-x(12)-x(13)-x(14)-x(15)-...
%                 x(16)-x(17)-x(18)-x(19)-x(20)-x(21)-x(22)-x(23)-x(24)-x(25)-2.0)*udata.^26;
            
% fun = @(x,udata)(x(1)*udata.^2 + x(2)*udata + 2*x(3))./(udata.^3 + x(4)*udata.^2 + x(5)*udata + x(3));

capacity = -sum(current);
SoC(1) = 1;
for i = 1:N_end-1
    SoC(i+1) = SoC(i)+current(i)/capacity;
end

x0 = 10*ones(1,5);
x_guess = fmincon(@lsq,x0,[],[],[1 1 Vmin-Vh -Vh -Vh],[Vh],[],[]);
% x_guess = lsqcurvefit(fun,x0,SoC,V)
fun = @(x,udata)(x_guess(1)*udata.^2 + x_guess(2)*udata + Vmin*x_guess(3))./...
    (udata.^3 + x_guess(4)*udata.^2 + x_guess(5)*udata + x_guess(3));
figure;plot(SoC,V,'b-','linewidth',1.5);
hold on;plot(SoC,fun(x_guess,SoC),'r:','linewidth',3);
set(gca,'Units','normalized','Position',[.15 .2 .75 .7],...
'FontUnits','points','FontWeight','normal','FontSize',18,'FontName','Times')
axis([0,1,2.4,4.3]); xlabel('SoC'); ylabel('OCV')
legend('Measured OCV','Predicted OCV','Location','SouthEast');

