close all;clear all; clc
global SoC V

data = readmatrix("..\Dataset\NDC\(sam)0.08A_discharge.csv");

global Vmin
Vmin = 2.5;
for ct = 1:1e6
    if data(ct,5)<Vmin
        break
    end
end

data = data(2:ct-1,:);
voltage = data(:,5);current = data(:,6); Vh = voltage(1);

N_end = length(voltage);

V = voltage(1:N_end)';

capacity = -sum(current);
SoC(1) = 1;
for i = 1:N_end-1
    SoC(i+1) = SoC(i)+current(i)/capacity;
end

x0 = 10*ones(1,5);
x_guess = fmincon(@lsq,x0,[],[],[1 1 Vmin-Vh -Vh -Vh],[Vh],[],[]);

fun = @(x,udata)(x_guess(1)*udata.^2 + x_guess(2)*udata + Vmin*x_guess(3))./...
    (udata.^3 + x_guess(4)*udata.^2 + x_guess(5)*udata + x_guess(3));
figure;plot(SoC,V,'b-','linewidth',1.5);
hold on;plot(SoC,fun(x_guess,SoC),'r:','linewidth',3);
set(gca,'Units','normalized','Position',[.15 .2 .75 .7],...
'FontUnits','points','FontWeight','normal','FontSize',18,'FontName','Times')
axis([0,1,2.4,4.3]); xlabel('SoC'); ylabel('OCV')
legend('Measured OCV','Predicted OCV','Location','SouthEast');

