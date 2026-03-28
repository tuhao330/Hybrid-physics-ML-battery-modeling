clear all;close all;clc;

global Ic Vh Vsh t_end Ntotal N_end Cb Cs Rb  Qt Vmin

data = readmatrix("..\Dataset\NDC\(sam)2.5A_discharge.csv");

Vmin = 2.5;

data = data(2:end,:);%
Qt = 9.0509e+03; % Vmin=2.5V
voltage = data(:,5);current = data(:,6); Vh = voltage(1); % column vector
N_start = 1;
N_end = length(voltage);
%% problem setting
I = current(N_start:N_end)'; Ic = -2.5;% the constant current, negative
V = voltage(N_start:N_end)';  Vsh = 1; % bar{Vs}
t = N_start:1:N_end;
t_end = t(end); 
uu = [t;V]; % input in fmincon function

x0 = .01*ones(1,9);
lb = [0.005    0.005    0.01     0.05       1       0.01    1          0.001        1/800]; % lower bound
ub = [0.2       0.2     0.09      0.35       15      0.12    15        0.03         1/10]; % uppper bound

poly = [-3.0338    9.6609    0.1148   -1.1439    2.6867]*1e4; % SoC-OCV polynominal
%% parameter identification
fun_error = @(x)((poly(1)*(Vsh+1/Qt*Ic*uu(1,:)+Ic*x(1)*(1-exp(-uu(1,:)*x(2)))).^2 +...
    poly(2)*(Vsh+1/Qt*Ic*uu(1,:)+Ic*x(1)*(1-exp(-uu(1,:)*x(2)))) +...
    Vmin*poly(3))./...
    ((Vsh+1/Qt*Ic*uu(1,:)+Ic*x(1)*(1-exp(-uu(1,:)*x(2)))).^3 +...
    poly(4)*(Vsh+1/Qt*Ic*uu(1,:)+Ic*x(1)*(1-exp(-uu(1,:)*x(2)))).^2 +...
    poly(5)*(Vsh+1/Qt*Ic*uu(1,:)+Ic*x(1)*(1-exp(-uu(1,:)*x(2)))) +...
    poly(3))+...
    Ic*x(3)+Ic*x(4)*exp(-x(5)*(1+1/Qt*Ic*uu(1,:)/Vsh))+Ic*x(6)*exp(-x(7)*(-1/Qt*Ic*uu(1,:)/Vsh))+...
    Ic*x(8)*(1-exp(-uu(1,:)*x(9)))-uu(2,:))*...
    ((poly(1)*(Vsh+1/Qt*Ic*uu(1,:)+Ic*x(1)*(1-exp(-uu(1,:)*x(2)))).^2 +...
    poly(2)*(Vsh+1/Qt*Ic*uu(1,:)+Ic*x(1)*(1-exp(-uu(1,:)*x(2)))) +...
    Vmin*poly(3))./...
    ((Vsh+1/Qt*Ic*uu(1,:)+Ic*x(1)*(1-exp(-uu(1,:)*x(2)))).^3 +...
    poly(4)*(Vsh+1/Qt*Ic*uu(1,:)+Ic*x(1)*(1-exp(-uu(1,:)*x(2)))).^2 +...
    poly(5)*(Vsh+1/Qt*Ic*uu(1,:)+Ic*x(1)*(1-exp(-uu(1,:)*x(2)))) +...
    poly(3))+...
    Ic*x(3)+Ic*x(4)*exp(-x(5)*(1+1/Qt*Ic*uu(1,:)/Vsh))+Ic*x(6)*exp(-x(7)*(-1/Qt*Ic*uu(1,:)/Vsh))+...
    Ic*x(8)*(1-exp(-uu(1,:)*x(9)))-uu(2,:))';

x_guess = fmincon(fun_error,x0,[],[],[],[],lb,ub);
%% print parameter identification results
Cb = fun_RC_parameter(x_guess)*[1;0;0];
Cs = fun_RC_parameter(x_guess)*[0;1;0];
Rb = fun_RC_parameter(x_guess)*[0;0;1];
SoC_end = fun_SoC(x_guess,t); SoC_end(end); % SoC at the end of discharge
R1 = x_guess(8);
R1C1 = 1/x_guess(9);
sprintf('Q=%0.5s,Cb=%0.5s,Cs=%0.5s,Rb=%0.5s,R1=%0.5s,R1C1=%0.5s',...
        num2str(Cb/3600+Cs/3600),num2str(Cb),num2str(Cs),...
        num2str(Rb),num2str(R1),num2str(R1C1))
%% fitting result
figure;
plot([0 t],[Vh V],'b-','LineWidth',1.5); hold on
plot([0 t],[Vh fun_terminal_voltage(x_guess,t)],'r:','LineWidth',3);  
set(gca,'Units','normalized','Position',[.15 .2 .75 .7],...
'FontUnits','points','FontWeight','normal','FontSize',18,'FontName','Times')
axis([0,t_end,2.4,4.3]); xlabel('Time (s)'); ylabel('Voltage (V)')
legend('Measured terminal voltage','Predicted terminal voltage');
%% plot: SoC vs R0
funR0 = @(x,udata)x(3)+x(4)*exp(-x(5)*udata)+x(6)*exp(-x(7)*(1-udata));
figure;plot(0:0.01:1,funR0(x_guess,0:0.01:1),'b-','LineWidth',1.5);
set(gca,'Units','normalized','Position',[.15 .2 .75 .7],...
'FontUnits','points','FontWeight','normal','FontSize',18,'FontName','Times')
xlabel('SoC');ylabel('$R_0$ ($\Omega$)','interpreter','latex')
%% plot: time vs SoC
figure;plot(t,fun_SoC(x_guess,t),'b-'); set(gca,'FontSize',18,'FontName','Times');
xlabel('Time (s)');ylabel('SoC')
%% plot: Vb and Vs
figure;plot(t,fun_Vb(x_guess,t),'b-'); hold on; plot(t,fun_Vs(x_guess,t),'r-'); set(gca,'FontSize',18,'FontName','Times');
legend('Vb','Vs')