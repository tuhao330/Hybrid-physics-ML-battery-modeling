%% Single Particle Model w/ Electrolyte & Temperature
%   Published December 18, 2016 by Professor Scott Moura
%   Energy, Controls, and Applications Lab (eCAL)
%   University of California, Berkeley
%   http://ecal.berkeley.edu/

%   Code based on publications
%   Battery State Estimation for a Single Particle Model with Electrolyte Dynamics 
%   S. J. Moura, F. Bribiesca Argomedo, R. Klein, A. Mirtabatabaei, M. Krstic 
%   IEEE Transactions on Control System Technology, to appear 
%   DOI: 10.1109/TCST.2016.2571663

%   Optimal Charging of Batteries via a Single Particle Model with Electrolyte and Thermal Dynamics 
%   H. Perez, X. Hu, S. J. Moura
%   2016 American Control Conference
%   DOI: 10.1109/ACC.2016.7525538

clear;
clc;
close all;

% global flag

disp('Single Particle Model w/ Electrolyte & Temperature (SPMeT)')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

% DFN = csvread('DFN/cycle0/DFN_1C_cyc0.csv');
% ttotal = length(DFN)-1;
%% Electrochemical Model Parameters
% Load Lithium Cobolt Oxide Params, adopted from DUALFOIL
run param/params_LCO
% n_Li_s3 = readmatrix('n_Li_s_250.csv');
% % delta_sei3 = readmatrix('delta_sei_250.csv');
% n_Li_s3 = 2.5;% n_Li_s3(150);
% delta_sei0 = delta_sei3(250);
delta_sei0 = 0;
%% Input charge/discharge Current Data %%
% % Current | Positive <=> Discharge, Negative <=> Charge

% Calculate C-rate in terms of [A/m^2] using low/high voltage cutoffs
[cn_low,cp_low] = init_cs(p,p.volt_min);
[cn_high,cp_high] = init_cs(p,p.volt_max);
Delta_cn = cn_high-cn_low;
Delta_cp = cp_low-cp_high;
p.OneC = min(p.epsilon_s_n*p.L_n*Delta_cn*p.Faraday/3600, p.epsilon_s_p*p.L_p*Delta_cp*p.Faraday/3600);
% p.OneC = 39.05; % 39.05 (nLis=2.5)
% p.OneC = 33.9923;
% cn_low = 235.4512;
% cn_high = 1.6766e+04;
% p.OneC = 26.5824;
%%%%%%%%%%%%%%% MANUAL INPUT WITH C-RATE %%%%%%%%%%%%%%%%%%%%%%%%%
p.delta_t = 1;
t = 0:p.delta_t:(260); %1C 3425 %10C 283
% cycle 0: 106638; cycle 250:
% 4.3V-3.1V: 1C 3427; 5C 622; 8C 368
% 4.1V-3.1V: 1C 3423; 5C 613; 8C 363
%     Charge: 1/30 106730(n=2.5);
I1 = 7*p.OneC*ones(1,60);
I2 = 6.5*p.OneC*ones(1,200);
I = [0 I1 I2];
%%%%%%%%%%%%%%% DYNAMIC CHARGE/DISCHARGE CYCLES FROM EXPERIMENTS %%%%%%%%%%%%%%%
% load('input-data/UDDS');
% load('data/UDDSx2_batt_ObsData');
% load('data/UDDS_data_Oct_26_2015_Sample_05sec');
% load('data/LA92x2_batt_ObsData');
% I = -current_exp'/p.Area*20;
% t = time_exp';

% Current = importdata('data/C_rate_US06.txt');
% current_exp=Current(:,1);
% I = current_exp/p.Area*p.OneC;
% t = 0:1:length(I)-1;
% t = t';

% % load('data/US06x3_batt_ObsData');
% I = readmatrix('Current/HWFET_cut_sim.csv');
% I = [0; I(1:1096)];
% t = 0:1:length(I)-1;
% t = t';
% p.delta_t = t(2)-t(1);

% Data structure with time,current, initial condition
data.time = t;
data.cur = I;
NT = length(t);

%% Preallocation & Initial Conditions

%%% Finite difference for spherical particle
% p.Nr = 30; % 100 Make this very large so it closely approximates the true model
p.Nr = 100;
Nr = p.Nr;
p.delta_r_n = 1/p.Nr;
p.delta_r_p = 1/p.Nr;
r_vec = (0:p.delta_r_n:1)';
r_vecx = r_vec(2:end-1);

% Finite difference points along x-coordinate
% p.Nxn = 10;
% p.Nxs = 5;
% p.Nxp = 10;
p.Nxn = 70;
p.Nxs = 35;
p.Nxp = 70;
p.Nx = p.Nxn+p.Nxs+p.Nxp;
Nx = p.Nx - 3;
x_vec_spme = linspace(0,1,Nx+4);

p.delta_x_n = 1 / p.Nxn;
p.delta_x_s = 1 / p.Nxs;
p.delta_x_p = 1 / p.Nxp;

% Output Discretization params
disp('Discretization Params:');
fprintf(1,'No. of FDM nodes in Anode | Separator | Cathode : %1.0f | %1.0f | %1.0f\n',p.Nxn,p.Nxs,p.Nxp);
fprintf(1,'No. of FDM nodes in Single Particles : %1.0f\n',p.Nr);
fprintf(1,'Time Step : %2.2f sec\n',p.delta_t);
disp(' ');


%%% INITIAL CONDITIONS
% Solid concentration
V0 = 4.1;
[csn0,csp0] = init_cs(p,V0);
% [csn00,csp00] = Equ_init(p,n_Li_s3,V0); % calculate equivalent IC
% [csn0,csp0] = Equ_init2(p,2.5,V0); % calculate equivalent IC
c_n0 = csn0 * ones(p.Nr-1,1);
c_p0 = csp0 * ones(p.Nr-1,1);

% Electrolyte concentration
ce0 = p.c_e*ones(Nx,1);

% Temperature
T10 = p.T_amb;
T20 = p.T_amb;

% SEI layer
% delta_sei0 = 0;
% delta_sei0 = 2.6660e-06;

disp('Initial Conditions:');
fprintf(1,'Voltage : %1.3f V\n',V0);
fprintf(1,'Normalized Solid Concentration in Anode | Cathode : %1.2f | %1.2f\n',csn0/p.c_s_n_max,csp0/p.c_s_p_max);
fprintf(1,'Initial SoC in Anode : %1.3f\n',(csn0-cn_low)/Delta_cn);
fprintf(1,'Electrolyte Concentration : %2.3f kmol/m^3\n',ce0(1)/1e3);
fprintf(1,'Temperature in Roll | Can : %3.2f K | %3.2f K \n',T10,T20);
fprintf(1,'SEI Layer in Anode : %2f um \n',delta_sei0*1e6);
disp(' ');

%% Generate Constant System Matrices

% Electrolyte concentration matrices
[M1n,M2n,M3n,M4n,M5n, M1s,M2s,M3s,M4s, M1p,M2p,M3p,M4p,M5p, C_ce] = c_e_mats(p);

p.ce.M1n = M1n;
p.ce.M2n = M2n;
p.ce.M3n = M3n;
p.ce.M4n = M4n;
p.ce.M5n = M5n;

p.ce.M1s = M1s;
p.ce.M2s = M2s;
p.ce.M3s = M3s;
p.ce.M4s = M4s;

p.ce.M1p = M1p;
p.ce.M2p = M2p;
p.ce.M3p = M3p;
p.ce.M4p = M4p;
p.ce.M5p = M5p;

p.ce.C = C_ce;

clear M1n M2n M3n M4n M5n M1s M2s M3s M4s M1p M2p M3p M4p M5p C_ce;

%% Simulate SPMeT Plant
tic;
disp('Simulating SPMeT Plant...');

% Initial Conditions
x0 = [c_n0; c_p0; ce0; T10; T20; delta_sei0];

% INTEGRATE !!!!
[t,x] = ode23s(@(t,x) ode_spmet(t,x,data,p),t,x0);

% Parse states
c_s_n = x(:,1:(p.Nr-1));
c_s_p = x(:,p.Nr : 2*(p.Nr-1));
c_ex = x(:,2*p.Nr-1 : 2*p.Nr-1+p.Nx-4);
T1 = x(:,end-2);
T2 = x(:,end-1);
delta_sei = x(:,end);


% Output Function %%%
V = zeros(NT,1);
V_spm = zeros(NT,1);
SOC_n = zeros(NT,1);
SOC_p = zeros(NT,1);
c_ss_n = zeros(NT,1);
c_ss_p = zeros(NT,1);
c_n = zeros(NT,p.Nr+1);
c_p = zeros(NT,p.Nr+1);
c_e = zeros(p.Nx+1,NT);
n_Li_s = zeros(NT,1);
j_s = zeros(NT,1);
% flag = 0;

for k = 1:NT
    
    % Compute outputs
    [~,V(k),V_spm(k),SOC_n(k),SOC_p(k),c_ss_n(k),c_ss_p(k),c_e(:,k),j_s(k)] = ...
        ode_spmet(t(k),x(k,:)',data,p);
%     if flag == 1 
%         break
%     end
    % Aggregate Solid concentrations
    c_n(k,:) = [c_s_n(k,1), c_s_n(k,:), c_ss_n(k)];
    c_p(k,:) = [c_s_p(k,1), c_s_p(k,:), c_ss_p(k)];
    
    % Total Moles of Lithium in Solid
    n_Li_s(k) = (3*p.epsilon_s_p*p.L_p*p.Area) * trapz(r_vec,r_vec.^2.*c_p(k,:)') ...
            + (3*p.epsilon_s_n*p.L_n*p.Area) * trapz(r_vec,r_vec.^2.*c_n(k,:)');
    
end

n_Li_loss = sum(j_s)*p.L_n*p.a_s_n;
% Output Elapsed time
simtime = toc;
fprintf(1,'Elapsed time: %4.1f sec or %2.2f min \n',simtime,simtime/60);

disp('To plots results, run...');
disp(' plot_spmet')
disp(' animate_spmet')

%%
% data = zeros(length(V_spm),6);
% data(:,1) = V_spm;
% data(:,2) = (SOC_n*p.c_s_n_max-cn_low)/Delta_cn; %bulk SoC
% data(:,3) = (c_ss_n-cn_low)/Delta_cn; %surface SoC
% % data(:,4) = (c_n(:,89)-cn_low)/Delta_cn; %surface SoC 1/10
% data(:,4) = T1; %Core temp
% data(:,5) = T2; %Surf temp
% data(:,6) = I;

% csvwrite("ThermalData/SPMT_Test.csv",data)

%% Data generation
% clear all; close all; clc
% run param/params_LCO
% n_Li_s3 = readmatrix('n_Li_s_250.csv');
% n_Li_s3 = n_Li_s3(200);
% data = readmatrix('SPM_LA92_cyc250');
% [csn00,csp00] = Equ_init(p,n_Li_s3,3.9); % calculate equivalent IC
% data(:,6) = ones(length(data),1)*round(csn00/1.6766e+04,3); % SoH
% csvwrite("SPM_LA92_cyc200_2.csv",data)

%%
% x = rescale(c_e(1,:)');
% y = rescale(data(:,2)-data(:,3));
% figure
% plot(rescale(x))
% hold on
% plot(rescale(y))
% 
% [c,lags] = xcorr(x,y,'normalized');
% figure
% stem(lags,c)