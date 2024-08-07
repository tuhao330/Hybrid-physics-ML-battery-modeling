%% Single Particle Model w/ Electrolyte
%   Published November 5, 2016 by Professor Scott Moura
%   Energy, Controls, and Applications Lab (eCAL)
%   University of California, Berkeley
%   http://ecal.berkeley.edu/

%   Code based on publication
%   Battery State Estimation for a Single Particle Model with Electrolyte Dynamics 
%   S. J. Moura, F. Bribiesca Argomedo, R. Klein, A. Mirtabatabaei, M. Krstic 
%   IEEE Transactions on Control System Technology, to appear 
%   DOI: 10.1109/TCST.2016.2571663


clear;
clc;
close all;

disp('Single Particle Model w/ Electrolyte (SPMe)')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

%% Electrochemical Model Parameters
% Load Lithium Cobolt Oxide Params, adopted from DUALFOIL
run param/params_LCO

%% Input charge/discharge Current Data %%
% % Current | Positive <=> Discharge, Negative <=> Charge

% Calculate C-rate in terms of [A/m^2] using low/high voltage cutoffs
[cn_low,cp_low] = init_cs(p,p.volt_min);
[cn_high,cp_high] = init_cs(p,p.volt_max);
Delta_cn = cn_high-cn_low;
Delta_cp = cp_low-cp_high;
p.OneC = min(p.epsilon_s_n*p.L_n*Delta_cn*p.Faraday/3600, p.epsilon_s_p*p.L_p*Delta_cp*p.Faraday/3600);

%%%%%%%%%%%%%%% MANUAL INPUT WITH C-RATE %%%%%%%%%%%%%%%%%%%%%%%%%
p.delta_t = 1;
t = 0:p.delta_t:(120);
I = 2*p.OneC*ones(size(t));


%%%%%%%%%%%%%%% DYNAMIC CHARGE/DISCHARGE CYCLES FROM EXPERIMENTS %%%%%%%%%%%%%%%
% load('input-data/UDDS');
% 
% I = -current_exp'/p.Area*10;
% t = time_exp';
% p.delta_t = t(2)-t(1);


% Data structure with time,current, initial condition
data.time = t;
data.cur = I;
NT = length(t);

%% Preallocation & Initial Conditions

%%% Finite difference for spherical particle
p.Nr = 30; % 100 Make this very large so it closely approximates the true model
Nr = p.Nr;
p.delta_r_n = 1/p.Nr;
p.delta_r_p = 1/p.Nr;
r_vec = (0:p.delta_r_n:1)';
r_vecx = r_vec(2:end-1);

% Finite difference points along x-coordinate
p.Nxn = 10;
p.Nxs = 5;
p.Nxp = 10;
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
V0 = 3.8;
[csn0,csp0] = init_cs(p,V0);
c_n0 = csn0 * ones(p.Nr-1,1);
c_p0 = csp0 * ones(p.Nr-1,1);

% Electrolyte concentration
ce0 = p.c_e;
data.c_e0 = ce0*ones(Nx,1);

% Temperature
T0 = p.T_amb;

disp('Initial Conditions:');
fprintf(1,'Voltage : %1.3f V\n',V0);
fprintf(1,'Normalized Solid Concentration in Anode | Cathode : %1.2f | %1.2f\n',csn0/p.c_s_n_max,csp0/p.c_s_p_max);
fprintf(1,'Electrolyte Concentration : %2.3f kmol/m^3\n',ce0/1e3);
fprintf(1,'Temperature : %3.2f Kelvin\n',T0);
disp(' ');

%% ISOTHERMAL assumption: Set temperature-dependent parameters as constants

% Solid phase diffusivity
p.D_s_n = p.D_s_n0 * exp(p.E.Dsn/p.R*(1/p.T_ref - 1/T0));
p.D_s_p = p.D_s_n0 * exp(p.E.Dsp/p.R*(1/p.T_ref - 1/T0));

% Kinetic reaction rate
p.k_n = p.k_n0 * exp(p.E.kn/p.R*(1/p.T_ref - 1/T0));
p.k_p = p.k_p0 * exp(p.E.kp/p.R*(1/p.T_ref - 1/T0));

%% Generate SPMe System Matrices (A,C,B,D)
% Construct (A,B) matrices for solid-phase Li diffusion
[A_n,A_p,B_n,B_p,C_n,C_p,D_n,D_p,MN_mats] = spm_plant_obs_mats(p);

sys_n = ss(A_n,B_n,C_n,D_n);
sys_p = ss(A_p,B_p,C_p,D_p);

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

%% Simulate SPMe Plant
tic;
disp('Simulating SPMe Plant...');

%%% Solid-phase Diffusion %%%
% Spherical Diffusion Simulation
[c_ss_n,t,c_nx] = lsim(sys_n,I,t,c_n0);     % Anode Particle
[c_ss_p,t,c_px] = lsim(sys_p,I,t,c_p0);     % Cathode Particle

% Aggregate states
c_n = [c_nx(:,1), c_nx, c_ss_n];
c_p = [c_px(:,1), c_px, c_ss_p];


%%% Electrolyte-phase Diffusion %%%
% Simulate electrolyte dynamics for a given input & initial conditions
c_e = electrolyte_scott(data,p);

% Electrolyte concentration endpoints
ce0n = c_e(1,:);
cens = c_e(p.Nxn+1,:);
cesp = c_e(p.Nxn+p.Nxs+1,:);
ce0p = c_e(end,:);

% Average electrolyte concentrations
cen_bar = mean(c_e(1:p.Nxn+1,:));
ces_bar = mean(c_e((p.Nxn+1):(p.Nxn+p.Nxs+1),:));
cep_bar = mean(c_e((p.Nxn+p.Nxs+1):(p.Nxn+p.Nxs+p.Nxp+1),:));

%%% Voltage Output Function %%%
V_noVCE = zeros(size(t));
V_electrolyteCond = zeros(size(t));
V_electrolytePolar = zeros(size(t));
V = zeros(size(t));
V_spm = zeros(size(t));
SOC = zeros(size(t));
n_Li_s = zeros(size(t));

for k = 1:NT
    
    % SPMe Voltage w/o electrolyte concentration term
    V_noVCE(k) = nonlinearSPMOutputVoltage_Scott(p,c_ss_n(k),c_ss_p(k),cen_bar(k),ces_bar(k),cep_bar(k),I(k));
    
    % Overpotentials due to electrolyte subsystem
    kap_n = electrolyteCond(cen_bar(k));
    kap_s = electrolyteCond(ces_bar(k));
    kap_p = electrolyteCond(cep_bar(k));
    
    kap_n_eff = kap_n * p.epsilon_e_n.^(p.brug);
    kap_s_eff = kap_s * p.epsilon_e_s.^(p.brug);
    kap_p_eff = kap_p * p.epsilon_e_p.^(p.brug);
    
    % Activity coefficient
    dfca_n = electrolyteAct(cen_bar(k),T0,p);
    dfca_s = electrolyteAct(ces_bar(k),T0,p);
    dfca_p = electrolyteAct(cep_bar(k),T0,p);
    
    if( (ce0p(k) < 0) || ce0n(k) < 0)
        error('Error: The electrolyte concentration became negative, which is non-physical. You probably tried a C-rate that is too high. Try reducing your C-rate. (This occured because the SPMe does not capture a self-limiting process that prevents c_e from becoming negative.)');
    end
    
    % Overpotential due to electrolyte conductivity
    V_electrolyteCond(k) = (p.L_n/(2*kap_n_eff) + 2*p.L_s/(2*kap_s_eff) + p.L_p/(2*kap_p_eff))*I(k); ...
        
    % Overpotential due to electrolyte polarization
    V_electrolytePolar(k) = (2*p.R*p.T_amb)/(p.Faraday) * (1-p.t_plus)* ...
            ( (1+dfca_n) * (log(cens(k)) - log(ce0n(k))) ...
             +(1+dfca_s) * (log(cesp(k)) - log(cens(k))) ...
             +(1+dfca_p) * (log(ce0p(k)) - log(cesp(k))));
    
    % Add 'em up!
    V(k) = V_noVCE(k) + V_electrolyteCond(k) + V_electrolytePolar(k);
    
    % SPM Voltage
    V_spm(k) = nonlinearSPMOutputVoltage_Scott(p,c_ss_n(k),c_ss_p(k),p.c_e,p.c_e,p.c_e,I(k));
    
    % State-of-Charge (Bulk)
    SOC_n(k) = 3/p.c_s_n_max * trapz(r_vec,r_vec.^2.*c_n(k,:)');
    SOC_p(k) = 3/p.c_s_p_max * trapz(r_vec,r_vec.^2.*c_p(k,:)');
    
    % Total Moles of Lithium in Solid
    n_Li_s(k) = (3*p.epsilon_s_p*p.L_p*p.Area) * trapz(r_vec,r_vec.^2.*c_p(k,:)') ...
            + (3*p.epsilon_s_n*p.L_n*p.Area) * trapz(r_vec,r_vec.^2.*c_n(k,:)');
    
end

% Output Elapsed time
simtime = toc;
fprintf(1,'Elapsed time: %4.1f sec or %2.2f min \n',simtime,simtime/60);

disp('To plots results, run...');
disp(' plot_spme')
disp(' animate_spme')


