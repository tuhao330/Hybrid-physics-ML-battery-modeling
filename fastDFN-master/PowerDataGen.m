clear all; close all; clc
global n_crate
n_crate = 8;
dfndata = readmatrix('ThermalData\DFN_10C.csv');
spmtdata = readmatrix('ThermalData\SPMT_10C.csv');
lng = length(spmtdata);
crate = 33.9923;

engy = zeros(lng,1);

for startpt = 1:lng
    
    I = [ spmtdata(1:startpt-1,6); n_crate*crate*ones(10000,1)];
    engy(startpt) = dfn_scott_func(I,startpt);
    writematrix([spmtdata engy*crate/3600],'ThermalData\SPMT_10C_PW8C.csv')
    
end

%% save data