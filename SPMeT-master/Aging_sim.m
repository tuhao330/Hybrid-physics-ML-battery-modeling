clear all; close all; clc

global flag

cycle = 1;
n_Li_s = zeros(cycle,1);
delta_sei = zeros(cycle,1);
idxx = zeros(cycle,1);

n_Li_s0 = 2.5;
delta_sei0 = 0;

% n_Li_s0 = readmatrix('n_Li_s_250.csv');
% n_Li_s0 = n_Li_s0(1);
% delta_sei0 = readmatrix('delta_sei_250.csv');
% delta_sei0 = delta_sei0(1);

for i = 1:cycle
disp(" ")
disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
disp("Charging Cycle: " + num2str(i))
flag = 0;
otpt = spmetfunc(n_Li_s0, delta_sei0);
% otpt(1) = round(otpt(1),4);
% otpt(2) = round(otpt(2)*1e6,4)/1e6;
n_Li_s(i) = otpt(1);
delta_sei(i) = otpt(2);
n_Li_s0 = otpt(1);
delta_sei0 = otpt(2);
% idxx(i) = otpt(3:end)';
end
Li_loss = [0 ; otpt(3:end)'];%Lithium loss every second
Total_loss = zeros(length(Li_loss),1);
for ct = 1:length(Li_loss)
    Total_loss(ct) = sum(Li_loss(1:ct));
end
time = 0:1:length(Li_loss)-1;
plot(time,-Total_loss)
xlabel('Charging Time(s)')
ylabel('Li-ion loss(mol)')
writematrix(-Total_loss,'loss_1_20.csv')
% writematrix(n_Li_s,'n_Li_s_300.csv')
% writematrix(delta_sei,'delta_sei_300.csv')
% %
% n_Li_s1 = readmatrix('n_Li_s_73.csv');
% delta_sei1 = readmatrix('delta_sei_73.csv');
% n_Li_s2 = readmatrix('n_Li_s_195.csv');
% delta_sei2 = readmatrix('delta_sei_195.csv');
% n_Li_s3 = readmatrix('n_Li_s_250.csv');
% delta_sei3 = readmatrix('delta_sei_250.csv');
% n_Li_s = [n_Li_s1;n_Li_s2;n_Li_s3];
% delta_sei = [delta_sei1;delta_sei2;delta_sei3];
% writematrix(n_Li_s,'n_Li_s_250.csv')
% writematrix(delta_sei,'delta_sei_250.csv')

%% 
clear all; close all; clc

L1 = readmatrix('loss_1_30.csv');
L2 = readmatrix('loss_1_20.csv');
L3 = readmatrix('loss_1_10.csv');
m = min(length(L1),length(L2));
m = min(m,length(L3));
plot(0:1:m-1,L1(1:m))
hold on
plot(0:1:m-1,L2(1:m))
hold on
plot(0:1:m-1,L3(1:m))
legend("1/30C","1/20C","1/10C")
xlabel('Charging Time(s)')
ylabel('Li-ion loss(mol)')