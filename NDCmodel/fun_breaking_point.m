% judge the breaking point of the terminal voltage
function y = fun_breaking_point(udata)
global Ntotal N_end
%% breaking point: voltage reduces to 2.5 V
for i = 4:Ntotal
    if udata(i) < udata(5)+0.05 && udata(i+1) > udata(5)+0.05
%     if udata(i) < udata(600)+0.05 && udata(i+1) > udata(600)+0.05
        N_end = i;
    end
end
y = N_end;
%% breaking point: voltage reduces to 3 V
% k = find(abs(udata-2.5)<0.005); % to find the sequence that close to 3
% N_end = k(1);
% y = N_end;
end
