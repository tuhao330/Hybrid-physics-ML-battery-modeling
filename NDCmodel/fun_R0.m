% compute R0
function y = fun_R0(x,udata)
global Ic Vsh  Qt
y = x(3)+x(4)*exp(-x(5)*(1+1/Qt*Ic*udata/Vsh))+x(6)*exp(-x(7)*(-1/Qt*Ic*udata/Vsh));
end