% compute SoC
function y = fun_SoC(x,udata)
global Ic Vsh Qt
y = 1+Ic*udata/Qt/Vsh;
end