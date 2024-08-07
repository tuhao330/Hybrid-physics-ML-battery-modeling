% compute Vs
function y = fun_Vs(x,udata)
global Ic Vsh Cb Rb Qt
y = Vsh+1/Qt*Ic*udata+Cb*Rb*Cb*1/Qt*1/Qt*Ic*(1-exp(-udata*x(2)));
end
