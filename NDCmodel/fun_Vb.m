% compute Vb
function y = fun_Vb(x,udata)
global Ic Vsh  Cb Cs Rb Qt
y = Vsh+1/Qt*Ic*udata-Cs*Rb*Cb*1/Qt*1/Qt*Ic*(1-exp(-udata*x(2)));
end