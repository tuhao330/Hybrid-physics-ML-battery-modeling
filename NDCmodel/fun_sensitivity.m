% compute sensitivity
function y = fun_sensitivity(x,udata)
global Vsh Ic Qt
poly = [3.2  2.5898   -9.0032   18.8674  -17.8199  6.3254];
voltage_s = Vsh+1/Qt*Ic*udata+Ic*x(1)*(1-exp(-udata*x(2)));%fun_Vs(x,udata);
SoC = 1+Ic*udata/Qt/Vsh;
y(:,1) = (poly(2)+2*poly(3)*voltage_s+3*poly(4)*voltage_s.^2+4*poly(5)*voltage_s.^3+5*poly(6)*voltage_s.^4).*...
         (Ic*(1-exp(-udata*x(2)))); %  % column of x(1)
y(:,2) = (poly(2)+2*poly(3)*voltage_s+3*poly(4)*voltage_s.^2+4*poly(5)*voltage_s.^3+5*poly(6)*voltage_s.^4).*...
         (Ic*x(1)*exp(-udata*x(2)).*udata); % column of x(2) 
y(:,3) = Ic*ones(length(SoC),1); % column of beta0
y(:,4) = Ic*exp(-x(5)*SoC); % column of beta1
y(:,5) = Ic*x(4)*exp(-x(5)*SoC).*(-SoC); % column of beta2
y(:,6) = Ic*exp(-x(7)*(1-SoC)); % column of beta3
y(:,7) = Ic*x(6)*exp(-x(7)*(1-SoC)).*(SoC-1); % column of beta4
y(:,8) = Ic*(1-exp(-udata*x(9))); % column of 
y(:,9) = Ic*x(8)*exp(-udata*x(9)).*udata; % column of 
end




