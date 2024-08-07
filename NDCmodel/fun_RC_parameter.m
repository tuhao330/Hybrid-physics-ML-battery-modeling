% compute Cb, Cs, and Rb
function y = fun_RC_parameter(x)
global Qt
a = 1/Qt/x(1)/x(2);
y(1) = Qt/(a+1); % Cb
y(2) = a*Qt/(a+1); % Cs
y(3) = Qt/x(2)/y(1)/y(2); % Rb  
end

