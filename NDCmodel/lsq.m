function sum = lsq(x)
global SoC V Vmin
sum = 0;
for i = 1:length(SoC)
    func = (x(1)*SoC(i)^2 + x(2)*SoC(i) + Vmin*x(3))/(SoC(i)^3 + x(4)*SoC(i)^2 + x(5)*SoC(i) + x(3));
    sum = sum + (func - V(i))^2;
end
end