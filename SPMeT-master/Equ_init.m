function [csn0,csp0] = Equ_init(p,n_Li_s,V0)
p.n_Li_s = n_Li_s;
[csn0,~] = init_cs(p,V0);
p.n_Li_s = 2.5;
csp0 = (p.n_Li_s - p.epsilon_s_n * p.L_n * p.Area * csn0) / (p.epsilon_s_p * p.L_p * p.Area);
end