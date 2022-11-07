function f3 = f_hmhr_bg_rhor(tr, rhor, const, para)
% 201809271446
r = 8.3144598;
tc = const(1);
rhoc = const(3);
f3 = tr * ModelEos_MP_CO2.f_p_a(tr, 1, const, para) / rhoc / rhor^2 / r / tc;