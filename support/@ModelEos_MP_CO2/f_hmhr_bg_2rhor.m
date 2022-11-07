function f14 = f_hmhr_bg_2rhor(tr, rhor, const, para)
% 201809271448
r = 8.3144598;
tc = const(1);
rhoc = const(3);
f14 = - 2 * tr * ModelEos_MP_CO2.f_p_a(tr, 1, const, para) / rhoc / rhor^3 / r / tc;