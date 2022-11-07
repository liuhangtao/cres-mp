function f37 = f_hmhr_bg_rhor_tr(tr, rhor, const, para)
% 201809272238
r = 8.3144598;
tc = const(1);
rhoc = const(3);
f37 = (ModelEos_MP_CO2.f_p_a(tr, 1, const, para) + tr * ModelEos_MP_CO2.f_p_a_tr(tr, 1, const, para)) / ...
    rhoc / r / tc / rhor^2;