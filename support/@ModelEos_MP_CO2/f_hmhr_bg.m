function f22 = f_hmhr_bg(tr, rhor, const, para)
% 201809271502
r = 8.3144598;
tc = const(1);
rhoc = const(3);
f22 = ModelEos_MP_CO2.f_hmhr_res(tr, 1, para) - tr * ModelEos_MP_CO2.f_p_a(tr, 1, const, para) / rhoc / ...
    r / tc * (1 / rhor - 1) + ModelEos_MP_CO2.f_hmhr_id(tr, rhor) - log(rhor);