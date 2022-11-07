function f7 = f_p_a_tr(tr, rhor, const, para)
% 201809271453
r = 8.3144598;
tc = const(1);
rhoc = const(3);
f7 = - r * tc * rhoc * rhor / tr^2 * (1 + rhor * ...
    ModelEos_MP_CO2.f_hmhr_res_rhor(tr, rhor, para) - rhor * tr * ...
    ModelEos_MP_CO2.f_hmhr_res_rhor_tr(tr, rhor, para));