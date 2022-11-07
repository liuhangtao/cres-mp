function f31 = f_hmhr_bg_2tr(tr, rhor, const, para)
%201809272238
r = 8.3144598;
tc = const(1);
rhoc = const(3);
hmhr_2tr = ModelEos_MP_CO2.f_hmhr_res_2tr(tr, 1, para) + ModelEos_MP_CO2.f_hmhr_id_2tr(tr);
f31 = hmhr_2tr - (1 / rhor - 1) / rhoc / r / tc * (2 * ...
    ModelEos_MP_CO2.f_p_a_tr(tr, 1, const, para) + tr * ModelEos_MP_CO2.f_p_a_2tr(tr, 1, const, para));