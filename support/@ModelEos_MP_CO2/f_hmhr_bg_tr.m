function f31 = f_hmhr_bg_tr(tr, rhor, const, para)
% 202104110131
r = 8.3144598;
tc = const(1);
rhoc = const(3);
hmhr_tr = ModelEos_MP_CO2.f_hmhr_res_tr(tr, 1, para) + ModelEos_MP_CO2.f_hmhr_id_tr(tr);
f31 = hmhr_tr - (1 / rhor - 1) / rhoc / r / tc * (ModelEos_MP_CO2.f_p_a(tr, 1, const, para) + tr * ModelEos_MP_CO2.f_p_a_tr(tr, 1, const, para));