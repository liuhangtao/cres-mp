function f35 = f_p_t(t, rho, const, para)
% 201809272237
r = 8.3144598;
tc = const(1);
rhoc = const(3);
tr = tc / t;
rhor = rho / rhoc;
hmhr_rhor = ModelEos_MP_CO2.f_hmhr_cr_rhor(tr, rhor, const, para) + ...
    ModelEos_MP_CO2.f_hmhr_bg_rhor(tr, rhor, const, para);
hmhr_rhor_tr_cr = ModelEos_MP_CO2.f_hmhr_cr_rhor_tr(tr, rhor, const, para);
hmhr_rhor_tr_bg = ModelEos_MP_CO2.f_hmhr_bg_rhor_tr(tr, rhor, const, para);
hmhr_rhor_tr = hmhr_rhor_tr_cr + hmhr_rhor_tr_bg;
f35 = rhoc * r * rhor^2 * (hmhr_rhor - tr * hmhr_rhor_tr);