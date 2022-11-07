function f6 = f_p_a(tr, rhor, const, para)
% 201809271447
r = 8.3144598;
tc = const(1);
rhoc = const(3);
f6 = rhoc * r * tc * rhor / tr * (1 + rhor * ...
    ModelEos_MP_CO2.f_hmhr_res_rhor(tr, rhor, para));