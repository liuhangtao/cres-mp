function f1 = f_p(t, rho, const, para)
% 201809271446
r = 8.3144598;
tc = const(1);
rhoc = const(3);
tr = tc / t;
rhor = rho / rhoc;
hmhr_rhor = ModelEos_MP_CO2.f_hmhr_cr_rhor(tr, rhor, const, para) + ...
    ModelEos_MP_CO2.f_hmhr_bg_rhor(tr, rhor, const, para);
f1 = rhor^2 * rhoc * r * tc / tr * hmhr_rhor;