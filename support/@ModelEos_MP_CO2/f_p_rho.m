function f12 = f_p_rho(t, rho, const, para)
% 201904161621
%rho1 = rho + 1e-5;
%rho2 = rho - 1e-5;
%f12_ref = (ModelEos_MP_CO2.f_p(t, rho1, const, para) - ModelEos_MP_CO2.f_p(t, rho2, const, para)) / 2e-5;

r = 8.3144598;
tc = const(1);
rhoc = const(3);
tr = tc / t;
rhor = rho / rhoc;
f12 = r * tc * rhor / tr * ...
    (2 * (ModelEos_MP_CO2.f_hmhr_cr_rhor(tr, rhor, const, para) + ...
    ModelEos_MP_CO2.f_hmhr_bg_rhor(tr, rhor, const, para)) + ...
    rhor * (ModelEos_MP_CO2.f_hmhr_cr_2rhor(tr, rhor, const, para) + ...
    ModelEos_MP_CO2.f_hmhr_bg_2rhor(tr, rhor, const, para)));
end