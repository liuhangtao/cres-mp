function f21 = f_hmhr_cr(tr, rhor, const, para)
% 201809271503
r = 8.3144598;
tc = const(1);
rhoc = const(3);
alpha = const(4);
beta = const(5);
gama = const(6);
deta1 = const(7);
n_y = const(8);
v1 = para(1);
xishu = para(2);
d1 = para(3);
b = para(4);
m0 = para(5);
gi = para(6);
rhor_s = (1 / rhor - 1) * (1 + v1 * exp(- xishu * (1 / rhor - 1))) + ...
    d1 * (1 / tr - 1);
prs = 4 * (b * sqrt(rhor_s^2) / m0)^(1 / beta) + 2 * (1 / tr - 1);
q = sqrt((prs + sqrt(prs^2 + 12 * (1 / tr - 1)^2)) / 6 / gi);
y1 = 0;
for i = 1: n_y
    y1 = y1 + q^i;
end
y = (y1 / (1 + y1))^(2 * deta1);
if not(or(y < 0, y > 0))
    y = 1;
end
ka_cr = (1 / tr - 1) * y^(- alpha / 2 / deta1);
tr_cr = 1 / (1 + ka_cr);
ph_cr = (1 / rhor - 1) * y^((gama - 2 * beta) / 4 / deta1);
rhor_cr = 1 / (1 + ph_cr);
f21 = ModelEos_MP_CO2.f_hmhr_res(tr_cr, rhor_cr, para) - ModelEos_MP_CO2.f_hmhr_res(tr_cr, 1, para) + ...
    log(rhor_cr) + tr_cr * ModelEos_MP_CO2.f_p_a(tr_cr, 1, const, para) / rhoc / r / ...
    tc * (1 / rhor_cr - 1) - ModelEos_MP_CO2.f_kernel(tr, rhor, const, para);