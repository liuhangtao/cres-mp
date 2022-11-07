function f2 = f_hmhr_cr_rhor(tr, rhor, const, para)
% 201809271458
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



try

rhor_s = (1 / rhor - 1) * (1 + v1 * exp(- xishu * (1 / rhor - 1))) + ...
    d1 * (1 / tr - 1);

catch
    
    disp('debug here.');
    
end






rhor_s_rhor = - 1 / rhor^2 + (xishu / rhor^3 - (xishu + 1) / rhor^2) * ...
    v1 * exp(- xishu * (1 / rhor - 1));
prs = 4 * (b * sqrt(rhor_s^2) / m0)^(1 / beta) + 2 * (1 / tr - 1);
prs_rhor = 4 * b / beta / m0 * rhor_s / sqrt(rhor_s^2) * rhor_s_rhor * ...
    (b * sqrt(rhor_s^2) / m0)^(1 / beta - 1);
q = sqrt((prs + sqrt(prs^2 + 12 * (1 / tr - 1)^2)) / 6 / gi);
q_rhor = prs_rhor / 12 / gi / q * ...
    (1 + prs / sqrt(prs^2 + 12 * (1 / tr - 1)^2));
y1 = 0;
y2 = 0;
for i = 1: n_y
    y1 = y1 + q^i;
    y2 = y2 + i * q^(i - 1);
end
y = (y1 / (1 + y1))^(2 * deta1);
y_rhor = 2 * deta1 * y1^(2 * deta1 - 1) / (1 + y1)^(2 * deta1 + 1) * ...
    y2 * q_rhor;
if not(or(y < 0, y > 0))
    y = 1;
end
if not(or(y_rhor < 0, y_rhor > 0))
    y_rhor = 0;
end
ka_cr = (1 / tr - 1) * y^(- alpha / 2 / deta1);
tr_cr = 1 / (1 + ka_cr);
ph_cr = (1 / rhor - 1) * y^((gama - 2 * beta) / 4 / deta1);
rhor_cr = 1 / (1 + ph_cr);
ka_cr_rhor = - alpha / 2 / deta1 * (1 / tr - 1) * ...
    y^(- alpha / 2 / deta1 - 1) * y_rhor;
tr_cr_rhor = - ka_cr_rhor / (1 + ka_cr)^2;
ph_cr_rhor = - y^((gama - 2 * beta) / 4 / deta1) / rhor^2 + ...
    (gama - 2 * beta) / 4 / deta1 * (1 / rhor - 1) * ...
    y^((gama - 2 * beta) / 4 / deta1 - 1) * y_rhor;
rhor_cr_rhor = - ph_cr_rhor / (1 + ph_cr)^2;
f2 = ModelEos_MP_CO2.f_hmhr_res_rhor(tr_cr, rhor_cr, para) * rhor_cr_rhor + ...
    (ModelEos_MP_CO2.f_hmhr_res_tr(tr_cr, rhor_cr, para) - ...
    ModelEos_MP_CO2.f_hmhr_res_tr(tr_cr, 1, para)) * tr_cr_rhor - tr_cr * ...
    ModelEos_MP_CO2.f_p_a(tr_cr, 1, const, para) / rhoc / r / tc / rhor_cr^2 * ...
    rhor_cr_rhor + tr_cr * ModelEos_MP_CO2.f_p_a_tr(tr_cr, 1, const, para) / rhoc / ...
    r / tc * (1 / rhor_cr - 1) * tr_cr_rhor + ...
    ModelEos_MP_CO2.f_p_a(tr_cr, 1, const, para) / rhoc / r / tc * (1 / rhor_cr - 1) * ...
    tr_cr_rhor + rhor_cr_rhor / rhor_cr - ...
    ModelEos_MP_CO2.f_kernel_rhor(tr, rhor, const, para);