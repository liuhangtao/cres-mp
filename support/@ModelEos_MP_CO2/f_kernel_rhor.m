function f8 = f_kernel_rhor(tr, rhor, const, para)
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
a20 = para(7);
a21 = para(11);
a12 = para(12);
rhor_s = (1 / rhor - 1) * (1 + v1 * exp(- xishu * (1 / rhor - 1))) + ...
    d1 * (1 / tr - 1);
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
f8 = - alpha / 2 / deta1 * a20 * (1 / tr - 1)^2 * ...
    y^(- alpha / deta1 - 1) * y_rhor - (alpha - deta1) / 2 / deta1 * ...
    a21 * (1 / tr - 1)^2 * y^(- alpha / deta1) * y_rhor - 2 * a12 * ...
    (1 / tr - 1) * (rhor - 1) * ...
    (y^((gama - alpha - 2 * beta) / 2 / deta1) - 1) - ...
    (gama - alpha - 2 * beta) / 2 / deta1 * a12 * (1 / tr - 1) * ...
    (rhor - 1)^2 * y^((gama - alpha - 2 * beta) / 2 / deta1 - 1) * y_rhor;