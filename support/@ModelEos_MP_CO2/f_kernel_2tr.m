function f32 = f_kernel_2tr(tr, rhor, const, para)
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
rhor_s_tr = - d1 / tr^2;
rhor_s_2tr = 2 * d1 / tr^3;
prs = 4 * (b * sqrt(rhor_s^2) / m0)^(1 / beta) + 2 * (1 / tr - 1);
prs_tr = 4 * b / beta / m0 * (b * sqrt(rhor_s^2) / m0)^(1 / beta - 1) * ...
    rhor_s / sqrt(rhor_s^2) * rhor_s_tr - 2 / tr^2;
prs_2tr = 4 * b^2 / beta / m0^2 * (1 / beta - 1) * ...
    (b * sqrt(rhor_s^2) / m0)^(1 / beta - 2) * rhor_s_tr^2 + 4 * b / ...
    beta / m0 * (b * sqrt(rhor_s^2) / m0)^(1 / beta - 1) * rhor_s / ...
    sqrt(rhor_s^2) * rhor_s_2tr + 4 / tr^3;
q = sqrt((prs + sqrt(prs^2 + 12 * (1 / tr - 1)^2)) / 6 / gi);
q_tr = (prs_tr + (prs * prs_tr + 12 * (1 / tr^2 - 1 / tr^3)) / ...
    sqrt(prs^2 + 12 * (1 / tr - 1)^2)) / 12 / q / gi;
q_2tr = (prs_2tr - (prs * prs_tr + 12 * (1 / tr^2 - 1 / tr^3))^2 / ...
    (prs^2 + 12 * (1 / tr - 1)^2)^1.5 + (prs * prs_2tr + ...
    prs_tr^2 + 12 * (3 / tr^4 - 2 / tr^3)) / ...
    sqrt(prs^2 + 12 * (1 / tr - 1)^2)) / 12 / q / gi - q_tr^2 / q;
y1 = 0;
y2 = 0;
y3 = 0;
for i = 1: n_y
    y1 = y1 + q^i;
    y2 = y2 + i * q^(i - 1);
    y3 = y3 + i * (i - 1) * q^(i - 2);
end
y = (y1 / (1 + y1))^(2 * deta1);
y_tr = 2 * deta1 * y1^(2 * deta1 - 1) / (1 + y1)^(2 * deta1 + 1) * ...
    y2 * q_tr;
y_2tr = 2 * deta1 * y2 * y1^(2 * deta1 - 1) / ...
    (1 + y1)^(2 * deta1 + 1) * q_2tr + 2 * deta1 * y3 * ...
    y1^(2 * deta1 - 1) / (1 + y1)^(2 * deta1 + 1) * q_tr^2 + 2 * ...
    deta1 * (2 * deta1 - 1) * y2^2 * y1^(2 * deta1 - 2) / ...
    (1 + y1)^(2 * deta1 + 1) * q_tr^2 - 2 * deta1 * (2 * deta1 + 1) * ...
    y2^2 * y1^(2 * deta1 - 1) / (1 + y1)^(2 * deta1 + 2) * q_tr^2;
if not(or(y < 0, y > 0))
    y = 1;
end
if not(or(y_tr < 0, y_tr > 0))
    y_tr = 0;
end
if not(or(y_2tr < 0, y_2tr > 0))
    y_2tr = 0;
end
f32 = a20 * (3 / tr^4 - 2 / tr^3) * (y^(- alpha / deta1) - 1) - 2 * ...
    a20 * (1 / tr^2 - 1 / tr^3) * alpha / deta1 * ...
    y^(- alpha / deta1 - 1) * y_tr + a20 / 2 * (1 / tr - 1)^2 * alpha / ...
    deta1 * (alpha / deta1 + 1) * y^(- alpha / deta1 - 2) * y_tr^2 - ...
    a20 / 2 * (1 / tr - 1)^2 * alpha / deta1 * ...
    y^(- alpha / deta1 - 1) * y_2tr + a21 * (3 / tr^4 - 2 / tr^3) * ...
    (y^(- (alpha - deta1) / deta1) - 1) - 2 * (alpha - deta1) / deta1 * ...
    a21 * (1 / tr^2 - 1 / tr^3) * y^(- alpha / deta1) * y_tr + ...
    (alpha - deta1) / 2 / deta1 * alpha / deta1 * a21 * ...
    (1 / tr - 1)^2 * y^(- alpha / deta1 - 1) * y_tr^2 - ...
    (alpha - deta1) / 2 / deta1 * a21 * (1 / tr - 1)^2 * ...
    y^(- alpha / deta1) * y_2tr - 2 * a12 / tr^3 * (rhor - 1)^2 * ...
    (y^((gama - alpha - 2 * beta) / 2 / deta1) - 1) + ...
    (gama - alpha - 2 * beta) / deta1 * a12 / tr^2 * (rhor - 1)^2 * ...
    y^((gama - alpha - 2 * beta) / 2 / deta1 - 1) * y_tr - ...
    (gama - alpha - 2 * beta) / 2 / deta1 * ...
    ((gama - alpha - 2 * beta) / 2 / deta1 - 1) * a12 * (1 / tr - 1) * ...
    (rhor - 1)^2 * y^((gama - alpha - 2 * beta) / 2 / deta1 - 2) * ...
    y_tr^2 - (gama - alpha - 2 * beta) / 2 / deta1 * a12 * ...
    (1 / tr - 1) * (rhor - 1)^2 * ...
    y^((gama - alpha - 2 * beta) / 2 / deta1 - 1) * y_2tr;