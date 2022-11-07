function f24 = f_kernel(tr, rhor, const, para)
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
f24 = a20 / 2 * (1 / tr - 1)^2 * (y^(- alpha / deta1) - 1) + a21 / 2 * ...
    (1 / tr - 1)^2 * (y^(- (alpha - deta1) / deta1) - 1) - a12 * ...
    (1 / tr - 1) * (rhor - 1)^2 * ...
    (y^((gama - alpha - 2 * beta) / 2 / deta1) - 1);