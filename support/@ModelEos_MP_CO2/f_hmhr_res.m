function f23 = f_hmhr_res(tr, rhor, para)
n(40) = para(8);
n(41) = - para(8);
n(42) = para(9);
n(43) = - para(9);
n(44) = para(10);
n(45) = - para(10);
d(40) = 0;
d(41) = 0;
d(42) = 1;
d(43) = 1;
d(44) = 2;
d(45) = 2;
tt(40) = 0;
tt(41) = 1;
tt(42) = 0;
tt(43) = 1;
tt(44) = 0;
tt(45) = 1;
gama(40) = 1.15;
gama(41) = 1.15;
gama(42) = 1.15;
gama(43) = 1.15;
gama(44) = 1.15;
gama(45) = 1.15;
alpha(40) = 20;
alpha(41) = 20;
alpha(42) = 20;
alpha(43) = 20;
alpha(44) = 20;
alpha(45) = 20;
beta(40) = 250;
beta(41) = 250;
beta(42) = 250;
beta(43) = 250;
beta(44) = 250;
beta(45) = 250;
epsilon(40) = 1;
epsilon(41) = 1;
epsilon(42) = 1;
epsilon(43) = 1;
epsilon(44) = 1;
epsilon(45) = 1;
n(1) = 0.38856823203161;
n(2) = 0.2938547594274 * 10;
n(3) = - 0.55867188534934 * 10;
n(4) = - 0.76753199592477;
n(5) = 0.31729005580416;
n(6) = 0.54803315897767;
n(7) = 0.12279411220335;
n(8) = 0.2165896154322 * 10;
n(9) = 0.15841735109724 * 10;
n(10) = - 0.23132705405503;
n(11) = 0.58116916431436 / 10;
n(12) = - 0.55369137205382;
n(13) = 0.48946615909422;
n(14) = - 0.24275739843501 / 10;
n(15) = 0.62494790501678 / 10;
n(16) = - 0.12175860225246;
n(17) = - 0.37055685270086;
n(18) = - 0.16775879700426 / 10;
n(19) = - 0.11960736637987;
n(20) = - 0.45619362508778 / 10;
n(21) = 0.35612789270346 / 10;
n(22) = - 0.74427727132052 / 100;
n(23) = - 0.17395704902432 / 100;
n(24) = - 0.21810121289527 / 10;
n(25) = 0.24332166559236 / 10;
n(26) = - 0.37440133423463 / 10;
n(27) = 0.14338715756878;
n(28) = - 0.13491969083286;
n(29) = - 0.2315122505348 / 10;
n(30) = 0.12363125492901 / 10;
n(31) = 0.2105832197294 / 100;
n(32) = -0.33958519026368 / 1000;
n(33) = 0.55993651771592 / 100;
n(34) = - 0.30335118055646 / 1000;
n(35) = -0.2136548868832 * 1000;
n(36) = 0.26641569149272 * 1e5;
n(37) = -0.24027212204557 * 1e5;
n(38) = -0.28341603423999 * 1000;
n(39) = 0.21247284400179 * 1000;
d(1) = 1;
d(2) = 1;
d(3) = 1;
d(4) = 1;
d(5) = 2;
d(6) = 2;
d(7) = 3;
d(8) = 1;
d(9) = 2;
d(10) = 4;
d(11) = 5;
d(12) = 5;
d(13) = 5;
d(14) = 6;
d(15) = 6;
d(16) = 6;
d(17) = 1;
d(18) = 1;
d(19) = 4;
d(20) = 4;
d(21) = 4;
d(22) = 7;
d(23) = 8;
d(24) = 2;
d(25) = 3;
d(26) = 3;
d(27) = 5;
d(28) = 5;
d(29) = 6;
d(30) = 7;
d(31) = 8;
d(32) = 10;
d(33) = 4;
d(34) = 8;
d(35) = 2;
d(36) = 2;
d(37) = 2;
d(38) = 3;
d(39) = 3;
tt(1) = 0;
tt(2) = 0.75;
tt(3) = 1;
tt(4) = 2;
tt(5) = 0.75;
tt(6) = 2;
tt(7) = 0.75;
tt(8) = 1.5;
tt(9) = 1.5;
tt(10) = 2.5;
tt(11) = 0;
tt(12) = 1.5;
tt(13) = 2;
tt(14) = 0;
tt(15) = 1;
tt(16) = 2;
tt(17) = 3;
tt(18) = 6;
tt(19) = 3;
tt(20) = 6;
tt(21) = 8;
tt(22) = 6;
tt(23) = 0;
tt(24) = 7;
tt(25) = 12;
tt(26) = 16;
tt(27) = 22;
tt(28) = 24;
tt(29) = 16;
tt(30) = 24;
tt(31) = 8;
tt(32) = 2;
tt(33) = 28;
tt(34) = 14;
tt(35) = 1;
tt(36) = 0;
tt(37) = 1;
tt(38) = 3;
tt(39) = 3;
c(8) = 1;
c(9) = 1;
c(10) = 1;
c(11) = 1;
c(12) = 1;
c(13) = 1;
c(14) = 1;
c(15) = 1;
c(16) = 1;
c(17) = 2;
c(18) = 2;
c(19) = 2;
c(20) = 2;
c(21) = 2;
c(22) = 2;
c(23) = 2;
c(24) = 3;
c(25) = 3;
c(26) = 3;
c(27) = 4;
c(28) = 4;
c(29) = 4;
c(30) = 4;
c(31) = 4;
c(32) = 4;
c(33) = 5;
c(34) = 6;
alpha(35) = 25;
alpha(36) = 25;
alpha(37) = 25;
alpha(38) = 15;
alpha(39) = 20;
beta(35) = 325;
beta(36) = 300;
beta(37) = 300;
beta(38) = 275;
beta(39) = 275;
gama(35) = 1.16;
gama(36) = 1.19;
gama(37) = 1.19;
gama(38) = 1.25;
gama(39) = 1.22;
epsilon(35) = 1;
epsilon(36) = 1;
epsilon(37) = 1;
epsilon(38) = 1;
epsilon(39) = 1;
f23 = 0;
for i = 1:7
    f23 = f23 + n(i) * rhor^d(i) * tr^tt(i);
end
for i = 8:34
    f23 = f23 + n(i) * rhor^d(i) * tr^tt(i) * exp(- rhor^c(i));
end
for i = 35: 45
    f23 = f23 + n(i) * rhor^d(i) * tr^tt(i) * exp(- alpha(i) * ...
        (rhor - epsilon(i))^2 - beta(i) * (tr - gama(i))^2);
end