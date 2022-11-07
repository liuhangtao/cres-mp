function f33 = f_hmhr_id_tr(tr)
a0 = zeros(8, 1);
s0 = zeros(8, 1);
a0(1) = 8.37304456;
a0(2) = -3.70454304;
a0(3) = 2.5;
a0(4) = 1.99427042;
a0(5) = 0.62105248;
a0(6) = 0.41195293;
a0(7) = 1.04028922;
a0(8) = 0.08327678;
s0(4) = 3.15163;
s0(5) = 6.1119;
s0(6) = 6.77708;
s0(7) = 11.32384;
s0(8) = 27.08792;
f33 = a0(2) + a0(3)/tr; 
for i=4:8
    f33 = f33  - (a0(i) * s0(i) * exp( ( - s0(i)*tr ) ) ) / ( exp( ( - s0(i) * tr) ) - 1);
end

