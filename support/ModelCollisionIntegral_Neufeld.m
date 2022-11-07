classdef ModelCollisionIntegral_Neufeld < Model

    % ref: Philip D. Neufeld; J. Chem. Phys. 57, 1100 (1972); https://doi.org/10.1063/1.1678363

    properties (Constant)
        para = [
            1.06036, 0.15610, 0.19300, 0.47635, 1.03587, 1.52996, 1.76474, 3.89411,  0,      0,        0,       0;        % (1, 1) % 1
            1.00220, 0.15530, 0.16105, 0.72751, 0.86125, 2.06848, 1.95162, 4.84492,  0,      0,        0,       0;        % (1, 2) % 2
            0.96573, 0.15611, 0.44067, 1.52420, 2.38981, 5.08063, 0,       0,       -5.373, 19.2866,  -1.30775, 6.58711;  % (1, 3) % 3
            0.93447, 0.15578, 0.39478, 1.85761, 2.45988, 6.15727, 0,       0,        4.246, 12.9880,  -1.36399, 3.33290;  % (1, 4) % 4
            0.90972, 0.15565, 0.35967, 2.18528, 2.45169, 7.17936, 0,       0,       -3.814,  9.38191,  0.14025, 9.93802;  % (1, 5) % 5
            0.88928, 0.15562, 0.33305, 2.51303, 2.36298, 8.11690, 0,       0,       -4.649,  9.86928,  0.12851, 9.82414;  % (1, 6) % 6
            0.87208, 0.15568, 0.36583, 3.01399, 2.70659, 9.92310, 0,       0,       -4.902, 10.2274,   0.12306, 9.97712;  % (1, 7) % 7
            1.16145, 0.14874, 0.52487, 0.77320, 2.16178, 2.43787, 0,       0,       -6.435, 18.0323,  -0.76830, 7.27371;  % (2, 2) % 8
            1.11521, 0.14796, 0.44844, 0.99548, 2.30009, 3.06031, 0,       0,        4.565, 38.5868,  -0.69403, 2.56375;  % (2, 3) % 9
            1.08228, 0.14807, 0.47128, 1.31596, 2.42738, 3.90018, 0,       0,       -5.623,  3.08449,  0.28271, 3.22871;  % (2, 4) % 10
            1.05581, 0.14822, 0.51203, 1.67007, 2.57317, 4.85939, 0,       0,       -7.120,  4.71210,  0.21730, 4.73530;  % (2, 5) % 11
            1.03358, 0.14834, 0.53928, 2.01942, 2.72350, 5.84817, 0,       0,       -8.576,  7.66012,  0.15493, 7.60110;  % (2, 6) % 12
            1.05567, 0.14980, 0.30887, 0.86437, 1.35766, 2.44123, 1.29030, 5.55734,  2.339, 57.7757,  -1.08980, 6.94750;  % (3, 3) % 13
            1.02621, 0.15050, 0.55381, 1.40070, 2.06176, 4.26234, 0,       0,        5.227, 11.3331,  -0.82090, 3.87185;  % (3, 4) % 14
            0.99958, 0.15029, 0.50441, 1.64304, 2.06947, 4.87712, 0,       0,       -5.184,  3.45031,  0.26821, 3.73348;  % (3, 5) % 15
            1.12007, 0.14578, 0.53347, 1.11986, 2.28803, 3.27567, 0,       0,        7.427, 21.0480,  -0.28759, 6.69149]; % (4, 4) % 16
        %   A,       B,       C,       D,       E,       F,       G,       H,        R*1E4,  S,        W,       P
        %   1        2        3        4        5        6        7        8         9      10        11       12
    end

    methods (Static)
        function [val0, valS] = Omega(l, s, Ts, pd)
            if nargin == 3
                pd = 0;
            end
            i_row_list = [
                 1,  2,  3,  4,  5,  6,  7;
                -5,  8,  9, 10, 11, 12, -5;
                -5, -5, 13, 14, 15, -5, -5;
                -5, -5, -5, 16, -5, -5, -5];
            i_row = i_row_list(l, s);
            if i_row == -5
                errordlg('ModelCollisionIntegral: Invalid order of the collision integral.');
            end

            switch pd
            case 0
                A = ModelCollisionIntegral_Neufeld.para(i_row,  1);
                B = ModelCollisionIntegral_Neufeld.para(i_row,  2);
                C = ModelCollisionIntegral_Neufeld.para(i_row,  3);
                D = ModelCollisionIntegral_Neufeld.para(i_row,  4);
                E = ModelCollisionIntegral_Neufeld.para(i_row,  5);
                F = ModelCollisionIntegral_Neufeld.para(i_row,  6);
                G = ModelCollisionIntegral_Neufeld.para(i_row,  7);
                H = ModelCollisionIntegral_Neufeld.para(i_row,  8);
                R = 1E-04 * ModelCollisionIntegral_Neufeld.para(i_row,  9);
                S = ModelCollisionIntegral_Neufeld.para(i_row, 10);
                W = ModelCollisionIntegral_Neufeld.para(i_row, 11);
                P = ModelCollisionIntegral_Neufeld.para(i_row, 12);
                val0 = A * Ts.^(- B) + C * exp(- D * Ts) + E * exp(- F * Ts) + G * exp(- H * Ts);
                valS = val0 + R * Ts.^B .* sin(S * Ts.^W - P);
            case 1
                [O_l_s_val0,  O_l_s_valS]  = ModelCollisionIntegral_Neufeld.Omega(l, s,     Ts, 0);
                [O_l_sp_val0, O_l_sp_valS] = ModelCollisionIntegral_Neufeld.Omega(l, s + 1, Ts, 0);
                val0 = (s + 2) ./ Ts .* (O_l_sp_val0 - O_l_s_val0);
                valS = (s + 2) ./ Ts .* (O_l_sp_valS - O_l_s_valS);
            case 2
                [O_l_s_val0,   O_l_s_valS]   = ModelCollisionIntegral_Neufeld.Omega(l, s,     Ts, 0);
                [O_l_sp_val0,  O_l_sp_valS]  = ModelCollisionIntegral_Neufeld.Omega(l, s + 1, Ts, 0);
                [O_l_spp_val0, O_l_spp_valS] = ModelCollisionIntegral_Neufeld.Omega(l, s + 2, Ts, 0);
                val0 = (s + 2) * (s + 3) ./ Ts .* (O_l_spp_val0 - 2 * O_l_sp_val0 + O_l_s_val0);
                valS = (s + 2) * (s + 3) ./ Ts .* (O_l_spp_valS - 2 * O_l_sp_valS + O_l_s_valS);
            end
        end

    end
end