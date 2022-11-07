function [D1, D2, P1, P2] = spinodal(obj, varargin)
% To calculate the extremum of pressure on the isotherm line.
%  -  P_minimum on isothermal
%  -- D1 [mol/m3]
%  -- P1 [Pa]
%  -  P_maximum on isothermal
%  -- D2 [mol/m3]
%  -- P2 [Pa]
% format:
%  - spinodal(point)
%  - spinodal(T)
%  - spinodal(T, X)

    switch class(varargin{1})
    case 'Point' % spinodal(point)
        point = varargin{1};
        T = point.T;
        X = point.X;

    case 'double'
        switch nargin
        case 2 % spinodal(T)
            T = varargin{1};
            X = 1;
        case 3 % spinodal(T, X)
            T = varargin{1};
            X = varargin{2};
        end
        
        point = Point('MTRCODE', obj, 'MC', obj.name, 'T', T, 'X', X);
        obj.getPointPara(point);
    end

    if obj.mc_num == 1
        Dcr = obj.fluid.Dcr;
    else
        Dcr = (X * obj.fluid.Vcr.^(2/3)) / (X * obj.fluid.Vcr.^(5/3)); 
        % Chueh-Prausnitz method. ref: Tong Jingshan (1982), eq. (4-5.8)
    end
    
    try
    if (obj.eos.P(point, 'D', [], Dcr, []) >= 0) % debug here for mixture
        D1 = Dcr; 
        D2 = Dcr;
        P1 = min(obj.fluid.Pcr); 
        P2 = max(obj.fluid.Pcr);
        return;
    end
    catch
        errordlg('Material.spinodal: debug here at line 50.');
        pause
    end

    eps0 = 1E-02;
    epsP = 10000;
    Rg   = Phys.Rg;

    if isa(obj.eos, 'Eos_SRK')
        b = point.tdy.b;
        a = point.tdy.a;
        c = point.tdy.c;
    else
        b0 = (2^(1/3)-1)/3 * Rg * obj.fluid.Tcr ./ obj.fluid.Pcr;
        b = X * b0;
    end

    % Newton Itt
    D01 = Dcr * (1 + 0.01 * b * Dcr); % V01 = Vcr - 0.01 * b;
    D02 = Dcr * (1 - 0.01 * b * Dcr); % V02 = Vcr + 0.01 * b;
    
    if isa(obj.eos, 'Eos_SRK')
        % to solve rho by setting partial(p)/partial(rho) = 0
        D0 = roots([- b^3, 0, b * (b * Rg*T/a + 3), 2 * (b * Rg*T/a - 1), Rg*T/a]); % [mol/m^3]
        D0 = real(D0(abs(imag(D0)) <= eps0));
        if ~isempty(D0)
            D0 = D0(D0 < 1 / b);
            for i = 1 : length(D0)
                P_DD = obj.eos.P(point, 'DD', [], D0(i), []);
                if abs(P_DD) > eps0
                    if P_DD > 0
                        if D0(i) > Dcr
                            D01 = D0(i);
                        end
                    else
                        if D0(i) > Dcr
                            D02 = D0(i);
                        end
                    end
                end
            end
        end
    end

    itt_count = 0;
    D1 = D01;
    P_D = obj.eos.P(point, 'D', [], D1, []);
    while abs(P_D) > eps0
        D1 = D1 - P_D / obj.eos.P(point, 'DD', [], D1, []);
        if (D1 < Dcr) || (D1 > 1/b) || (abs(imag(D1)) > eps0)
            D1 = -1;
            break;
        end

        itt_count = itt_count + 1;
        if itt_count > 100
            D1 = -1;
            break;
        end
    end

    itt_count = 0;
    D2 = D02;
    P_D = obj.eos.P(point, 'D', [], D2, []);
    while abs(P_D) > eps0
        D2 = D2 - P_D / obj.eos.P(point, 'DD', [], D2, []);
        if (D2 < eps0) || (D2 > Dcr) || (abs(imag(D2)) > eps0)
            D2 = -1;
            break;
        end
        
        itt_count = itt_count + 1;
        if itt_count > 100
            D2 = -1;
            break;
        end
    end

    % dichotomy 
    if D1 < 0
        D11 = 0.96 / b; % V11 = 1.05 * b;
        P11 = obj.eos.P(point, 'D', [], D11, []);
        D12 = Dcr + 1e-07;  %V12 = Vcr - 2e-10 * b;
        P12 = obj.eos.P(point, 'D', [], D12, []);
        if P11 * P12 > 0
            errordlg('Debug here! 11232040');
        end
        D1  = (D11 + D12)/2;
        P1  = obj.eos.P(point, 'D', [], D1, []);
        itt_count = 0;
        while abs(P1) > eps0
            if P11 * P1 > 0
                D11 = D1;
                P11 = P1;
            else
                D12 = D1;
            end
            D1 = (D11 + D12)/2;
            if (D1 <= Dcr) || (D1 > 1/b) || (abs(imag(D1)) > eps0)
                D1 = -1;
                break;
            end

            itt_count = itt_count + 1;
            if itt_count > 200
                D1 = -1;
                break;
            end
            
            P1 = obj.eos.P(point, 'D', [], D1, []);
        end
    end

    if D2 < 0
        D21 = Dcr - 1e-07; %V21 = Vcr + 2e-10 * b;
        P21 = obj.eos.P(point, 'D', [], D21, []);
        D22 = epsP / Rg / T; %V22 = Rg * T / epsP;
        P22 = epsP; %P22 = obj.eos.P(point, V22, 'V');
        if P21 * P22 > 0
            errordlg('Debug here! 11232041');
        end
        D2  = (D21 + D22)/2;
        P2  = obj.eos.P(point, 'D', [], D2, []);
        itt_count = 0;
        while abs(P2) > eps0
            if P21 * P2 > 0
                D21 = D2;
                P21 = P2;
            else
                D22 = D2;
            end
            D2 = (D21 + D22)/2;
            if (D2 < eps0) || (D2 >= Dcr) || (abs(imag(D2)) > eps0)
                D2 = -1;
                break;
            end
            
            itt_count = itt_count + 1;
            if itt_count > 200
                D2 = -1;
                break;
            end

            P2 = obj.eos.P(point, 'D', [], D2, []);
        end
    end

    % isotherm search, not finished!!!
    if D1 < 0
        
        
        % errordlg('unfinished!!!');
        
        
        D = 0.96 / b;

        step_D = 0.01 / b;
        while obj.eos.P(point, '0', [], D - step_D, []) < obj.eos.P(point, '0', [], D, [])
            D = D - step_D;
            if D < Dcr
                break;
            end
        end
        step_D = 0.01 * step_D;
        while obj.eos.P(point, '0', [], D + step_D, []) < obj.eos.P(point, '0', [], D, [])
            D = D + step_D;
            if D > 0.96 * b
                break;
            end
        end
        step_D = 0.01 * step_D;
        while obj.eos.P(point, '0', [], D - step_D, []) < obj.eos.P(point, '0', [], D, []) % very slow for 300 K. debug here.
            D = D - step_D;
            if D < Dcr
                D = D + step_D; % only for the last attempt.
                break;
            end
        end

        D1 = D;
    end

    if D2 < 0
        
        
        % errordlg('unfinished!!!');
        
        
        
        D = Dcr;
        step_D = 0.05 / b;
        while obj.eos.P(point, '0', [], D - step_D, []) < obj.eos.P(point, '0', [], D, [])
            D = D - step_D;
            if D < 10000 / Rg / T % specific volume of ideal gas under 10 kPa 
                break;
            end
        end
        step_D = 0.10 * step_D;
        while obj.eos.P(point, '0', [], D + step_D, []) < obj.eos.P(point, '0', [], D, [])
            D = D + step_D;
            if D > Dcr
                break;
            end
            % P = obj.eos.P(point, V);
        end
        step_D = 0.02 * step_D; %slow
        while obj.eos.P(point, '0', [], D - step_D, []) < obj.eos.P(point, '0', [], D, [])
            D = D - step_D;
            if D < 10000 / Rg / T
                break;
            end
        end
        step_D = 0.01 * step_D;
        while obj.eos.P(point, '0', [], D + step_D, []) < obj.eos.P(point, '0', [], D, [])
            D = D + step_D;
            if D > Dcr
                D = Dcr - step_D; % only for the last attempt.
                break;
            end
            % p = obj.eos.P(point, V);
        end    
        
        D2 = D;
    end

    P1 = obj.eos.P(point, '0', [], D1, []);
    P2 = obj.eos.P(point, '0', [], D2, []);

end
