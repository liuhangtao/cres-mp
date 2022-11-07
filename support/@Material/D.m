function val = D(obj, varargin)
    % molar density [mol/m3]
    % format:
    %  - D(point)
    %  - D(point, phase)
    %  - D(point, phase, T, P, X)
    %  - D(       phase, T, P, X)

    flag_point_edit = 0;
    if isa(varargin{1}, 'Point')
        point = varargin{1};
        switch nargin
        case 2 % D(point)
            flag_point_edit = 1;
            phase = point.phase;
            T = point.T;
            P = point.P;
            X = point.X;
        case 3 % D(point, phase)
            phase = varargin{2};
            if strcmp(phase, 'L&G')
                flag_point_edit = 1;
            end
            T = point.T;
            P = point.P;
            X = point.X;
        case 6 % D(point, phase, T, P, X)
            phase = varargin{2};
            T = varargin{3};
            P = varargin{4};
            X = varargin{5};
            if isempty(phase)
                phase = point.phase;
            end
            if isempty(T)
                T = point.T;
            end
            if isempty(P)
                P = point.P;
            end
            if isempty(X)
                X = point.X;
            end
        end
    else % D(phase, T, P, X)
        phase = varargin{1};
        T = varargin{2};
        P = varargin{3};
        X = varargin{4};
        point = Point('T', T, 'P', P, 'X', X);
        obj.getPointPara(point);
    end

    if P == 0
        val = 0; 
        if flag_point_edit == 1
            point.D = 0;
            point.V = Inf;
            point.Z = 1;
        end
        return
    end

    if isprop(obj.eos, 'spc')
        val = obj.eos.D(phase, T, P, X);
        if flag_point_edit
            point.D = val;
            point.V = 1 / val;
            obj.Z(point);
        end
        return
    end

    if isprop(obj.eos, 'is_cubic')
        if obj.eos.is_cubic
            if flag_point_edit
                val = obj.eos.D(point);
            else
                val = obj.eos.D(point, phase, T, P, X);
            end
            if prod(val >= 0)
                if flag_point_edit
                    point.D = val;
                    point.V = 1 / val;
                    obj.Z(point);
                end
                return
            end
        end
    end

    switch phase
    case {'L', 'Lsat', 'L, G'}
        flag_phase = 1;
    case {'G', 'Gsat', 'G, L'}
        flag_phase = 2;
    case 'SC'
        flag_phase = 0;
    otherwise
        flag_phase = [1, 2];
    end

    eps0 = 2E-8;   % digit limit = 1.3859e-08 %1.5179e-10
    epsP = 10000;   % [Pa]
    Rg   = Phys.Rg; % [J/mol/K] 

    b0   = (2^(1/3)-1)/3 * Rg * obj.fluid.Tcr ./ obj.fluid.Pcr;
    b    = X * b0;

    val = zeros(1, length(flag_phase));
    
    for i_phase = 1 : length(flag_phase)

        switch flag_phase(i_phase)
        case 1 % liquid
            D0 = 0.999 / b; % VM0 = 1.01 * b;
        case {2, 0} % gas
            D0 = 1.001 * P / Rg / T; % VM0 = 0.999 * Rg * T / P;
        end

        % Newton Itt
        %if (T > 0.95 * min(obj.fluid.Tcr)) && (T < max(obj.fluid.Tcr))
        %    Di = -101;
        %else
        itt_count = 0;
        Di = D0;
        P_cal = obj.eos.P(point, '0', [], Di, []);
        while abs(P_cal - P) > P * eps0
            Di = Di - (P_cal - P) / obj.eos.P(point, 'D', [], Di, []);

            if (Di > 1/b) || (Di < eps0) || (abs(imag(Di)) > eps0)
                Di = -101;
                break;
            end

            itt_count = itt_count + 1;
            if itt_count > 100
                Di = -101;
                break;
            end

            P_cal = obj.eos.P(point, '0', [], Di, []);
            %disp([Di, P_cal]);
        end
        %end
        
        % dichotomy 
        if (Di < 0)
            if flag_phase(i_phase) == 1 % liquid
                D1 = 0.999 / b;
                P1 = obj.eos.P(point, '0', [], D1, []);
                [D2, ~, P2, ~] = obj.spinodal(point);
                if (P1 - P)*(P2 - P) > 0
                    flag_phase(i_phase) = 2;
                end
            end
            if flag_phase(i_phase) == 2 % gas
                [~, D1, ~, P1] = obj.spinodal(point);
                D2 = epsP / Rg / T; % density of ideal gas under epsP
                P2 = epsP;
                if (P1 - P)*(P2 - P) > 0
                    flag_phase(i_phase) = 0;
                end
            end
            if flag_phase(i_phase) == 0 % supercritical
                D1 = 0.999 / b;
                P1 = obj.eos.P(point, '0', [], D1, []);
                D2 = epsP / Rg / T; % density of ideal gas under epsP
                P2 = epsP;
                if (P1 - P)*(P2 - P) > 0
                    Di = -1;
                end
            end

            if Di ~= -1
                D3 = (D1 + D2) / 2;
                P3 = obj.eos.P(point, '0', [], D3, []);
                while P3 < 0
                    D2 = D3;
                    % P2 = P3;
                    D3 = (D1 + D2) / 2;
                    P3 = obj.eos.P(point, '0', [], D3, []);
                end
                
                count = 0;
                while abs(P3 - P) > eps0 * P %[Pa]
                    if (P1 - P)*(P3 - P) > 0
                        D1 = D3;
                        P1 = P3;
                    else
                        D2 = D3;
                    end
                    D3 = (D1 + D2) / 2;
                    %disp( abs(P3 - P)/P );
                    P3 = obj.eos.P(point, '0', [], D3, []);
                    count = count + 1;
                    if count > 200
                        point.logError('Fail to solve density.');
                        break;
                    end
                end
                if P3 < 0
                    point.logError('Material.D: Fail to solve density.');
                    disp(['Material.D: Fail to solve density, point no.',num2str(point.id)]);
                    %pause;
                end
                Di = D3;
            end
            
        end

        val(i_phase) = Di;

    end

    if flag_point_edit
        point.D = val;
        point.V = 1 ./ val;
        obj.Z(point);
        if sum(val < 0)
            point.logError('Fail to solve density.');
            point.status_eos = -1;
        end
    end

end