function VM = V(obj, varargin)
    % molar volume [m^3/mol]
    % format:
    %  - V(point)
    %  - V(point, P, phase)
    %  - V(point, P, X, phase)
    %  - V(T, P)
    %  - V(T, P, phase)
    %  - V(T, P, X)
    %  - V(T, P, X, phase)

    flag_point_edit = 0;
    switch class(varargin{1})
    case 'Point'
        point = varargin{1};
        T = point.T;
        switch nargin
        case 2
            P = point.P;
            X = point.X;
            phase = point.phase;
            flag_point_edit = 1;
        case 4
            P = varargin{2};
            X = point.X;
            phase = varargin{3};
        case 5
            P = varargin{2};
            X = varargin{3};
            phase = varargin{4};
        end

    case 'double'
        T = varargin{1};
        P = varargin{2};
        switch nargin
        case 3
            X = 1;
            phase = 'FL';
        case 4
            X = 1;
            phase = varargin{3};
        case 5
            X = varargin{3};
            phase = varargin{4};
        end

        point = Point('T',T,'P',P,'X',X);
        obj.getPointPara(point);
        
    end

    if isa(obj.eos, 'ModelEos_SRK')
        if flag_point_edit
            VM = obj.eos.V(point);
        else
            VM = obj.eos.V(point, P, X, phase);
        end
        if prod(VM > 0)
            return
        end
    end

    eps0 = 2E-10; % digit limit = 1.5179e-10
    epsP = 10000; % [Pa]
    Rg = Phys.Rg; % [J/mol/K]

    if P == 0
        VM = Inf; 
        if flag_point_edit == 1
            point.V = VM;
        end
        return
    end

    flag_phase = [1, 2];
    switch phase
    case {'L', 'Lsat', 'L, G'}
        flag_phase = 1;
    case {'G', 'Gsat', 'G, L'}
        flag_phase = 2;
    case 'SC'
        flag_phase = 0;
    end
    
    b0 = (2^(1/3)-1)/3 * Rg * obj.fluid.Tcr ./ obj.fluid.Pcr;
    b = X * b0;

    VM = zeros(1, length(flag_phase));
    
    for i_phase = 1 : length(flag_phase)

        switch flag_phase(i_phase)
        case 1 % liquid
            VM0 = 1.01 * b;
            dVM = 0.01 * b;
        case {2, 0} % gas
            VM0 = 0.999 * Rg * T / P;
            %dVM = eps * VM0; %very slow
            dVM = 0.01 * b;
        end

        % Newton Itt
        VMi = VM0;
        count = 0;
        while abs(obj.eos.P(point, VMi) - P) > P * eps0
            VMi = VMi - (obj.eos.P(point, VMi) - P) / obj.eos.P(point, VMi, 'V');
            count = count + 1;

            if (VMi > 1/eps0) || (VMi < b) || (abs(imag(VMi)) > eps0)
                VMi = -101;
                break;
            end
            if count > 100
                VMi = -101;
                break;
            end
        end
        
        % dichotomy 
        if (VMi < 0)
            if flag_phase == 1 % liquid
                V1 = 1.01 * b;
                V2 = obj.spinodal(point);
                P1 = obj.eos.P(point, V1);
                P2 = obj.eos.P(point, V2);
                if (P1 - P)*(P2 - P) > 0
                    flag_phase = 2;
                end
            end
            if flag_phase == 2 % gas
                [~, V1] = obj.spinodal(point);
                V2 = Rg * T / epsP; % specific volume of ideal gas under epsP
                P1 = obj.eos.P(point, V1);
                P2 = epsP;
                if (P1 - P)*(P2 - P) > 0
                    flag_phase = 0;
                end
            end
            if flag_phase == 0 % supercritical
                V1 = 1.01 * b;
                V2 = Rg * T / epsP; % specific volume of ideal gas under epsP
                P1 = obj.eos.P(point, V1);
                P2 = epsP;
                if (P1 - P)*(P2 - P) > 0
                    VMi = -1;
                end
            end

            if VMi ~= -1
                V3 = 2 * V1 * V2 / (V1 + V2);
                P3 = obj.eos.P(point, V3);
                while P3 < 0
                    V2 = V3;
                    P2 = P3;
                    V3 = 2 * V1 * V2 / (V1 + V2);
                    P3 = obj.eos.P(point, V3);
                end
                
                % count = 0;
                while abs(P3 - P) > eps * P %[Pa]
                    if (P1 - P)*(P3 - P) > 0
                        V1 = V3;
                        P1 = P3;
                    else
                        V2 = V3;
                    end
                    V3 = (V1 + V2)/2;
                    % disp( abs(P3 - P)/P );
                    P3 = obj.eos.P(point, V3);
                    % count = count + 1;
                end
                if P3 < 0
                    pause;
                end
                VMi = V3;

            end
            
        end

        VM(i_phase) = VMi;

    end

    if flag_point_edit
        point.V = VM;
        if sum(VM < 0)
            point.status_eos = -1;
        end
    end

end