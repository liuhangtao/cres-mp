function p_sat = Psat(obj, point) 
% saturated pressure at a specific temperature [Pa]

     if ~ismethod(obj.eos, 'Ares')
        p_sat = obj.eos.Psat(point);
        return
    end
    
    T = point.T;

    if (T > obj.fluid.Tcr)
        p_sat = -5;
        return;
    elseif (T == obj.fluid.Tcr)
        p_sat = obj.fluid.Pcr;
        return
    elseif (T < obj.fluid.Tlower)
        p_sat = -3;
        return;
    end

    epslK = 1E-10;
    epsPr = 1E-4;
    epsP  = 10000;
    
    % *** initial guess ***
    % Wilson initial estimate
    p_sat = obj.fluid.Pcr * exp(5.373 * (1 + obj.fluid.AF) * (1 - obj.fluid.Tcr / T));
    
    % debug here
    % [~, ~, P1, P2] = obj.spinodal(point);
    % p_sat = (max(0, P1) + P2)/2;
    
    % debug here
    % p_sat = 1000 * refpropm('P','T',T,'Q', 0, obj.fluid.file); %[Pa]
    
    % debug here
    % p_sat = T / Tcr * Pcr;

    % *** Newton Itt ***
    dP = min(1, 0.01 * p_sat); %[Pa]
    %p_step = p_sat;
    lK = log(K(obj, point, p_sat));

    itt_count = 0;

    while abs(lK) > epslK
    %while abs(p_step/p_sat) > epsPr

        lKp = log(K(obj, point, p_sat + dP/2));
        lKm = log(K(obj, point, p_sat - dP/2));
        p_step = lK * dP / (lKp - lKm);
        p_sat = p_sat - p_step;

        if (p_sat < 0)||(p_sat > obj.fluid.Pcr)
            p_sat = -1;
            break;
        end
        
        [K_factor, flag_diverge] = K(obj, point, p_sat);
        if flag_diverge == -1
            p_sat = -1;
            break;
        else
            lK = log(K_factor);
        end

        itt_count = itt_count + 1;
        
        if itt_count > 200
            %if abs(lK) > epslKr
            %    p_sat = -1;
            %end
            break;
        end
        
        dP = min(1, 0.001 * p_sat); %[Pa]
    end

    % *** dichotomy ***

    if p_sat == -1

        [~, ~, P1, P2] = obj.spinodal(T);
        if P1 < 0
            P1 = epsP;
        end
        P3 = (P1 + P2)/2;

        % if lK > 0, steady gas; if lK < 0, steady liquid

        lK1 = log(K(obj, point, P1));
        lK2 = log(K(obj, point, P2));
        lK3 = log(K(obj, point, P3));

        if lK1 < 0
            P1 = 2000;
            lK1 = log(K(obj, point, P1));
            if lK1 < 0
                return
            end
        end
        if lK2 > 0
            return
        end

        count = 0;
        
        while abs(lK3) > epslK
            if lK3 < 0
                P2 = P3;
                %lK2 = lK3;
            else
                P1 = P3;
                %lK1 = lK3;
            end
            P3 = (P1 + P2)/2;

            lK3 = log(K(obj, point, P3));
            
            count = count + 1;
            if count > 200
                break;
            end
        end
        p_sat = P3;
    end
end

function [K_val, flag_diverge] = K(obj, point, P)
    epslKr = 1E-10;

    point.P = P; % point is edited here.
    rho = obj.D(point, 'L&G');
    %V_liq = obj.V(point, P, 'L');
    %V_vap = obj.V(point, P, 'G');
    fc = obj.FC(point, rho);
    K_val = fc(1) / fc(2);

    if (abs(log(K_val)) <= epslKr) && (abs(rho(2) - rho(1))/rho(2) <= epslKr)
        flag_diverge = -1;
    else
        flag_diverge = 0;
    end
end