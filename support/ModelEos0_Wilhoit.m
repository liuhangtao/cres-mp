classdef ModelEos0_Wilhoit < Model

    properties
        % id % fluid id
        range % temperature range of the cp0 correlation
        a0
        a1
        a2
        a3
        a4
        a5
        a6
        a7
        Tref  % refrence temperature for representation of enthalpy and entropy
        I     % integration constant for calculation of enthalpy increments
        J     % integration constan for calculation of entropies
        s00   % entropy of the ideal gas at T = 298.15 K and p = 100 kPa
    end

    methods
        function obj = ModelEos0_Wilhoit(fluid, data)
            if nargin == 2
                obj.id = fluid.id;
                data = cell2mat(data{1});
                if size(data, 2) == 14
                    obj.range = [min(data(1:2)), max(data(1:2))];
                    obj.a0    = data(3);
                    obj.a1    = data(4) * 1e5;
                    obj.a2    = data(5);
                    obj.a3    = data(6);
                    obj.a4    = data(7);
                    obj.a5    = data(8) * 1e6;
                    obj.a6    = data(9);
                    obj.a7    = data(10);
                    obj.Tref  = data(11);
                    obj.I     = data(12);
                    obj.J     = data(13);
                    obj.s00   = data(14);
                else
                    errordlg('ModelEos0_Wilhoit data: Wrong input arguments.');
                end     
            else
                errordlg('ModelEos0_Wilhoit: Wrong input arguments.');
            end
        end
        
        function val = CP0(obj, point)
            T = [point.T]';
            if or(min(T) < obj.range(1), max(T) > obj.range(2))
                disp(strcat('Fluid #', num2str(obj.id),'out of the range of cp0 correlation! T/K: ', num2str([min(T), max(T)])));
            end
            y = (T > obj.a7).*(T - obj.a7)./(T + obj.a6);
            val = Phys.Rg * (obj.a0 + obj.a1./T.^2 .* exp(-obj.a2./T) + obj.a3 * y.^2 + (obj.a4 - obj.a5./(T - obj.a7).^2).*y.^8); 
            point.tdy.CP0 = val;
        end

        function val = H0(obj, point)
            T = [point.T]';
            x = exp(-obj.a2 ./ T);
            y = (T <= obj.a7).*(T - obj.a7)./(T + obj.a6);
            h = (obj.a6 + obj.a7) * ((2 * obj.a3 + 8 * obj.a4) * log(1 - y) + (obj.a3 * (1 + 1./(1 - y)) + obj.a4 * (7 + 1./(1 - y))).*y + ...
                obj.a4 * (3 * y.^2 + 5/3 * y.^3 + y.^4 + 3/5 * y.^5 + 1/3 * y.^6) + 1/7 * (obj.a4 - obj.a5/(obj.a6 + obj.a7)^2) * y.^7);
            val = Phys.Rg * T .* (obj.a0 + obj.a1 * x./(obj.a2 * T) + (obj.I + h)./ T);  
        end

        function val = S0(obj, point)
            % s0 is related to (T, p) or (T, rho)
            % here, s0 is calculated according to (T, rho)
            
            T = [point.T]';
            x = exp(-obj.a2 ./ T);
            y = (T <= obj.a7).*(T - obj.a7)./(T + obj.a6);
            z = T ./ (T + obj.a6) * (1 + obj.a6/obj.a7);
            tmp_ele = ((obj.a4 * obj.a7^2 - obj.a5)/obj.a6^2 * (- obj.a7 / obj.a6) .^ (6 - (1:7)) - obj.a4) .* y.^(1:7) ./(1:7);
            tmp_sum = sum(tmp_ele, 2);
            s = (obj.a3 + (obj.a4*obj.a7^2 - obj.a5)/obj.a6^2 * (obj.a7/obj.a6)^4) * (obj.a7/obj.a6)^2 * log(z) + ...
                (obj.a3 + obj.a4) * log((T + obj.a6)/(obj.a6 + obj.a7)) + ...
                tmp_sum - ...
                (obj.a3 * (obj.a6 + obj.a7) / obj.a6 + obj.a5 * y.^6 / 7 / obj.a7 / (obj.a6 + obj.a7)) .* y;
            
            % s0 at T = point.T, p = 100 kPa: 
            s0_100 = Phys.Rg * (obj.J + obj.a0 * log(T) + obj.a1/obj.a2^2 * (1 + obj.a2 ./ T) * x + s);

            s0_TP = s0_100 - Phys.Rg * log(point.P / 100e3);
            val = s0_TP - Phys.Rg * log(point.Z);
        end

    end
    
end