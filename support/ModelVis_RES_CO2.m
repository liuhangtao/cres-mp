classdef ModelVis_RES_CO2 < Model

    properties
        % mc_num

        % id

        % for residual entropy scaling
        coeff
        xi

        % for critical enhancement
        eos

        qC1
        qD1
        Gamma0
        kexi0
        Trefr

        % setup   % reserved
    end

    methods
        function obj = ModelVis_RES_CO2(varargin)
            % compound:
            %  - ModelVis_RES_CO2(fluid, eos, {res_para; cr_para}, res_coeff)
            % mixture:
            %  - ModelVis_RES_CO2([component_vis], eos, [])
            % note for model:
            %  * res_coeff = 
            %        [a1, a2, a3, a4, a5]
            % note for compound:
            %  * res_para =
            %        xi
            %  * cr_para = 
            %        [qC1, qD1, Gamma0, kexi0, Trefr]

            switch class(varargin{1})
            case 'Fluid'
                fluid = varargin{1};
                obj.xi    = cell2mat(varargin{3}{1});
                cr_para   = cell2mat(varargin{3}{2});
                obj.coeff = cell2mat(varargin{4}{1});
                
                obj.mc_num = fluid.mc_num;

                obj.id = fluid.id;

                obj.eos    = varargin{2};

                obj.qC1    = cr_para(:, 1);
                obj.qD1    = cr_para(:, 2);
                obj.Gamma0 = cr_para(:, 3);
                obj.kexi0  = cr_para(:, 4);
                obj.Trefr  = cr_para(:, 5);

            case 'ModelVis_RES_CO2'
                component_vis = varargin{1};

                obj.mc_num = length(component_vis);
                
                obj.id = [component_vis.id]';

                obj.coeff = component_vis(1).coeff;
                obj.xi    = [component_vis.xi]';

                obj.eos   = varargin{2};
                
                obj.qC1   = [component_vis.qC1]';
                obj.qD1   = [component_vis.qD1]';
                %obj.Gamma0, determined by EoS
                obj.kexi0 = [component_vis.kexi0]';
                obj.Trefr = max([component_vis.Trefr]);
            end
            
        end
        
        function [Vis_val, Viss_val] = Vis(obj, point)
            a1 = obj.coeff(1);
            a2 = obj.coeff(2);
            a3 = obj.coeff(3);
            a4 = obj.coeff(4);
            a5 = obj.coeff(5);

            Vis0   = point.vis.vis0;
            Viscrs = obj.Viscrs(point);
            sx     = - point.tdy.Sres / (point.X * obj.xi);

            Viss_val = exp( a1 * sx + a2 * sx^2 + a3 * sx^3 + a4 * (exp(- a5 * sx) - 1) ); 
            Visbg_val = Vis0 * Viss_val;
            Vis_val = Visbg_val * Viscrs;
            
            point.vis.visbg = Visbg_val;
            point.vis.vis = Vis_val;
            point.vis.res_kexi = point.X * obj.xi;
        end

        function val = Viscrs(obj, point)
            if point.vis.viscrs > 0
                val = point.vis.viscrs;
                return
            end
            
            qC = 1 / (point.X * obj.qC1); % [nm]
            qD = 1 / (point.X * obj.qD1); % [nm]
            kexi = point.vis.cr_kexi;

            psiD = acos( ( 1 + (qD * kexi)^2 )^(-0.5) );
            omega = sqrt(abs((qC*kexi - 1)/(qC*kexi + 1))) * tan(psiD / 2);
            if (qC * kexi < 1)
                L = 2 * atan(abs(omega));
            else
                L = log((1 + omega)/(1 - omega));
            end

            H = 1/12 * sin(3 * psiD) ...
                 - (1/4/qC/kexi) * sin(2 * psiD) ...
                 + (1/qC/kexi)^2 * (1 - 5/4*(qC * kexi)^2) * sin(psiD) ...
                 - (1/qC/kexi)^3 * ((1 - 3/2*(qC * kexi)^2) * psiD ...
                                    - (abs((qC*kexi)^2 - 1))^(3/2) * L);

            %z = 0.06667;
            z = 0.063; % ref: Luettmer-Strathmann 1995
            
            val = exp(z * H);

            point.vis.viscrs = val;
        end

    end

    methods (Static)
        function setPointModel(point)
            [point.vis] = deal(PropVis_RES_CO2());
        end
    end

end