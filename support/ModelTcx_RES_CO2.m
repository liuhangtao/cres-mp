classdef ModelTcx_RES_CO2 < Model

    properties
        % mc_num

        % id

        % for residual entropy scaling
        coeff
        zeta

        % for critical enhancement
        eos
        vis

        qD1
        Gamma0
        kexi0
        Trefr

        % setup   % reserved
    end

    methods
        function obj = ModelTcx_RES_CO2(varargin)
            % compound:
            %  - ModelTcx_RES_CO2(fluid, eos, {res_para; cr_para}, res_coeff)
            % mixture:
            %  - ModelTcx_RES_CO2([component_tcx], eos, [])
            % note for model:
            %  * res_coeff = 
            %        [a1, a2, a3, a4, a5]
            % note for compound:
            %  * res_para =
            %        zeta
            %  * cr_para = 
            %        [qD1, Gamma0, kexi0, Trefr]

            switch class(varargin{1})
            case 'Fluid'
                fluid = varargin{1};
                obj.zeta = cell2mat(varargin{3}{1});
                cr_para = cell2mat(varargin{3}{2});
                obj.coeff = cell2mat(varargin{4}{1});
                
                obj.mc_num = fluid.mc_num;

                obj.id = fluid.id;

                obj.eos    = varargin{2};

                obj.qD1    = cr_para(:, 1);
                obj.Gamma0 = cr_para(:, 2);
                obj.kexi0  = cr_para(:, 3);
                obj.Trefr  = cr_para(:, 4);

            case 'ModelTcx_RES_CO2'
                component_tcx = varargin{1};

                obj.mc_num = length(component_tcx);
                
                obj.id = [component_tcx.id]';

                obj.coeff = component_tcx(1).coeff;
                obj.zeta  = [component_tcx.zeta]';

                obj.eos   = varargin{2};
                
                obj.qD1   = [component_tcx.qD1]';
                %obj.Gamma0, determined by EoS
                obj.kexi0 = [component_tcx.kexi0]';
                obj.Trefr = max([component_tcx.Trefr]);
            end
            
        end
        
        function [Tcx_val, Tcxs_val] = Tcx(obj, point)
            a1 = obj.coeff(1);
            a2 = obj.coeff(2);
            a3 = obj.coeff(3);
            a4 = obj.coeff(4);
            a5 = obj.coeff(5);

            Tcxc   = point.tcx.tcx0tr + point.tcx.tcx0int * exp(point.tdy.Sres / Phys.Rg);
            Tcxcr = obj.Tcxcr(point);
            sx    = - point.tdy.Sres / (point.X * obj.zeta);

            Tcxs_val = exp( a1 * sx^(1/3) + a2 * sx^(2/3) + a3 * sx + a4 *  (exp(- a5 * sx) - 1) ); 
            Tcxbg_val = Tcxc * Tcxs_val;
            Tcx_val = Tcxbg_val + Tcxcr;
            
            point.tcx.tcxc = Tcxc;
            point.tcx.tcxbg = Tcxbg_val;
            point.tcx.tcx = Tcx_val;
            point.tcx.res_zeta = point.X * obj.zeta;
        end



        function val = Tcxcr(obj, point)
            % ref: Int J Thremophys (2013) 34: 191-212
            if point.tcx.tcxcr > 0
                val = point.tcx.tcxcr;
                return
            end

            % eq.(36)
            RD = 1.02;

            qD = 1 / (point.X * obj.qD1); % [nm]

            D = point.D; % [mol/m3]
            Dr = D / point.Dcr;
            CP = point.tdy.CP;
            kapa1 = 1 / point.tdy.Cratio;
            vis = point.vis.vis;
            kexi = point.tcx.cr_kexi;

            y = qD * kexi;
            Yval = 2 / pi * (((1 - kapa1)*atan(y) + kapa1 * y) - (1 - exp(-1 / (1/y + y^2 /(3*Dr^2)))));
            %Yval = 2 / pi / y * (((1 - kapa1)*atan(y) + kapa1 * y) - (1 - exp(-1 / (1/y + y^2 /(3*Dr^2))))); % add /y to denominator,
            %20220321 -- it is wrong, 20220325

            if point.mc_num == 1
                val = D * CP * RD * Phys.kB * point.T / (6 * pi * vis * 1E-09*kexi) * Yval; %[W/(m K)]
            else
                chempot_T = point.tcx.chempot_T;
                chempot_n = point.tcx.chempot_n;
                alpha_bg = point.tcx.alpha_bg;
                beta_bg = point.tcx.beta_bg;
                d_gamma_bg = point.tcx.d_gamma_bg;
                alpha_cr = Phys.kB * point.T * D * (1E-09*qD) / (6 * pi * vis) / chempot_n * Yval;
                val = D * CP * (1E-09*qD) * Phys.kB * point.T / (6 * pi * vis) * Yval ...
                    - point.T / (alpha_bg + alpha_cr) * (beta_bg + chempot_T * alpha_bg)^2 ...
                    + point.T * chempot_T^2 * alpha_bg + 2 * point.T * chempot_T * beta_bg + d_gamma_bg; %[W/(m K)] eq.(10.57)
            end

            point.tcx.tcxcr = val;
        end 

    end

    methods (Static)
        function setPointModel(point)
            [point.tcx] = deal(PropTcx_RES_CO2());
        end
    end
    
end