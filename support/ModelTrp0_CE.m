classdef ModelTrp0_CE < Model

    properties
        % mc_num

        % id
        sigmaij
        epsilonij
        deltaij
        alpha0
        mij
        eos0
        fint_coeff

        setup_collision_integral = {'Nfd0'}
        setup_LJpara_mix         = 0
        setup_delta_mix          = 0
        setup_Tcx_fint           = 0
    end

    methods
        function obj = ModelTrp0_CE(varargin)
            % compound:
            %  - ModelTrp0_CE(fluid, eos0)
            %  - ModelTrp0_CE(fluid, eos0, {LJ_para; polar_para; fint_coeff})
            %  - ModelTrp0_CE(fluid, eos0, {LJ_para; polar_para; fint_coeff}, {code_collision_integral; code_pure; code_mix})
            % mixture:
            %  - ModelTrp0_CE(fluid, eos0)
            %  - ModelTrp0_CE([component_trp], eos0)
            %  - ModelTrp0_CE([component_trp], eos0, {binary_para})
            %  - ModelTrp0_CE([component_trp], eos0, {binary_para}, {code_mix})
            % ----------------------------------------------------------------------------------------------
            % note for compound:
            %  * LJ_para = 
            %        []
            %        [sigma0, epsilon0]
            %  * polar_para = 
            %        []
            %        [delta0, alpha0]
            %  * fint_coeff = 
            %        []
            %        [c0, c1, c2, p0, p1, p2]
            %  * code_collision_integral options of the collision integral model:
            %        'Nfd0' -- Neufeld (1972), sine term not included
            %        'NfdS' -- Neufeld (1972), sine term included
            %        'Brk'  -- Brokaw (1969), for polar fluid
            %        'Mck'  -- Monchick (1961), for polar fluid, correlated by Xiang (2002)
            %        'LN'   -- Liu (2021), correlated for noble gas transport property using ab initio data
            %  * code_pure = [0, 0]
            %    code_pure(1) parameters sigma amd epsilon in the LJ potential function:
            %        0 -- calculated from fluid basic properties
            %        1 -- calculated from fluid basic properties
            %        2 -- read from excel file
            %    code_pure(2) f_int in the internal term of dilute gas thermal conductivity:
            %        0 -- calculated from collision integral
            %        1 -- calculated from collision integral
            %        2 -- calculated by a polynomial, coefficients read from excel file
            % ----------------------------------------------------------------------------------------------
            % note for mixture:
            %  * binary_para = 
            %        {[sigmaij], [epsilonij]}
            %  * code_mix = 0
            %    code_mix(1) parameters sigma amd epsilon of mixtures in the LJ potential function:
            %        0 -- calculated by mixing rule using pure LJ potential parameters
            %        1 -- calculated by mixing rule using pure LJ potential parameters
            %        2 -- read from excel file

            switch class(varargin{1})
            case 'Fluid'
                fluid = varargin{1};
                obj.mc_num = fluid.mc_num;

                obj.id = fluid.id;

                switch nargin
                case 2 % ModelTrp0_CE(fluid, eos0)
                    code_pure = [0, 0];
                    code_mix  = 0;
                    code_collision_integral = {'Nfd0'};
                    polar_para = [];
                case 3 % ModelTrp0_CE(fluid, eos0, {LJ_para; polar_para; fint_coeff})
                    LJ_para = cell2mat(varargin{3}{1});
                    polar_para = cell2mat(varargin{3}{2});
                    fint_coeff = cell2mat(varargin{3}{3});
                    if isempty(LJ_para)
                        code_pure = [0, 0];
                        code_mix  = 0;
                    else
                        if isempty(fint_coeff)
                            code_pure = [2, 0];
                        else
                            if sum(fint_coeff) == 0
                                code_pure = [2, 0];
                            else
                                code_pure = [2, 2];
                            end
                        end
                        code_mix  = -1;
                    end
                    code_collision_integral = {'Nfd0'};
                case 4 % ModelTrp0_CE(fluid, eos0, {LJ_para; polar_para; fint_coeff}, {code_collision_integral; code_pure; code_mix})
                    LJ_para = cell2mat(varargin{3}{1});
                    polar_para = cell2mat(varargin{3}{2});
                    fint_coeff = cell2mat(varargin{3}{3});
                    if numel(varargin{4}) == 3
                        code_collision_integral = varargin{4}{1};
                        code_pure = cell2mat(varargin{4}{2});
                        code_mix  = cell2mat(varargin{4}{3});
                        if code_pure(2) ~= 0
                            if isempty(fint_coeff)
                                code_pure(2) = 0;
                            else
                                if sum(fint_coeff) == 0
                                    code_pure(2) = 0;
                                end
                            end
                        end
                    else
                        if isempty(LJ_para)
                            code_pure = [0, 0];
                            code_mix  = 0;
                        else
                            if isempty(fint_coeff)
                                code_pure = [2, 0];
                            else
                                if sum(fint_coeff) == 0
                                    code_pure = [2, 0];
                                else
                                    code_pure = [2, 2];
                                end
                            end
                            code_mix  = -1;
                        end
                        code_collision_integral = {'Nfd0'};
                    end
                end
                
                obj.eos0 = varargin{2};
                obj.setup_collision_integral = code_collision_integral; 
                obj.setup_Tcx_fint = code_pure(2);
                obj.setup_LJpara_mix = code_mix;

                switch code_pure(1) % parameters sigma amd epsilon in the LJ potential function
                case {0, 1} % calculated from fluid basic properties
                    obj.epsilonij = fluid.Tcr/1.2593; %[K]
                    obj.sigmaij   = 0.809E-09./fluid.Dcr.^(1/3); %[m]
                case 2
                    obj.sigmaij   = 1E-09 * LJ_para(:, 1); %[m]
                    obj.epsilonij = LJ_para(:, 2); %[K]
                end
                obj.mij = fluid.MW / Phys.NA;

                if isempty(polar_para)
                    obj.deltaij = 1e-49 * 0.5 * fluid.DPM.^2 ./ obj.epsilonij / Phys.kB ./ obj.sigmaij.^3;
                    obj.alpha0 = zeros(obj.mc_num, 1);
                else
                    obj.deltaij = polar_para(:, 1);
                    obj.alpha0 = polar_para(:, 2);
                end

                if code_pure(2) == 2 % f_int in the internal term of dilute gas thermal conductivity
                    obj.fint_coeff = fint_coeff;
                end

                if obj.mc_num > 1
                    if obj.setup_LJpara_mix == 0
                        obj.sigmaij = diag(obj.sigmaij);
                        obj.epsilonij = diag(obj.epsilonij);
                        obj.mij = diag(obj.mij);
                        for i = 1 : obj.mc_num
                            for j = (i + 1) : obj.mc_num
                                if (obj.deltaij(i) * obj.deltaij(j) == 0) && (obj.deltaij(i) + obj.deltaij(j) > 0)
                                    if obj.deltaij(i) == 0
                                        i_nonpolar = i;
                                        i_polar = j;
                                    else
                                        i_nonpolar = j;
                                        i_polar = i;
                                    end
                                    kexi = 1 + 0.25 * obj.alpha0(i_nonpolar)/obj.sigma0(i_nonpolar)^3 * obj.DPM(i_polar)/sqrt(obj.epsilon0(i_polar) * obj.sigma0(i_polar)^3) * sqrt(obj.epsilon0(i_polar) / obj.epsilon0(i_nonpolar));
                                else
                                    kexi = 1;
                                end
                                obj.sigmaij(i, j) = (obj.sigmaij(i, i) + obj.sigmaij(j, j))/2 * kexi^(-1/6);
                                obj.epsilonij(i, j) = sqrt(obj.epsilonij(i, i) * obj.epsilonij(j, j)) * kexi^2;
                                obj.mij(i, j) = 2 * obj.mij(i, i) * obj.mij(j, j) / (obj.mij(i, i) + obj.mij(j, j));
                            end
                        end
                        obj.sigmaij = obj.sigmaij + tril(obj.sigmaij', -1);
                        obj.epsilonij = obj.epsilonij + tril(obj.epsilonij', -1);
                        obj.mij = obj.mij + tril(obj.mij', -1);
                    else
                        errordlg('Debug ModelTrp0_CE.');
                    end

                    if obj.setup_delta_mix == 0
                        obj.deltaij = diag(obj.deltaij);
                        for i = 1 : obj.mc_num
                            for j = (i + 1) : obj.mc_num
                                obj.deltaij(i, j) = sqrt(obj.deltaij(i, i) * obj.deltaij(j, j)) * (sqrt(obj.sigmaij(i, i) * obj.sigmaij(j, j)) / obj.sigmaij(i, j))^3;
                            end
                        end
                        obj.deltaij = obj.deltaij + tril(obj.deltaij', -1);
                    else
                        errordlg('Debug ModelTrp0_CE.');
                    end
                end

            case 'ModelTrp0_CE'
                component_trp = varargin{1};
                obj.mc_num    = length(component_trp);
                
                obj.id         = [component_trp.id]';
                obj.fint_coeff = cat(1, component_trp.fint_coeff);

                obj.eos0 = varargin{2};
                obj.setup_collision_integral = cat(1, component_trp.setup_collision_integral);
                obj.setup_collision_integral = obj.setup_collision_integral(1);
                obj.setup_Tcx_fint = cat(1, component_trp.setup_Tcx_fint);
                obj.setup_Tcx_fint = obj.setup_Tcx_fint(1);
                % debug here.

                switch nargin
                case 2 % ModelTrp0_CE([component_trp], eos0)
                    if obj.setup_LJpara_mix == -1
                        obj.setup_LJpara_mix = 0;
                    end
                case 3 % ModelTrp0_CE([component_trp], eos0, {binary_para})
                    binary_para = varargin{3}{1};
                    if obj.setup_LJpara_mix == -1
                        if isempty(binary_para)
                            obj.setup_LJpara_mix = 0;
                        else
                            obj.setup_LJpara_mix = 2;
                        end
                    end
                case 4 % ModelTrp0_CE([component_trp], eos0, {binary_para}, {code_mix})
                    binary_para = varargin{3}{1};
                    obj.setup_LJpara_mix = cell2mat(varargin{4}{1});
                end

                switch obj.setup_LJpara_mix
                case {0, 1}
                    obj.sigmaij = diag([component_trp.sigmaij]);
                    obj.epsilonij = diag([component_trp.epsilonij]);
                    obj.mij = diag([component_trp.mij]);
                    for i = 1 : obj.mc_num
                        for j = (i + 1) : obj.mc_num
                            obj.sigmaij(i, j) = (obj.sigmaij(i, i) + obj.sigmaij(j, j))/2;
                            obj.epsilonij(i, j) = sqrt(obj.epsilonij(i, i) * obj.epsilonij(j, j));
                            obj.mij(i, j) = 2 * obj.mij(i, i) * obj.mij(j, j) / (obj.mij(i, i) + obj.mij(j, j));
                        end
                    end
                    obj.sigmaij = obj.sigmaij + tril(obj.sigmaij', -1);
                    obj.epsilonij = obj.epsilonij + tril(obj.epsilonij', -1);
                    obj.mij = obj.mij + tril(obj.mij', -1);
                case 2
                    obj.sigmaij = binary_para{1};
                    obj.epsilonij = binary_para{2};
                    obj.mij = diag([component_trp.mij]);
                    for i = 1 : obj.mc_num
                        for j = (i + 1) : obj.mc_num
                            obj.mij(i, j) = 2 * obj.mij(i, i) * obj.mij(j, j) / (obj.mij(i, i) + obj.mij(j, j));
                        end
                    end
                    obj.mij = obj.mij + tril(obj.mij', -1);
                otherwise
                    errordlg('ModelTrp0_CE code_mix: Wrong input arguments.');
                end

                switch obj.setup_delta_mix
                case {0, 1}
                    obj.deltaij = diag([component_trp.deltaij]);
                    for i = 1 : obj.mc_num
                        for j = (i + 1) : obj.mc_num
                            obj.deltaij(i, j) = sqrt(obj.deltaij(i, i) * obj.deltaij(j, j)) * (sqrt(obj.sigmaij(i, i) * obj.sigmaij(j, j)) / obj.sigmaij(i, j))^3;
                        end
                    end
                    obj.deltaij = obj.deltaij + tril(obj.deltaij', -1);
                case 2
                    obj.deltaij = binary_para{3};
                end

            otherwise
                errordlg('ModelTrp0_CE: Wrong input arguments.');
                
            end
            
        end
        
        function getPointPara(obj, point)
            % No point parameter needed.
        end
        
        function Vis0_val = Vis0(obj, point, code)
            Tsij  = [point.T] ./ [obj.epsilonij];
            O22ij = obj.Omega(2, 2, Tsij);
            Visij = 5/16 * sqrt(Phys.kB * [point.T] .* [obj.mij] / pi) ./ [obj.sigmaij].^2 ./O22ij;
            if nargin == 2
                code = 0;
            end
            if code == 1
                Vis0_val = Visij;
                return
            end

            if obj.mc_num == 1
                Vis0_val = Visij'; % in vertical
            % elseif obj.mc_num == 2
            %     O11ij = obj.Omega(1, 1, Tsij);
            %     A = O22ij(1, 2) ./ O11ij(1, 2);
            %     Xeta = sum(point.X' * point.X ./ Visij, 'all');
            %     Yeta = 3/5 * A * ( ...
            %             point.X(1)^2 / Visij(1, 1) * obj.mij(1, 1)/obj.mij(2, 2) + ...
            %             point.X(2)^2 / Visij(2, 2) * obj.mij(2, 2)/obj.mij(1, 1) + ...
            %             2 * point.X(1) * point.X(2) * ...
            %                 (obj.mij(1, 1)+obj.mij(2, 2))^2 /obj.mij(1, 1)/obj.mij(2, 2)/4 * ...
            %                 Visij(1, 2) / Visij(1, 1) / Visij(2, 2) );
            %     Zeta = 3/5 * A * ( ...
            %             point.X(1)^2 * obj.mij(1, 1)/obj.mij(2, 2) + ...
            %             point.X(2)^2 * obj.mij(2, 2)/obj.mij(1, 1) + ...
            %             2 * point.X(1) * point.X(2) * ...
            %                 ((obj.mij(1, 1)+obj.mij(2, 2))^2 /obj.mij(1, 1)/obj.mij(2, 2)/4 * ...
            %                 Visij(1, 2) * (1 / Visij(1, 1) + 1 / Visij(2, 2)) - 1) );
            %     Vis0_val = (1 + Zeta)/(Xeta + Yeta);
            else
                O11ij = obj.Omega(1, 1, Tsij);
                Aij = O22ij ./ O11ij;
                Hij = zeros(obj.mc_num);
                for ii = 1 : obj.mc_num
                    Hii_tmp = zeros(obj.mc_num, 1);
                    for kk = 1 : obj.mc_num
                        if kk ~= ii
                            Hii_tmp(kk) = 2 * point.X(ii) * point.X(kk) / Visij(ii, kk) * ...
                                obj.mij(ii, ii) * obj.mij(kk, kk) / (obj.mij(ii, ii) + obj.mij(kk, kk))^2 * ...
                                (5/3 /Aij(ii, kk) + obj.mij(kk, kk) / obj.mij(ii, ii));
                        end
                    end
                    Hij(ii, ii) = point.X(ii)^2 / Visij(ii, ii) + sum(Hii_tmp);
                    for jj = (ii + 1) : obj.mc_num
                        Hij(ii, jj) = - 2 * point.X(ii) * point.X(jj) / Visij(ii, jj) * ...
                            obj.mij(ii, ii) * obj.mij(jj, jj) / (obj.mij(ii, ii) + obj.mij(jj, jj))^2 * ...
                            (5/3 /Aij(ii, jj) - 1);
                    end
                end
                Hij = Hij + tril(Hij', -1);
                Vis0_val = - det(Bsc.eliminateZero([Hij, [point.X]'; [point.X], 0])) / det(Bsc.eliminateZero(Hij));
            end
            
            point.structPoint('VIS0', num2cell(Vis0_val));
            % point.vis.vis0 = Vis0_val;
        end
        
        function varargout = Tcx0(obj, point, code)
            if nargin == 2
                code = '0';
            end
            
            Tsij  = [point.T] ./ [obj.epsilonij];
            O22ij = obj.Omega(2, 2, Tsij);

            % translational contribution
            Tcx0trij = 75/64 * sqrt(Phys.kB^3 * [point.T] ./ [obj.mij] / pi) ./ [obj.sigmaij].^2 ./O22ij;

            % internal contribution
            switch obj.setup_Tcx_fint
            case {0, 1}
                fint = 1.2 * O22ij ./ obj.Omega(1, 1, Tsij);
            case 2
                fint = obj.Fint(Tsij);
            end
            if fint == 0
                Tcx0intj = zeros(size(point));
            else
                if obj.mc_num == 1
                    Tcx0intj = obj.Vis0(point) .* (obj.eos0.CP0(point) / Phys.NA - 2.5 * Phys.kB) .* fint ./ obj.mij;
                else
                    Tcx0intj = diag(obj.Vis0(point, 1)) .* (obj.eos0.CP0(point, 1) / Phys.NA - 2.5 * Phys.kB) .* diag(fint) ./diag(obj.mij);
                end
            end

            % translational contribution
            O22ij_ref = zeros(1, length(point));
            if isprop(point(1).si, 'propname_ref')
                if ~Bsc.isblank(point(1).si.propname_ref)
                    if strcmp(point(1).si.propname_ref, 'TCX')
                        tmp_si = [point.si];
                        tcx0_ref = [tmp_si.val_ref];
                        O22ij_ref = 75/64 * sqrt(Phys.kB^3 * [point.T] ./ [obj.mij] / pi) ./ [obj.sigmaij].^2 ./ tcx0_ref;
                    end 
                end
            end

            % pure compound: result output
            if obj.mc_num == 1
                switch code
                case '0'
                    Tcx0_val = Tcx0trij' + Tcx0intj';
                    varargout = {Tcx0_val, Tcx0trij', Tcx0intj'};
                    point.structPoint('TCX0', num2cell(Tcx0_val)); % point.tcx.tcx0 = Tcx0_val
                    point.structPoint('TCX0TR', num2cell(Tcx0trij')); % point.tcx.tcx0tr = Tcx0trij
                    point.structPoint('TCX0INT', num2cell(Tcx0intj')); % point.tcx.tcx0int = Tcx0intj
                    point.structPoint('TCX0CETS', num2cell(Tsij'));
                    point.structPoint('TCX0CECOLITG22', num2cell(O22ij_ref'));
                case 'tr'
                    varargout = {Tcx0trij'};
                    point.structPoint('TCX0TR', num2cell(Tcx0trij')); % point.tcx.tcx0tr = Tcx0trij
                    point.structPoint('TCX0CETS', num2cell(Tsij'));
                    point.structPoint('TCX0CECOLITG22', num2cell(O22ij_ref'));
                case 'int'
                    varargout = {Tcx0intj'};
                    point.structPoint('TCX0INT', num2cell(Tcx0intj')); % point.tcx.tcx0int = Tcx0intj
                end
                return
            end

            O11ij = obj.Omega(1, 1, Tsij);
            O12ij = obj.Omega(1, 2, Tsij);
            O13ij = obj.Omega(1, 3, Tsij);
            Aij = O22ij ./ O11ij;
            Cij = O12ij ./ O11ij;
            Bij = (5 * O12ij - 4 * O13ij) ./ O11ij;

            Deltaij = zeros(obj.mc_num);
            for ii = 1 : obj.mc_num
                for jj = (ii + 1) : obj.mc_num
                    if obj.mij(ii, ii) >= obj.mij(jj, jj)
                        X1 = point.X(ii) / (point.X(ii) + point.X(jj));
                        delta_c = obj.mij(jj, jj) / obj.mij(ii, ii);
                        delta_a = sqrt(2)/8 / (1 + 1.8 * delta_c)^2 * O11ij(ii, jj) / O22ij(jj, jj);
                    else
                        X1 = point.X(jj) / (point.X(ii) + point.X(jj));
                        delta_c = obj.mij(ii, ii) / obj.mij(jj, jj);
                        delta_a = sqrt(2)/8 / (1 + 1.8 * delta_c)^2 * O11ij(ii, jj) / O22ij(ii, ii);
                    end
                    delta_b = 10 * delta_a * (1 + 1.8 * delta_c + 3 * delta_c ^ 2) - 1;
                    Deltaij(ii, jj) = 1.3 * (6 * Cij(ii, jj) - 5)^2 * delta_a * X1 / (1 + delta_b * X1);
                end
            end
            Deltaij = Deltaij + tril(Deltaij', -1);

            % mixture: translational contribution
            Lij = zeros(obj.mc_num);
            for ii = 1 : obj.mc_num
                Lii_tmp = zeros(obj.mc_num, 1);
                for kk = 1 : obj.mc_num
                    if kk ~= ii
                        Lii_tmp(kk) = 2 * point.X(ii) * point.X(kk) * (...
                            15/2 * obj.mij(ii, ii)^2 + 25/4 * obj.mij(kk, kk)^2 ...
                            - 3 * obj.mij(kk, kk)^2 * Bij(ii, kk) ...
                            + 4 * obj.mij(ii, ii) * obj.mij(kk, kk) * Aij(ii, kk)) ...
                            / (obj.mij(ii, ii) + obj.mij(kk, kk))^2 ...
                            / Aij(ii, kk) / Tcx0trij(ii, kk) / (1 + Deltaij(ii, kk));
                    end
                end
                Lij(ii, ii) = -4 * point.X(ii)^2 / Tcx0trij(ii, ii) -sum(Lii_tmp);
                for jj = (ii + 1) : obj.mc_num
                    Lij(ii, jj) = 2 * point.X(ii) * point.X(jj) * obj.mij(ii, ii) * obj.mij(jj, jj) ...
                        / (obj.mij(ii, ii) + obj.mij(jj, jj))^2 / Aij(ii, jj) / Tcx0trij(ii, jj) / (1 + Deltaij(ii, jj)) ...
                        * (55/4 - 3 * Bij(ii, jj) - 4 * Aij(ii, jj));
                end
            end
            Lij = Lij + tril(Lij', -1);
            Tcx0tr = 4 * det(Bsc.eliminateZero([Lij, [point.X]'; [point.X], 0])) / det(Bsc.eliminateZero(Lij));

            if strcmp(code, 'tr')
                varargout = {Tcx0tr};
                point.structPoint('TCX0TR', num2cell(Tcx0tr')); % point.tcx.tcx0tr = Tcx0tr_val
                return
            end

            % mixture: internal contribution
            % ref: Vesovic (2001), Eq. (8), typographical error fixed
            tmp_Tcx0int = zeros(obj.mc_num, 1);
            for uu = 1 : obj.mc_num
                tmp_Tcx0int(uu) = Tcx0intj(uu) * point.X(uu) / sum([point.X] * Aij(uu, uu) ./ Aij(uu, :) * Tcx0trij(uu, uu) ./ Tcx0trij(uu, :));
            end
            Tcx0int = sum(tmp_Tcx0int);

            switch code
            case '0'
                Tcx0_val = Tcx0tr' + Tcx0int';
                varargout = {Tcx0_val, Tcx0tr', Tcx0int'};
                point.structPoint('TCX0', num2cell(Tcx0_val)); % point.tcx.tcx0 = Tcx0_val
                point.structPoint('TCX0TR', num2cell(Tcx0tr')); % point.tcx.tcx0tr = Tcx0tr_val
                point.structPoint('TCX0INT', num2cell(Tcx0int')); % point.tcx.tcx0int = Tcx0int_val
            case 'int'
                varargout = {Tcx0int'};
                point.structPoint('TCX0INT', num2cell(Tcx0int')); % point.tcx.tcx0int = Tcx0int_val
            end
        end

        function [O22_val, Ts_val] = refO22FromVis0(obj, varargin)
            % convert {vis0 v.s. T} data to {Omega22 v.s. Ts}
            switch nargin
            case 2
                point = varargin{1};
                if strcmp(point(1).si.propname_ref, 'VIS')
                    T = [point.T]';
                    tmp_si = [point.si];
                    vis0 = [tmp_si.val_ref]';
                else
                    errordlg('Wrong point.propname_ref.');
                end
            case 3
                T = varargin{1};
                vis0 = varargin{2};
            end
            O22_val = 5/16 * sqrt(Phys.kB * T .* obj.mij / pi) ./ obj.sigmaij.^2 ./vis0;
            Ts_val  = T ./ obj.epsilonij;
        end

        function [fint_val, Ts_val] = refFintFromTcx0(obj, varargin)
            % convert {tcx0 v.s. T} data to {fint v.s. Ts}
            switch nargin
            case 2
                point = varargin{1};
                if strcmp(point(1).si.propname_ref, 'TCX')
                    T = [point.T]';
                    tmp_si = [point.si];
                    tcx0 = [tmp_si.val_ref]';
                else
                    errordlg('Wrong point.propname_ref.');
                end
            case 3
                T = varargin{1};
                tcx0 = varargin{2};
            end

            fint_val = (tcx0 - obj.Tcx0(point, 'tr'))./obj.Vis0(point) .* obj.mij ./ (obj.eos0.CP0(point) / Phys.NA - 2.5 * Phys.kB);
            Ts_val  = T ./ obj.epsilonij;
        end

        function val = Omega(obj, l, s, Ts, pd)
            if nargin == 4
                pd = 0;
            end
            switch obj.setup_collision_integral{1}
            case 'Nfd0'
                val = ModelCollisionIntegral_Neufeld.Omega(l, s, Ts, pd);
            case 'NfdS'
                [~, val] = ModelCollisionIntegral_Neufeld.Omega(l, s, Ts, pd);
            case 'Brk'
                val = ModelCollisionIntegral_Brokaw.Omega(l, s, Ts, obj.deltaij, pd);
            case 'Mck'
                val = ModelCollisionIntegral_Monchick.Omega(l, s, Ts, obj.deltaij, pd);
            case 'LN'
                val = ModelCollisionIntegral_Liu_noble.Omega(l, s, Ts, pd);
            end
        end
        
        function val = Fint(obj, Ts)
            if obj.mc_num == 1
                val = obj.fint_coeff(1) * Ts.^ obj.fint_coeff(4) ...
                    + obj.fint_coeff(2) * Ts.^ obj.fint_coeff(5) ...
                    + obj.fint_coeff(3) * Ts.^ obj.fint_coeff(6);
            else
                val = obj.fint_coeff(:, 1) .* diag(Ts).^ obj.fint_coeff(:, 4) ...
                    + obj.fint_coeff(:, 2) .* diag(Ts).^ obj.fint_coeff(:, 5) ...
                    + obj.fint_coeff(:, 3) .* diag(Ts).^ obj.fint_coeff(:, 6);
                val = diag(val);
            end
        end

    end
end