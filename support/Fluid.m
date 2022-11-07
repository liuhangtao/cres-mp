classdef Fluid < handle

    properties
        mc_num
        id
        name
        file
        MW0
        Tb
        Ttr
        Tcr
        Pcr
        Dcr
        Vcr
        Zcr
        DPM % dipole moment
        AF % acentric factor
        Tlower
        Tupper
        Pupper
        
    end

    methods
        function obj = Fluid(varargin)
            % compound: 
            %  - Fluid(id, {data})
            % mixture: 
            %  - Fluid(component_fluids)
            %  * note: 'component_fluids' is an object array. 
            %          e.g., [fluid_1; fluid_2; ...]

            switch nargin
            case 2
                % for compound (pure fluid) ---------------------------------------------------
                obj.mc_num = 1;
                obj.id = varargin{1};
                data = varargin{2}{1};
                if size(data, 2) == 30
                    obj.name   = data{:, 1};
                    obj.file   = data{:, 5}(1 : (end - 4) ); % delete the last 4 characters ('.FLD')
                    obj.MW0    = 1e-3 * data{:, 8}; % [kg/mol]
                    obj.Tb     = data{:, 9}; % [K]
                    obj.Ttr    = data{:, 10}; % [K]
                    obj.Tcr    = data{:, 11}; % [K]
                    obj.Pcr    = 1e6 * data{:, 12}; % [Pa]
                    obj.Dcr    = 1e3 * data{:, 13}; % [mol/m3]
                    obj.Vcr    = 1 ./ obj.Dcr; % [m3/mol]
                    obj.Zcr    = obj.Pcr * obj.Vcr / Phys.Rg / obj.Tcr;
                    obj.DPM    = data{:, 14};
                    obj.AF     = data{:, 15};
                    obj.Tlower = data{:, 19}; % [K]
                    obj.Tupper = data{:, 20}; % [K]
                    obj.Pupper = 1e6 * data{:, 21}; % [Pa]
                else
                    errordlg('Fluid data: Wrong input arguments.');
                end

            case 1
                % for mixture -----------------------------------------------------------------
                component_fluids = varargin{1};
                obj.mc_num = length(component_fluids);
                obj.id     = [component_fluids.id]';
                obj.name   = {component_fluids.name}';
                obj.file   = {component_fluids.file}';
                obj.MW0    = [component_fluids.MW0]';
                obj.Tb     = [component_fluids.Tb]';
                obj.Ttr    = [component_fluids.Ttr]';
                obj.Tcr    = [component_fluids.Tcr]';
                obj.Pcr    = [component_fluids.Pcr]';
                obj.Dcr    = [component_fluids.Dcr]';
                obj.Vcr    = [component_fluids.Vcr]';
                obj.Zcr    = [component_fluids.Zcr]';
                obj.DPM    = [component_fluids.DPM]';
                obj.AF     = [component_fluids.AF]';
                obj.Tlower = [component_fluids.Tlower]';
                obj.Tupper = [component_fluids.Tupper]';
                obj.Pupper = [component_fluids.Pupper]';

            otherwise
                errordlg('Fluid: Wrong input arguments.');
            end
        end
        
        function point = reduceState(obj, point)
            point.Tr0_b = point.T ./ obj.Tb';
            point.Tr0   = point.T ./ obj.Tcr';
            point.Pr0   = point.P ./ obj.Pcr';
            point.Dr0   = point.D ./ obj.Dcr';
        end
        
        function val = cvtMolar2MassX(obj, var)
            switch class(var)
            case 'point'
                X = var.X;
            case 'double'
                X = var;
            end
            val = X .* obj.MW0' /obj.MW(X);
            dev = 1 - sum(val);
            while dev~= 0
                val(obj.mc_num) = val(obj.mc_num) + dev;
                dev = 1 - sum(val);
            end
        end

        function val = MW(obj, varargin) 
            % molar weight [g/mol][not SI]
            %  - MW()
            %  - MW(point)
            %  - MW(X)
        
            switch nargin
            case 1 % MW()
                val = obj.MW0;
            case 2
                if isa(varargin{1}, 'Point') % MW(point)
                    point = varargin{1};
                    if Bsc.isblank(point.X)
                        point.logError('Fluid.MW: missing X.');
                        X = zeros(1, point.mc_num);
                        X(1) = 1;
                    else
                        X = point.X;
                    end
                    if X(1) == 1
                        val = obj.MW0;
                    else
                        val = X * obj.MW0;
                    end
                    point.MW = val;
                else % MW(X)
                    X = varargin{1};
                    val = X * obj.MW0;
                end
            end
            
        end

    end
    
end