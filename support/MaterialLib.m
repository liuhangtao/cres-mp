classdef MaterialLib < handle

    properties
        lib_size  % number of compounds in the library  @double
        pure      % pure Material() list                @Material
        binary    % binary Material() list              @Material
        setup     % reserved
    end
    
    methods
        function obj = MaterialLib()
            obj.lib_size = 1;
            obj.pure     = Material();
            obj.binary   = Material();
        end
        
        function varargout = calProp(obj, prop_req, varargin)
            % examples:
            % -  P = lib1.calProp({'VLE'},'MC',{'R32','R1234yf'},'PHASE','Gsat','T',280,'X',{0.5,0.5})
            % -  [D, CP, Cratio, W] = lib1.calProp({'D','CP','Cratio','W'},'MC',{'R134a'},'PHASE','L','T',280,'P',0.5e6)
            
            varargin = [{'MTRCODE', obj.pure(1)}, varargin];
            point = Point(varargin);
            obj.calPointProp(prop_req, point);
            
            varargout = cell(size(prop_req));

            for i_prop = 1 : length(prop_req)
                switch upper(prop_req{i_prop})
                case 'D'
                    varargout{i_prop} = [point.D]';
                case 'V'
                    varargout{i_prop} = [point.V]';
                case 'P'
                    varargout{i_prop} = [point.P]';
                case 'VLE'
                    varargout{i_prop} = [point.P]';
                case 'CP'
                    tdy = [point.tdy]';
                    varargout{i_prop} = [tdy.CP]';
                case 'CV'
                    tdy = [point.tdy]';
                    varargout{i_prop} = [tdy.CV]';
                case 'CRATIO'
                    tdy = [point.tdy]';
                    varargout{i_prop} = [tdy.Cratio]';
                case 'W'
                    tdy = [point.tdy]';
                    varargout{i_prop} = [tdy.W]';
                case 'SRES'
                    tdy = [point.tdy]';
                    varargout{i_prop} = [tdy.Sres]';
                case 'VIS0'
                    vis = [point.vis]';
                    varargout{i_prop} = [vis.vis0]';
                case 'TCX0'
                    tcx = [point.tcx]';
                    varargout{i_prop} = {[tcx.tcx0tr]', [tcx.tcx0int]', [tcx.tcx0]'};
                case 'VIS'
                    vis = [point.vis]';
                    varargout{i_prop} = [vis.vis]';
                case 'TCX'
                    tcx = [point.tcx]';
                    varargout{i_prop} = [tcx.tcx]';
                end
            end
        end

        function calPointProp(obj, prop_req, point)
            mtr = obj.getMaterial(point);
            mtr.calProp(prop_req, point);
        end

        function mc_id = getMcid(obj, val_input)
            % To get the material id in the MaterialLib by the name.

            switch class(val_input)
            case 'cell'
                mc_name = val_input;
                mc_num = length(mc_name);
            case 'Point'
                mc_name = val_input.mc_name;
                mc_num  = val_input.mc_num;
            case 'char'
                mc_name = {val_input};
                mc_num  = 1;
            end

            mc_id = zeros(1, mc_num);
            for i = 1 : mc_num
                obj_material = findobj(obj.pure, 'name', mc_name{i});
                if isempty(obj_material)
                    tmp_pure = [obj.pure];
                    obj_material = findobj([tmp_pure.fluid], 'file', mc_name{i});
                    if isempty(obj_material)
                        errordlg(strcat('Wrong component name: ', mc_name{i}));
                    end
                end
                mc_id(i) = obj_material.id;
            end

            if isa(val_input, 'Point')
                val_input.mc_id = mc_id;
                if mc_num == 1
                    val_input.mc_name = mc_name{1};
                end
            end

        end

        function mtr = getMaterial(obj, point)
            % To get the material in the MaterialLib by mc_id.

            if isempty(point.mc_id)
                obj.getMcid(point);
            end

            switch point.mc_num
            case 1
                mtr = obj.pure(point.mc_id);
            case 2
                mtr = obj.binary(point.mc_id(1), point.mc_id(2));
            otherwise
                errordlg('Invalid mc_num.');
            end
            

        end

    end

end