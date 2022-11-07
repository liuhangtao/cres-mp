classdef Point < matlab.mixin.Copyable

    properties
        status      = 0
        status_eos  = 0
        status_eos0 = 0
        status_vis  = 0
        status_vis0 = 0
        status_tcx  = 0
        status_tcx0 = 0

        % basic info
        id
        mc_num
        mc_name
        mc_id

        phase = 'FL'
        X     = 1
        Xliq
        T
        P
        D
        V
        Z
        MW
        D_mass

        % critical val
        Tcr
        Pcr
        Dcr

        % reduced val
        Tr_b
        Tr
        Pr
        Dr

        % reduced val for each component (dimention: row_num = 1, col_num = mc_num)
        Tr0_b
        Tr0
        Pr0
        Dr0

        % otpt obj
        tdy    % thermodynamic properties, including ideal gas specific heat
        vis    % viscosity, including dilute gas viscosity
        tcx    % thermal conductivity, including dilute gas thermal conductivity

        point_type
        si     % user defined supporting information

        error_num = 0
        error_log

        setup  % reserved
    end

    properties (Constant)
        mc_name_hyphen = ' + '
    end
    
    methods
        function obj = Point(varargin)
            if nargin == 0
                obj.tdy = PropTdy;
                obj.vis = PropVis;
                obj.tcx = PropTcx;
            elseif nargin > 1
                obj = Point();
                obj.structPointSet(varargin);
                obj.standardizePoint;
                if mod(nargin, 2) == 1
                    obj.setup = varargin(nargin);
                end
            elseif nargin == 1
                % debug here !!!
                obj = Point();
                obj.structPointSet(varargin{1});
                obj.standardizePoint;
            end
        end

        function standardizePoint(obj)
            if ischar(obj.mc_name)
                obj.mc_name = {obj.mc_name};
                obj.mc_num = 1;
            else
                for ii = length(obj.mc_name) : -1 : 1
                    if ismissing(obj.mc_name{ii})
                        obj.mc_name(ii) = [];
                    else
                        break
                    end
                end
                obj.mc_num = length(obj.mc_name);
            end
            
            if iscell(obj.X)
                if or(iscell(obj.X{1}), isempty(obj.X{1}))
                    obj.X = [];
                    obj.logError('Point.standardizePoint: Invalid format for mole fraction.');
                    return
                else
                    tmp_X = zeros(1, obj.mc_num);
                    for ii = 1 : obj.mc_num
                        tmp_X(ii) = obj.X{ii};
                    end
                    obj.X = tmp_X;
                end
            end
                        
            %else
                if length(obj.X) > obj.mc_num
                    obj.X = obj.X(1 : obj.mc_num);
                end
                if abs(sum([0, obj.X]) - 1) > 1e-6
                    obj.logError('Point.standardizePoint: Molar fraction must sum to 1.');
                end
            %end
            
        end

        function structPointSet(obj, varargin)
            if nargin == 2
                varargin = varargin{1};
                num_prop = floor((length(varargin))/2);
            else
                num_prop = floor((nargin - 1)/2);
            end
            for i = 1 : num_prop
                obj.structPoint(char(varargin(i*2-1)), varargin(i*2));
            end
        end

        function structPoint(obj, prop_name, prop_val, prop_unit)
            if nargin == 3
                prop_unit = 1;
            else
                if or(ismissing(prop_unit), isempty(prop_unit))
                    prop_unit = 1;
                end
            end
            
            if or(numel(prop_val) == 1, all(size(prop_val) == size(obj)))
                if isa(prop_val{1, 1}, 'double')
                    % prop_val = prop_unit * cell2mat(prop_val);
                    % error when prop_val including missing data
                    % adapted to the following code:
                    prop_val_tbl = cell2table(prop_val);
                    prop_val = num2cell(prop_unit * prop_val_tbl{:,:});
                end
            elseif all(~mod(size(prop_val), size(obj)))
                if isa(prop_val{1, 1}, 'double')
                    % prop_val = prop_unit * cell2mat(prop_val);
                    % error when prop_val including missing data
                    % adapted to the following code:
                    prop_val_tbl = cell2table(prop_val);
                    prop_val = mat2cell(prop_unit * prop_val_tbl{:,:}, zeros(size(obj, 1), 1) + size(prop_val, 1)/size(obj, 1), zeros(size(obj, 2), 1) + size(prop_val, 2)/size(obj, 2));
                else
                    val_size = size(prop_val) ./ size(obj);
                    prop_val_tmp = cell(size(obj));
                    for i = 1 : size(prop_val_tmp, 1)
                        for j = 1 : size(prop_val_tmp, 2)                        
                            prop_val_tmp{i, j} = prop_val(...
                                (i - 1) * val_size(1) + 1 : i * val_size(1), ...
                                (j - 1) * val_size(2) + 1 : j * val_size(2));
                        end
                    end
                    prop_val = prop_val_tmp;
                end
            else
                errordlg('Point.structPoint: Wrong prop_val dimension.');
                return
            end

            switch upper(prop_name)
            case 'ID'
                [obj.id] = deal(prop_val{:});
            case 'MC'
                [obj.mc_name] = deal(prop_val{:});
            case 'MCID'
                [obj.mc_id] = deal(prop_val{:});
            case 'MCNUM'
                [obj.mc_num] = deal(prop_val{:});
            case 'PHASE'
                [obj.phase] = deal(prop_val{:});
            case {'X', 'XREF'}
                [obj.X] = deal(prop_val{:});
            case 'T'
                [obj.T] = deal(prop_val{:});
            case {'P', 'PREF'}
                [obj.P] = deal(prop_val{:});
            case 'D'
                [obj.D] = deal(prop_val{:});
            case 'V'
                [obj.V] = deal(prop_val{:});
            case 'Z'
                [obj.Z] = deal(prop_val{:});
            case 'MW'
                [obj.MW] = deal(prop_val{:});
            case 'DM'
                [obj.D_mass] = deal(prop_val{:});
            case 'TCR'
                [obj.Tcr] = deal(prop_val{:});
            case 'PCR'
                [obj.Pcr] = deal(prop_val{:});
            case 'DCR'
                [obj.Dcr] = deal(prop_val{:});
            case 'TR'
                [obj.Tr] = deal(prop_val{:});
            case 'PR'
                [obj.Pr] = deal(prop_val{:});
            case 'DR'
                [obj.Dr] = deal(prop_val{:});
            case 'SRES'
                tdy_obj = [obj.tdy]';
                [tdy_obj.Sres] = deal(prop_val{:});
                % useless: tdy_cell = num2cell(tdy_obj); [obj.tdy] = deal(tdy_cell{:});
                tdy_cell = num2cell(tdy_obj); [obj.tdy] = deal(tdy_cell{:});
            case 'CP0'
                tdy_obj = [obj.tdy]';
                [tdy_obj.CP0] = deal(prop_val{:});
                tdy_cell = num2cell(tdy_obj); [obj.tdy] = deal(tdy_cell{:});
            case 'CP'
                tdy_obj = [obj.tdy]';
                [tdy_obj.CP] = deal(prop_val{:});
                tdy_cell = num2cell(tdy_obj); [obj.tdy] = deal(tdy_cell{:});
            case 'CV0'
                tdy_obj = [obj.tdy]';
                [tdy_obj.CV0] = deal(prop_val{:});
                tdy_cell = num2cell(tdy_obj); [obj.tdy] = deal(tdy_cell{:});
            case 'CV'
                tdy_obj = [obj.tdy]';
                [tdy_obj.CV] = deal(prop_val{:});
                tdy_cell = num2cell(tdy_obj); [obj.tdy] = deal(tdy_cell{:});
            case 'CRATIO'
                tdy_obj = [obj.tdy]';
                [tdy_obj.Cratio] = deal(prop_val{:});
                tdy_cell = num2cell(tdy_obj); [obj.tdy] = deal(tdy_cell{:});    
            case 'VIS0'
                vis_obj = [obj.vis]';
                [vis_obj.vis0] = deal(prop_val{:});
                vis_cell = num2cell(vis_obj); [obj.vis] = deal(vis_cell{:});
            case 'VIS'
                vis_obj = [obj.vis]';
                [vis_obj.vis] = deal(prop_val{:});
                vis_cell = num2cell(vis_obj); [obj.vis] = deal(vis_cell{:});
            case 'TCX0'
                tcx_obj = [obj.tcx]';
                [tcx_obj.tcx0] = deal(prop_val{:});
                tcx_cell = num2cell(tcx_obj); [obj.tcx] = deal(tcx_cell{:});
            case 'TCX0TR'
                tcx_obj = [obj.tcx]';
                [tcx_obj.tcx0tr] = deal(prop_val{:});
                tcx_cell = num2cell(tcx_obj); [obj.tcx] = deal(tcx_cell{:});
            case 'TCX0INT'
                tcx_obj = [obj.tcx]';
                [tcx_obj.tcx0int] = deal(prop_val{:});
                tcx_cell = num2cell(tcx_obj); [obj.tcx] = deal(tcx_cell{:});
            case 'TCX'
                tcx_obj = [obj.tcx]';
                [tcx_obj.tcx] = deal(prop_val{:});
                tcx_cell = num2cell(tcx_obj); [obj.tcx] = deal(tcx_cell{:});
            case 'TCX0CETS'
                tcx_obj = [obj.tcx]';
                [tcx_obj.ce_ts_tcx] = deal(prop_val{:});
                tcx_cell = num2cell(tcx_obj); [obj.tcx] = deal(tcx_cell{:});
            case 'TCX0CECOLITG22'
                tcx_obj = [obj.tcx]';
                [tcx_obj.ce_colitg22_tcx] = deal(prop_val{:});
                tcx_cell = num2cell(tcx_obj); [obj.tcx] = deal(tcx_cell{:});    
            case 'MTRCODE'
                prop_val{:}.setPointModel(obj);
            case 'POINTTYPE'
                [obj.point_type] = deal(prop_val{:});
                si_cell = obj.setSI;
                [obj.si] = deal(si_cell{:});
            case 'IDM'
                si_obj = [obj.si]';
                [si_obj.idm] = deal(prop_val{:});
                si_cell = num2cell(si_obj); [obj.si] = deal(si_cell{:});
            case 'FILTER'
                si_obj = [obj.si]';
                [si_obj.filter] = deal(prop_val{:});
                si_cell = num2cell(si_obj); [obj.si] = deal(si_cell{:});
            case 'AUTHORFULL'
                si_obj = [obj.si]';
                [si_obj.authorfull] = deal(prop_val{:});
                si_cell = num2cell(si_obj); [obj.si] = deal(si_cell{:});
            case 'AUTHOR'
                si_obj = [obj.si]';
                [si_obj.author] = deal(prop_val{:});
                si_cell = num2cell(si_obj); [obj.si] = deal(si_cell{:});
            case 'YEAR'
                si_obj = [obj.si]';
                [si_obj.year] = deal(prop_val{:});
                si_cell = num2cell(si_obj); [obj.si] = deal(si_cell{:});
            case 'REFLIT'
                si_obj = [obj.si]';
                [si_obj.reflit] = deal(prop_val{:});
                si_cell = num2cell(si_obj); [obj.si] = deal(si_cell{:});
            case 'REFTITLE'
                si_obj = [obj.si]';
                [si_obj.reftitle] = deal(prop_val{:});
                si_cell = num2cell(si_obj); [obj.si] = deal(si_cell{:});
            case 'PROPREF'
                si_obj = [obj.si]';
                [si_obj.propname_ref] = deal(prop_val{:});
                si_cell = num2cell(si_obj); [obj.si] = deal(si_cell{:});
            case 'VALREF'
                si_obj = [obj.si]';
                [si_obj.val_ref] = deal(prop_val{:});
                si_cell = num2cell(si_obj); [obj.si] = deal(si_cell{:});
            case 'VALDEV'
                si_obj = [obj.si]';
                [si_obj.val_dev] = deal(prop_val{:});
                si_cell = num2cell(si_obj); [obj.si] = deal(si_cell{:});
            case 'VALUNCRTNTY'
                si_obj = [obj.si]';
                [si_obj.val_uncertainty] = deal(prop_val{:});
                si_cell = num2cell(si_obj); [obj.si] = deal(si_cell{:});
            case 'STATUS'
                [obj.status] = deal(prop_val{:});
            end
        end
        
        function si_cell = setSI(obj)
            si_cell = cell(size(obj));
            for i = 1 : size(obj, 1)
                for j = 1 : size(obj, 2)
                    switch upper(obj(i, j).point_type)
                    case 'REFDATA'
                        si_cell{i, j} = PropSI_RefData();
                    end
                end
            end
        end

        function logError(obj, error_info)
            obj.error_num = obj.error_num + 1;
            obj.error_log = strcat(obj.error_log, error_info, {32});
        end

   % end

   % methods (Static)
        function mc_name_str = getMcName(obj)
            if obj.mc_num == 1
                mc_name_str = obj.mc_name;
            else
                mc_name_str = obj.mc_name{1};
                for i_mc = 2 : obj.mc_num
                    mc_name_str = [mc_name_str, obj.mc_name_hyphen, obj.mc_name{i_mc}];
                end
            end
        end
    end
    
    methods (Access = protected)
        function newobj = copyElement(obj)
            newobj = copyElement@matlab.mixin.Copyable(obj);
            % newobj.si = copy(obj.si);
        end
    end

end