function varargout = calProp(obj, prop_req, point)

%    if isa(point.si, 'PropSI_RefData')    
%        if point.si.filter < 0
%            varargout = num2cell(point.si.filter * ones(size(prop_req)));
%            return
%        end
%    end

    if 0 %point.T < max([obj.fluid.Tlower; 180])
        point.logError('Unsupported range.');
        return
    end 
    
    obj.getPointPara(point);
    if ~Bsc.isblank(point.D_mass)
        if Bsc.isblank(point.D)
            obj.cvtMass2MolarD(point);
        end
    end
    obj.calTherm(point);

    varargout = cell(size(prop_req));

    for i_prop = 1 : length(prop_req)
        switch upper(prop_req{i_prop})
        case 'D'
            varargout{i_prop} = point.D;
        case 'V'
            varargout{i_prop} = point.V;
        case 'P'
            varargout{i_prop} = point.P;    
        case 'VLE'
            varargout{i_prop} = point.P;
        case 'CP'
            varargout{i_prop} = obj.CP(point);
        case 'CV'
            varargout{i_prop} = obj.CV(point);
        case 'CRATIO'
            varargout{i_prop} = obj.Cratio(point);
        case 'W'
            varargout{i_prop} = obj.W(point);
        case 'SRES'
            varargout{i_prop} = obj.Sres(point);
        case 'VIS0'
            varargout{i_prop} = obj.vis0.Vis0(point);
        case 'TCX0'
            varargout{i_prop} = obj.tcx0.Tcx0(point);
        case 'VIS'
            if isa(obj.vis, 'ModelVis_RES_CO2')
                obj.Trp_RES_CO2_pre(point, 'VIS');
            end
            varargout{i_prop} = obj.vis.Vis(point);
        case 'TCX'
            if isa(obj.tcx, 'ModelTcx_RES_CO2')
                obj.Trp_RES_CO2_pre(point);
            end
            varargout{i_prop} = obj.tcx.Tcx(point);
        end
    end

end