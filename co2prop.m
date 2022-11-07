function varargout = co2prop(propReq,propTyp1,prop1,propTyp2,prop2)
    addpath support;
    lib_CO2 = MaterialLib_CO2;
    prop_req_seq = {'CP','Cratio','SRES','VIS0','VIS','TCX0','TCX'};
    [~, ~, ~, ~, vis, ~, tcx] = lib_CO2.calProp(prop_req_seq,'MC',{'CO2'},propTyp1,prop1,propTyp2,prop2);
    switch upper(propReq)
    case 'V'
        varargout = vis;
    case 'L'
        varargout = tcx;
    case 'VL'
        varargout = {vis, tcx};
    case 'LV'
        varargout = {tcx, vis};
    otherwise
        error('Unknown property type requested.');
    end
end
    