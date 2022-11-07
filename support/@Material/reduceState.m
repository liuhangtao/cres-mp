function reduceState(obj, point)
    obj.fluid.reduceState(point);
    if ismethod(obj.eos, 'reduceState') 
        obj.eos.reduceState(point);
    else
        if point.mc_num > 1
            point.logError('Mixture critical parameters are estimated.');
        end
        Tcr = point.X * obj.fluid.Tcr;
        Dcr = 1 / (point.X * (1 ./ obj.fluid.Dcr));
        Pcr = obj.eos.P('0', Tcr, Dcr, point.X);
        point.Tcr = Tcr;
        point.Pcr = Pcr;
        point.Dcr = Dcr;
        point.Tr  = point.T / Tcr;
        point.Pr  = point.P / Pcr;
        point.Dr  = point.D / Dcr;
    end
end
