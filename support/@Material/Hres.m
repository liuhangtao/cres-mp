function val = Hres(obj, point)
    % residual enthalpy [J/mol]
    val = obj.eos.Ares(point, '0') - T * obj.eos.Ares(point, 'T') + Phys.Rg * point.T * (point.Z - 1);
end