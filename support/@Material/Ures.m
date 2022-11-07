function val = Ures(obj, point)
    % residual internal energy [J/mol]
    val = obj.eos.Ares(point, '0') - point.T * obj.eos.Ares(point, 'T');
end