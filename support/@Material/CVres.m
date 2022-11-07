function val = CVres(obj, point)
    % residual specific heat capacity at constant volume [J/mol/K]
    val = - point.T * obj.eos.Ares(point, 'TT');
end