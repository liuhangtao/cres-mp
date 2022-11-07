function val = Sres(obj, point)
    % residual entropy [J/mol/K]
    if isprop(obj.eos, 'spc')
        val = obj.eos.Sres(point);
    else
        val = - obj.eos.Ares(point, 'T');
    end

    if isprop(point.tdy, 'Sres')
        point.tdy.Sres = val;
    end
end