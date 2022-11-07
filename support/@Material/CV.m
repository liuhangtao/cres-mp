function val = CV(obj, point)
    % specific heat capacity at constant volume [J/mol/K]

    if isprop(obj.eos, 'spc')
        val = obj.eos.CV(point);
    else
        val = obj.CVres(point) + obj.eos0.CP0(point) - Phys.Rg;
        if isprop(point.tdy, 'CV')
            point.tdy.CV = val;
        end
    end
end