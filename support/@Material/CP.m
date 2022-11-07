function val = CP(obj, point)
    % specific heat capacity at constant pressure [J/mol/K]
    
    if isprop(obj.eos, 'spc')
        val = obj.eos.CP(point);
    else
        val = obj.CPres(point) + obj.eos0.CP0(point);
        if isprop(point.tdy, 'CP')
            if val < 0
                point.logError('CP < 0.');
            end
            point.tdy.CP = val;
        end
    end

end