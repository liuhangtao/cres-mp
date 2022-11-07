function val = Cratio(obj, point)
    % specific heat capacity at constant pressure [J/mol/K]
    
    
    % debug here
    obj.CV(point);
    
    

    val = point.tdy.CP/point.tdy.CV;
    if isprop(point.tdy, 'Cratio')
        point.tdy.Cratio = val;
    end
end