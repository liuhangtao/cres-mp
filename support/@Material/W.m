function val = W(obj, point)
    % sound speed [m/s]
    val = sqrt(point.tdy.CP / point.tdy.CV * obj.eos.P(point, 'D') / point.MW);
    if isprop(point.tdy, 'W')
        point.tdy.W = val;
    end
end