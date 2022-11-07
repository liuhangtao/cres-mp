function val = CPres(obj, point)
    % specific heat capacity at constant pressure [J/mol/K]

    if point.D == 0
        val = 0;
    else
        P_T =   obj.eos.P(point, 'T');
        P_V = - obj.eos.P(point, 'D') * point.D^2;
        CV_res = - point.T * obj.eos.Ares(point, 'TT');
        val = CV_res - point.T * P_T^2 / P_V - Phys.Rg;
    end

end