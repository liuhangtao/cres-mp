function val = Z(obj, point)
    % compression factor [dimensionless]
    val = point.P * point.V / Phys.Rg / point.T;
    point.Z = val;
end