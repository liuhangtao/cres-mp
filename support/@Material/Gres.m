function val = Gres(obj, point)
    % residual Gibbs function [J/mol]
    val = obj.eos.Ares(point, '0') + Phys.Rg * point.T * (point.Z - 1);
end