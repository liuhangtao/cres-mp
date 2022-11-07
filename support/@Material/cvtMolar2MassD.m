function val = cvtMolar2MassD(obj, varargin)
    % convert molar density [mol/m3] to mass density [kg/m3]
    % format:
    %   - cvtMolar2MassD(point)
    %   - cvtMolar2MassD(D)
    %   - cvtMolar2MassD(D, X)

    switch class(varargin{1})
    case 'point'
        point = varargin{1};
        val = point.MW * point.D;
        point.D_mass = val;
    case 'double'
        D = varargin{1};
        switch nargin
        case 2
            val = obj.fluid.MW0 * D;
        case 3
            X = varargin{3};
            val = obj.fluid.MW(X) * D;
        end
    end

end