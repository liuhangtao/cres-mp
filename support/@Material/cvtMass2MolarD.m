function val = cvtMass2MolarD(obj, varargin)
    % convert mass density [kg/m3] to molar density [mol/m3]
    % format:
    %   - cvtMass2MolarD(point)
    %   - cvtMass2MolarD(D_mass)
    %   - cvtMass2MolarD(D_mass, X)

    switch class(varargin{1})
    case 'Point'
        point = varargin{1};
        val = point.D_mass / point.MW;
        point.D = val;
    case 'double'
        D_mass = varargin{1};
        switch nargin
        case 2
            val = D_mass / obj.fluid.MW0;
        case 3
            X = varargin{3};
            val = D_mass / obj.fluid.MW(X);
        end
    end

end