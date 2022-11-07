function val = P(varargin)
    % To calculate the pressure [Pa] and its partial derivatives from temperature and density.
    % format:
    %  - P(point)
    %  - P(point, pd)
    %  - P(point, pd, T, D, X)
    %  - P(       pd, T, D, X)
    % -- D can be a matrix.

    flag_point_edit = 0;
    if isa(varargin{1}, 'Point')
        flag_point_para = 1; 
        point = varargin{1};
        switch nargin
        case 1 % P(point)
            flag_point_edit = 1;
            pd = '0';
            T = point.T;
            D = point.D;
            X = point.X;
        case 2 % P(point, pd)
            pd = varargin{2};
            T = point.T;
            D = point.D;
            X = point.X;
        case 5 % P(point, pd, T, D, X)
            pd = varargin{2};
            T = varargin{3};
            D = varargin{4};
            X = varargin{5};
            if isempty(pd)
                pd = '0';
            end
            if isempty(T)
                T = point.T;
            else
                flag_point_para = 0;
            end
            if isempty(D)
                D = point.D;
            end
            if isempty(X)
                X = point.X;
            else
                flag_point_para = 0;
            end
        end

    else % P(pd, T, D, X)
        pd = varargin{1};
        T = varargin{2};
        D = varargin{3};
        X = varargin{4};

        if D < 0
            val = D;
            return
        end

    end
    
    if length(D) > 1
        disp('debug here.');
    end

    const = ModelEos_MP_CO2.const;
    para = ModelEos_MP_CO2.para;

    switch pd
    case '0' 
        val = ModelEos_MP_CO2.f_p(T, D, const, para);
        if flag_point_edit
            point.P = val;
        end
    case 'T'
        val = ModelEos_MP_CO2.f_p_t(T, D, const, para);
    case 'D' 
        val = ModelEos_MP_CO2.f_p_rho(T, D, const, para);
    case 'DD'
        dD = 1e-5;
        val = (ModelEos_MP_CO2.f_p(T, D+dD, const, para) - 2*ModelEos_MP_CO2.f_p(T, D, const, para) + ModelEos_MP_CO2.f_p(T, D-dD, const, para)) / dD^2;
    case 'n'
        val = ModelEos_MP_CO2.f_p(T, D, const, para);
        disp('ModelEos_MP_CO2: partial P / partial n unsupported.');
    end
end