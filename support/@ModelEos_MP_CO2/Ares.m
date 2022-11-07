function val = Ares(varargin)
    % To calculate the residual Helmholtz function [J/mol] and its partial derivatives from temperature and density.
    % res = real - ideal
    % format:
    %  - Ares(point)
    %  - Ares(point, pd)
    %  - Ares(point, pd, T, D, X)
    %  - Ares(       pd, T, D, X)
    % -- D can be a matrix.

    % flag_point_edit = 0;
    if isa(varargin{1}, 'Point')
        flag_point_para = 1;  
        point = varargin{1};
        switch nargin
        case 1 % Ares(point)
            % flag_point_edit = 1;
            pd = '0';
            T = point.T;
            D = point.D;
            X = point.X;
        case 2 % Ares(point, pd)
            pd = varargin{2};
            T = point.T;
            D = point.D;
            X = point.X;
        case 5 % Ares(point, pd, T, D, X)
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

    if D == 0
        val = 0;
        return
    end

    const = ModelEos_MP_CO2.const;
    para = ModelEos_MP_CO2.para;
    
    Rg = Phys.Rg;
    tr = ModelEos_MP_CO2.fluid.Tcr / T;
    rhor = D / ModelEos_MP_CO2.fluid.Dcr;

    switch pd
    case '0'
        hmhr_res = ModelEos_MP_CO2.f_hmhr_cr(tr, rhor, const, para) + ModelEos_MP_CO2.f_hmhr_bg(tr, rhor, const, para) - ModelEos_MP_CO2.f_hmhr_id(tr, rhor);
        val = Rg * T * hmhr_res;
    case 'T'
        hmhr_res = ModelEos_MP_CO2.f_hmhr_cr(tr, rhor, const, para) + ModelEos_MP_CO2.f_hmhr_bg(tr, rhor, const, para) - ModelEos_MP_CO2.f_hmhr_id(tr, rhor);
        hmhr_res_tr = ModelEos_MP_CO2.f_hmhr_cr_tr(tr, rhor, const, para) + ModelEos_MP_CO2.f_hmhr_bg_tr(tr, rhor, const, para) - ModelEos_MP_CO2.f_hmhr_id_tr(tr);
        val = Rg * hmhr_res - Rg * tr * hmhr_res_tr;
    case 'TT'
        dT = 1e-7;
        Tp = T + dT;
        Tm = T - dT;
        val = (ModelEos_MP_CO2.Ares('T',Tp,D,1) - ModelEos_MP_CO2.Ares('T',Tm,D,1)) /2 /dT;
    case 'n'
        val = ModelEos_MP_CO2.Ares('0',T,D,1);
        disp('ModelEos_MP_CO2: partial Ares / partial n unsupported.');
    end
end