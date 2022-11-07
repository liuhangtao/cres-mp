function val = FC(obj, varargin)
% fugacity coefficient (i.e. {phi}) [dimensionless]
% format:
%  - FC(point)
%  - FC(point, D)
%  - FC(point, D, X)

% for pure fluid: (f   /P     ) = exp( Ares(T, D) / RT - ln(Z) + (Z - 1) );
% for mixture   : (f_i /P /X_i) = exp( Ares(T, D, X, 'n') / RT - ln(Z) ).
% res = real - ideal;
% Ares(..., 'n') = partial(Ares) / partial(n_i).

    switch nargin
    case 2 % FC(point)
        point = varargin{1};
        T = point.T;
        D = point.D;
        Z = point.Z;
        X = point.X;
        mc_num = point.mc_num;
    case 3 % FC(point, D)
        point = varargin{1};
        T = point.T;
        P = point.P;
        D = varargin{2};
        Z = P ./ D / Phys.Rg / T;
        X = point.X;
        mc_num = point.mc_num;
    case 4 % FC(point, D, X)
        point = varargin{1};
        T = point.T;
        P = point.P;
        D = varargin{2};
        Z = P ./ D / Phys.Rg / T;
        X = varargin{3};
        mc_num = point.mc_num;
    end

    Rg = Phys.Rg;

    if mc_num == 1
        Ares = zeros(size(D));
        for ii = 1 : length(D)
            Ares(ii) = obj.eos.Ares(point, '0', [], D(ii), []);
        end
        % *****************************************************
        % (f /P) = exp( Ares(T, D) / RT - ln(Z) + (Z - 1) )
        % -----------------------------------------------------
        val = exp(Ares / Rg / T + (Z - 1)) ./ Z;
        % -----------------------------------------------------
    else
        Aresn = zeros(obj.mc_num, length(D));
        for ii = 1 : length(D)
            Aresn(:, ii) = obj.eos.Ares(point, 'n', [], D(ii), X);
        end
        % *************************************************************************
        % (f_i /P /x_i) = exp( partial(Ares(T, D, X)) / partial(n_i) / RT - ln(Z) )
        % -------------------------------------------------------------------------
        val = exp(Aresn / Rg / T) ./ Z; % dimention: r = mc_num, c = phase_num
        % -------------------------------------------------------------------------
    end

end
