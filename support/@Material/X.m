function val = X(obj, varargin)
    % format:
    %  - X(point)
    %  - X(point, T, D, P)
    %  - X(       T, D, P)

    flag_point_edit = 0;
    if isa(varargin{1}, 'Point')
        point = varargin{1};
        switch nargin
        case 2 % X(point)
            flag_point_edit = 1;
            T = point.T;
            P = point.P;
            D = point.D;
        case 5 % D(point, T, D, P)
            T = varargin{2};
            D = varargin{3};
            P = varargin{4};
            if isempty(T)
                T = point.T;
            end
            if isempty(P)
                P = point.P;
            end
            if isempty(D)
                D = point.D;
            end
        end
    else % D(T, D, P, X)
        T = varargin{1};
        D = varargin{2};
        P = varargin{3};
        X = varargin{4};
        point = Point('T', T, 'D', D, 'P', P, 'X', X);
        obj.getPointPara(point);
    end

    %if Bsc.isblank(point.X)
    %    x3 = 0.5;
    %else
    %    x3 = 1 - point.X(1);
    %end
    %P3 = obj.eos.P(point, '0', [], [], [1-x3,x3]);

    eps0 = 2e-8;
    itt_count = 0;
    xi = 1 - point.X(1);
    P_cal = obj.eos.P(point, '0', [], [], [1-xi,xi]);
    while abs(P_cal - P) > P * eps0
        P_n = (1 - xi) * obj.eos.P(point, 'n', [], [], [1-xi,xi]);
        xi = xi - (P_cal - P) / P_n(2);
        if xi < 0
            xi = eps0;
        end
        if xi > 1
            xi = 1 - eps0;
        end
        
        itt_count = itt_count + 1;
        if itt_count > 100
            xi = -101;
            point.logError('point.X: not converge.');
            break;
        end
        
        P_cal = obj.eos.P(point, '0', [], [], [1-xi,xi]);
        %disp([Di, P_cal]);
    end
    
%     while abs(P3 - P) > eps0 * P %[Pa]
%         if (P1 - P)*(P3 - P) > 0
%             x1 = x3;
%             P1 = P3;
%         else
%             x2 = x3;
%         end
%         x3 = (x1 + x2) / 2;
%         %disp( abs(P3 - P)/P );
%         P3 = obj.eos.P(point, '0', [], [], [1-x3,x3]);
%         count = count + 1;
%         if count > 200
%             point.logError('Fail to solve molar fraction.');
%             break;
%         end
%     end
%     if P3 < 0
%         pause;
%     end
    val = [1-xi,xi];
    
    if flag_point_edit
        point.X = val;
    end
end