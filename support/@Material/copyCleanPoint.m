function pt = copyCleanPoint(obj, point, T, D, P, X)
    % T, D, P, X < 0 : clean
    % T, D, P, X = [] : reserve
    % T, D, P, X > 0 : change to new value

    pt = copy(point);

    if ~isempty(T)
        if T < 0
            pt.T = [];
        else
            pt.T = T;
        end
    end
    if ~isempty(D)
        if D < 0
            pt.D = [];
        else
            pt.D = D;
        end
    end
    if ~isempty(P)
        if P < 0
            pt.P = [];
        else
            pt.P = P;
        end
    end
    if ~isempty(X)
        if X(1) < 0
            pt.X = obj.X(pt);
        else
            pt.X = X;
        end
    end

    obj.getPointPara(pt);
    obj.calTherm(pt);
end