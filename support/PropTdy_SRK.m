classdef PropTdy_SRK < PropTdy

    properties
        b
        a
        at
        att
        c
        ct
        ctt
    end
    
    methods
        function obj = PropTdy_SRK(eos_SRK, point)
            switch nargin
            case 0
            case 2
                T = point.T;
                X = point.X;

                obj.b                    = eos_SRK.b(   X   );
                [obj.att, obj.at, obj.a] = eos_SRK.a(T, X, 2);
                [obj.ctt, obj.ct, obj.c] = eos_SRK.c(T, X, 2);
            end
        end

    end

end