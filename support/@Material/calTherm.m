function calTherm(obj, point)
    
    obj.getPointPara(point);

    if Bsc.isblank(point.X)
        obj.X(point);
    end

    if Bsc.isblank(point.phase)
        point.phase = 'FL';
    end
    switch point.phase
    case {'L', 'G', 'SC'}
        if Bsc.isblank(point.D)
            if Bsc.isblank(point.P)
                point.logError('Incomplete information about the thermodynamic state.');
                point.P = -1;
                point.D = -1;
            elseif point.P < 0
                point.logError('Incomplete information about the thermodynamic state.');
                point.P = -1;
                point.D = -1;
            else
                obj.D(point);
            end
        elseif Bsc.isblank(point.P)
            obj.eos.P(point);
        end
    case {'Lsat', 'L, G', 'Gsat', 'G, L'}
        if Bsc.isblank(point.P) || point.P < 0
            if obj.mc_num == 1
                point.P = obj.Psat(point);
                %obj.D(point);
                switch point.phase
                case {'Lsat', 'L, G'}
                    point.D = point.D(1);
                    point.V = point.V(1);
                    point.Z = point.Z(1);
                case {'Gsat', 'G, L'}
                    point.D = point.D(2);
                    point.V = point.V(2);
                    point.Z = point.Z(2);
                end
            else
                obj.VLE(point);
            end
        end
    otherwise
        % if ~isempty(point.D)
        %     if (point.D > obj.fluid.Dcr)
        %         point.phase = 'L';
        %     else
        %         point.phase = 'G';
        %     end
        % else
            obj.phase(point);
            %return
        % end
    end
    
    obj.reduceState(point);

end