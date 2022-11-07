function [phase_flag, rho] = phase(obj, varargin)
% To determine the phase
% format:
% - [phase, D] = phase(point)
% - [phase, D] = phase(T, P)
% - [phase, D] = phase(T, P, X)

    flag_point_edit = 0;
    switch class(varargin{1})
    case 'Point'
        point = varargin{1};
        flag_point_edit = 1;
    case 'double'
        T = varargin{1};
        P = varargin{2};
        switch nargin
        case 3
            X = 1;
        case 4
            X = varargin{3};
        end
        point = Point('T',T, 'P',P, 'X',X);
        obj.getPointPara(point);
    end

    if point.P == 0
        rho = 0;
        phase_flag = 'G';
    else
        if point.mc_num == 1
            if point.T < obj.fluid.Tcr
                if ~Bsc.isblank(point.D)
                    if (point.D > obj.fluid.Dcr)
                        point.phase = 'L';
                    else
                        point.phase = 'G';
                    end
                    point.P = obj.eos.P(point);
                    point.V = 1 / point.D;
                    obj.Z(point);
                    return
                    % debug here.
                end
                if point.P < obj.fluid.Pcr
                    rho = obj.D(point, 'L&G');
                    if length(rho) ~= 2
                        point.logError('Material.phase: Invalid density.');
                    end
                    fc = obj.FC(point, rho);
                    if length(fc) < 2
                        %errordlg('Material.phase: length(fc)<2.')
                        %pause
                        disp('Material.phase: length(fc)<2.');
                        point.logError('Material.phase: length(fc)<2.');
                    end
                    if fc(1) < fc(2)
                        phase_flag = 'L';
                        rho = rho(1);
                    elseif fc(1) >fc(2)
                        phase_flag = 'G';
                        rho = rho(2);
                    else
                        if rho(1) > obj.fluid.Dcr
                            phase_flag = 'L';
                        else
                            phase_flag = 'G';
                        end
                        rho = rho(1);
                        point.logError('phase: 2 identical real roots for density.');
                        % phase_flag = 'sat';
                        % disp('Material.phase(): debug here! (11042308)');
                        % pause;
                    end
                else
                    phase_flag = 'L';
                    rho = obj.D(point, phase_flag);
                end
            else
                if ~Bsc.isblank(point.D)
                    point.P = obj.eos.P(point);
                    if point.P < obj.fluid.Pcr
                        phase_flag = 'G';
                    else
                        phase_flag = 'SC';
                    end
                    point.phase = phase_flag;
                    point.V = 1 / point.D;
                    obj.Z(point);
                    return
                    % debug here.
                end
                if point.P < obj.fluid.Pcr
                    phase_flag = 'G';
                    rho = obj.D(point, phase_flag);
                else
                    phase_flag = 'SC';
                    rho = obj.D(point, phase_flag);
                end
            end
        else % for mixture
            if ~Bsc.isblank(point.D)
                if (point.D > obj.fluid.Dcr)
                    point.phase = 'L';
                else
                    point.phase = 'G';
                end
                point.P = obj.eos.P(point);
                point.V = 1 / point.D;
                obj.Z(point);
                return
                % debug here.
            end
            rho = obj.D(point);
            rho = rho(rho > 0);
            switch length(rho)
            case 0
                point.logError(strcat('phase: fail to solve density at T =', {32}, point.T, {32}, 'K, P =', {32}, point.P/1e6, {32}, 'MPa.'));
                rho = -1;
            case 1
                if point.T > min(obj.fluid.Tcr)
                    if point.P < max(obj.fluid.Pcr)
                        phase_flag = 'G';
                    else
                        phase_flag = 'SC'; % not exactly.
                    end
                end                
            case 2
                fc = obj.FC(point, rho);
                %FC_liq = obj.FC(point, V_liq);
                %FC_vap = obj.FC(point, V_vap);
                %if prod(FC_liq < FCG_vap)
                if prod(fc(:, 1) < fc(:, 2))
                    phase_flag = 'L';
                    rho = rho(1);
                elseif prod(fc(:, 1) > fc(:, 2))
                    phase_flag = 'G';
                    rho = rho(2);
                else
                    phase_flag = 'sat'; % debug here!
                    disp('Material.phase(): debug here! (11042309)');
                    %pause;
                    rho = rho(1);
                    point.logError('Material.phase: Fail-11042309.');
                end
            end
        end
    end

    if flag_point_edit
        point.phase = phase_flag;
        point.D = rho;
        point.V = 1 / rho;
        obj.Z(point);
        if rho == -1
            point.status_eos = -1;
        end
    end

end