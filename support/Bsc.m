classdef Bsc < handle
properties (Constant)

    fluid_data_default = {
        'Carbon dioxide', ...
        [], ...
        'Carbon dioxide', ...
        'R-744', ...
        'CO2.FLD', ...
        'CO2', ...
        'CO2', ...
        44.0098, ...
        194.6860, ...
        216.5920, ...
        304.1282, ...
        7.3773, ...
        10.6249, ...
        0, ...
        0.2239, ...
        0, ...
        [], ...
        'other', ...
        216.5920, ...
        2000, ...
        800, ...
        77, ...
        '77c8bee0', ...
        'IIR', ...
        1013, ...
        1, ...
        [], ...
        'A1', ...
        'CURLTUGMZLYLDI-UHFFFAOYSA-N', ...
        '1S/CO2/c2-1-3' ...
    };

    point_set_code_data_default = {
        'ID',      [],   1, 0, 'A',   {};
        'IDM',     [],   1, 0, 'B',   {}; 
        'FILTER',  [],   1, 0, 'C',   {};
        'FLAG',    [],   1, 0, 'D',   {};
        'MC',      [],   1, 0, 'E:F', {};
        'MCID',    [],   0, 0, {},    {};
        'MCNUM',   [],   0, 0, 'H',   {};
        'XREF',    [],   1, 0, 'I:J', {};
        'X',       [],   0, 0, {},    {};
        'AUTHOR',  [],   1, 0, 'M',   {};
        'YEAR',    [],   1, 0, 'N',   {};
        'PHASE',   [],   1, 0, 'V',   {};
        'T',       [],   1, 0, 'W',   {};
        'PREF',    1e+3, 1, 0, 'X',   {};
        'P',       1e+3, 0, 0, {},    {};
        'PROPREF', [],   0, 0, {},    {};
        'VALREF',  [],   1, 0, 'AB',  {};        
    }

end
methods (Static)
    
    function l = isblank(var)
        if isempty(var)
            l = 1;
        else
            if sum(ismissing(var))
                l = 1;
            else
                l = 0;
            end
        end
    end

    function num =  convertColToNum(col)
    % column in excel to number, e.g., column 'AB' - 28
        base = 64; % base = abs('A') - 1;
        col = upper(col);
        num = 0;
        j = 0;
        for i = length(col) : -1 : 1
            num = num + (abs(col(i)) - base) * 26 ^ j;
            j = j + 1;
        end
    end

    function col =  convertNumToCol(num)
    % number to column in excel, e.g., 28 - column 'AB'

        if num > 26 * 27
            errordlg('convertNumToCol error: num out of range!');
        end

        num_1 = fix(num/26);
        num_0 = mod(num, 26);

        if (num_0 == 0)
            num_0 = 26;
            num_1 = num_1 - 1;
        end

        base = 64; % base = abs('A') - 1;
        col_0 = char(base + num_0);
        if (num_1 == 0)
            col_1 = '';
        else
            col_1 = char(base + num_1);
        end

        col = strcat(col_1, col_0);

    end
    
    function arry = eliminateZero(arry, code)
        if nargin == 1
            code = 0;
        end
        if code ~= 1
            arry(all(arry == 0, 2), :) = [];
        end
        if code ~= 2
            arry(:, all(arry == 0, 1)) = [];
        end
    end
    
    function [col1, col_upper, col_lower] = convertRange(colA, num)
        % To convert the column ranges numbered in alphabat to column numbers.
        % e.g.:
        % [col1, col_upper, col_lower] = convertRange('A:C, E, I:AA', 3)
        % col1 = {[1, 3], [5, 5], [9, 27]}
        % col_upper = 27
        % col_lower = 1

        if isa(colA, 'char')
            colA = strsplit(colA, {',', '，', ' '});
            if nargin == 2
                if numel(colA) ~= num
                    errordlg('Bsc.convertRange: Wrong column info.');
                end
            end
        end
        col1 = cell(size(colA));
        % col_upper = 1;
        % col_lower = 1;
        for i = 1 : size(col1, 1)
            for j = 1 : size(col1, 2)
                if or(isempty(colA{i, j}), ismissing(colA{i, j}))
                    col1{i, j} = [NaN, NaN];
                elseif strcmp(colA{i, j}(1), '#')
                    col1{i, j} = [NaN, NaN];
                else
                    tmp = strsplit(colA{i, j}, {':', '：'});
                    switch length(tmp)
                    case 1
                        col1{i, j} = [Bsc.convertColToNum(tmp{1}), Bsc.convertColToNum(tmp{1})];
                    case 2
                        col1{i, j} = [Bsc.convertColToNum(tmp{1}), Bsc.convertColToNum(tmp{2})];
                    otherwise
                        errordlg(strcat('MaterialCode: invalid range argument. ',colA{i, j}));
                    end
                    % if col1{i, j}(2) > col_upper
                    %     col_upper = col1{i, j}(2);
                    % end
                    % if col1{i, j}(1) < col_lower
                    %     col_lower = col1{i, j}(1);
                    % end
                end
            end
        end
        tmp_mat = cell2mat(col1);
        if isempty(tmp_mat)
            col_lower = NaN;
            col_upper = NaN;
        else
            [col_lower, col_upper] = bounds(tmp_mat, 'all');
        end
    end

    function [Xp, Xm, np, nm] = PartialMolar(X, ii)
        Xp = X; Xm = X;
        Xstep = max(1e-8 * X(ii), 1e-12);
        Xp(ii) = min(X(ii) + Xstep, 1);
        Xm(ii) = max(X(ii) - Xstep, 0);
        np = sum(Xp);
        nm = sum(Xm);
        Xp = Xp / np;
        Xm = Xm / nm;
    end

end
end