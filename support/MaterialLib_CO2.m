classdef MaterialLib_CO2 < MaterialLib

    methods
        function obj = MaterialLib_CO2()
            obj.lib_size = 1;
            obj.pure     = Material_CO2;
        end
    end

end