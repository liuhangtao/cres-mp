classdef Material_CO2 < Material

    methods

        function obj = Material_CO2()
            data = Para_CO2;

            obj.mc_num = 1;
            obj.id     = 1;
            obj.name   = 'CO2';
            obj.fluid  = Fluid(obj.id, data.fluid{1});

            obj.eos  = ModelEos_MP_CO2(obj.fluid, data.eos{2}, data.eos{1});
            obj.eos0 = ModelEos0_Wilhoit(obj.fluid, data.eos0{2});
            obj.vis  = ModelVis_RES_CO2(obj.fluid, obj.eos, data.vis{2}, data.vis{1});
            obj.vis0 = ModelTrp0_CE(obj.fluid, obj.eos0, data.vis0{2}, data.vis0{1});
            obj.tcx  = ModelTcx_RES_CO2(obj.fluid, obj.eos, data.tcx{2}, data.tcx{1});
            obj.tcx0 = ModelTrp0_CE(obj.fluid, obj.eos0, data.tcx0{2}, data.tcx0{1});

        end

    end

end