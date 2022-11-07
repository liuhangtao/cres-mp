classdef Material < matlab.mixin.Copyable

    properties
        mc_num
        id
        name
        fluid
        eos
        eos0
        vis
        vis0
        tcx
        tcx0
        setup  % reserved
    end
    
    methods

        varargout = calProp(obj, prop_req, point);

        calTherm(obj, point);
        
        reduceState(obj, point);

        function setPointModel(obj, point)
            obj.eos.setPointModel(point);
            obj.vis.setPointModel(point);
            obj.tcx.setPointModel(point);
        end 
        
        function getPointPara(obj, point)
            if ~or(isempty(point.T), isempty(point.X))
                obj.eos.getParaT(point);
                obj.eos0.getParaT(point);
                obj.vis.getParaT(point);
                obj.vis0.getParaT(point);
                obj.tcx.getParaT(point);
                obj.tcx0.getParaT(point);
            end
            if ~isempty(point.X)
                obj.fluid.MW(point);
            end
        end

        val = D(obj, varargin);
            % molar density [mol/m3]

        val = X(obj, varargin);
            
        val = Z(obj, point);
            % compression factor [dimensionless]
            
        val = cvtMolar2MassD(obj, varargin);
            % convert molar volume [mol/m3] to mass density [kg/m3]
            
        val = cvtMass2MolarD(obj, varargin);
            % convert mass density [kg/m3] to molar volume [mol/m3]
            
        val = Gres(obj, point);
            % residual Gibbs function [J/mol]

        val = Sres(obj, point);
            % residual entropy [J/mol/K]

        val = Ures(obj, point);
            % residual internal energy [J/mol]

        val = Hres(obj, point);
            % residual enthalpy [J/mol]

        val = CVres(obj, point);
            % residual specific heat capacity at constant volume [J/mol/K]

        val = CPres(obj, point);
            % specific heat capacity at constant pressure [J/mol/K]
            
        val = CV(obj, point)
            % specific heat capacity at constant volume [J/mol/K]

        val = CP(obj, point)
            % specific heat capacity at constant pressure [J/mol/K]

        val = Cratio(obj, point)
            % specific heat capacity at constant pressure [J/mol/K]    

        val = W(obj, point);
            % sound speed [m/s]

        val = FC(obj, varargin);
            % fugacity coefficient [dimensionless]

        [D1, D2, P1, P2] = spinodal(obj, varargin);
            % To calculate the extremum of pressure on the isotherm line. 

        [phase_flag, D] = phase(obj, varargin);
            % To determine the phase
            
        p_sat = Psat(obj, point);
            % saturated pressure at a specific temperature [Pa]

        val = Trp_RES_CO2_pre(obj, varargin);
            % preparation for the crossover RES model of CO2
            
        pt = copyCleanPoint(obj, point, T, D, P, X);
    end
    
    methods (Access = protected)
        function newobj = copyElement(obj)
            newobj = copyElement@matlab.mixin.Copyable(obj);
            newobj.eos = copy(obj.eos);
            newobj.eos0 = copy(obj.eos0);
            newobj.vis = copy(obj.vis);
            newobj.vis0 = copy(obj.vis0);
            newobj.tcx = copy(obj.tcx);
            newobj.tcx0 = copy(obj.tcx0);
        end
    end

end