classdef Model < matlab.mixin.Copyable

    properties
        mc_num
        id
        setup   % reserved

    end

    methods
        function obj = Model(varargin)
            % compound:
            %  - Model()
            %  - Model(fluid)
            % mixture:
            %  - Model(fluid)
            %  - Model([component_mdl])

            if nargin == 0
                obj.mc_num = 1;
                obj.id = 1;
                return
            end
            
            switch class(varargin{1})
            case 'Fluid'
                fluid = varargin{1};
                obj.mc_num = fluid.mc_num;
                obj.id     = fluid.id;

                % switch nargin
                % case 1
                %     code_pure = 0;
                % case 2
                %     code_pure = varargin{2};
                % end

                % switch code_pure
                % ...

            case 'Model'
                component_mdl = varargin{1};
                obj.mc_num = length(component_mdl);
                obj.id     = [component_mdl.id]';

                % obj.name   = {component_mdl.name}';
                % obj.setup  = cat(1, component_mdl.setup);
                % ...

            otherwise
                errordlg('Model: Wrong input arguments.')

            end
                
        end
        
        function getParaT(obj, point)
            
        end

%         function [Tcr, Pcr, Dcr] = reduceState(obj, vararg)
%             if isa(vararg, 'Point')
%                 Tcr = vararg.X * obj.fluid.Tcr;
%                 Pcr = vararg.X * obj.fluid.Pcr;
%                 Dcr = vararg.X * obj.fluid.Dcr;
%                 %Zcr = Pcr / Dcr / Phys.Rg / Tcr;
%                 vararg.Tcr = obj.fluid.Tcr;
%                 vararg.Pcr = obj.fluid.Pcr;
%                 vararg.Dcr = obj.fluid.Dcr;
%                 vararg.Tr  = vararg.T / Tcr;
%                 vararg.Pr  = vararg.P / Pcr;
%                 vararg.Dr  = vararg.D / Dcr;
%             else
%                 Tcr = obj.fluid.Tcr;
%                 Dcr = obj.fluid.Dcr;
%                 Pcr = obj.fluid.Pcr;
%             end
%         end
        
    end
    
    methods (Static)
        function setPointModel(point)
            % e.g.: [point.tdy] = deal(PropTdy_SRK());
        end
        
    end
    
end