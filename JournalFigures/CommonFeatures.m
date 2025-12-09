classdef CommonFeatures
    %COMMONFEATURES Summary of this class goes here
    %   Detailed explanation goes here
    
    properties

    end
    
    methods (Static)
        function text_flame_behaviour = Label_FlameBehaviour5Rows()
            text_flame_behaviour = {'Pilot ignition','Jets intersection','Main ignition','First HRR peak','Quasi-steady'};
        end
        function text_flame_behaviour = Color_Gradient()
            text_flame_behaviour = {'#000000','#001C7D','#0039FF','#7291FF'};
        end
        function time_eoi_pilot = TimeEOI_Pilot()
            time_eoi_pilot = (50-11)/37.5;
        end

        function time_eoi_main = TimeEOI_Main()
            time_eoi_main = (511-11)/37.5;
        end
        function [text_x_percent, text_y_percent] = TextPosition3Columns()
            text_x_percent = 0.025;
            text_y_percent = 0.86;
        end
    end
end

