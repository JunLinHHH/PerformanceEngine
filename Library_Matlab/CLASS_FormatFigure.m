classdef CLASS_FormatFigure
    %CLASS_JOURNALFORMATFIGURE Summary of this class goes here
    %   Detailed explanation goes here
% objFig = CLASS_FormatFigure;
%         objFig.PositionAxis = [60 60 270 270];
%         objFig.PositionFigure = [50 50 340 340];
% objFig.SingleYAxis(HandleFigure);
    properties
        PositionAxis = [50 50 280 140];
        PositionFigure = [50 50 340 200];
    end
    
    methods
        function this = SingleYAxis(this, HandleFigure)
            set(HandleFigure,'color','white');
            idx = CLASS_FormatFigure.GetAxesHandle(HandleFigure);
            ChildrenObject = HandleFigure.Children;
            HandleAxis = ChildrenObject(idx);
            set(HandleAxis,'box','on');
            grid(HandleAxis,'on');
            set(HandleAxis,'fontsize',9,'fontname','helvetica');
            set(HandleAxis,'unit','pixel','position',this.PositionAxis)
            set(HandleFigure,'unit','pixel','position',this.PositionFigure)
        end
    end
    methods (Static)
        function idx = GetAxesHandle(HandleFigure)
            % Check which children is the 'axes' object
            ChildrenObject = HandleFigure.Children;
            idx = 0;
            for ii = 1:length(ChildrenObject)
                thisObject = ChildrenObject(ii);
                if strcmp(get(thisObject,'type'),'axes')
                    idx = ii;
                end
            end
        end
        
    end
end

