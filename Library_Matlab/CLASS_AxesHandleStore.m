classdef CLASS_AxesHandleStore

    properties
        MarginLeft = 50             % left white margin
        MarginBottom = 50        % bottom white margin
        MarginRight = 10
        MarginTop = 10
        
        GapRow = -1                 % gap between rows
        GapColumn = 3            % gap between columns
        RowNumber = 5             % how many rows
        ColumnNumber = 5        % how many columns
        AxesWidth = 100            % width of the axes
        AxesHeight = 100            % height if the axes

        FontSize = 9;
    end
    
    properties (Access = private)
        AxesHandleMatrix           % store axes handles
        FigureHandle
    end
    
    methods
        function obj = ConstuctAxes(obj)
            f = figure;
            obj.FigureHandle = f; 
            % construct axes
            axesMatrix = gobjects(obj.RowNumber, obj.ColumnNumber);
            for ii = 1:1:obj.ColumnNumber
                    for jj = 1:1:obj.RowNumber
                            axPosition1 = obj.MarginLeft + (obj.GapColumn + obj.AxesWidth) * (ii - 1);
                            axPosition2 = obj.MarginBottom + (obj.GapRow + obj.AxesHeight) * (obj.RowNumber - jj);
                            axPosition3 = obj.AxesWidth;
                            axPosition4 = obj.AxesHeight;
                            axPosition = [axPosition1, axPosition2, axPosition3, axPosition4];
                            currentAxis = axes('Unit','Pixel','Position',axPosition); box(currentAxis,"on"); set(currentAxis,'FontSize', obj.FontSize) 
                            axesMatrix(jj,ii) = currentAxis;
                    end
            end
            obj.AxesHandleMatrix = axesMatrix;
            
            % resize figure handle based on axes
            figureWidth =  obj.MarginLeft + (obj.AxesWidth + obj.GapColumn ) * obj.ColumnNumber - obj.GapColumn + obj.MarginRight;
            figureHeight = obj.MarginBottom + (obj.AxesHeight + obj.GapRow ) * obj.RowNumber - obj.GapRow + obj.MarginTop;
            set(f,'Unit','Pixel','Position',[5, 50, figureWidth, figureHeight])
            set(f,'color','white');
            
            if figureWidth >  3.4811e+03 % Max figure width
                widthAllAxes = 3.4811e+03 - obj.MarginLeft + obj.GapColumn - obj.MarginRight;
                amountAxes = floor(widthAllAxes / (obj.AxesWidth + obj.GapColumn));
                warning('Figure width exceeds maximum, only %1.0d/%1.0d axes are shown', amountAxes, obj.ColumnNumber);
            end
        end
        %

        %%
        function obj = ConstuctAxesUnevenGap(obj)
            if length(obj.GapRow) ~= obj.RowNumber
                warning('Gap array row is different from Axes size')
                return
            end

            if length(obj.GapColumn) ~= obj.ColumnNumber
                warning('Gap array collumn is different from Axes size')
                return
            end
            GapColumnArray = obj.GapColumn;
            GapRowArray = obj.GapRow;

            f = figure;
            obj.FigureHandle = f;
            % construct axes
            axesMatrix = gobjects(obj.RowNumber, obj.ColumnNumber);
            for ii = 1:1:obj.ColumnNumber
                    for jj = 1:1:obj.RowNumber
                            axPosition1 = obj.MarginLeft + sum(GapColumnArray(1 : (ii - 1))) + obj.AxesWidth * (ii - 1);
                            axPosition2 = obj.MarginBottom + sum(GapRowArray(1 : (obj.RowNumber - jj))) + obj.AxesHeight * (obj.RowNumber - jj);
                            axPosition3 = obj.AxesWidth;
                            axPosition4 = obj.AxesHeight;
                            axPosition = [axPosition1, axPosition2, axPosition3, axPosition4];
                            currentAxis = axes('Unit','Pixel','Position',axPosition); box(currentAxis,"on"); set(currentAxis,'FontSize', obj.FontSize) 
                            axesMatrix(jj,ii) = currentAxis;
                    end
            end
            obj.AxesHandleMatrix = axesMatrix;
            
            % resize figure handle based on axes
            figureWidth =  obj.MarginLeft + obj.AxesWidth  * obj.ColumnNumber  + obj.MarginRight + sum(GapColumnArray(1:end-1));
            figureHeight = obj.MarginBottom + obj.AxesHeight * obj.RowNumber  + obj.MarginTop + sum(GapRowArray(1:end-1));
            set(f,'Unit','Pixel','Position',[5, 50, figureWidth, figureHeight])
            set(f,'color','white');
            
            if figureWidth >  3.4811e+03
                widthAllAxes = 3.4811e+03 - obj.MarginLeft + obj.GapColumn - obj.MarginRight;
                amountAxes = floor(widthAllAxes / (obj.AxesWidth + obj.GapColumn));
                warning('Figure width exceeds maximum, only %1.0d/%1.0d axes are shown', amountAxes, obj.ColumnNumber);
            end
        end
        %%
        function axesHandleMatrix = GetAxesHandleMatrix(obj)
             % access the axes handles
             axesHandleMatrix = obj.AxesHandleMatrix;
        end

        function figureHandle = GetFigureHandle(obj)
             % access the axes handles
             figureHandle = obj.FigureHandle;
        end
    end
end

