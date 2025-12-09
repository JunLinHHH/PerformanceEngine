% Generate video matrix for all runs at the same test conditions 

% Input:  DataTable (T_DataAll), DataSheet(.xlsx file)
% Output: Video .avi
% Requires load('MethaneAutoignitionDataSet.mat');

% If .mat exist
clearvars; close all; clc; 
Configs = Configuration();
Loaded = load(Configs.DataObjReadDir);

obj = Loaded.obj;
T_DataAll = obj.DataMatrix;
T_Condition_ROI = obj.ConditionMatrix;
%% Video configuration
% Inputs (Change accordingly)
ConditionSet = [1];
SaveConfig.FileOutDir = Configs.VideoOutput;
SaveConfig.DoSave = 1;

VideoConfig.FrameRange = 1:600;  %

VideoConfig.StepInMM = 10;
VideoConfig.PixelDownstream = 750;

VideoConfig.RatioResize = 0.3;
VideoConfig.RowNumber = 6;
VideoConfig.ColumnNumber = 6;
%% Build storage sturcture
for ii = ConditionSet 
    ConditionTag = Class_TableManagement.ConditionTagGenerator(T_Condition_ROI, obj, ii);

    runsForThisCondition = T_Condition_ROI.RunNumbers{ii};
    ConditionSet_NonDiscarded = T_DataAll(runsForThisCondition, :);

    SaveConfig.SaveName = fullfile(Configs.VideoOutput, ConditionTag);
    fprintf('CONDITION: #%02d\n',ii);
    fprintf('Processing %02d: %s\n',ii,ConditionTag)

    [FrameCapture, DoSave, SaveName] = FUNC_MatrixVideo(ConditionSet_NonDiscarded,VideoConfig, SaveConfig);

    if SaveConfig.DoSave == 1
        CLASS_Utilis.SaveFrameDataToVideo(FrameCapture,SaveName,15)
        fprintf('Video Saved\n')
    end

    close all
end

fprintf('Script >> %s << is finished\n', mfilename());

%%
function [FrameCapture, SaveName] = FUNC_MatrixVideo(inputTable,VideoConfig, SaveConfig)
    if isempty(inputTable)
        fprintf("No runs for this condition\n")
        return
    end
    % Inputs unpacking 
    % Video configurations
    FrameRange = VideoConfig.FrameRange;  %
    StepInMM = VideoConfig.StepInMM;
    PixelDownstream = VideoConfig.PixelDownstream;
    RatioResize = VideoConfig.RatioResize;
    RowNumber = VideoConfig.RowNumber;
    ColumnNumber = VideoConfig.ColumnNumber;    
    % Save configurations
    SaveName = SaveConfig.SaveName;
    %%
    % Store images to cells
    ImageTifCellTotal = cell(height(inputTable),1);

    fprintf('Reading Image ... \nCurrent File:      ')
    for ii = 1:height(inputTable)
        fprintf('\b\b\b\b\b%04d',ii)

        % Read from table
        thisTable = inputTable(ii,:);
        ImageDir = thisTable.File_Image;
        PixelResolution = thisTable.Pixel_res;
        NozzleX = thisTable.Nozzle_X;
        FPms =thisTable.FPS/1000;

        try
            ImageTif = CLASS_MultiPageTif(ImageDir{1});
            ImageTif = ImageTif.StoreTifToCell(FrameRange);
            ImageTifCellSingleRun = ImageTif.GetTifStoreCell();         % Get tif cell for each .tif file 
        catch
            ImageTifCellSingleRun = {};
            fprintf('\n')
            warning('Tif file loading error: %s',ImageDir{1})
            if isempty(ImageDir)
                warning('Tif file not Found')
            end
        end

        % Crop nozzle
        if NozzleX~= 0 && ~isempty(ImageTifCellSingleRun)
            for jj = 1:numel(ImageTifCellSingleRun)
                temp = ImageTifCellSingleRun{jj};
                temp(:,1:NozzleX) = [];
                temp(:,PixelDownstream+1:end) = [];
                ImageTifCellSingleRun{jj} = temp;
            end
        end
    
        ImageTifCellTotal{ii} = ImageTifCellSingleRun; % store in the total cell (size = file number)
        fprintf('\n')
    end
    
    % Resize (Optional)
    ImageHeight = [];
    ImageWidth= [];
    if RatioResize ~= 1
        [ImageHeight,ImageWidth] = size(temp);
        ImageWidth = ImageWidth - NozzleX;      % imageTIF.ImageWidth;
        
        ImageSample = ones(ImageHeight,ImageWidth);
        ImageSampleResize =  imresize(ImageSample,RatioResize);
        ImageHeightResize = size(ImageSampleResize,1);
        ImageWidthResize = size(ImageSampleResize,2);
    else
        ImageHeightResize = ImageHeight;
        ImageWidthResize = ImageWidth;
    end
    
    %% Build axes matrix
    AxesHandles = CLASS_AxesHandleStore;
    AxesHandles.RowNumber = RowNumber;     % make this flexibe 
    AxesHandles.ColumnNumber = ColumnNumber;  % make this flexibe 
    AxesHandles.MarginRight = 10;
    AxesHandles.MarginTop = 40;
    AxesHandles.GapRow = -1;                 % gap between rows
    AxesHandles.GapColumn = 3;            % gap between columns
    AxesHandles.AxesWidth = ImageWidthResize;
    AxesHandles.AxesHeight = ImageHeightResize;
    
    AxesHandles = AxesHandles.ConstuctAxes();
    AxesMatrix = AxesHandles.GetAxesHandleMatrix();
    
    %% Add images (each axes for one file updating)
    if numel(ImageTifCellTotal) < RowNumber*ColumnNumber
        EmptyCellAdd = RowNumber*ColumnNumber - numel(ImageTifCellTotal);
        while EmptyCellAdd > 0
            ImageTifCellTotal{numel(ImageTifCellTotal)+1} = [];
            EmptyCellAdd = EmptyCellAdd - 1;
        end
    elseif numel(ImageTifCellTotal) > RowNumber*ColumnNumber
        warning('Axes are not enough for all images')
        return
    end

    ImageTifCellTotal = reshape(ImageTifCellTotal,AxesHandles.RowNumber,AxesHandles.ColumnNumber);   % depends on axes [m x n]
    FrameCapture = struct('cdata',[], 'colormap',[]);
    LabelAxis = AxesMatrix(end,1);
    
    ReferenceAxis = AxesMatrix(1,1);
    set(ReferenceAxis,'units','pixels');
    Position = get(ReferenceAxis,'position');
    Position(2) = Position(2) + Position(4) + 10;
    Position(4) = 20;
    
    HandleTime = createtextbox(gcf, Position);
    
    fprintf('Image processing ... \nCurrent Frame:     ')
    for kk = 1:length(FrameRange)
        fprintf('\b\b\b\b%04d',kk)
        
        HandleTime.String = sprintf('Time after trigger [ms]: %0.2f', FrameRange(kk)/FPms); % FPS and Injection delay
        
        for ii = 1:size(AxesMatrix,2)
            for jj = 1:size(AxesMatrix,1)
            
                this_image_cell = ImageTifCellTotal{jj, ii};
                this_axis = AxesMatrix(jj,ii);
    
                if ~isempty(this_image_cell)
                    this_image= ImageTifCellTotal{jj, ii}{kk};
                    imshow(this_image,[],'Parent',this_axis);
                    this_axis.Visible = 'on';
                    if ~(ii == 1 && jj == size(AxesMatrix,1))
                        % Hide axis label if not the left-bottom axis
                        FUNCTION_AxisRelabel(this_axis, ImageWidth,PixelResolution,StepInMM,'axial','x',0);
                        FUNCTION_AxisRelabel(this_axis, ImageHeight,PixelResolution,StepInMM,'radial','y');
                        this_axis.XTickLabel = [];
                        this_axis.YTickLabel = [];
                    else
                        FUNCTION_AxisRelabel(LabelAxis, ImageWidth,PixelResolution,StepInMM,'axial','x',0);
                        FUNCTION_AxisRelabel(LabelAxis, ImageHeight,PixelResolution,StepInMM,'radial','y');
    
                        LabelAxis.XLabel.String = 'Axial [mm]';
                        LabelAxis.YLabel.String = 'Radial [mm]';
                    end
                    text_x_percent = 0.02;
                    text_y_percent = 0.86;
                    this_text = CLASS_Utilis.InsertFigureText(this_axis, text_x_percent ,text_y_percent, sprintf('%04d', FrameRange(kk)));
                    set(this_text,'Color',[1,1,1])
                    set(this_text,'BackgroundColor','#373A40')
                    set(this_text,'FontName','Helvetica')
                    set(this_text,'Fontsize',9)
                else
                    % Hide axis if no image
                    axis off
                    this_axis.Visible = 'off';
                end  
            end
        end
        
        FrameCapture(kk) = getframe(gcf);
    end
        fprintf('\n')
        
    %%
end


%%
function [pos,label] = FUNCTION_AxisRelabel(axis, varargin) 
    % Input:
    % 1 size (double: image row or col)
    % 2 pres (double: pixel/mm)
    % 3 step (double: how many mm between two ticks)
    % 4 type (string: 'axial' distance or 'radial')
    % 5 axis (string: scaling to which axis, 'x' or 'y')
    % 6 start (double: default 5) offset of first value shown on scale
    
    % Size is after trim
    size = varargin{1};
    pres = varargin{2};
    step = varargin{3};
    type = varargin{4};
    axisXY = varargin{5};
    
    if nargin > 6
        start = varargin{6};
    else
        start = 0;
    end
    
    if strcmp(type,'axial')
               
                pos = 1:round(step/pres):size;
                axial_base = (1:length((pos))) - 1;   % for getting 0 as begining 
                label = axial_base * step + start;
                        
                    if strcmp(axisXY,'x')
                        set(axis, 'XTick', pos);
                        set(axis, 'XTickLabel', label);
                    elseif strcmp(axisXY,'y')
                        set(axis, 'YTick', pos);
                        set(axis, 'YTickLabel', label);
                    end    
    
    elseif strcmp(type,'radial')
            
                center = size/2;
                base = 0; % base is center postion (in pixel)
                
                % step/pres gives pixel count for one step
                % loop until out of range
                while (center + base * (step/pres) ) < max(size)        
                    base = base + 1;       
                end
                
                pos = ( -(base-1) : 1 :(base-1)) * (step/pres) + center;
                label = ((-(base-1):1:(base-1))*step);  % remove abs if negative is required
                
                    if strcmp(axisXY,'x')
                        set(axis, 'XTick', pos);
                        set(axis, 'XTickLabel', label);
                    elseif strcmp(axisXY,'y')
                        set(axis, 'YTick', pos);
                        set(axis, 'YTickLabel', label);
                    end    
    
    end
end

%%
function HandleTime = createtextbox(Figure,Position)
    %CREATETEXTBOX(figure1)
    %  FIGURE1:  annotation figure
    
    %  Auto-generated by MATLAB on 13-Jun-2021 01:01:51
    
    % Create textbox
    HandleTime = annotation(Figure,'textbox',...
        [0.3 0.9 0.300 0.045],...
        'HorizontalAlignment','center',...
        'VerticalAlignment','middle',...
        'String','Time aSOI [ms]:',...
        'FitBoxToText','off');
    set(HandleTime,'units','pixels',...
        'Position',Position)
end


