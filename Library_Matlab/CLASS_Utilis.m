classdef CLASS_Utilis
    % Class for utility
    properties
    end
    
    methods (Static)
        % UtilisClass can only access static methods
        % Instance of UtilisClass can access all methods
        
%%
        function ImageShifted = ImageBitShift(imInput, bitShift)
            % Needs Image Processing Tool Box
            % Adjust image display range
            % need to install "Image Processing Toolbox"
            % n = 1 means intensity of each pixel X 2
            % n = 2 means intensity of each pixel X 4
            % uint8 is easier to get pure white (255 max)
            % uint16 is recommended
            ImageShifted = imadjust(imInput,[0 1/(2^bitShift)],[0 1]);
        end
        
%%
        function NewHandleTickLabel = ConvertAxesLabel(handleTickLabel, ratioOffset, pmOffset)
            % Input: axes ticklabel
            % ratioOffest: original label x n
            % pmOffset: original label +/- n
            NewHandleTickLabel = cell(numel(handleTickLabel),1);
            for ii = 1:numel(handleTickLabel)
                NewHandleTickLabel{ii} = str2double(handleTickLabel{ii}) * ratioOffset + pmOffset;
            end
        end
        
%%
        function SaveFrameDataToVideo(frameCapture, saveDirectory, FPS)
            % frameCapture: generate from getFrame, format frameCapture = struct('cdata',[], 'colormap',[]);
            % saveDirectory: string fullfile name
            % FPS: int
            vOutput = VideoWriter(saveDirectory,'Uncompressed AVI');
            vOutput.FrameRate = FPS;
            open(vOutput)
            writeVideo(vOutput,frameCapture);
            close(vOutput)
        end
        
%%
function SaveFigureToPDF(saveDirectory,varargin)
            % %Save figure to PDF
            % if isempty(varargin)
            %     h = gcf;
            % else
            %     h = varargin{1};
            % end
            % set(h,'Units','Inches');
            % pos = get(h,'Position');
            % set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
            % set(h, 'Renderer', 'painters');
            % print(h,saveDirectory,'-dpdf')
            % print(h,saveDirectory,'-dpng')
            % savefig(gcf,saveDirectory)
            %Save figure to PDF
            if isempty(varargin)
                h = gcf;
            else
                h = varargin{1};
            end
            set(h,'Units','centimeters');
            pos = get(h,'Position');
            set(h,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[ceil(pos(3)), ceil(pos(4))])
            fprintf('Paper size in cm [%f, %f]\n', ceil(pos(3)), ceil(pos(4)))
            %set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[3.1250, 5.5208])
            set(h, 'Renderer', 'painter');
            print(h,saveDirectory,'-dpdf')
            print(h,saveDirectory,'-dpng')
            savefig(gcf,saveDirectory)
        end
%%
        function Value = MathClamp(x, minVal, maxVal)
            Value = min(max(x,minVal),maxVal);
        end
%%
        function InsertText = InsertFigureText(axisHandle, Xpercentage, Ypercentage, TextHere)
            % plot start from bottom left / imshow start from top left
            % Input:
            % axisHandle: current axis
            % Xpercentage: (double)
            % Ypercentage: (double)
            % TextHere: (String)
            XLim = get(axisHandle,'XLim');
            YLim = get(axisHandle,'YLim');
            
            XRange = max(XLim) - min(XLim);
            YRange = max(YLim) - min(YLim);
            X = XRange*Xpercentage;
            Y = YRange*Ypercentage;
            
            if isnumeric(TextHere)
                TextHere = sprintf('%0.2f',TextHere);
            end
            
            InsertText  = text('Parent',axisHandle,'Position',[min(XLim)+X, min(YLim)+Y],'String',TextHere);
        end
        
%%
        function LineNearestPoint(axisHandle, dataX,dataY,thres,color)
            % Input:
            % axisHandle: current axis
            % dataX: vector
            % dataY: vector same size as data X
            % thres: no connection beyond this thres
            % color: line color
            num_points = length(dataX);
            cur_X = dataX(1);
            cur_Y = dataY(1);
            dataX(1) = []; % clear storage, since it stores in cur_X now
            dataY(1) = []; % clear storage, since it stores in cur_Y now
            for ii = 1:num_points
                
                list_distance = sqrt((dataX-cur_X).^2 +  (dataY-cur_Y).^2);
                [distance_min,idx_min] = min(list_distance);
                if distance_min > thres
                    continue
                end
                
                tar_X = dataX(idx_min);
                tar_Y = dataY(idx_min);
                
                plot(axisHandle,[cur_X,tar_X], [cur_Y,tar_Y], color ,'linewidth',0.5,'HandleVisibility','off')
                drawnow
                cur_X = tar_X;
                cur_Y = tar_Y;
                dataX(idx_min) = [];
                dataY(idx_min) = [];
                
                if isempty(dataX)
                    return
                end
            end
        end
        
%%
        function myColorMap = CustomizedColorMap()
            % Load customised colomap
            if ~isfile('CustomizedColorMap.mat')
                error('Colormap .mat fiel not exist')
            end
            cMap = load('CustomizedColorMap.mat');
            myColorMap = cMap.myColorMap;
        end

        %%
        function dataTDMS = ReadTdmsData(FullFilename)
            dataTDMS = [];
            if ~isfile(FullFilename)
                warning("TDMS file is invalid");
                return
            end

            % addpath 'C:\Users\markz\Desktop\ImageProcess' % For func_file_filter
            % addpath 'G:\Fork\DataProcess_CVCC\Library_Matlab\tdms_package\Version_2p5_Final'
            % addpath 'G:\Fork\DataProcess_CVCC\Library_Matlab\tdms_package\Version_2p5_Final\v2p5\tdmsSubfunctions'
            dataTDMS = TDMS_getStruct(FullFilename);          % read data
        end

%%
        function allFileInfo = FileTypeFilter(file_pattern,varargin)
            % Syntax: allFileInfo = func_file_filter('\**\*.sif')
            % Input:
            % file_pattern: string
            % Varargin:
            % {1} 'filename' / 'filecreate' (default): decide the file date 
            % {2} topfolder specified
            
            % >>>>> specify topfolder / use 'desktop' as default
            if length(varargin) == 1
                    topLevelFolder = varargin{1};
            else
                    % Open uigetdir with start path as starting
                    start_path =  'C:\**\Desktop';
                    topLevelFolder = uigetdir(start_path);
            end
            
            % >>>>> Counter cancel selection
            if topLevelFolder == 0
                    allFileInfo = [];
                    fprintf('Selection cancelled \n')
                    return;
            else
                    fprintf('The top level folder is "%s".\n', topLevelFolder);
            end
            
            % >>>>> Specify the file pattern
            % filePattern = [topLevelFolder,file_pattern]; % e.g. filePattern = [topLevelFolder,'\**\*.abc'];
            filePattern = fullfile(topLevelFolder,file_pattern); % changed 20200723
            allFileInfo = dir(filePattern);
            
            % >>>>> Remove any folders with same names as target files
            isFolder = [allFileInfo.isdir];  % isdir: 1 is folder
            allFileInfo(isFolder) = [];
            
            % Get a cell array of strings (not using as output now)
            listOfFolderNames = unique({allFileInfo.folder});
            numberOfFolders = length(listOfFolderNames);
            fprintf('The total number of folders to look in is %d.\n', numberOfFolders);
            
            % Get a cell array of base filename strings (not using as output now)
            listOfFileNames = {allFileInfo.name};
            totalNumberOfFiles = length(listOfFileNames);
            fprintf('The total number of files in those %d folders is %d.\n', numberOfFolders, totalNumberOfFiles);
            
            % Display all files in those folders.
            totalNumberOfFiles = length(allFileInfo);
            
            if totalNumberOfFiles > 0 
                % Get run ID (format: YYYYMMDD_XX)
	            for ii = 1 : totalNumberOfFiles
                    % >>>>> Get date
                    if any(strcmpi(varargin,'filename'))      
                        file_date_index = strfind(allFileInfo(ii).name,'2019');
                        file_date = allFileInfo(ii).name(file_date_index : file_date_index+7);
                    elseif any(strcmpi(varargin,'filecreate'))   
                        file_date = yyyymmdd(datetime(allFileInfo(ii).date,'InputFormat','dd-MMM-yyyy HH:mm:ss'));  
                    else
                        file_date = yyyymmdd(datetime(allFileInfo(ii).date,'InputFormat','dd-MMM-yyyy HH:mm:ss'));  
                    end
                    
                    % >>>>> Get run number if file name contains 'run_'
                    if contains(allFileInfo(ii).name,'run_','IgnoreCase',true)
                        run_idx = strfind(allFileInfo(ii).name,'run_') + length('run_');
                        run_num = allFileInfo(ii).name(run_idx : run_idx+1);
                        allFileInfo(ii).runID = [num2str(file_date),'_',run_num];
                    end
                    
		            % thisFolder = allFileInfo(k).folder;
		            thisBaseFileName = allFileInfo(ii).name;
		            %fullFileName = fullfile(thisFolder, thisBaseFileName);
		            [~, baseNameNoExt, ~] = fileparts(thisBaseFileName); % syntax: [filepath,name,ext] = fileparts(file)
		            fprintf('%s\n', baseNameNoExt);
	            end
            else
	            fprintf('Folder %s has no files in it.\n', topLevelFolder);
            end
            
            fprintf('\nDone looking in all %d folders \nFound %d files in the %d folders.\n', numberOfFolders, totalNumberOfFiles, numberOfFolders);
        end
%%
% CLASS_Utilis.AxisRelabel(LabelAxis, ImageWidth,PixelResolution,StepInMM,'axial','x',0);
% CLASS_Utilis.AxisRelabel(LabelAxis, ImageHeight,PixelResolution,StepInMM,'radial','y');
        function [pos,label] = AxisRelabel(axis, varargin) 
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
        function cMask = CircularMask(pres,imageRow,imageCol,radiusOffset100,offsetXY)
            % offset100 is added for 100 mm window to fit the raw image
            % offsetXY  is vector [X,Y], "+" = center moves left/up
            if nargin == 4 || isempty(offsetXY)
                centerX = imageCol/2;
                centerY = imageRow/2;
            else
                X = offsetXY(:,1);
                Y = offsetXY(:,2);
                centerX = imageCol/2 - X; % since Theoretical - true = offset
                centerY = imageRow/2 - Y;
            end
            
            [columnsInImage,rowsInImage] = meshgrid(1:imageCol, 1:imageRow);
            radius = ((100+radiusOffset100)/2)/pres;
            cMask = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= radius.^2; % compute a map of radius value
        end

        %% Plot pressure from T data all
        function PlotPressureFromTable(InputData, ii, ratio, offset)
            dir_tdms = InputData{ii,'File_Pressure'};
            dir_tdms = dir_tdms{1};
            data = CLASS_Utilis.ReadTdmsData(dir_tdms);
            f = figure; plot(data.Untitled.Pressure__bar_.data*ratio + offset);
            f.Name = sprintf('Run number %2d', ii);
            set(f,'color','white')
            grid minor
        end
    end
end

