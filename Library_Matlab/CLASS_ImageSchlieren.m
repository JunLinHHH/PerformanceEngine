classdef CLASS_ImageSchlieren
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ImageCellArray = []        % Raw image (stores in cell array)
        FactorBS = 6            % Bitshift factor
        SizeFilter = 5;
        FactorOpen_1 = 50       % Factor for "bwareaopen" 1
        FactorOpen_2 = 100      % Factor for "bwareaopen" 2
        FactorDE = 5            % Factor for erode and dilate
        FilterKernal            % Standard deviation filter window

        ImCompacted
    end
    
    properties (Access = private)
        JetBoundary = []
        JetBWImage = []
    end
    
    methods
        %% constructor
        function this = CLASS_ImageSchlieren()
            this.FilterKernal = ones(this.SizeFilter, this.SizeFilter);

            imCompacted.imRaw = [];
            imCompacted.imStd = [];
            imCompacted.imStdFiltered = [];
            imCompacted.imbinarized = [];
            imCompacted.imbinarized_Opened_1 = [];
            imCompacted.imbinarized_Opened_2 = [];
            imCompacted.imEroded = [];
            imCompacted.imDilatedComplementedFilled = [];
            
            this.ImCompacted = imCompacted;
        end

        %% Getter
        function [jetBoundary, jetBW] = GetBoundary(this)
            jetBW = this.JetBWImage;
            jetBoundary = this.JetBoundary;
        end
        
        %% Setter
        function this = SetBoundary(this, jetBoundary)
            this.JetBoundary = jetBoundary; 
        end

        %%
        function this = JetBoundaryProcess(this, imageCellArray)
            this.ImageCellArray = imageCellArray;
            if ~iscell(imageCellArray)
                error("Input data is not valid (type must be 'cell')")
            elseif numel(imageCellArray) ~= 3
                error("Input cell must include 3 concecutive frames")
            end

            % Get image size
            imTemp = this.ImageCellArray{1};
            imWidth = size(imTemp,2);
            imHeight =  size(imTemp,1);
            imMatrix = zeros(imHeight, imWidth, 3); % Image index i-1, i, i+1

            % Image index i-1, i, i+1
            for ii = 1:numel(this.ImageCellArray)
                    imMatrix(:,:,ii) = CLASS_Utilis.ImageBitShift(this.ImageCellArray{ii},this.FactorBS);
            end

            % Processing jet
            % Standard deviation of each image pixel in the 'time' of frame dimension (1 is weight, 3 is dim)
            imStd = std(imMatrix,1,3);                                          
            imStdFiltered = stdfilt(imStd,this.FilterKernal);    

            imStdFiltered = imStdFiltered/max(imStdFiltered(:)); % Normalise % Output 1

            % Ostu Level
            levelBinarize = graythresh(imStdFiltered);

            % Binarizing and processing (Constant can be modified)
            imbinarized = imbinarize(imStdFiltered, levelBinarize);                         % !0.12 works well for non-background   
            imbinarized_Opened_1 = bwareaopen(imbinarized,this.FactorOpen_1);               % Counter small white pixels outside of jet  % Output 2
            imComplemented = imcomplement(imbinarized_Opened_1);
            imbinarized_Opened_2 = bwareaopen(imComplemented,this.FactorOpen_2);            % Counter small white pixels inside of jet    % Output 3

            imEroded = bwmorph( imbinarized_Opened_2,'erode',this.FactorDE);
            imEroded(1:this.FactorDE,:) = true;
            imEroded(end-this.FactorDE+1:end,:) = true;
            imEroded = ~imfill(~imEroded, 'holes');   % Output 4

            % Counter white pixels inside jet after Erode 
            % However, it will prevent HTC detection
            imDilated = bwmorph(imEroded,'dilate',this.FactorDE);
            imDilatedComplemented = imcomplement(imDilated);
            imDilatedComplementedFilled = imfill(imDilatedComplemented,'holes'); % Output 5

            trueRegions = bwconncomp(imDilatedComplementedFilled);
            numPixels = cellfun(@numel,trueRegions.PixelIdxList);
            [~,idx] = max(numPixels);

            imOutput =  zeros(imHeight, imWidth);
            if ~isempty(idx)
                imOutput(trueRegions.PixelIdxList{idx}) = 1;
            else
                imOutput(round(imHeight/2), 1) = 1;                     % Prevent clash
            end

            % Boundary detection & Get coordinatesimage 
            jetBoundary = bwboundaries(imOutput);
            [~,idxMaxBoundary] = max(cellfun('size',jetBoundary,1));
            jetBoundary = jetBoundary{idxMaxBoundary};
            
            % Output
            this.JetBoundary = jetBoundary;
            this.JetBWImage = imOutput;
            
            % For debug
            imCompacted.imRaw = this.ImageCellArray{1}; % after BS
            imCompacted.imStd = imStd;
            imCompacted.imStdFiltered = imStdFiltered;
            imCompacted.imbinarized = imbinarized;
            imCompacted.imbinarized_Opened_1 = imbinarized_Opened_1;
            imCompacted.imbinarized_Opened_2 = imbinarized_Opened_2;
            imCompacted.imEroded = imEroded;
            imCompacted.imDilatedComplementedFilled = imDilatedComplementedFilled;
            this.ImCompacted = imCompacted;
        end
        
        %% Plot
        function ImShow(this, axisHandle)
            if isempty(this.JetBoundary)
                warning('Jet boundary is empty')
                return
            end

            jetBoundary = this.JetBoundary;
            if isempty(axisHandle)
                figure;
                imshow(this.ImageCellArray{1},[]); hold on;
                plot(jetBoundary(:,2), jetBoundary(:,1), 'linewidth',1.5,'color','r');
                drawnow
            else
                imshow(this.ImageCellArray{1},[],'Parent',axisHandle); hold on;
                plot(axisHandle, jetBoundary(:,2), jetBoundary(:,1), 'linewidth',1.5,'color','r');
                drawnow
            end
        end
        
    end
end

