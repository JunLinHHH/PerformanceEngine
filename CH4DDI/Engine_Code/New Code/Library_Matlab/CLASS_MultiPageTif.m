classdef CLASS_MultiPageTif
    % construct a cell contains selected frames from a single .tif file
    properties
        ImageFileName
        ImageWidth
        ImageHeight
        MaxFrame
    end
    
    properties (Access = private)
        TifStoreCell
        TifStoreCell_SubtractPrevious
        IntensityProfile
    end
    
    methods
        %% Constructor
        function this = CLASS_MultiPageTif(tifFileName)
            this.ImageFileName = tifFileName;
            if ~isfile(tifFileName)
                error("Selected file is not a valid .tif file")
            end
            warnSts = warning('off', 'imageio:tifftagsread:nextIfdPointerOutOfRange');
            warnSts.state = 'off';
            imageInfo = imfinfo(tifFileName);

            this.MaxFrame = length(imageInfo);
            this.ImageWidth = imageInfo(1).Width;
            this.ImageHeight = imageInfo(1).Height;
        end
        %% Store tif to cell (arg1: frame range) 
        % Raw image
        function this = StoreTifToCell(this, target_frames)
            % read multipage tif files to a cell
            % frameRange: int array
            % dim: int
                % 1: column vector
                % 2: row vector
            warnSts = warning('off', 'imageio:tifftagsread:nextIfdPointerOutOfRange');
            warnSts.state = 'off';
            imageInfo = imfinfo(this.ImageFileName);
            fileFrameMaximum = length(imageInfo);

            this.ImageWidth = imageInfo(1).Width;
            this.ImageHeight = imageInfo(1).Height;
            
            if max(target_frames) > this.MaxFrame
                error("Selected frame exceeds the file frame range (%d)",...
                    fileFrameMaximum);
            elseif min(target_frames) <=0
                error("Selected frame must be positve integers");
            end

            % Storage
            tifStoreCell = cell(length(target_frames),1); 
            
            % Read .tif into a column cell
            for ii = 1:length(target_frames)
                    imgFrame = imread(this.ImageFileName, target_frames(ii));
                    if length(size(imgFrame)) == 3
                        imgFrame = rgb2gray(imgFrame);
                    end
                    tifStoreCell{ii} = imgFrame;
            end
            this.TifStoreCell = tifStoreCell;
        end

        % Background subtracted (if frame 1, subtract by itself)
        function this = StoreTifToCell_SubtractPrevious(this, target_frames)
            warnSts = warning('off', 'imageio:tifftagsread:nextIfdPointerOutOfRange');
            warnSts.state = 'off';
            imageInfo = imfinfo(this.ImageFileName);
            fileFrameMaximum = length(imageInfo);

            this.ImageWidth = imageInfo(1).Width;
            this.ImageHeight = imageInfo(1).Height;
            
            if max(target_frames) > this.MaxFrame
                error("Selected frame exceeds the file frame range (%d)",...
                    fileFrameMaximum);
            elseif min(target_frames) <=0
                error("Selected frame must be positve integers");
            end
            
            % Storage
            tifStoreCell_SubtractPrevious = cell(length(target_frames),1); 
            
            % Read .tif into a column cell
            for ii = 1:length(target_frames)
                    this_frame = target_frames(ii);
                    if this_frame ~= 1
                        img_frame = imread(this.ImageFileName, this_frame);
                        img_frame_previous = imread(this.ImageFileName, this_frame - 1);
                    else
                        img_frame = imread(this.ImageFileName, this_frame);
                        img_frame_previous = img_frame;
                    end

                    if length(size(img_frame)) == 3
                        img_frame = rgb2gray(img_frame);
                        img_frame_previous = rgb2gray(img_frame_previous);
                    end

                    img_frame_subtraction = img_frame - img_frame_previous;

                    tifStoreCell_SubtractPrevious{ii} = img_frame_subtraction;
            end
            this.TifStoreCell_SubtractPrevious = tifStoreCell_SubtractPrevious;
        end

        %% Get 1D array intensity change with respect to frame
        function this = CalculateIntensityProfile(this)
            maxFrame = this.MaxFrame;
            intensityProfile = zeros(maxFrame, 1);
            for ii = 1:length(intensityProfile)
                img_this_frame = imread(this.ImageFileName, ii);
                [img_height, img_width] = size(img_this_frame);
                img_this_frame = img_this_frame(1:round(img_height/2), 1:(img_width/2));
                intensityProfile(ii) = sum(sum(img_this_frame));
            end
            this.IntensityProfile = intensityProfile;
        end

        %% Getter
        function tifStoreCell = GetTifStoreCell(this)
            tifStoreCell = this.TifStoreCell;
        end

        function tifStoreCell_SubtractPrevious = GetTifStoreCell_SubtractPrevious(this)
            tifStoreCell_SubtractPrevious = this.TifStoreCell_SubtractPrevious;
        end

        function intensityProfile = GetIntensityProfile(this)
            intensityProfile = this.IntensityProfile;
        end        
    end
end

