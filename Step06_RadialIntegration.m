%% -- Load data
% If .mat exist
clearvars; close all; clc; 
Configs = Configuration();
Loaded = load(Configs.DataObjReadDir);

obj = Loaded.obj;
T_DataAll = obj.DataMatrix;
T_Condition_ROI = obj.ConditionMatrix;

%%
T_Condition_ROI.FlameRecession = cell(height(T_Condition_ROI), 1);
T_Condition_ROI.FlamePenetration = cell(height(T_Condition_ROI), 1);
%% Target Runs and Frames
SelectedConditions = 10:20; % 1:height(T_Condition_ROI);
MaxFrame = 600;
MaxLengthPixel = 700; % how long is the axis distance (pixel)
EdgeOffset = 50; % 600 frame is EOI, cut last 50 frame to prevent wrong detection at edge
bPlot = 0;
for  ss = SelectedConditions

    TargetRuns = T_Condition_ROI{ss,'RunNumbers'}{1};
    TargetFrames = (1:MaxFrame)';
    TargetFrames = repmat(TargetFrames,1,length(TargetRuns));
    
    if length(TargetRuns) ~= size(TargetFrames,2)
        fprintf('Frame column number does not match run number\m')
        return
    end
    
    row_number_store = size(TargetFrames,1);
    column_number_store = length(TargetRuns);
    
    fprintf('Processing condition #%03d: \n', ss)
    
    %% Storage
    % For all runs at this conditions
    store_sum_all_runs = zeros(length(TargetFrames)-EdgeOffset, MaxLengthPixel);
    store_recession = cell(length(TargetRuns),1);
    store_penetration = cell(length(TargetRuns),1);

    for rr = 1:length(TargetRuns)
        
        % load image and hrr arrys
        run_number = TargetRuns(rr);
        fprintf('Processing run #%03d: \n', run_number)
        this_frames = TargetFrames(:, rr);
        img_dir =  T_DataAll(T_DataAll.Number == run_number, :).File_Image{1};
    
        % load constant parameters (run-based)
        nozzle_position_x = T_DataAll(T_DataAll.Number == run_number, :).Nozzle_X;
        % pixel_resolution = T_DataAll(T_DataAll.Number == run_number, :).Pixel_res;
        % injection_delay_ms = T_DataAll(T_DataAll.Number == run_number, :).Inj_delay;
        % FPS = T_DataAll(T_DataAll.Number == run_number, :).FPS;
        %     time_scale_ms = (this_frames-1)/FPS/1000  - injection_delay_ms; 
    
        % Read image
        obj_img = CLASS_MultiPageTif(img_dir);
        obj_img = obj_img.StoreTifToCell(this_frames);
        img_loaded = obj_img.GetTifStoreCell();
        
        % Crop image
        for jj = 1:numel(img_loaded)
            temp = img_loaded{jj};
            temp(:,1:nozzle_position_x) = [];
            temp(:, MaxLengthPixel+1:end) = []; % keep 
            % temp = CLASS_Utilis.ImageBitShift(temp,1); % !!!!!!!!!!!!Bitshift
            img_loaded{jj} = temp;
        end
    
        img_size = size(temp);
    
        % Storage (for each run)
        store_RI_img = zeros(length(this_frames),img_size(2));
        store_RI_imBgSubtract = zeros(length(this_frames),img_size(2));
        store_RI_imTempStdFiltered = zeros(length(this_frames),img_size(2));
        store_RI_imSpatialStdFiltered = zeros(length(this_frames),img_size(2));
        
        % Mask - normal distribution
        valve_std = 15;
        desired_max_value = 0.03;
        y = GenerateNormalDistribution(valve_std, size(temp,1), desired_max_value, bPlot)';
        mask_normal_dis =  repmat(y,1,img_size(2));
        FilterKernal = ones(5,5);
        
        % Processing each frame
        fprintf('\tCurrent Frame:      ')
        for ii = 1:numel(img_loaded)
            fprintf('\b\b\b\b%04d',ii)
            if ii == 1 || ii == numel(img_loaded)
                continue
            end
            
            imMatrix = zeros(img_size(1), img_size(2), 3);
        
            mask_apply =  mask_normal_dis;
            % Masking
            imMatrix(:,:,1) = double(img_loaded{ii-1}) .*mask_apply;
            imMatrix(:,:,2) = double(img_loaded{ii}) .*mask_apply;  % Current frame
            imMatrix(:,:,3) = double(img_loaded{ii+1}) .*mask_apply; 
    
            imBgSubtract = imMatrix(:,:,2) - imMatrix(:,:,1);
    
            imTempStdFiltered = std(imMatrix,1,3);                                          
            imSpatialStdFiltered = stdfilt(imTempStdFiltered,FilterKernal);    
    
            img_sum = sum(img_loaded{ii},1);  % Raw image 
            imBgSubtract_sum = sum(imBgSubtract,1); % Bg subtraction
            imTempStdFiltered_sum = sum(imTempStdFiltered,1); % Temp std
            imSpatialStdFiltered_sum = sum(imSpatialStdFiltered,1);  % Spatial std
            
            % Storage
            store_RI_img(ii, :) = img_sum;
            store_RI_imBgSubtract(ii, :) = imBgSubtract_sum;
            store_RI_imTempStdFiltered(ii, :) = imTempStdFiltered_sum;
            store_RI_imSpatialStdFiltered(ii, :) = imSpatialStdFiltered_sum;
        end
        fprintf('\n')
        %%
        % f1 = figure; f1.Name = 'store_RI_img';
        % imshow(store_RI_img, [])
        % axis on
        % f2 = figure; f2.Name = 'store_RI_imBgSubtract';
        % imshow(store_RI_imBgSubtract, [])
        % axis on
        % f3 = figure; f3.Name = 'store_RI_imTempStdFiltered';
        % imshow(store_RI_imTempStdFiltered, [])
        % axis on
        % f4 = figure; f4.Name = 'store_RI_imSpatialStdFiltered';
        % imshow(store_RI_imSpatialStdFiltered, [])
        % axis on
        %%
        im_binary = imbinarize(uint8(store_RI_imSpatialStdFiltered));
        im_binary(:, 1:100) = 1; % Cut near nozzle side
        im_binary(end-EdgeOffset+1:end, :) = [];
        se = strel('disk',1);
        im_close = imclose(im_binary, se);
        im_open = bwareaopen(~im_close, 1000);
        im_open = imfill(im_open,'holes');
    
        %keyboard
        % Get largest region
        trueRegions = bwconncomp(im_open);
        numPixels = cellfun(@numel,trueRegions.PixelIdxList);
        [~,idx] = max(numPixels);
    
        im_final_BW =  zeros(length(TargetFrames)-EdgeOffset, MaxLengthPixel);
        if ~isempty(idx)
            im_final_BW(trueRegions.PixelIdxList{idx}) = 1;
        else
            im_final_BW(round(imHeight/2), 1) = 1;                     % Prevent clash
        end
        
        [time_recession_frame, data_recession_pixel] = GetFlameBase(im_final_BW,bPlot);
        [time_penetration_frame, data_penetration_pixel] = GetFlamePenetration(im_final_BW, bPlot);
        if bPlot    
            figure;hold on;
            plot(time_recession_frame, smoothdata(data_recession_pixel,'gaussian',15),'r-'); 
            plot(time_penetration_frame, data_penetration_pixel,'k-');
            figure; imshow(im_final_BW,[])
            axis on
        end

        % keyboard

        store_sum_all_runs = store_sum_all_runs + im_open;
        store_recession{rr} = [time_recession_frame', data_recession_pixel'];
        store_penetration{rr} = [time_penetration_frame', data_penetration_pixel'];    
    end
    T_Condition_ROI.FlameRecession{ss} = store_recession;
    T_Condition_ROI.FlamePenetration{ss} = store_penetration;
    obj.ConditionMatrix= T_Condition_ROI;

end
%%
keyboard
DataObjSaveDir = Configs.DataObjSaveDir;     
save(DataObjSaveDir,"obj")
fprintf("Saving finished\n")
%%
function y = GenerateNormalDistribution(valve_std, length, desired_max_value, bPlot)
    std_val = valve_std;
    scale_min = 1;
    scale_max = length;
    mean_val = (scale_max - scale_min) / 2;
    % valve_std = 30; % Specify the standard deviation as an argument
    
    % Generate the scale
    x = scale_min:scale_max;
    
    % Generate the normal distribution
    y = (1 / (std_val * sqrt(2 * pi))) * exp(-0.5 * ((x - mean_val) / std_val) .^ 2);
    scaling_factor = desired_max_value / max(y);
    y = y * scaling_factor;
    % Plot the normal distribution
    if bPlot
        figure;
        plot(x, y, 'LineWidth', 2);
        xlabel('Scale');
        ylabel('Probability Density');
        title(['Normal Distribution from ' num2str(scale_min) ' to ' num2str(scale_max)]);
        grid on;
    end
end

function [time_frame, flamebase_pixel] = GetFlameBase(im_input,bPlot)
    ref_scale = 1: size(im_input,1); % y scale (frame)
    data = zeros(size(ref_scale));
    for ii = ref_scale
        this_line = im_input(ii, :);
        idx_true = find(this_line ==1,1, 'first');
        if isempty(idx_true)
            data(ii) = 1;
        else
            data(ii) = idx_true;
        end
    end

    time_frame = ref_scale;
    flamebase_pixel = data;

    if bPlot
        plot(time_frame,flamebase_pixel)
    end
end

function [time_frame_transfer, flamepenetration_pixel_transfer] = GetFlamePenetration(im_input,bPlot)
    ref_scale = 1: size(im_input,2); % x scale (axial distance pixel)
    data = zeros(size(ref_scale));
    for ii = ref_scale
        this_line = im_input(:, ii);
        idx_true = find(this_line ==1,1, 'first');
        if isempty(idx_true)
            data(ii) = size(im_input,1);
        else
            data(ii) = idx_true;
        end
    end

    time_frame = data;
    flamepenetration_pixel = ref_scale;

    if bPlot
        plot(time_frame,flamepenetration_pixel)
    end

    % Axis transfer 
    [unique_time, ~, idx] = unique(time_frame);
    v_max = accumarray(idx, flamepenetration_pixel, [], @max);
    time_frame_transfer = 1:size(im_input,1);
    flamepenetration_pixel_transfer = interp1(unique_time, v_max, time_frame_transfer);
    if bPlot
        plot(time_frame,flamepenetration_pixel,'kx')
        hold on;
        plot(time_frame_transfer,flamepenetration_pixel_transfer,'ro')
    end
end


% function [x_min, y_min] = GetFlameBase(im_input,bPlot)
%     edges = edge((im_input), 'Canny');
%     [x,y] = find(edges);
% 
%     % Find the unique x values and their corresponding indices
%     [uniqueX, ~, idx] = unique(x);
% 
%     % Use accumarray to find the minimum y for each unique x
%     minY = accumarray(idx, y, [], @min);
% 
%     % Combine the unique x values with their corresponding minimum y values
%     minEdgeData = [uniqueX, minY];
%     % Separate the result into x and y vectors
%     x_min = minEdgeData(:, 1);
%     y_min = minEdgeData(:, 2);
%     if bPlot 
%         figure; hold on;
%         plot(x, y, 'g.', 'MarkerSize', 1);
%         plot(x_min, y_min, 'ro', 'MarkerSize', 1);
%     end
% end
% 
% function [y_min, x_min] = GetFlamePenetration(im_input, bPlot)
%     edges = edge((im_input), 'Canny');
%     [y,x] = find(edges);
% 
%     % Find the unique x values and their corresponding indices
%     [uniqueX, ~, idx] = unique(x);
% 
%     % Use accumarray to find the minimum y for each unique x
%     minY = accumarray(idx, y, [], @min);
% 
%     % Combine the unique x values with their corresponding minimum y values
%     minEdgeData = [uniqueX, minY];
%     % Separate the result into x and y vectors
%     x_min = minEdgeData(:, 1);
%     y_min = minEdgeData(:, 2);
%     if bPlot 
%         figure; hold on;
%         plot(y, x, 'g.', 'MarkerSize', 1);
%         plot(y_min, x_min, 'ro', 'MarkerSize', 1);
%     end
% end