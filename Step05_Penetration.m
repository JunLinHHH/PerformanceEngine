%% -- Load data
% If .mat exist
clearvars; close all; clc; 
Configs = Configuration();
Loaded = load(Configs.DataObjReadDir);

obj = Loaded.obj;
T_DataAll = obj.DataMatrix;
T_Condition_ROI = obj.ConditionMatrix;

%% Target Runs and Frames
SelectedConditions = [10:20];
MaxFrame = 600;
MaxLengthPixel = 700; % how long is the axis distance (pixel)
for ss = SelectedConditions
    TargetRuns = T_Condition_ROI{ss,'RunNumbers'}{1};
    TargetFrames = (1:300)';
    TargetFrames = repmat(TargetFrames,1,length(TargetRuns));
    
    if length(TargetRuns) ~= size(TargetFrames,2)
        fprintf('Frame column number does not match run number\m')
        return
    end
    
    row_number_store = size(TargetFrames,1);
    column_number_store = length(TargetRuns);
    
    fprintf('Processing condition #%03d: \n', ss)
    %% Storage
    image_raw_store = cell(row_number_store, column_number_store);
    boundary_store_trans = cell(row_number_store, column_number_store);
    boundary_store_black = cell(row_number_store, column_number_store);
    array_time_scale_frame = zeros(size(boundary_store_black));
    %%
    for rr = 1:length(TargetRuns)
        
        % load image and hrr arrys
        run_number = TargetRuns(rr);
        fprintf('...Processing %03d/%03d, run: %03d\n', rr, length(TargetRuns), TargetRuns(rr));
        selected_frames = TargetFrames(:, rr);
        img_dir =  T_DataAll(T_DataAll.Number == run_number, :).File_Image{1};
    
        % load constant parameters (run-based)
        nozzle_position_x = T_DataAll(T_DataAll.Number == run_number, :).Nozzle_X;
        pixel_resolution = T_DataAll(T_DataAll.Number == run_number, :).Pixel_res;
            time_scale_frame = TargetFrames(:, rr); 
        fprintf('\tCurrent Frame:      ')
        for ff = 1:size(TargetFrames,1)
            thisFrame = selected_frames(ff);
            
            fprintf('\b\b\b\b%04d',ff)
            if thisFrame == 1 || thisFrame == max(max(TargetFrames))
                % Skip first & last since they have no 3 frame comparison
                continue
            end
    
            thisFrame_expanded =  GetThreeConsecutiveFrames(thisFrame);
            
            obj_img = CLASS_MultiPageTif(img_dir);
            obj_img = obj_img.StoreTifToCell(thisFrame_expanded);
            image_frame3consec = obj_img.GetTifStoreCell();
    
            % Cut nozzle (before boundary processing)
            for jj = 1:numel(image_frame3consec)
                temp = image_frame3consec{jj};
                temp(:,1:nozzle_position_x) = [];
                temp(:, MaxLengthPixel+1:end) = []; % keep 
                image_frame3consec{jj} = temp;
            end
            [im_height_cut, im_width_cut] = size(temp);
            % Store image
            image_raw_store{ff,rr} = image_frame3consec{2};
    
            % Store boundary
            obj_schlieren_trans = CLASS_ImageSchlieren();
            obj_schlieren_trans.FactorBS = 1; % can change by runs
            obj_schlieren_trans.FactorDE = 4;
            obj_schlieren_trans = obj_schlieren_trans.JetBoundaryProcess(image_frame3consec);
            boundary_frame_trans = obj_schlieren_trans.GetBoundary();
            boundary_store_trans{ff,rr} = boundary_frame_trans;
            
            obj_schlieren_black = CLASS_ImageSchlieren();
            % 0 for non-reacting
            % 2 for 10, 18, 19, 20 pilot-main
            % 3 for 13, 14 ,15 ,16, 17,  11, 12 main pilot, T var, O2 var
            if ismember(ss, 29)
                FactorBS = 0;
                SizeFilter = 5;
            elseif ismember(ss, [10,18,19,20])
                FactorBS = 3;
                SizeFilter = 21;
            elseif ss == 3
                FactorBS = 3;
                SizeFilter = 15;
            else
                FactorBS = 3;
                SizeFilter = 21;
            end
            obj_schlieren_black.FactorBS = FactorBS;  
            obj_schlieren_black.SizeFilter = SizeFilter;
            obj_schlieren_black.FactorOpen_1 = 50;
            obj_schlieren_black.FactorOpen_2 = 100;
            obj_schlieren_black.FactorDE = 5;
            obj_schlieren_black = obj_schlieren_black.JetBoundaryProcess(image_frame3consec);
            boundary_frame_black = obj_schlieren_black.GetBoundary();
            boundary_store_black{ff,rr} = boundary_frame_black;
        end
        fprintf('\n')
        array_time_scale_frame(:,rr) = time_scale_frame;
    end
    
    [array_penetration, array_recession] = Boundary_to_Penetration(boundary_store_black);
    
    %% Plot
    % % Penetration
    % figure; hold on
    % for ii = 1:size(array_penetration, 2)
    %     plot(array_time_scale_frame(:,ii), array_penetration(:,ii) * pixel_resolution)
    % end
    % % Recession
    % figure; hold on
    % for ii = 1:size(array_recession, 2)
    %     plot(array_time_scale_frame(:,ii), array_recession(:,ii) * pixel_resolution)
    % end
    % % Penetration & Recession
    % fff = figure; hold on; fff.Name = sprintf('set: %03d', ss);
    % for ii = 1:size(array_recession, 2)
    %     plot(array_time_scale_frame(:,ii), array_penetration(:,ii) * pixel_resolution)
    %     plot(array_time_scale_frame(:,ii), array_recession(:,ii) * pixel_resolution)
    % end
    % mean_array_penetration = mean(array_penetration,2);
    % mean_array_recession = mean(array_recession,2);
    % plot(array_time_scale_frame(:,ii), mean_array_penetration * pixel_resolution,'LineWidth',2,'Color','b');
    % plot(array_time_scale_frame(:,ii), mean_array_recession * pixel_resolution,'LineWidth',2,'Color','r');
    % ylim([-10,100])
    %% Output
    PenetrationRecession.array_time_scale_frame = array_time_scale_frame;
    PenetrationRecession.array_penetration = array_penetration;
    PenetrationRecession.array_recession = array_recession;
    PenetrationRecession.runs = TargetRuns;
    %% Save
    keyboard
    close all
    save_name = sprintf('Penetration_set_%03d', ss);
    fullfile_output = fullfile(Configs.MatDataOutput, save_name);
    save(fullfile_output,'PenetrationRecession');
    fprintf('Save completed\n')
end
%%
function three_frames = GetThreeConsecutiveFrames(frame_number)
    previous_frame = frame_number - 1;
    next_frame = frame_number + 1;
    three_frames = sort([previous_frame; frame_number; next_frame]);
end

function [array_penetration, array_recession] = Boundary_to_Penetration(boundary_store_black)
    array_penetration = zeros(size(boundary_store_black));
    array_recession = zeros(size(boundary_store_black));

    for ii = 1: size(boundary_store_black, 2)
        for jj = 1: size(boundary_store_black, 1)
            this_boundary = boundary_store_black(jj, ii);
            this_boundary = this_boundary{1};
            if isempty(this_boundary)
                array_penetration(jj, ii) = 0;
                array_recession(jj, ii) = 0;
            else
                pixel_penetration = max(this_boundary(:,2));
                array_penetration(jj, ii) = pixel_penetration;

                pixel_recession = min(this_boundary(:,2));
                array_recession(jj, ii) = pixel_recession;
            end
        end
    end
end

