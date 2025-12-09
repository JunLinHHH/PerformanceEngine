% Extract centroid of black circle marker from 4-panel figure
% Focus on detecting only the black circle marker centroid using RGB color detection
clear; clc; close all;

%% Configuration - UPDATE THIS PATH
img_path_full = 'C:\Users\jun-y\OneDrive - UNSW\CH4DDI\CH4DDI vs H2DDI CVCC\Patrick2022_handled\Patrick2022_Figure10.png';

% Subplot names (will be used for saving)
subplot_names = {'H-0.07ms-D', 'D-0.93ms-H', 'D-1.93ms-H', 'D-2.93ms-H'};

% Load FULL image
img = imread(img_path_full);
[height, width, ~] = size(img);

% Get image path for saving
[img_path, ~, ~] = fileparts(img_path_full);

% Axis ranges for each subplot (all have same range)
x_min_data = 0;
x_max_data = 5;
y_min_data = 0;
y_max_data = 95;

%% Check if boundaries file exists and ask user
boundary_file = fullfile(img_path, 'plot_boundaries_Patrick2022_Figure10.mat');
load_existing = false;

if exist(boundary_file, 'file')
    % Create a simple dialog to ask user
    answer = questdlg('Previous plot boundaries found. Do you want to load them?', ...
        'Load Existing Boundaries', ...
        'Yes', 'No', 'Yes');
    
    if strcmp(answer, 'Yes')
        load(boundary_file, 'plot_boundaries', 'subplot_names');
        load_existing = true;
        fprintf('Loaded existing boundaries from: %s\n', boundary_file);
    end
end

%% Storage for boundaries (in FULL image coordinates)
if ~load_existing
    plot_boundaries = cell(4, 1);
end

%% Interactive boundary selection - show FULL image
if ~load_existing
    fprintf('\n========================================\n');
    fprintf('BLACK CIRCLE CENTROID DETECTION\n');
    fprintf('========================================\n');
    fprintf('You will see the FULL 4-panel figure.\n');
    fprintf('For each subplot, click 4 points to define the plot boundary:\n');
    fprintf('  1. TOP-LEFT corner of the plot area\n');
    fprintf('  2. TOP-RIGHT corner of the plot area\n');
    fprintf('  3. BOTTOM-RIGHT corner of the plot area\n');
    fprintf('  4. BOTTOM-LEFT corner of the plot area\n');
    fprintf('========================================\n\n');
else
    fprintf('\n========================================\n');
    fprintf('Using previously saved boundaries...\n');
    fprintf('========================================\n\n');
end

if ~load_existing
    for sp = 1:4
    fprintf('\n=== Subplot %d: %s ===\n', sp, subplot_names{sp});
    
    % Create figure showing FULL image
    fig = figure('Name', sprintf('Select Boundary - %s (Subplot %d/4)', subplot_names{sp}, sp), ...
                 'Position', [50, 50, 1400, 800]);
    imshow(img);
    hold on;
    
    % Show previously selected boundaries
    for i = 1:sp-1
        bounds = plot_boundaries{i};
        rectangle('Position', [bounds(1), bounds(3), ...
            bounds(2)-bounds(1), bounds(4)-bounds(3)], ...
            'EdgeColor', 'g', 'LineWidth', 2, 'LineStyle', '--');
        text(bounds(1), bounds(3)-10, sprintf('Subplot %d (done)', i), ...
            'Color', 'g', 'FontSize', 10, 'FontWeight', 'bold');
    end
    
    title(sprintf('Click 4 corners: 1=TOP-LEFT, 2=TOP-RIGHT, 3=BOTTOM-RIGHT, 4=BOTTOM-LEFT\nSubplot %d/%d: %s', ...
        sp, 4, subplot_names{sp}), 'Interpreter', 'none', 'FontSize', 14);
    hold off;
    
    % Get 4 corner points
    fprintf('Click the 4 corners of the plot area (in order)...\n');
    [x_clicks, y_clicks] = ginput(4);
    
    % Extract boundary from clicks
    x_min_pixel = min(x_clicks);
    x_max_pixel = max(x_clicks);
    y_min_pixel = min(y_clicks);
    y_max_pixel = max(y_clicks);
    
    % Store boundaries (in full image coordinates)
    plot_boundaries{sp} = [x_min_pixel, x_max_pixel, y_min_pixel, y_max_pixel];
    
    % Show selected region
    hold on;
    rectangle('Position', [x_min_pixel, y_min_pixel, ...
        x_max_pixel-x_min_pixel, y_max_pixel-y_min_pixel], ...
        'EdgeColor', 'r', 'LineWidth', 3);
    plot(x_clicks, y_clicks, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    text(x_min_pixel, y_min_pixel-20, sprintf('Subplot %d', sp), ...
        'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', 'white');
    hold off;
    
    fprintf('Selected boundary: x[%.1f, %.1f], y[%.1f, %.1f]\n', ...
        x_min_pixel, x_max_pixel, y_min_pixel, y_max_pixel);
    fprintf('Plot size: %.0f x %.0f pixels\n', ...
        x_max_pixel - x_min_pixel, y_max_pixel - y_min_pixel);
    
    % Ask user to confirm
    answer = questdlg('Is this boundary correct?', ...
        'Confirm Boundary', ...
        'Yes', 'No, redo', 'Yes');
    
    if strcmp(answer, 'No, redo')
        sp = sp - 1;  % Redo this subplot
        close(fig);
        continue;
    end
    
    pause(0.5);  % Brief pause before next subplot
    close(fig);
end

% Show final result with all boundaries
figure('Name', 'All Boundaries Selected', 'Position', [50, 50, 1400, 800]);
imshow(img);
hold on;
for sp = 1:4
    bounds = plot_boundaries{sp};
    rectangle('Position', [bounds(1), bounds(3), ...
        bounds(2)-bounds(1), bounds(4)-bounds(3)], ...
        'EdgeColor', 'g', 'LineWidth', 3);
    text(bounds(1)+10, bounds(3)+20, sprintf('%d: %s', sp, subplot_names{sp}), ...
        'Color', 'g', 'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', 'white');
end
if load_existing
    title('Loaded Plot Boundaries', 'FontSize', 14);
else
    title('All Plot Boundaries Selected', 'FontSize', 14);
end
hold off;
end

fprintf('\n========================================\n');
fprintf('Boundary selection complete!\n');
fprintf('Now processing black circle detection...\n');
fprintf('========================================\n\n');

%% Process each subplot with selected boundaries - DETECT BLUE FIRST, THEN BLACK CIRCLE
circle_centroids = cell(4, 1);
detection_results = cell(4, 1);

for sp = 1:4
    fprintf('\n=== Processing subplot %d: %s ===\n', sp, subplot_names{sp});
    
    % Get boundary in full image coordinates
    x_min_pixel_full = plot_boundaries{sp}(1);
    x_max_pixel_full = plot_boundaries{sp}(2);
    y_min_pixel_full = plot_boundaries{sp}(3);
    y_max_pixel_full = plot_boundaries{sp}(4);
    
    % Extract subplot region from FULL image
    y_min_idx = max(1, round(y_min_pixel_full));
    y_max_idx = min(height, round(y_max_pixel_full));
    x_min_idx = max(1, round(x_min_pixel_full));
    x_max_idx = min(width, round(x_max_pixel_full));
    
    subplot_img = img(y_min_idx:y_max_idx, x_min_idx:x_max_idx, :);
    
    fprintf('Extracted region: x[%d:%d], y[%d:%d]\n', x_min_idx, x_max_idx, y_min_idx, y_max_idx);
    fprintf('Subplot size: %d x %d pixels\n', size(subplot_img, 2), size(subplot_img, 1));
    
    % Now work in subplot coordinates (relative to extracted region)
    x_min_pixel = 1;
    x_max_pixel = size(subplot_img, 2);
    y_min_pixel = 1;
    y_max_pixel = size(subplot_img, 1);
    
    % Create plot mask (entire subplot region)
    plot_mask = true(size(subplot_img, 1), size(subplot_img, 2));
    
    %% STEP 1: DETECT BLUE MARKERS FIRST
    fprintf('Step 1: Detecting blue markers...\n');
    
    % Extract blue markers
    blue_mask = (subplot_img(:,:,3) > 150) & (subplot_img(:,:,1) < 150) & (subplot_img(:,:,2) < 150);
    blue_mask = blue_mask & plot_mask;
    blue_mask = bwareaopen(blue_mask, 5);
    
    % Find blue marker centroids
    cc_blue = bwconncomp(blue_mask);
    stats_blue = regionprops(cc_blue, 'Area', 'Centroid', 'PixelIdxList');
    
    blue_centroids = [];
    blue_exclusion_mask = false(size(blue_mask));
    
    for i = 1:length(stats_blue)
        if stats_blue(i).Area > 5 && stats_blue(i).Area < 300
            blue_centroids(end+1, :) = stats_blue(i).Centroid;
            % Create exclusion zone around blue marker (dilate)
            temp_mask = false(size(blue_mask));
            temp_mask(stats_blue(i).PixelIdxList) = true;
            blue_exclusion_mask = blue_exclusion_mask | imdilate(temp_mask, strel('disk', 15));
        end
    end
    
    fprintf('Found %d blue markers\n', size(blue_centroids, 1));
    
    %% STEP 2: BLACK CIRCLE DETECTION - Exclude blue markers and use RGB filtering
    fprintf('Step 2: Detecting black circle marker (excluding blue markers, filtering by color)...\n');
    
    % Method 1: RGB-based black detection
    % Black should have low values in ALL RGB channels
    red_channel = subplot_img(:,:,1);
    green_channel = subplot_img(:,:,2);
    blue_channel = subplot_img(:,:,3);
    
    % Black mask: low in all RGB channels (not red, not green, not blue)
    black_rgb_mask = (red_channel < 80) & (green_channel < 80) & (blue_channel < 80);
    
    % Also exclude red markers: high red, low green and blue
    red_mask = (red_channel > 150) & (green_channel < 100) & (blue_channel < 100);
    red_exclusion_mask = imdilate(red_mask, strel('disk', 10));
    
    % Also exclude green markers: high green, low red and blue  
    green_mask = (green_channel > 150) & (red_channel < 100) & (blue_channel < 100);
    green_exclusion_mask = imdilate(green_mask, strel('disk', 10));
    
    % Combine all exclusions
    total_exclusion_mask = blue_exclusion_mask | red_exclusion_mask | green_exclusion_mask;
    
    % Final black mask
    black_mask = black_rgb_mask & plot_mask & ~total_exclusion_mask;
    
    % Remove noise (very small components)
    black_mask = bwareaopen(black_mask, 20);
    
    fprintf('  Red exclusion regions: %d pixels\n', sum(red_exclusion_mask(:)));
    fprintf('  Green exclusion regions: %d pixels\n', sum(green_exclusion_mask(:)));
    fprintf('  Blue exclusion regions: %d pixels\n', sum(blue_exclusion_mask(:)));
    fprintf('  Remaining black pixels: %d\n', sum(black_mask(:)));
    
    % Find connected components
    cc_black = bwconncomp(black_mask);
    stats_black = regionprops(cc_black, 'Area', 'Circularity', 'Eccentricity', ...
                              'Centroid', 'BoundingBox', 'PixelIdxList');
    
    % Filter for circular objects (potential circles)
    circle_candidates = [];
    
    for i = 1:length(stats_black)
        area = stats_black(i).Area;
        
        % Area filter - looking for medium-sized objects
        if area >= 50 && area <= 800
            % Check circularity and eccentricity for round objects
            is_circular = false;
            
            if isfield(stats_black(i), 'Circularity') && ~isempty(stats_black(i).Circularity)
                if stats_black(i).Circularity > 0.6  % High circularity
                    is_circular = true;
                end
            end
            
            if isfield(stats_black(i), 'Eccentricity') && ~isempty(stats_black(i).Eccentricity)
                if stats_black(i).Eccentricity < 0.5  % Low eccentricity (round)
                    is_circular = true;
                end
            end
            
            % Check aspect ratio of bounding box
            bbox = stats_black(i).BoundingBox;
            aspect_ratio = bbox(4) / bbox(3);  % height/width
            if aspect_ratio > 0.7 && aspect_ratio < 1.3  % Nearly square
                is_circular = true;
            end
            
            if is_circular
                circle_candidates(end+1) = i;
            end
        end
    end
    
    fprintf('Found %d circle candidates\n', length(circle_candidates));
    
    % Select the best circle candidate (largest area among candidates)
    if ~isempty(circle_candidates)
        areas = [stats_black(circle_candidates).Area];
        [~, best_idx] = max(areas);
        best_circle_idx = circle_candidates(best_idx);
        
        % Get centroid of the best circle
        centroid_px = stats_black(best_circle_idx).Centroid;
        circle_area = stats_black(best_circle_idx).Area;
        
        % Convert centroid to data coordinates
        [centroid_x_data, centroid_y_data] = pixel_to_data(centroid_px(1), centroid_px(2), ...
            x_min_pixel, x_max_pixel, y_min_pixel, y_max_pixel, ...
            x_min_data, x_max_data, y_min_data, y_max_data);
        
        % Store centroid information
        circle_centroids{sp} = struct('x_pixel', centroid_px(1), 'y_pixel', centroid_px(2), ...
                                     'x_data', centroid_x_data, 'y_data', centroid_y_data, ...
                                     'area', circle_area);
        
        fprintf('Black circle detected:\n');
        fprintf('  Centroid (pixels): (%.1f, %.1f)\n', centroid_px(1), centroid_px(2));
        fprintf('  Centroid (data): (%.3f, %.3f)\n', centroid_x_data, centroid_y_data);
        fprintf('  Area: %.0f pixels\n', circle_area);
        
        % Create visualization mask for the detected circle
        circle_mask = false(size(black_mask));
        circle_mask(stats_black(best_circle_idx).PixelIdxList) = true;
        
    else
        fprintf('WARNING: No black circle detected in subplot %d\n', sp);
        circle_centroids{sp} = struct('x_pixel', NaN, 'y_pixel', NaN, ...
                                     'x_data', NaN, 'y_data', NaN, ...
                                     'area', 0);
        circle_mask = false(size(black_mask));
    end
    
    %% STEP 3: CALCULATE TIME DIFFERENCES
    fprintf('Step 3: Calculating time differences...\n');
    
    time_differences = [];
    blue_data_coords = [];
    
    if ~isempty(blue_centroids) && ~isnan(circle_centroids{sp}.x_data)
        % Convert blue marker centroids to data coordinates
        for b = 1:size(blue_centroids, 1)
            [blue_x_data, blue_y_data] = pixel_to_data(blue_centroids(b,1), blue_centroids(b,2), ...
                x_min_pixel, x_max_pixel, y_min_pixel, y_max_pixel, ...
                x_min_data, x_max_data, y_min_data, y_max_data);
            blue_data_coords(end+1, :) = [blue_x_data, blue_y_data];
            
            % Calculate time difference (x-value difference)
            time_diff = abs(blue_x_data - circle_centroids{sp}.x_data);
            time_differences(end+1) = time_diff;
            
            fprintf('  Blue marker %d: (%.3f, %.3f), Time diff: %.3f ms\n', ...
                    b, blue_x_data, blue_y_data, time_diff);
        end
        
        if ~isempty(time_differences)
            fprintf('  Minimum time difference: %.3f ms\n', min(time_differences));
            fprintf('  Maximum time difference: %.3f ms\n', max(time_differences));
        end
    else
        fprintf('  Cannot calculate time differences (missing markers)\n');
    end
    
    % Store all detection results
    detection_results{sp} = struct('circle_centroid', circle_centroids{sp}, ...
                                  'blue_centroids', blue_data_coords, ...
                                  'time_differences', time_differences);
    
    %% DON'T SAVE - Just display results
    fprintf('Results (not saved):\n');
    fprintf('  Circle centroid: (%.3f, %.3f)\n', circle_centroids{sp}.x_data, circle_centroids{sp}.y_data);
    if ~isempty(blue_data_coords)
        fprintf('  Blue markers: %d found\n', size(blue_data_coords, 1));
        fprintf('  Time differences: [%s] ms\n', num2str(time_differences, '%.3f '));
    end
    
    %% Visualization for this subplot
    figure('Position', [100 + sp*50, 100, 1400, 900]);
    
    subplot(3,4,1);
    imshow(subplot_img);
    title(['Original: ' subplot_names{sp}], 'Interpreter', 'none');
    
    subplot(3,4,2);
    imshow(blue_mask);
    title('Blue Marker Mask');
    
    subplot(3,4,3);
    imshow(red_mask);
    title('Red Marker Mask');
    
    subplot(3,4,4);
    imshow(green_mask);
    title('Green Marker Mask');
    
    subplot(3,4,5);
    imshow(blue_exclusion_mask);
    title('Blue Exclusion Zone');
    
    subplot(3,4,6);
    imshow(red_exclusion_mask);
    title('Red Exclusion Zone');
    
    subplot(3,4,7);
    imshow(total_exclusion_mask);
    title('All Exclusion Zones');
    
    subplot(3,4,8);
    imshow(black_mask);
    title('Black Mask (Color Filtered)');
    
    subplot(3,4,9);
    imshow(circle_mask);
    title('Detected Black Circle');
    
    subplot(3,4,10);
    imshow(subplot_img);
    hold on;
    % Show blue markers
    if ~isempty(blue_centroids)
        plot(blue_centroids(:,1), blue_centroids(:,2), 'bo', 'MarkerSize', 12, 'LineWidth', 2);
        plot(blue_centroids(:,1), blue_centroids(:,2), 'b+', 'MarkerSize', 15, 'LineWidth', 2);
    end
    % Show detected red markers
    cc_red = bwconncomp(red_mask);
    if cc_red.NumObjects > 0
        stats_red = regionprops(cc_red, 'Centroid');
        red_centroids = cat(1, stats_red.Centroid);
        plot(red_centroids(:,1), red_centroids(:,2), 'ro', 'MarkerSize', 12, 'LineWidth', 2);
        plot(red_centroids(:,1), red_centroids(:,2), 'r+', 'MarkerSize', 15, 'LineWidth', 2);
    end
    % Show black circle
    if ~isnan(circle_centroids{sp}.x_pixel)
        [r, c] = find(circle_mask);
        plot(c, r, 'k.', 'MarkerSize', 2);
        plot(circle_centroids{sp}.x_pixel, circle_centroids{sp}.y_pixel, ...
             'ko', 'MarkerSize', 15, 'LineWidth', 3, 'MarkerFaceColor', 'white');
        plot(circle_centroids{sp}.x_pixel, circle_centroids{sp}.y_pixel, ...
             'k+', 'MarkerSize', 20, 'LineWidth', 3);
    end
    title('All Detected Markers');
    hold off;
    
    subplot(3,4,[11,12]);
    hold on; grid on; box on;
    % Plot blue markers
    if ~isempty(blue_data_coords)
        plot(blue_data_coords(:,1), blue_data_coords(:,2), 'bo', 'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', 'Blue Markers');
        for b = 1:size(blue_data_coords, 1)
            text(blue_data_coords(b,1) + 0.05, blue_data_coords(b,2) + 2, ...
                 sprintf('B%d: %.3f', b, blue_data_coords(b,1)), ...
                 'FontSize', 8, 'BackgroundColor', 'cyan');
        end
    end
    % Plot black circle
    if ~isnan(circle_centroids{sp}.x_data)
        plot(circle_centroids{sp}.x_data, circle_centroids{sp}.y_data, ...
             'ko', 'MarkerSize', 15, 'LineWidth', 3, 'MarkerFaceColor', 'white', 'DisplayName', 'Black Circle');
        text(circle_centroids{sp}.x_data + 0.05, circle_centroids{sp}.y_data + 2, ...
             sprintf('BC: %.3f', circle_centroids{sp}.x_data), ...
             'FontSize', 8, 'BackgroundColor', 'white');
        
        % Draw time difference lines
        if ~isempty(blue_data_coords)
            for b = 1:size(blue_data_coords, 1)
                plot([blue_data_coords(b,1), circle_centroids{sp}.x_data], ...
                     [blue_data_coords(b,2), circle_centroids{sp}.y_data], ...
                     'k--', 'LineWidth', 1);
                mid_x = (blue_data_coords(b,1) + circle_centroids{sp}.x_data) / 2;
                mid_y = (blue_data_coords(b,2) + circle_centroids{sp}.y_data) / 2;
                text(mid_x, mid_y, sprintf('Δt=%.3f', time_differences(b)), ...
                     'FontSize', 8, 'BackgroundColor', 'yellow', 'Rotation', 15);
            end
        end
    end
    xlabel('Time aSOI [ms]');
    ylabel('Distance from nozzle [mm]');
    title('Markers and Time Differences');
    xlim([0 5]); ylim([0 100]);
    legend('Location', 'best');
    hold off;
end

fprintf('\n=== All subplots processed ===\n');

%% Summary of all detected markers and time differences
fprintf('\n========================================\n');
fprintf('SUMMARY: MARKERS AND TIME DIFFERENCES\n');
fprintf('========================================\n');

for sp = 1:4
    fprintf('\nSubplot %d (%s):\n', sp, subplot_names{sp});
    
    if ~isnan(circle_centroids{sp}.x_data)
        fprintf('  Black circle: (%.3f, %.3f)\n', ...
                circle_centroids{sp}.x_data, circle_centroids{sp}.y_data);
    else
        fprintf('  Black circle: NOT DETECTED\n');
    end
    
    if isfield(detection_results{sp}, 'blue_centroids') && ~isempty(detection_results{sp}.blue_centroids)
        blue_coords = detection_results{sp}.blue_centroids;
        time_diffs = detection_results{sp}.time_differences;
        
        fprintf('  Blue markers: %d found\n', size(blue_coords, 1));
        for b = 1:size(blue_coords, 1)
            fprintf('    Marker %d: (%.3f, %.3f), Time diff: %.3f ms\n', ...
                    b, blue_coords(b,1), blue_coords(b,2), time_diffs(b));
        end
        
        if ~isempty(time_diffs)
            fprintf('  Min time difference: %.3f ms\n', min(time_diffs));
            fprintf('  Max time difference: %.3f ms\n', max(time_diffs));
        end
    else
        fprintf('  Blue markers: NONE DETECTED\n');
    end
end

fprintf('\n=== Analysis complete (no files saved) ===\n');

% Save boundaries for future use (but not the detection results)
boundary_file = fullfile(img_path, 'plot_boundaries_Patricl2022_Figure10.mat');
save(boundary_file, 'plot_boundaries', 'subplot_names');
fprintf('Boundaries saved to: %s\n', boundary_file);

%% Create summary plot showing both markers and time differences
figure('Position', [200, 200, 1000, 600]);

subplot(1,2,1);
hold on; grid on; box on;
for sp = 1:4
    % Plot black circles
    if ~isnan(circle_centroids{sp}.x_data)
        plot(circle_centroids{sp}.x_data, circle_centroids{sp}.y_data, 'ko', 'MarkerSize', 12, 'LineWidth', 2);
        text(circle_centroids{sp}.x_data + 0.05, circle_centroids{sp}.y_data + 1, ...
             sprintf('BC%d', sp), 'FontSize', 10, 'BackgroundColor', 'white');
    end
    
    % Plot blue markers
    if isfield(detection_results{sp}, 'blue_centroids') && ~isempty(detection_results{sp}.blue_centroids)
        blue_coords = detection_results{sp}.blue_centroids;
        plot(blue_coords(:,1), blue_coords(:,2), 'bo', 'MarkerSize', 10, 'LineWidth', 2);
        for b = 1:size(blue_coords, 1)
            text(blue_coords(b,1) + 0.05, blue_coords(b,2) + 1, ...
                 sprintf('B%d-%d', sp, b), 'FontSize', 8, 'BackgroundColor', 'cyan');
        end
    end
end
xlabel('Time aSOI [ms]');
ylabel('Distance from nozzle [mm]');
title('All Detected Markers');
xlim([0 5]); ylim([0 100]);
legend({'Black Circles', 'Blue Markers'}, 'Location', 'best');

subplot(1,2,2);
hold on; grid on; box on;
all_time_diffs = [];
subplot_labels = {};
for sp = 1:4
    if isfield(detection_results{sp}, 'time_differences') && ~isempty(detection_results{sp}.time_differences)
        time_diffs = detection_results{sp}.time_differences;
        x_pos = sp * ones(size(time_diffs));
        plot(x_pos, time_diffs, 'bo', 'MarkerSize', 8, 'LineWidth', 2);
        all_time_diffs = [all_time_diffs, time_diffs];
        subplot_labels{sp} = sprintf('SP%d', sp);
        
        % Add values as text
        for i = 1:length(time_diffs)
            text(sp + 0.1, time_diffs(i), sprintf('%.3f', time_diffs(i)), ...
                 'FontSize', 8, 'BackgroundColor', 'white');
        end
    else
        subplot_labels{sp} = sprintf('SP%d', sp);
    end
end
xlabel('Subplot');
ylabel('Time Difference [ms]');
title('Time Differences Between Blue and Black Markers');
xlim([0.5 4.5]);
if ~isempty(all_time_diffs)
    ylim([0 max(all_time_diffs) * 1.1]);
end
set(gca, 'XTick', 1:4, 'XTickLabel', subplot_labels);
grid on;

%% Functions
function [x_actual, y_actual] = pixel_to_data(x_pixel, y_pixel, ...
    x_min_pix, x_max_pix, y_min_pix, y_max_pix, ...
    x_min_dat, x_max_dat, y_min_dat, y_max_dat)
    x_actual = x_min_dat + (x_pixel - x_min_pix) .* ...
        (x_max_dat - x_min_dat) / (x_max_pix - x_min_pix);
    y_actual = y_max_dat - (y_pixel - y_min_pix) .* ...
        (y_max_dat - y_min_dat) / (y_max_pix - y_min_pix);
end