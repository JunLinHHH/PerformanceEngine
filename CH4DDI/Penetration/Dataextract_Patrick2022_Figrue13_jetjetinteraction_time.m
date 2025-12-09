% Extract centroid of blue and red markers from 4-panel figure
% Focus on detecting only blue and red markers and calculating time differences
clear; clc; close all;

%% Configuration - UPDATE THIS PATH
img_path_full = 'C:\Users\jun-y\OneDrive - UNSW\CH4DDI\CH4DDI vs H2DDI CVCC\Patrick2022_handled\Patrick2022_Figure13.png';

% Subplot names (will be used for saving)
subplot_names = {'H-0.07ms-D', 'H-1.07ms-D', 'H-2.07ms-D', 'H-3.07ms-D'};

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
boundary_file = fullfile(img_path, 'plot_boundaries_Patrick2022_Figure13.mat');
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
    fprintf('BLUE AND RED MARKER DETECTION\n');
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
fprintf('Now processing blue and red marker detection...\n');
fprintf('========================================\n\n');

%% Process each subplot with selected boundaries - DETECT BLUE AND RED MARKERS
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
    
    %% STEP 1: DETECT BLUE MARKERS
    fprintf('Step 1: Detecting blue markers...\n');
    
    % Extract blue markers - focus on color difference rather than absolute values
    % For high-contrast images, look for pixels where blue is significantly higher
    blue_diff = double(subplot_img(:,:,3)) - double(subplot_img(:,:,1)); % Blue - Red
    blue_diff2 = double(subplot_img(:,:,3)) - double(subplot_img(:,:,2)); % Blue - Green
    
    blue_mask = (subplot_img(:,:,3) > 100) & ...           % Some blue content
                (blue_diff > 50) & ...                     % Blue much higher than red
                (blue_diff2 > 30) & ...                    % Blue higher than green
                (subplot_img(:,:,1) < 200) & ...           % Not pure white
                (subplot_img(:,:,2) < 200);                % Not pure white
    blue_mask = blue_mask & plot_mask;
    
    % Remove small noise and keep only substantial regions
    blue_mask = bwareaopen(blue_mask, 5);
    
    % Additional morphological operations to clean up
    se = strel('disk', 1);
    blue_mask = imclose(blue_mask, se);
    blue_mask = imfill(blue_mask, 'holes');
    
    % Find blue marker centroids with shape criteria
    cc_blue = bwconncomp(blue_mask);
    stats_blue = regionprops(cc_blue, 'Area', 'Centroid', 'PixelIdxList', 'Solidity', 'Eccentricity', 'BoundingBox');
    
    blue_centroids = [];
    
    for i = 1:length(stats_blue)
        % Criteria to ensure we get actual markers:
        % - Area between 8 and 1000 pixels (very flexible for different marker sizes)
        % - Reasonable solidity (not too stringent)
        % - Not too elongated
        % - Aspect ratio check for circular shapes
        bbox = stats_blue(i).BoundingBox;
        aspect_ratio = bbox(4) / bbox(3); % height/width
        
        if stats_blue(i).Area > 8 && stats_blue(i).Area < 1000 && ...
           stats_blue(i).Solidity > 0.3 && stats_blue(i).Eccentricity < 0.95 && ...
           aspect_ratio > 0.2 && aspect_ratio < 5.0
            centroid = stats_blue(i).Centroid;
            blue_centroids = [blue_centroids; centroid];
        end
    end
    
    fprintf('Found %d blue markers\n', size(blue_centroids, 1));
    
    % Debug: Show some color statistics
    if size(blue_centroids, 1) == 0
        fprintf('Debug - Blue detection: Checking color ranges in image...\n');
        fprintf('  Blue channel: min=%.1f, max=%.1f, mean=%.1f\n', ...
                min(subplot_img(:,:,3),[],'all'), max(subplot_img(:,:,3),[],'all'), mean(subplot_img(:,:,3),'all'));
        fprintf('  Red channel: min=%.1f, max=%.1f, mean=%.1f\n', ...
                min(subplot_img(:,:,1),[],'all'), max(subplot_img(:,:,1),[],'all'), mean(subplot_img(:,:,1),'all'));
        fprintf('  Green channel: min=%.1f, max=%.1f, mean=%.1f\n', ...
                min(subplot_img(:,:,2),[],'all'), max(subplot_img(:,:,2),[],'all'), mean(subplot_img(:,:,2),'all'));
        fprintf('  Blue mask pixels: %d\n', sum(blue_mask(:)));
    end
    
    %% STEP 2: DETECT RED MARKERS
    fprintf('Step 2: Detecting red markers...\n');
    
    % Extract red markers - focus on color difference rather than absolute values
    % For high-contrast images, look for pixels where red is significantly higher
    red_diff = double(subplot_img(:,:,1)) - double(subplot_img(:,:,2)); % Red - Green
    red_diff2 = double(subplot_img(:,:,1)) - double(subplot_img(:,:,3)); % Red - Blue
    
    red_mask = (subplot_img(:,:,1) > 100) & ...            % Some red content
               (red_diff > 30) & ...                       % Red higher than green
               (red_diff2 > 30) & ...                      % Red higher than blue
               (subplot_img(:,:,2) < 200) & ...            % Not pure white
               (subplot_img(:,:,3) < 200);                 % Not pure white
    red_mask = red_mask & plot_mask;
    
    % Remove small noise and keep only substantial regions
    red_mask = bwareaopen(red_mask, 5);
    
    % Additional morphological operations to clean up and ensure we get compact regions
    se = strel('disk', 1);
    red_mask = imclose(red_mask, se);
    red_mask = imfill(red_mask, 'holes');
    
    % Find red marker centroids with shape criteria
    cc_red = bwconncomp(red_mask);
    stats_red = regionprops(cc_red, 'Area', 'Centroid', 'PixelIdxList', 'Solidity', 'Eccentricity', 'BoundingBox');
    
    red_centroids = [];
    
    for i = 1:length(stats_red)
        % Criteria to avoid dashed lines while detecting markers:
        % - Area between 8 and 1000 pixels (very flexible for different marker sizes)
        % - Reasonable solidity (compact shape, not elongated lines)
        % - Not too elongated (eccentricity)
        % - Aspect ratio check for circular shapes
        bbox = stats_red(i).BoundingBox;
        aspect_ratio = bbox(4) / bbox(3); % height/width
        
        if stats_red(i).Area > 8 && stats_red(i).Area < 1000 && ...
           stats_red(i).Solidity > 0.3 && stats_red(i).Eccentricity < 0.95 && ...
           aspect_ratio > 0.2 && aspect_ratio < 5.0
            centroid = stats_red(i).Centroid;
            red_centroids = [red_centroids; centroid];
        end
    end
    
    fprintf('Found %d red markers\n', size(red_centroids, 1));
    
    % Debug: Show some color statistics for red detection
    if size(red_centroids, 1) == 0
        fprintf('Debug - Red detection: Checking color ranges and shapes...\n');
        fprintf('  Red mask pixels: %d\n', sum(red_mask(:)));
        fprintf('  Connected components found: %d\n', cc_red.NumObjects);
        if cc_red.NumObjects > 0
            for j = 1:min(3, cc_red.NumObjects) % Show first 3 components
                bbox = stats_red(j).BoundingBox;
                aspect_ratio = bbox(4) / bbox(3);
                fprintf('    Component %d: Area=%.1f, Solidity=%.2f, Eccentricity=%.2f, AspectRatio=%.2f\n', ...
                        j, stats_red(j).Area, stats_red(j).Solidity, stats_red(j).Eccentricity, aspect_ratio);
            end
        end
        % Show a sample of red pixels if any exist
        red_pixels = subplot_img(:,:,1) > 100;
        fprintf('  Potential red pixels (R>100): %d\n', sum(red_pixels(:)));
    end
    
    %% Convert pixel coordinates to data coordinates
    blue_data_coords = [];
    red_data_coords = [];
    
    if ~isempty(blue_centroids)
        for i = 1:size(blue_centroids, 1)
            [x_data, y_data] = pixel_to_data(blue_centroids(i,1), blue_centroids(i,2), ...
                x_min_pixel, x_max_pixel, y_min_pixel, y_max_pixel, ...
                x_min_data, x_max_data, y_min_data, y_max_data);
            blue_data_coords = [blue_data_coords; x_data, y_data];
        end
    end
    
    if ~isempty(red_centroids)
        for i = 1:size(red_centroids, 1)
            [x_data, y_data] = pixel_to_data(red_centroids(i,1), red_centroids(i,2), ...
                x_min_pixel, x_max_pixel, y_min_pixel, y_max_pixel, ...
                x_min_data, x_max_data, y_min_data, y_max_data);
            red_data_coords = [red_data_coords; x_data, y_data];
        end
    end
    
    %% Calculate time differences between blue and red markers
    time_differences = [];
    marker_pairs = [];
    
    if ~isempty(blue_data_coords) && ~isempty(red_data_coords)
        % Calculate all pairwise time differences
        for b = 1:size(blue_data_coords, 1)
            for r = 1:size(red_data_coords, 1)
                time_diff = abs(blue_data_coords(b,1) - red_data_coords(r,1));
                time_differences = [time_differences; time_diff];
                marker_pairs = [marker_pairs; b, r];
            end
        end
        
        fprintf('Calculated %d time differences between blue and red markers\n', length(time_differences));
        if ~isempty(time_differences)
            fprintf('Min time difference: %.3f ms\n', min(time_differences));
            fprintf('Max time difference: %.3f ms\n', max(time_differences));
            fprintf('Mean time difference: %.3f ms\n', mean(time_differences));
        end
    else
        fprintf('Cannot calculate time differences - missing blue or red markers\n');
    end
    
    %% Store results
    detection_results{sp}.blue_centroids_pixel = blue_centroids;
    detection_results{sp}.red_centroids_pixel = red_centroids;
    detection_results{sp}.blue_centroids_data = blue_data_coords;
    detection_results{sp}.red_centroids_data = red_data_coords;
    detection_results{sp}.time_differences = time_differences;
    detection_results{sp}.marker_pairs = marker_pairs;
    
    %% Visualization for this subplot
    figure('Position', [100 + sp*50, 100 + sp*50, 1400, 1000]);
    sgtitle(sprintf('Subplot %d: %s - Blue and Red Marker Detection', sp, subplot_names{sp}), ...
            'FontSize', 16, 'FontWeight', 'bold');
    
    % Original image
    subplot(3,4,1);
    imshow(subplot_img);
    title('Original Subplot');
    
    % Blue mask
    subplot(3,4,2);
    imshow(blue_mask);
    title('Blue Marker Mask');
    
    % Red mask
    subplot(3,4,3);
    imshow(red_mask);
    title('Red Marker Mask');
    
    % Combined visualization
    subplot(3,4,4);
    imshow(subplot_img);
    hold on;
    % Show blue markers
    if ~isempty(blue_centroids)
        plot(blue_centroids(:,1), blue_centroids(:,2), 'bo', 'MarkerSize', 12, 'LineWidth', 2);
        plot(blue_centroids(:,1), blue_centroids(:,2), 'b+', 'MarkerSize', 15, 'LineWidth', 2);
    end
    % Show red markers
    if ~isempty(red_centroids)
        plot(red_centroids(:,1), red_centroids(:,2), 'ro', 'MarkerSize', 12, 'LineWidth', 2);
        plot(red_centroids(:,1), red_centroids(:,2), 'r+', 'MarkerSize', 15, 'LineWidth', 2);
    end
    title('Detected Markers');
    hold off;
    
    % Data coordinate plot
    subplot(3,4,[5,6,7,8]);
    hold on; grid on; box on;
    
    % Plot blue markers
    if ~isempty(blue_data_coords)
        plot(blue_data_coords(:,1), blue_data_coords(:,2), 'bo', 'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', 'Blue Markers');
        for b = 1:size(blue_data_coords, 1)
            text(blue_data_coords(b,1) + 0.05, blue_data_coords(b,2) + 2, ...
                 sprintf('B%d: %.3f ms', b, blue_data_coords(b,1)), ...
                 'FontSize', 8, 'BackgroundColor', 'cyan');
        end
    end
    
    % Plot red markers
    if ~isempty(red_data_coords)
        plot(red_data_coords(:,1), red_data_coords(:,2), 'ro', 'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', 'Red Markers');
        for r = 1:size(red_data_coords, 1)
            text(red_data_coords(r,1) + 0.05, red_data_coords(r,2) + 2, ...
                 sprintf('R%d: %.3f ms', r, red_data_coords(r,1)), ...
                 'FontSize', 8, 'BackgroundColor', 'white');
        end
    end
    
    % Draw connection lines and time differences
    if ~isempty(time_differences) && ~isempty(marker_pairs)
        for i = 1:length(time_differences)
            b_idx = marker_pairs(i,1);
            r_idx = marker_pairs(i,2);
            plot([blue_data_coords(b_idx,1), red_data_coords(r_idx,1)], ...
                 [blue_data_coords(b_idx,2), red_data_coords(r_idx,2)], ...
                 'k--', 'LineWidth', 1);
            mid_x = (blue_data_coords(b_idx,1) + red_data_coords(r_idx,1)) / 2;
            mid_y = (blue_data_coords(b_idx,2) + red_data_coords(r_idx,2)) / 2;
            text(mid_x, mid_y, sprintf('Δt=%.3f ms', time_differences(i)), ...
                 'FontSize', 8, 'BackgroundColor', 'yellow', 'Rotation', 15);
        end
    end
    
    xlabel('Time aSOI [ms]');
    ylabel('Distance from nozzle [mm]');
    title('Markers and Time Differences in Data Coordinates');
    xlim([0 5]); ylim([0 100]);
    legend('Location', 'best');
    hold off;
    
    % Time difference histogram
    subplot(3,4,[9,10]);
    if ~isempty(time_differences)
        histogram(time_differences, 'BinWidth', 0.1, 'FaceColor', 'yellow', 'EdgeColor', 'black');
        xlabel('Time Difference [ms]');
        ylabel('Count');
        title('Distribution of Time Differences');
        grid on;
    else
        text(0.5, 0.5, 'No time differences calculated', 'HorizontalAlignment', 'center');
        title('No Time Differences');
    end
    
    % Statistics
    subplot(3,4,[11,12]);
    axis off;
    stats_text = {};
    stats_text{end+1} = sprintf('Blue markers detected: %d', size(blue_centroids, 1));
    stats_text{end+1} = sprintf('Red markers detected: %d', size(red_centroids, 1));
    stats_text{end+1} = sprintf('Time differences calculated: %d', length(time_differences));
    if ~isempty(time_differences)
        stats_text{end+1} = sprintf('Min time diff: %.3f ms', min(time_differences));
        stats_text{end+1} = sprintf('Max time diff: %.3f ms', max(time_differences));
        stats_text{end+1} = sprintf('Mean time diff: %.3f ms', mean(time_differences));
        stats_text{end+1} = sprintf('Std time diff: %.3f ms', std(time_differences));
    end
    
    text(0.1, 0.9, stats_text, 'FontSize', 12, 'VerticalAlignment', 'top', ...
         'FontWeight', 'bold', 'BackgroundColor', 'white');
    title('Detection Statistics', 'FontWeight', 'bold');
end

fprintf('\n=== All subplots processed ===\n');

%% Summary of all detected markers and time differences
fprintf('\n========================================\n');
fprintf('SUMMARY: BLUE AND RED MARKERS WITH TIME DIFFERENCES\n');
fprintf('========================================\n');

for sp = 1:4
    fprintf('\nSubplot %d (%s):\n', sp, subplot_names{sp});
    
    blue_coords = detection_results{sp}.blue_centroids_data;
    red_coords = detection_results{sp}.red_centroids_data;
    time_diffs = detection_results{sp}.time_differences;
    marker_pairs = detection_results{sp}.marker_pairs;
    
    if ~isempty(blue_coords)
        fprintf('  Blue markers: %d found\n', size(blue_coords, 1));
        for b = 1:size(blue_coords, 1)
            fprintf('    Blue %d: (%.3f ms, %.3f mm)\n', ...
                    b, blue_coords(b,1), blue_coords(b,2));
        end
    else
        fprintf('  Blue markers: NONE DETECTED\n');
    end
    
    if ~isempty(red_coords)
        fprintf('  Red markers: %d found\n', size(red_coords, 1));
        for r = 1:size(red_coords, 1)
            fprintf('    Red %d: (%.3f ms, %.3f mm)\n', ...
                    r, red_coords(r,1), red_coords(r,2));
        end
    else
        fprintf('  Red markers: NONE DETECTED\n');
    end
    
    if ~isempty(time_diffs)
        fprintf('  Time differences: %d calculated\n', length(time_diffs));
        for i = 1:length(time_diffs)
            b_idx = marker_pairs(i,1);
            r_idx = marker_pairs(i,2);
            fprintf('    Blue%d-Red%d: %.3f ms\n', b_idx, r_idx, time_diffs(i));
        end
        fprintf('  Min: %.3f ms, Max: %.3f ms, Mean: %.3f ms\n', ...
                min(time_diffs), max(time_diffs), mean(time_diffs));
    else
        fprintf('  Time differences: NONE CALCULATED\n');
    end
end

fprintf('\n=== Analysis complete ===\n');

% Save boundaries for future use
boundary_file = fullfile(img_path, 'plot_boundaries_Patrick2022_Figure13.mat');
save(boundary_file, 'plot_boundaries', 'subplot_names');
fprintf('Boundaries saved to: %s\n', boundary_file);

%% Create comprehensive summary plot
figure('Position', [200, 200, 1400, 800]);
sgtitle('Summary: Blue and Red Marker Detection with Time Differences', 'FontSize', 16, 'FontWeight', 'bold');

% All markers plot
subplot(2,3,1);
hold on; grid on; box on;
colors = {'r', 'g', 'b', 'm'};
for sp = 1:4
    % Plot blue markers
    blue_coords = detection_results{sp}.blue_centroids_data;
    if ~isempty(blue_coords)
        plot(blue_coords(:,1), blue_coords(:,2), 'bo', 'MarkerSize', 10, 'LineWidth', 2);
        for b = 1:size(blue_coords, 1)
            text(blue_coords(b,1) + 0.05, blue_coords(b,2) + 1, ...
                 sprintf('B%d-%d', sp, b), 'FontSize', 8, 'BackgroundColor', 'cyan');
        end
    end
    
    % Plot red markers
    red_coords = detection_results{sp}.red_centroids_data;
    if ~isempty(red_coords)
        plot(red_coords(:,1), red_coords(:,2), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
        for r = 1:size(red_coords, 1)
            text(red_coords(r,1) + 0.05, red_coords(r,2) + 1, ...
                 sprintf('R%d-%d', sp, r), 'FontSize', 8, 'BackgroundColor', 'white');
        end
    end
end
xlabel('Time aSOI [ms]');
ylabel('Distance from nozzle [mm]');
title('All Detected Markers');
xlim([0 5]); ylim([0 100]);
legend({'Blue Markers', 'Red Markers'}, 'Location', 'best');

% Time differences by subplot
subplot(2,3,2);
hold on; grid on; box on;
all_time_diffs = [];
for sp = 1:4
    time_diffs = detection_results{sp}.time_differences;
    if ~isempty(time_diffs)
        % Ensure time_diffs is a row vector for consistent concatenation
        time_diffs = time_diffs(:)'; % Convert to row vector
        x_pos = sp * ones(size(time_diffs));
        plot(x_pos, time_diffs, 'o', 'Color', colors{sp}, 'MarkerSize', 8, 'LineWidth', 2);
        all_time_diffs = [all_time_diffs, time_diffs]; % Now both are row vectors
        
        % Add values as text
        for i = 1:length(time_diffs)
            text(sp + 0.1, time_diffs(i), sprintf('%.3f', time_diffs(i)), ...
                 'FontSize', 8, 'BackgroundColor', 'white');
        end
    end
end
xlabel('Subplot');
ylabel('Time Difference [ms]');
title('Time Differences by Subplot');
xlim([0.5 4.5]);
if ~isempty(all_time_diffs)
    ylim([0 max(all_time_diffs) * 1.1]);
end
set(gca, 'XTick', 1:4, 'XTickLabel', subplot_names);

% Overall time difference histogram
subplot(2,3,3);
if ~isempty(all_time_diffs)
    histogram(all_time_diffs, 'BinWidth', 0.1, 'FaceColor', 'yellow', 'EdgeColor', 'black');
    xlabel('Time Difference [ms]');
    ylabel('Count');
    title('Overall Time Difference Distribution');
    grid on;
else
    text(0.5, 0.5, 'No time differences calculated', 'HorizontalAlignment', 'center');
    title('No Time Differences');
end

% Summary statistics
subplot(2,3,[4,5,6]);
axis off;
summary_text = {};
summary_text{end+1} = 'DETECTION SUMMARY:';
summary_text{end+1} = '';

total_blue = 0;
total_red = 0;
total_time_diffs = 0;

for sp = 1:4
    blue_count = size(detection_results{sp}.blue_centroids_data, 1);
    red_count = size(detection_results{sp}.red_centroids_data, 1);
    time_diff_count = length(detection_results{sp}.time_differences);
    
    total_blue = total_blue + blue_count;
    total_red = total_red + red_count;
    total_time_diffs = total_time_diffs + time_diff_count;
    
    summary_text{end+1} = sprintf('Subplot %d (%s):', sp, subplot_names{sp});
    summary_text{end+1} = sprintf('  Blue: %d, Red: %d, Time diffs: %d', blue_count, red_count, time_diff_count);
end

summary_text{end+1} = '';
summary_text{end+1} = sprintf('TOTALS: Blue: %d, Red: %d, Time diffs: %d', total_blue, total_red, total_time_diffs);

if ~isempty(all_time_diffs)
    summary_text{end+1} = '';
    summary_text{end+1} = 'TIME DIFFERENCE STATISTICS:';
    summary_text{end+1} = sprintf('  Min: %.3f ms', min(all_time_diffs));
    summary_text{end+1} = sprintf('  Max: %.3f ms', max(all_time_diffs));
    summary_text{end+1} = sprintf('  Mean: %.3f ms', mean(all_time_diffs));
    summary_text{end+1} = sprintf('  Std: %.3f ms', std(all_time_diffs));
end

text(0.1, 0.9, summary_text, 'FontSize', 12, 'VerticalAlignment', 'top', ...
     'FontWeight', 'bold', 'BackgroundColor', 'white');
title('Summary Statistics', 'FontWeight', 'bold', 'FontSize', 14);

%% Functions
function [x_actual, y_actual] = pixel_to_data(x_pixel, y_pixel, ...
    x_min_pix, x_max_pix, y_min_pix, y_max_pix, ...
    x_min_dat, x_max_dat, y_min_dat, y_max_dat)
    x_actual = x_min_dat + (x_pixel - x_min_pix) .* ...
        (x_max_dat - x_min_dat) / (x_max_pix - x_min_pix);
    y_actual = y_max_dat - (y_pixel - y_min_pix) .* ...
        (y_max_dat - y_min_dat) / (y_max_pix - y_min_pix);
end

