% Extract centroid of black circle marker from 4-panel figure
% Focus on detecting only the black circle marker centroid
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
x_max_data = 4;
y_min_data = 0;
y_max_data = 95;

%% Storage for boundaries (in FULL image coordinates)
plot_boundaries = cell(4, 1);

%% Interactive boundary selection - show FULL image
fprintf('\n========================================\n');
fprintf('INTERACTIVE BOUNDARY SELECTION\n');
fprintf('========================================\n');
fprintf('You will see the FULL 4-panel figure.\n');
fprintf('For each subplot, click 4 points to define the plot boundary:\n');
fprintf('  1. TOP-LEFT corner of the plot area\n');
fprintf('  2. TOP-RIGHT corner of the plot area\n');
fprintf('  3. BOTTOM-RIGHT corner of the plot area\n');
fprintf('  4. BOTTOM-LEFT corner of the plot area\n');
fprintf('========================================\n\n');

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
title('All Plot Boundaries Selected', 'FontSize', 14);
hold off;

fprintf('\n========================================\n');
fprintf('Boundary selection complete!\n');
fprintf('Now processing data extraction...\n');
fprintf('========================================\n\n');

%% Process each subplot with selected boundaries
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
    % The entire extracted region IS the plot area
    x_min_pixel = 1;
    x_max_pixel = size(subplot_img, 2);
    y_min_pixel = 1;
    y_max_pixel = size(subplot_img, 1);
    
    % Create plot mask (entire subplot region)
    plot_mask = true(size(subplot_img, 1), size(subplot_img, 2));
    
    %% Extract BLUE MARKER
    blue_mask = (subplot_img(:,:,3) > 150) & (subplot_img(:,:,1) < 150) & (subplot_img(:,:,2) < 150);
    blue_mask = blue_mask & plot_mask;
    blue_mask = bwareaopen(blue_mask, 5);
    
    % Find connected components
    cc_blue = bwconncomp(blue_mask);
    stats_blue = regionprops(cc_blue, 'Area', 'Centroid', 'PixelIdxList');
    
    % Extract centroids of blue markers
    blue_marker_mask = false(size(blue_mask));
    x_blue_markers = [];
    y_blue_markers = [];
    
    for i = 1:length(stats_blue)
        if stats_blue(i).Area > 5 && stats_blue(i).Area < 300
            blue_marker_mask(stats_blue(i).PixelIdxList) = true;
            x_blue_markers(end+1) = stats_blue(i).Centroid(1);
            y_blue_markers(end+1) = stats_blue(i).Centroid(2);
        end
    end
    
    fprintf('Found %d blue markers\n', length(x_blue_markers));
    
    %% Extract GREEN SOLID LINE
    green_mask = (subplot_img(:,:,2) > 150) & (subplot_img(:,:,1) < 150) & (subplot_img(:,:,3) < 150);
    green_mask = green_mask & plot_mask;
    green_mask = bwareaopen(green_mask, 10);
    
    cc_green = bwconncomp(green_mask);
    stats_green = regionprops(cc_green, 'Area', 'PixelIdxList');
    
    green_solid_mask = false(size(green_mask));
    for i = 1:length(stats_green)
        if stats_green(i).Area > 100
            green_solid_mask(stats_green(i).PixelIdxList) = true;
        end
    end
    green_solid_mask = imclose(green_solid_mask, strel('disk', 3));
    
    %% Extract BLACK DOTTED LINE - Improved method
    gray = rgb2gray(subplot_img);
    black_mask = gray < 100;
    black_mask = black_mask & plot_mask;
    
    % Remove the green line from black mask
    green_removed = imdilate(green_solid_mask, strel('disk', 5));
    black_mask = black_mask & ~green_removed;
    
    % Clean up noise
    black_mask = bwareaopen(black_mask, 3);
    
    % Separate dotted line from circle markers
    cc_black = bwconncomp(black_mask);
    stats_black = regionprops(cc_black, 'Area', 'Eccentricity', 'PixelIdxList', 'BoundingBox');
    
    black_dotted_mask = false(size(black_mask));
    
    for i = 1:length(stats_black)
        area = stats_black(i).Area;
        
        if area >= 5 && area <= 150
            if isfield(stats_black(i), 'Eccentricity') && stats_black(i).Eccentricity > 0.5
                black_dotted_mask(stats_black(i).PixelIdxList) = true;
            elseif area < 50
                black_dotted_mask(stats_black(i).PixelIdxList) = true;
            end
        end
    end
    
    % Connect nearby dots
    black_dotted_mask = imclose(black_dotted_mask, strel('line', 15, 0));
    black_dotted_mask = imclose(black_dotted_mask, strel('line', 15, 45));
    black_dotted_mask = imclose(black_dotted_mask, strel('line', 15, 90));
    black_dotted_mask = imclose(black_dotted_mask, strel('disk', 3));
    
    %% Extract BLACK CIRCLE (marker)
    black_mask_full = gray < 100;
    black_mask_full = black_mask_full & plot_mask;
    black_mask_full = black_mask_full & ~green_removed;
    
    % Remove the dotted line
    dotted_expanded = imdilate(black_dotted_mask, strel('disk', 3));
    black_mask_full = black_mask_full & ~dotted_expanded;
    
    cc_circles = bwconncomp(black_mask_full);
    stats_circles = regionprops(cc_circles, 'Area', 'Circularity', 'Eccentricity', 'PixelIdxList');
    
    black_circle_mask = false(size(black_mask));
    for i = 1:length(stats_circles)
        area = stats_circles(i).Area;
        
        if area > 50 && area < 500
            if isfield(stats_circles(i), 'Eccentricity')
                if stats_circles(i).Eccentricity < 0.7
                    black_circle_mask(stats_circles(i).PixelIdxList) = true;
                end
            else
                black_circle_mask(stats_circles(i).PixelIdxList) = true;
            end
        end
    end
    
    %% Extract data from masks
    [x_green_px, y_green_px] = extract_line_raw(green_solid_mask);
    [x_black_dot_px, y_black_dot_px] = extract_line_raw(black_dotted_mask);
    [x_circle_px, y_circle_px] = extract_line_raw(black_circle_mask);
    
    % Blue markers - use centroids directly
    x_blue_px = x_blue_markers';
    y_blue_px = y_blue_markers';
    
    %% Convert to data coordinates
    [x_green_data, y_green_data] = pixel_to_data(x_green_px, y_green_px, ...
        x_min_pixel, x_max_pixel, y_min_pixel, y_max_pixel, ...
        x_min_data, x_max_data, y_min_data, y_max_data);
    
    [x_black_dotted_data, y_black_dotted_data] = pixel_to_data(x_black_dot_px, y_black_dot_px, ...
        x_min_pixel, x_max_pixel, y_min_pixel, y_max_pixel, ...
        x_min_data, x_max_data, y_min_data, y_max_data);
    
    [x_circle_data, y_circle_data] = pixel_to_data(x_circle_px, y_circle_px, ...
        x_min_pixel, x_max_pixel, y_min_pixel, y_max_pixel, ...
        x_min_data, x_max_data, y_min_data, y_max_data);
    
    [x_blue_data, y_blue_data] = pixel_to_data(x_blue_px, y_blue_px, ...
        x_min_pixel, x_max_pixel, y_min_pixel, y_max_pixel, ...
        x_min_data, x_max_data, y_min_data, y_max_data);
    
    %% Save to MAT file
    output_file = fullfile(img_path, [subplot_names{sp} '_Patrick2022_Figure13_extracted_data.mat']);
    save(output_file, 'x_green_data', 'y_green_data', ...
        'x_black_dotted_data', 'y_black_dotted_data', ...
        'x_circle_data', 'y_circle_data', ...
        'x_blue_data', 'y_blue_data');
    
    fprintf('Saved: %s\n', output_file);
    fprintf('  Green line: %d points\n', length(x_green_data));
    fprintf('  Black dotted: %d points\n', length(x_black_dotted_data));
    fprintf('  Circle: %d points\n', length(x_circle_data));
    fprintf('  Blue markers: %d points\n', length(x_blue_data));
    
    %% Visualization for this subplot
    figure('Position', [100 + sp*50, 100, 1200, 800]);
    
    subplot(2,4,1);
    imshow(subplot_img);
    title(['Original: ' subplot_names{sp}], 'Interpreter', 'none');
    
    subplot(2,4,2);
    imshow(blue_marker_mask);
    title('Blue Marker Mask');
    
    subplot(2,4,3);
    imshow(green_solid_mask);
    title('Green Solid Line Mask');
    
    subplot(2,4,4);
    imshow(black_dotted_mask);
    title('Black Dotted Line Mask');
    
    subplot(2,4,5);
    imshow(black_circle_mask);
    title('Circle Marker Mask');
    
    subplot(2,4,6);
    % Show all masks overlaid
    imshow(subplot_img);
    hold on;
    [r, c] = find(blue_marker_mask);
    plot(c, r, 'b.', 'MarkerSize', 3);
    [r, c] = find(green_solid_mask);
    plot(c, r, 'g.', 'MarkerSize', 2);
    [r, c] = find(black_dotted_mask);
    plot(c, r, 'k.', 'MarkerSize', 2);
    [r, c] = find(black_circle_mask);
    plot(c, r, 'r.', 'MarkerSize', 3);
    title('All Masks Overlay');
    hold off;
    
    subplot(2,4,[7,8]);
    hold on; grid on; box on;
    plot(x_green_data, y_green_data, 'g-', 'LineWidth', 2, 'DisplayName', 'Green solid');
    plot(x_black_dotted_data, y_black_dotted_data, 'k:', 'LineWidth', 2, 'DisplayName', 'Black dotted');
    plot(x_circle_data, y_circle_data, 'ko', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Circle marker');
    plot(x_blue_data, y_blue_data, 'b*', 'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', 'Blue marker');
    xlabel('Time aSOI [ms]');
    ylabel('Distance from nozzle [mm]');
    title('Extracted Data');
    xlim([0 5]); ylim([0 100]);
    legend('Location', 'best');
    hold off;
end

fprintf('\n=== All subplots processed ===\n');

% Save boundaries for future use
boundary_file = fullfile(img_path, 'plot_boundaries_Patrick2022_Figure13.mat');
save(boundary_file, 'plot_boundaries', 'subplot_names');
fprintf('\nBoundaries saved to: %s\n', boundary_file);
fprintf('You can load these boundaries next time to skip manual selection.\n');

%% Functions
function [x_data, y_data] = extract_line_raw(mask)
    [rows, cols] = find(mask);
    if isempty(rows)
        x_data = [];
        y_data = [];
        return;
    end
    x_unique = unique(cols);
    y_median = zeros(size(x_unique));
    for i = 1:length(x_unique)
        y_median(i) = median(rows(cols == x_unique(i)));
    end
    [x_data, sort_idx] = sort(x_unique);
    y_data = y_median(sort_idx);
end

function [x_actual, y_actual] = pixel_to_data(x_pixel, y_pixel, ...
    x_min_pix, x_max_pix, y_min_pix, y_max_pix, ...
    x_min_dat, x_max_dat, y_min_dat, y_max_dat)
    x_actual = x_min_dat + (x_pixel - x_min_pix) .* ...
        (x_max_dat - x_min_dat) / (x_max_pix - x_min_pix);
    y_actual = y_max_dat - (y_pixel - y_min_pix) .* ...
        (y_max_dat - y_min_dat) / (y_max_pix - y_min_pix);
end