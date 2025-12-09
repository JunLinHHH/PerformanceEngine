% Complete data extraction from plot image - Older MATLAB Compatible
% Step 1: Find plot boundary, then extract data only within that region
% All functions are at the end of the script

clear; clc; close all;

% Read the image
img = imread('C:\Users\jun-y\OneDrive - UNSW\CH4DDI\CH4 vs H2\Penetration\Source\H2_1060k_200bar_Paul.png');
[height, width, ~] = size(img);

% Extract the filename without extension for saving
[img_path, img_name, ~] = fileparts('C:\Users\jun-y\OneDrive - UNSW\CH4DDI\CH4 vs H2\Penetration\Source\H2_1060k_200bar_Paul.png');
output_filename = fullfile(img_path, [img_name '_extracted_data.mat']);

%% STEP 0: CALIBRATION HELPER (Run this first to find pixel coordinates)
% Uncomment the lines below to find the exact pixel coordinates
% figure; imshow(img); 
% title('Click on plot corners: 1) Bottom-left, 2) Bottom-right, 3) Top-left, 4) Top-right');
% [x_cal, y_cal] = ginput(4);
% fprintf('Calibration coordinates:\n');
% fprintf('x_min_pixel = %.0f;  %% Bottom-left and top-left x\n', min(x_cal));
% fprintf('x_max_pixel = %.0f;  %% Bottom-right and top-right x\n', max(x_cal));
% fprintf('y_max_pixel = %.0f;  %% Bottom-left and bottom-right y\n', max(y_cal));
% fprintf('y_min_pixel = %.0f;  %% Top-left and top-right y\n', min(y_cal));

%% STEP 1: Detect plot boundary (rectangular frame)
fprintf('=== STEP 1: Detecting plot boundary ===\n');

% Method 1: Try automatic detection based on the plot frame
gray = rgb2gray(img);
% Look for the black frame of the plot
frame_mask = gray < 50;  % Dark pixels (black frame)

% Find the largest rectangular region
[rows, cols] = find(~frame_mask);  % Find non-frame pixels (inside plot)
if ~isempty(rows) && ~isempty(cols)
    x_min_pixel = min(cols);
    x_max_pixel = max(cols);
    y_min_pixel = min(rows);
    y_max_pixel = max(rows);
    
    % Validate the detection
    plot_width = x_max_pixel - x_min_pixel;
    plot_height = y_max_pixel - y_min_pixel;
    
    % If detection seems reasonable, use it
    if plot_width > 200 && plot_height > 200 && plot_width < width && plot_height < height
        fprintf('Automatic boundary detection successful!\n');
    else
        fprintf('Automatic detection gave unreasonable values, using manual defaults.\n');
        x_min_pixel = 113;
        x_max_pixel = 336;
        y_min_pixel = 66;
        y_max_pixel = 287;
    end
else
    fprintf('Automatic detection failed, using manual defaults.\n');
    x_min_pixel = 113;
    x_max_pixel = 336;
    y_min_pixel = 66;
    y_max_pixel = 287;
end

fprintf('Plot boundary:\n');
fprintf('  x_min_pixel = %.0f\n', x_min_pixel);
fprintf('  x_max_pixel = %.0f\n', x_max_pixel);
fprintf('  y_min_pixel = %.0f\n', y_min_pixel);
fprintf('  y_max_pixel = %.0f\n', y_max_pixel);
fprintf('  Plot width: %.0f pixels, Plot height: %.0f pixels\n', x_max_pixel-x_min_pixel, y_max_pixel-y_min_pixel);

% Create plot boundary mask
plot_mask = false(height, width);
plot_mask(round(y_min_pixel):round(y_max_pixel), round(x_min_pixel):round(x_max_pixel)) = true;

% Visualize the detected boundary
figure('Position', [100, 100, 1000, 400]);
subplot(1,2,1);
imshow(img);
hold on;
rectangle('Position', [x_min_pixel, y_min_pixel, x_max_pixel-x_min_pixel, y_max_pixel-y_min_pixel], ...
    'EdgeColor', 'g', 'LineWidth', 3);
title('Detected Plot Boundary', 'FontSize', 14, 'FontWeight', 'bold');
hold off;

subplot(1,2,2);
imshow(plot_mask);
title('Plot Region Mask', 'FontSize', 14, 'FontWeight', 'bold');

%% STEP 2: Manual calibration for axis values
% Actual data ranges
x_min_data = 0;
x_max_data = 7;
y_min_data = 0;
y_max_data = 80;

fprintf('\n=== Axis Calibration ===\n');
fprintf('X-axis: %d to %d (positions)\n', x_min_data, x_max_data);
fprintf('Y-axis: %d to %d mm (penetration)\n', y_min_data, y_max_data);

%% STEP 3: Extract each line by color WITHIN plot boundary only

% RED LINES - extract all red within plot boundary
red_mask_all = (img(:,:,1) > 180) & (img(:,:,2) < 100) & (img(:,:,3) < 100);
red_mask_all = red_mask_all & plot_mask;  % Apply plot boundary mask
red_mask_all = bwareaopen(red_mask_all, 5);

% Separate red lines by pattern:
% 1. Solid line (large continuous areas)
% 2. Dash-dotted line (medium segments) - ONLY detect in LEFT HALF (x <= 3)
% 3. Dots only (small isolated points)

% Create a mask for left half of plot (x <= 3 in data space, which is x <= midpoint in pixels)
x_midpoint_pixel = x_min_pixel + (x_max_pixel - x_min_pixel) * (3.0 / 7.0);
left_half_mask = false(size(red_mask_all));
left_half_mask(:, 1:round(x_midpoint_pixel)) = true;

cc_red = bwconncomp(red_mask_all);
stats_red = regionprops(cc_red, 'Area', 'PixelIdxList', 'Centroid');

red_solid_mask = false(size(red_mask_all));
red_dashdot_mask = false(size(red_mask_all));
red_dots_mask = false(size(red_mask_all));

for i = 1:length(stats_red)
    area = stats_red(i).Area;
    centroid_x = stats_red(i).Centroid(1);  % x-coordinate of centroid
    
    % Check if component is in left half for dash-dot detection
    in_left_half = centroid_x <= x_midpoint_pixel;
    
    % Solid line: very large connected components
    if area > 800
        red_solid_mask(stats_red(i).PixelIdxList) = true;
    % Dash-dotted line: medium segments AND in left half only
    elseif area > 30 && area <= 800 && in_left_half
        red_dashdot_mask(stats_red(i).PixelIdxList) = true;
    % Dots only: small isolated points OR medium segments in right half
    else
        red_dots_mask(stats_red(i).PixelIdxList) = true;
    end
end

% Connect nearby segments for continuous lines
red_solid_mask = imclose(red_solid_mask, strel('disk', 4));
red_dashdot_mask = imclose(red_dashdot_mask, strel('disk', 3));
% Don't connect dots - keep them as individual points
red_dots_mask = bwareaopen(red_dots_mask, 3);

% BLUE LINES - solid line and dots only, within plot boundary
blue_mask_all = (img(:,:,3) > 180) & (img(:,:,1) < 100) & (img(:,:,2) < 100);
blue_mask_all = blue_mask_all & plot_mask;  % Apply plot boundary mask
blue_mask_all = bwareaopen(blue_mask_all, 5);

cc_blue = bwconncomp(blue_mask_all);
stats_blue = regionprops(cc_blue, 'Area', 'PixelIdxList');

blue_solid_mask = false(size(blue_mask_all));
blue_dots_mask = false(size(blue_mask_all));

for i = 1:length(stats_blue)
    area = stats_blue(i).Area;
    % Solid line: large connected areas
    if area > 800
        blue_solid_mask(stats_blue(i).PixelIdxList) = true;
    % Dots: small isolated points
    else
        blue_dots_mask(stats_blue(i).PixelIdxList) = true;
    end
end

blue_solid_mask = imclose(blue_solid_mask, strel('disk', 4));
blue_dots_mask = bwareaopen(blue_dots_mask, 3);

% BLACK LINES - solid line and dots only, within plot boundary
% More lenient threshold for black detection
black_mask_all = (img(:,:,1) < 80) & (img(:,:,2) < 80) & (img(:,:,3) < 80);
black_mask_all = black_mask_all & plot_mask;  % Apply plot boundary mask

% Remove legend area from black mask
legend_x_start = round(x_min_pixel + (x_max_pixel - x_min_pixel) * 0.6);
legend_y_start = round(y_min_pixel + (y_max_pixel - y_min_pixel) * 0.1);
legend_y_end = round(y_min_pixel + (y_max_pixel - y_min_pixel) * 0.4);
legend_x_end = round(x_max_pixel - 10);
black_mask_all(legend_y_start:legend_y_end, legend_x_start:legend_x_end) = false;

% CRITICAL: Remove bottom boundary of plot (axis line)
% Remove bottom 3% of plot area to avoid detecting the x-axis line (reduced from 5%)
bottom_boundary_cutoff = round(y_max_pixel - (y_max_pixel - y_min_pixel) * 0.03);
black_mask_all(bottom_boundary_cutoff:end, :) = false;

% Also remove other boundaries (left, right, top) - reduced from 5% to 3%
top_boundary_cutoff = round(y_min_pixel + (y_max_pixel - y_min_pixel) * 0.03);
left_boundary_cutoff = round(x_min_pixel + (x_max_pixel - x_min_pixel) * 0.03);
right_boundary_cutoff = round(x_max_pixel - (x_max_pixel - x_min_pixel) * 0.03);
black_mask_all(1:top_boundary_cutoff, :) = false;
black_mask_all(:, 1:left_boundary_cutoff) = false;
black_mask_all(:, right_boundary_cutoff:end) = false;

% NEW: Only detect black lines where x <= 3.5 (left half of the plot)
x_cutoff_pixel = x_min_pixel + (x_max_pixel - x_min_pixel) * (3.5 / 7.0);
black_mask_all(:, round(x_cutoff_pixel):end) = false;

black_mask_all = bwareaopen(black_mask_all, 3);  % Reduced from 5 to 3

cc_black = bwconncomp(black_mask_all);
stats_black = regionprops(cc_black, 'Area', 'PixelIdxList', 'Centroid');

black_solid_mask = false(size(black_mask_all));
black_dots_mask = false(size(black_mask_all));

% Calculate y position corresponding to y_data = 30
% In pixel space: higher y_pixel = lower y_data (inverted)
% y_data = 30 means lower part of plot
y_cutoff_pixel = y_max_pixel - (y_max_pixel - y_min_pixel) * (30.0 / 80.0);

for i = 1:length(stats_black)
    area = stats_black(i).Area;
    centroid_y = stats_black(i).Centroid(2);  % y-coordinate of centroid
    
    % Check if component is above y=30 threshold (y_pixel < y_cutoff)
    above_threshold = centroid_y < y_cutoff_pixel;
    
    % Solid line: large connected areas AND above y=30
    if area > 400 && above_threshold
        black_solid_mask(stats_black(i).PixelIdxList) = true;
    % Dots: everything else (small areas or below threshold)
    else
        black_dots_mask(stats_black(i).PixelIdxList) = true;
    end
end

black_solid_mask = imclose(black_solid_mask, strel('disk', 4));
black_dots_mask = bwareaopen(black_dots_mask, 3);

%% STEP 4: Extract line data WITHOUT smoothing (raw data)

[x_red_solid_px, y_red_solid_px] = extract_line_raw(red_solid_mask);
[x_red_dashdot_px, y_red_dashdot_px] = extract_line_raw(red_dashdot_mask);
[x_red_dots_px, y_red_dots_px] = extract_line_raw(red_dots_mask);

[x_blue_solid_px, y_blue_solid_px] = extract_line_raw(blue_solid_mask);
[x_blue_dots_px, y_blue_dots_px] = extract_line_raw(blue_dots_mask);

[x_black_solid_px, y_black_solid_px] = extract_line_raw(black_solid_mask);
[x_black_dots_px, y_black_dots_px] = extract_line_raw(black_dots_mask);

% Debug: Print extraction info
fprintf('\n=== Extraction Debug Info ===\n');
fprintf('Red solid: %d pixels extracted\n', length(x_red_solid_px));
fprintf('Red dash-dot: %d pixels extracted\n', length(x_red_dashdot_px));
fprintf('Red dots: %d pixels extracted\n', length(x_red_dots_px));
fprintf('Blue solid: %d pixels extracted\n', length(x_blue_solid_px));
fprintf('Blue dots: %d pixels extracted\n', length(x_blue_dots_px));
fprintf('Black solid: %d pixels extracted\n', length(x_black_solid_px));
fprintf('Black dots: %d pixels extracted\n', length(x_black_dots_px));

%% STEP 5: Convert pixel coordinates to actual data coordinates

[x_red_solid_data, y_red_solid_data] = pixel_to_data(x_red_solid_px, y_red_solid_px, ...
    x_min_pixel, x_max_pixel, y_min_pixel, y_max_pixel, ...
    x_min_data, x_max_data, y_min_data, y_max_data);

[x_red_dashdot_data, y_red_dashdot_data] = pixel_to_data(x_red_dashdot_px, y_red_dashdot_px, ...
    x_min_pixel, x_max_pixel, y_min_pixel, y_max_pixel, ...
    x_min_data, x_max_data, y_min_data, y_max_data);

[x_red_dots_data, y_red_dots_data] = pixel_to_data(x_red_dots_px, y_red_dots_px, ...
    x_min_pixel, x_max_pixel, y_min_pixel, y_max_pixel, ...
    x_min_data, x_max_data, y_min_data, y_max_data);

[x_blue_solid_data, y_blue_solid_data] = pixel_to_data(x_blue_solid_px, y_blue_solid_px, ...
    x_min_pixel, x_max_pixel, y_min_pixel, y_max_pixel, ...
    x_min_data, x_max_data, y_min_data, y_max_data);

[x_blue_dots_data, y_blue_dots_data] = pixel_to_data(x_blue_dots_px, y_blue_dots_px, ...
    x_min_pixel, x_max_pixel, y_min_pixel, y_max_pixel, ...
    x_min_data, x_max_data, y_min_data, y_max_data);

[x_black_solid_data, y_black_solid_data] = pixel_to_data(x_black_solid_px, y_black_solid_px, ...
    x_min_pixel, x_max_pixel, y_min_pixel, y_max_pixel, ...
    x_min_data, x_max_data, y_min_data, y_max_data);

[x_black_dots_data, y_black_dots_data] = pixel_to_data(x_black_dots_px, y_black_dots_px, ...
    x_min_pixel, x_max_pixel, y_min_pixel, y_max_pixel, ...
    x_min_data, x_max_data, y_min_data, y_max_data);

%% STEP 6: Remove resampling - keep raw data only
% Skip resampling - data is already converted to actual coordinates

%% STEP 7: Save extracted data to MAT file

% Save all raw data to a single MAT file with the same name as the image
save(output_filename, ...
    'x_red_solid_data', 'y_red_solid_data', ...
    'x_red_dashdot_data', 'y_red_dashdot_data', ...
    'x_red_dots_data', 'y_red_dots_data', ...
    'x_blue_solid_data', 'y_blue_solid_data', ...
    'x_blue_dots_data', 'y_blue_dots_data', ...
    'x_black_solid_data', 'y_black_solid_data', ...
    'x_black_dots_data', 'y_black_dots_data');

fprintf('\n=== DATA SAVED ===\n');
fprintf('Output file: %s\n', output_filename);
fprintf('All raw extracted data saved to MAT file.\n');

%% STEP 8: Reproduce the figure with raw data and save as PNG
fig = figure('Position', [100, 100, 900, 700]);
hold on; grid on; box on;

% Plot extracted data with correct line styles
% Red: solid line, dash-dot line, dots only
plot(x_red_solid_data, y_red_solid_data, 'r-', 'LineWidth', 2.5);
plot(x_red_dashdot_data, y_red_dashdot_data, 'r-.', 'LineWidth', 2);
plot(x_red_dots_data, y_red_dots_data, 'r.', 'MarkerSize', 12);

% Blue: solid line, dots only
plot(x_blue_solid_data, y_blue_solid_data, 'b-', 'LineWidth', 2.5);
plot(x_blue_dots_data, y_blue_dots_data, 'b.', 'MarkerSize', 12);

% Black: solid line, dots only
plot(x_black_solid_data, y_black_solid_data, 'k-', 'LineWidth', 2.5);
plot(x_black_dots_data, y_black_dots_data, 'k.', 'MarkerSize', 12);

xlabel('Position (discrete index)', 'FontSize', 13);
ylabel('Jet/flame penetration (mm)', 'FontSize', 13);
xlim([0 7]);
ylim([0 80]);
set(gca, 'XTick', 0:1:7);
set(gca, 'YTick', 0:10:80);
legend('Red solid', 'Red dash-dot', 'Red dots', ...
    'Blue solid', 'Blue dots', ...
    'Black solid', 'Black dots', 'Location', 'best');
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
title('Extracted Raw Data', 'FontSize', 14);
hold off;

% Save figure as PNG
output_png = fullfile(img_path, [img_name '_extracted_data.png']);
saveas(fig, output_png);
fprintf('Figure saved to: %s\n', output_png);

%% STEP 9: Display detailed comparison (keep plotting uniform data for visualization)
figure('Position', [100, 100, 1800, 1000]);

% Original image
subplot(3,4,1);
imshow(img);
title('Original Image', 'FontSize', 13, 'FontWeight', 'bold');

% Red solid mask
subplot(3,4,2);
imshow(red_solid_mask);
title('Red Solid Line Mask', 'FontSize', 11);

% Red dash-dot mask
subplot(3,4,3);
imshow(red_dashdot_mask);
title('Red Dash-Dot Line Mask', 'FontSize', 11);

% Red dots mask
subplot(3,4,4);
imshow(red_dots_mask);
title('Red Dots Mask', 'FontSize', 11);

% Blue solid mask
subplot(3,4,5);
imshow(blue_solid_mask);
title('Blue Solid Line Mask', 'FontSize', 11);

% Blue dots mask
subplot(3,4,6);
imshow(blue_dots_mask);
title('Blue Dots Mask', 'FontSize', 11);

% Black solid mask
subplot(3,4,7);
imshow(black_solid_mask);
title('Black Solid Line Mask', 'FontSize', 11);

% Black dots mask
subplot(3,4,8);
imshow(black_dots_mask);
title('Black Dots Mask', 'FontSize', 11);

% Overlay all red masks on original
subplot(3,4,9);
imshow(img); hold on;
[r_rows, r_cols] = find(red_solid_mask);
plot(r_cols, r_rows, 'r.', 'MarkerSize', 1);
[r_rows, r_cols] = find(red_dashdot_mask);
plot(r_cols, r_rows, 'r.', 'MarkerSize', 1);
[r_rows, r_cols] = find(red_dots_mask);
plot(r_cols, r_rows, 'r.', 'MarkerSize', 1);
title('All Red Lines Overlay', 'FontSize', 11);
hold off;

% Final extracted plot - Red lines
subplot(3,4,10);
hold on; grid on; box on;
plot(x_red_solid_data, y_red_solid_data, 'r-', 'LineWidth', 2);
plot(x_red_dashdot_data, y_red_dashdot_data, 'r-.', 'LineWidth', 2);
plot(x_red_dots_data, y_red_dots_data, 'r.', 'MarkerSize', 12);
xlabel('Position', 'FontSize', 11);
ylabel('Penetration (mm)', 'FontSize', 11);
xlim([0 7]); ylim([0 80]);
set(gca, 'XTick', 0:1:7);
title('Red Lines Extracted', 'FontSize', 11);
legend('Solid', 'Dash-dot', 'Dots', 'Location', 'best', 'FontSize', 9);
hold off;

% Final extracted plot - Blue lines
subplot(3,4,11);
hold on; grid on; box on;
plot(x_blue_solid_data, y_blue_solid_data, 'b-', 'LineWidth', 2);
plot(x_blue_dots_data, y_blue_dots_data, 'b.', 'MarkerSize', 12);
xlabel('Position', 'FontSize', 11);
ylabel('Penetration (mm)', 'FontSize', 11);
xlim([0 7]); ylim([0 80]);
set(gca, 'XTick', 0:1:7);
title('Blue Lines Extracted', 'FontSize', 11);
legend('Solid', 'Dots', 'Location', 'best', 'FontSize', 9);
hold off;

% Final extracted plot - Black lines
subplot(3,4,12);
hold on; grid on; box on;
plot(x_black_solid_data, y_black_solid_data, 'k-', 'LineWidth', 2);
plot(x_black_dots_data, y_black_dots_data, 'k.', 'MarkerSize', 12);
xlabel('Position', 'FontSize', 11);
ylabel('Penetration (mm)', 'FontSize', 11);
xlim([0 7]); ylim([0 80]);
set(gca, 'XTick', 0:1:7);
title('Black Lines Extracted', 'FontSize', 11);
legend('Solid', 'Dots', 'Location', 'best', 'FontSize', 9);
hold off;

%% STEP 10: Remove quality check section - not needed for raw data

%% Display summary
fprintf('\n=== EXTRACTION COMPLETE ===\n');
fprintf('RED LINES:\n');
fprintf('  Solid: %d raw points\n', length(x_red_solid_data));
fprintf('  Dash-dot: %d raw points\n', length(x_red_dashdot_data));
fprintf('  Dots: %d raw points\n', length(x_red_dots_data));
fprintf('\nBLUE LINES:\n');
fprintf('  Solid: %d raw points\n', length(x_blue_solid_data));
fprintf('  Dots: %d raw points\n', length(x_blue_dots_data));
fprintf('\nBLACK LINES:\n');
fprintf('  Solid: %d raw points\n', length(x_black_solid_data));
fprintf('  Dots: %d raw points\n', length(x_black_dots_data));
fprintf('\nData saved to: %s\n', output_filename);
fprintf('==========================\n');


%% ========== FUNCTIONS (must be at end for older MATLAB) ==========

function [x_data, y_data] = extract_line_raw(mask)
    % Extract raw line data without any smoothing
    [rows, cols] = find(mask);
    
    if isempty(rows)
        x_data = [];
        y_data = [];
        return;
    end
    
    % Group by column (x-coordinate) and take median
    x_unique = unique(cols);
    y_median = zeros(size(x_unique));
    for i = 1:length(x_unique)
        y_median(i) = median(rows(cols == x_unique(i)));
    end
    
    % Sort by x
    [x_data, sort_idx] = sort(x_unique);
    y_data = y_median(sort_idx);
end

function [x_data, y_data] = extract_and_smooth_line(mask, smooth_window)
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
    
    if length(y_data) > smooth_window
        y_data = medfilt1(y_data, smooth_window);
    end
    
    if length(y_data) > 5
        window = ones(smooth_window, 1) / smooth_window;
        y_data = conv(y_data, window, 'same');
    end
end

function [x_actual, y_actual] = pixel_to_data(x_pixel, y_pixel, ...
    x_min_pix, x_max_pix, y_min_pix, y_max_pix, ...
    x_min_dat, x_max_dat, y_min_dat, y_max_dat)
    
    x_actual = x_min_dat + (x_pixel - x_min_pix) .* ...
        (x_max_dat - x_min_dat) / (x_max_pix - x_min_pix);
    
    y_actual = y_max_dat - (y_pixel - y_min_pix) .* ...
        (y_max_dat - y_min_dat) / (y_max_pix - y_min_pix);
end

function [x_uniform, y_uniform] = resample_curve(x_data, y_data, x_min, x_max, num_points)
    if isempty(x_data)
        x_uniform = [];
        y_uniform = [];
        return;
    end
    
    x_uniform = linspace(x_min, x_max, num_points);
    y_uniform = interp1(x_data, y_data, x_uniform, 'pchip', 'extrap');
    
    valid_range = (x_uniform >= min(x_data)) & (x_uniform <= max(x_data));
    x_uniform = x_uniform(valid_range);
    y_uniform = y_uniform(valid_range);
end