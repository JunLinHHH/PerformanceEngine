% Data extraction for new plot figure
% X-axis: 0 to 10 (0,2,4,6,8,10)
% Y-axis: 0 to 140 (0,20,40,60,80,100,120,140)
% Line styles: Inert jet (dotted), Flame recession (dash-dot), Flame penetration (solid)
% Colors: Black, Blue, Red

clear; clc; close all;

%% STEP 0: CALIBRATION HELPER
% Uncomment to manually calibrate
% figure; imshow(img); 
% title('Click 4 corners: 1)Bottom-left 2)Bottom-right 3)Top-left 4)Top-right');
% [x_cal, y_cal] = ginput(4);
% x_min_pixel = round(min(x_cal));
% x_max_pixel = round(max(x_cal));
% y_max_pixel = round(max(y_cal));
% y_min_pixel = round(min(y_cal));
% close;

%% Load image - UPDATE THIS PATH
img_path_full = 'C:\Users\jun-y\OneDrive - UNSW\CH4DDI\CH4 vs H2\Penetration\Source\CH4_200bar_Bob_r1060_black1200_blue1140.png';
img = imread(img_path_full);
[height, width, ~] = size(img);

% Extract filename for saving
[img_path, img_name, ~] = fileparts(img_path_full);
output_filename = fullfile(img_path, [img_name '_extracted_data.mat']);

%% STEP 1: Detect plot boundary
fprintf('=== STEP 1: Detecting plot boundary ===\n');

% Automatic detection based on non-frame pixels
gray = rgb2gray(img);
frame_mask = gray < 50;
[rows, cols] = find(~frame_mask);

if ~isempty(rows) && ~isempty(cols)
    x_min_pixel = min(cols);
    x_max_pixel = max(cols);
    y_min_pixel = min(rows);
    y_max_pixel = max(rows);
    
    plot_width = x_max_pixel - x_min_pixel;
    plot_height = y_max_pixel - y_min_pixel;
    
    if plot_width > 200 && plot_height > 200
        fprintf('Automatic boundary detection successful!\n');
    else
        fprintf('Using manual defaults.\n');
        x_min_pixel = 100;
        x_max_pixel = 900;
        y_min_pixel = 100;
        y_max_pixel = 900;
    end
else
    x_min_pixel = 100;
    x_max_pixel = 900;
    y_min_pixel = 100;
    y_max_pixel = 900;
end

fprintf('Plot boundary:\n');
fprintf('  x: [%d, %d], y: [%d, %d]\n', x_min_pixel, x_max_pixel, y_min_pixel, y_max_pixel);

% Create plot boundary mask
plot_mask = false(height, width);
plot_mask(round(y_min_pixel):round(y_max_pixel), round(x_min_pixel):round(x_max_pixel)) = true;

%% STEP 2: Axis calibration
x_min_data = 0;
x_max_data = 11;
y_min_data = 0;
y_max_data = 140;

fprintf('\n=== Axis Calibration ===\n');
fprintf('X-axis: %d to %d\n', x_min_data, x_max_data);
fprintf('Y-axis: %d to %d\n', y_min_data, y_max_data);

%% STEP 3: Extract lines by color and style

% BLACK LINES - Apply spatial constraints: x<=5, dash-dot x>=2.5, dotted x<=2.5, dash-dot y<=30, solid y>=30
black_mask_all = (img(:,:,1) < 80) & (img(:,:,2) < 80) & (img(:,:,3) < 80);
black_mask_all = black_mask_all & plot_mask;

% Remove legend and boundaries
legend_x_start = round(x_min_pixel + (x_max_pixel - x_min_pixel) * 0.5);
legend_y_start = round(y_min_pixel + (y_max_pixel - y_min_pixel) * 0.05);
legend_y_end = round(y_min_pixel + (y_max_pixel - y_min_pixel) * 0.25);
black_mask_all(legend_y_start:legend_y_end, legend_x_start:end) = false;

% Remove boundaries
margin = 0.03;
black_mask_all(1:round(y_min_pixel + (y_max_pixel-y_min_pixel)*margin), :) = false;
black_mask_all(round(y_max_pixel - (y_max_pixel-y_min_pixel)*margin):end, :) = false;
black_mask_all(:, 1:round(x_min_pixel + (x_max_pixel-x_min_pixel)*margin)) = false;
black_mask_all(:, round(x_max_pixel - (x_max_pixel-x_min_pixel)*margin):end) = false;

% Black: only detect where x <= 5
x_cutoff_black = x_min_pixel + (x_max_pixel - x_min_pixel) * (5.0 / 10.0);
black_mask_all(:, round(x_cutoff_black):end) = false;

black_mask_all = bwareaopen(black_mask_all, 3);

% Classify black lines by size AND position
cc_black = bwconncomp(black_mask_all);
stats_black = regionprops(cc_black, 'Area', 'PixelIdxList', 'Centroid');

black_solid_mask = false(size(black_mask_all));
black_dashdot_mask = false(size(black_mask_all));
black_dotted_mask = false(size(black_mask_all));

% Calculate pixel positions for thresholds
y_30_pixel = y_max_pixel - (y_max_pixel - y_min_pixel) * (30.0 / 140.0);
x_2_5_pixel = x_min_pixel + (x_max_pixel - x_min_pixel) * (2.8 / 10.0);

for i = 1:length(stats_black)
    area = stats_black(i).Area;
    centroid_x = stats_black(i).Centroid(1);
    centroid_y = stats_black(i).Centroid(2);
    
    % Black solid: area > 500 AND y >= 30
    if area > 500 && centroid_y < y_30_pixel
        black_solid_mask(stats_black(i).PixelIdxList) = true;
    % Black dash-dot: area 30-500 AND y <= 30 AND x >= 2.5
    elseif area > 30 && area <= 500 && centroid_y >= y_30_pixel && centroid_x >= x_2_5_pixel
        black_dashdot_mask(stats_black(i).PixelIdxList) = true;
    % Black dotted: small area AND x <= 2.5
    elseif centroid_x <= x_2_5_pixel
        black_dotted_mask(stats_black(i).PixelIdxList) = true;
    end
end

black_solid_mask = imclose(black_solid_mask, strel('disk', 4));
black_dashdot_mask = imclose(black_dashdot_mask, strel('disk', 3));

% BLUE LINES - dash-dot y<=45, solid y>=45
blue_mask_all = (img(:,:,3) > 150) & (img(:,:,1) < 100) & (img(:,:,2) < 150);
blue_mask_all = blue_mask_all & plot_mask;
blue_mask_all = bwareaopen(blue_mask_all, 3);

cc_blue = bwconncomp(blue_mask_all);
stats_blue = regionprops(cc_blue, 'Area', 'PixelIdxList', 'Centroid');

blue_solid_mask = false(size(blue_mask_all));
blue_dashdot_mask = false(size(blue_mask_all));
blue_dotted_mask = false(size(blue_mask_all));

% Calculate pixel position for y = 45
y_45_pixel = y_max_pixel - (y_max_pixel - y_min_pixel) * (45.0 / 140.0);

for i = 1:length(stats_blue)
    area = stats_blue(i).Area;
    centroid_y = stats_blue(i).Centroid(2);
    
    % Blue solid: area > 500 AND y >= 45
    if area > 500 && centroid_y < y_45_pixel
        blue_solid_mask(stats_blue(i).PixelIdxList) = true;
    % Blue dash-dot: area 30-500 AND y <= 45
    elseif area > 30 && centroid_y >= y_45_pixel
        blue_dashdot_mask(stats_blue(i).PixelIdxList) = true;
    % Blue dotted: small areas
    else
        blue_dotted_mask(stats_blue(i).PixelIdxList) = true;
    end
end

blue_solid_mask = imclose(blue_solid_mask, strel('disk', 4));
blue_dashdot_mask = imclose(blue_dashdot_mask, strel('disk', 3));

% RED LINES - dash-dot y<=63, solid y>=53
red_mask_all = (img(:,:,1) > 150) & (img(:,:,2) < 100) & (img(:,:,3) < 100);
red_mask_all = red_mask_all & plot_mask;
red_mask_all = bwareaopen(red_mask_all, 3);

cc_red = bwconncomp(red_mask_all);
stats_red = regionprops(cc_red, 'Area', 'PixelIdxList', 'Centroid');

red_solid_mask = false(size(red_mask_all));
red_dashdot_mask = false(size(red_mask_all));
red_dotted_mask = false(size(red_mask_all));

% Calculate pixel positions for y = 53 and y = 63
y_53_pixel = y_max_pixel - (y_max_pixel - y_min_pixel) * (53.0 / 140.0);
y_63_pixel = y_max_pixel - (y_max_pixel - y_min_pixel) * (63.0 / 140.0);

for i = 1:length(stats_red)
    area = stats_red(i).Area;
    centroid_y = stats_red(i).Centroid(2);
    
    % Red solid: area > 500 AND y >= 53
    if area > 500 && centroid_y < y_53_pixel
        red_solid_mask(stats_red(i).PixelIdxList) = true;
    % Red dash-dot: area 30-500 AND y <= 63
    elseif area > 30 && centroid_y >= y_63_pixel
        red_dashdot_mask(stats_red(i).PixelIdxList) = true;
    % Red dotted: small areas
    else
        red_dotted_mask(stats_red(i).PixelIdxList) = true;
    end
end

red_solid_mask = imclose(red_solid_mask, strel('disk', 4));
red_dashdot_mask = imclose(red_dashdot_mask, strel('disk', 3));

%% STEP 4: Extract raw data
[x_black_solid_px, y_black_solid_px] = extract_line_raw(black_solid_mask);
[x_black_dashdot_px, y_black_dashdot_px] = extract_line_raw(black_dashdot_mask);
[x_black_dotted_px, y_black_dotted_px] = extract_line_raw(black_dotted_mask);

[x_blue_solid_px, y_blue_solid_px] = extract_line_raw(blue_solid_mask);
[x_blue_dashdot_px, y_blue_dashdot_px] = extract_line_raw(blue_dashdot_mask);
[x_blue_dotted_px, y_blue_dotted_px] = extract_line_raw(blue_dotted_mask);

[x_red_solid_px, y_red_solid_px] = extract_line_raw(red_solid_mask);
[x_red_dashdot_px, y_red_dashdot_px] = extract_line_raw(red_dashdot_mask);
[x_red_dotted_px, y_red_dotted_px] = extract_line_raw(red_dotted_mask);

%% STEP 5: Convert to data coordinates
[x_black_solid_data, y_black_solid_data] = pixel_to_data(x_black_solid_px, y_black_solid_px, ...
    x_min_pixel, x_max_pixel, y_min_pixel, y_max_pixel, x_min_data, x_max_data, y_min_data, y_max_data);
[x_black_dashdot_data, y_black_dashdot_data] = pixel_to_data(x_black_dashdot_px, y_black_dashdot_px, ...
    x_min_pixel, x_max_pixel, y_min_pixel, y_max_pixel, x_min_data, x_max_data, y_min_data, y_max_data);
[x_black_dotted_data, y_black_dotted_data] = pixel_to_data(x_black_dotted_px, y_black_dotted_px, ...
    x_min_pixel, x_max_pixel, y_min_pixel, y_max_pixel, x_min_data, x_max_data, y_min_data, y_max_data);

[x_blue_solid_data, y_blue_solid_data] = pixel_to_data(x_blue_solid_px, y_blue_solid_px, ...
    x_min_pixel, x_max_pixel, y_min_pixel, y_max_pixel, x_min_data, x_max_data, y_min_data, y_max_data);
[x_blue_dashdot_data, y_blue_dashdot_data] = pixel_to_data(x_blue_dashdot_px, y_blue_dashdot_px, ...
    x_min_pixel, x_max_pixel, y_min_pixel, y_max_pixel, x_min_data, x_max_data, y_min_data, y_max_data);
[x_blue_dotted_data, y_blue_dotted_data] = pixel_to_data(x_blue_dotted_px, y_blue_dotted_px, ...
    x_min_pixel, x_max_pixel, y_min_pixel, y_max_pixel, x_min_data, x_max_data, y_min_data, y_max_data);

[x_red_solid_data, y_red_solid_data] = pixel_to_data(x_red_solid_px, y_red_solid_px, ...
    x_min_pixel, x_max_pixel, y_min_pixel, y_max_pixel, x_min_data, x_max_data, y_min_data, y_max_data);
[x_red_dashdot_data, y_red_dashdot_data] = pixel_to_data(x_red_dashdot_px, y_red_dashdot_px, ...
    x_min_pixel, x_max_pixel, y_min_pixel, y_max_pixel, x_min_data, x_max_data, y_min_data, y_max_data);
[x_red_dotted_data, y_red_dotted_data] = pixel_to_data(x_red_dotted_px, y_red_dotted_px, ...
    x_min_pixel, x_max_pixel, y_min_pixel, y_max_pixel, x_min_data, x_max_data, y_min_data, y_max_data);

%% STEP 6: Save to MAT file
save(output_filename, ...
    'x_black_solid_data', 'y_black_solid_data', ...
    'x_black_dashdot_data', 'y_black_dashdot_data', ...
    'x_black_dotted_data', 'y_black_dotted_data', ...
    'x_blue_solid_data', 'y_blue_solid_data', ...
    'x_blue_dashdot_data', 'y_blue_dashdot_data', ...
    'x_blue_dotted_data', 'y_blue_dotted_data', ...
    'x_red_solid_data', 'y_red_solid_data', ...
    'x_red_dashdot_data', 'y_red_dashdot_data', ...
    'x_red_dotted_data', 'y_red_dotted_data');

fprintf('\n=== DATA SAVED ===\n');
fprintf('Output file: %s\n', output_filename);

%% STEP 7: Plot and save figure
fig = figure('Position', [100, 100, 900, 700]);
hold on; grid on; box on;

% Black lines: Inert jet (dotted), Flame recession (dash-dot), Flame penetration (solid)
plot(x_black_solid_data, y_black_solid_data, 'k-', 'LineWidth', 2.5);
plot(x_black_dashdot_data, y_black_dashdot_data, 'k-.', 'LineWidth', 2);
plot(x_black_dotted_data, y_black_dotted_data, 'k:', 'LineWidth', 2);

% Blue lines
plot(x_blue_solid_data, y_blue_solid_data, 'b-', 'LineWidth', 2.5);
plot(x_blue_dashdot_data, y_blue_dashdot_data, 'b-.', 'LineWidth', 2);
plot(x_blue_dotted_data, y_blue_dotted_data, 'b:', 'LineWidth', 2);

% Red lines
plot(x_red_solid_data, y_red_solid_data, 'r-', 'LineWidth', 2.5);
plot(x_red_dashdot_data, y_red_dashdot_data, 'r-.', 'LineWidth', 2);
plot(x_red_dotted_data, y_red_dotted_data, 'r:', 'LineWidth', 2);

xlabel('X', 'FontSize', 13);
ylabel('Y', 'FontSize', 13);
xlim([0 10]);
ylim([0 140]);
set(gca, 'XTick', 0:2:10);
set(gca, 'YTick', 0:20:140);
legend('Black solid', 'Black dash-dot', 'Black dotted', ...
    'Blue solid', 'Blue dash-dot', 'Blue dotted', ...
    'Red solid', 'Red dash-dot', 'Red dotted', 'Location', 'best');
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
title('Extracted Raw Data', 'FontSize', 14);
hold off;

% Save PNG
output_png = fullfile(img_path, [img_name '_extracted_data.png']);
saveas(fig, output_png);
fprintf('Figure saved to: %s\n', output_png);

%% STEP 8: Display mask visualization
figure('Position', [100, 100, 1800, 1200]);

% Original image
subplot(4,3,1);
imshow(img);
title('Original Image', 'FontSize', 11, 'FontWeight', 'bold');

% Black masks
subplot(4,3,2);
imshow(black_solid_mask);
title('Black Solid Mask', 'FontSize', 10);

subplot(4,3,3);
imshow(black_dashdot_mask);
title('Black Dash-Dot Mask', 'FontSize', 10);

subplot(4,3,4);
imshow(black_dotted_mask);
title('Black Dotted Mask', 'FontSize', 10);

% Blue masks
subplot(4,3,5);
imshow(blue_solid_mask);
title('Blue Solid Mask', 'FontSize', 10);

subplot(4,3,6);
imshow(blue_dashdot_mask);
title('Blue Dash-Dot Mask', 'FontSize', 10);

subplot(4,3,7);
imshow(blue_dotted_mask);
title('Blue Dotted Mask', 'FontSize', 10);

% Red masks
subplot(4,3,8);
imshow(red_solid_mask);
title('Red Solid Mask', 'FontSize', 10);

subplot(4,3,9);
imshow(red_dashdot_mask);
title('Red Dash-Dot Mask', 'FontSize', 10);

subplot(4,3,10);
imshow(red_dotted_mask);
title('Red Dotted Mask', 'FontSize', 10);

% Extracted plots by color
subplot(4,3,11);
hold on; grid on; box on;
plot(x_black_solid_data, y_black_solid_data, 'k-', 'LineWidth', 2);
plot(x_black_dashdot_data, y_black_dashdot_data, 'k-.', 'LineWidth', 2);
plot(x_black_dotted_data, y_black_dotted_data, 'k:', 'LineWidth', 2);
xlabel('X', 'FontSize', 10);
ylabel('Y', 'FontSize', 10);
xlim([0 10]); ylim([0 140]);
title('Black Lines Extracted', 'FontSize', 10);
legend('Solid', 'Dash-dot', 'Dotted', 'Location', 'best', 'FontSize', 8);
hold off;

subplot(4,3,12);
hold on; grid on; box on;
plot(x_blue_solid_data, y_blue_solid_data, 'b-', 'LineWidth', 2);
plot(x_blue_dashdot_data, y_blue_dashdot_data, 'b-.', 'LineWidth', 2);
plot(x_blue_dotted_data, y_blue_dotted_data, 'b:', 'LineWidth', 2);
xlabel('X', 'FontSize', 10);
ylabel('Y', 'FontSize', 10);
xlim([0 10]); ylim([0 140]);
title('Blue Lines Extracted', 'FontSize', 10);
legend('Solid', 'Dash-dot', 'Dotted', 'Location', 'best', 'FontSize', 8);
hold off;

% Red lines extracted (in 4th row, spanning position)
axes('Position', [0.4, 0.05, 0.25, 0.18]);
hold on; grid on; box on;
plot(x_red_solid_data, y_red_solid_data, 'r-', 'LineWidth', 2);
plot(x_red_dashdot_data, y_red_dashdot_data, 'r-.', 'LineWidth', 2);
plot(x_red_dotted_data, y_red_dotted_data, 'r:', 'LineWidth', 2);
xlabel('X', 'FontSize', 10);
ylabel('Y', 'FontSize', 10);
xlim([0 10]); ylim([0 140]);
title('Red Lines Extracted', 'FontSize', 10);
legend('Solid', 'Dash-dot', 'Dotted', 'Location', 'best', 'FontSize', 8);
hold off;

fprintf('\n=== COMPLETE ===\n');

%% ========== FUNCTIONS ==========

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