%% Complete Black Dotted Line Detection with Manual Calibration
% This code detects dotted lines and uses manual calibration for accurate coordinates
clear all, close all, clc

%% Configuration - Set your save location here
save_location = 'C:\Users\jun-y\OneDrive - UNSW\CH4DDI\Results\';  % Change this to your desired save path
output_filename = 'H2_jet_penetration.mat';

%% Step 1: Load the image
img = imread('C:\Users\jun-y\OneDrive - UNSW\CH4DDI\CH4DDI vs H2DDI CVCC\Patrick2022_handled\Patrick2022_Jetpenetration_detect .png');

[img_height, img_width, channels] = size(img);

if channels == 3
    gray_img = rgb2gray(img);
else
    gray_img = img;
end

%% Step 2: Detect all black/dark pixels
threshold = 220;
black_mask = gray_img < threshold;
black_mask = bwareaopen(black_mask, 2);

%% Step 3: Find all dot-like components
CC = bwconncomp(black_mask, 8);
stats = regionprops(CC, 'Area', 'BoundingBox', 'Centroid', 'PixelList', 'MajorAxisLength', 'MinorAxisLength');

fprintf('Total components found: %d\n', CC.NumObjects);

% Filter for dot-like components
all_dots = [];
dot_indices = [];
solid_components = [];

for i = 1:CC.NumObjects
    area = stats(i).Area;
    bbox = stats(i).BoundingBox;
    width = bbox(3);
    height = bbox(4);
    
    % Aspect ratio
    aspect_ratio = max(width, height) / max(min(width, height), 1);
    
    if area >= 3 && area <= 200 && aspect_ratio <= 3
        all_dots = [all_dots; stats(i).Centroid];
        dot_indices = [dot_indices; i];
    elseif area > 200 || aspect_ratio > 5
        solid_components = [solid_components; i];
    else
        perimeter = 2 * (width + height);
        compactness = (perimeter^2) / (4 * pi * area);
        
        if compactness < 2.5 && area <= 150
            all_dots = [all_dots; stats(i).Centroid];
            dot_indices = [dot_indices; i];
        else
            solid_components = [solid_components; i];
        end
    end
end

fprintf('Initial dot candidates: %d\n', size(all_dots, 1));

%% Step 4: Track the curved dotted line
if ~isempty(all_dots)
    [sorted_x, sort_idx] = sort(all_dots(:,1));
    sorted_dots = all_dots(sort_idx, :);
    
    fprintf('Total dots available for tracking: %d\n', size(sorted_dots, 1));
    
    curve_dots = [];
    used_dots = false(size(sorted_dots, 1), 1);
    
    start_x = min(sorted_dots(:,1));
    start_candidates = sorted_dots(sorted_dots(:,1) < start_x + 30, :);
    
    if ~isempty(start_candidates)
        best_curve = [];
        best_length = 0;
        
        for start_attempt = 1:min(5, size(start_candidates, 1))
            fprintf('Trying starting point %d\n', start_attempt);
            
            temp_used = false(size(sorted_dots, 1), 1);
            
            if start_attempt == 1
                [~, start_idx] = min(start_candidates(:,2));
            elseif start_attempt == 2
                [~, start_idx] = max(start_candidates(:,2));
            elseif start_attempt == 3
                start_idx = round(size(start_candidates, 1)/2);
            else
                start_idx = start_attempt - 2;
            end
            
            current_dot = start_candidates(start_idx, :);
            
            start_sorted_idx = 0;
            for i = 1:size(sorted_dots, 1)
                if norm(sorted_dots(i,:) - current_dot) < 1
                    temp_used(i) = true;
                    start_sorted_idx = i;
                    break;
                end
            end
            
            if start_sorted_idx == 0
                continue;
            end
            
            temp_curve = current_dot;
            max_dist = 60;
            iteration = 0;
            max_iterations = 200;
            
            while iteration < max_iterations
                iteration = iteration + 1;
                
                candidates = [];
                candidate_indices = [];
                
                for i = 1:size(sorted_dots, 1)
                    if ~temp_used(i)
                        dot = sorted_dots(i, :);
                        dist = norm(dot - current_dot);
                        
                        if dist < max_dist && dist > 0.5
                            candidates = [candidates; dot];
                            candidate_indices = [candidate_indices; i];
                        end
                    end
                end
                
                if isempty(candidates)
                    break;
                end
                
                if size(temp_curve, 1) >= 2
                    prev_direction = temp_curve(end,:) - temp_curve(end-1,:);
                    
                    best_score = inf;
                    best_idx = 1;
                    
                    for j = 1:size(candidates, 1)
                        new_direction = candidates(j,:) - current_dot;
                        dist = norm(candidates(j,:) - current_dot);
                        
                        dist_score = dist / max_dist;
                        
                        angle_score = 0;
                        if norm(prev_direction) > 1e-10 && norm(new_direction) > 1e-10
                            dot_product = prev_direction(1) * new_direction(1) + prev_direction(2) * new_direction(2);
                            magnitude_product = norm(prev_direction) * norm(new_direction);
                            cos_angle = dot_product / magnitude_product;
                            cos_angle = max(-1, min(1, cos_angle));
                            angle = acos(cos_angle);
                            if isfinite(angle)
                                angle_score = angle / pi;
                            end
                        end
                        
                        x_movement = candidates(j,1) - current_dot(1);
                        if x_movement >= 0
                            x_score = 0;
                        else
                            x_score = abs(x_movement) / max_dist;
                        end
                        
                        score = 0.4 * dist_score + 0.3 * angle_score + 0.3 * x_score;
                        
                        if score < best_score
                            best_score = score;
                            best_idx = j;
                        end
                    end
                else
                    dists = sqrt(sum((candidates - repmat(current_dot, size(candidates,1), 1)).^2, 2));
                    [~, best_idx] = min(dists);
                end
                
                current_dot = candidates(best_idx, :);
                temp_curve = [temp_curve; current_dot];
                temp_used(candidate_indices(best_idx)) = true;
            end
            
            if size(temp_curve, 1) > best_length
                best_curve = temp_curve;
                best_length = size(temp_curve, 1);
                used_dots = temp_used;
            end
            
            fprintf('Attempt %d: Found %d dots\n', start_attempt, size(temp_curve, 1));
        end
        
        curve_dots = best_curve;
    end
    
    % Extend curve
    if ~isempty(curve_dots)
        fprintf('Attempting to extend curve...\n');
        
        max_extend_dist = 80;
        extended = true;
        extend_iterations = 0;
        
        while extended && extend_iterations < 50
            extend_iterations = extend_iterations + 1;
            extended = false;
            
            current_dot = curve_dots(end, :);
            
            for i = 1:size(sorted_dots, 1)
                if ~used_dots(i)
                    dot = sorted_dots(i, :);
                    dist = norm(dot - current_dot);
                    
                    if dist < max_extend_dist && dist > 0.5
                        x_movement = dot(1) - current_dot(1);
                        
                        if x_movement > -20
                            curve_dots = [curve_dots; dot];
                            used_dots(i) = true;
                            extended = true;
                            fprintf('Extended curve with dot at (%.1f, %.1f)\n', dot(1), dot(2));
                            break;
                        end
                    end
                end
            end
        end
    end
    
    if ~isempty(curve_dots)
        final_x = curve_dots(:, 1);
        final_y = curve_dots(:, 2);
    else
        final_x = [];
        final_y = [];
    end
else
    final_x = [];
    final_y = [];
end

fprintf('Dots following curve: %d\n', length(final_x));

%% Step 5: Display detection results
figure('Position', [50, 50, 1200, 800]);

subplot(2,3,1);
imshow(img);
title('Original Image');

subplot(2,3,2);
imshow(black_mask);
title('All Black Pixels');

subplot(2,3,3);
imshow(img);
hold on;
if ~isempty(final_x)
    plot(final_x, final_y, 'r.', 'MarkerSize', 12);
    plot(final_x, final_y, 'r-', 'LineWidth', 1);
end
title('Detected Dotted Line');

%% Step 6: Manual Calibration
% Show image for manual calibration
subplot(2,3,4:6);
imshow(img);
title('Manual Calibration - Click 4 Points');

% Instructions
annotation('textbox', [0.02, 0.45, 0.6, 0.05], ...
    'String', 'Click 4 points in order: (0,0), (6,0), (6,95), (0,95)', ...
    'FontSize', 12, 'BackgroundColor', 'yellow');

% Collect 4 calibration points
calibration_points = zeros(4, 2);
point_labels = {'(0,0) - Bottom Left', '(6,0) - Bottom Right', ...
               '(6,95) - Top Right', '(0,95) - Top Left'};

for i = 1:4
    fprintf('Click point %d: %s\n', i, point_labels{i});
    title(sprintf('Click point %d: %s', i, point_labels{i}));
    
    [x, y] = ginput(1);
    calibration_points(i, :) = [x, y];
    
    hold on;
    plot(x, y, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    text(x+10, y, sprintf('P%d', i), 'Color', 'red', 'FontSize', 12, 'FontWeight', 'bold');
end

%% Step 7: Apply manual calibration
if ~isempty(final_x)
    fprintf('\nApplying manual calibration...\n');
    
    % Extract calibration points
    p1 = calibration_points(1, :); % (0,0)
    p2 = calibration_points(2, :); % (6,0)
    p3 = calibration_points(3, :); % (6,95)
    p4 = calibration_points(4, :); % (0,95)
    
    x_data = zeros(size(final_x));
    y_data = zeros(size(final_y));
    
    for i = 1:length(final_x)
        px = final_x(i);
        py = final_y(i);
        
        % Bilinear interpolation
        x1 = p1(1); y1 = p1(2); % (0,0)
        x2 = p2(1); y2 = p2(2); % (6,0)
        x4 = p4(1); y4 = p4(2); % (0,95)
        
        % Simple rectangular mapping
        left_x = x1;
        right_x = x2;
        width_pixels = right_x - left_x;
        
        bottom_y = y1;
        top_y = y4;
        height_pixels = top_y - bottom_y;
        
        % Convert to normalized coordinates [0,1]
        norm_x = (px - left_x) / width_pixels;
        norm_y = (py - bottom_y) / height_pixels;
        
        % Map to data coordinates
        x_data(i) = norm_x * 6.0;           % 0 to 6
        y_data(i) = norm_y * 95.0;          % 0 to 95 (no flip needed - the normalization already handles it)
    end
    
    % Add (0,0) as the starting point for jet penetration
    x_data = [0; x_data(:)];
    y_data = [0; y_data(:)];
    
    fprintf('Added (0,0) starting point to data\n');
    
    % Display results
    figure('Position', [100, 100, 800, 600]);
    plot(x_data, y_data, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
    xlim([0 6]);
    ylim([0 100]);
    set(gca, 'XTick', 0:1:6);
    set(gca, 'YTick', [0 10 20 30 40 50 60 70 80 90 95]);
    grid on;
    xlabel('X');
    ylabel('Y');
    title('Manually Calibrated Dotted Line Data');
    
    % Display all data points
    fprintf('\n========================================\n');
    fprintf('Final Calibrated Results:\n');
    fprintf('========================================\n');
    fprintf('Total points: %d\n', length(x_data));
    fprintf('X range: [%.3f, %.3f]\n', min(x_data), max(x_data));
    fprintf('Y range: [%.3f, %.3f]\n', min(y_data), max(y_data));
    
    fprintf('\nAll extracted data points (X, Y):\n');
    for i = 1:length(x_data)
        fprintf('  Point %d: (%.4f, %.4f)\n', i, x_data(i), y_data(i));
    end
    
    % Save results to specified location
    if ~exist(save_location, 'dir')
        mkdir(save_location);
        fprintf('Created directory: %s\n', save_location);
    end
    
    full_save_path = fullfile(save_location, output_filename);
    save(full_save_path, 'x_data', 'y_data');
    fprintf('\nData saved to: %s\n', full_save_path);
    fprintf('Variables saved: x_data, y_data\n');
else
    fprintf('No dots detected to calibrate\n');
end