% Compare CH4DDI and H2DDI data - Modified for new data structures
clear; clc; close all;

% Load both MAT files
file1 = 'C:\Users\jun-y\OneDrive - UNSW\CH4DDI\CH4DDI vs H2DDI CVCC\Mark2024_handled\Mark2024_Figure_extractedData_D_1_03ms_M.mat';
file2 = 'C:\Users\jun-y\OneDrive - UNSW\CH4DDI\CH4DDI vs H2DDI CVCC\Patrick2022_handled\D-0.93ms-H_Patrick2022_Figure10_extracted_data.mat';

% Load the data files
data_ch4 = load(file1);
data_h2 = load(file2);

fprintf('=== Data loaded successfully ===\n');

% Extract CH4 data from structure
if isfield(data_ch4, 'data_ch4_D_1_03ms_M')
    ch4_data = data_ch4.data_ch4_D_1_03ms_M;
elseif ~isempty(fieldnames(data_ch4))
    fields = fieldnames(data_ch4);
    ch4_data = data_ch4.(fields{1});
else
    error('Cannot find data structure in CH4 file');
end 

% Extract H2 data - check if it's already in the right format or needs extraction
if isfield(data_h2, 'x_black_dotted_data')
    % Data is already in the expected format
    h2_data = data_h2;
elseif isfield(data_h2, 'data_h2')
    h2_data = data_h2.data_h2;
elseif ~isempty(fieldnames(data_h2))
    % Try to extract the first field
    fields = fieldnames(data_h2);
    if isstruct(data_h2.(fields{1}))
        h2_data = data_h2.(fields{1});
    else
        h2_data = data_h2;
    end
else
    h2_data = data_h2;
end

% Extract CH4 data fields
x_ch4_nonreacting = [];
y_ch4_nonreacting = [];
x_ch4_penetration = [];
y_ch4_penetration = [];
x_ch4_recession = [];
y_ch4_recession = [];

if isfield(ch4_data, 'NonReactingJetPenetration') && isstruct(ch4_data.NonReactingJetPenetration)
    fields_nr = fieldnames(ch4_data.NonReactingJetPenetration);
    if length(fields_nr) >= 2
        x_ch4_nonreacting = ch4_data.NonReactingJetPenetration.(fields_nr{1});
        y_ch4_nonreacting = ch4_data.NonReactingJetPenetration.(fields_nr{2});
    end
end

if isfield(ch4_data, 'Penetration_average') && isstruct(ch4_data.Penetration_average)
    fields_pen = fieldnames(ch4_data.Penetration_average);
    if length(fields_pen) >= 2
        x_ch4_penetration = ch4_data.Penetration_average.(fields_pen{1});
        y_ch4_penetration = ch4_data.Penetration_average.(fields_pen{2});
    end
end

if isfield(ch4_data, 'Recession_average') && isstruct(ch4_data.Recession_average)
    fields_rec = fieldnames(ch4_data.Recession_average);
    if length(fields_rec) >= 2
        x_ch4_recession = ch4_data.Recession_average.(fields_rec{1});
        y_ch4_recession = ch4_data.Recession_average.(fields_rec{2});
    end
end

% Define colors
% H2: Different shades of blue
h2_light_blue = [0.3 0.6 1];    % Light blue for inert jet
h2_medium_blue = [0 0.4 0.8];   % Medium blue for penetration
h2_dark_blue = [0 0.2 0.6];     % Dark blue for other data
h2_cyan = [0 0.8 0.8];          % Cyan for circle data

% CH4: Different shades of orange
ch4_dark_orange = [0.8 0.3 0];    % Dark orange for inert jet
ch4_orange = [1 0.5 0];           % Medium orange for penetration
ch4_light_orange = [1 0.65 0.2];  % Light orange for recession

%% Figure 1: Original data comparison
figure('Position', [50, 50, 1000, 700]);
hold on; grid on; box on;

% Plot H2 Data with new format
% Black dotted data - jet penetration (inert jet) - SMOOTH BEFORE x=1.5
if isfield(h2_data, 'x_black_dotted_data') && isfield(h2_data, 'y_black_dotted_data')
    if ~isempty(h2_data.x_black_dotted_data)
        % Create smoothed version for x < 1.5
        x_h2_inert = h2_data.x_black_dotted_data;
        y_h2_inert = h2_data.y_black_dotted_data;
        
        % Find indices where x < 1.5
        idx_early = x_h2_inert < 1.5;
        
        % Apply smoothing to the early part
        y_h2_inert_smooth = y_h2_inert;
        if any(idx_early)
            % Apply median filter and smoothing to early portion
            y_early_smooth = medfilt1(y_h2_inert(idx_early), 5);
            y_early_smooth = smooth(y_early_smooth, 10, 'loess');
            y_h2_inert_smooth(idx_early) = y_early_smooth;
        end
        
        plot(x_h2_inert, y_h2_inert_smooth, '-', ...
            'Color', h2_light_blue, 'LineWidth', 2.5, 'DisplayName', 'H2: inert jet penetration');
    end
end

% Green data - penetration average
if isfield(h2_data, 'x_green_data') && isfield(h2_data, 'y_green_data')
    if ~isempty(h2_data.x_green_data)
        plot(h2_data.x_green_data, h2_data.y_green_data, '-', ...
            'Color', h2_medium_blue, 'LineWidth', 2.5, 'DisplayName', 'H2: flame penetration average');
    end
end

% Blue data 
if isfield(h2_data, 'x_blue_data') && isfield(h2_data, 'y_blue_data')
    if ~isempty(h2_data.x_blue_data) && ~isscalar(h2_data.x_blue_data)
        plot(h2_data.x_blue_data, h2_data.y_blue_data, '-', ...
            'Color', h2_dark_blue, 'LineWidth', 2.5, 'DisplayName', 'H2: flame penetration');
    end
end

% Circle data (recession) - REMOVED as per request
% Not plotting H2 recession points

% Plot CH4 Data
if ~isempty(x_ch4_nonreacting) && ~isempty(y_ch4_nonreacting)
    plot(x_ch4_nonreacting, y_ch4_nonreacting, '-', ...
        'Color', ch4_dark_orange, 'LineWidth', 2.5, 'DisplayName', 'CH4: inert jet');
end

if ~isempty(x_ch4_penetration) && ~isempty(y_ch4_penetration)
    plot(x_ch4_penetration, y_ch4_penetration, '-', ...
        'Color', ch4_orange, 'LineWidth', 2.5, 'DisplayName', 'CH4: flame penetration average');
end

if ~isempty(x_ch4_recession) && ~isempty(y_ch4_recession)
    plot(x_ch4_recession, y_ch4_recession, '.', ...
        'Color', ch4_light_orange, 'MarkerSize', 15, 'DisplayName', 'CH4: flame recession');
end

xlabel('Time aSOI [ms]', 'FontSize', 13);
ylabel('Jet/flame penetration [mm]', 'FontSize', 13);
title('CH4 vs H2 Comparison', 'FontSize', 14, 'FontWeight', 'bold');
xlim([0 3.5]);  % Changed from [0 10] to [0 3.5]
ylim([0 80]);
set(gca, 'XTick', 0:0.5:3.5);  % Adjusted tick marks for new xlim
set(gca, 'YTick', 0:20:80);
set(gcf, 'Color', 'w');
legend('Location', 'best', 'FontSize', 10, 'Box', 'on');
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
hold off;

%% Compute velocities
% Function to compute derivative and apply smoothing
compute_velocity = @(x, y) gradient(y) ./ gradient(x);
smooth_data = @(v, window) movmean(v, window, 'omitnan');

% Smoothing window sizes
smooth_window = 20;
smooth_window2 = 50;
smooth_window3 = 30;

% H2 velocities
v_h2_black_dotted = [];
x_h2_black_dotted_shifted = [];
v_h2_green = [];
x_h2_green_shifted = [];
v_h2_blue = [];
x_h2_blue_shifted = [];

if isfield(h2_data, 'x_black_dotted_data') && isfield(h2_data, 'y_black_dotted_data')
    if ~isempty(h2_data.x_black_dotted_data)
        v_h2_black_dotted = compute_velocity(h2_data.x_black_dotted_data, h2_data.y_black_dotted_data);
        v_h2_black_dotted = smooth_data(v_h2_black_dotted, smooth_window2);
        x_h2_black_dotted_shifted = h2_data.x_black_dotted_data - min(h2_data.x_black_dotted_data);
    end
end

if isfield(h2_data, 'x_green_data') && isfield(h2_data, 'y_green_data')
    if ~isempty(h2_data.x_green_data)
        v_h2_green = compute_velocity(h2_data.x_green_data, h2_data.y_green_data);
        v_h2_green = smooth_data(v_h2_green, smooth_window);
        x_h2_green_shifted = h2_data.x_green_data - h2_data.x_blue_data - 0.373;
    end
end

if isfield(h2_data, 'x_blue_data') && isfield(h2_data, 'y_blue_data')
    if ~isempty(h2_data.x_blue_data) && ~isscalar(h2_data.x_blue_data)
        v_h2_blue = compute_velocity(h2_data.x_blue_data, h2_data.y_blue_data);
        v_h2_blue = smooth_data(v_h2_blue, smooth_window);
        x_h2_blue_shifted = h2_data.x_blue_data - h2_data.x_blue_data;
    end
end

% CH4 velocities
v_ch4_penetration = [];
x_ch4_penetration_shifted = [];
v_ch4_nonreacting = [];
x_ch4_nonreacting_shifted = [];

if ~isempty(x_ch4_penetration) && ~isempty(y_ch4_penetration)
    v_ch4_penetration = compute_velocity(x_ch4_penetration, y_ch4_penetration);
    v_ch4_penetration = smooth_data(v_ch4_penetration, smooth_window3);
    x_ch4_penetration_shifted = x_ch4_penetration - min(x_ch4_penetration)- 1.8667; %% Hard coded !!!;
end

if ~isempty(x_ch4_nonreacting) && ~isempty(y_ch4_nonreacting)
    v_ch4_nonreacting = compute_velocity(x_ch4_nonreacting, y_ch4_nonreacting);
    v_ch4_nonreacting = smooth_data(v_ch4_nonreacting, smooth_window);
    x_ch4_nonreacting_shifted = x_ch4_nonreacting - min(x_ch4_nonreacting);
end

%% Figure 2: Velocity comparison - focusing on penetration rates
figure('Position', [100, 100, 1000, 700]);
hold on; grid on; box on;

% H2 velocities - focusing on penetration (green data)
if ~isempty(v_h2_green)
    plot(x_h2_green_shifted, v_h2_green, '-', ...
        'Color', h2_medium_blue, 'LineWidth', 2.5, 'DisplayName', 'H2: flame penetration average');
end

if ~isempty(v_h2_blue)
    plot(x_h2_blue_shifted, v_h2_blue, '-', ...
        'Color', h2_dark_blue, 'LineWidth', 2.5, 'DisplayName', 'H2: flame penetration');
end

% Optionally include inert jet velocity for comparison
% if ~isempty(v_h2_black_dotted)
%     plot(x_h2_black_dotted_shifted, v_h2_black_dotted, '--', ...
%         'Color', h2_light_blue, 'LineWidth', 2, 'DisplayName', 'H2: inert jet velocity');
% end

% CH4 velocities
if ~isempty(v_ch4_penetration)
    plot(x_ch4_penetration_shifted, v_ch4_penetration, '-', ...
        'Color', ch4_orange, 'LineWidth', 2.5, 'DisplayName', 'CH4: flame penetration average');
end

% Optionally include CH4 inert jet velocity
% if ~isempty(v_ch4_nonreacting)
%     plot(x_ch4_nonreacting_shifted, v_ch4_nonreacting, '--', ...
%         'Color', ch4_dark_orange, 'LineWidth', 2, 'DisplayName', 'CH4: inert jet velocity');
% end

xlabel('Time after start of main fuel ignition events [ms]', 'FontSize', 13);
ylabel('Penetration rate [mm/ms]', 'FontSize', 13);
title('CH4 vs H2 Penetration Rate Comparison', 'FontSize', 14, 'FontWeight', 'bold');
xlim([0 2]);
ylim([0 150]);
set(gca, 'XTick', 0:0.25:2);
legend('Location', 'best', 'FontSize', 10, 'Box', 'on');
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
set(gcf, 'Color', 'w');
grid on;
hold off;

%% Display data structure information
fprintf('\n=== Data Structure Summary ===\n');

% H2 data fields
fprintf('\nH2 Data Fields Found:\n');
h2_fields = fieldnames(h2_data);
for i = 1:length(h2_fields)
    if ~isempty(h2_data.(h2_fields{i}))
        fprintf('  - %s: %s\n', h2_fields{i}, mat2str(size(h2_data.(h2_fields{i}))));
    end
end

% CH4 data fields
fprintf('\nCH4 Data Fields Found:\n');
ch4_fields = fieldnames(ch4_data);
for i = 1:length(ch4_fields)
    fprintf('  - %s\n', ch4_fields{i});
end

fprintf('\n=== Plotting complete ===\n');
fprintf('Figure 1: Original data (penetration vs time)\n');
fprintf('Figure 2: Smoothed velocity (dy/dx) - all lines shifted to start at x=0\n');
fprintf('Smoothing windows: %d, %d, %d points\n', smooth_window, smooth_window2, smooth_window3);