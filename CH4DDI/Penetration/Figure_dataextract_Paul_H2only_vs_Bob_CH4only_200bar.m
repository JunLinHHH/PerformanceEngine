% Compare CH4 and H2 data - Selected lines only
clear; clc; close all;

% Load both MAT files
file1 = 'C:\Users\jun-y\OneDrive - UNSW\CH4DDI\CH4 vs H2\Penetration\Source\CH4_200bar_Bob_r1060_black1200_blue1140_extracted_data.mat';
file2 = 'C:\Users\jun-y\OneDrive - UNSW\CH4DDI\CH4 vs H2\Penetration\Source\H2_1060k_200bar_Paul_extracted_data.mat';

data_ch4 = load(file1);
data_h2 = load(file2);

fprintf('=== Data loaded successfully ===\n');

% Define colors
% H2: Different shades of blue
h2_light_blue = [0.3 0.6 1];    % Light blue for red dash-dot
h2_medium_blue = [0 0.4 0.8];   % Medium blue for blue lines
h2_dark_blue = [0 0.2 0.6];     % Dark blue for black lines

% CH4: Different shades of orange
ch4_dark_orange = [0.8 0.3 0];    % Dark orange for black dotted
ch4_orange = [1 0.5 0];           % Medium orange for red solid
ch4_light_orange = [1 0.65 0.2];  % Light orange for red dotted

%% Figure 1: Original data comparison
figure('Position', [50, 50, 1000, 700]);
hold on; grid on; box on;

% H2 Data (all in blue shades)
if isfield(data_h2, 'x_red_dashdot_data') && ~isempty(data_h2.x_red_dashdot_data)
    plot(data_h2.x_red_dashdot_data, data_h2.y_red_dashdot_data, '-', ...
        'Color', h2_light_blue, 'LineWidth', 2.5, 'DisplayName', 'H2: inert jet');
end

if isfield(data_h2, 'x_blue_solid_data') && ~isempty(data_h2.x_blue_solid_data)
    plot(data_h2.x_blue_solid_data, data_h2.y_blue_solid_data, '-', ...
        'Color', h2_medium_blue, 'LineWidth', 2.5, 'DisplayName', 'H2: flame penetration (multi-kernel)');
end
if isfield(data_h2, 'x_blue_dots_data') && ~isempty(data_h2.x_blue_dots_data)
    plot(data_h2.x_blue_dots_data, data_h2.y_blue_dots_data, '.', ...
        'Color', h2_medium_blue, 'MarkerSize', 15, 'DisplayName', 'H2: flame recession (multi-kernel)');
end

if isfield(data_h2, 'x_black_solid_data') && ~isempty(data_h2.x_black_solid_data)
    plot(data_h2.x_black_solid_data, data_h2.y_black_solid_data, '-', ...
        'Color', h2_dark_blue, 'LineWidth', 2.5, 'DisplayName', 'H2: flame penetration (single-kernel)');
end
if isfield(data_h2, 'x_black_dots_data') && ~isempty(data_h2.x_black_dots_data)
    plot(data_h2.x_black_dots_data, data_h2.y_black_dots_data, '.', ...
        'Color', h2_dark_blue, 'MarkerSize', 15, 'DisplayName', 'H2: flame recession (single-kernel)');
end

% CH4 Data (different shades of orange)
% if isfield(data_ch4, 'x_black_dotted_data') && ~isempty(data_ch4.x_black_dotted_data)
%     plot(data_ch4.x_black_dotted_data, data_ch4.y_black_dotted_data, '-', ...
%         'Color', ch4_dark_orange, 'LineWidth', 2.5, 'DisplayName', 'CH4: inert jet');
% end

if isfield(data_ch4, 'x_red_solid_data') && ~isempty(data_ch4.x_red_solid_data)
    plot(data_ch4.x_red_solid_data, data_ch4.y_red_solid_data, '-', ...
        'Color', ch4_orange, 'LineWidth', 2.5, 'DisplayName', 'CH4: flame penetration (1060 K)');
end

if isfield(data_ch4, 'x_red_dashdot_data') && ~isempty(data_ch4.x_red_dashdot_data)
    plot(data_ch4.x_red_dashdot_data, data_ch4.y_red_dashdot_data, '-.', ...
        'Color', ch4_light_orange, 'LineWidth', 2.5, 'DisplayName', 'CH4: flame recession (1060 K)');
end

xlabel('Time aSOI [ms]', 'FontSize', 13);
ylabel('Jet/flame penetration [mm]', 'FontSize', 13);
title('CH4 vs H2 Comparison', 'FontSize', 14, 'FontWeight', 'bold');
xlim([0 8.5]);
ylim([0 80]);
set(gca, 'XTick', 0:1:10);
set(gca, 'YTick', 0:20:80);
set(gcf, 'Color', 'w');   % only outside area → white
legend('Location', 'best', 'FontSize', 10, 'Box', 'on');
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
hold off;

%% Compute dy/dx (velocity) for each line and shift to start at x=0

% Option 1 Function to compute derivative and apply smoothing
compute_velocity = @(x, y) deal(gradient(y) ./ gradient(x));
smooth_data = @(v, window) movmean(v, window, 'omitnan');

% Smoothing window size
% smooth_window = 20; % Adjust this value for more/less smoothing
% smooth_window2 = 50;
% smooth_window3 = 30;
smooth_window = 5; % Adjust this value for more/less smoothing
smooth_window2 = 5;
smooth_window3 = 2;


% H2 velocities
if isfield(data_h2, 'x_red_dashdot_data') && ~isempty(data_h2.x_red_dashdot_data)
    v_h2_red = compute_velocity(data_h2.x_red_dashdot_data, data_h2.y_red_dashdot_data);
    v_h2_red = smooth_data(v_h2_red, smooth_window2);
    x_h2_red_shifted = data_h2.x_red_dashdot_data - min(data_h2.x_red_dashdot_data);
else
    v_h2_red = []; x_h2_red_shifted = [];
end

if isfield(data_h2, 'x_blue_solid_data') && ~isempty(data_h2.x_blue_solid_data)
    v_h2_blue_solid = compute_velocity(data_h2.x_blue_solid_data, data_h2.y_blue_solid_data);
    v_h2_blue_solid = smooth_data(v_h2_blue_solid, smooth_window);
    x_h2_blue_solid_shifted = data_h2.x_blue_solid_data - min(data_h2.x_blue_solid_data);
else
    v_h2_blue_solid = []; x_h2_blue_solid_shifted = [];
end

if isfield(data_h2, 'x_blue_dots_data') && ~isempty(data_h2.x_blue_dots_data)
    v_h2_blue_dots = compute_velocity(data_h2.x_blue_dots_data, data_h2.y_blue_dots_data);
    v_h2_blue_dots = smooth_data(v_h2_blue_dots, smooth_window);
    x_h2_blue_dots_shifted = data_h2.x_blue_dots_data - min(data_h2.x_blue_dots_data);
else
    v_h2_blue_dots = []; x_h2_blue_dots_shifted = [];
end

if isfield(data_h2, 'x_black_solid_data') && ~isempty(data_h2.x_black_solid_data)
    v_h2_black_solid = compute_velocity(data_h2.x_black_solid_data, data_h2.y_black_solid_data);
    v_h2_black_solid = smooth_data(v_h2_black_solid, smooth_window);
    x_h2_black_solid_shifted = data_h2.x_black_solid_data - min(data_h2.x_black_solid_data);
else
    v_h2_black_solid = []; x_h2_black_solid_shifted = [];
end

if isfield(data_h2, 'x_black_dots_data') && ~isempty(data_h2.x_black_dots_data)
    v_h2_black_dots = compute_velocity(data_h2.x_black_dots_data, data_h2.y_black_dots_data);
    v_h2_black_dots = smooth_data(v_h2_black_dots, smooth_window);
    x_h2_black_dots_shifted = data_h2.x_black_dots_data - min(data_h2.x_black_dots_data);
else
    v_h2_black_dots = []; x_h2_black_dots_shifted = [];
end

% CH4 velocities
if isfield(data_ch4, 'x_black_dotted_data') && ~isempty(data_ch4.x_black_dotted_data)
    v_ch4_black = compute_velocity(data_ch4.x_black_dotted_data, data_ch4.y_black_dotted_data);
    v_ch4_black = smooth_data(v_ch4_black, smooth_window);
    x_ch4_black_shifted = data_ch4.x_black_dotted_data - min(data_ch4.x_black_dotted_data);
else
    v_ch4_black = []; x_ch4_black_shifted = [];
end

if isfield(data_ch4, 'x_red_solid_data') && ~isempty(data_ch4.x_red_solid_data)
    v_ch4_red_solid = compute_velocity(data_ch4.x_red_solid_data, data_ch4.y_red_solid_data);
    v_ch4_red_solid = smooth_data(v_ch4_red_solid, smooth_window3);
    x_ch4_red_solid_shifted = data_ch4.x_red_solid_data - min(data_ch4.x_red_solid_data);
else
    v_ch4_red_solid = []; x_ch4_red_solid_shifted = [];
end

if isfield(data_ch4, 'x_red_dotted_data') && ~isempty(data_ch4.x_red_dotted_data)
    v_ch4_red_dots = compute_velocity(data_ch4.x_red_dotted_data, data_ch4.y_red_dotted_data);
    v_ch4_red_dots = smooth_data(v_ch4_red_dots, smooth_window3);
    x_ch4_red_dots_shifted = data_ch4.x_red_dotted_data - min(data_ch4.x_red_dotted_data);
else
    v_ch4_red_dots = []; x_ch4_red_dots_shifted = [];
end

%% Figure 2: Velocity (dy/dx) comparison - all starting from x=0
figure('Position', [100, 100, 1000, 700]);
hold on; grid on; box on;

% H2 velocities
% if ~isempty(v_h2_red)
%     plot(x_h2_red_shifted, v_h2_red, '-', ...
%         'Color', h2_light_blue, 'LineWidth', 2.5, 'DisplayName', 'H2: inert jet');
% end

if ~isempty(v_h2_blue_solid)
    plot(x_h2_blue_solid_shifted, v_h2_blue_solid, '-', ...
        'Color', h2_medium_blue, 'LineWidth', 2.5, 'DisplayName', 'H2: flame penetration (multi-kernel)');
end
% if ~isempty(v_h2_blue_dots)
%     plot(x_h2_blue_dots_shifted, v_h2_blue_dots, '.', ...
%         'Color', h2_medium_blue, 'MarkerSize', 15, 'DisplayName', 'H2: flame recession (multi-kernel)');
% end

if ~isempty(v_h2_black_solid)
    plot(x_h2_black_solid_shifted, v_h2_black_solid, '-', ...
        'Color', h2_dark_blue, 'LineWidth', 2.5, 'DisplayName', 'H2: flame penetration (single-kernel)');
end

% if ~isempty(v_h2_black_dots)
%     plot(x_h2_black_dots_shifted, v_h2_black_dots, '.', ...
%         'Color', h2_dark_blue, 'MarkerSize', 15, 'DisplayName', 'H2: flame recession (single-kernel)');
% end

% CH4 velocities
% if ~isempty(v_ch4_black)
%     plot(x_ch4_black_shifted, v_ch4_black, '-', ...
%         'Color', ch4_dark_orange, 'LineWidth', 2.5, 'DisplayName', 'CH4: inert jet');
% end

if ~isempty(v_ch4_red_solid)
    plot(x_ch4_red_solid_shifted, v_ch4_red_solid, '-', ...
        'Color', ch4_orange, 'LineWidth', 2.5, 'DisplayName', 'CH4: flame penetration (1060 K)');
end

% if ~isempty(v_ch4_red_dots)
%     plot(x_ch4_red_dots_shifted, v_ch4_red_dots, '.', ...
%         'Color', ch4_light_orange, 'MarkerSize', 15, 'DisplayName', 'CH4: flame recession');
% end

xlabel('Time after start of ignition events [ms]', 'FontSize', 13);
ylabel('Penentration rate [mm/ms]', 'FontSize', 13);
title('CH4 vs H2 Penetration Rate Comparison', 'FontSize', 14, 'FontWeight', 'bold');
xlim([0 0.75]);
ylim([0 150]);
set(gca, 'XTick', 0:0.25:0.75);
legend('Location', 'best', 'FontSize', 10, 'Box', 'on');
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
set(gcf, 'Color', 'w');   % only outside area → white
grid on;
hold off;

fprintf('\n=== Plotting complete ===\n');
fprintf('Figure 1: Original data (penetration vs time)\n');
fprintf('Figure 2: Smoothed velocity (dy/dx) - all lines shifted to start at x=0\n');
fprintf('Smoothing window: %d points\n', smooth_window);