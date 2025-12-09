% Compare CH4DDI and H2DDI data - Modified for new data structures
clear; clc; close all;

% Load both MAT files
file1 = 'C:\Users\jun-y\OneDrive - UNSW\CH4DDI\CH4DDI vs H2DDI CVCC\Mark2024_handled\Mark2024_Figure_extractedData_D_1_03ms_M.mat';
file2 = 'C:\Users\jun-y\OneDrive - UNSW\CH4DDI\Results\H2_jet_penetration.mat';  % Updated path for new H2 data

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

% Extract H2 data - NEW FORMAT with x_data and y_data
h2_data = struct();
if isfield(data_h2, 'x_data') && isfield(data_h2, 'y_data')
    % New format: x_data and y_data directly
    h2_data.x_black_dotted_data = data_h2.x_data;
    h2_data.y_black_dotted_data = data_h2.y_data;
    fprintf('H2 data loaded: %d points\n', length(data_h2.x_data));
else
    error('Cannot find x_data and y_data in H2 file');
end

% Initialize empty arrays for other H2 data (since new format only has jet penetration)
h2_data.x_green_data = [];
h2_data.y_green_data = [];
h2_data.x_blue_data = [];
h2_data.y_blue_data = [];

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

%% Figure 1: Jet penetration comparison only (no flame data)
figure('Position', [50, 50, 1000, 700]);
hold on; grid on; box on;

% Plot H2 jet penetration data (already starts from 0)
x_h2_jet = [];
y_h2_jet = [];
if isfield(h2_data, 'x_black_dotted_data') && isfield(h2_data, 'y_black_dotted_data')
    if ~isempty(h2_data.x_black_dotted_data)
        x_h2_jet = h2_data.x_black_dotted_data;
        y_h2_jet = h2_data.y_black_dotted_data;
        
        plot(x_h2_jet, y_h2_jet, '-', ...
            'Color', h2_light_blue, 'LineWidth', 2.5, 'DisplayName', 'H2: jet penetration');
    end
end

% Process CH4 inert jet data - find shift point and apply
x_ch4_jet_shifted = [];
y_ch4_jet_shifted = [];
if ~isempty(x_ch4_nonreacting) && ~isempty(y_ch4_nonreacting)
    % Find the last point where y = 0 from the start
    shift_idx = 1;
    for i = 1:length(y_ch4_nonreacting)
        if y_ch4_nonreacting(i) == 0
            shift_idx = i;  % Keep updating to find the LAST zero point
        elseif y_ch4_nonreacting(i) > 0
            break;  % Stop when we hit the first non-zero value
        end
    end
    
    % Shift time to make the point AFTER the last zero as t=0
    if shift_idx < length(y_ch4_nonreacting)
        shift_idx = shift_idx + 1;  % Move to first non-zero point
    end
    
    % Shift both time and penetration to start from (0,0)
    time_shift = x_ch4_nonreacting(shift_idx);
    y_shift = y_ch4_nonreacting(shift_idx);
    
    x_ch4_jet_shifted = x_ch4_nonreacting - time_shift;
    y_ch4_jet_shifted = y_ch4_nonreacting - y_shift;  % Also shift y-values to start from 0
    
    fprintf('CH4 jet shifted by %.3f ms and %.1f mm (from point %d)\n', time_shift, y_shift, shift_idx);
    
    plot(x_ch4_jet_shifted, y_ch4_jet_shifted, '-', ...
        'Color', ch4_dark_orange, 'LineWidth', 2.5, 'DisplayName', 'CH4: inert jet (shifted)');
end

xlabel('Time aSOI [ms]', 'FontSize', 13);
ylabel('Jet penetration [mm]', 'FontSize', 13);
title('H2 vs CH4 Jet Penetration Comparison', 'FontSize', 14, 'FontWeight', 'bold');
xlim([0 3.5]);
ylim([0 80]);
set(gca, 'XTick', 0:0.5:3.5);
set(gca, 'YTick', 0:20:80);
set(gcf, 'Color', 'w');
legend('Location', 'best', 'FontSize', 10, 'Box', 'on');
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
hold off;

%% Figure 2: Calculate and plot difference between curves
figure('Position', [100, 100, 1000, 700]);

% Calculate difference between H2 and CH4 jet penetration
if ~isempty(x_h2_jet) && ~isempty(x_ch4_jet_shifted)
    % Find common time range
    t_min = max(min(x_h2_jet), min(x_ch4_jet_shifted));
    t_max = min(max(x_h2_jet), max(x_ch4_jet_shifted));
    
    % Create common time vector
    t_common = linspace(t_min, t_max, 1000);
    
    % Interpolate both curves to common time points
    y_h2_interp = interp1(x_h2_jet, y_h2_jet, t_common, 'linear', 'extrap');
    y_ch4_interp = interp1(x_ch4_jet_shifted, y_ch4_jet_shifted, t_common, 'linear', 'extrap');
    
    % Calculate difference (H2 - CH4)
    y_difference = y_h2_interp - y_ch4_interp;
    
    subplot(2,1,1);
    hold on; grid on; box on;
    plot(x_h2_jet, y_h2_jet, '-', 'Color', h2_light_blue, 'LineWidth', 2.5, 'DisplayName', 'H2: jet penetration');
    plot(x_ch4_jet_shifted, y_ch4_jet_shifted, '-', 'Color', ch4_dark_orange, 'LineWidth', 2.5, 'DisplayName', 'CH4: inert jet (shifted)');
    xlabel('Time [ms]', 'FontSize', 12);
    ylabel('Jet penetration [mm]', 'FontSize', 12);
    title('H2 vs CH4 Jet Penetration (Aligned)', 'FontSize', 13, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 10);
    xlim([t_min, t_max]);
    set(gca, 'FontSize', 11);
    
    subplot(2,1,2);
    hold on; grid on; box on;
    plot(t_common, y_difference, '-', 'Color', [0.8 0.2 0.8], 'LineWidth', 2.5, 'DisplayName', 'Difference (H2 - CH4)');
    plot([t_min t_max], [0 0], 'k--', 'LineWidth', 1, 'DisplayName', 'Zero line');
    xlabel('Time [ms]', 'FontSize', 12);
    ylabel('Penetration difference [mm]', 'FontSize', 12);
    title('Difference in Jet Penetration (H2 - CH4)', 'FontSize', 13, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 10);
    xlim([t_min, t_max]);
    set(gca, 'FontSize', 11);
    
    % Print statistics
    fprintf('\n=== Penetration Difference Analysis ===\n');
    fprintf('Time range for comparison: %.3f to %.3f ms\n', t_min, t_max);
    fprintf('Mean difference (H2 - CH4): %.2f mm\n', mean(y_difference));
    fprintf('Max difference: %.2f mm at t = %.3f ms\n', max(y_difference), t_common(y_difference == max(y_difference)));
    fprintf('Min difference: %.2f mm at t = %.3f ms\n', min(y_difference), t_common(y_difference == min(y_difference)));
    fprintf('RMS difference: %.2f mm\n', sqrt(mean(y_difference.^2)));
    
    % Determine which fuel penetrates faster on average
    if mean(y_difference) > 0
        fprintf('Result: H2 penetrates %.2f mm further on average\n', abs(mean(y_difference)));
    else
        fprintf('Result: CH4 penetrates %.2f mm further on average\n', abs(mean(y_difference)));
    end
else
    fprintf('Cannot calculate difference - missing data\n');
end

set(gcf, 'Color', 'w');

%% Display data structure information
fprintf('\n=== Data Structure Summary ===\n');

% H2 data fields
fprintf('\nH2 Data Fields Found:\n');
fprintf('  - x_data: %s (jet penetration)\n', mat2str(size(data_h2.x_data)));
fprintf('  - y_data: %s (jet penetration)\n', mat2str(size(data_h2.y_data)));

% CH4 data fields
fprintf('\nCH4 Data Fields Found:\n');
ch4_fields = fieldnames(ch4_data);
for i = 1:length(ch4_fields)
    fprintf('  - %s\n', ch4_fields{i});
end

fprintf('\n=== Plotting complete ===\n');
fprintf('Figure 1: Original data (penetration vs time) - no smoothing\n');
fprintf('Figure 2: Raw velocity (dy/dx) - original data with no smoothing or shifting\n');
fprintf('Note: H2 data now contains only jet penetration data from new detection method\n');
