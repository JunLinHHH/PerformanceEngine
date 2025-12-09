% Compare CH4DDI and H2DDI data - Multi-case analysis with structured data
clear; clc; close all;

%% Define file paths for all cases
% CH4 files
ch4_files = {
    'C:\Users\jun-y\OneDrive - UNSW\CH4DDI\CH4DDI vs H2DDI CVCC\Mark2024_handled\Mark2024_Figure_extractedData_D_0_03ms_M.mat';
    'C:\Users\jun-y\OneDrive - UNSW\CH4DDI\CH4DDI vs H2DDI CVCC\Mark2024_handled\Mark2024_Figure_extractedData_D_1_03ms_M.mat';
    'C:\Users\jun-y\OneDrive - UNSW\CH4DDI\CH4DDI vs H2DDI CVCC\Mark2024_handled\Mark2024_Figure_extractedData_D_2_03ms_M.mat';
    'C:\Users\jun-y\OneDrive - UNSW\CH4DDI\CH4DDI vs H2DDI CVCC\Mark2024_handled\Mark2024_Figure_extractedData_D_3_03ms_M.mat'
};

% H2 files
h2_files = {
    'C:\Users\jun-y\OneDrive - UNSW\CH4DDI\CH4DDI vs H2DDI CVCC\Patrick2022_handled\H-0.07ms-D_Patrick2022_Figure10_extracted_data.mat';
    'C:\Users\jun-y\OneDrive - UNSW\CH4DDI\CH4DDI vs H2DDI CVCC\Patrick2022_handled\D-0.93ms-H_Patrick2022_Figure10_extracted_data.mat';
    'C:\Users\jun-y\OneDrive - UNSW\CH4DDI\CH4DDI vs H2DDI CVCC\Patrick2022_handled\D-1.93ms-H_Patrick2022_Figure10_extracted_data.mat';
    'C:\Users\jun-y\OneDrive - UNSW\CH4DDI\CH4DDI vs H2DDI CVCC\Patrick2022_handled\D-2.93ms-H_Patrick2022_Figure10_extracted_data.mat'
};

% Additional H2 jet penetration data
h2_jet_penetration_file = 'C:\Users\jun-y\OneDrive - UNSW\CH4DDI\CH4DDI vs H2DDI CVCC\Patrick2022_handled\H2_jet_penetration_Patrick2022Figure13.mat';

% Case names for legends
ch4_case_names = {
    'CH4: D-0.03ms-M';
    'CH4: D-1.03ms-M';
    'CH4: D-2.03ms-M';
    'CH4: D-3.03ms-M'
};

h2_case_names = {
    'H2: H-0.07ms-D';
    'H2: D-0.93ms-H';
    'H2: D-1.93ms-H';
    'H2: D-2.93ms-H'
};

%% Load and organize all data into structures
fprintf('=== Loading all data files ===\n');

% Initialize structured data
CH4_data = struct();
H2_data = struct();

% Load CH4 data
for i = 1:length(ch4_files)
    fprintf('Loading CH4 case %d: %s\n', i, ch4_case_names{i});
    try
        loaded_data = load(ch4_files{i});
        CH4_data.(['case_' num2str(i)]) = extract_data_from_structure(loaded_data, 'ch4');
        CH4_data.(['case_' num2str(i)]).case_name = ch4_case_names{i};
        CH4_data.(['case_' num2str(i)]).file_path = ch4_files{i};
        fprintf('  ✓ CH4 case %d loaded successfully\n', i);
    catch ME
        fprintf('  ✗ Error loading CH4 case %d: %s\n', i, ME.message);
        CH4_data.(['case_' num2str(i)]) = struct();
    end
end

% Load H2 data
for i = 1:length(h2_files)
    fprintf('Loading H2 case %d: %s\n', i, h2_case_names{i});
    try
        loaded_data = load(h2_files{i});
        H2_data.(['case_' num2str(i)]) = extract_data_from_structure(loaded_data, 'h2');
        H2_data.(['case_' num2str(i)]).case_name = h2_case_names{i};
        H2_data.(['case_' num2str(i)]).file_path = h2_files{i};
        fprintf('  ✓ H2 case %d loaded successfully\n', i);
    catch ME
        fprintf('  ✗ Error loading H2 case %d: %s\n', i, ME.message);
        H2_data.(['case_' num2str(i)]) = struct();
    end
end

% Load additional H2 jet penetration data
fprintf('\nLoading additional H2 jet penetration data...\n');
try
    h2_jet_data = load(h2_jet_penetration_file);
    fprintf('  ✓ H2 jet penetration data loaded successfully\n');
catch ME
    fprintf('  ✗ Error loading H2 jet penetration data: %s\n', ME.message);
    h2_jet_data = struct();
end

%% Define color schemes for multiple cases
% CH4: Different shades of orange/red
ch4_colors = [
    0.8, 0.2, 0.0;    % Dark red-orange (case 1)
    1.0, 0.4, 0.0;    % Orange (case 2)
    1.0, 0.6, 0.2;    % Light orange (case 3)
    1.0, 0.8, 0.4     % Very light orange (case 4)
];

% H2: Different shades of blue
h2_colors = [
    0.0, 0.2, 0.6;    % Dark blue (case 1)
    0.0, 0.4, 0.8;    % Medium blue (case 2)
    0.3, 0.6, 1.0;    % Light blue (case 3)
    0.6, 0.8, 1.0     % Very light blue (case 4)
];

% Line styles for different data types
line_styles = {'-', '-.', '--', ':'};

%% Figure 1: Original data comparison - All cases
figure('Position', [50, 50, 1200, 800]);
hold on; grid on; box on;

% Plot H2 Data for all cases
for i = 1:4
    case_name = ['case_' num2str(i)];
    if isfield(H2_data, case_name) && ~isempty(fieldnames(H2_data.(case_name)))
        current_h2 = H2_data.(case_name);
        
        % Inert jet penetration (black dotted equivalent)
        if ~isempty(current_h2.x_black_dotted) && ~isempty(current_h2.y_black_dotted)
            % Apply smoothing for x < 1.5
            x_h2_inert = current_h2.x_black_dotted;
            y_h2_inert = current_h2.y_black_dotted;
            
            idx_early = x_h2_inert < 1.5;
            y_h2_inert_smooth = y_h2_inert;
            if any(idx_early)
                y_early_smooth = medfilt1(y_h2_inert(idx_early), 5);
                y_early_smooth = smooth(y_early_smooth, 10, 'loess');
                y_h2_inert_smooth(idx_early) = y_early_smooth;
            end
            
            plot(x_h2_inert, y_h2_inert_smooth, line_styles{1}, ...
                'Color', h2_colors(i,:), 'LineWidth', 2.5, ...
                'DisplayName', [current_h2.case_name ': inert jet']);
        end
        
        % Flame penetration average (green equivalent)
        if ~isempty(current_h2.x_green) && ~isempty(current_h2.y_green)
            plot(current_h2.x_green, current_h2.y_green, line_styles{2}, ...
                'Color', h2_colors(i,:), 'LineWidth', 2.5, ...
                'DisplayName', [current_h2.case_name ': flame penetration avg']);
        end
        
        % Flame penetration (blue equivalent)
        if ~isempty(current_h2.x_blue) && ~isempty(current_h2.y_blue) && ~isscalar(current_h2.x_blue)
            plot(current_h2.x_blue, current_h2.y_blue, line_styles{3}, ...
                'Color', h2_colors(i,:), 'LineWidth', 2.5, ...
                'DisplayName', [current_h2.case_name ': flame penetration']);
        end
    end
end

% Plot CH4 Data for all cases
for i = 1:4
    case_name = ['case_' num2str(i)];
    if isfield(CH4_data, case_name) && ~isempty(fieldnames(CH4_data.(case_name)))
        current_ch4 = CH4_data.(case_name);
        
        % Inert jet
        if ~isempty(current_ch4.x_nonreacting) && ~isempty(current_ch4.y_nonreacting)
            plot(current_ch4.x_nonreacting, current_ch4.y_nonreacting, line_styles{1}, ...
                'Color', ch4_colors(i,:), 'LineWidth', 2.5, ...
                'DisplayName', [current_ch4.case_name ': inert jet']);
        end
        
        % Flame penetration average
        if ~isempty(current_ch4.x_penetration) && ~isempty(current_ch4.y_penetration)
            plot(current_ch4.x_penetration, current_ch4.y_penetration, line_styles{2}, ...
                'Color', ch4_colors(i,:), 'LineWidth', 2.5, ...
                'DisplayName', [current_ch4.case_name ': flame penetration avg']);
        end
        
        % Flame recession
        if ~isempty(current_ch4.x_recession) && ~isempty(current_ch4.y_recession)
            plot(current_ch4.x_recession, current_ch4.y_recession, '.', ...
                'Color', ch4_colors(i,:), 'MarkerSize', 15, ...
                'DisplayName', [current_ch4.case_name ': flame recession']);
        end
    end
end

xlabel('Time aSOI [ms]', 'FontSize', 13);
ylabel('Jet/flame penetration [mm]', 'FontSize', 13);
title('CH4 vs H2 Multi-Case Comparison', 'FontSize', 14, 'FontWeight', 'bold');
xlim([0 3.5]);
ylim([0 80]);
set(gca, 'XTick', 0:0.5:3.5);
set(gca, 'YTick', 0:20:80);
set(gcf, 'Color', 'w');
legend('Location', 'eastoutside', 'FontSize', 9, 'Box', 'on');
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
hold off;

%% Figure 2: Velocity comparison - All cases
figure('Position', [100, 100, 1200, 800]);
hold on; grid on; box on;

% Smoothing parameters
smooth_window = 10;
smooth_window_h2 = 50;
smooth_window_ch4 = 20;

% Manual adjustments for fine-tuning (all zero for now, but available for adjustment)
ch4_manual_adjustments = [-1, 0, 0, 0]; % Fine-tuning adjustments for each CH4 case
h2_manual_adjustments = [-0.63-0.12, -1.19, -2.3, -3.43];  % Fine-tuning adjustments for each H2 case

% Plot H2 velocities for all cases
for i = 1:4
    case_name = ['case_' num2str(i)];
    if isfield(H2_data, case_name) && ~isempty(fieldnames(H2_data.(case_name)))
        current_h2 = H2_data.(case_name);
        
        % Flame penetration average velocity (green equivalent)
        if ~isempty(current_h2.x_green) && ~isempty(current_h2.y_green)
            [v_h2, x_h2_shifted] = compute_case_velocity(current_h2.x_green, current_h2.y_green, smooth_window, h2_manual_adjustments(i), 'h2');
            if ~isempty(v_h2)
                plot(x_h2_shifted, v_h2, line_styles{1}, ...
                    'Color', h2_colors(i,:), 'LineWidth', 2.5, ...
                    'DisplayName', [current_h2.case_name]);
            end
        end
        
        % Flame penetration velocity (blue equivalent)
        if ~isempty(current_h2.x_blue) && ~isempty(current_h2.y_blue) && ~isscalar(current_h2.x_blue)
            [v_h2_blue, x_h2_blue_shifted] = compute_case_velocity(current_h2.x_blue, current_h2.y_blue, smooth_window, h2_manual_adjustments(i), 'h2');
            if ~isempty(v_h2_blue)
                plot(x_h2_blue_shifted, v_h2_blue, line_styles{1}, ...
                    'Color', h2_colors(i,:), 'LineWidth', 2.5, ...
                    'DisplayName', [current_h2.case_name]);
            end
        end
    end
end

% Plot additional H2 jet penetration velocity
if ~isempty(fieldnames(h2_jet_data))
    % Extract jet penetration data based on the actual structure (x_data, y_data)
    if isfield(h2_jet_data, 'x_data') && isfield(h2_jet_data, 'y_data')
        x_jet = h2_jet_data.x_data;
        y_jet = h2_jet_data.y_data;
        if ~isempty(x_jet) && ~isempty(y_jet)
            [v_h2_jet, x_h2_jet_shifted] = compute_case_velocity(x_jet, y_jet, smooth_window, 0, 'h2');
            if ~isempty(v_h2_jet)
                plot(x_h2_jet_shifted, v_h2_jet, '--', ...
                    'Color', [0.2, 0.2, 0.2], 'LineWidth', 3, ...
                    'DisplayName', 'H2: Jet Penetration (Patrick2022 Fig13)');
            end
        end
    end
end

% Plot CH4 velocities for all cases
for i = 1:4
    case_name = ['case_' num2str(i)];
    if isfield(CH4_data, case_name) && ~isempty(fieldnames(CH4_data.(case_name)))
        current_ch4 = CH4_data.(case_name);
        
        % Flame penetration velocity
        if ~isempty(current_ch4.x_penetration) && ~isempty(current_ch4.y_penetration)
            [v_ch4, x_ch4_shifted] = compute_case_velocity(current_ch4.x_penetration, current_ch4.y_penetration, smooth_window_ch4, ch4_manual_adjustments(i), 'ch4');
            if ~isempty(v_ch4)
                plot(x_ch4_shifted, v_ch4, line_styles{1}, ...
                    'Color', ch4_colors(i,:), 'LineWidth', 2.5, ...
                    'DisplayName', [current_ch4.case_name]);
            end
        end
        
        % CH4 non-reacting jet velocity
        if ~isempty(current_ch4.x_nonreacting) && ~isempty(current_ch4.y_nonreacting)
            [v_ch4_nonreacting, x_ch4_nonreacting_shifted] = compute_case_velocity(current_ch4.x_nonreacting, current_ch4.y_nonreacting, smooth_window_ch4, ch4_manual_adjustments(i), 'ch4');
            if ~isempty(v_ch4_nonreacting)
                plot(x_ch4_nonreacting_shifted, v_ch4_nonreacting, ':', ...
                    'Color', ch4_colors(i,:), 'LineWidth', 2, ...
                    'DisplayName', [current_ch4.case_name ': inert jet']);
            end
        end
    end
end

xlabel('Time after start of main fuel ignition events [ms]', 'FontSize', 13);
ylabel('Penetration rate [mm/ms]', 'FontSize', 13);
title('CH4 vs H2 Multi-Case Penetration Rate Comparison', 'FontSize', 14, 'FontWeight', 'bold');
xlim([0 1]);
ylim([5 90]);
set(gca, 'XTick', 0:0.25:2);
legend('Location', 'north', 'FontSize', 9, 'Box', 'on', 'NumColumns', 2);
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
set(gcf, 'Color', 'w');
grid on;
hold off;

%% Figure 3: Complete Penetration Rate Comparison (No shifting)
figure('Position', [150, 150, 1200, 800]);
hold on; grid on; box on;

% Plot H2 velocities for all cases (no shifting)
for i = 1:4
    case_name = ['case_' num2str(i)];
    if isfield(H2_data, case_name) && ~isempty(fieldnames(H2_data.(case_name)))
        current_h2 = H2_data.(case_name);
        
        % Flame penetration average velocity (green equivalent) - no shifting
        if ~isempty(current_h2.x_green) && ~isempty(current_h2.y_green)
            velocity_h2 = gradient(current_h2.y_green) ./ gradient(current_h2.x_green);
            velocity_h2 = movmean(velocity_h2, smooth_window, 'omitnan');
            plot(current_h2.x_green, velocity_h2, line_styles{1}, ...
                'Color', h2_colors(i,:), 'LineWidth', 2.5, ...
                'DisplayName', [current_h2.case_name]);
        end
        
        % Flame penetration velocity (blue equivalent) - no shifting
        if ~isempty(current_h2.x_blue) && ~isempty(current_h2.y_blue) && ~isscalar(current_h2.x_blue)
            velocity_h2_blue = gradient(current_h2.y_blue) ./ gradient(current_h2.x_blue);
            velocity_h2_blue = movmean(velocity_h2_blue, smooth_window, 'omitnan');
            plot(current_h2.x_blue, velocity_h2_blue, line_styles{1}, ...
                'Color', h2_colors(i,:), 'LineWidth', 2.5, ...
                'DisplayName', [current_h2.case_name]);
        end
    end
end

% Plot additional H2 jet penetration velocity (no shifting)
if ~isempty(fieldnames(h2_jet_data))
    % Extract jet penetration data based on the actual structure (x_data, y_data)
    if isfield(h2_jet_data, 'x_data') && isfield(h2_jet_data, 'y_data')
        x_jet = h2_jet_data.x_data;
        y_jet = h2_jet_data.y_data;
        if ~isempty(x_jet) && ~isempty(y_jet)
            % Compute velocity without any shifting
            velocity_jet = gradient(y_jet) ./ gradient(x_jet);
            velocity_jet = movmean(velocity_jet, smooth_window, 'omitnan');
            
            plot(x_jet, velocity_jet, '--', ...
                'Color', [0.2, 0.2, 0.2], 'LineWidth', 3, ...
                'DisplayName', 'H2: Jet Penetration (Patrick2022 Fig13)');
        end
    end
end

% Plot CH4 velocities for all cases (no shifting)
for i = 1:4
    case_name = ['case_' num2str(i)];
    if isfield(CH4_data, case_name) && ~isempty(fieldnames(CH4_data.(case_name)))
        current_ch4 = CH4_data.(case_name);
        
        % Flame penetration velocity - no shifting
        if ~isempty(current_ch4.x_penetration) && ~isempty(current_ch4.y_penetration)
            velocity_ch4 = gradient(current_ch4.y_penetration) ./ gradient(current_ch4.x_penetration);
            velocity_ch4 = movmean(velocity_ch4, smooth_window_ch4, 'omitnan');
            plot(current_ch4.x_penetration, velocity_ch4, line_styles{1}, ...
                'Color', ch4_colors(i,:), 'LineWidth', 2.5, ...
                'DisplayName', [current_ch4.case_name]);
        end
        
        % CH4 non-reacting jet velocity - no shifting
        if ~isempty(current_ch4.x_nonreacting) && ~isempty(current_ch4.y_nonreacting)
            velocity_ch4_nr = gradient(current_ch4.y_nonreacting) ./ gradient(current_ch4.x_nonreacting);
            velocity_ch4_nr = movmean(velocity_ch4_nr, smooth_window_ch4, 'omitnan');
            
            plot(current_ch4.x_nonreacting, velocity_ch4_nr, ':', ...
                'Color', ch4_colors(i,:), 'LineWidth', 2, ...
                'DisplayName', [current_ch4.case_name ': inert jet']);
        end
    end
end

xlabel('Time aSOI [ms]', 'FontSize', 13);
ylabel('Penetration rate [mm/ms]', 'FontSize', 13);
title('Complete Penetration Rate Comparison (Original Time Scale)', 'FontSize', 14, 'FontWeight', 'bold');
xlim([0 6.5]);
% ylim([0 150]);
set(gca, 'XTick', 0:0.5:6.5);
legend('Location', 'north', 'FontSize', 9, 'Box', 'on', 'NumColumns', 2);
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
set(gcf, 'Color', 'w');
grid on;
hold off;

%% Figure 4: CH4 Only Penetration Rate Analysis
figure('Position', [200, 200, 800, 400]);
hold on; grid on; box on;

% Plot CH4 velocities for all cases (no shifting)
for i = 1:4
    case_name = ['case_' num2str(i)];
    if isfield(CH4_data, case_name) && ~isempty(fieldnames(CH4_data.(case_name)))
        current_ch4 = CH4_data.(case_name);
        
        % Flame penetration velocity - no shifting
        if ~isempty(current_ch4.x_penetration) && ~isempty(current_ch4.y_penetration)
            velocity_ch4 = gradient(current_ch4.y_penetration) ./ gradient(current_ch4.x_penetration);
            velocity_ch4 = movmean(velocity_ch4, smooth_window_ch4, 'omitnan');
            plot(current_ch4.x_penetration, velocity_ch4, line_styles{1}, ...
                'Color', ch4_colors(i,:), 'LineWidth', 3, ...
                'DisplayName', [current_ch4.case_name]);
        end
        
        % CH4 non-reacting jet velocity - cut before peak
        if ~isempty(current_ch4.x_nonreacting) && ~isempty(current_ch4.y_nonreacting)
            velocity_ch4_nr = gradient(current_ch4.y_nonreacting) ./ gradient(current_ch4.x_nonreacting);
            velocity_ch4_nr = movmean(velocity_ch4_nr, smooth_window_ch4, 'omitnan');
            
            % Find the peak (maximum velocity)
            [~, peak_idx] = max(velocity_ch4_nr);
            
            % Cut everything before the peak: set to NaN
            velocity_ch4_nr_cut = velocity_ch4_nr;
            velocity_ch4_nr_cut(1:peak_idx-1) = NaN;
            
            plot(current_ch4.x_nonreacting, velocity_ch4_nr_cut, '--', ...
                'Color', ch4_colors(i,:), 'LineWidth', 2, ...
                'DisplayName', [current_ch4.case_name ': inert jet']);
        end
    end
end

xlabel('Time aSOI [ms]', 'FontSize', 13);
ylabel('Penetration rate [mm/ms]', 'FontSize', 13);
xlim([0 5]);
ylim([0 100]);
set(gca, 'XTick', 0:0.5:5);
set(gca, 'YTick', 0:20:100);
legend('Location', 'northeast', 'FontSize', 10, 'Box', 'on');
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
set(gcf, 'Color', 'w');
grid on;
hold off;

%% Figure 5: H2 Only Penetration Rate Analysis
figure('Position', [200, 200, 800, 400]);
hold on; grid on; box on;

% Time shifts for H2 jet penetration lines
h2_jet_time_shifts = [0, 0.93, 1.93, 2.93];

% Plot H2 velocities for all cases (no shifting)
for i = 1:4
    case_name = ['case_' num2str(i)];
    if isfield(H2_data, case_name) && ~isempty(fieldnames(H2_data.(case_name)))
        current_h2 = H2_data.(case_name);
        
        % Flame penetration average velocity (green equivalent) - no shifting
        if ~isempty(current_h2.x_green) && ~isempty(current_h2.y_green)
            velocity_h2 = gradient(current_h2.y_green) ./ gradient(current_h2.x_green);
            velocity_h2 = movmean(velocity_h2, smooth_window_h2, 'omitnan');
            plot(current_h2.x_green, velocity_h2, line_styles{1}, ...
                'Color', h2_colors(i,:), 'LineWidth', 3, ...
                'DisplayName', [current_h2.case_name]);
        end
        
        % Flame penetration velocity (blue equivalent) - no shifting
        if ~isempty(current_h2.x_blue) && ~isempty(current_h2.y_blue) && ~isscalar(current_h2.x_blue)
            velocity_h2_blue = gradient(current_h2.y_blue) ./ gradient(current_h2.x_blue);
            velocity_h2_blue = movmean(velocity_h2_blue, smooth_window_h2, 'omitnan');
            plot(current_h2.x_blue, velocity_h2_blue, line_styles{1}, ...
                'Color', h2_colors(i,:), 'LineWidth', 3, ...
                'DisplayName', [current_h2.case_name]);
        end
    end
end

% Plot 4 H2 jet penetration lines with specific time shifts and matching colors
if ~isempty(fieldnames(h2_jet_data))
    if isfield(h2_jet_data, 'x_data') && isfield(h2_jet_data, 'y_data')
        x_jet = h2_jet_data.x_data;
        y_jet = h2_jet_data.y_data;
        if ~isempty(x_jet) && ~isempty(y_jet)
            % Compute velocity without any shifting first
            velocity_jet = gradient(y_jet) ./ gradient(x_jet);
            velocity_jet = movmean(velocity_jet, smooth_window, 'omitnan');
            
            % Plot 4 jet penetration lines with different time shifts
            for i = 1:4
                x_jet_shifted = x_jet + h2_jet_time_shifts(i);
                plot(x_jet_shifted, velocity_jet, '--', ...
                    'Color', h2_colors(i,:), 'LineWidth', 2, ...
                    'DisplayName', sprintf('H2 Inert Jet: +%.2fms shift', h2_jet_time_shifts(i)));
            end
        end
    end
end

xlabel('Time aSOI [ms]', 'FontSize', 13);
ylabel('Penetration rate [mm/ms]', 'FontSize', 13);
%title('H2 Penetration Rate Analysis - All Cases (No Shift)', 'FontSize', 14, 'FontWeight', 'bold');
xlim([0 5]);
% ylim([0 150]);
set(gca, 'XTick', 0:0.5:6.5);
legend('Location', 'northeast', 'FontSize', 9, 'Box', 'on');
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
set(gcf, 'Color', 'w');
grid on;
hold off;

%% Figure 6: Combined CH4 and H2 Penetration Rate Analysis (Shifted to Zero)
figure('Position', [200, 200, 1200, 800]);
hold on; grid on; box on;

ch4_shift = [0.69, 1.82, 2.91, 4.33]; % Cut points for each H2 case
% Plot CH4 velocities for all cases (shifted to zero)
for i = 1:4
    case_name = ['case_' num2str(i)];
    if isfield(CH4_data, case_name) && ~isempty(fieldnames(CH4_data.(case_name)))
        current_ch4 = CH4_data.(case_name);
        
        % Flame penetration velocity - shifted to zero
        if ~isempty(current_ch4.x_penetration) && ~isempty(current_ch4.y_penetration)
            velocity_ch4 = gradient(current_ch4.y_penetration) ./ gradient(current_ch4.x_penetration);
            velocity_ch4 = movmean(velocity_ch4, smooth_window_ch4, 'omitnan');
            
            % Shift to zero: subtract the first x value
            x_ch4_shifted = current_ch4.x_penetration - ch4_shift(i);
            
            plot(x_ch4_shifted, velocity_ch4, line_styles{i}, ...
                'Color', ch4_colors(3,:), 'LineWidth', 3, ...
                'DisplayName', [current_ch4.case_name]);
        end
    end
end

% Plot H2 velocities for all cases (shifted to zero)
h2_cut_points = [0.69, 1.19, 2.3, 3.43]; % Cut points for each H2 case

for i = 1:4
    case_name = ['case_' num2str(i)];
    if isfield(H2_data, case_name) && ~isempty(fieldnames(H2_data.(case_name)))
        current_h2 = H2_data.(case_name);
        
        % Flame penetration average velocity (green equivalent) - shifted to zero
        if ~isempty(current_h2.x_green) && ~isempty(current_h2.y_green)
            velocity_h2 = gradient(current_h2.y_green) ./ gradient(current_h2.x_green);
            velocity_h2 = movmean(velocity_h2, smooth_window_h2, 'omitnan');
            
            % Cut the line: set points before cut_point to NaN
            x_h2_green_cut = current_h2.x_green;
            velocity_h2_cut = velocity_h2;
            velocity_h2_cut(x_h2_green_cut < h2_cut_points(i)) = NaN;
            
            % Shift to zero: subtract the cut point value
            x_h2_green_shifted = x_h2_green_cut - h2_cut_points(i);
            
            plot(x_h2_green_shifted, velocity_h2_cut, line_styles{i}, ...
                'Color', h2_colors(2,:), 'LineWidth', 3, ...
                'DisplayName', [current_h2.case_name]);
        end
        
        % Flame penetration velocity (blue equivalent) - shifted to zero
        if ~isempty(current_h2.x_blue) && ~isempty(current_h2.y_blue) && ~isscalar(current_h2.x_blue)
            velocity_h2_blue = gradient(current_h2.y_blue) ./ gradient(current_h2.x_blue);
            velocity_h2_blue = movmean(velocity_h2_blue, smooth_window_h2, 'omitnan');
            
            % Cut the line: set points before cut_point to NaN
            x_h2_blue_cut = current_h2.x_blue;
            velocity_h2_blue_cut = velocity_h2_blue;
            velocity_h2_blue_cut(x_h2_blue_cut < h2_cut_points(i)) = NaN;
            
            % Shift to zero: subtract the cut point value
            x_h2_blue_shifted = x_h2_blue_cut - h2_cut_points(i);
            
            plot(x_h2_blue_shifted, velocity_h2_blue_cut, line_styles{i}, ...
                'Color', h2_colors(2,:), 'LineWidth', 3, ...
                'DisplayName', [current_h2.case_name]);
        end
    end
end

xlabel('Time after main fuel ignition [ms]', 'FontSize', 13);
ylabel('Penetration rate [mm/ms]', 'FontSize', 13);
xlim([0 2]);
ylim([0 60]);
set(gca, 'XTick', 0:0.25:5);
set(gca, 'YTick', 0:10:100);
legend('Location', 'northeast', 'FontSize', 9, 'Box', 'on', 'NumColumns', 2);
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
set(gcf, 'Color', 'w');
grid on;
hold off;

%% Display summary information
fprintf('\n=== Multi-Case Data Loading Summary ===\n');

fprintf('\nCH4 Cases Loaded:\n');
for i = 1:4
    case_name = ['case_' num2str(i)];
    if isfield(CH4_data, case_name) && ~isempty(fieldnames(CH4_data.(case_name)))
        fprintf('  ✓ Case %d: %s\n', i, CH4_data.(case_name).case_name);
    else
        fprintf('  ✗ Case %d: Failed to load\n', i);
    end
end

fprintf('\nH2 Cases Loaded:\n');
for i = 1:4
    case_name = ['case_' num2str(i)];
    if isfield(H2_data, case_name) && ~isempty(fieldnames(H2_data.(case_name)))
        fprintf('  ✓ Case %d: %s\n', i, H2_data.(case_name).case_name);
    else
        fprintf('  ✗ Case %d: Failed to load\n', i);
    end
end

fprintf('\n=== Plotting Complete ===\n');
fprintf('Figure 1: Multi-case penetration vs time comparison\n');
fprintf('Figure 2: Multi-case penetration rate comparison (time-aligned)\n');
fprintf('Figure 3: Complete penetration rate comparison (original time scale)\n');
fprintf('Figure 4: CH4 only penetration rate analysis\n');
fprintf('Figure 5: H2 only penetration rate analysis\n');
fprintf('Total cases: %d CH4 + %d H2 = %d pairs\n', 4, 4, 4);

% Save the organized data structures for future use
fprintf('\nSaving organized data structures...\n');
save('CH4_H2_MultiCase_Data.mat', 'CH4_data', 'H2_data', 'h2_jet_data', 'ch4_case_names', 'h2_case_names');
fprintf('Data saved to: CH4_H2_MultiCase_Data.mat\n');

%% Velocity computation function (moved to end as requested)
function [velocity, x_shifted] = compute_case_velocity(x_data, y_data, smooth_window, manual_adjustment, data_type, reference_x_data)
    if isempty(x_data) || isempty(y_data)
        velocity = [];
        x_shifted = [];
        return;
    end
    
    % Determine zero reference based on data type
    if nargin >= 5 && strcmp(data_type, 'h2')
        % For H2: no automatic shifting, only manual adjustment
        zero_reference = 0;  % No automatic shift
    else
        % For CH4: use first non-NaN point in y_data to determine zero reference
        valid_indices = ~isnan(y_data);
        if ~any(valid_indices)
            velocity = [];
            x_shifted = [];
            return;
        end
        first_valid_idx = find(valid_indices, 1, 'first');
        zero_reference = x_data(first_valid_idx);
    end
    
    % Compute velocity
    velocity = gradient(y_data) ./ gradient(x_data);
    
    % Apply smoothing
    velocity = movmean(velocity, smooth_window, 'omitnan');
    
    % Apply shift based on zero reference + manual adjustment
    if nargin < 4
        manual_adjustment = 0;
    end
    
    if strcmp(data_type, 'h2')
        % For H2: keep original x_data, only apply manual adjustment
        x_shifted = x_data + manual_adjustment;
    else
        % For CH4: apply automatic shift + manual adjustment
        x_shifted = x_data - zero_reference + manual_adjustment;
    end
end

%% Function to extract data from loaded structure
function extracted_data = extract_data_from_structure(loaded_data, data_type)
    % data_type: 'ch4' or 'h2'
    extracted_data = struct();
    
    if strcmp(data_type, 'ch4')
        % Handle CH4 data structure
        if isfield(loaded_data, 'data_ch4_D_1_03ms_M')
            ch4_data = loaded_data.data_ch4_D_1_03ms_M;
        elseif ~isempty(fieldnames(loaded_data))
            fields = fieldnames(loaded_data);
            ch4_data = loaded_data.(fields{1});
        else
            error('Cannot find data structure in CH4 file');
        end
        
        % Extract CH4 specific fields
        extracted_data.x_nonreacting = [];
        extracted_data.y_nonreacting = [];
        extracted_data.x_penetration = [];
        extracted_data.y_penetration = [];
        extracted_data.x_recession = [];
        extracted_data.y_recession = [];
        
        if isfield(ch4_data, 'NonReactingJetPenetration') && isstruct(ch4_data.NonReactingJetPenetration)
            fields_nr = fieldnames(ch4_data.NonReactingJetPenetration);
            if length(fields_nr) >= 2
                extracted_data.x_nonreacting = ch4_data.NonReactingJetPenetration.(fields_nr{1});
                extracted_data.y_nonreacting = ch4_data.NonReactingJetPenetration.(fields_nr{2});
            end
        end
        
        if isfield(ch4_data, 'Penetration_average') && isstruct(ch4_data.Penetration_average)
            fields_pen = fieldnames(ch4_data.Penetration_average);
            if length(fields_pen) >= 2
                extracted_data.x_penetration = ch4_data.Penetration_average.(fields_pen{1});
                extracted_data.y_penetration = ch4_data.Penetration_average.(fields_pen{2});
            end
        end
        
        if isfield(ch4_data, 'Recession_average') && isstruct(ch4_data.Recession_average)
            fields_rec = fieldnames(ch4_data.Recession_average);
            if length(fields_rec) >= 2
                extracted_data.x_recession = ch4_data.Recession_average.(fields_rec{1});
                extracted_data.y_recession = ch4_data.Recession_average.(fields_rec{2});
            end
        end
        
    elseif strcmp(data_type, 'h2')
        % Handle H2 data structure
        if isfield(loaded_data, 'x_black_dotted_data')
            h2_data = loaded_data;
        elseif isfield(loaded_data, 'data_h2')
            h2_data = loaded_data.data_h2;
        elseif ~isempty(fieldnames(loaded_data))
            fields = fieldnames(loaded_data);
            if isstruct(loaded_data.(fields{1}))
                h2_data = loaded_data.(fields{1});
            else
                h2_data = loaded_data;
            end
        else
            h2_data = loaded_data;
        end
        
        % Extract H2 specific fields
        extracted_data.x_black_dotted = [];
        extracted_data.y_black_dotted = [];
        extracted_data.x_green = [];
        extracted_data.y_green = [];
        extracted_data.x_blue = [];
        extracted_data.y_blue = [];
        extracted_data.x_circle = [];
        extracted_data.y_circle = [];
        
        if isfield(h2_data, 'x_black_dotted_data') && isfield(h2_data, 'y_black_dotted_data')
            extracted_data.x_black_dotted = h2_data.x_black_dotted_data;
            extracted_data.y_black_dotted = h2_data.y_black_dotted_data;
        end
        
        if isfield(h2_data, 'x_green_data') && isfield(h2_data, 'y_green_data')
            extracted_data.x_green = h2_data.x_green_data;
            extracted_data.y_green = h2_data.y_green_data;
        end
        
        if isfield(h2_data, 'x_blue_data') && isfield(h2_data, 'y_blue_data')
            extracted_data.x_blue = h2_data.x_blue_data;
            extracted_data.y_blue = h2_data.y_blue_data;
        end
        
        if isfield(h2_data, 'x_circle_data') && isfield(h2_data, 'y_circle_data')
            extracted_data.x_circle = h2_data.x_circle_data;
            extracted_data.y_circle = h2_data.y_circle_data;
        end
    end
end

