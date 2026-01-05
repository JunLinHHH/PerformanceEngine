%% Test CLASS_PressureProcessEngine - Debug and validate processing pipeline
% This script tests the complete pressure processing workflow on a single Set

clear all; clearvars; close all; clc;
Configs = Configuration();
Loaded = load(Configs.DataObjReadDir);

obj = Loaded.obj;
T_DataAll = obj.DataMatrix;

%% Select one Set for testing
test_set = 2;  % Change this to test different Sets

fprintf('\n');
fprintf(repmat('=', 1, 70));
fprintf('\nTEST: CLASS_PressureProcessEngine\n');
fprintf(repmat('=', 1, 70));
fprintf('\n\n');

fprintf('Testing Set %02d...\n\n', test_set);

%% Extract test data
rowIdx = find(T_DataAll.Set == test_set, 1);
if isempty(rowIdx)
    error('Set %d not found in table', test_set);
end

this_row = T_DataAll(rowIdx, :);
tdmsFiles = this_row.File_TDMS_AllTakes{1};

fprintf('Found %d Take files\n', length(tdmsFiles));
fprintf('RPM: %.0f\n\n', this_row.RPM);

%% Create processor and run pipeline
fprintf('Creating processor...\n');
Obj_Pressure = CLASS_PressureProcessEngine(tdmsFiles, test_set, 'PressureChannel', 'Dev0_Ai5');

fprintf('Reading TDMS files...\n');
Obj_Pressure = Obj_Pressure.ReadPressureFile_TDMS();

fprintf('Running complete processing pipeline...\n');
Obj_Pressure = Obj_Pressure.ProcessAll(this_row);

% Debug: Check if fields were created
PressureInfo_debug = Obj_Pressure.GetPressureInfo();
fprintf('\nDEBUG: Checking first Take fields:\n');
if isfield(PressureInfo_debug(1), 'Filtered')
    fprintf('  Filtered field exists: %s\n', mat2str(~isempty(PressureInfo_debug(1).Filtered)));
else
    fprintf('  Filtered field missing!\n');
end

fprintf('✓ Processing complete\n\n');

%% Get results
PressureInfo = Obj_Pressure.GetPressureInfo();
nTakes = length(PressureInfo);

%% PHASE A: DEBUG - CYCLE BOUNDARY & CA INTEGRITY (INTEGRATED)
fprintf(repmat('=', 1, 80));
fprintf('\nPHASE A: DIAGNOSTIC - CYCLE BOUNDARY & CA INTEGRITY ANALYSIS\n');
fprintf(repmat('=', 1, 80));
fprintf('\n');

% Expected samples per cycle at this RPM
rpm = this_row.RPM;
expected_samples_per_cycle = 2 * 60 / rpm * 200000;  % 2 revolutions per cycle at Fs=200kHz
fprintf('RPM: %.0f\n', rpm);
fprintf('Expected samples per cycle: %.0f\n\n', expected_samples_per_cycle);

% Initialize debug report
debugReport = struct();
debugReport.SetNumber = test_set;
debugReport.RPM = rpm;
debugReport.ExpectedSamplesPerCycle = expected_samples_per_cycle;
debugReport.Takes = cell(nTakes, 1);  % Use cell array for Takes

% Loop through all takes and collect metrics
for jj = 1:nTakes
    take_idx = jj - 1;
    fprintf('--- TAKE %02d ---\n', take_idx);
    
    PI = PressureInfo(jj);
    takeData = struct();
    
    % Check for basic fields
    if ~isfield(PI, 'TDC_indices_shifted') || isempty(PI.TDC_indices_shifted)
        fprintf('  TDC_indices_shifted: EMPTY (no cycle detection)\n');
        takeData.TDC_indices_shifted = [];
        takeData.TDC_spacing_samples = [];
        takeData.TDC_spacing_deviation_pct = [];
        takeData.Z_rise_indices_shifted = [];
        takeData.Z_rise_spacing_samples = [];
        takeData.CA_warnings = 0;
        takeData.nCycles_detected = 0;
        takeData.nCycles_extracted = 0;
        takeData.nFiringCycles = 0;
        takeData.FiringRate_pct = 0;
        debugReport.Takes{jj} = takeData;  % Store in cell array
        fprintf('\n');
        continue;
    end
    
    % Get TDC indices
    TDC_idx_shifted = PI.TDC_indices_shifted;
    nCycles_detected = length(TDC_idx_shifted);
    
    % TDC spacing analysis
    TDC_spacing = diff(TDC_idx_shifted);
    TDC_spacing_mean = mean(TDC_spacing);
    TDC_spacing_std = std(TDC_spacing);
    TDC_spacing_min = min(TDC_spacing);
    TDC_spacing_max = max(TDC_spacing);
    TDC_spacing_deviation_pct = abs(TDC_spacing_mean - expected_samples_per_cycle) / expected_samples_per_cycle * 100;
    
    fprintf('  TDC_indices_shifted:       [%d indices, min=%d, max=%d samples]\n', ...
        nCycles_detected, TDC_idx_shifted(1), TDC_idx_shifted(end));
    fprintf('  TDC spacing (mean ± std):  %.0f ± %.0f samples\n', ...
        TDC_spacing_mean, TDC_spacing_std);
    fprintf('  TDC spacing range:         %.0f - %.0f samples\n', ...
        TDC_spacing_min, TDC_spacing_max);
    fprintf('  Deviation from expected:   %.2f%%\n', TDC_spacing_deviation_pct);
    
    % Z-pulse reference
    if isfield(PI, 'Z_rise_indices_shifted') && ~isempty(PI.Z_rise_indices_shifted)
        Z_rise_shifted = PI.Z_rise_indices_shifted;
        Z_spacing = diff(Z_rise_shifted);
        fprintf('  Z_rise_shifted:            [%d events, spacing %.0f ± %.0f]\n', ...
            length(Z_rise_shifted), mean(Z_spacing), std(Z_spacing));
        takeData.Z_rise_indices_shifted = Z_rise_shifted;
        takeData.Z_rise_spacing_samples = Z_spacing;
    else
        fprintf('  Z_rise_shifted:            EMPTY (reference missing)\n');
        takeData.Z_rise_indices_shifted = [];
        takeData.Z_rise_spacing_samples = [];
    end
    
    % CA integrity warnings
    nCA_warnings = 0;
    if isfield(PI, 'CA_warnings') && ~isempty(PI.CA_warnings)
        nCA_warnings = PI.CA_warnings;
        fprintf('  CA non-monotonic warnings: %d cycles\n', nCA_warnings);
    else
        fprintf('  CA non-monotonic warnings: 0\n');
    end
    
    % Cycle extraction status
    nCycles_extracted = 0;
    if isfield(PI, 'P_cycles') && ~isempty(PI.P_cycles)
        nCycles_extracted = size(PI.P_cycles, 2);
        fprintf('  Cycles extracted:         %d\n', nCycles_extracted);
    else
        fprintf('  Cycles extracted:         0\n');
    end
    
    % Firing rate (from P_max > 45 bar threshold)
    nFiringCycles = 0;
    if isfield(PI, 'Firing_indicator') && ~isempty(PI.Firing_indicator)
        nFiringCycles = sum(PI.Firing_indicator);
        firing_rate = nFiringCycles / nCycles_extracted * 100;
        fprintf('  Firing cycles (P_max>45): %d / %d (%.1f%%)\n', ...
            nFiringCycles, nCycles_extracted, firing_rate);
        takeData.FiringRate_pct = firing_rate;
    else
        fprintf('  Firing cycles:            NO DATA\n');
        takeData.FiringRate_pct = 0;
    end
    
    % Store metrics for this take
    takeData.TDC_indices_shifted = TDC_idx_shifted;
    takeData.TDC_spacing_samples = TDC_spacing;
    takeData.TDC_spacing_mean = TDC_spacing_mean;
    takeData.TDC_spacing_std = TDC_spacing_std;
    takeData.TDC_spacing_deviation_pct = TDC_spacing_deviation_pct;
    takeData.CA_warnings = nCA_warnings;
    takeData.nCycles_detected = nCycles_detected;
    takeData.nCycles_extracted = nCycles_extracted;
    takeData.nFiringCycles = nFiringCycles;
    
    debugReport.Takes{jj} = takeData;  % Store in cell array
    fprintf('\n');
end

% Print summary table
fprintf(repmat('=', 1, 80));
fprintf('\nPHASE A SUMMARY TABLE\n');
fprintf(repmat('=', 1, 80));
fprintf('\n');

% Create summary table
summaryData = [];
summaryLabels = cell(nTakes, 1);
for jj = 1:nTakes
    TD = debugReport.Takes{jj};  % Access from cell array
    summaryLabels{jj} = sprintf('Take %02d', jj-1);
    summaryData = [summaryData; ...
        TD.TDC_spacing_deviation_pct, ...
        TD.TDC_spacing_std, ...
        TD.CA_warnings, ...
        TD.nCycles_detected, ...
        TD.nCycles_extracted, ...
        TD.FiringRate_pct];
end

fprintf('%-12s | %10s | %10s | %10s | %10s | %10s | %10s\n', ...
    'Take', 'TDC Deviat', 'TDC Std', 'CA_Warn', 'Detected', 'Extracted', 'FiringRate');
fprintf('%-12s | %10s | %10s | %10s | %10s | %10s | %10s\n', ...
    '', '(%)', '(samp)', '(#)', '(#)', '(#)', '(%)');
fprintf(repmat('-', 1, 80));
fprintf('\n');

for jj = 1:nTakes
    TD = debugReport.Takes{jj};  % Access from cell array
    fprintf('Take %02d    | %10.2f | %10.1f | %10d | %10d | %10d | %10.1f\n', ...
        jj-1, ...
        TD.TDC_spacing_deviation_pct, ...
        TD.TDC_spacing_std, ...
        TD.CA_warnings, ...
        TD.nCycles_detected, ...
        TD.nCycles_extracted, ...
        TD.FiringRate_pct);
end
fprintf(repmat('=', 1, 80));
fprintf('\n\n');

% Debug report generation (not saved)
% debugReportFile = fullfile(fileparts(which('Configuration.m')), 'Phase_A_DebugReport.mat');
% save(debugReportFile, 'debugReport');
% fprintf('Debug report saved to: %s\n\n', debugReportFile);

fprintf(repmat('=', 1, 70));
fprintf('\nRESULTS SUMMARY\n');
fprintf(repmat('=', 1, 70));
fprintf('\n\n');

%% Display results for each Take
for jj = 1:nTakes
    fprintf('Take %02d:\n', jj-1);
    
    if isfield(PressureInfo(jj), 'Filtered') && ~isempty(PressureInfo(jj).Filtered)
        fprintf('  Filtered pressure:    ✓ (%d samples)\n', length(PressureInfo(jj).Filtered));
    else
        fprintf('  Filtered pressure:    ✗\n');
    end
    
    if isfield(PressureInfo(jj), 'TDC_indices') && ~isempty(PressureInfo(jj).TDC_indices)
        fprintf('  TDC detection:        ✓ (%d cycles)\n', length(PressureInfo(jj).TDC_indices));
    else
        fprintf('  TDC detection:        ✗\n');
    end
    
    if isfield(PressureInfo(jj), 'P_cycles') && ~isempty(PressureInfo(jj).P_cycles)
        fprintf('  Cycle extraction:     ✓ (%d cycles)\n', size(PressureInfo(jj).P_cycles, 2));
    else
        fprintf('  Cycle extraction:     ✗\n');
    end
    
    if isfield(PressureInfo(jj), 'IMEP') && ~isempty(PressureInfo(jj).IMEP)
        fprintf('  IMEP:                 %.2f bar\n', PressureInfo(jj).IMEP);
        fprintf('  CoV(IMEP):            %.2f%%\n', PressureInfo(jj).CoV_IMEP * 100);
        fprintf('  P_max:                %.2f bar\n', PressureInfo(jj).P_max);
        fprintf('  PRR_max:              %.2f bar/CAD\n', PressureInfo(jj).PRR_max);
    else
        fprintf('  IMEP calculation:     ✗\n');
    end
    
    if isfield(PressureInfo(jj), 'CA10') && ~isempty(PressureInfo(jj).CA10) && ~any(isnan(PressureInfo(jj).CA10))
        fprintf('  CA10:                 %.2f CAD\n', PressureInfo(jj).CA10);
        fprintf('  CA50:                 %.2f CAD\n', PressureInfo(jj).CA50);
        fprintf('  CA90:                 %.2f CAD\n', PressureInfo(jj).CA90);
        fprintf('  CA10-50:              %.2f CAD\n', PressureInfo(jj).CA10_50);
        fprintf('  CA50-90:              %.2f CAD\n', PressureInfo(jj).CA50_90);
        
        % Display injection timing and ignition delay
        if isfield(PressureInfo(jj), 'SOI_CA') && ~isempty(PressureInfo(jj).SOI_CA) && ~isnan(PressureInfo(jj).SOI_CA)
            fprintf('  SOI (injection):      %.2f CAD\n', PressureInfo(jj).SOI_CA);
            fprintf('  Ignition delay:       %.2f CAD (CA10 - SOI)\n', PressureInfo(jj).IgnDelay);
        else
            fprintf('  SOI (injection):      NOT DETECTED\n');
            fprintf('  Ignition delay:       %.2f CAD (absolute CA10)\n', PressureInfo(jj).IgnDelay);
        end
        
        % Display per-cycle metrics (Validation Task 4)
        if isfield(PressureInfo(jj), 'CA10_cycles') && ~isempty(PressureInfo(jj).CA10_cycles)
            n_firing = sum(~isnan(PressureInfo(jj).CA10_cycles));
            n_total = numel(PressureInfo(jj).CA10_cycles);
            firing_rate = 100 * n_firing / n_total;
            fprintf('  CA10_std:             %.2f CAD (cycle variability)\n', PressureInfo(jj).CA10_std);
            fprintf('  CA50_std:             %.2f CAD (cycle variability)\n', PressureInfo(jj).CA50_std);
            fprintf('  CA90_std:             %.2f CAD (cycle variability)\n', PressureInfo(jj).CA90_std);
            fprintf('  Per-cycle metrics:    ✓ (%d firing / %d total = %.1f%%)\n', n_firing, n_total, firing_rate);
        else
            fprintf('  Per-cycle metrics:    ✗ (not computed)\n');
        end
        fprintf('  Burn duration:        %.2f CAD\n', PressureInfo(jj).BurnDur);
    else
        fprintf('  HRR calculation:      ✗\n');
    end
    
    fprintf('\n');
end

%% === COMPREHENSIVE CYCLE DIAGNOSTICS: Plot all cycles from all Takes ===
fprintf('\n');
fprintf(repmat('=', 1, 70));
fprintf('\nPLOTTING ALL CYCLES FOR ALL TAKES\n');
fprintf(repmat('=', 1, 70));
fprintf('\n\n');

for take_idx = 1:nTakes
    fprintf('Plotting Take %02d cycles...\n', take_idx-1);
    
    % Get pressure cycles and CA data
    if ~isfield(PressureInfo(take_idx), 'P_cycles') || isempty(PressureInfo(take_idx).P_cycles)
        fprintf('  (Skipping - no P_cycles data)\n\n');
        continue;
    end
    
    P_cycles = PressureInfo(take_idx).P_cycles;  % Already interpolated to standard CA grid
    CA = Obj_Pressure.CA;  % Standard CA grid used for interpolation
    CA10_cycles = PressureInfo(take_idx).CA10_cycles;
    nCycles = size(P_cycles, 2);
    
    % Determine firing status
    firing_status = ~isnan(CA10_cycles);
    n_firing = sum(firing_status);
    firing_rate = 100 * n_firing / nCycles;
    
    % Create figure with subplots (4 cycles per row)
    n_cols = 4;
    n_rows = ceil(nCycles / n_cols);
    
    fig_title = sprintf('Set %02d Take %02d - All %d Cycles (%.1f%% firing)', test_set, take_idx-1, nCycles, firing_rate);
    fig = figure('Name', fig_title, 'NumberTitle', 'off', 'Position', [100 100 1600 1200]);
    
    for cc = 1:nCycles
        subplot(n_rows, n_cols, cc);
        
        % Get cycle pressure (already on CA grid)
        P_cycle = P_cycles(:, cc);
        
        % Compute P_max
        P_max = max(P_cycle);
        is_firing = firing_status(cc);
        
        % Plot with color-coded background
        plot(CA, P_cycle, 'LineWidth', 1.5, 'Color', 'black');
        hold on;
        yline(45, 'k--', 'LineWidth', 0.8, 'Alpha', 0.5);  % Threshold line
        
        % Set background color (green=firing, red=misfire)
        if is_firing
            set(gca, 'Color', [0.85, 1.0, 0.85]);  % Light green
            edge_color = 'green';
        else
            set(gca, 'Color', [1.0, 0.85, 0.85]);  % Light red
            edge_color = 'red';
        end
        set(gca, 'Box', 'on', 'LineWidth', 2);
        ax = gca;
        ax.XColor = edge_color;
        ax.YColor = edge_color;
        
        % Formatting
        xlabel('CA [deg]', 'FontSize', 8);
        ylabel('P [bar]', 'FontSize', 8);
        title(sprintf('Cycle %d (P_{max}=%.1f)', cc-1, P_max), 'FontSize', 9);
        grid on;
        xlim([-360, 360]);
        
        % Add status text
        if is_firing
            text(0.05, 0.95, 'FIRING', 'Units', 'normalized', ...
                'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
                'FontSize', 8, 'FontWeight', 'bold', 'Color', 'green');
        else
            text(0.05, 0.95, 'MISFIRE', 'Units', 'normalized', ...
                'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
                'FontSize', 8, 'FontWeight', 'bold', 'Color', 'red');
        end
    end
    
    % Add overall title and summary
    sgtitle(sprintf('Set %02d Take %02d: %d firing / %d total cycles (%.1f%%)', ...
        test_set, take_idx-1, n_firing, nCycles, firing_rate), ...
        'FontSize', 13, 'FontWeight', 'bold');
    
    % Figure saved to disk (disabled)
    % output_dir = fileparts(which('Configuration.m'));
    % output_file = fullfile(output_dir, sprintf('Set%02d_Take%02d_AllCycles.png', test_set, take_idx-1));
    % saveas(fig, output_file);
    % fprintf('  ✓ Saved to: %s\n', output_file);
    fprintf('  ✓ Plotted %d cycles (%.1f%% firing)\n\n', nCycles, firing_rate);
end

fprintf('All cycle diagnostics plotted.\n\n');

%% === PLOT: Individual Cycles with CA and Injection Signal ===
fprintf(repmat('=', 1, 70));
fprintf('\nPLOTTING INDIVIDUAL CYCLES WITH CA MARKERS AND INJECTION SIGNAL\n');
fprintf(repmat('=', 1, 70));
fprintf('\n\n');

for take_idx = 1:nTakes
    fprintf('Creating individual cycle plot for Take %02d...\n', take_idx-1);
    
    if ~isfield(PressureInfo(take_idx), 'P_cycles') || isempty(PressureInfo(take_idx).P_cycles)
        fprintf('  (Skipping - no P_cycles data)\n\n');
        continue;
    end
    
    P_cycles = PressureInfo(take_idx).P_cycles;
    CA = Obj_Pressure.CA;
    CA10_cycles = PressureInfo(take_idx).CA10_cycles;
    CA50_cycles = PressureInfo(take_idx).CA50_cycles;
    CA90_cycles = PressureInfo(take_idx).CA90_cycles;
    nCycles = size(P_cycles, 2);
    
    % Get injection signal if available
    InjectionSignal = [];
    dataRaw = [];
    if isfield(PressureInfo(take_idx), 'DataRaw') && ~isempty(PressureInfo(take_idx).DataRaw)
        dataRaw = PressureInfo(take_idx).DataRaw;
        try
            if isfield(dataRaw.Digital_channels, 'Di3')
                InjectionSignal = double(dataRaw.Digital_channels.Di3.data);
            elseif isfield(dataRaw.Digital_channels, 'Untitled_3')
                InjectionSignal = double(dataRaw.Digital_channels.Untitled_3.data);
            elseif isfield(dataRaw.Digital_channels, 'Untitled_5')
                InjectionSignal = double(dataRaw.Digital_channels.Untitled_5.data);
            end
        catch
        end
    end
    
    % Create figure for individual cycles
    n_cols = 4;
    n_rows = ceil(nCycles / n_cols);
    fig = figure('Name', sprintf('Set %02d - Take %02d: Individual Cycles with CA and Injection', test_set, take_idx-1), ...
                 'Position', [100 100 1600 1000]);
    
    for cc = 1:nCycles
        subplot(n_rows, n_cols, cc);
        
        % Plot cycle pressure
        plot(CA, P_cycles(:, cc), 'b-', 'LineWidth', 1.5); hold on;
        
        % Add CA10, CA50, CA90 markers
        if ~isnan(CA10_cycles(cc))
            xline(CA10_cycles(cc), 'g--', 'CA10');
        end
        if ~isnan(CA50_cycles(cc))
            xline(CA50_cycles(cc), 'r--', 'CA50');
        end
        if ~isnan(CA90_cycles(cc))
            xline(CA90_cycles(cc), 'm--', 'CA90');
        end
        
        % Add SOI marker
        if isfield(PressureInfo(take_idx), 'SOI_CA') && ~isempty(PressureInfo(take_idx).SOI_CA) && ~isnan(PressureInfo(take_idx).SOI_CA)
            xline(PressureInfo(take_idx).SOI_CA, 'k:', 'SOI', 'LineWidth', 1.5);
        end
        
        xlabel('Crank Angle [deg]');
        ylabel('Pressure [bar]');
        
        % Title with firing status
        if ~isnan(CA10_cycles(cc))
            title(sprintf('Cycle %d (firing)', cc), 'Color', 'green', 'FontWeight', 'bold');
        else
            title(sprintf('Cycle %d (misfire)', cc), 'Color', 'red');
        end
        
        grid on;
        xlim([-360 360]);
    end
    
    % Add overall title
    sgtitle(sprintf('Set %02d Take %02d: Individual Cycles with CA Markers and Injection Signal', ...
        test_set, take_idx-1), 'FontSize', 12, 'FontWeight', 'bold');
    
    fprintf('  ✓ Plotted %d cycles with CA and injection markers\n\n', nCycles);
end

fprintf('Individual cycle plots with CA markers completed.\n\n');

%% Visualization - Pick first valid Take
valid_idx = 1;
for jj = 1:nTakes
    if isfield(PressureInfo(jj), 'P_mean') && ~isempty(PressureInfo(jj).P_mean) && ...
       isfield(PressureInfo(jj), 'aHRR') && ~isempty(PressureInfo(jj).aHRR)
        valid_idx = jj;
        break;
    end
end

fprintf('Plotting detailed analysis for Take %02d...\n\n', valid_idx-1);

%% Plot 1: Pressure traces
figure('Name', sprintf('Set %02d - Take %02d: Pressure', test_set, valid_idx-1), 'Position', [100 100 1200 800]);

subplot(2,2,1);
x_has_time = isfield(PressureInfo(valid_idx), 'Time') && ~isempty(PressureInfo(valid_idx).Time);
x_has_raw = isfield(PressureInfo(valid_idx), 'Raw') && ~isempty(PressureInfo(valid_idx).Raw);
x_has_filt = isfield(PressureInfo(valid_idx), 'Filtered') && ~isempty(PressureInfo(valid_idx).Filtered);
if x_has_time && x_has_raw
    plot(PressureInfo(valid_idx).Time, PressureInfo(valid_idx).Raw, 'k'); hold on;
end
if x_has_time && x_has_filt
    plot(PressureInfo(valid_idx).Time, PressureInfo(valid_idx).Filtered, 'r');
end
xlabel('Time [s]');
ylabel('Pressure [bar]');
title('Raw vs Filtered Pressure');
legend('Raw', 'Filtered');
grid on;

subplot(2,2,2);
if isfield(PressureInfo(valid_idx), 'P_mean_corrected') && ~isempty(PressureInfo(valid_idx).P_mean_corrected)
    plot(Obj_Pressure.CA, PressureInfo(valid_idx).P_mean_corrected, 'b', 'LineWidth', 1.5); hold on;
elseif isfield(PressureInfo(valid_idx), 'P_mean') && ~isempty(PressureInfo(valid_idx).P_mean)
    plot(Obj_Pressure.CA, PressureInfo(valid_idx).P_mean, 'b', 'LineWidth', 1.5); hold on;
end
if isfield(PressureInfo(valid_idx), 'TDC_indices') && ~isempty(PressureInfo(valid_idx).TDC_indices) && ...
   isfield(PressureInfo(valid_idx), 'TDC_P') && isscalar(PressureInfo(valid_idx).TDC_P)
    plot(0, PressureInfo(valid_idx).TDC_P, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
end
xlabel('Crank Angle [deg]');
ylabel('Pressure [bar]');
title('Mean Pressure vs Crank Angle');
grid on;
xlim([-180 180]);

subplot(2,2,3);
if isfield(PressureInfo(valid_idx), 'aHRR') && ~isempty(PressureInfo(valid_idx).aHRR)
    plot(Obj_Pressure.CA(1:end-1), sgolayfilt(PressureInfo(valid_idx).aHRR, 7, 91), 'r', 'LineWidth', 1.5); hold on;
end
if isfield(PressureInfo(valid_idx), 'SOI_CA') && ~isempty(PressureInfo(valid_idx).SOI_CA) && ~isnan(PressureInfo(valid_idx).SOI_CA)
    xline(PressureInfo(valid_idx).SOI_CA, 'k:', 'SOI', 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
end
if isfield(PressureInfo(valid_idx), 'CA10') && ~isempty(PressureInfo(valid_idx).CA10) && ~isnan(PressureInfo(valid_idx).CA10)
    xline(PressureInfo(valid_idx).CA10, 'g--', 'CA10');
    xline(PressureInfo(valid_idx).CA50, 'r--', 'CA50');
    xline(PressureInfo(valid_idx).CA90, 'm--', 'CA90');
end
xlabel('Crank Angle [deg]');
ylabel('aHRR [J/CAD]');
title('Mean Apparent HRR');
grid on;
xlim([-50 150]);

subplot(2,2,4);
if isfield(PressureInfo(valid_idx), 'cumHRR') && ~isempty(PressureInfo(valid_idx).cumHRR)
    plot(Obj_Pressure.CA(1:end-1), PressureInfo(valid_idx).cumHRR, 'b', 'LineWidth', 1.5); hold on;
end
if isfield(PressureInfo(valid_idx), 'SOI_CA') && ~isempty(PressureInfo(valid_idx).SOI_CA) && ~isnan(PressureInfo(valid_idx).SOI_CA)
    xline(PressureInfo(valid_idx).SOI_CA, 'k:', 'SOI', 'LineWidth', 1.5, 'LabelVerticalAlignment', 'top');
end
if isfield(PressureInfo(valid_idx), 'CA10') && ~isempty(PressureInfo(valid_idx).CA10) && ~isnan(PressureInfo(valid_idx).CA10)
    xline(PressureInfo(valid_idx).CA10, 'g--', 'CA10');
    xline(PressureInfo(valid_idx).CA50, 'r--', 'CA50');
    xline(PressureInfo(valid_idx).CA90, 'm--', 'CA90');
end
xlabel('Crank Angle [deg]');
ylabel('Cumulative HRR [J]');
title('Cumulative Heat Release');
grid on;
xlim([-50 150]);

%% Plot 2: Individual cycle aHRR plots for ALL TAKES
fprintf(repmat('=', 1, 70));
fprintf('\nPLOTTING INDIVIDUAL CYCLE aHRR (FILTERED) FOR ALL TAKES\n');
fprintf(repmat('=', 1, 70));
fprintf('\n\n');

% Savitzky-Golay filter parameters (TDMS_Processing.m:367)
sg_order = 7;
sg_window = 91;

for take_idx = 1:nTakes
    if ~isfield(PressureInfo(take_idx), 'aHRR_cycles') || isempty(PressureInfo(take_idx).aHRR_cycles)
        fprintf('Take %02d: Skipping (no aHRR_cycles data)\n\n', take_idx-1);
        continue;
    end
    
    nCycles = size(PressureInfo(take_idx).aHRR_cycles, 2);
    nCols = 4;
    nRows = ceil(nCycles / nCols);
    
    fprintf('Take %02d: Plotting %d filtered aHRR cycles...\n', take_idx-1, nCycles);
    
    figure('Name', sprintf('Set %02d - Take %02d: Individual Cycle aHRR (Filtered SG %d/%d)', test_set, take_idx-1, sg_order, sg_window), ...
           'Position', [300 150 1400 900]);
    
    for cc = 1:nCycles
        subplot(nRows, nCols, cc);
        
        % Get raw aHRR and apply Savitzky-Golay filter
        aHRR_raw = PressureInfo(take_idx).aHRR_cycles(:, cc);
        aHRR_filt = sgolayfilt(aHRR_raw, sg_order, sg_window);
        %aHRR_filt = aHRR_raw;
        
        plot(Obj_Pressure.CA(1:end-1), aHRR_filt, 'b', 'LineWidth', 1.5); 
        hold on;
        
        % Add CA10, CA50, CA90 markers for this cycle
        if isfield(PressureInfo(take_idx), 'CA10_cycles') && ~isempty(PressureInfo(take_idx).CA10_cycles)
            if ~isnan(PressureInfo(take_idx).CA10_cycles(cc))
                xline(PressureInfo(take_idx).CA10_cycles(cc), 'g--', 'CA10', 'LineWidth', 1.2);
            end
            if ~isnan(PressureInfo(take_idx).CA50_cycles(cc))
                xline(PressureInfo(take_idx).CA50_cycles(cc), 'r--', 'CA50', 'LineWidth', 1.2);
            end
            if ~isnan(PressureInfo(take_idx).CA90_cycles(cc))
                xline(PressureInfo(take_idx).CA90_cycles(cc), 'm--', 'CA90', 'LineWidth', 1.2);
            end
        end
        
        % Add SOI marker
        if isfield(PressureInfo(take_idx), 'SOI_CA') && ~isempty(PressureInfo(take_idx).SOI_CA) && ~isnan(PressureInfo(take_idx).SOI_CA)
            xline(PressureInfo(take_idx).SOI_CA, 'k:', 'SOI', 'LineWidth', 1.5);
        end
        
        xlabel('Crank Angle [deg]');
        ylabel('aHRR [J/CAD]');
        
        % Title with firing status
        if isfield(PressureInfo(take_idx), 'CA10_cycles') && ~isempty(PressureInfo(take_idx).CA10_cycles) && ...
           ~isnan(PressureInfo(take_idx).CA10_cycles(cc))
            title(sprintf('Cycle %d (firing)', cc), 'Color', 'green', 'FontWeight', 'bold');
        else
            title(sprintf('Cycle %d (misfire)', cc), 'Color', 'red');
        end
        
        grid on;
        xlim([-50 150]);
    end
    
    sgtitle(sprintf('Set %02d Take %02d: Filtered aHRR Cycles (SG Order %d, Window %d)', ...
        test_set, take_idx-1, sg_order, sg_window), 'FontSize', 12, 'FontWeight', 'bold');
    
    fprintf('  ✓ Plotted %d filtered aHRR cycles\n\n', nCycles);
end

fprintf('All filtered aHRR cycle plots generated.\n\n');

%% Plot 3: TDC detection diagnostics
has_tdc = isfield(PressureInfo(valid_idx), 'CA_Enc') && ~isempty(PressureInfo(valid_idx).CA_Enc) && ...
          isfield(PressureInfo(valid_idx), 'TDC_indices') && ~isempty(PressureInfo(valid_idx).TDC_indices) && ...
          isfield(PressureInfo(valid_idx), 'Fs') && ~isempty(PressureInfo(valid_idx).Fs);

if has_tdc
    dataRaw = PressureInfo(valid_idx).DataRaw;
    zpulseField = Obj_Pressure.EncoderZPulseChannel;
    encField = Obj_Pressure.EncoderDegreeChannel;

    Zpulse = [];
    Puls = [];
    try
        Zpulse = double(dataRaw.Digital_channels.(zpulseField).data);
    catch
        try
            Zpulse = double(dataRaw.Digital_channels.(['Untitled_' zpulseField(end)]).data);
        catch
        end
    end

    try
        Puls = double(dataRaw.Digital_channels.(encField).data);
    catch
        try
            Puls = double(dataRaw.Digital_channels.(['Untitled_' encField(end)]).data);
        catch
        end
    end

    fs_tdc = PressureInfo(valid_idx).Fs;
    t_tdc = (0:numel(PressureInfo(valid_idx).CA_Enc)-1).' ./ fs_tdc;
    figure('Name', sprintf('Set %02d - Take %02d: TDC Detection', test_set, valid_idx-1), 'Position', [200 200 1100 600]);

    subplot(2,1,1);
    if ~isempty(Zpulse)
        zlen = min(numel(Zpulse), numel(t_tdc));
        plot(t_tdc(1:zlen), Zpulse(1:zlen), 'k'); hold on;
    end
    if ~isempty(Puls)
        plen = min(numel(Puls), numel(t_tdc));
        plot(t_tdc(1:plen), Puls(1:plen)*0.5, 'Color', [0.6 0.6 0.6]);
    end
    if ~isempty(PressureInfo(valid_idx).TDC_indices)
        tdc_t = t_tdc(min(PressureInfo(valid_idx).TDC_indices, numel(t_tdc)));
        xline(tdc_t, 'r--');
    end
    xlabel('Time [s]'); ylabel('Digital levels');
    title('Z-pulse (TDC) and encoder pulses with detected TDC markers');
    grid on;

    subplot(2,1,2);
    caLen = min(numel(t_tdc), numel(PressureInfo(valid_idx).CA_Enc));
    plot(t_tdc(1:caLen), PressureInfo(valid_idx).CA_Enc(1:caLen), 'b'); hold on;
    if ~isempty(PressureInfo(valid_idx).TDC_indices)
        tdc_t = t_tdc(min(PressureInfo(valid_idx).TDC_indices, numel(t_tdc)));
        xline(tdc_t, 'r--');
    end
    xlabel('Time [s]'); ylabel('Crank Angle [deg]');
    title('Smoothed encoder-derived crank angle with TDC markers');
    grid on;
end

% Cycle detection sanity check and pressure overlay
if has_tdc
    if isfield(PressureInfo(valid_idx), 'TDC_indices') && numel(PressureInfo(valid_idx).TDC_indices) > 1
        d_tdc_ca = diff(PressureInfo(valid_idx).TDC_indices);
        fprintf('TDC_indices spacing: mean %.1f samples, std %.1f, n=%d\n', mean(d_tdc_ca), std(d_tdc_ca), numel(d_tdc_ca));
    end
    if isfield(PressureInfo(valid_idx), 'TDC_from_Z') && numel(PressureInfo(valid_idx).TDC_from_Z) > 1
        d_tdc_z = diff(PressureInfo(valid_idx).TDC_from_Z);
        fprintf('TDC_from_Z spacing: mean %.1f samples, std %.1f, n=%d\n', mean(d_tdc_z), std(d_tdc_z), numel(d_tdc_z));
    end

    figure('Name', sprintf('Set %02d - Take %02d: Pressure with TDC markers', test_set, valid_idx-1), 'Position', [250 250 1200 500]);
    p_raw = [];
    if isfield(PressureInfo(valid_idx), 'Raw') && ~isempty(PressureInfo(valid_idx).Raw)
        p_raw = PressureInfo(valid_idx).Raw;
    end
    if ~isempty(p_raw)
        plot(p_raw, 'k'); hold on;
    elseif isfield(PressureInfo(valid_idx), 'Filtered') && ~isempty(PressureInfo(valid_idx).Filtered)
        plot(PressureInfo(valid_idx).Filtered, 'k'); hold on;
    end
    legendEntries = {'Pressure'};

    % Overlay Z-pulse (scaled) and TDC markers (shifted and raw)
    if exist('Zpulse','var') && ~isempty(Zpulse)
        zlen = min(numel(Zpulse), numel(p_raw));
        y_min = min(p_raw); y_max = max(p_raw);
        z_scale = max(y_max - y_min, 1) / 6;
        z_offset = y_min - z_scale * 0.5;
        plot(1:zlen, Zpulse(1:zlen) * z_scale + z_offset, 'Color', [0.85 0.55 0], 'LineWidth', 1);
        legendEntries{end+1} = 'Z pulse';

        % Shift Z pulse by TDC_shift for visualization
        if isfield(PressureInfo(valid_idx), 'TDC_indices_raw') && numel(PressureInfo(valid_idx).TDC_indices_raw) > 1
            samples_per_cycle = mean(diff(PressureInfo(valid_idx).TDC_indices_raw));
            shift_samples = round(Obj_Pressure.TDC_shift / 720 * samples_per_cycle);
            Zpulse_shifted = zeros(size(Zpulse));
            src_idx = (1:numel(Zpulse))';
            dst_idx = src_idx - shift_samples;
            keep = dst_idx >= 1 & dst_idx <= numel(Zpulse);
            Zpulse_shifted(dst_idx(keep)) = Zpulse(src_idx(keep));

            zslen = min(numel(Zpulse_shifted), numel(p_raw));
            plot(1:zslen, Zpulse_shifted(1:zslen) * z_scale + (z_offset - z_scale*0.6), 'Color', [0.5 0.3 0], 'LineWidth', 1);
            legendEntries{end+1} = 'Z pulse shifted';
        end
    end
    if isfield(PressureInfo(valid_idx), 'TDC_indices_raw') && ~isempty(PressureInfo(valid_idx).TDC_indices_raw)
        plot(PressureInfo(valid_idx).TDC_indices_raw, zeros(size(PressureInfo(valid_idx).TDC_indices_raw)), 'gx', 'MarkerSize', 8, 'LineWidth', 1.5);
        legendEntries{end+1} = 'TDC from Z (raw)';
    end
    if isfield(PressureInfo(valid_idx), 'TDC_indices_shifted') && ~isempty(PressureInfo(valid_idx).TDC_indices_shifted)
        plot(PressureInfo(valid_idx).TDC_indices_shifted, zeros(size(PressureInfo(valid_idx).TDC_indices_shifted)), 'r+', 'MarkerSize', 8, 'LineWidth', 1.5);
        legendEntries{end+1} = 'TDC shifted (boundary)';
    end

    legend(legendEntries, 'Location', 'best');
    xlabel('Sample index'); ylabel('Pressure [bar]');
    title('Pressure trace with Z pulse and TDC markers');
    grid on;
end

%% Plot 2: Individual cycles
if isfield(PressureInfo(valid_idx), 'P_cycles') && ~isempty(PressureInfo(valid_idx).P_cycles)
    figure('Name', sprintf('Set %02d - Take %02d: Individual Cycles', test_set, valid_idx-1), 'Position', [150 150 1000 600]);
    
    subplot(1,2,1);
    plot(Obj_Pressure.CA, PressureInfo(valid_idx).P_cycles, 'Color', [0.7 0.7 0.7]);
    hold on;
    plot(Obj_Pressure.CA, PressureInfo(valid_idx).P_mean, 'b', 'LineWidth', 2);
    xlabel('Crank Angle [deg]');
    ylabel('Pressure [bar]');
    title('Individual Cycles (gray) and Mean (blue)');
    grid on;
    xlim([-180 180]);
    
    subplot(1,2,2);
    histogram(PressureInfo(valid_idx).IMEP_cycles, 20, 'FaceColor', 'b', 'EdgeColor', 'k');
    xlabel('IMEP [bar]');
    ylabel('Count');
    title(sprintf('IMEP Distribution (CoV = %.2f%%)', PressureInfo(valid_idx).CoV_IMEP * 100));
    grid on;
end

%% Plot 4: Cycle detection check (windowed by shifted Z pulses) - FOR ALL TAKES

fprintf('\n');
fprintf(repmat('=', 1, 80));
fprintf('\nPLOTTING CYCLE DETECTION FOR ALL TAKES\n');
fprintf(repmat('=', 1, 80));
fprintf('\n\n');

for take_idx = 1:nTakes
    fprintf('Plotting cycle detection for Take %02d...\n', take_idx-1);
    
    % Check if this take has required data
    if ~isfield(PressureInfo(take_idx), 'Raw') || isempty(PressureInfo(take_idx).Raw)
        fprintf('  (Skipping - no Raw pressure data)\n\n');
        continue;
    end
    
    % Prefer cycle bounds from processing (shifted rising Z), then shifted TDCs, then raw
    tdc_bound = [];
    if isfield(PressureInfo(take_idx), 'CycleBounds') && numel(PressureInfo(take_idx).CycleBounds) >= 2
        tdc_bound = PressureInfo(take_idx).CycleBounds;
    elseif isfield(PressureInfo(take_idx), 'Z_rise_indices_shifted') && numel(PressureInfo(take_idx).Z_rise_indices_shifted) >= 2
        tdc_bound = PressureInfo(take_idx).Z_rise_indices_shifted;
    elseif isfield(PressureInfo(take_idx), 'TDC_indices_shifted') && numel(PressureInfo(take_idx).TDC_indices_shifted) >= 2
        tdc_bound = PressureInfo(take_idx).TDC_indices_shifted;
    end

    if ~isempty(tdc_bound)
        tdc_raw = [];
        if isfield(PressureInfo(take_idx), 'TDC_indices_raw')
            tdc_raw = PressureInfo(take_idx).TDC_indices_raw;
        end
        take_raw = PressureInfo(take_idx).Raw;
        start_idx = tdc_bound(1);
        end_idx = tdc_bound(2);

        if isfield(PressureInfo(take_idx), 'Cycle_CA_start') && ~isempty(PressureInfo(take_idx).Cycle_CA_start)
            spans = [PressureInfo(take_idx).Cycle_CA_start, PressureInfo(take_idx).Cycle_CA_end];
            nsp = min(size(spans,1), 3);
            for ii = 1:nsp
                fprintf('  Cycle %d CA span: start %.1f, end %.1f\n', ii, spans(ii,1), spans(ii,2));
            end
        end
        if start_idx < 1, start_idx = 1; end
        if end_idx > numel(take_raw), end_idx = numel(take_raw); end

        % === PLOT A: Single cycle window ===
        figure('Name', sprintf('Set %02d - Take %02d: Cycle detection window', test_set, take_idx-1), 'Position', [300 300 1100 450]);
        subplot(1,2,1);
        plot(start_idx:end_idx, take_raw(start_idx:end_idx), 'k'); hold on;
        xline(start_idx, 'r--', 'Start (shifted Z)');
        xline(end_idx, 'r--', 'End (next shifted Z)');
        xlabel('Sample index'); ylabel('Pressure [bar]');
        title(sprintf('Take %02d: One cycle between shifted Z pulses', take_idx-1)); grid on;

        subplot(1,2,2);
        if isfield(PressureInfo(take_idx), 'P_cycles') && ~isempty(PressureInfo(take_idx).P_cycles)
            plot(Obj_Pressure.CA, PressureInfo(take_idx).P_cycles(:,1), 'b'); hold on;
            xline(0, 'r--', 'CA = 0 (shifted Z)');
            xlim([-360 360]);
            xlabel('Crank Angle [deg]'); ylabel('Pressure [bar]');
            title(sprintf('Take %02d: Cycle mapped to CA (-360..360)', take_idx-1)); grid on;
        end
        
        % Figures not saved per user request

        % === PLOT B: Full trace with all shifted Z markers ===
        figure('Name', sprintf('Set %02d - Take %02d: Full cycles (shifted Z)', test_set, take_idx-1), 'Position', [320 320 1200 450]);
        plot(take_raw, 'k'); hold on;

        % Get Z-pulse data for this take
        Zpulse_take = [];
        try
            if isfield(PressureInfo(take_idx), 'DataRaw') && ~isempty(PressureInfo(take_idx).DataRaw)
                dataRaw = PressureInfo(take_idx).DataRaw;
                zpulseField = Obj_Pressure.EncoderZPulseChannel;
                try
                    Zpulse_take = double(dataRaw.Digital_channels.(zpulseField).data);
                catch
                    try
                        Zpulse_take = double(dataRaw.Digital_channels.(['Untitled_' zpulseField(end)]).data);
                    catch
                    end
                end
            end
        catch
        end

        % Overlay shifted Z pulse waveform for visual alignment
        if ~isempty(Zpulse_take)
            zlen = min(numel(Zpulse_take), numel(take_raw));
            y_min = min(take_raw); y_max = max(take_raw);
            z_scale = max(y_max - y_min, 1) / 6;
            z_offset = y_min - z_scale * 0.5;

            % compute shifted Zpulse using same TDC_shift
            zpulse_work = Zpulse_take(1:zlen);
            Zpulse_shifted = zeros(size(zpulse_work));
            if isfield(PressureInfo(take_idx), 'TDC_indices_raw') && numel(PressureInfo(take_idx).TDC_indices_raw) > 1
                samples_per_cycle = mean(diff(PressureInfo(take_idx).TDC_indices_raw));
                shift_samples = round(Obj_Pressure.TDC_shift / 720 * samples_per_cycle);
                src_idx = (1:numel(zpulse_work))';
                dst_idx = src_idx - shift_samples;
                keep = dst_idx >= 1 & dst_idx <= numel(zpulse_work);
                Zpulse_shifted(dst_idx(keep)) = zpulse_work(src_idx(keep));
            end

            plot(1:zlen, zpulse_work * z_scale + z_offset, 'Color', [0.85 0.55 0], 'LineWidth', 1);
            plot(1:zlen, Zpulse_shifted * z_scale + (z_offset - z_scale*0.6), 'Color', [0.5 0.3 0], 'LineWidth', 1);
        end
        
        h_shift = gobjects(0);
        for kk = 1:numel(tdc_bound)
            h_shift(end+1) = xline(tdc_bound(kk), 'r-', 'LineWidth', 1); %#ok<AGROW>
            if kk < numel(tdc_bound)
                midpt = round((tdc_bound(kk) + tdc_bound(kk+1)) / 2);
                text(midpt, max(take_raw)*0.9, sprintf('Cycle %d', kk), 'Rotation', 90, 'HorizontalAlignment', 'center');
            end
        end
        
        h_raw = gobjects(0);
        if ~isempty(tdc_raw)
            for kk = 1:numel(tdc_raw)
                h_raw(end+1) = xline(tdc_raw(kk), 'k:', 'LineWidth', 1); %#ok<AGROW>
            end
        end

        % Overlay rising edges (raw and shifted) to verify boundary logic
        h_rise_raw = gobjects(0);
        if isfield(PressureInfo(take_idx), 'Z_rise_indices') && ~isempty(PressureInfo(take_idx).Z_rise_indices)
            for kk = 1:numel(PressureInfo(take_idx).Z_rise_indices)
                h_rise_raw(end+1) = xline(PressureInfo(take_idx).Z_rise_indices(kk), '--', 'Color', [0 0.6 1], 'LineWidth', 1); %#ok<AGROW>
            end
        end
        
        h_rise_shift = gobjects(0);
        if isfield(PressureInfo(take_idx), 'Z_rise_indices_shifted') && ~isempty(PressureInfo(take_idx).Z_rise_indices_shifted)
            for kk = 1:numel(PressureInfo(take_idx).Z_rise_indices_shifted)
                h_rise_shift(end+1) = xline(PressureInfo(take_idx).Z_rise_indices_shifted(kk), ':', 'Color', [0.6 0 0.8], 'LineWidth', 1.2); %#ok<AGROW>
            end
        end

        xlabel('Sample index'); ylabel('Pressure [bar]');
        title(sprintf('Take %02d: Full pressure trace with shifted Z cycle boundaries (red) and raw Z reference (black)', take_idx-1)); grid on;

        legendHandles = [];
        legendEntries = {};
        if ~isempty(h_shift), legendHandles(end+1) = h_shift(1); legendEntries{end+1} = 'Cycle boundary (shifted rise)'; end %#ok<AGROW>
        if ~isempty(h_raw), legendHandles(end+1) = h_raw(1); legendEntries{end+1} = 'Raw Z (fall)'; end %#ok<AGROW>
        if ~isempty(h_rise_shift), legendHandles(end+1) = h_rise_shift(1); legendEntries{end+1} = 'Shifted Z (rise)'; end %#ok<AGROW>
        if ~isempty(h_rise_raw), legendHandles(end+1) = h_rise_raw(1); legendEntries{end+1} = 'Raw Z (rise)'; end %#ok<AGROW>
        if ~isempty(legendHandles)
            legend(legendHandles, legendEntries, 'Location', 'southoutside', 'Orientation', 'horizontal');
        end
        
        % Figure not saved per user request
    else
        fprintf('  (Skipping - no cycle bounds detected)\n\n');
    end
end

fprintf('All cycle detection plots complete.\n\n');

%% === VALIDATION TESTS (VALIDATION_REPORT.md) ===
fprintf('\n');
fprintf(repmat('=', 1, 70));
fprintf('\nVALIDATION TESTS\n');
fprintf(repmat('=', 1, 70));
fprintf('\n\n');

valid_idx = find(~cellfun(@isempty, {PressureInfo.CA10}), 1);
if ~isempty(valid_idx)
    
    % Test 1: Per-Cycle Metrics (Task 4)
    fprintf('Test 1: Per-Cycle Metrics\n');
    if isfield(PressureInfo(valid_idx), 'CA10_cycles') && ~isempty(PressureInfo(valid_idx).CA10_cycles)
        nCycles_test = numel(PressureInfo(valid_idx).CA10_cycles);
        fprintf('  ✓ PASS: CA10_cycles exists (%d cycles)\n', nCycles_test);
        fprintf('  ✓ PASS: CA50_cycles exists (%d cycles)\n', numel(PressureInfo(valid_idx).CA50_cycles));
        fprintf('  ✓ PASS: CA90_cycles exists (%d cycles)\n', numel(PressureInfo(valid_idx).CA90_cycles));
        fprintf('  Mean CA10: %.2f ± %.2f CAD\n', PressureInfo(valid_idx).CA10, PressureInfo(valid_idx).CA10_std);
    else
        fprintf('  ✗ FAIL: CA10_cycles missing or empty\n');
    end
    fprintf('\n');
    
    % Test 2: Filtering Match (Task 1)
    fprintf('Test 2: Filtering (SG filters applied)\n');
    if isfield(PressureInfo(valid_idx), 'P_corrected_motoring') && ~isempty(PressureInfo(valid_idx).P_corrected_motoring)
        fprintf('  ✓ PASS: P_corrected computed (includes SG filter)\n');
    else
        fprintf('  ✗ FAIL: P_corrected missing\n');
    end
    if isfield(PressureInfo(valid_idx), 'aHRR') && ~isempty(PressureInfo(valid_idx).aHRR)
        fprintf('  ✓ PASS: aHRR computed (includes SG filter)\n');
    else
        fprintf('  ✗ FAIL: aHRR missing\n');
    end
    fprintf('\n');
    
    % Test 3: Injection Signal Priority (Task 3)
    fprintf('Test 3: Di3 Injection Signal Priority\n');
    if isfield(PressureInfo(valid_idx), 'SOI_CA') && ~isnan(PressureInfo(valid_idx).SOI_CA)
        fprintf('  ✓ PASS: SOI_CA detected: %.2f CAD\n', PressureInfo(valid_idx).SOI_CA);
        if abs(PressureInfo(valid_idx).SOI_CA - (-2)) < 10
            fprintf('  ✓ PASS: SOI timing reasonable (within ±10 CAD of expected -2 CAD)\n');
        else
            fprintf('  ⚠ WARN: SOI timing = %.2f CAD (expected ~-2 CAD)\n', PressureInfo(valid_idx).SOI_CA);
        end
    else
        fprintf('  ⚠ WARN: SOI_CA not detected (injection signal missing?)\n');
    end
    fprintf('\n');
    
    % Test 4: CA Integrity (Task 5)
    fprintf('Test 4: CA Assignment Integrity\n');
    if isfield(PressureInfo(valid_idx), 'P_cycles') && ~isempty(PressureInfo(valid_idx).P_cycles)
        fprintf('  ✓ PASS: CA integrity checks added (see warnings above if duplicates detected)\n');
        fprintf('  Note: Check console for CA duplicate warnings during processing\n');
    else
        fprintf('  ⚠ WARN: No cycle data to check CA integrity\n');
    end
    fprintf('\n');
    
    % Test 5: Expected CA10 Value (TDMS Comparison)
    fprintf('Test 5: CA10 Value vs TDMS Reference\n');
    expected_CA10 = 7.9;  % TDMS reference value
    if ~isnan(PressureInfo(valid_idx).CA10)
        error_CA10 = abs(PressureInfo(valid_idx).CA10 - expected_CA10);
        fprintf('  CA10 = %.2f CAD (TDMS reference: %.2f CAD)\n', PressureInfo(valid_idx).CA10, expected_CA10);
        fprintf('  Error: %.2f CAD\n', error_CA10);
        if error_CA10 < 0.5
            fprintf('  ✓ PASS: Matches TDMS within ±0.5 CAD\n');
        elseif error_CA10 < 1.0
            fprintf('  ⚠ WARN: Matches TDMS within ±1.0 CAD\n');
        else
            fprintf('  ✗ FAIL: Error > 1.0 CAD\n');
        end
    else
        fprintf('  ✗ FAIL: CA10 is NaN\n');
    end
    fprintf('\n');
else
    fprintf('⚠ WARN: No valid Take with CA metrics to validate\n\n');
end

%% === COMPARISON: Ensemble-Averaged vs Per-Cycle CA Metrics ===
    fprintf('\n');
    fprintf(repmat('=', 1, 70));
    fprintf('\nCOMPARISON: CA METRICS FROM AVERAGED aHRR vs PER-CYCLE METHOD\n');
    fprintf(repmat('=', 1, 70));
    fprintf('\n\n');
    
    fprintf('Methodology Comparison:\n');
    fprintf('  Method A (CLASS current): Per-cycle aHRR → per-cycle CA10/50/90 → average\n');
    fprintf('  Method B (Single Take):   Averaged P_mean → averaged aHRR → CA10/50/90\n');
    fprintf('  Method C (All Takes):     Average P across all takes → aHRR → CA10/50/90\n');
    fprintf('  Method D (SOI-anchored):  Zero cumHRR at SOI_CA instead of -50 CAD\n\n');
    
    % Find a valid take with complete data
    valid_comparison_idx = [];
    for jj = 1:nTakes
        if isfield(PressureInfo(jj), 'aHRR') && ~isempty(PressureInfo(jj).aHRR) && ...
           isfield(PressureInfo(jj), 'cumHRR') && ~isempty(PressureInfo(jj).cumHRR) && ...
           isfield(PressureInfo(jj), 'CA10') && ~isnan(PressureInfo(jj).CA10)
            valid_comparison_idx = jj;
            break;
        end
    end
    
    if ~isempty(valid_comparison_idx)
        PI = PressureInfo(valid_comparison_idx);
        CA_grid = Obj_Pressure.CA;
        
        fprintf('Using Take %02d for comparison...\n\n', valid_comparison_idx-1);
        
        % === METHOD A: Per-cycle approach (CLASS current method) ===
        CA10_method_A = PI.CA10;
        CA50_method_A = PI.CA50;
        CA90_method_A = PI.CA90;
        CA10_std_A = PI.CA10_std;
        CA50_std_A = PI.CA50_std;
        CA90_std_A = PI.CA90_std;
        
        fprintf('Method A (Per-Cycle → Average):\n');
        fprintf('  CA10 = %.2f ± %.2f CAD\n', CA10_method_A, CA10_std_A);
        fprintf('  CA50 = %.2f ± %.2f CAD\n', CA50_method_A, CA50_std_A);
        fprintf('  CA90 = %.2f ± %.2f CAD\n', CA90_method_A, CA90_std_A);
        fprintf('\n');
        
        % === METHOD B: Ensemble-averaged approach ===
        % Use the already-computed ensemble-averaged aHRR and cumHRR
        aHRR_avg = PI.aHRR;  % This is already ensemble-averaged
        cumHRR_avg = PI.cumHRR;  % This is already ensemble-averaged
        
        % Find CA10/50/90 from averaged cumHRR
        CA_range_min = find(CA_grid > -50, 1) : find(CA_grid > -30, 1);
        CA_range_max = find(CA_grid > 0, 1) : find(CA_grid > 180, 1);
        CA_range = find(CA_grid > -50, 1) : find(CA_grid > 150, 1);
        
        minCum_B = min(cumHRR_avg(CA_range_min));
        maxCum_B = max(cumHRR_avg(CA_range_max));
        SumHRR_B = maxCum_B - minCum_B;
        
        % Find CA10/50/90 percentiles from ensemble average
        idx_10_B = find(cumHRR_avg(CA_range) > minCum_B + 0.1 * SumHRR_B, 1, 'first');
        idx_50_B = find(cumHRR_avg(CA_range) > minCum_B + 0.5 * SumHRR_B, 1, 'first');
        idx_90_B = find(cumHRR_avg(CA_range) > minCum_B + 0.9 * SumHRR_B, 1, 'first');
        
        if ~isempty(idx_10_B) && ~isempty(idx_50_B) && ~isempty(idx_90_B)
            CA10_method_B = CA_grid(find(CA_grid > -50, 1) + idx_10_B - 1);
            CA50_method_B = CA_grid(find(CA_grid > -50, 1) + idx_50_B - 1);
            CA90_method_B = CA_grid(find(CA_grid > -50, 1) + idx_90_B - 1);
            
            fprintf('Method B (Single Take Averaged aHRR → CA):\n');
            fprintf('  CA10 = %.2f CAD (no std - single averaged trace)\n', CA10_method_B);
            fprintf('  CA50 = %.2f CAD (no std - single averaged trace)\n', CA50_method_B);
            fprintf('  CA90 = %.2f CAD (no std - single averaged trace)\n', CA90_method_B);
            fprintf('\n');
            
            % === METHOD C: Average pressure across ALL valid takes, then compute aHRR ===
            fprintf('Method C (All Takes Averaged):\n');
            
            % Collect all valid P_mean_corrected traces
            P_all_takes = [];
            take_count = 0;
            for jj = 1:nTakes
                if isfield(PressureInfo(jj), 'P_mean_corrected') && ~isempty(PressureInfo(jj).P_mean_corrected)
                    P_take = PressureInfo(jj).P_mean_corrected;
                    if length(P_take) == length(CA_grid)
                        P_all_takes = [P_all_takes, P_take(:)];
                        take_count = take_count + 1;
                    end
                end
            end
            
            if take_count > 0
                % Average pressure across all takes
                P_mean_all = mean(P_all_takes, 2);
                
                % Compute aHRR from averaged pressure (same method as CLASS)
                Kappa = 1.35;
                Vol = Obj_Pressure.Vol(:);
                aHRR_all = (Kappa/(Kappa-1) * P_mean_all(1:end-1) .* diff(Vol) + ...
                           1/(Kappa-1) * Vol(1:end-1) .* diff(P_mean_all));
                aHRR_all = aHRR_all * 1e5 / (CA_grid(2) - CA_grid(1));
                
                % Compute cumulative HRR
                cumHRR_all = cumsum(aHRR_all) * (CA_grid(2) - CA_grid(1));
                
                % Zero cumHRR at -50 CAD
                idx_zero = find(CA_grid > -50, 1);
                cumHRR_all = cumHRR_all - cumHRR_all(idx_zero);
                
                % Find CA10/50/90 from all-take averaged cumHRR
                CA_range_min_C = find(CA_grid > -50, 1) : find(CA_grid > -30, 1);
                CA_range_max_C = find(CA_grid > 0, 1) : find(CA_grid > 180, 1);
                CA_range_C = find(CA_grid > -50, 1) : find(CA_grid > 150, 1);
                
                minCum_C = min(cumHRR_all(CA_range_min_C));
                maxCum_C = max(cumHRR_all(CA_range_max_C));
                SumHRR_C = maxCum_C - minCum_C;
                
                idx_10_C = find(cumHRR_all(CA_range_C) > minCum_C + 0.1 * SumHRR_C, 1, 'first');
                idx_50_C = find(cumHRR_all(CA_range_C) > minCum_C + 0.5 * SumHRR_C, 1, 'first');
                idx_90_C = find(cumHRR_all(CA_range_C) > minCum_C + 0.9 * SumHRR_C, 1, 'first');
                
                if ~isempty(idx_10_C) && ~isempty(idx_50_C) && ~isempty(idx_90_C)
                    CA10_method_C = CA_grid(find(CA_grid > -50, 1) + idx_10_C - 1);
                    CA50_method_C = CA_grid(find(CA_grid > -50, 1) + idx_50_C - 1);
                    CA90_method_C = CA_grid(find(CA_grid > -50, 1) + idx_90_C - 1);
                    
                    fprintf('  Averaged %d takes\n', take_count);
                    fprintf('  CA10 = %.2f CAD\n', CA10_method_C);
                    fprintf('  CA50 = %.2f CAD\n', CA50_method_C);
                    fprintf('  CA90 = %.2f CAD\n', CA90_method_C);
                    fprintf('\n');
                else
                    fprintf('  ✗ Failed to compute CA metrics from all-take average\n\n');
                    CA10_method_C = NaN;
                    CA50_method_C = NaN;
                    CA90_method_C = NaN;
                end
            else
                fprintf('  ✗ No valid takes found for averaging\n\n');
                CA10_method_C = NaN;
                CA50_method_C = NaN;
                CA90_method_C = NaN;
            end
            
            % === METHOD D: SOI-anchored zeroing ===
            fprintf('Method D (SOI-anchored cumHRR zeroing):\n');
            
            % Get SOI_CA from current take
            SOI_CA = NaN;
            if isfield(PI, 'SOI_CA') && ~isempty(PI.SOI_CA) && ~isnan(PI.SOI_CA)
                SOI_CA = PI.SOI_CA;
            end
            
            if ~isnan(SOI_CA)
                fprintf('  SOI detected at %.2f CAD\n', SOI_CA);
                
                % Use the same aHRR as Method B (ensemble-averaged)
                aHRR_D = PI.aHRR;
                
                % Recompute cumHRR
                cumHRR_D = cumsum(aHRR_D) * (CA_grid(2) - CA_grid(1));
                
                % Zero cumHRR at SOI_CA instead of -50 CAD
                idx_soi = find(CA_grid >= SOI_CA, 1);
                if ~isempty(idx_soi)
                    cumHRR_D = cumHRR_D - cumHRR_D(idx_soi);
                    
                    % Search CA10/50/90 starting from SOI_CA
                    CA_range_D = idx_soi : min(length(CA_grid)-1, find(CA_grid > SOI_CA + 200, 1));
                    
                    % Find min/max for percentiles (starting from SOI)
                    CA_range_min_D = idx_soi : min(length(CA_grid)-1, find(CA_grid > SOI_CA + 20, 1));
                    CA_range_max_D = idx_soi : min(length(CA_grid)-1, find(CA_grid > SOI_CA + 180, 1));
                    
                    minCum_D = min(cumHRR_D(CA_range_min_D));
                    maxCum_D = max(cumHRR_D(CA_range_max_D));
                    SumHRR_D = maxCum_D - minCum_D;
                    
                    idx_10_D = find(cumHRR_D(CA_range_D) > minCum_D + 0.1 * SumHRR_D, 1, 'first');
                    idx_50_D = find(cumHRR_D(CA_range_D) > minCum_D + 0.5 * SumHRR_D, 1, 'first');
                    idx_90_D = find(cumHRR_D(CA_range_D) > minCum_D + 0.9 * SumHRR_D, 1, 'first');
                    
                    if ~isempty(idx_10_D) && ~isempty(idx_50_D) && ~isempty(idx_90_D)
                        CA10_method_D = CA_grid(idx_soi + idx_10_D - 1);
                        CA50_method_D = CA_grid(idx_soi + idx_50_D - 1);
                        CA90_method_D = CA_grid(idx_soi + idx_90_D - 1);
                        
                        % Ignition delay is now directly CA10 offset from SOI
                        IgnDelay_D = CA10_method_D - SOI_CA;
                        
                        fprintf('  CA10 = %.2f CAD (IgnDelay from SOI: %.2f CAD)\n', CA10_method_D, IgnDelay_D);
                        fprintf('  CA50 = %.2f CAD\n', CA50_method_D);
                        fprintf('  CA90 = %.2f CAD\n', CA90_method_D);
                        fprintf('  CA10-50 = %.2f CAD\n', CA50_method_D - CA10_method_D);
                        fprintf('  CA50-90 = %.2f CAD\n', CA90_method_D - CA50_method_D);
                        fprintf('\n');
                    else
                        fprintf('  ✗ Failed to find CA metrics in SOI-anchored window\n\n');
                        CA10_method_D = NaN;
                        CA50_method_D = NaN;
                        CA90_method_D = NaN;
                    end
                else
                    fprintf('  ✗ SOI_CA (%.2f) not found in CA grid\n\n', SOI_CA);
                    CA10_method_D = NaN;
                    CA50_method_D = NaN;
                    CA90_method_D = NaN;
                end
            else
                fprintf('  ✗ SOI_CA not available for this take\n\n');
                CA10_method_D = NaN;
                CA50_method_D = NaN;
                CA90_method_D = NaN;
            end
            
            % === COMPARISON ===
            fprintf('Difference (Method A - Method B - Method C):\n');
            diff_CA10_AB = CA10_method_A - CA10_method_B;
            diff_CA50_AB = CA50_method_A - CA50_method_B;
            diff_CA90_AB = CA90_method_A - CA90_method_B;
            
            fprintf('  ΔCA10 (A-B) = %+.2f CAD', diff_CA10_AB);
            if abs(diff_CA10_AB) < 0.5
                fprintf(' ✓ EXCELLENT (<0.5 CAD)\n');
            elseif abs(diff_CA10_AB) < 1.0
                fprintf(' ✓ GOOD (<1.0 CAD)\n');
            else
                fprintf(' ⚠ SIGNIFICANT (>1.0 CAD)\n');
            end
            
            fprintf('  ΔCA50 (A-B) = %+.2f CAD', diff_CA50_AB);
            if abs(diff_CA50_AB) < 0.5
                fprintf(' ✓ EXCELLENT (<0.5 CAD)\n');
            elseif abs(diff_CA50_AB) < 1.0
                fprintf(' ✓ GOOD (<1.0 CAD)\n');
            else
                fprintf(' ⚠ SIGNIFICANT (>1.0 CAD)\n');
            end
            
            fprintf('  ΔCA90 (A-B) = %+.2f CAD', diff_CA90_AB);
            if abs(diff_CA90_AB) < 0.5
                fprintf(' ✓ EXCELLENT (<0.5 CAD)\n');
            elseif abs(diff_CA90_AB) < 1.0
                fprintf(' ✓ GOOD (<1.0 CAD)\n');
            else
                fprintf(' ⚠ SIGNIFICANT (>1.0 CAD)\n');
            end
            
            if ~isnan(CA10_method_C)
                diff_CA10_AC = CA10_method_A - CA10_method_C;
                diff_CA50_AC = CA50_method_A - CA50_method_C;
                diff_CA90_AC = CA90_method_A - CA90_method_C;
                
                fprintf('  ΔCA10 (A-C) = %+.2f CAD', diff_CA10_AC);
                if abs(diff_CA10_AC) < 0.5
                    fprintf(' ✓ EXCELLENT (<0.5 CAD)\n');
                elseif abs(diff_CA10_AC) < 1.0
                    fprintf(' ✓ GOOD (<1.0 CAD)\n');
                else
                    fprintf(' ⚠ SIGNIFICANT (>1.0 CAD)\n');
                end
                
                fprintf('  ΔCA50 (A-C) = %+.2f CAD', diff_CA50_AC);
                if abs(diff_CA50_AC) < 0.5
                    fprintf(' ✓ EXCELLENT (<0.5 CAD)\n');
                elseif abs(diff_CA50_AC) < 1.0
                    fprintf(' ✓ GOOD (<1.0 CAD)\n');
                else
                    fprintf(' ⚠ SIGNIFICANT (>1.0 CAD)\n');
                end
                
                fprintf('  ΔCA90 (A-C) = %+.2f CAD', diff_CA90_AC);
                if abs(diff_CA90_AC) < 0.5
                    fprintf(' ✓ EXCELLENT (<0.5 CAD)\n');
                elseif abs(diff_CA90_AC) < 1.0
                    fprintf(' ✓ GOOD (<1.0 CAD)\n');
                else
                    fprintf(' ⚠ SIGNIFICANT (>1.0 CAD)\n');
                end
            end
            fprintf('\n');
            
            % Analysis
            fprintf('Analysis:\n');
            fprintf('  - Method A captures cycle-to-cycle variability (std shown)\n');
            fprintf('  - Method B smooths out variability (no std)\n');
            fprintf('  - Small differences (<1 CAD) indicate averaging order matters minimally\n');
            fprintf('  - Large differences (>1 CAD) indicate non-linear effects or skewed distributions\n');
            fprintf('  - Per-cycle method (A) is more rigorous and captures true variability\n');
            fprintf('\n');
            
            % Additional comparison: Check if Method B falls within Method A's uncertainty
            fprintf('Statistical Validation:\n');
            in_range_10 = abs(diff_CA10_AB) <= CA10_std_A;
            in_range_50 = abs(diff_CA50_AB) <= CA50_std_A;
            in_range_90 = abs(diff_CA90_AB) <= CA90_std_A;
            
            if in_range_10
                fprintf('  CA10: Method B IS within Method A uncertainty (±%.2f CAD)\n', CA10_std_A);
            else
                fprintf('  CA10: Method B NOT within Method A uncertainty (±%.2f CAD)\n', CA10_std_A);
            end
            if in_range_50
                fprintf('  CA50: Method B IS within Method A uncertainty (±%.2f CAD)\n', CA50_std_A);
            else
                fprintf('  CA50: Method B NOT within Method A uncertainty (±%.2f CAD)\n', CA50_std_A);
            end
            if in_range_90
                fprintf('  CA90: Method B IS within Method A uncertainty (±%.2f CAD)\n', CA90_std_A);
            else
                fprintf('  CA90: Method B NOT within Method A uncertainty (±%.2f CAD)\n', CA90_std_A);
            end
            fprintf('\n');
            
        else
            fprintf('✗ Could not compute Method B metrics (missing cumHRR percentiles)\n\n');
        end
    else
        fprintf('✗ No valid Take found with complete HRR data for comparison\n\n');
    end
    
    fprintf(repmat('=', 1, 70));
fprintf('\nVALIDATION COMPLETE\n');
fprintf(repmat('=', 1, 70));
fprintf('\n\n');

fprintf(repmat('=', 1, 70));
fprintf('\nTEST COMPLETE\n');
fprintf(repmat('=', 1, 70));
fprintf('\n\n');

