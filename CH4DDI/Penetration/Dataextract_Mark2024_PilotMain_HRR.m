clear all; clc; close all

% Load the figure
fig = openfig('I:\CH4DDI\Mark‘s paper\JournalFigures\Figure_HRR_DualFuel_Pilot_Main.fig');

% Get all axes
allAxes = findall(fig, 'Type', 'axes');
fprintf('========== TOTAL AXES FOUND: %d ==========\n', length(allAxes));

% Sort axes by position (LEFT to RIGHT instead of top to bottom)
axesPositions = zeros(length(allAxes), 4);
for i = 1:length(allAxes)
    axesPositions(i, :) = get(allAxes(i), 'Position');
end

% Sort by X position (ascending - left first)
[~, sortIdx] = sort(axesPositions(:, 1), 'ascend');
sortedAxes = allAxes(sortIdx);

% Identify subplots (exclude legend axes - they usually have very small positions)
plotAxes = [];
for i = 1:length(sortedAxes)
    pos = get(sortedAxes(i), 'Position');
    % Filter out legend axes (usually very small)
    if pos(3) > 0.1 && pos(4) > 0.1  % width and height > 0.1
        plotAxes(end+1) = sortedAxes(i);
    end
end

fprintf('Main plot axes found: %d\n', length(plotAxes));

if length(plotAxes) >= 2
    leftSubplot = plotAxes(1);     % Left subplot (was top - with shaded areas)
    rightSubplot = plotAxes(2);    % Right subplot (was bottom - clean average lines)
else
    error('Expected at least 2 subplots, found %d', length(plotAxes));
end

fprintf('\n========== ANALYZING LEFT SUBPLOT (with variations) ==========\n');
leftChildren = get(leftSubplot, 'Children');
fprintf('Left subplot children: %d\n', length(leftChildren));

% Separate patch objects (shaded areas) and line objects
leftPatches = findall(leftSubplot, 'Type', 'patch');
leftLines = findall(leftSubplot, 'Type', 'line');

fprintf('Left subplot patches (shaded areas): %d\n', length(leftPatches));
fprintf('Left subplot lines: %d\n', length(leftLines));

fprintf('\n========== ANALYZING RIGHT SUBPLOT (average lines) ==========\n');
rightChildren = get(rightSubplot, 'Children');
fprintf('Right subplot children: %d\n', length(rightChildren));

rightLines = findall(rightSubplot, 'Type', 'line');
fprintf('Right subplot lines: %d\n', length(rightLines));

% Extract average lines from RIGHT subplot (was bottom subplot)
fprintf('\n========== EXTRACTING AVERAGE LINES (Right Subplot) ==========\n');
averageLineData = struct();

for i = 1:length(rightLines)
    currentLine = rightLines(i);
    
    % Get line properties
    xData = get(currentLine, 'XData');
    yData = get(currentLine, 'YData');
    lineStyle = get(currentLine, 'LineStyle');
    lineWidth = get(currentLine, 'LineWidth');
    color = get(currentLine, 'Color');
    displayName = get(currentLine, 'DisplayName');
    
    % Store in structure
    averageLineData(i).xData = xData;
    averageLineData(i).yData = yData;
    averageLineData(i).lineStyle = lineStyle;
    averageLineData(i).lineWidth = lineWidth;
    averageLineData(i).color = color;
    averageLineData(i).displayName = displayName;
    averageLineData(i).dataPoints = length(xData);
    
    fprintf('\n--- Average Line %d ---\n', i);
    fprintf('Display Name: %s\n', displayName);
    fprintf('Data Points: %d\n', length(xData));
    fprintf('X Range: [%.3f, %.3f]\n', min(xData), max(xData));
    fprintf('Y Range: [%.3f, %.3f]\n', min(yData), max(yData));
    fprintf('Color: [%.3f, %.3f, %.3f]\n', color(1), color(2), color(3));
end

% Extract actual colors from the figure and create color mapping
fprintf('\n========== CREATING COLOR MAPPING FROM ORIGINAL FIGURE ==========\n');
originalColorMap = containers.Map();
for i = 1:length(averageLineData)
    if ~isempty(averageLineData(i).displayName)
        conditionName = averageLineData(i).displayName;
        actualColor = averageLineData(i).color;
        originalColorMap(conditionName) = actualColor;
        fprintf('Condition: %s -> Color: [%.3f, %.3f, %.3f]\n', ...
                conditionName, actualColor(1), actualColor(2), actualColor(3));
    end
end

% Extract shaded areas from LEFT subplot (was top subplot)
fprintf('\n========== EXTRACTING SHADED AREAS (Left Subplot) ==========\n');
shadedAreaData = struct();

for i = 1:length(leftPatches)
    currentPatch = leftPatches(i);
    
    % Get patch properties
    xData = get(currentPatch, 'XData');
    yData = get(currentPatch, 'YData');
    faceColor = get(currentPatch, 'FaceColor');
    edgeColor = get(currentPatch, 'EdgeColor');
    faceAlpha = get(currentPatch, 'FaceAlpha');
    displayName = get(currentPatch, 'DisplayName');
    
    % Store in structure
    shadedAreaData(i).xData = xData;
    shadedAreaData(i).yData = yData;
    shadedAreaData(i).faceColor = faceColor;
    shadedAreaData(i).edgeColor = edgeColor;
    shadedAreaData(i).faceAlpha = faceAlpha;
    shadedAreaData(i).displayName = displayName;
    shadedAreaData(i).dataPoints = length(xData);
    
    fprintf('\n--- Shaded Area %d ---\n', i);
    fprintf('Display Name: %s\n', displayName);
    fprintf('Data Points: %d\n', length(xData));
    if ~isempty(xData) && ~isempty(yData)
        fprintf('X Range: [%.3f, %.3f]\n', min(xData), max(xData));
        fprintf('Y Range: [%.3f, %.3f]\n', min(yData), max(yData));
    end
    fprintf('Face Color: [%.3f, %.3f, %.3f]\n', faceColor(1), faceColor(2), faceColor(3));
    fprintf('Face Alpha: %.3f\n', faceAlpha);
end

% Also extract mean lines from LEFT subplot for comparison (was top subplot)
fprintf('\n========== EXTRACTING MEAN LINES (Left Subplot) ==========\n');
leftMeanLineData = struct();

for i = 1:length(leftLines)
    currentLine = leftLines(i);
    
    % Get line properties
    xData = get(currentLine, 'XData');
    yData = get(currentLine, 'YData');
    color = get(currentLine, 'Color');
    displayName = get(currentLine, 'DisplayName');
    
    % Store in structure
    leftMeanLineData(i).xData = xData;
    leftMeanLineData(i).yData = yData;
    leftMeanLineData(i).color = color;
    leftMeanLineData(i).displayName = displayName;
    
    fprintf('\n--- Left Mean Line %d ---\n', i);
    fprintf('Display Name: %s\n', displayName);
    fprintf('Data Points: %d\n', length(xData));
    if ~isempty(xData) && ~isempty(yData)
        fprintf('X Range: [%.3f, %.3f]\n', min(xData), max(xData));
        fprintf('Y Range: [%.3f, %.3f]\n', min(yData), max(yData));
    end
end

fprintf('\n========== DATA EXTRACTION COMPLETE ==========\n');
fprintf('Average lines extracted: %d\n', length(averageLineData));
fprintf('Shaded areas extracted: %d\n', length(shadedAreaData));
fprintf('Left mean lines extracted: %d\n', length(leftMeanLineData));

% =========================================================================
% CA10/CA50/CA90 COMBUSTION ANALYSIS
% =========================================================================

fprintf('\n========== CA10/CA50/CA90 COMBUSTION ANALYSIS ==========\n');

% Check if data exists from previous extraction
if ~exist('averageLineData', 'var') || ~exist('shadedAreaData', 'var')
    error('Please run the data extraction script first to load averageLineData and shadedAreaData');
end

% Initialize results structure
combustionResults = struct();

fprintf('\n========== PROCESSING AVERAGE LINES ==========\n');

% Process each average line (RIGHT subplot data - was bottom subplot)
for lineIdx = 1:length(averageLineData)
    if isempty(averageLineData(lineIdx).xData)
        continue;
    end
    
    condition = averageLineData(lineIdx).displayName;
    time = averageLineData(lineIdx).xData;
    hrr = averageLineData(lineIdx).yData;
    
    fprintf('\n--- Processing: %s ---\n', condition);
    fprintf('Original data: %d points, Time range: [%.3f, %.3f] ms\n', ...
            length(time), min(time), max(time));
    
    % Limit analysis to 0-20 ms time window and remove negative HRR values
    timeWindowIdx = time >= 0 & time <= 20;
    positiveHRRIdx = hrr >= 0;
    validIdx = timeWindowIdx & positiveHRRIdx;
    time_valid = time(validIdx);
    hrr_valid = hrr(validIdx);
    
    fprintf('Filtered data (0-20ms, HRR>0): %d points, Time range: [%.3f, %.3f] ms\n', ...
            length(time_valid), min(time_valid), max(time_valid));
    fprintf('HRR range: [%.3f, %.3f] J/ms\n', min(hrr_valid), max(hrr_valid));
    
    if length(time_valid) < 2
        fprintf('Warning: Insufficient positive HRR data in 0-20ms range for %s\n', condition);
        continue;
    end
    
    % Calculate cumulative heat release using trapezoidal integration (0-20ms, HRR>0 only)
    cumulative_heat = cumtrapz(time_valid, hrr_valid);
    total_heat = cumulative_heat(end);
    
    fprintf('Total heat released (0-20ms): %.2f J\n', total_heat);
    
    % Calculate fractional heat release
    fractional_heat = cumulative_heat / total_heat;
    
    % Find CA10 (10% heat release)
    ca10_idx = find(fractional_heat >= 0.10, 1, 'first');
    if ~isempty(ca10_idx)
        % Linear interpolation for more accurate CA10
        if ca10_idx > 1
            x1 = fractional_heat(ca10_idx-1);
            x2 = fractional_heat(ca10_idx);
            t1 = time_valid(ca10_idx-1);
            t2 = time_valid(ca10_idx);
            ca10_time = t1 + (0.10 - x1) * (t2 - t1) / (x2 - x1);
        else
            ca10_time = time_valid(ca10_idx);
        end
    else
        ca10_time = NaN;
        fprintf('Warning: Could not find CA10 for %s\n', condition);
    end
    
    % Find CA50 (50% heat release)
    ca50_idx = find(fractional_heat >= 0.50, 1, 'first');
    if ~isempty(ca50_idx)
        % Linear interpolation for more accurate CA50
        if ca50_idx > 1
            x1 = fractional_heat(ca50_idx-1);
            x2 = fractional_heat(ca50_idx);
            t1 = time_valid(ca50_idx-1);
            t2 = time_valid(ca50_idx);
            ca50_time = t1 + (0.50 - x1) * (t2 - t1) / (x2 - x1);
        else
            ca50_time = time_valid(ca50_idx);
        end
    else
        ca50_time = NaN;
        fprintf('Warning: Could not find CA50 for %s\n', condition);
    end
    
    % Find CA90 (90% heat release)
    ca90_idx = find(fractional_heat >= 0.90, 1, 'first');
    if ~isempty(ca90_idx)
        % Linear interpolation for more accurate CA90
        if ca90_idx > 1
            x1 = fractional_heat(ca90_idx-1);
            x2 = fractional_heat(ca90_idx);
            t1 = time_valid(ca90_idx-1);
            t2 = time_valid(ca90_idx);
            ca90_time = t1 + (0.90 - x1) * (t2 - t1) / (x2 - x1);
        else
            ca90_time = time_valid(ca90_idx);
        end
    else
        ca90_time = NaN;
        fprintf('Warning: Could not find CA90 for %s\n', condition);
    end
    
    % Use CA90 as burn duration
    burn_duration = ca90_time;
    
    % Store results
    combustionResults(lineIdx).condition = condition;
    combustionResults(lineIdx).total_heat = total_heat;
    combustionResults(lineIdx).ca10_time = ca10_time;
    combustionResults(lineIdx).ca50_time = ca50_time;
    combustionResults(lineIdx).ca90_time = ca90_time;
    combustionResults(lineIdx).burn_duration = burn_duration;
    combustionResults(lineIdx).time_valid = time_valid;
    combustionResults(lineIdx).hrr_valid = hrr_valid;
    combustionResults(lineIdx).cumulative_heat = cumulative_heat;
    combustionResults(lineIdx).fractional_heat = fractional_heat;
    
    % Display results
    fprintf('--- RESULTS (based on 0-20ms integration, HRR>0 only) ---\n');
    fprintf('Total Heat Released: %.2f J\n', total_heat);
    fprintf('CA10 Time: %.3f ms\n', ca10_time);
    fprintf('CA50 Time: %.3f ms\n', ca50_time);
    fprintf('CA90 Time: %.3f ms\n', ca90_time);
    fprintf('Burn Duration (CA90): %.3f ms\n', burn_duration);
end

fprintf('\n========== ANALYZING VARIATIONS FROM SHADED AREAS ==========\n');

% For variation analysis, we need to extract upper and lower bounds from patches
% and calculate CA10/CA50/CA90 for each bound to get standard deviations

variationResults = struct();

% Group patches by color/condition (assuming patches come in pairs - upper and lower bounds)
fprintf('Total patches found: %d\n', length(shadedAreaData));

% We'll try to match patches to conditions by color similarity
for lineIdx = 1:length(averageLineData)
    if isempty(averageLineData(lineIdx).xData)
        continue;
    end
    
    condition = averageLineData(lineIdx).displayName;
    lineColor = averageLineData(lineIdx).color;
    
    fprintf('\n--- Analyzing variations for: %s ---\n', condition);
    fprintf('Line color: [%.3f, %.3f, %.3f]\n', lineColor(1), lineColor(2), lineColor(3));
    
    % Find matching patches by color similarity (tolerance for color matching)
    colorTolerance = 0.3;
    matchingPatches = [];
    
    for patchIdx = 1:length(shadedAreaData)
        if isempty(shadedAreaData(patchIdx).xData)
            continue;
        end
        
        patchColor = shadedAreaData(patchIdx).faceColor;
        colorDiff = sqrt(sum((lineColor - patchColor).^2));
        
        if colorDiff < colorTolerance
            matchingPatches(end+1) = patchIdx;
            fprintf('  Matched patch %d (color diff: %.3f)\n', patchIdx, colorDiff);
        end
    end
    
    if isempty(matchingPatches)
        fprintf('  No matching patches found for %s\n', condition);
        continue;
    end
    
    % Extract boundary data from all matching patches
    allBoundaryData = [];
    for patchIdx = matchingPatches
        patchX = shadedAreaData(patchIdx).xData;
        patchY = shadedAreaData(patchIdx).yData;
        
        % Remove NaN values
        validIdx = ~isnan(patchX) & ~isnan(patchY);
        patchX = patchX(validIdx);
        patchY = patchY(validIdx);
        
        if ~isempty(patchX)
            allBoundaryData = [allBoundaryData; patchX(:), patchY(:)];
        end
    end
    
    if isempty(allBoundaryData)
        fprintf('  No valid boundary data for %s\n', condition);
        continue;
    end
    
    % Sort boundary data by time
    [~, sortIdx] = sort(allBoundaryData(:, 1));
    allBoundaryData = allBoundaryData(sortIdx, :);
    
    % Extract upper and lower bounds
    % Strategy: for each unique time point, find min and max HRR values
    uniqueTimes = unique(allBoundaryData(:, 1));
    upperBound = zeros(size(uniqueTimes));
    lowerBound = zeros(size(uniqueTimes));
    
    for i = 1:length(uniqueTimes)
        timeIdx = allBoundaryData(:, 1) == uniqueTimes(i);
        hrrValues = allBoundaryData(timeIdx, 2);
        upperBound(i) = max(hrrValues);
        lowerBound(i) = min(hrrValues);
    end
    
    % Calculate CA10/CA50/CA90 for upper and lower bounds
    bounds = {upperBound, lowerBound};
    boundNames = {'Upper', 'Lower'};
    boundResults = struct();
    
    for boundIdx = 1:2
        currentBound = bounds{boundIdx};
        boundName = boundNames{boundIdx};
        
        % Limit to 0-20 ms time window and remove negative values (same as main analysis)
        timeWindowIdx = uniqueTimes >= 0 & uniqueTimes <= 20;
        positiveHRRIdx = currentBound >= 0;
        validIdx = timeWindowIdx & positiveHRRIdx;
        time_bound = uniqueTimes(validIdx);
        hrr_bound = currentBound(validIdx);
        
        if length(time_bound) < 2
            continue;
        end
        
        % Calculate cumulative heat release
        cumulative_heat_bound = cumtrapz(time_bound, hrr_bound);
        total_heat_bound = cumulative_heat_bound(end);
        fractional_heat_bound = cumulative_heat_bound / total_heat_bound;
        
        % Find CA10
        ca10_idx = find(fractional_heat_bound >= 0.10, 1, 'first');
        if ~isempty(ca10_idx) && ca10_idx > 1
            x1 = fractional_heat_bound(ca10_idx-1);
            x2 = fractional_heat_bound(ca10_idx);
            t1 = time_bound(ca10_idx-1);
            t2 = time_bound(ca10_idx);
            ca10_bound = t1 + (0.10 - x1) * (t2 - t1) / (x2 - x1);
        elseif ~isempty(ca10_idx)
            ca10_bound = time_bound(ca10_idx);
        else
            ca10_bound = NaN;
        end
        
        % Find CA50
        ca50_idx = find(fractional_heat_bound >= 0.50, 1, 'first');
        if ~isempty(ca50_idx) && ca50_idx > 1
            x1 = fractional_heat_bound(ca50_idx-1);
            x2 = fractional_heat_bound(ca50_idx);
            t1 = time_bound(ca50_idx-1);
            t2 = time_bound(ca50_idx);
            ca50_bound = t1 + (0.50 - x1) * (t2 - t1) / (x2 - x1);
        elseif ~isempty(ca50_idx)
            ca50_bound = time_bound(ca50_idx);
        else
            ca50_bound = NaN;
        end
        
        % Find CA90
        ca90_idx = find(fractional_heat_bound >= 0.90, 1, 'first');
        if ~isempty(ca90_idx) && ca90_idx > 1
            x1 = fractional_heat_bound(ca90_idx-1);
            x2 = fractional_heat_bound(ca90_idx);
            t1 = time_bound(ca90_idx-1);
            t2 = time_bound(ca90_idx);
            ca90_bound = t1 + (0.90 - x1) * (t2 - t1) / (x2 - x1);
        elseif ~isempty(ca90_idx)
            ca90_bound = time_bound(ca90_idx);
        else
            ca90_bound = NaN;
        end
        
        % Use CA90 as burn duration
        burn_duration_bound = ca90_bound;
        
        boundResults(boundIdx).boundName = boundName;
        boundResults(boundIdx).ca10 = ca10_bound;
        boundResults(boundIdx).ca50 = ca50_bound;
        boundResults(boundIdx).ca90 = ca90_bound;
        boundResults(boundIdx).burn_duration = burn_duration_bound;
        
        fprintf('  %s bound - CA10: %.3f ms, CA50: %.3f ms, CA90: %.3f ms, Duration: %.3f ms\n', ...
                boundName, ca10_bound, ca50_bound, ca90_bound, burn_duration_bound);
    end
    
    % Calculate standard deviations
    if length(boundResults) >= 2
        ca10_values = [boundResults.ca10];
        ca50_values = [boundResults.ca50];
        ca90_values = [boundResults.ca90];
        duration_values = [boundResults.burn_duration];
        
        % Remove NaN values before calculating std
        ca10_valid = ca10_values(~isnan(ca10_values));
        ca50_valid = ca50_values(~isnan(ca50_values));
        ca90_valid = ca90_values(~isnan(ca90_values));
        duration_valid = duration_values(~isnan(duration_values));
        
        if length(ca10_valid) >= 2
            ca10_std = std(ca10_valid);
        else
            ca10_std = NaN;
        end
        
        if length(ca50_valid) >= 2
            ca50_std = std(ca50_valid);
        else
            ca50_std = NaN;
        end
        
        if length(ca90_valid) >= 2
            ca90_std = std(ca90_valid);
        else
            ca90_std = NaN;
        end
        
        if length(duration_valid) >= 2
            duration_std = std(duration_valid);
        else
            duration_std = NaN;
        end
        
        variationResults(lineIdx).condition = condition;
        variationResults(lineIdx).ca10_std = ca10_std;
        variationResults(lineIdx).ca50_std = ca50_std;
        variationResults(lineIdx).ca90_std = ca90_std;
        variationResults(lineIdx).burn_duration_std = duration_std;
        variationResults(lineIdx).boundResults = boundResults;
        
        fprintf('  Standard Deviations:\n');
        fprintf('    CA10 std: %.3f ms\n', ca10_std);
        fprintf('    CA50 std: %.3f ms\n', ca50_std);
        fprintf('    CA90 std: %.3f ms\n', ca90_std);
        fprintf('    Burn Duration (CA90) std: %.3f ms\n', duration_std);
    end
end

fprintf('\n========== SUMMARY TABLE ==========\n');
fprintf('%-15s | %8s | %8s | %8s | %8s | %8s | %8s | %8s | %8s\n', ...
        'Condition', 'CA10', 'CA10_std', 'CA50', 'CA50_std', 'CA90', 'CA90_std', 'Duration', 'Dur_std');
fprintf('%-15s-|-%8s-|-%8s-|-%8s-|-%8s-|-%8s-|-%8s-|-%8s-|-%8s\n', ...
        repmat('-', 1, 15), repmat('-', 1, 8), repmat('-', 1, 8), repmat('-', 1, 8), ...
        repmat('-', 1, 8), repmat('-', 1, 8), repmat('-', 1, 8), repmat('-', 1, 8), repmat('-', 1, 8));

for lineIdx = 1:length(combustionResults)
    if isempty(combustionResults(lineIdx).condition)
        continue;
    end
    
    condition = combustionResults(lineIdx).condition;
    ca10_mean = combustionResults(lineIdx).ca10_time;
    ca50_mean = combustionResults(lineIdx).ca50_time;
    ca90_mean = combustionResults(lineIdx).ca90_time;
    duration_mean = combustionResults(lineIdx).burn_duration;
    
    % Find corresponding variation data
    ca10_std = NaN;
    ca50_std = NaN;
    ca90_std = NaN;
    duration_std = NaN;
    
    for varIdx = 1:length(variationResults)
        if strcmp(variationResults(varIdx).condition, condition)
            ca10_std = variationResults(varIdx).ca10_std;
            ca50_std = variationResults(varIdx).ca50_std;
            ca90_std = variationResults(varIdx).ca90_std;
            duration_std = variationResults(varIdx).burn_duration_std;
            break;
        end
    end
    
    fprintf('%-15s | %8.3f | %8.3f | %8.3f | %8.3f | %8.3f | %8.3f | %8.3f | %8.3f\n', ...
            condition, ca10_mean, ca10_std, ca50_mean, ca50_std, ca90_mean, ca90_std, duration_mean, duration_std);
end

% Create results table
conditions = {};
ca10_means = [];
ca10_stds = [];
ca50_means = [];
ca50_stds = [];
ca90_means = [];
ca90_stds = [];
duration_means = [];
duration_stds = [];

for lineIdx = 1:length(combustionResults)
    if isempty(combustionResults(lineIdx).condition)
        continue;
    end
    
    conditions{end+1} = combustionResults(lineIdx).condition;
    ca10_means(end+1) = combustionResults(lineIdx).ca10_time;
    ca50_means(end+1) = combustionResults(lineIdx).ca50_time;
    ca90_means(end+1) = combustionResults(lineIdx).ca90_time;
    duration_means(end+1) = combustionResults(lineIdx).burn_duration;
    
    % Find corresponding variation data
    ca10_std = NaN;
    ca50_std = NaN;
    ca90_std = NaN;
    duration_std = NaN;
    
    for varIdx = 1:length(variationResults)
        if strcmp(variationResults(varIdx).condition, combustionResults(lineIdx).condition)
            ca10_std = variationResults(varIdx).ca10_std;
            ca50_std = variationResults(varIdx).ca50_std;
            ca90_std = variationResults(varIdx).ca90_std;
            duration_std = variationResults(varIdx).burn_duration_std;
            break;
        end
    end
    
    ca10_stds(end+1) = ca10_std;
    ca50_stds(end+1) = ca50_std;
    ca90_stds(end+1) = ca90_std;
    duration_stds(end+1) = duration_std;
end

% Create figure with white background
figure('Name', 'Combustion Analysis Results', 'Position', [100, 100, 1200, 800], 'Color', 'white');

% Plot 1: CA10, CA50, CA90 comparison with CA metrics on x-axis
subplot(2, 2, 1);

% Prepare data for plotting - transpose to have CA metrics on x-axis
ca_metrics = {'CA10', 'CA50', 'CA90'};
x_ca = 1:3;  % x-axis positions for CA10, CA50, CA90

% Plot lines for each condition
hold on;
for condIdx = 1:length(conditions)
    % Get CA values for this condition
    ca_values = [ca10_means(condIdx), ca50_means(condIdx), ca90_means(condIdx)];
    ca_stds = [ca10_stds(condIdx), ca50_stds(condIdx), ca90_stds(condIdx)];
    
    % Use actual color from original figure
    lineColor = getConditionColorFromMap(conditions{condIdx}, originalColorMap);
    
    % Plot line (all solid lines, no markers)
    plot(x_ca, ca_values, '-', 'Color', lineColor, 'LineWidth', 1, ...
         'DisplayName', conditions{condIdx});
    
    % Add error bars
    errorbar(x_ca, ca_values, ca_stds, 'Color', lineColor, 'LineStyle', 'none', 'LineWidth', 1);
end

set(gca, 'XTick', x_ca, 'XTickLabel', ca_metrics);
xlim([0.8 3.2])
ylabel('Time [ms]');
title('Combuistion Phasing Comparison');
%legend('show', 'Location', 'best');
grid on;
set(gca, 'Color', 'white');
% Add black boundary
set(gca, 'XColor', 'black', 'YColor', 'black', 'LineWidth', 1);
box on;

% Plot 2: Burn Duration (CA90)
subplot(2, 2, 2);
x_pos = 1:length(conditions);

% Plot each condition with its color
hold on;
for condIdx = 1:length(conditions)
    % Use actual color from original figure
    lineColor = getConditionColorFromMap(conditions{condIdx}, originalColorMap);
    
    plot(x_pos(condIdx), duration_means(condIdx), 'o', 'Color', lineColor, ...
         'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', lineColor, ...
         'DisplayName', conditions{condIdx});
    
    % Add error bar
    errorbar(x_pos(condIdx), duration_means(condIdx), duration_stds(condIdx), ...
             'Color', lineColor, 'LineStyle', 'none', 'LineWidth', 1);
end

set(gca, 'XTick', x_pos, 'XTickLabel', conditions, 'XTickLabelRotation', 0);
ylabel('Duration [ms]');
title('Burn Duration (CA90)');
%legend('show', 'Location', 'best');
grid on;
set(gca, 'Color', 'white');
% Add black boundary
set(gca, 'XColor', 'black', 'YColor', 'black', 'LineWidth', 1);
box on;

% Plot 3: HRR curves with CA markers (using original figure color scheme)
subplot(2, 1, 2);
hold on;

% Sort conditions in the order shown in legend
conditionOrder = {};
conditionColors = {};
orderedResults = {};

for condIdx = 1:length(combustionResults)
    if ~isempty(combustionResults(condIdx).condition)
        conditionOrder{end+1} = combustionResults(condIdx).condition;
        orderedResults{end+1} = combustionResults(condIdx);
        % Use actual color from original figure
        conditionColors{end+1} = getConditionColorFromMap(combustionResults(condIdx).condition, originalColorMap);
    end
end

% Flip the order
conditionOrder = flip(conditionOrder);
orderedResults = flip(orderedResults);
conditionColors = flip(conditionColors);

% Plot in order
for plotIdx = 1:length(orderedResults)
    result = orderedResults{plotIdx};
    lineColor = conditionColors{plotIdx};
    
    time_data = result.time_valid;
    hrr_data = result.hrr_valid;
    ca10_time = result.ca10_time;
    ca50_time = result.ca50_time;
    ca90_time = result.ca90_time;
    condition = result.condition;
    
    % Get original data from averageLineData (before any filtering)
    for origIdx = 1:length(averageLineData)
        if strcmp(averageLineData(origIdx).displayName, condition)
            original_time = averageLineData(origIdx).xData;
            original_hrr = averageLineData(origIdx).yData;
            break;
        end
    end
    plot(original_time, original_hrr, 'Color', lineColor, 'LineWidth', 1.5, 'DisplayName', condition);
    
    % Mark CA10, CA50, and CA90 with consistent markers
    if ~isnan(ca10_time)
        ca10_hrr = interp1(time_data, hrr_data, ca10_time);
        plot(ca10_time, ca10_hrr, 'o', 'Color', lineColor, 'MarkerSize', 8, ...
             'MarkerFaceColor', lineColor, 'MarkerEdgeColor', 'r', 'LineWidth', 1);
    end
    
    if ~isnan(ca50_time)
        ca50_hrr = interp1(time_data, hrr_data, ca50_time);
        plot(ca50_time, ca50_hrr, 's', 'Color', lineColor, 'MarkerSize', 8, ...
             'MarkerFaceColor', lineColor, 'MarkerEdgeColor', 'r', 'LineWidth', 1);
    end
    
    if ~isnan(ca90_time)
        ca90_hrr = interp1(time_data, hrr_data, ca90_time);
        plot(ca90_time, ca90_hrr, '^', 'Color', lineColor, 'MarkerSize', 8, ...
             'MarkerFaceColor', lineColor, 'MarkerEdgeColor', 'r', 'LineWidth', 1);
    end
end

xlabel('Time relative to main fuel SOI [ms]');
ylabel('HRR [J/ms]');
title('HRR Curves with Combustion Phasing Markers');

% Add marker explanation
text(0.02, 0.95, 'Markers: ○ CA10, □ CA50, △ CA90', 'Units', 'normalized', ...
     'FontSize', 10, 'BackgroundColor', 'white');

% Create legend with proper order and symbols matching original
% Create legend lines without markers
legendHandles = [];
for plotIdx = 1:length(orderedResults)
    condition = orderedResults{plotIdx}.condition;
    lineColor = conditionColors{plotIdx};
    
    % Create invisible line just for legend (solid line, no markers)
    h = plot(NaN, NaN, '-', 'Color', lineColor, 'LineWidth', 2, 'DisplayName', condition);
    legendHandles(end+1) = h;
end

legend(legendHandles, conditionOrder, 'Location', 'best');
ylim([-50 600])
xlim([-0.5 20])
grid on;
set(gca, 'Color', 'white');
% Add black boundary
set(gca, 'XColor', 'black', 'YColor', 'black', 'LineWidth', 1);
box on;

% Set white background for entire figure
set(gcf, 'Color', 'white');

fprintf('Visualization created with white background and original color scheme\n');

fprintf('\n========== ANALYSIS COMPLETE ==========\n');
fprintf('Results summary:\n');
fprintf('- CA10, CA50, CA90, and burn duration (CA90) calculated for each condition\n');
fprintf('- Standard deviations calculated from shaded area variations\n');
fprintf('- File saving disabled per user request\n');
fprintf('- Visualization plots created with white background and original color scheme\n');
fprintf('- Line plots with error bars instead of bar charts\n');
fprintf('- CA90 used as burn duration metric\n');
fprintf('- Colors extracted from original loaded figure\n');

% Helper function to get color from original figure (placed at the back)
function color = getConditionColorFromMap(conditionName, colorMap)
    if isKey(colorMap, conditionName)
        color = colorMap(conditionName);
    else
        % Default fallback color if condition not found
        color = [0, 0, 0];  % black
        fprintf('Warning: Color not found for condition %s, using black\n', conditionName);
    end
end
