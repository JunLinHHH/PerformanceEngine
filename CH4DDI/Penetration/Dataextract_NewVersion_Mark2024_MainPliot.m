clear all; clc; close all
% Load the figure
fig = openfig('C:\Users\jun-y\OneDrive - UNSW\CH4DDI\Mark2024 CH4DDI CVCC\JournalFigures\Figure_PenetrationRecession_Main_Pilot.fig');

% Get all axes (excluding legend axes)
allAxes = findall(fig, 'Type', 'axes');

fprintf('========== TOTAL AXES FOUND: %d ==========\n', length(allAxes));

% Filter out legend axes
nonLegendAxes = [];
for i = 1:length(allAxes)
    if ~strcmp(get(allAxes(i), 'Tag'), 'legend')
        nonLegendAxes = [nonLegendAxes; allAxes(i)];
    end
end

fprintf('Non-legend axes: %d\n', length(nonLegendAxes));

% Get positions and sort by Y (descending) then X (ascending)
positions = zeros(length(nonLegendAxes), 4);
for i = 1:length(nonLegendAxes)
    positions(i, :) = nonLegendAxes(i).Position;
end

[~, sortIdx] = sortrows(positions, [-2, 1]);
allAxes = nonLegendAxes(sortIdx);

fprintf('\n========== IDENTIFYING D-CASE AXES ==========\n');
fprintf('Based on the analysis, the 4 D-case subplots are at the BOTTOM ROW\n');

% Find axes with Y position around 50 (bottom row) and have data
dCaseAxes = [];
for i = 1:length(allAxes)
    yPos = allAxes(i).Position(2);
    lines = findall(allAxes(i), 'Type', 'line');
    numLinesWithData = 0;
    
    % Count lines with more than 2 points
    for j = 1:length(lines)
        if length(lines(j).XData) > 10
            numLinesWithData = numLinesWithData + 1;
        end
    end
    
    % Bottom row axes with actual data
    if yPos < 100 && numLinesWithData >= 3
        dCaseAxes = [dCaseAxes, i];
        fprintf('Found D-case at axes %d: Y=%.1f, Lines with data=%d\n', ...
            i, yPos, numLinesWithData);
    end
end

fprintf('\nD-case axes indices: %s\n', mat2str(dCaseAxes));

% Define subplot names and get actual subplot titles from the figure
subplotNames = {'D_0_03ms_M', 'D_1_03ms_M', 'D_2_03ms_M', 'D_3_03ms_M'};
subplotTitles = {'D-0.03ms-M', 'D-1.03ms-M', 'D-2.03ms-M', 'D-3.03ms-M'};

actualSubplotNames = cell(4, 1);
for i = 1:min(4, length(dCaseAxes))
    axIdx = dCaseAxes(i);
    currentAx = allAxes(axIdx);
    titleObj = get(currentAx, 'Title');
    actualTitle = get(titleObj, 'String');
    if ~isempty(actualTitle)
        % Clean up the title for use as field name
        cleanTitle = strrep(actualTitle, '-', '_');
        cleanTitle = strrep(cleanTitle, '.', '_');
        cleanTitle = strrep(cleanTitle, ' ', '_');
        cleanTitle = matlab.lang.makeValidName(cleanTitle);
        actualSubplotNames{i} = cleanTitle;
        fprintf('Subplot %d actual title: %s -> %s\n', i, actualTitle, cleanTitle);
    else
        actualSubplotNames{i} = subplotNames{i}; % fallback
    end
end

% Use actual names if available
if ~any(cellfun(@isempty, actualSubplotNames))
    subplotNames = actualSubplotNames;
end

% Define save folder
saveFolder = 'C:\Users\jun-y\OneDrive - UNSW\CH4DDI\CH4DDI vs H2DDI CVCC\Mark2024_handled';

if ~exist(saveFolder, 'dir')
    mkdir(saveFolder);
end

% Define color-to-name mapping
colorMap = {
    struct('rgb', [0.043, 0.702, 0.255], 'name', 'NonReactingJetPenetration', 'desc', 'Non-reacting jet penetration');
    struct('rgb', [0.596, 0.000, 0.000], 'name', 'Recession_average', 'desc', 'Recession (average)');
    struct('rgb', [0.094, 0.271, 0.800], 'name', 'Penetration_average', 'desc', 'Penetration (average)');
};

% Initialize data storage with proper field names
allData = struct();
for i = 1:4
    fieldName = matlab.lang.makeValidName(subplotNames{i});
    allData.(fieldName) = struct();
end

fprintf('\n========== PROCESSING SUBPLOTS ==========\n');

% Process each D-case (ENHANCED TO INCLUDE MARKERS AND ALL LINES)
for dataIdx = 1:min(4, length(dCaseAxes))
    axIdx = dCaseAxes(dataIdx);
    currentAx = allAxes(axIdx);
    
    % Get the current subplot field name
    currentSubplotName = matlab.lang.makeValidName(subplotNames{dataIdx});
    
    fprintf('\n--- Processing D-case %d: Axes index %d (%s) ---\n', ...
        dataIdx, axIdx, currentSubplotName);
    
    % Get all lines and patches
    lines = findall(currentAx, 'Type', 'line');
    patches = findall(currentAx, 'Type', 'patch');
    
    fprintf('  Total lines: %d, Total patches: %d\n', length(lines), length(patches));
    
    % Process lines (ENHANCED VERSION)
    lineCount = 0;
    markerCount = 0;
    
    for i = 1:length(lines)
        numPoints = length(lines(i).XData);
        lineStyle = lines(i).LineStyle;
        lineColor = lines(i).Color;
        marker = lines(i).Marker;
        markerSize = lines(i).MarkerSize;
        
        % FIRST: Process single-point markers (the key addition!)
        if numPoints == 1 && ~strcmp(marker, 'none') && ~isnan(lines(i).XData) && ~isnan(lines(i).YData)
            markerCount = markerCount + 1;
            
            % Determine marker type based on the actual marker symbol
            if strcmp(marker, '*')
                fieldName = 'IgnitionDelay_Main';
                displayName = 'Ignition delay (main)';
            elseif strcmp(marker, 'o')
                fieldName = 'JetInteraction';
                displayName = 'Jet interaction';
            else
                fieldName = sprintf('Marker_%s', marker);
                displayName = sprintf('Marker %s', marker);
            end
            
            % Handle duplicates for markers
            originalFieldName = fieldName;
            if isfield(allData.(currentSubplotName), fieldName)
                count = 2;
                while isfield(allData.(currentSubplotName), sprintf('%s_%d', originalFieldName, count))
                    count = count + 1;
                end
                fieldName = sprintf('%s_%d', originalFieldName, count);
            end
            
            % Store marker data
            allData.(currentSubplotName).(fieldName).XData = lines(i).XData;
            allData.(currentSubplotName).(fieldName).YData = lines(i).YData;
            allData.(currentSubplotName).(fieldName).Color = lines(i).Color;
            allData.(currentSubplotName).(fieldName).LineStyle = lines(i).LineStyle;
            allData.(currentSubplotName).(fieldName).LineWidth = lines(i).LineWidth;
            allData.(currentSubplotName).(fieldName).Marker = lines(i).Marker;
            allData.(currentSubplotName).(fieldName).MarkerSize = lines(i).MarkerSize;
            
            % Handle marker colors safely
            try
                allData.(currentSubplotName).(fieldName).MarkerFaceColor = lines(i).MarkerFaceColor;
            catch
                allData.(currentSubplotName).(fieldName).MarkerFaceColor = lineColor;
            end
            
            try
                allData.(currentSubplotName).(fieldName).MarkerEdgeColor = lines(i).MarkerEdgeColor;
            catch
                allData.(currentSubplotName).(fieldName).MarkerEdgeColor = lineColor;
            end
            
            allData.(currentSubplotName).(fieldName).DisplayName = displayName;
            allData.(currentSubplotName).(fieldName).Type = 'marker';
            
            fprintf('  [M%d] Added Marker: %s at [%.2f, %.2f]\n', ...
                markerCount, fieldName, lines(i).XData, lines(i).YData);
            
            continue; % Skip to next line
        end
        
        % SECOND: Process continuous lines (both solid AND dashed)
        % Skip short lines (legend entries)
        if numPoints <= 2
            continue;
        end
        
        % Skip lines with all NaN data
        if all(isnan(lines(i).YData))
            continue;
        end
        
        % Process both solid and dashed lines, but classify correctly
        if strcmp(lineStyle, '-') || strcmp(lineStyle, '--')
            lineCount = lineCount + 1;
            
            % Match color and determine if it's average or individual data
            [baseFieldName, dist, baseDisplayName] = findClosestColor(lineColor, colorMap);
            
            % For dashed lines, append suffix to indicate they're not averages
            if strcmp(lineStyle, '--')
                if strcmp(baseFieldName, 'Penetration_average')
                    fieldName = 'Penetration_individual';
                    displayName = 'Penetration (individual)';
                elseif strcmp(baseFieldName, 'Recession_average')
                    fieldName = 'Recession_individual';
                    displayName = 'Recession (individual)';
                else
                    fieldName = sprintf('%s_dashed', baseFieldName);
                    displayName = sprintf('%s (dashed)', baseDisplayName);
                end
            else
                % Solid lines keep their original names
                fieldName = baseFieldName;
                displayName = baseDisplayName;
            end
            
            if dist > 0.3
                fieldName = sprintf('UnknownLine_RGB_%d_%d_%d_%s', ...
                    round(lineColor(1)*255), round(lineColor(2)*255), round(lineColor(3)*255), lineStyle);
                displayName = sprintf('Unknown line (%s)', lineStyle);
            end
            
            % Handle duplicates
            originalFieldName = fieldName;
            if isfield(allData.(currentSubplotName), fieldName)
                count = 2;
                while isfield(allData.(currentSubplotName), sprintf('%s_%d', originalFieldName, count))
                    count = count + 1;
                end
                fieldName = sprintf('%s_%d', originalFieldName, count);
            end
            
            % Store data
            allData.(currentSubplotName).(fieldName).XData = lines(i).XData;
            allData.(currentSubplotName).(fieldName).YData = lines(i).YData;
            allData.(currentSubplotName).(fieldName).Color = lines(i).Color;
            allData.(currentSubplotName).(fieldName).LineStyle = lines(i).LineStyle;
            allData.(currentSubplotName).(fieldName).LineWidth = lines(i).LineWidth;
            allData.(currentSubplotName).(fieldName).Marker = lines(i).Marker;
            allData.(currentSubplotName).(fieldName).MarkerSize = lines(i).MarkerSize;
            allData.(currentSubplotName).(fieldName).DisplayName = displayName;
            allData.(currentSubplotName).(fieldName).Type = 'line';
            
            fprintf('  [L%d] Added: %s (%d points, Style=%s, Color=[%.3f %.3f %.3f], dist=%.3f)\n', ...
                lineCount, fieldName, numPoints, lineStyle, lineColor(1), lineColor(2), lineColor(3), dist);
        end
    end
    
    % Process patches
    if ~isempty(patches)
        allData.(currentSubplotName).PilotInjectionDuration.XData = patches(1).XData;
        allData.(currentSubplotName).PilotInjectionDuration.YData = patches(1).YData;
        allData.(currentSubplotName).PilotInjectionDuration.FaceColor = patches(1).FaceColor;
        allData.(currentSubplotName).PilotInjectionDuration.EdgeColor = patches(1).EdgeColor;
        allData.(currentSubplotName).PilotInjectionDuration.FaceAlpha = patches(1).FaceAlpha;
        allData.(currentSubplotName).PilotInjectionDuration.DisplayName = 'Pilot injection duration';
        allData.(currentSubplotName).PilotInjectionDuration.Type = 'patch';
        fprintf('  Added: PilotInjectionDuration (patch)\n');
    end
    
    fprintf('  TOTAL STORED: %d lines, %d markers, %d patches\n', ...
        lineCount, markerCount, ~isempty(patches));
end

% =============== CALCULATE DIFFERENCES ===============
fprintf('\n========== CALCULATING DIFFERENCES ==========\n');

subplotFieldNames = fieldnames(allData);
for dataIdx = 1:length(subplotFieldNames)
    currentSubplotName = subplotFieldNames{dataIdx};
    currentData = allData.(currentSubplotName);
    fieldNames = fieldnames(currentData);
    
    % Find jet interaction and ignition delay data
    jetInteractionField = '';
    ignitionDelayField = '';
    
    for i = 1:length(fieldNames)
        if isfield(currentData.(fieldNames{i}), 'Type') && strcmp(currentData.(fieldNames{i}).Type, 'marker')
            if contains(fieldNames{i}, 'JetInteraction')
                jetInteractionField = fieldNames{i};
            elseif contains(fieldNames{i}, 'IgnitionDelay')
                ignitionDelayField = fieldNames{i};
            end
        end
    end
    
    fprintf('Subplot %s:\n', currentSubplotName);
    fprintf('  Jet Interaction field: %s\n', jetInteractionField);
    fprintf('  Ignition Delay field: %s\n', ignitionDelayField);
    
    if ~isempty(jetInteractionField) && ~isempty(ignitionDelayField)
        % Get the data
        jetInteractionTime = currentData.(jetInteractionField).XData;
        jetInteractionDistance = currentData.(jetInteractionField).YData;
        ignitionDelayTime = currentData.(ignitionDelayField).XData;
        ignitionDelayDistance = currentData.(ignitionDelayField).YData;
        
        % Calculate differences (simple since these are single points)
        timeDiff = ignitionDelayTime - jetInteractionTime;
        distanceDiff = ignitionDelayDistance - jetInteractionDistance;
        
        % Store the difference data
        allData.(currentSubplotName).Diff_jetinteraction_ID.TimeDifference = timeDiff;
        allData.(currentSubplotName).Diff_jetinteraction_ID.DistanceDifference = distanceDiff;
        allData.(currentSubplotName).Diff_jetinteraction_ID.AverageTimeDifference = timeDiff;
        allData.(currentSubplotName).Diff_jetinteraction_ID.AverageDistanceDifference = distanceDiff;
        allData.(currentSubplotName).Diff_jetinteraction_ID.JetInteractionTime = jetInteractionTime;
        allData.(currentSubplotName).Diff_jetinteraction_ID.JetInteractionDistance = jetInteractionDistance;
        allData.(currentSubplotName).Diff_jetinteraction_ID.IgnitionDelayTime = ignitionDelayTime;
        allData.(currentSubplotName).Diff_jetinteraction_ID.IgnitionDelayDistance = ignitionDelayDistance;
        allData.(currentSubplotName).Diff_jetinteraction_ID.DisplayName = 'Time difference (ID - JI)';
        allData.(currentSubplotName).Diff_jetinteraction_ID.Type = 'difference';
        
        fprintf('  Time difference: %.3f ms\n', timeDiff);
        fprintf('  Distance difference: %.3f mm\n', distanceDiff);
        fprintf('  ✓ Difference calculated and stored\n');
    else
        fprintf('  ✗ Missing marker data - cannot calculate difference\n');
    end
    fprintf('\n');
end

% =============== CREATE DIFFERENCE PLOT ===============
fprintf('========== CREATING DIFFERENCE ANALYSIS PLOT ==========\n');
figDiff = figure('Position', [200, 200, 1200, 400]);

% Summary comparison plot
subplot(1, 2, 1);
hold on; grid on; box on;

avgTimeDiffs = [];
subplotFieldNames = fieldnames(allData);

for i = 1:length(subplotFieldNames)
    currentSubplotName = subplotFieldNames{i};
    if isfield(allData.(currentSubplotName), 'Diff_jetinteraction_ID') && ~isempty(allData.(currentSubplotName).Diff_jetinteraction_ID.TimeDifference)
        avgTimeDiffs(i) = allData.(currentSubplotName).Diff_jetinteraction_ID.TimeDifference;
    else
        avgTimeDiffs(i) = NaN;
    end
end

x = 1:length(subplotFieldNames);
validData = ~isnan(avgTimeDiffs);

if any(validData)
    b = bar(x(validData), avgTimeDiffs(validData), 'FaceColor', [0.3 0.6 0.9], ...
            'EdgeColor', 'black', 'LineWidth', 1.5);
    
    for i = find(validData)
        text(i, avgTimeDiffs(i) + 0.05*sign(avgTimeDiffs(i)), ...
             sprintf('%.2f ms', avgTimeDiffs(i)), ...
             'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
    end
    
    xlabel('Subplot', 'FontSize', 12);
    ylabel('Time Difference [ms]', 'FontSize', 12);
    title('Time Difference: Ignition Delay - Jet Interaction', ...
          'FontSize', 12, 'FontWeight', 'bold');
    
    set(gca, 'XTick', 1:length(subplotFieldNames), 'XTickLabel', strrep(subplotFieldNames, '_', '-'));
    yline(0, '--k', 'LineWidth', 1);
    xlim([0.5 length(subplotFieldNames)+0.5]);
end

hold off;

% Plot timing comparison
subplot(1, 2, 2);
hold on; grid on; box on;

if any(validData)
    jetTimes = [];
    ignTimes = [];
    
    for i = 1:length(subplotFieldNames)
        currentSubplotName = subplotFieldNames{i};
        if isfield(allData.(currentSubplotName), 'Diff_jetinteraction_ID')
            jetTimes(i) = allData.(currentSubplotName).Diff_jetinteraction_ID.JetInteractionTime;
            ignTimes(i) = allData.(currentSubplotName).Diff_jetinteraction_ID.IgnitionDelayTime;
        else
            jetTimes(i) = NaN;
            ignTimes(i) = NaN;
        end
    end
    
    plot(x(validData), jetTimes(validData), 'o-', 'MarkerSize', 8, 'LineWidth', 2, ...
         'Color', [0 0.7 0], 'MarkerFaceColor', [0 0.7 0], 'DisplayName', 'Jet Interaction');
    
    plot(x(validData), ignTimes(validData), '*-', 'MarkerSize', 10, 'LineWidth', 2, ...
         'Color', [0.8 0 0], 'MarkerFaceColor', [0.8 0 0], 'DisplayName', 'Ignition Delay');
    
    xlabel('Subplot', 'FontSize', 12);
    ylabel('Time [ms]', 'FontSize', 12);
    title('Timing Comparison', 'FontSize', 12, 'FontWeight', 'bold');
    
    set(gca, 'XTick', 1:length(subplotFieldNames), 'XTickLabel', strrep(subplotFieldNames, '_', '-'));
    legend('Location', 'best');
    xlim([0.5 length(subplotFieldNames)+0.5]);
end

hold off;

sgtitle('Marker Analysis Results', 'FontSize', 14, 'FontWeight', 'bold');

% Verification
fprintf('\n========== DATA VERIFICATION ==========\n');
subplotFieldNames = fieldnames(allData);
for i = 1:length(subplotFieldNames)
    currentSubplotName = subplotFieldNames{i};
    fprintf('allData.%s: %d fields\n', currentSubplotName, length(fieldnames(allData.(currentSubplotName))));
    if ~isempty(fieldnames(allData.(currentSubplotName)))
        fprintf('  %s\n', strjoin(fieldnames(allData.(currentSubplotName)), ', '));
    end
end

% Create verification plot
fprintf('\n========== CREATING VERIFICATION PLOT ==========\n');
figVerify = figure('Position', [100, 100, 1400, 800]);

subplotFieldNames = fieldnames(allData);
numSubplots = length(subplotFieldNames);

% Determine subplot layout
if numSubplots <= 4
    rows = 2; cols = 2;
else
    rows = ceil(sqrt(numSubplots)); cols = ceil(numSubplots/rows);
end

for plotIdx = 1:numSubplots
    subplot(rows, cols, plotIdx);
    hold on; grid on; box on;
    
    currentSubplotName = subplotFieldNames{plotIdx};
    currentData = allData.(currentSubplotName);
    
    if isempty(fieldnames(currentData))
        title(sprintf('%s - NO DATA', strrep(currentSubplotName, '_', '-')));
        continue;
    end
    
    fieldNames = fieldnames(currentData);
    
    % Plot patches first
    patchHandles = [];
    for i = 1:length(fieldNames)
        if contains(fieldNames{i}, 'PilotInjectionDuration')
            h = patch(currentData.(fieldNames{i}).XData, currentData.(fieldNames{i}).YData, ...
                  currentData.(fieldNames{i}).FaceColor, ...
                  'EdgeColor', currentData.(fieldNames{i}).EdgeColor, ...
                  'FaceAlpha', currentData.(fieldNames{i}).FaceAlpha);
            patchHandles = [patchHandles, h];
        end
    end
    
    % Plot lines and markers and collect handles for legend
    lineHandles = [];
    legendEntries = {};
    
    for i = 1:length(fieldNames)
        if ~contains(fieldNames{i}, 'PilotInjectionDuration') && ~contains(fieldNames{i}, 'Diff_jetinteraction_ID')
            if isfield(currentData.(fieldNames{i}), 'Type')
                if strcmp(currentData.(fieldNames{i}).Type, 'marker')
                    % Plot marker
                    h = plot(currentData.(fieldNames{i}).XData, currentData.(fieldNames{i}).YData, ...
                         'Color', currentData.(fieldNames{i}).Color, ...
                         'Marker', currentData.(fieldNames{i}).Marker, ...
                         'MarkerSize', currentData.(fieldNames{i}).MarkerSize, ...
                         'LineStyle', 'none', ...
                         'MarkerFaceColor', currentData.(fieldNames{i}).MarkerFaceColor, ...
                         'MarkerEdgeColor', currentData.(fieldNames{i}).MarkerEdgeColor);
                else
                    % Plot line
                    h = plot(currentData.(fieldNames{i}).XData, currentData.(fieldNames{i}).YData, ...
                         'Color', currentData.(fieldNames{i}).Color, ...
                         'LineStyle', currentData.(fieldNames{i}).LineStyle, ...
                         'LineWidth', 1.5);
                end
            else
                % Backward compatibility
                h = plot(currentData.(fieldNames{i}).XData, currentData.(fieldNames{i}).YData, ...
                     'Color', currentData.(fieldNames{i}).Color, ...
                     'LineWidth', 1.5);
            end
            
            lineHandles = [lineHandles, h];
            
            % Use the stored display name
            if isfield(currentData.(fieldNames{i}), 'DisplayName')
                legendEntries{end+1} = currentData.(fieldNames{i}).DisplayName;
            else
                legendEntries{end+1} = strrep(fieldNames{i}, '_', ' ');
            end
        end
    end
    
    xlabel('Time [ms]', 'FontSize', 10);
    ylabel('Axial distance [mm]', 'FontSize', 10);
    title(strrep(currentSubplotName, '_', '-'), 'FontSize', 11, 'FontWeight', 'bold');
    xlim([0 15]);
    ylim([0 100]);
    
    % Create legend with proper handles
    if ~isempty(lineHandles)
        if ~isempty(patchHandles)
            legend([patchHandles, lineHandles], ['Pilot injection duration', legendEntries], ...
                'Location', 'southeast', 'FontSize', 8);
        else
            legend(lineHandles, legendEntries, 'Location', 'southeast', 'FontSize', 8);
        end
    end
    
    hold off;
end

sgtitle('Verification: Extracted Data (Lines + Markers) - Check if this matches the original', 'FontSize', 14, 'FontWeight', 'bold');

set(fig, 'Visible', 'on');

% Confirmation dialog
answer = questdlg('Does the extracted data match the original figure?', ...
    'Confirm Data Extraction', 'Yes - Save the data', 'No - Cancel', 'Yes - Save the data');

if strcmp(answer, 'Yes - Save the data')
    fprintf('\n========== SAVING DATA ==========\n');
    
    % Save the complete data structure using actual subplot names
    CompleteExtractedData = allData;  % allData is already properly structured
    
    % Save the complete data structure to the specified path
    saveFilePath = fullfile(saveFolder, 'Mark2024_Figure_CompleteExtractedData.mat');
    save(saveFilePath, 'CompleteExtractedData');
    
    fprintf('\n========== COMPLETE DATA STRUCTURE SAVED ==========\n');
    fprintf('File: %s\n', saveFilePath);
    fprintf('\nData Structure:\n');
    fprintf('CompleteExtractedData:\n');
    
    subplotFieldNames = fieldnames(CompleteExtractedData);
    for i = 1:length(subplotFieldNames)
        fprintf('  %s:\n', subplotFieldNames{i});
        dataFields = fieldnames(CompleteExtractedData.(subplotFieldNames{i}));
        for j = 1:length(dataFields)
            if isfield(CompleteExtractedData.(subplotFieldNames{i}).(dataFields{j}), 'DisplayName')
                displayName = CompleteExtractedData.(subplotFieldNames{i}).(dataFields{j}).DisplayName;
                fprintf('    %s (%s)\n', dataFields{j}, displayName);
            else
                fprintf('    %s\n', dataFields{j});
            end
        end
    end
    
    if ~isempty(subplotFieldNames)
        firstSubplot = subplotFieldNames{1};
        firstField = fieldnames(CompleteExtractedData.(firstSubplot));
        if ~isempty(firstField)
            fprintf('\nAccess example: CompleteExtractedData.%s.%s.XData\n', ...
                firstSubplot, firstField{1});
        end
    end
    
    msgbox(sprintf('Complete data structure saved!\n\nFile: Mark2024_Figure_CompleteExtractedData.mat\nLocation: %s\n\nStructure: CompleteExtractedData.SubplotName.LineType.data', saveFolder), 'Success');
else
    fprintf('\n========== SAVE CANCELLED ==========\n');
    msgbox('Save cancelled.', 'Cancelled', 'warn');
end

function [name, distance, displayName] = findClosestColor(testColor, colorMap)
    minDist = Inf;
    name = 'Unknown';
    displayName = 'Unknown';
    for i = 1:length(colorMap)
        dist = norm(testColor - colorMap{i}.rgb);
        if dist < minDist
            minDist = dist;
            name = colorMap{i}.name;
            displayName = colorMap{i}.desc;
        end
    end
    distance = minDist;
end