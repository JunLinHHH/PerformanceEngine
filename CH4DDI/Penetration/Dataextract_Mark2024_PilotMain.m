clear all; clc; close all
% Load the figure
fig = openfig('C:\Users\jun-y\OneDrive - UNSW\CH4DDI\Mark2024 CH4DDI CVCC\JournalFigures\Figure_PenetrationRecession_Pilot_Main.fig');

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

% Define subplot names
subplotNames = {'D_0_03ms_M', 'D_1_03ms_M', 'D_2_03ms_M', 'D_3_03ms_M'};
subplotTitles = {'D-0.03ms-M', 'D-1.03ms-M', 'D-2.03ms-M', 'D-3.03ms-M'};

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

% Initialize data storage
allData = cell(4, 1);
for i = 1:4
    allData{i} = struct();
end

fprintf('\n========== PROCESSING SUBPLOTS ==========\n');

% Process each D-case
for dataIdx = 1:min(4, length(dCaseAxes))
    axIdx = dCaseAxes(dataIdx);
    currentAx = allAxes(axIdx);
    
    fprintf('\n--- Processing D-case %d: Axes index %d (%s) ---\n', ...
        dataIdx, axIdx, subplotNames{dataIdx});
    
    % Get all lines and patches
    lines = findall(currentAx, 'Type', 'line');
    patches = findall(currentAx, 'Type', 'patch');
    
    fprintf('  Total lines: %d, Total patches: %d\n', length(lines), length(patches));
    
    % Process lines
    lineCount = 0;
    for i = 1:length(lines)
        numPoints = length(lines(i).XData);
        lineStyle = lines(i).LineStyle;
        lineColor = lines(i).Color;
        
        % Skip short lines (legend entries)
        if numPoints <= 2
            continue;
        end
        
        % Only solid lines
        if ~strcmp(lineStyle, '-')
            continue;
        end
        
        lineCount = lineCount + 1;
        
        % Match color
        [fieldName, dist, displayName] = findClosestColor(lineColor, colorMap);
        
        if dist > 0.3
            fieldName = sprintf('UnknownLine_RGB_%d_%d_%d', ...
                round(lineColor(1)*255), round(lineColor(2)*255), round(lineColor(3)*255));
            displayName = fieldName;
        end
        
        % Handle duplicates
        originalFieldName = fieldName;
        if isfield(allData{dataIdx}, fieldName)
            count = 2;
            while isfield(allData{dataIdx}, sprintf('%s_%d', originalFieldName, count))
                count = count + 1;
            end
            fieldName = sprintf('%s_%d', originalFieldName, count);
        end
        
        % Store data
        allData{dataIdx}.(fieldName).XData = lines(i).XData;
        allData{dataIdx}.(fieldName).YData = lines(i).YData;
        allData{dataIdx}.(fieldName).Color = lines(i).Color;
        allData{dataIdx}.(fieldName).LineStyle = lines(i).LineStyle;
        allData{dataIdx}.(fieldName).LineWidth = lines(i).LineWidth;
        allData{dataIdx}.(fieldName).Marker = lines(i).Marker;
        allData{dataIdx}.(fieldName).MarkerSize = lines(i).MarkerSize;
        allData{dataIdx}.(fieldName).DisplayName = displayName; % Store display name
        
        fprintf('  [%d] Added: %s (%d points, Color=[%.3f %.3f %.3f], dist=%.3f)\n', ...
            lineCount, fieldName, numPoints, lineColor(1), lineColor(2), lineColor(3), dist);
    end
    
    % Process patches
    if ~isempty(patches)
        allData{dataIdx}.PilotInjectionDuration.XData = patches(1).XData;
        allData{dataIdx}.PilotInjectionDuration.YData = patches(1).YData;
        allData{dataIdx}.PilotInjectionDuration.FaceColor = patches(1).FaceColor;
        allData{dataIdx}.PilotInjectionDuration.EdgeColor = patches(1).EdgeColor;
        allData{dataIdx}.PilotInjectionDuration.FaceAlpha = patches(1).FaceAlpha;
        allData{dataIdx}.PilotInjectionDuration.DisplayName = 'Pilot injection duration';
        fprintf('  Added: PilotInjectionDuration (patch)\n');
    end
    
    fprintf('  TOTAL STORED FIELDS: %d\n', length(fieldnames(allData{dataIdx})));
end

% Verification
fprintf('\n========== DATA VERIFICATION ==========\n');
for i = 1:4
    fprintf('allData{%d} (%s): %d fields\n', i, subplotNames{i}, length(fieldnames(allData{i})));
    if ~isempty(fieldnames(allData{i}))
        fprintf('  %s\n', strjoin(fieldnames(allData{i}), ', '));
    end
end

% Create verification plot
fprintf('\n========== CREATING VERIFICATION PLOT ==========\n');
figVerify = figure('Position', [100, 100, 1400, 800]);

for plotIdx = 1:4
    subplot(2, 2, plotIdx);
    hold on; grid on; box on;
    
    currentData = allData{plotIdx};
    
    if isempty(fieldnames(currentData))
        title(sprintf('%s - NO DATA', subplotTitles{plotIdx}));
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
    
    % Plot lines and collect handles for legend
    lineHandles = [];
    legendEntries = {};
    
    for i = 1:length(fieldNames)
        if ~contains(fieldNames{i}, 'PilotInjectionDuration')
            h = plot(currentData.(fieldNames{i}).XData, currentData.(fieldNames{i}).YData, ...
                 'Color', currentData.(fieldNames{i}).Color, ...
                 'LineWidth', 1.5);
            
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
    title(subplotTitles{plotIdx}, 'FontSize', 11, 'FontWeight', 'bold');
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

sgtitle('Verification: Extracted Data - Check if this matches the original', 'FontSize', 14, 'FontWeight', 'bold');

set(fig, 'Visible', 'on');

% Confirmation dialog
answer = questdlg('Does the extracted data match the original figure?', ...
    'Confirm Data Extraction', 'Yes - Save the data', 'No - Cancel', 'Yes - Save the data');

if strcmp(answer, 'Yes - Save the data')
    fprintf('\n========== SAVING DATA ==========\n');
    
    D_0_03ms_M = allData{1};
    save(fullfile(saveFolder, 'Mark2024_Figure_extractedData_D_0_03ms_M.mat'), 'D_0_03ms_M');
    fprintf('  ✓ Saved: D_0_03ms_M\n');
    
    D_1_03ms_M = allData{2};
    save(fullfile(saveFolder, 'Mark2024_Figure_extractedData_D_1_03ms_M.mat'), 'D_1_03ms_M');
    fprintf('  ✓ Saved: D_1_03ms_M\n');
    
    D_2_03ms_M = allData{3};
    save(fullfile(saveFolder, 'Mark2024_Figure_extractedData_D_2_03ms_M.mat'), 'D_2_03ms_M');
    fprintf('  ✓ Saved: D_2_03ms_M\n');
    
    D_3_03ms_M = allData{4};
    save(fullfile(saveFolder, 'Mark2024_Figure_extractedData_D_3_03ms_M.mat'), 'D_3_03ms_M');
    fprintf('  ✓ Saved: D_3_03ms_M\n');
    
    fprintf('\n========== ALL FILES SAVED ==========\n');
    fprintf('Location: %s\n', saveFolder);
    msgbox('Data successfully saved!', 'Success');
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