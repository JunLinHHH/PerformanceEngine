% Generate Statistical Analysis for:
    % Pressure
    % Heat release rate (HRR)
    % Ignition delay
    
% Input:   DataTable (T_DataAll), DataSheet(.xlsx file)
% Output: .png/.pdf image

%clearvars; close all; clc
%% Process Data
%% -- Load data
% If .mat exist
clearvars; close all; clc; 
Configs = Configuration();
Loaded = load(Configs.DataObjReadDir);

obj = Loaded.obj;
T_DataAll = obj.DataMatrix;
T_Condition_ROI = obj.ConditionMatrix;
%% -- Process Mean and STD
T_Condition_ROI.Time_Scale                     = cell(height(T_Condition_ROI), 1);
T_Condition_ROI.P_Mean                          = cell(height(T_Condition_ROI), 1);
T_Condition_ROI.P_Std                           = cell(height(T_Condition_ROI), 1);
T_Condition_ROI.P_Indie                         = cell(height(T_Condition_ROI), 1);

T_Condition_ROI.HRR_Mean                        = cell(height(T_Condition_ROI), 1);
T_Condition_ROI.HRR_Std                         = cell(height(T_Condition_ROI), 1);
T_Condition_ROI.HRR_Indie                       = cell(height(T_Condition_ROI), 1);

T_Condition_ROI.PressureIgnitionDelay_Mean      = cell(height(T_Condition_ROI), 1);
T_Condition_ROI.PressureIgnitionDelay_Std       = cell(height(T_Condition_ROI), 1);
T_Condition_ROI.PressureIgnitionDelay_Indie     = cell(height(T_Condition_ROI), 1);

T_Condition_ROI.IntensityIgnitionDelay_Mean     = cell(height(T_Condition_ROI), 1);
T_Condition_ROI.IntensityIgnitionDelay_Std      = cell(height(T_Condition_ROI), 1);
T_Condition_ROI.IntensityIgnitionDelay_Indie    = cell(height(T_Condition_ROI), 1);

T_Condition_ROI.HRRIgnitionDelay_Mean           = cell(height(T_Condition_ROI), 1);
T_Condition_ROI.HRRIgnitionDelay_Std            = cell(height(T_Condition_ROI), 1);
T_Condition_ROI.HRRIgnitionDelay_Indie          = cell(height(T_Condition_ROI), 1);

T_Condition_ROI.ShiftedHRR_PressureID_Mean      = cell(height(T_Condition_ROI), 1);
T_Condition_ROI.ShiftedHRR_PressureID_Std       = cell(height(T_Condition_ROI), 1);
T_Condition_ROI.ShiftedHRR_PressureID_Indie     = cell(height(T_Condition_ROI), 1);

T_Condition_ROI.ShiftedHRR_IntensityID_Mean     = cell(height(T_Condition_ROI), 1);
T_Condition_ROI.ShiftedHRR_IntensityID_Std      = cell(height(T_Condition_ROI), 1);
T_Condition_ROI.ShiftedHRR_IntensityID_Indie    = cell(height(T_Condition_ROI), 1);

T_Condition_ROI.Photodiode_Mean                 = cell(height(T_Condition_ROI), 1);
T_Condition_ROI.Photodiode_Std                  = cell(height(T_Condition_ROI), 1);
T_Condition_ROI.Photodiode_Indie                = cell(height(T_Condition_ROI), 1);

T_Condition_ROI.ShiftedPD_Mean                  = cell(height(T_Condition_ROI), 1);
T_Condition_ROI.ShiftedPD_Std                   = cell(height(T_Condition_ROI), 1);
T_Condition_ROI.ShiftedPD_Indie                 = cell(height(T_Condition_ROI), 1);


%% Get Mean and Std for Pressure, HRR and Ignition delay

for ii = 1: height(T_Condition_ROI) 

    ConditionTag = Class_TableManagement.ConditionTagGenerator(T_Condition_ROI, obj, ii);
    runsForThisConsition =  T_Condition_ROI.RunNumbers{ii};

    fprintf('\nCONDITION: #%02d\n',ii);
    fprintf([ConditionTag,'\n']);
    fprintf('\t\tValid %03d runs\n', length(runsForThisConsition));
    RowsForThisCondition = T_DataAll(runsForThisConsition, :);  % Non discarded

    % Calculate Mean and STD
    if ~isempty(RowsForThisCondition)
        % Time scale
        for rr = 1:height(RowsForThisCondition)
             temp_time_scale = RowsForThisCondition.TimeScale{rr};
             if ~isempty(temp_time_scale)
                    T_Condition_ROI.Time_Scale{ii} = temp_time_scale*1000;  % unit: ms
             end
        end

        % Pressure
        fprintf('Calculating >>> Pressure <<< mean and std\n')
        [T_Condition_ROI.P_Mean{ii}, T_Condition_ROI.P_Std{ii}] = ...
            StatisticFromCell(RowsForThisCondition.P_Corrected); % RowsForThisCondition.P_Corrected: cell array
        T_Condition_ROI.P_Indie{ii} = RowsForThisCondition.P_Corrected;
    
        % HRR
        fprintf('Calculating >>> HRR <<< mean and std\n')
        [T_Condition_ROI.HRR_Mean{ii}, T_Condition_ROI.HRR_Std{ii}] = ...
            StatisticFromCell(RowsForThisCondition.HRR);
        T_Condition_ROI.HRR_Indie{ii} = RowsForThisCondition.HRR;
    
        % Pressure Ignition Delay
        fprintf('Calculating >>> Pressure Ignition Delay <<< mean and std\n')
        [T_Condition_ROI.PressureIgnitionDelay_Mean{ii}, T_Condition_ROI.PressureIgnitionDelay_Std{ii}] = ...
            StatisticFromCell(RowsForThisCondition.IgnitionDelay_Pressure);
        T_Condition_ROI.PressureIgnitionDelay_Indie{ii} = RowsForThisCondition.IgnitionDelay_Pressure;
    
        % HRR Ignition Delay
        fprintf('Calculating >>> HRR Ignition Delay <<< mean and std\n')
        [T_Condition_ROI.HRRIgnitionDelay_Mean{ii}, T_Condition_ROI.HRRIgnitionDelay_Std{ii}] = ...
            StatisticFromCell(RowsForThisCondition.IgnitionDelay_HRR);
        T_Condition_ROI.HRRIgnitionDelay_Indie{ii} = RowsForThisCondition.IgnitionDelay_HRR;
    
        % Intensity Ignition delay
        fprintf('Calculating >>> Intensity Ignition Delay <<< mean and std\n')
        [T_Condition_ROI.IntensityIgnitionDelay_Mean{ii}, T_Condition_ROI.IntensityIgnitionDelay_Std{ii}] = ...
            StatisticFromCell(RowsForThisCondition.IgnitionDelay_Intensity);
        T_Condition_ROI.IntensityIgnitionDelay_Indie{ii} = RowsForThisCondition.IgnitionDelay_Intensity;
    
        % Photodiode
        fprintf('Calculating >>> Photodiode <<< mean and std\n')
        [T_Condition_ROI.Photodiode_Mean{ii}, T_Condition_ROI.Photodiode_Std{ii}] = ...
            StatisticFromCell(RowsForThisCondition.PD_Signal);
        T_Condition_ROI.Photodiode_Indie{ii} = RowsForThisCondition.PD_Signal;
    
        % Shift HRR
        fprintf('Calculating >>> Shift HRR (Pressure ID) <<< mean and std\n')
        RowInput = T_Condition_ROI(ii, :);

        ShiftedHRR_PressureID_Indie = ShiftHRR_PressureIgnitionDelay(RowInput);

        [T_Condition_ROI.ShiftedHRR_PressureID_Mean{ii}, T_Condition_ROI.ShiftedHRR_PressureID_Std{ii}] = ...
            StatisticFromCell(ShiftedHRR_PressureID_Indie);
        T_Condition_ROI.ShiftedHRR_PressureID_Indie{ii} = ShiftedHRR_PressureID_Indie;
        
        fprintf('Calculating >>> Shift HRR (Intensity ID) <<< mean and std\n')
        ShiftedHRR_IntensityID_Indie = ShiftHRR_IntensityIgnitionDelay(RowInput);

        [T_Condition_ROI.ShiftedHRR_IntensityID_Mean{ii}, T_Condition_ROI.ShiftedHRR_IntensityID_Std{ii}] = ...
            StatisticFromCell(ShiftedHRR_IntensityID_Indie);
        T_Condition_ROI.ShiftedHRR_IntensityID_Indie{ii} = ShiftedHRR_IntensityID_Indie;

        % Shift PD
        fprintf('Calculating >>> Shift PD (Pressure ID) <<< mean and std\n')
        Shifted_PD_Indie = ShiftPD(RowInput);
        [T_Condition_ROI.ShiftedPD_Mean{ii}, T_Condition_ROI.ShiftedPD_Std{ii}] = ...
            StatisticFromCell(Shifted_PD_Indie);
        T_Condition_ROI.ShiftedPD_Indie{ii} = Shifted_PD_Indie;        
    end
    pause(0.5)
    close all

end

% Overwrite_ID_Mean = dataTable(:, 15);
% Overwrite_ID_Std = dataTable(:, 16);

%% Plot and save figures
% Three secionts below
% Can run any of them after the previous steps
% Can show run number 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ¯\_( ͡° ͜ʖ ͡°)_/¯✊   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ¯\_( ͡° ͜ʖ ͡°)_/¯✊   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ¯\_( ͡° ͜ʖ ͡°)_/¯✊   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ¯\_( ͡° ͜ʖ ͡°)_/¯✊   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -- Shared properties 
% FigureOutput = fullfile(pwd,'E:\Outputs\figure');
FigureOutput = Configs.CommonFiguresOutput;
Color = {'k', [0 0.4470 0.7410],[0.6350 0.0780 0.1840],[0.9290 0.6940 0.1250]};
LineWidth = 2;
if ~isempty(RowsForThisCondition)
    TimeScale = RowsForThisCondition.('TimeScale'){1}*1000;  % unit: ms
else
    error('The sub set is empty, cannot get TimeScale reference')
end
TimeScaleClosedLoop = [TimeScale,fliplr(TimeScale)]; 

%%
keyboard
obj.ConditionMatrix = T_Condition_ROI;
DataObjSaveDir = Configs.DataObjSaveDir;     
save(DataObjSaveDir,"obj")
fprintf("Saving finished\n")
%% Internal functions
%% -- Click callback
function ButtonDownCallBack(hObj,~,ColorNonActive) % src is ~

    thisFig = gcf;
    SelectType = thisFig.SelectionType;

    if strcmp(SelectType,'normal')
        fprintf('Clicked data is #%s\n',hObj.XDataSource)
        hObj.Color = [1,0,0];
    elseif strcmp(SelectType,'alt')
        hObj.Color = ColorNonActive;
    end
end

%% -- Calculate mean and std from a cell
function [MeanData,StdData] = StatisticFromCell(C_Input)
    MeanData = [];
    StdData = [];
    % Check if the input cell is empty
    Flags = ~cellfun(@isempty, C_Input);

    if sum(Flags) == 0
        fprintf("- - - Input cell is empty, return empty\n")
        return
    end
    
    % Find the data size
    idx = find(Flags == 1,1,'first');
    SumData = zeros(numel(C_Input), length(C_Input{idx}));
    for ii = 1:numel(C_Input)
        if ~isempty(C_Input{ii})
            SumData(ii,:) = C_Input{ii};
        else
            SumData(ii,:) = NaN;
            % fprintf('Empty...\n')
        end
    end
    MeanData = mean(SumData,1,'omitnan');
    StdData = std(SumData,1,'omitnan');

%     f = figure; grid on; hold on; box on; set(f,'color','white')
%     for ii = 1:numel(C_Input)
%         plot(C_Input{ii});
%     end
%     plot(MeanData,'color','k','LineWidth',3,'LineStyle','--','Marker','o','MarkerSize',1);
end

%% -- Shift HRR
function newHRR_Indie = ShiftHRR_PressureIgnitionDelay(DataPassIn)
    SampingRate = 200000/1000; % sampling per ms
    IgnitionDelay_Mean = DataPassIn.PressureIgnitionDelay_Mean{1};
    IgnitionDelay_Indie = DataPassIn.PressureIgnitionDelay_Indie{1};
    HRR_Indie = DataPassIn.HRR_Indie{1};
    newHRR_Indie = cell(size(HRR_Indie));
    RunAmount = length(IgnitionDelay_Indie);
    if  isempty(IgnitionDelay_Mean)
            fprintf("- - - Input IgnitionDelay_Mean is empty, return empty\n")
            return
    end
    % Circshift HRR
    for jj = 1:RunAmount
            thisIgnitionDelay = IgnitionDelay_Indie{jj};
            thisHRR_Indie = HRR_Indie{jj};
            if  isempty(thisIgnitionDelay) || isempty(thisHRR_Indie)
                continue
            end  
            ID_difference = floor((IgnitionDelay_Mean - thisIgnitionDelay) * SampingRate);
            newHRR = circshift(thisHRR_Indie, ID_difference);   
            newHRR_Indie{jj} = newHRR;
    end
end

function newHRR_Indie = ShiftHRR_IntensityIgnitionDelay(DataPassIn)
    SampingRate = 200000/1000; % sampling per ms
    IgnitionDelay_Mean = DataPassIn.IntensityIgnitionDelay_Mean{1};
    IgnitionDelay_Indie = DataPassIn.IntensityIgnitionDelay_Indie{1};
    HRR_Indie = DataPassIn.HRR_Indie{1};
    newHRR_Indie = cell(size(HRR_Indie));
    RunAmount = length(IgnitionDelay_Indie);
    if  isempty(IgnitionDelay_Mean)
            fprintf("- - - Input IgnitionDelay_Mean is empty, return empty\n")
            return
    end
    % Circshift HRR
    for jj = 1:RunAmount
            thisIgnitionDelay = IgnitionDelay_Indie{jj};
            thisHRR_Indie = HRR_Indie{jj};
            if  isempty(thisIgnitionDelay) || isempty(thisHRR_Indie)
                continue
            end  
            ID_difference = floor((IgnitionDelay_Mean - thisIgnitionDelay) * SampingRate);
            newHRR = circshift(thisHRR_Indie, ID_difference);   
            newHRR_Indie{jj} = newHRR;
    end
end

%% -- Shift PD
function newPD_Indie = ShiftPD(DataPassIn)
    SampingRate = 200000/1000; % sampling per ms
    IgnitionDelay_Mean = DataPassIn.PressureIgnitionDelay_Mean{1};
    IgnitionDelay_Indie = DataPassIn.PressureIgnitionDelay_Indie{1};
    PD_Indie = DataPassIn.Photodiode_Indie{1};
    newPD_Indie = cell(size(PD_Indie));
    RunAmount = length(IgnitionDelay_Indie);
    if  isempty(IgnitionDelay_Mean)
            fprintf("- - - Input IgnitionDelay_Mean is empty, return empty\n")
            return
    end
    % Circshift PD
    for jj = 1:RunAmount
            thisIgnitionDelay = IgnitionDelay_Indie{jj};
            thisPD_Indie = PD_Indie{jj};
            if  isempty(thisIgnitionDelay) || isempty(thisPD_Indie)
                continue
            end  
            ID_difference = floor((IgnitionDelay_Mean - thisIgnitionDelay) * SampingRate);
            newPD = circshift(thisPD_Indie, ID_difference);   
            newPD_Indie{jj} = newPD;
    end
end

%% -- CreateRectangle For 4 elements
function CreateRectangle_4Elements(HandleFigure, L_position, offset)
    %CREATERECTANGLE(figure1)
    %  FIGURE1:  annotation figure

    %  Auto-generated by MATLAB on 08-Jul-2021 16:15:43

    % Create rectangle
    color_shade = {'r','b','k','g'};
    Step_height = 12.5;

    %color_shade = fliplr(color_shade);
    L_position = L_position + offset;

    L_base = L_position(2);
    stepY = 14;
    for ii = 1:length(color_shade)

        L_position(2) = L_base + stepY*(ii-1);
        L_position(4) = Step_height;
        
        annotation(HandleFigure,'rectangle',...
            'unit','pixel','position',...
            L_position,...
            'LineStyle','none',...
            'FaceColor', color_shade{ii},...
            'FaceAlpha',0.2);
    end
end
%% -- CreateRectangle For 3 elements
function CreateRectangle_3Element(HandleFigure, L_position, offset)
    %CREATERECTANGLE(figure1)
    %  FIGURE1:  annotation figure

    %  Auto-generated by MATLAB on 08-Jul-2021 16:15:43

    % Create rectangle
    color_shade = {'r','b','k'};

    %color_shade = fliplr(color_shade);
    L_position = L_position + offset;

    L_base = L_position(2);
    stepY = 14;
    for ii = 1:length(color_shade)

        L_position(2) = L_base + stepY*(ii-1);
        
        annotation(HandleFigure,'rectangle',...
            'unit','pixel','position',...
            L_position,...
            'LineStyle','none',...
            'FaceColor', color_shade{ii},...
            'FaceAlpha',0.2);
    end
end