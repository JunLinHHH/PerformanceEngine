% Generate video matrix for all runs at the same test conditions 

% Input:  DataTable (T_DataAll), DataSheet(.xlsx file)
% Output: Video .avi
% Requires load('MethaneAutoignitionDataSet.mat');

% If .mat exist
clearvars; close all; clc; 
Configs = Configuration();
Loaded = load(Configs.DataObjReadDir);

obj = Loaded.obj;
T_DataAll = obj.DataMatrix;
T_Condition_ROI = obj.ConditionMatrix;

%% Process pressure signal (.tdms)
% Process all data (discarded run will be skipped by using try-catch) 
% Pressure data processing for selected rows

run_selected = 1:height(T_DataAll);
bShowPlot = 0;

T_DataAll = Table_PressureProcess(T_DataAll, run_selected, bShowPlot);
T_DataAll = Table_PhotodiodeProcess(T_DataAll, run_selected, bShowPlot);
%T_DataAll = Table_IntensityProfile(T_DataAll, run_selected, bShowPlot);
T_DataAll = Table_HRRIgnitionDelay(T_DataAll, run_selected);
obj = obj.SetDataMatrix(T_DataAll); % Put modified Table to replace the one in object
fprintf("Process finished\n")

%%
% obj= obj.RecoverRun([191, 192]);
% T_DataAll = obj.DataMatrix;
% T_Condition_ROI = obj.ConditionMatrix;
%%
keyboard
DataObjSaveDir = Configs.DataObjSaveDir;     
save(DataObjSaveDir,"obj")
fprintf("Saving finished\n")
%% Table processing functions
% Pressure
function T_DataAll = Table_PressureProcess(T_DataAll, run_selected, plotCheck)
    for ii = run_selected
        try
            fprintf('# Processing (Pressure) index %03d ...\n', ii);
            this_row = T_DataAll(T_DataAll.Number==ii,:);
            output_row = this_row;
            dir_tdms = this_row.File_Pressure{1};
            Obj_Pressure = CLASS_PressureProcess();

            % Injection Delay
            Obj_Pressure.InjectionDelay = this_row.Inj_delay;

            % Toggle checking plot 
                Obj_Pressure.bPlotIgnitionDelay = plotCheck;
                Obj_Pressure.bPlotPressureCorrected = plotCheck;
            
            Obj_Pressure = Obj_Pressure.PressureProcess(dir_tdms);
            PressureInfo = Obj_Pressure.GetPressureInfo;
    
            % Check if the table ambient conditions are the same as the ones in log file
            if ~isempty(PressureInfo.TriggerPressureSet) % This means .log reading above failed 
                try
                    assert(PressureInfo.TriggerPressureSet == this_row.('P_trig'));
                    assert(PressureInfo.O2Percentage == this_row.('O2'));
                catch
                    warning('Condition assertion failed')
                    fprintf("Log file pressure: %d\n",PressureInfo.TriggerPressureSet)
                    fprintf("Table pressure: %d\n", this_row.('P_trig'))
                    fprintf("Log file O2 percentage: %d\n",PressureInfo.O2Percentage)
                    fprintf("Table O2 percentage: %d\n", this_row.('O2'))                    
                end
            end
            
            output_row.TimeScale = {PressureInfo.TimeScale};
            output_row.TimeScale = {PressureInfo.TimeScale};
            output_row.P_Raw = {PressureInfo.Pressure};
            output_row.P_Corrected = {PressureInfo.P_Corrected};
            output_row.HRR = {PressureInfo.HRR};
            output_row.IgnitionDelay_Pressure = {PressureInfo.IgnitionDelay};
            output_row.P_TrigTrue = {PressureInfo.TriggerPressureTrue};
            output_row.Trigger_Signal = {PressureInfo.Trigger_Signal};

            T_DataAll(T_DataAll.Number==ii,:) = output_row;
        catch
            warning('# Error for the run with index %03d ...', ii);
        end
    end
end

% Photodiode
function T_DataAll = Table_PhotodiodeProcess(T_DataAll, run_selected, plotCheck)
    for ii = run_selected
        try
            fprintf('# Processing (Photodiode) index %03d ...\n', ii);
            this_row = T_DataAll(T_DataAll.Number==ii,:);
            output_row = this_row;
            dir_tdms = this_row.File_Photodiode{1};
            Gain_PD = this_row.Gain(1);
            % Injection Delay
   
            Obj_Photodiode = CLASS_PhotodiodeProcess();
            Obj_Photodiode.Gain_PD = Gain_PD;
            Obj_Photodiode.InjectionDelay = this_row.Inj_delay;      
            % Toggle checking plot 
                % NOT Implemented
                % Obj_Photodiode.bPlotIgnitionDelay = plotCheck;
                % Obj_Photodiode.bPlotPressureCorrected = plotCheck;

            Obj_Photodiode = Obj_Photodiode.PhotodiodeProcess(dir_tdms);
            PhotodiodeInfo = Obj_Photodiode.GetPhotodiodeInfo;
            
            output_row.PD_TimeScale = {PhotodiodeInfo.TimeScale};
            output_row.PD_Signal = {PhotodiodeInfo.Photodiode};
            T_DataAll(T_DataAll.Number==ii,:) = output_row;
        catch
            warning('# Error for the run with index %03d ...', ii);
        end
    end
end

% Intensity
function T_DataAll = Table_IntensityProfile(T_DataAll, run_selected, checkPlot)
    for ii = run_selected
            this_row = T_DataAll(T_DataAll.Number==ii,:);
            
            fprintf('# Processing (Schlieren Intensity) index %03d ...\n', ii);
            O2 = this_row.O2;
            P_trigger = this_row.P_trig;
            Type = this_row.Type;
            Type= Type{1};
            
            if O2 == 21 && P_trigger == 52 && strcmp(Type, 'AutoIgnitionHeptane')
                threshold = 10000;
            elseif O2 == 15 && P_trigger == 52 && strcmp(Type, 'AutoIgnitionHeptane')
                threshold = 5000;
            else
                threshold = 30000;
            end

            output_row = IgnitionDelayFromIntensity(this_row, 'ByCriticalPoint', threshold, checkPlot);
            T_DataAll(T_DataAll.Number==ii,:) = output_row;
    end
end

% HRR
function T_DataAll = Table_HRRIgnitionDelay(T_DataAll, run_selected)
    for ii = run_selected
            HRR_threshold_percent = 0.3;
            fprintf('# Processing (HRR ignition delay) index %03d ...\n', ii);
            this_row = T_DataAll(T_DataAll.Number == ii,:);
            if strcmp(this_row.Type,  'DualFuel')
                HRR_threshold_percent = 0.17;
            elseif strcmp(this_row.Type,  'AutoIgnitionHeptane') && this_row.P_trig == 52 && this_row.O2 == 21
                HRR_threshold_percent = 0.7;
            elseif strcmp(this_row.Type,  'AutoIgnitionHeptane') && this_row.P_trig == 63
                HRR_threshold_percent = 0.3;
            elseif strcmp(this_row.Type,  'AutoIgnitionHeptane') && this_row.P_trig == 44
                HRR_threshold_percent = 0.9;     
            elseif strcmp(this_row.Type,  'AutoIgnitionHeptane') && this_row.O2 == 15
                HRR_threshold_percent = 0.8;  
            elseif strcmp(this_row.Type,  'AutoIgnitionHeptane') && this_row.O2 == 10
                HRR_threshold_percent = 0.55;  
            end
            output_row = IgnitionDelayFromHRR(this_row, HRR_threshold_percent);
            T_DataAll(T_DataAll.Number==ii,:) = output_row;
    end
end

%% Process for each individuial row
% Process for each individuial file
function output_row = IgnitionDelayFromIntensity(input_row, method, threshold, checkPlot)
        output_row = input_row;
        imageDir = input_row.File_Image;
        FPS = input_row.FPS;
        FPms = FPS/1000;
        injectionDelay = input_row.Inj_delay;

        obj = CLASS_MultiPageTif( imageDir{1});
        obj = obj.CalculateIntensityProfile();

        
        intensityProfile = obj.GetIntensityProfile();
        % method  - by critical point
        if strcmp(method,'ByCriticalPoint')
            [~ ,index_min] = min(intensityProfile);
            [logical_local_max,~] = islocalmax(intensityProfile,'MinProminence',threshold); % manual check this threshold
            xxx = 1:length(intensityProfile);
            index_local_max = xxx(logical_local_max);

            % FOR CHECK
            figure; plot(xxx/FPms,intensityProfile, xxx(index_local_max)/FPms,intensityProfile(index_local_max),'ro')

            index_local_max(index_local_max>index_min) = [];
            index_ignition = index_local_max(end-1);
        elseif strcmp(method,'ByThreshold')
            intensityProfile_diff = -diff(intensityProfile) / max(-diff(intensityProfile));
            index_ignition = find(intensityProfile_diff>threshold,1,'first'); 
        end
        %
        ignDelay_Intensity = index_ignition/FPms - injectionDelay; % !!!! in ms
        obj = obj.StoreTifToCell(index_ignition:index_ignition+2);
        imageFrame = obj.GetTifStoreCell();
        instantImage = imageFrame{1};

        output_row.Intensity_Profile = {intensityProfile};        
        output_row.IgnitionDelay_Intensity = {ignDelay_Intensity};
        output_row.Frame_At_Ignition = {index_ignition};

        this_HRR = input_row.HRR{1};
        this_time_HRR_sec = input_row.TimeScale{1};

        if checkPlot
            obj_ax = CLASS_AxesHandleStore();
            obj_ax.RowNumber = 4;
            obj_ax.ColumnNumber = 1;
            obj_ax.AxesWidth = 800/2;            % width of the axes
            obj_ax.AxesHeight = 300/2;            % height if the axes
            obj_ax.GapRow = 30;            % gap between columns
            obj_ax.MarginTop = 30;
            obj_ax = obj_ax.ConstuctAxes();
            handles_ax = obj_ax.GetAxesHandleMatrix();
            handles_fig = obj_ax.GetFigureHandle();
            [~,file_name,~] = fileparts(imageDir{1});
            handles_fig.Name =file_name;
            hold(handles_ax(1),'on');
            
            plot(handles_ax(1), ((1:length(intensityProfile))-1), intensityProfile);
            xline(handles_ax(1),ignDelay_Intensity)
            xline(handles_ax(1),index_ignition,'r')
            % set(handles_ax(1),'YLim',[0,1.1]);
            set(handles_ax(1),'XLim',[0,50]);
            % plot(handles_ax(1),this_time_HRR_sec*1000, this_HRR/max(this_HRR));
            plot(handles_ax(1), xxx(index_local_max),intensityProfile(index_local_max),'ro')
            
            imshow(instantImage,'Parent',handles_ax(2));
            imshow(imageFrame{2},'Parent',handles_ax(3));
            imshow(imageFrame{3},'Parent',handles_ax(4));
            CLASS_Utilis.InsertFigureText(handles_ax(1), 0.5, 0.8, sprintf(('Frame: %04d, ign: %0.2f ms'),index_ignition, ignDelay_Intensity));
            
            % For save
            % [~,bb,~] = fileparts(imageDir);
            % print(gcf,sprintf(('E:/Outputs/figure/IntensityPlot/IntensityPlot_%s'),bb),'-dpng')
            % close all
        end
end

function output_row = IgnitionDelayFromHRR(input_row, HRR_threshold_percent)
        output_row = input_row;
        ignitiondelay_hrr = nan;
        this_HRR = input_row.HRR{1};
        if isempty(this_HRR)
            fprintf('Empty HRR\n')
            return
        end
        this_time_HRR_sec = input_row.TimeScale{1};
        this_time_HRR_ms = this_time_HRR_sec*1000;
        [peakvalue, peakIndex] = max(this_HRR);
        
        % Step 2: Find the first index to the left of the peak where the value is smaller than 0
        leftIndex = -1; % Initialize with an invalid index
        for i = peakIndex-1:-1:1
            if this_HRR(i) < peakvalue * HRR_threshold_percent
                leftIndex = i;
                ignitiondelay_hrr = this_time_HRR_ms(leftIndex);
                break;
            end
        end
        
        if leftIndex == -1
            disp('No value less than 0 found to the left of the peak');
        end
        output_row.IgnitionDelay_HRR = {ignitiondelay_hrr};
end