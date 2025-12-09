% Read datasheet .xlsx to construct a datatable
% Check amount of total runs/non-discard runs
% Get directories for raw files (image/signal)
% Process pressure signal 
% Store data from -10 to 20 ms at aSOI
% Master/ Primary Branches 
%%%%% add something random
%% Read table
clearvars; close all; clc;
% addpath 'Library_Matlab\'
DataSheetDirectory   = 'E:\20240306_H2_DoubleInjection.xlsx';
obj = Class_TableManagement(DataSheetDirectory);
    obj.ConditionMatchColumn = 5:16;
T_DataAll = obj.GetDataMatrix();
T_DataAll = Class_TableManagement.ExtendTableWithDouble(T_DataAll, {'TimeScale', 'P_Raw','P_Corrected', 'HRR', 'IgnitionDelay', 'P_TrigTrue', 'PD_TimeScale', 'PD_Signal'}); % add cols (double) for store data 
    % Consider have something like "SO"
        % Alternative inheritence
        
%% Track raw files directories and add into table:
T_DataAll = Class_TableManagement.ExtendTableWithChar(T_DataAll, {'File_Image', 'File_Pressure','File_Log', 'File_Photodiode'}); % add cols (char) for store file directory
T_DataAll = Class_TableManagement.LocateRawFileDirectory(T_DataAll, 'F:\H2_Double');  % Change disk letter here

%% Process pressure signal (.tdms)
% Process all data (discarded run will be skipped by using try-catch) 
% Pressure data processing for selected rows
IndexStart = 1;
IndexEnd = 666;
T_DataAll = FUNCTION_PressureProcess(T_DataAll, IndexStart, IndexEnd,0);
%T_DataAll = FUNCTION_PhotodiodeProcess(T_DataAll, IndexStart, IndexEnd,1);
obj = obj.SetDataMatrix(T_DataAll);
save("Outputs\OBJ_Table","obj")

%% Internal functions
% Process pressure data from Table
function T_DataAll = FUNCTION_PressureProcess_Table(T_DataAll, IndexStart, IndexEnd, plotCheck)
    for ii = IndexStart:IndexEnd
        try
            fprintf('# Processing index %03d ...\n', ii);
            thisRow = T_DataAll(ii,:);
    
            Obj_Pressure = CLASS_PressureProcess();

            % Injection Delay
            Obj_Pressure.InjectionDelay = thisRow.Inj_delay;

            % Toggle checking plot 
                Obj_Pressure.bPlotIgnitionDelay = plotCheck;
                Obj_Pressure.bPlotPressureCorrected = plotCheck;

            Obj_Pressure = Obj_Pressure.FUNCTION_PressureProcess(thisRow);
            PressureInfo = Obj_Pressure.GetPressureInfo;
    
            % Check if the table ambient conditions are the same as the ones in log file
            if ~isempty(PressureInfo.TriggerPressureSet) % This means .log reading above failed 
                try
                    assert(PressureInfo.TriggerPressureSet == thisRow.('P_trig'));
                    assert(PressureInfo.O2Percentage == thisRow.('O2'));
                catch
                    warning('Condition assertion failed')
                end
            end
            
            T_DataAll(ii, 'TimeScale') = {{PressureInfo.TimeScale}};
            T_DataAll(ii, 'P_Raw') = {{PressureInfo.Pressure}};
            T_DataAll(ii, 'P_Corrected') = {{PressureInfo.P_Corrected}};
            T_DataAll(ii, 'HRR') = {{PressureInfo.HRR}};
            T_DataAll(ii, 'IgnitionDelay') = {{PressureInfo.IgnitionDelay}};
            T_DataAll(ii, 'P_TrigTrue') = {{PressureInfo.TriggerPressureTrue}};
        catch
            warning('# Error for the run with index %03d ...', ii);
        end
    end
end

%%
function T_DataAll = FUNCTION_PhotodiodeProcess_Table(T_DataAll, IndexStart, IndexEnd, plotCheck)
    for ii = IndexStart:IndexEnd
        try
            fprintf('# Processing index %03d ...\n', ii);
            thisRow = T_DataAll(ii,:);
            Gain_PD = thisRow.Gain(1);
            Obj_Photodiode = CLASS_PhotodiodeProcess();
            Obj_Photodiode.Gain_PD = Gain_PD;
            % Toggle checking plot 
                Obj_Photodiode.bPlotIgnitionDelay = plotCheck;
                Obj_Photodiode.bPlotPressureCorrected = plotCheck;

            Obj_Photodiode = Obj_Photodiode.FUNCTION_PhotodiodeProcess(thisRow);
            PhotodiodeInfo = Obj_Photodiode.GetPhotodiodeInfo;
            
            T_DataAll(ii, 'PD_TimeScale') = {{PhotodiodeInfo.TimeScale}};
            T_DataAll(ii, 'PD_Signal') = {{PhotodiodeInfo.Photodiode}};
        catch
            warning('# Error for the run with index %03d ...', ii);
        end
    end
end
