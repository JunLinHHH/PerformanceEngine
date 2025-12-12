% STEP 2 

clearvars; close all; clc; 
Configs = Configuration();
Loaded = load(Configs.DataObjReadDir);

obj = Loaded.obj;
T_DataAll = obj.DataMatrix;

set_selected = 1:height(T_DataAll);

T_DataAll = Table_PressureProcessEngine(T_DataAll, set_selected);
obj = obj.SetDataMatrix(T_DataAll); % Put modified Table to replace the one in object


fprintf("Process finished\n")
%%
function T_DataAll = Table_PressureProcessEngine(T_DataAll, set_selected)
    fprintf('\nProcessing %d Sets...\n', length(set_selected));
    
    for ii = set_selected
        try
            rowIdx = find(T_DataAll.Set == ii, 1);
            this_row = T_DataAll(rowIdx, :);
            output_row = this_row;
            fprintf('Set %02d: ', ii);
            
            tdmsFiles = this_row.File_TDMS_AllTakes{1};
            
            % Create processor
            Obj_Pressure = CLASS_PressureProcessEngine(tdmsFiles, ii);
            Obj_Pressure = Obj_Pressure.ReadPressureFile_TDMS();
            PressureInfo = Obj_Pressure.GetPressureInfo();
            
            % Store directly (like your example)
            output_row.P_Raw = {PressureInfo};
            % Later add more as you process:
            % output_row.TimeScale = {PressureInfo.TimeScale};
            % output_row.P_Corrected = {PressureInfo.P_Corrected};
            % output_row.HRR = {PressureInfo.HRR};
            
            T_DataAll(rowIdx, :) = output_row;
            
            fprintf('✓ %d Takes\n', length(PressureInfo));
            
        catch ME
            fprintf('✗ Error: %s\n', ME.message);
        end
    end
    
    fprintf('Done\n\n');
end