% Read datasheet .xlsx to construct a datatable
% Check amount of total runs/non-discard runs
% Get directories for raw files (image/signal)
% Process pressure signal 
% Store data from -10 to 20 ms at aSOI
%% Read table
clearvars; close all; clc;
Configs = Configuration();

DataSheetDirectory   = Configs.DataSheetDirectory;
ConditionMatchColumn =  Configs.ConditionMatchColumn;
DataTopfolder = Configs.DataTopfolder;
   
% Construct table 
obj = Class_TableManagement(DataSheetDirectory, ConditionMatchColumn);                                                                                             
T_DataAll = obj.DataMatrix;
T_ConditionMatrix = obj.ConditionMatrix;
% Assert if the fieldnames are the same
if ~isequal(fieldnames(obj.ConditionMatrix_ROI_Dependent), fieldnames(obj.DataMatrix_ROI_Dependent))
    fprintf('Excel Sheet1 and Sheet2 have different names')
    return
end

% Extend table to include expected post-processing results for each single run (i.e., each row)
T_DataAll = Class_TableManagement.ExtendTableWithDouble(T_DataAll, Configs.FieldNameExtended); 

% Track raw files directories and add into table:
T_DataAll = Class_TableManagement.LocateRawFileDirectory(T_DataAll, DataTopfolder);  % Change disk letter here
%
obj = obj.SetDataMatrix(T_DataAll); % Put modified Table to replace the one in object
%% Save
keyboard
DataObjSaveDir = Configs.DataObjSaveDir;     
save(DataObjSaveDir,"obj")
fprintf("Saving finished\n")