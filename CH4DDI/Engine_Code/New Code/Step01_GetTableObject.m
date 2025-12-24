%% STEP 1: Read Table and Locate Raw Files
% Read datasheet .xlsx to construct a datatable
% Check amount of total runs/non-discard runs
% Get directories for raw files (TDMS)
% 
% This script:
% 1. Reads Excel file and creates data matrix
% 2. Extends table with columns for post-processing
% 3. Locates all raw TDMS files (multiple Takes per Set)
% 4. Saves table object for later processing
%
% test vscode git push
% test vscode workspace
%% Read table
clearvars; close all; clc;

% Load configuration
Configs = Configuration();
DataSheetDirectory = Configs.DataSheetDirectory;
DataTopfolder = Configs.DataTopfolder;

%% Construct table
fprintf('\n');
fprintf(repmat('=', 1, 70));
fprintf('\nSTEP 1: READING DATA TABLE\n');
fprintf(repmat('=', 1, 70));
fprintf('\n\n');

fprintf('Reading Excel file: %s\n', DataSheetDirectory);
obj = CLASS_TableManagement_Engine(DataSheetDirectory, Configs);
T_DataAll = obj.DataMatrix;

fprintf('Data matrix loaded:\n');
fprintf('  - Total runs: %d\n', height(T_DataAll));
fprintf('  - Total columns: %d\n', width(T_DataAll));

%% Extend table to include expected post-processing results for each single run (i.e., each row)
fprintf('\nExtending table with processing columns...\n');
T_DataAll = CLASS_TableManagement_Engine.ExtendTableWithDouble(T_DataAll, Configs.FieldNameExtended);

fprintf('Extended columns added:\n');
for ii = 1:length(Configs.FieldNameExtended)
    fprintf('  - %s\n', Configs.FieldNameExtended{ii});
end

%% Track raw files directories and add into table
fprintf('\nLocating raw file directories...\n');
T_DataAll = CLASS_TableManagement_Engine.LocateRawFileDirectory(T_DataAll, DataTopfolder);

% Update object with modified table
obj = obj.SetDataMatrix(T_DataAll);

%% Display summary
fprintf('\n');
fprintf(repmat('=', 1, 70));
fprintf('\nSUMMARY\n');
fprintf(repmat('=', 1, 70));
fprintf('\n');

% Count files found
nSetsWithFiles = 0;
totalTakes = 0;
for ii = 1:height(T_DataAll)
    if ismember('File_TDMS_AllTakes', T_DataAll.Properties.VariableNames)
        if ~isempty(T_DataAll.File_TDMS_AllTakes{ii})
            nSetsWithFiles = nSetsWithFiles + 1;
            totalTakes = totalTakes + length(T_DataAll.File_TDMS_AllTakes{ii});
        end
    end
end

fprintf('Sets with TDMS files found: %d / %d\n', nSetsWithFiles, height(obj.GetValidRuns()));
fprintf('Total Take files found: %d\n', totalTakes);
fprintf(repmat('=', 1, 70));
fprintf('\n\n');

%% Save
keyboard
fprintf('Saving object...\n');
DataObjSaveDir = Configs.DataObjSaveDir;     
save(DataObjSaveDir,"obj")
DataObjSaveDir = Configs.DataObjSaveDir;
fprintf('Object saved to: %s\n',DataObjSaveDir);
fprintf('\n✓ STEP 1 COMPLETE\n\n');

fprintf('Next step: Process TDMS files (Step 2)\n');
