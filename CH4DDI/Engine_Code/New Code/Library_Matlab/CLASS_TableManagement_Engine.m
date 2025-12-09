classdef CLASS_TableManagement_Engine
    %CLASS_TABLEMANAGEMENT_ENGINE Simplified table management without condition matrix
    %   Manages data table for H2DDI experiments without condition matching
    
    properties
        DirectoryDateSheet = ''
        VariableNameRange = 'A2'  % Row with variable names
        DataRange = 'A3'          % First row of actual data
        DataSheetName = 'Sheet1'
        
        DataMatrix  % All experimental data
    end
    
    methods
        function obj = CLASS_TableManagement_Engine(DirectoryDateSheet, varargin)
            %Constructor - reads Excel file and initializes data matrix
            %
            % Usage:
            %   obj = Class_TableManagement_Engine(excelPath)
            %   obj = Class_TableManagement_Engine(excelPath, Configs)
            %
            % Inputs:
            %   DirectoryDateSheet - Path to Excel file
            %   Configs (optional) - Configuration struct with settings:
            %       .VariableNameRange (default: 'A2')
            %       .DataRange (default: 'A4')
            %       .DataSheetName (default: 'Sheet1')
            
            obj.DirectoryDateSheet = DirectoryDateSheet;
            
            % Handle optional Configs input
            if nargin > 1 && isstruct(varargin{1})
                Configs = varargin{1};
                
                % Use settings from Configs if available
                if isfield(Configs, 'VariableNameRange')
                    obj.VariableNameRange = Configs.VariableNameRange;
                end
                if isfield(Configs, 'DataRange')
                    obj.DataRange = Configs.DataRange;
                end
                if isfield(Configs, 'DataSheetName')
                    obj.DataSheetName = Configs.DataSheetName;
                end
            end
            
            % Configure table reading options
            TableOption = detectImportOptions(obj.DirectoryDateSheet, ...
                'Sheet', obj.DataSheetName, ...
                'VariableNamingRule', 'preserve');
            % Preserve variable names that are not valid MATLAB identifiers
            TableOption.VariableNamesRange = obj.VariableNameRange;
            TableOption.DataRange = obj.DataRange;
            
            % Read data matrix
            DataMatrix = readtable(obj.DirectoryDateSheet, TableOption);
            DataMatrix.Properties.VariableDescriptions = {};
            obj.DataMatrix = DataMatrix;
            
            % Display initial status
            obj = obj.DisplayDataMatrixStatus();
        end
        
        %% Set DataMatrix
        function obj = SetDataMatrix(obj, DataMatrix_Processed)
            %SETDATAMATRIX Update the data matrix with processed data
            obj.DataMatrix = DataMatrix_Processed;
        end
        
        %% Display Status
        function obj = DisplayDataMatrixStatus(obj)
            %DISPLAYDATAMATRIXSTATUS Display summary of data matrix status
            total_runs = height(obj.DataMatrix);
            
            % Check if 'Discard' column exists
            if ismember('Discard', obj.DataMatrix.Properties.VariableNames)
                discarded_runs = sum([obj.DataMatrix.Discard] == 1 | [obj.DataMatrix.Discard] == true);
                valid_runs = total_runs - discarded_runs;
                
                fprintf(['\n' repmat('=', 1, 60) '\n']);
                fprintf('DATA MATRIX STATUS\n');
                fprintf([repmat('=', 1, 60) '\n']);
                fprintf('Total runs (ALL):           %03d\n', total_runs);
                fprintf('Discarded runs:             %03d\n', discarded_runs);
                fprintf('Valid (non-discarded) runs: %03d\n', valid_runs);
                fprintf([repmat('=', 1, 60) '\n\n']);
            else
                fprintf(['\n' repmat('=', 1, 60) '\n']);
                fprintf('DATA MATRIX STATUS\n');
                fprintf([repmat('=', 1, 60) '\n']);
                fprintf('Total runs (ALL):           %03d\n', total_runs);
                fprintf('Note: No "Discard" column found in data matrix\n');
                fprintf([repmat('=', 1, 60) '\n\n']);
            end
        end
        
        %% Discard Run
        function obj = DiscardRun(obj, array_runNumber_for_discard)
            %DISCARDRUN Mark specified runs as discarded
            %   array_runNumber_for_discard: array of run numbers to discard
            
            dataMatrix = obj.DataMatrix;
            
            % Ensure Discard column exists
            if ~ismember('Discard', dataMatrix.Properties.VariableNames)
                warning('Discard column does not exist. Adding it now...');
                dataMatrix.Discard = false(height(dataMatrix), 1);
            end
            
            for this_runNumber_for_discard = array_runNumber_for_discard
                % Validate run number
                if this_runNumber_for_discard < 1 || this_runNumber_for_discard > height(dataMatrix)
                    warning('Run number %d is out of range (1-%d)', ...
                        this_runNumber_for_discard, height(dataMatrix));
                    continue;
                end
                
                % Mark as discarded
                dataMatrix(this_runNumber_for_discard, 'Discard') = {true};
                fprintf('Run %04d marked as DISCARDED\n', this_runNumber_for_discard);
            end
            
            obj.DataMatrix = dataMatrix;
            obj = obj.DisplayDataMatrixStatus();
        end
        
        %% Recover Run
        function obj = RecoverRun(obj, array_runNumber_for_recover)
            %RECOVERRUN Mark specified runs as non-discarded (recover them)
            %   array_runNumber_for_recover: array of run numbers to recover
            
            dataMatrix = obj.DataMatrix;
            
            % Ensure Discard column exists
            if ~ismember('Discard', dataMatrix.Properties.VariableNames)
                warning('Discard column does not exist. Cannot recover runs.');
                return;
            end
            
            for this_runNumber_for_recover = array_runNumber_for_recover
                % Validate run number
                if this_runNumber_for_recover < 1 || this_runNumber_for_recover > height(dataMatrix)
                    warning('Run number %d is out of range (1-%d)', ...
                        this_runNumber_for_recover, height(dataMatrix));
                    continue;
                end
                
                % Mark as non-discarded
                dataMatrix(this_runNumber_for_recover, 'Discard') = {false};
                fprintf('Run %04d marked as RECOVERED (non-discarded)\n', this_runNumber_for_recover);
            end
            
            obj.DataMatrix = dataMatrix;
            obj = obj.DisplayDataMatrixStatus();
        end
        
        %% Get Valid Runs
        function validDataMatrix = GetValidRuns(obj)
            %GETVALIDRUNS Return only non-discarded runs
            if ismember('Discard', obj.DataMatrix.Properties.VariableNames)
                validDataMatrix = obj.DataMatrix([obj.DataMatrix.Discard] == 0 | ...
                    [obj.DataMatrix.Discard] == false, :);
            else
                % If no Discard column, return all runs
                validDataMatrix = obj.DataMatrix;
            end
        end
        
        %% Get Discarded Runs
        function discardedDataMatrix = GetDiscardedRuns(obj)
            %GETDISCARDEDRUNS Return only discarded runs
            if ismember('Discard', obj.DataMatrix.Properties.VariableNames)
                discardedDataMatrix = obj.DataMatrix([obj.DataMatrix.Discard] == 1 | ...
                    [obj.DataMatrix.Discard] == true, :);
            else
                % If no Discard column, return empty table
                discardedDataMatrix = obj.DataMatrix([], :);
            end
        end
    end
    
    %% Static Methods
    methods (Static)
        function InputTable = ExtendTableWithChar(InputTable, columnNames)
            %EXTENDTABLEWITHCHAR Add new columns with empty char arrays
            %   columnNames: cell array of column names to add
            
            for ii = 1:length(columnNames)
                % Initialize as empty string cells
                InputTable.(columnNames{ii}) = repmat({''}, height(InputTable), 1);
            end
        end
        
        function InputTable = ExtendTableWithDouble(InputTable, columnNames)
            %EXTENDTABLEWITHDOUBLE Add new columns with empty cells for doubles
            %   columnNames: cell array of column names to add
            
            for ii = 1:length(columnNames)
                InputTable(:, columnNames{ii}) = cell(height(InputTable), 1);
            end
        end
        
        function InputTable = ExtendTableWithLogical(InputTable, columnNames, defaultValue)
            %EXTENDTABLEWITHLOGICAL Add new columns with logical values
            %   columnNames: cell array of column names to add
            %   defaultValue: default logical value (default: false)
            
            if nargin < 3
                defaultValue = false;
            end
            
            for ii = 1:length(columnNames)
                InputTable.(columnNames{ii}) = repmat(defaultValue, height(InputTable), 1);
            end
        end
        
        function InputTable = LocateRawFileDirectory(InputTable, TopFolderDirectory)
            %LOCATERAWFILEDIRECTORY Find and add raw file directories to table
            %   Locates .tdms files based on Date and Set number
            %   File naming pattern: _Set[XX]_[conditions]_Take[YY].tdms
            %   Stores all Takes for each Set as a cell array
            
            if isempty(InputTable)
                error("Data matrix is empty")
            end
            
            % Determine which column to use for set identification
            if ismember('Set', InputTable.Properties.VariableNames)
                setColumnName = 'Set';
            elseif ismember('Number', InputTable.Properties.VariableNames)
                setColumnName = 'Number';
            elseif ismember('Run', InputTable.Properties.VariableNames)
                setColumnName = 'Run';
            else
                error('Cannot find set/run number column (Set/Number/Run)');
            end
            
            % Check for Date column
            if ~ismember('Date', InputTable.Properties.VariableNames)
                error('Cannot find Date column in table');
            end
            
            fprintf('Locating raw .tdms files...\n');
            totalFilesFound = 0;
            
            for ii = 1:height(InputTable)
                % Get date and set number
                RunDate = InputTable{ii, 'Date'};
                SetNumber = InputTable{ii, setColumnName};
                
                % Convert to numeric if needed
                if iscell(RunDate)
                    RunDate = RunDate{1};
                end
                if iscell(SetNumber)
                    SetNumber = SetNumber{1};
                end
                
                % Validate inputs
                if ~isnumeric(RunDate)
                    warning('Row %03d: Date is not numeric (Date: %s)', ii, string(RunDate));
                    continue
                end
                if ~isnumeric(SetNumber)
                    warning('Row %03d: Set number is not numeric (Set: %s)', ii, string(SetNumber));
                    continue
                end
                
                % Get date folder
                FolderPattern = fullfile(TopFolderDirectory, num2str(RunDate));
                FolderByDate = dir(FolderPattern);
                
                if isempty(FolderByDate)
                    % Try with wildcard
                    FolderPattern = fullfile(TopFolderDirectory, strcat(num2str(RunDate), '*'));
                    FolderByDate = dir(FolderPattern);
                end
                
                if isempty(FolderByDate)
                    warning('Row %03d: Folder for date %d not found', ii, RunDate);
                    continue
                elseif length(FolderByDate) > 1
                    % Multiple folders found, use first one
                    warning('Row %03d: Multiple folders found for date %d, using first', ii, RunDate);
                end
                
                FileDirByDate = fullfile(FolderByDate(1).folder, FolderByDate(1).name);
                
                % Search for TDMS files in Pressure subfolder
                % Your structure: [Date]_H2DDI\[Date]_Pressure\_Set[XX]_*_Take*.tdms
                PressureFolder = dir(fullfile(FileDirByDate, '*Pressure'));
                
                if isempty(PressureFolder)
                    % Try searching directly in date folder if no Pressure subfolder
                    searchPath = FileDirByDate;
                    fprintf('  Row %03d: No Pressure subfolder found, searching in date folder\n', ii);
                else
                    % Use Pressure subfolder
                    searchPath = fullfile(PressureFolder(1).folder, PressureFolder(1).name);
                end
                
                % Search for TDMS files matching: _Set[XX]_*_Take*.tdms
                FilePattern = sprintf('_Set%02d_*_Take*.tdms', SetNumber);
                S_TdmsFiles = dir(fullfile(searchPath, FilePattern));
                
                if ~isempty(S_TdmsFiles)
                    % Sort files by Take number to ensure consistent ordering
                    fileNames = {S_TdmsFiles.name};
                    
                    % Extract Take numbers for sorting
                    takeNumbers = zeros(length(fileNames), 1);
                    for jj = 1:length(fileNames)
                        % Extract Take number using regexp
                        tokens = regexp(fileNames{jj}, '_Take(\d+)\.tdms', 'tokens');
                        if ~isempty(tokens)
                            takeNumbers(jj) = str2double(tokens{1}{1});
                        end
                    end
                    
                    % Sort by Take number
                    [~, sortIdx] = sort(takeNumbers);
                    S_TdmsFiles = S_TdmsFiles(sortIdx);
                    
                    % Create cell array with full paths to all Takes
                    allTakePaths = cell(length(S_TdmsFiles), 1);
                    for jj = 1:length(S_TdmsFiles)
                        allTakePaths{jj} = fullfile(S_TdmsFiles(jj).folder, S_TdmsFiles(jj).name);
                    end
                    
                    % Store in table as nested cell array
                    InputTable{ii, 'File_TDMS_AllTakes'} = {allTakePaths};
                    
                    totalFilesFound = totalFilesFound + length(S_TdmsFiles);
                    fprintf('  Row %03d (Set %02d, Date %d): Found %d Take(s)\n', ...
                        ii, SetNumber, RunDate, length(S_TdmsFiles));
                else
                    warning('Row %03d: No TDMS files found for Set %02d on date %d', ...
                        ii, SetNumber, RunDate);
                end
            end
            
            fprintf('Raw file location completed: %d total files found\n', totalFilesFound);
        end
        
        function summary = GetRunSummary(InputTable)
            %GETRUNSUMMARY Generate a summary of the data table
            %   Returns a struct with summary statistics
            
            summary = struct();
            summary.TotalRuns = height(InputTable);
            
            if ismember('Discard', InputTable.Properties.VariableNames)
                summary.DiscardedRuns = sum([InputTable.Discard] == 1 | [InputTable.Discard] == true);
                summary.ValidRuns = summary.TotalRuns - summary.DiscardedRuns;
            else
                summary.DiscardedRuns = 0;
                summary.ValidRuns = summary.TotalRuns;
            end
            
            % Add unique values for key columns if they exist
            if ismember('Date', InputTable.Properties.VariableNames)
                dates = InputTable.Date;
                if iscell(dates)
                    dates = cell2mat(dates);
                end
                summary.UniqueDates = unique(dates(~isnan(dates)));
                summary.NumberOfTestDays = length(summary.UniqueDates);
            end
            
            if ismember('Type', InputTable.Properties.VariableNames)
                summary.TestTypes = unique(InputTable.Type);
            end
        end
    end
end