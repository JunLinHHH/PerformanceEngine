classdef Class_TableManagement
    %CLASS_TABLEMANAGEMENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        DirectoryDateSheet = ''
        VariableNameRange =  'A2'
        DataRange = 'A3'
        DataSheetName = 'Sheet1'
        TestConditionSheetName = 'Sheet2'
        ConditionMatchColumn

        DataMatrix % All data
        ConditionMatrix
        DataMatrix_runsNotMeetConditionAfterDiscard
    end

    properties (Dependent)
        % Dependent on ConditionMatchColumn
        DataMatrix_ROI_Dependent
        ConditionMatrix_ROI_Dependent
    end

    methods
        function obj = Class_TableManagement(DirectoryDateSheet, ConditionMatchColumn)
            obj.DirectoryDateSheet = DirectoryDateSheet;

            TableOption = detectImportOptions(obj.DirectoryDateSheet, 'Sheet', obj.DataSheetName,'VariableNamingRule','preserve');
                    % Preserve variable names that are not valid MATLAB identifiers such as variable names that include spaces and non-ASCII characters.
                    TableOption.VariableNamesRange = obj.VariableNameRange;
                    TableOption.DataRange = obj.DataRange;

            DataMatrix = readtable(obj.DirectoryDateSheet,TableOption);
            DataMatrix.Properties.VariableDescriptions = {};
            obj.DataMatrix = DataMatrix;

            ConditionMatrix = readtable(obj.DirectoryDateSheet,'Sheet',obj.TestConditionSheetName,'VariableNamingRule','preserve');
            obj.ConditionMatrix = ConditionMatrix;
            obj.ConditionMatchColumn = ConditionMatchColumn;
            obj = obj.DataMatrixStatus; % added in the consturctor
        end
        %%
        % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Get
        function DataMatrix_ROI_Dependent = get.DataMatrix_ROI_Dependent(obj)
            DataMatrix_ROI_Dependent = obj.DataMatrix(:, obj.ConditionMatchColumn);
        end
        function ConditionMatrix_ROI_Dependent = get.ConditionMatrix_ROI_Dependent(obj)
            ConditionMatrix_ROI_Dependent = obj.ConditionMatrix(:, obj.ConditionMatchColumn);
        end        
        % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Get        
        function obj = SetDataMatrix(obj, DataMatrix_Processed)
            obj.DataMatrix = DataMatrix_Processed;
        end
        %%
        function obj = DataMatrixStatus(obj)
            % Crop ROI - ConditionMatchColumn
            ConditionMatrix_ROI = obj.ConditionMatrix(:, obj.ConditionMatchColumn); % Change depending on table style

            DataMatrix_ROI = obj.DataMatrix([obj.DataMatrix.Discard] == 0, :);    % Get non-discarded runs
            
            
            obj.ConditionMatrix.RunNumbers = num2cell(obj.ConditionMatrix.RunNumbers); % Type conversion (original is nan, double - one element)
            % Get run amounts (after discarding, i.e., valid) for each test condition
            % Get run amounts (after discarding, i.e., valid) for runs not
            % matching condition matrix
           
            runsNotMeetConditionAfterDiscard = ones(height(DataMatrix_ROI),1); % Flag all as not matching, change to 0 if matching during loop

            for ii = 1:height(ConditionMatrix_ROI)
                Amount = 0;
                RunNumbers = [];
                for jj = 1:height(DataMatrix_ROI)
                    this_DataMatrix = DataMatrix_ROI(jj,:);
                    this_DataMatrix_ROI = this_DataMatrix(:, obj.ConditionMatchColumn);
                    this_run_number = this_DataMatrix.Number;
                    if isequal(this_DataMatrix_ROI , ConditionMatrix_ROI(ii,:))
                        Amount = Amount + 1;
                        runsNotMeetConditionAfterDiscard(jj) = 0; % Switch flag if matching
                        RunNumbers(end+1) = this_run_number;
                    end
                end
                obj.ConditionMatrix(ii,'Amount') = {Amount}; % Amount is after discarding, i.e., valid runs
                obj.ConditionMatrix.RunNumbers{ii} = RunNumbers;
            end
            
            % Runs not matching any condition will remain as 1
            obj.DataMatrix_runsNotMeetConditionAfterDiscard = DataMatrix_ROI(runsNotMeetConditionAfterDiscard == 1, :);

            fprintf(['%03d total runs (ALL) \n' ...
                '%03d discarded runs among ALL \n' ...
                '%03d non-discarded runs among ALL \n' ...
                'After discarding: \n:'...
                '\t%03d non-macthing condition-matrix runs \n'],...
                height(obj.DataMatrix), ...
                height(obj.DataMatrix([obj.DataMatrix.Discard] == 1, :)), ...
                height(DataMatrix_ROI), ...
                sum(runsNotMeetConditionAfterDiscard))  
        end

        function obj = DiscardRun(obj, array_runNumber_for_discard)
            dataMatrix = obj.DataMatrix;
            conditionMatrix = obj.ConditionMatrix;
            ConditionMatrix_ROI = conditionMatrix(:, obj.ConditionMatchColumn); % Change depending on table style
            
            for this_runNumber_for_discard = array_runNumber_for_discard
                for ii = 1:height(ConditionMatrix_ROI)
                    this_DataMatrix = dataMatrix(this_runNumber_for_discard,:);
                    this_DataMatrix_ROI = this_DataMatrix(:, obj.ConditionMatchColumn);
                    if isequal(this_DataMatrix_ROI , ConditionMatrix_ROI(ii,:))
                        dataMatrix(this_runNumber_for_discard, 'Discard') = {true};

                        original_runNumbers = conditionMatrix.RunNumbers{ii};
                        new_runNumbers = original_runNumbers(original_runNumbers ~= this_runNumber_for_discard);
                        conditionMatrix{ii,'RunNumbers'} = {new_runNumbers};
                        obj.ConditionMatrix = conditionMatrix;
                        fprintf('Change run %04d as discard in Data Matrix\n', this_runNumber_for_discard)
                        obj.DataMatrix = dataMatrix;
                        fprintf('Change run %04d as discard in Condition Matrix\n', this_runNumber_for_discard)
                    end
                end
            end
        end

        function obj = RecoverRun(obj, array_runNumber_for_recover)
            dataMatrix = obj.DataMatrix;
            conditionMatrix = obj.ConditionMatrix;
            ConditionMatrix_ROI = conditionMatrix(:, obj.ConditionMatchColumn); % Change depending on table style

            for this_runNumber_for_recover = array_runNumber_for_recover
                for ii = 1:height(ConditionMatrix_ROI)
                    this_DataMatrix = dataMatrix(this_runNumber_for_recover,:);
                    this_DataMatrix_ROI = this_DataMatrix(:, obj.ConditionMatchColumn);
                    if isequal(this_DataMatrix_ROI , ConditionMatrix_ROI(ii,:))
                        dataMatrix(this_runNumber_for_recover, 'Discard') = {false};
                        
                        original_runNumbers = conditionMatrix.RunNumbers{ii};
                        new_runNumbers = [original_runNumbers, this_runNumber_for_recover];
                        conditionMatrix{ii,'RunNumbers'} = {new_runNumbers};
                        obj.ConditionMatrix = conditionMatrix;
                        fprintf('Change run %04d as NON-discard in Data Matrix\n', this_runNumber_for_recover)
                        obj.DataMatrix = dataMatrix;
                        fprintf('Change run %04d as NON-discard in Condition Matrix\n', this_runNumber_for_recover)
                    end
                end
            end
        end    
    end
    methods (Static)
        function InputTable = ExtendTableWithChar(InputTable, columnNames)
            % Loop through each column name and add it to the table
            for ii = 1:length(columnNames)
                % Initialize the new column as a 0x0 char (empty char array)
                % and add it to the table. Note: To strictly follow the 0x0 char requirement,
                % an adjustment is needed since MATLAB tables do not directly support this data type.
                % Instead, we can initialize as an empty cell array of chars or use strings.
                InputTable.(columnNames{ii}) = repmat({''}, height(InputTable), 1); % Using empty string cells
                
                % Alternatively, if using strings is acceptable, which is more common:
                % tbl.(columnNames{i}) = strings(height(tbl), 1); % Initialize as empty strings
            end
        end

        function InputTable = ExtendTableWithDouble(InputTable, columnNames)
            % Loop through each column name and add it to the table
            for ii = 1:length(columnNames)
                InputTable(:,columnNames{ii}) = cell(height(InputTable),1);
            end
        end
        %%
        function InputTable = LocateRawFileDirectory(InputTable, TopFolderDirectory)
        % Put files directories into the table by tracing:
            % Run Date + Number
        
            % Date+Run can generate a unique number, e.g. there will be only one run called Run#1 at 20211207
            % Make sure all signal files have been renamed by using the image file names 
            % (using SCRIPT_BatchRenameTdmsFiles)
            if isempty(InputTable)
                error("Data matrix is empty")
            end

            for ii = 1:height(InputTable)
                % Get run date and number
                RunDate = table2array(InputTable(ii,'Date'));
                RunNumber = table2array(InputTable(ii,'Run'));
                if ~isnumeric(RunDate)
                    warning('Run date read from table is not a number')
                    return 
                end
                if ~isnumeric(RunNumber)
                    warning('Run number read from table is not a number')
                    return 
                end    

                % Get date folder
                FolderPattern = fullfile(TopFolderDirectory,strcat(string(RunDate),'*'));
                FolderByDate = dir(FolderPattern);
                
                if isempty(FolderByDate)
                    % Need to check if the date is correct in Excel table if get warning
                    warning('Index %03d: Folder contains >> %d << is not FOUND',ii, RunDate)
                    continue
                elseif length(FolderByDate)~=1
                    warning('Index %03d: Folder contains >> %d << is not UNIQUE',ii, RunDate)
                    continue
                end
            
                FileDirByDate = fullfile(FolderByDate.folder, FolderByDate.name);
                
                % Get Image file name in the date folder
                S_Image = dir(strcat(FileDirByDate,'\*Tif\',string(RunDate),...
                    sprintf('*_run_%03d',RunNumber),'*.tif'));
                if ~isempty(S_Image)
                    Dir_Image = fullfile(S_Image.folder, S_Image.name);
                    InputTable(ii,'File_Image') = {Dir_Image};
                end
                
                % Get Pressure .tdms file name and .log file name in the date folder
                S_PressureFile = dir(strcat(FileDirByDate,'\*Pressure\',string(RunDate),...
                    sprintf('*_run_%03d',RunNumber),'*.tdms'));
                if ~isempty(S_PressureFile)
                    Dir_Pressure = fullfile(S_PressureFile.folder, S_PressureFile.name);
                    InputTable(ii,'File_Pressure') = {Dir_Pressure};
                end
                S_LogFile = dir(strcat(FileDirByDate,'\*Pressure\',string(RunDate),...
                    sprintf('*_run_%03d',RunNumber),'*.log'));
                if ~isempty(S_LogFile)
                    Dir_Log = fullfile(S_LogFile.folder, S_LogFile.name);
                    InputTable(ii,'File_Log') = {Dir_Log};
                end
            
                % Get Photodiode signal file name in the date folder
                S_PhotodiodeFile = dir(strcat(FileDirByDate,'\*Photodiode\',string(RunDate),...
                    sprintf('*_run_%03d',RunNumber),'*.tdms'));
                if ~isempty(S_PhotodiodeFile)
                    Dir_Photodiode = fullfile(S_PhotodiodeFile.folder, S_PhotodiodeFile.name);
                    InputTable(ii,'File_Photodiode') = {Dir_Photodiode};
                end
            end
        end
        
        function ConditionTag = ConditionTagGenerator(InputData, obj, ii)
            if istable(InputData)
                InputData = table2struct(InputData);
            end
            fieldNames = fieldnames(InputData);
            ConditionTag = '';
            fieldNamesRange = obj.ConditionMatchColumn;
            for jj = fieldNamesRange %!!!!!!!!
                value = InputData(ii).(fieldNames{jj});
                if jj == fieldNamesRange(end)
                    ConditionTag = sprintf('%s%s-%s', ConditionTag, fieldNames{jj}, string(value));        
                else
                    ConditionTag = sprintf('%s%s-%s-', ConditionTag, fieldNames{jj}, string(value)); 
                end
                % Append the field name and value to the all rows string with underscores
            end   
             ConditionTag =  sprintf('%02d-%s',ii,ConditionTag);
        end
    end
end

