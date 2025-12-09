%% Test Script - Verify Simple Table Management Setup
% This script tests the basic functionality of Class_SimpleTableManagement
% Run this to verify your installation is working correctly

clearvars; close all; clc;

fprintf('\n');
fprintf(repmat('=', 1, 70));
fprintf('\n');
fprintf('TEST SCRIPT: Class_TableManagement_Engine\n');
fprintf(repmat('=', 1, 70));
fprintf('\n\n');

testResults = struct();

%% Test 1: Check if Configuration exists
fprintf('TEST 1: Checking Configuration function...\n');
try
    which Configuration
    Configs = Configuration();
    fprintf('  ✓ Configuration function found and working\n');
    testResults.Config = 'PASS';
catch ME
    fprintf('  ✗ FAILED: %s\n', ME.message);
    fprintf('  → Make sure Configuration.m is in your MATLAB path\n');
    fprintf('  → Current folder: %s\n', pwd);
    testResults.Config = 'FAIL';
    fprintf('\n=== TEST ABORTED ===\n');
    fprintf('Please make sure Configuration.m is accessible and try again.\n');
    return;
end

%% Test 2: Check Excel File Exists
fprintf('\nTEST 2: Checking Excel File...\n');
try
    DataSheetDirectory = Configs.DataSheetDirectory;
    if isfile(DataSheetDirectory)
        fprintf('  ✓ Excel file found: %s\n', DataSheetDirectory);
        testResults.ExcelFile = 'PASS';
    else
        fprintf('  ✗ FAILED: Excel file not found\n');
        fprintf('  Path: %s\n', DataSheetDirectory);
        fprintf('  → Update Configs.DataSheetDirectory in Configuration.m\n');
        testResults.ExcelFile = 'FAIL';
        fprintf('\n=== TEST ABORTED ===\n');
        return;
    end
catch ME
    fprintf('  ✗ FAILED: %s\n', ME.message);
    testResults.ExcelFile = 'FAIL';
    return;
end

%% Test 3: Check if Class exists
fprintf('\nTEST 3: Checking Class_TableManagement_Engine...\n');
try
    which Class_TableManagement_Engine
    fprintf('  ✓ Class_TableManagement_Engine found\n');
    testResults.ClassExists = 'PASS';
catch ME
    fprintf('  ✗ FAILED: Class_TableManagement_Engine not found\n');
    fprintf('  → Make sure Class_TableManagement_Engine.m is in your MATLAB path\n');
    fprintf('  → Current folder: %s\n', pwd);
    testResults.ClassExists = 'FAIL';
    return;
end

%% Test 4: Create Object and Read Data
fprintf('\nTEST 4: Creating Object and Reading Excel Data...\n');
try
    obj = Class_TableManagement_Engine(DataSheetDirectory);
    fprintf('  ✓ Object created successfully\n');
    fprintf('  ✓ Data matrix loaded: %d rows × %d columns\n', ...
        height(obj.DataMatrix), width(obj.DataMatrix));
    testResults.ObjectCreation = 'PASS';
catch ME
    fprintf('  ✗ FAILED: %s\n', ME.message);
    fprintf('  → Check if Excel file format is correct\n');
    testResults.ObjectCreation = 'FAIL';
    return;
end

%% Test 5: Check Required Columns
fprintf('\nTEST 5: Checking Excel Columns...\n');
requiredCols = {'Date', 'Set'};
allPass = true;

for ii = 1:length(requiredCols)
    if ismember(requiredCols{ii}, obj.DataMatrix.Properties.VariableNames)
        fprintf('  ✓ Found column: %s\n', requiredCols{ii});
    else
        fprintf('  ✗ Missing column: %s\n', requiredCols{ii});
        allPass = false;
    end
end

if allPass
    testResults.RequiredColumns = 'PASS';
else
    testResults.RequiredColumns = 'FAIL';
    fprintf('  → Check Excel file row 2 has column names: Set, Date, etc.\n');
end

%% Test 6: Check FieldNameExtended Configuration
fprintf('\nTEST 6: Checking Configuration FieldNameExtended...\n');
try
    if isfield(Configs, 'FieldNameExtended')
        fprintf('  ✓ FieldNameExtended exists\n');
        
        % Check if it has File_TDMS_AllTakes
        if any(strcmp(Configs.FieldNameExtended, 'File_TDMS_AllTakes'))
            fprintf('  ✓ Contains ''File_TDMS_AllTakes'' (correct for multiple Takes)\n');
            testResults.FieldNames = 'PASS';
        else
            fprintf('  ⚠ WARNING: ''File_TDMS_AllTakes'' not found\n');
            fprintf('  Current fields: %s\n', strjoin(Configs.FieldNameExtended, ', '));
            fprintf('  → Update Configuration.m: Change ''File_Pressure'' to ''File_TDMS_AllTakes''\n');
            testResults.FieldNames = 'WARN';
        end
    else
        fprintf('  ✗ FieldNameExtended not found in Configuration\n');
        testResults.FieldNames = 'FAIL';
    end
catch ME
    fprintf('  ✗ FAILED: %s\n', ME.message);
    testResults.FieldNames = 'FAIL';
end

%% Test 7: Extend Table
fprintf('\nTEST 7: Testing Table Extension...\n');
try
    T = obj.DataMatrix;
    testCols = {'TestCol1', 'TestCol2'};
    T = Class_TableManagement_Engine.ExtendTableWithDouble(T, testCols);
    
    if ismember('TestCol1', T.Properties.VariableNames)
        fprintf('  ✓ Table extension works\n');
        testResults.TableExtension = 'PASS';
    else
        fprintf('  ✗ Column not added properly\n');
        testResults.TableExtension = 'FAIL';
    end
catch ME
    fprintf('  ✗ FAILED: %s\n', ME.message);
    testResults.TableExtension = 'FAIL';
end

%% Test 8: Get Valid/Discarded Runs
fprintf('\nTEST 8: Testing Get Valid/Discarded Runs...\n');
try
    validRuns = obj.GetValidRuns();
    discardedRuns = obj.GetDiscardedRuns();
    
    fprintf('  ✓ Valid runs: %d\n', height(validRuns));
    fprintf('  ✓ Discarded runs: %d\n', height(discardedRuns));
    testResults.GetRuns = 'PASS';
catch ME
    fprintf('  ✗ FAILED: %s\n', ME.message);
    testResults.GetRuns = 'FAIL';
end

%% Test 9: File Location (if raw data folder exists)
fprintf('\nTEST 9: Testing Raw File Location...\n');
try
    DataTopfolder = Configs.DataTopfolder;
    
    if isfolder(DataTopfolder)
        fprintf('  ✓ Data folder found: %s\n', DataTopfolder);
        
        % Try to locate files
        T = obj.DataMatrix;
        T = Class_TableManagement_Engine.LocateRawFileDirectory(T, DataTopfolder);
        
        % Check if File_TDMS_AllTakes column was created
        if ismember('File_TDMS_AllTakes', T.Properties.VariableNames)
            nSetsWithFiles = 0;
            totalTakes = 0;
            
            for ii = 1:height(T)
                if ~isempty(T.File_TDMS_AllTakes{ii})
                    nSetsWithFiles = nSetsWithFiles + 1;
                    totalTakes = totalTakes + length(T.File_TDMS_AllTakes{ii});
                end
            end
            
            if nSetsWithFiles > 0
                fprintf('  ✓ Found TDMS files for %d Set(s)\n', nSetsWithFiles);
                fprintf('  ✓ Total Take files found: %d\n', totalTakes);
                testResults.FileLocation = 'PASS';
            else
                fprintf('  ⚠ No files found\n');
                fprintf('  → Check folder structure: %s\\[Date]\\*.tdms\n', DataTopfolder);
                testResults.FileLocation = 'WARN';
            end
        else
            fprintf('  ✗ File_TDMS_AllTakes column not created\n');
            testResults.FileLocation = 'FAIL';
        end
    else
        fprintf('  ⊘ SKIP: Data folder not found: %s\n', DataTopfolder);
        fprintf('  → Update Configs.DataTopfolder in Configuration.m\n');
        testResults.FileLocation = 'SKIP';
    end
catch ME
    fprintf('  ✗ FAILED: %s\n', ME.message);
    testResults.FileLocation = 'FAIL';
end

%% Test 10: Discard and Recover (if data exists)
fprintf('\nTEST 10: Testing Discard and Recover...\n');
try
    if height(obj.DataMatrix) >= 1
        % Save original
        originalDiscard = obj.DataMatrix.Discard(1);
        
        % Test discard
        obj = obj.DiscardRun(1);
        if obj.DataMatrix.Discard(1) == true || obj.DataMatrix.Discard(1) == 1
            fprintf('  ✓ Discard function works\n');
            
            % Test recover
            obj = obj.RecoverRun(1);
            if obj.DataMatrix.Discard(1) == false || obj.DataMatrix.Discard(1) == 0
                fprintf('  ✓ Recover function works\n');
                testResults.DiscardRecover = 'PASS';
            else
                fprintf('  ✗ Recover failed\n');
                testResults.DiscardRecover = 'FAIL';
            end
        else
            fprintf('  ✗ Discard failed\n');
            testResults.DiscardRecover = 'FAIL';
        end
    else
        fprintf('  ⊘ SKIP: No data to test\n');
        testResults.DiscardRecover = 'SKIP';
    end
catch ME
    fprintf('  ✗ FAILED: %s\n', ME.message);
    testResults.DiscardRecover = 'FAIL';
end

%% Final Results
fprintf('\n');
fprintf(repmat('=', 1, 70));
fprintf('\n');
fprintf('TEST RESULTS SUMMARY\n');
fprintf(repmat('=', 1, 70));
fprintf('\n');

testNames = fieldnames(testResults);
passCount = 0;
failCount = 0;
warnCount = 0;
skipCount = 0;

for ii = 1:length(testNames)
    result = testResults.(testNames{ii});
    fprintf('%-25s: %s\n', testNames{ii}, result);
    
    if strcmp(result, 'PASS')
        passCount = passCount + 1;
    elseif strcmp(result, 'FAIL')
        failCount = failCount + 1;
    elseif strcmp(result, 'WARN')
        warnCount = warnCount + 1;
    elseif strcmp(result, 'SKIP')
        skipCount = skipCount + 1;
    end
end

fprintf(repmat('=', 1, 70));
fprintf('\n');
fprintf('Total: %d tests\n', length(testNames));
fprintf('Passed: %d | Failed: %d | Warnings: %d | Skipped: %d\n', ...
    passCount, failCount, warnCount, skipCount);
fprintf(repmat('=', 1, 70));
fprintf('\n');

%% Overall Result
if failCount == 0 && warnCount == 0
    fprintf('\n✓✓✓ ALL TESTS PASSED! ✓✓✓\n');
    fprintf('Your setup is working correctly.\n');
    fprintf('Next step: Run SCRIPT_ReadDataTable_Simple.m\n\n');
elseif failCount == 0 && warnCount > 0
    fprintf('\n⚠⚠⚠ TESTS PASSED WITH WARNINGS ⚠⚠⚠\n');
    fprintf('Basic functionality works but check warnings above.\n');
    fprintf('You can proceed but some features may not work optimally.\n\n');
else
    fprintf('\n✗✗✗ SOME TESTS FAILED ✗✗✗\n');
    fprintf('Please fix the failed tests before proceeding.\n\n');
    fprintf('Common fixes:\n');
    fprintf('  1. Update Configuration.m:\n');
    fprintf('     - Change ''File_Pressure'' to ''File_TDMS_AllTakes''\n');
    fprintf('  2. Check file paths are correct\n');
    fprintf('  3. Ensure Excel has ''Set'' and ''Date'' columns\n\n');
end

fprintf('Documentation:\n');
fprintf('  - QUICK_FIX.md - How to update Configuration.m\n');
fprintf('  - README.md - Full documentation\n');
fprintf('  - QUICK_REFERENCE_Takes.md - Working with multiple Takes\n\n');

