% 
% %% H2 LI
% function Configs = Configuration()
%     addpath 'Library_Matlab\'
%     addpath 'Library_Matlab\tdms_package\Version_2p5_Final\v2p5'
%     addpath 'Library_Matlab\tdms_package\Version_2p5_Final\v2p5\tdmsSubfunctions'
% 
%     % Excel sheet location
%     Configs.DataSheetDirectory   = 'E:\Paul_LaserIgnition\Table_H2_LI.xlsx';
%     % Condition match columns in Excel sheet
%     Configs.ConditionMatchColumn = 5:16;
%     % Top folder of data, e.g., "A_Raw_Data_8_CH4-DualFuel"
%     Configs.DataTopfolder = 'E:\Paul_LaserIgnition\Raw';
%     % Object save location
%     Configs.DataObjSaveDir = "Outputs\H2_LI";  
%     Configs.DataObjReadDir = "Outputs\H2_LI.mat";  
% 
%     Configs.FigureStatisticsOutput = 'Outputs\Statistics';
%     Configs.FigureConditionOutput = 'Outputs\Conditions';
%     Configs.FigureNamePattern = 'H2_LI';
% 
    % Configs.FieldNameExtended = {
    % 'File_Image', 'File_Pressure','File_Log', 'File_Photodiode',...  % Raw data directory
    % 'TimeScale',...            % Time near SOI (-10 to 20 ms centred at SOI (zero), in Class_PressureProcess)
    % 'P_Raw',...                  % Raw pressure near SOI 
    % 'P_Corrected',...         % Corrected pressure near SOI 
    % 'HRR',...                      % HRR near SOI 
    % 'IgnitionDelay',...
    % 'P_TrigTrue',...
    % 'PD_TimeScale',...
    % 'PD_Signal',...
    % 'IntensityProfile',...
    % 'IgnDelay_Intensity',...
    % 'Instant_Image',...
    % };
% 
%     % Discard runs
%     Configs.DiscardRuns = [2,3,4,5 , 15];  % PD saturated
% end

%% H2 Double
% function Configs = Configuration()
%     addpath 'Library_Matlab\'
%     addpath 'Library_Matlab\tdms_package\Version_2p5_Final\v2p5'
%     addpath 'Library_Matlab\tdms_package\Version_2p5_Final\v2p5\tdmsSubfunctions'
% 
%     Excel sheet location
%     Configs.DataSheetDirectory   = 'G:\My Drive\X_DataExcel\20240306_H2_DoubleInjection 1.xlsx';
%     Condition match columns in Excel sheet
%     Configs.ConditionMatchColumn = 5:14;
%     Top folder of data, e.g., "A_Raw_Data_8_CH4-DualFuel"
%     Configs.DataTopfolder = 'E:\Paul_LaserIgnition\Raw';
%     Object save location
%     Configs.DataObjSaveDir = "Outputs\H2_Double";  
%     Configs.DataObjReadDir = "Outputs\H2_Double.mat";  
% 
%     Configs.FigureStatisticsOutput = 'Outputs\Statistics';
%     Configs.FigureConditionOutput = 'Outputs\Conditions';
%     Configs.FigureNamePattern = 'H2_Double';
% 
    % Configs.FieldNameExtended = {
    % 'File_Image', 'File_Pressure','File_Log', 'File_Photodiode',...  % Raw data directory
    % 'TimeScale',...            % Time near SOI (-10 to 20 ms centred at SOI (zero), in Class_PressureProcess)
    % 'P_Raw',...                  % Raw pressure near SOI 
    % 'P_Corrected',...         % Corrected pressure near SOI 
    % 'HRR',...                      % HRR near SOI 
    % 'IgnitionDelay',...
    % 'P_TrigTrue',...
    % 'PD_TimeScale',...
    % 'PD_Signal',...
    % 'IntensityProfile',...
    % 'IgnDelay_Intensity',...
    % 'Instant_Image',...
    % };
% 
%     Discard runs
%     Configs.DiscardRuns = [2,3,4,5 , 15];  % PD saturated
% end


% %% CH4 Dual fuel
% function Configs = Configuration()
%     addpath 'Library_Matlab\'
%     addpath 'Library_Matlab\tdms_package\Version_2p5_Final\v2p5'
%     addpath 'Library_Matlab\tdms_package\Version_2p5_Final\v2p5\tdmsSubfunctions'
% 
%     % Excel sheet location
%     Configs.DataSheetDirectory   = 'G:\My Drive\X_DataExcel\Table_CH4_Diesel_DualFuel.xlsx';
%     % Condition match columns in Excel sheet
%     Configs.ConditionMatchColumn = 5:17;
%     % Top folder of data, e.g., "A_Raw_Data_8_CH4-DualFuel"
%     Configs.DataTopfolder = 'E:\A_Raw_Data_8_CH4-DualFuel';
%     % Object save location
%     Configs.DataObjSaveDir = "E:\A_Raw_Data_8_CH4-DualFuel\Outputs\CH4_n-heptane_DualFuel.mat";  
%     Configs.DataObjReadDir = "E:\A_Raw_Data_8_CH4-DualFuel\Outputs\CH4_n-heptane_DualFuel.mat";  
% 
%     Configs.CommonFiguresOutput = 'E:\A_Raw_Data_8_CH4-DualFuel\Outputs\CommonFigures';
%     Configs.VideoOutput = 'E:\A_Raw_Data_8_CH4-DualFuel\Outputs\Videos';
%     Configs.JournalFigureOutput = 'E:\A_Raw_Data_8_CH4-DualFuel\Outputs\JournalFigures';
%     Configs.MatDataOutput = 'E:\A_Raw_Data_8_CH4-DualFuel\Outputs\MatDataOutput';
%     Configs.CsvDataOutput = 'E:\A_Raw_Data_8_CH4-DualFuel\Outputs\CsvDataOutput';
%     Configs.FigureNamePattern = 'CH4_n-heptane_DualFuel';
% 
%     Configs.FieldNameExtended = {
%     'File_Image', 'File_Pressure','File_Log', 'File_Photodiode',...  % Raw data directory
%     'TimeScale',...            % Time near SOI (-10 to 20 ms centred at SOI (zero), in Class_PressureProcess)
%     'P_Raw',...                  % Raw pressure near SOI 
%     'P_Corrected',...         % Corrected pressure near SOI 
%     'HRR',...                      % HRR near SOI \
%     'P_TrigTrue',...
%     'Trigger_Signal',...
%     'PD_TimeScale',...
%     'PD_Signal',...
%     'Intensity_Profile',...
%     'Frame_At_Ignition',...
%     'IgnitionDelay_Intensity',...
%     'IgnitionDelay_Pressure',...
%     'IgnitionDelay_HRR'
%     };
% 
%     % Discard runs
%     Configs.DiscardRuns = [];
% end

%% Methanol Dual fuel
function Configs = Configuration()
    addpath '..\'
    addpath '..\Library_Matlab\'
    addpath '..\Library_Matlab\tdms_package\Version_2p5_Final\v2p5'
    addpath '..\Library_Matlab\tdms_package\Version_2p5_Final\v2p5\tdmsSubfunctions'

    % Excel sheet location
    Configs.DataSheetDirectory   = "E:\A_Raw_Data_10_Methanol\My_Run_Data.xlsx";
    % Condition match columns in Excel sheet
    Configs.ConditionMatchColumn = 5:19;
    % Top folder of data, e.g., "A_Raw_Data_8_CH4-DualFuel"
    Configs.DataTopfolder = "E:\A_Raw_Data_10_Methanol";
    % Object save location
    Configs.DataObjSaveDir = "E:\A_Raw_Data_10_Methanol\Outputs\Methanol_DualFuel";  
    Configs.DataObjReadDir = Configs.DataObjSaveDir;  

    Configs.CommonFiguresOutput = 'E:\A_Raw_Data_10_Methanol\Outputs\CommonFigures';
    Configs.VideoOutput = 'E:\A_Raw_Data_10_Methanol\Outputs\Videos';
    Configs.JournalFigureOutput = 'E:\A_Raw_Data_10_Methanol\Outputs\JournalFigures';
    Configs.MatDataOutput = 'E:\A_Raw_Data_10_Methanol\Outputs\MatDataOutput';
    Configs.CsvDataOutput = 'E:\A_Raw_Data_10_Methanol\Outputs\CsvDataOutput';
    Configs.FigureNamePattern = 'Methanol_n-heptane_DualFuel';

    Configs.FieldNameExtended = {
    'File_Image', 'File_Pressure','File_Log', 'File_Photodiode',...  % Raw data directory
    'TimeScale',...            % Time near SOI (-10 to 20 ms centred at SOI (zero), in Class_PressureProcess)
    'P_Raw',...                  % Raw pressure near SOI 
    'P_Corrected',...         % Corrected pressure near SOI 
    'HRR',...                      % HRR near SOI \
    'P_TrigTrue',...
    'Trigger_Signal',...
    'PD_TimeScale',...
    'PD_Signal',...
    'Intensity_Profile',...
    'Frame_At_Ignition',...
    'IgnitionDelay_Intensity',...
    'IgnitionDelay_Pressure',...
    'IgnitionDelay_HRR'
    };

    % Discard runs
    Configs.DiscardRuns = [];
end