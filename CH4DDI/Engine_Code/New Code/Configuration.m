%% Chirs H2DDI
function Configs = Configuration()
    addpath 'Library_Matlab\'
    addpath 'Library_Matlab\tdms_package\Version_2p5_Final\v2p5'
    addpath 'Library_Matlab\tdms_package\Version_2p5_Final\v2p5\tdmsSubfunctions'

    % Excel sheet location
    Configs.DataSheetDirectory   = "C:\Users\jun-y\OneDrive - UNSW\CH4DDI\Engine\H2DDI_Trial.xlsx";
    Configs.VariableNameRange = 'A2';   % Row with variable names
    Configs.DataRange = 'A4';           % First row of data  
    Configs.DataSheetName = 'Sheet1';   % Sheet name

    % Top folder of data, e.g., "A_Raw_Data_8_CH4-DualFuel"
    Configs.DataTopfolder = "I:\CH4DDI\Code_development\Chirs data";

    % Object save location
    Configs.DataObjSaveDir = "I:\CH4DDI\Code_development\Chirs data\Outputs\Trial20251208";  
    Configs.DataObjReadDir = Configs.DataObjSaveDir;  

    Configs.CommonFiguresOutput = 'I:\CH4DDI\Code_development\Chirs data\Outputs\CommonFigures';
    Configs.JournalFigureOutput = 'I:\CH4DDI\Code_development\Chirs data\Outputs\JournalFigures';
    Configs.MatDataOutput = 'I:\CH4DDI\Code_development\Chirs data\Outputs\MatDataOutput';
    Configs.CsvDataOutput = 'I:\CH4DDI\Code_development\Chirs data\CsvDataOutput';
    Configs.FigureNamePattern = 'Chirs_H2DDI';

    % Extended table object
    Configs.FieldNameExtended = {
    'File_TDMS_AllTakes',...  % Raw data directory
    'TimeScale',...           % Time near SOI 
    'P_Raw',...               % Raw pressure  
    'P_Corrected',...         % Corrected pressure 
    'aHRR',...
    'IMEP',...
    'IMEP_indie',...
    'CoV_ofIMEP',...
    'CA10_indie',...
    'CA50_indie',...
    'CA90_indie',...
    'ignDelay_indie',...
    'CA10_50_indie',...
    'CA50_90_indie',...
    'CA10_90_indie',...
    'PRR_indie',...
    'aHRR_indie',...
    'Total_HR_indie',...
    'P_peak_indie',...
    'indiEffe_indie',
    };

    % Discard runs
    Configs.DiscardRuns = [];
end

