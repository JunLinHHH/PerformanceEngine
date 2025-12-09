% Generate Statistical Analysis for:
    % Pressure
    % Heat release rate (HRR)
    % Ignition delay
    
% Input:   DataTable (T_DataAll), DataSheet(.xlsx file)
% Output: .png/.pdf image

%clearvars; close all; clc
%% Process Data
%% -- Load data
% If .mat exist
clearvars; close all; clc; 
Configs = Configuration();
Loaded = load(Configs.DataObjReadDir);

obj = Loaded.obj;
T_DataAll = obj.DataMatrix;
T_Condition_ROI = obj.ConditionMatrix;

%%
T_Condition_ROI.Time_Scale = cell(height(T_Condition_ROI), 1);
T_Condition_ROI.P_indie = cell(height(T_Condition_ROI), 1);
%% -- Process Mean and STD
S_StatisticalAnalysis = table2struct(T_Condition_ROI);

% Add fields for structure
S_StatisticalAnalysis(1).Time_Scale                     = {};
S_StatisticalAnalysis(1).P_Mean                     = {};
S_StatisticalAnalysis(1).P_Std                         = {};
S_StatisticalAnalysis(1).P_Indie                      = {};

S_StatisticalAnalysis(1).HRR_Mean                = {};
S_StatisticalAnalysis(1).HRR_Std                    = {};
S_StatisticalAnalysis(1).HRR_Indie                  = {};

S_StatisticalAnalysis(1).PressureIgnitionDelay_Mean  = {};
S_StatisticalAnalysis(1).PressureIgnitionDelay_Std      = {};
S_StatisticalAnalysis(1).PressureIgnitionDelay_Indie    = {};
S_StatisticalAnalysis(1).IntensityIgnitionDelay_Mean  = {};
S_StatisticalAnalysis(1).IntensityIgnitionDelay_Std      = {};
S_StatisticalAnalysis(1).IntensityIgnitionDelay_Indie    = {};
S_StatisticalAnalysis(1).HRRIgnitionDelay_Mean  = {};
S_StatisticalAnalysis(1).HRRIgnitionDelay_Std      = {};
S_StatisticalAnalysis(1).HRRIgnitionDelay_Indie    = {};

S_StatisticalAnalysis(1).ShiftedHRR_PressureID_Mean     = {};
S_StatisticalAnalysis(1).ShiftedHRR_PressureID_Std         = {};
S_StatisticalAnalysis(1).ShiftedHRR_PressureID_Indie       = {};
S_StatisticalAnalysis(1).ShiftedHRR_IntensityID_Mean     = {};
S_StatisticalAnalysis(1).ShiftedHRR_IntensityID_Std         = {};
S_StatisticalAnalysis(1).ShiftedHRR_IntensityID_Indie       = {};

S_StatisticalAnalysis(1).Photodiode_Mean      = {};
S_StatisticalAnalysis(1).Photodiode_Std         = {};
S_StatisticalAnalysis(1).Photodiode_Indie       = {};
S_StatisticalAnalysis(1).ShiftedPD_Mean     = {};
S_StatisticalAnalysis(1).ShiftedPD_Std         = {};
S_StatisticalAnalysis(1).ShiftedPD_Indie       = {};

%% Get Mean and Std for Pressure, HRR and Ignition delay

for ii = 1: height(T_Condition_ROI) 

    ConditionTag = Class_TableManagement.ConditionTagGenerator(T_Condition_ROI, obj, ii);
    % Runs fit this condition
    runsForThisConsition = table2array(T_Condition_ROI(ii,'RunNumbers'));
    runsForThisConsition = runsForThisConsition{1};
    
    
    fprintf('\nCONDITION: #%02d\n',ii);
    fprintf([ConditionTag,'\n']);
    fprintf('\t\tValid %03d runs\n', length(runsForThisConsition));
    ConditionSet = T_DataAll(runsForThisConsition, :);  % Non discarded

    % Calculate Mean and STD
    if ~isempty(ConditionSet)
        % TIme scale
        S_StatisticalAnalysis(ii).Time_Scale = ConditionSet.('TimeScale'){1}*1000;  % unit: ms

        % Pressure
        fprintf('Calculating >>> Pressure <<< mean and std\n')
        [S_StatisticalAnalysis(ii).P_Mean,S_StatisticalAnalysis(ii).P_Std] = ...
            StatisticFromCell(ConditionSet.P_Corrected); % ConditionSet_Valid.P_Corrected: cell arry
        S_StatisticalAnalysis(ii).P_Indie = ConditionSet.P_Corrected;

        % HRR
        fprintf('Calculating >>> HRR <<< mean and std\n')
        [S_StatisticalAnalysis(ii).HRR_Mean, S_StatisticalAnalysis(ii).HRR_Std] = ...
            StatisticFromCell(ConditionSet.HRR);
        S_StatisticalAnalysis(ii).HRR_Indie = ConditionSet.HRR;

        % Pressure Ignition Delay
        fprintf('Calculating >>> Pressure Ignition Delay <<< mean and std\n')
        [S_StatisticalAnalysis(ii).PressureIgnitionDelay_Mean, S_StatisticalAnalysis(ii).PressureIgnitionDelay_Std] = ...
            StatisticFromCell(ConditionSet.IgnitionDelay_Pressure);
        S_StatisticalAnalysis(ii).PressureIgnitionDelay_Indie = ConditionSet.IgnitionDelay_Pressure;

        % HRR Ignition Delay
        fprintf('Calculating >>> HRR Ignition Delay <<< mean and std\n')
        [S_StatisticalAnalysis(ii).HRRIgnitionDelay_Mean, S_StatisticalAnalysis(ii).HRRIgnitionDelay_Std] = ...
            StatisticFromCell(ConditionSet.IgnitionDelay_HRR);
        S_StatisticalAnalysis(ii).HRRIgnitionDelay_Indie = ConditionSet.IgnitionDelay_HRR;

        % Intensity Ignition delay
        % fprtinf('Calculating >>> Intensity Ignition Delay <<< mean and std\n')
        % [S_StatisticalAnalysis(ii).IntensityIgnitionDelay_Mean, S_StatisticalAnalysis(ii).IntensityIgnitionDelay_Std] = ...
        %    StatisticFromCell(ConditionSet_Valid.IgnDelay_Intensity);
        % S_StatisticalAnalysis(ii).IntensityIgnitionDelay_Indie = ConditionSet_NonDiscarded.IgnDelay_Intensity;

        % PD
        fprintf('Calculating >>> Photodiode <<< mean and std\n')
        [S_StatisticalAnalysis(ii).Photodiode_Mean, S_StatisticalAnalysis(ii).Photodiode_Std] = ...
            StatisticFromCell(ConditionSet.PD_Signal);
        S_StatisticalAnalysis(ii).Photodiode_Indie = ConditionSet.PD_Signal;

        % Shift HRR
        % fprtinf('Calculating >>> Shift HRR <<< mean and std\n')
        S_DataPassIn = S_StatisticalAnalysis(ii);
        % 
        % ShiftedHRR_PressureID_Indie = ShiftHRR_PressureIgnitionDelay(S_DataPassIn);
        % 
        % [S_StatisticalAnalysis(ii).ShiftedHRR_PressureID_Mean, S_StatisticalAnalysis(ii).ShiftedHRR_PressureID_Std] = ...
        %     StatisticFromCell(ShiftedHRR_PressureID_Indie);
        % S_StatisticalAnalysis(ii).ShiftedHRR_PressureID_Indie = ShiftedHRR_PressureID_Indie;
        % 
        % ShiftedHRR_IntensityID_Indie = ShiftHRR_IntensityIgnitionDelay(S_DataPassIn);
        % 
        % [S_StatisticalAnalysis(ii).ShiftedHRR_IntensityID_Mean, S_StatisticalAnalysis(ii).ShiftedHRR_IntensityID_Std] = ...
        %     StatisticFromCell(ShiftedHRR_IntensityID_Indie);
        % S_StatisticalAnalysis(ii).ShiftedHRR_IntensityID_Indie = ShiftedHRR_IntensityID_Indie; %% Fixed

        % Shift PD
        % newPD_Indie = ShiftPD(S_DataPassIn);
        % [S_StatisticalAnalysis(ii).ShiftedPD_Mean, S_StatisticalAnalysis(ii).ShiftedPD_Std] = ...
        %     StatisticFromCell(newPD_Indie);
        % S_StatisticalAnalysis(ii).ShiftedPD_Indie = newPD_Indie;        
    end
    pause(0.5)
    close all
end

% Overwrite_ID_Mean = dataTable(:, 15);
% Overwrite_ID_Std = dataTable(:, 16);

%% Plot and save figures
% Three secionts below
% Can run any of them after the previous steps
% Can show run number 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ¯\_( ͡° ͜ʖ ͡°)_/¯✊   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ¯\_( ͡° ͜ʖ ͡°)_/¯✊   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ¯\_( ͡° ͜ʖ ͡°)_/¯✊   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ¯\_( ͡° ͜ʖ ͡°)_/¯✊   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -- Shared properties 
% FigureOutput = fullfile(pwd,'E:\Outputs\figure');
FigureOutput = Configs.FigureStatisticsOutput;
Color = {'k', [0 0.4470 0.7410],[0.6350 0.0780 0.1840],[0.9290 0.6940 0.1250]};
LineWidth = 2;
if ~isempty(ConditionSet)
    TimeScale = ConditionSet.('TimeScale'){1}*1000;  % unit: ms
else
    error('The sub set is empty, cannot get TimeScale reference')
end
TimeScaleClosedLoop = [TimeScale,fliplr(TimeScale)]; 

keyboard

%% -- Section: Pressure
close all
ConditionTag = '';
for ii = 3%length(S_StatisticalAnalysis)
    f_Pressure = figure; grid on; hold on; box on; set(f_Pressure,'color','white')
    % Filename % Modifiled for DH2 test
    ConditionTag = Class_TableManagement.ConditionTagGenerator(S_StatisticalAnalysis, obj, ii);
    % Pressure
    Pressure_Mean = S_StatisticalAnalysis(ii).P_Mean;
    Pressure_Std = S_StatisticalAnalysis(ii).P_Std;
    DataInBetween = [Pressure_Mean - Pressure_Std,...
        fliplr(Pressure_Mean + Pressure_Std)];

    if isempty(Pressure_Mean)
        warning('Empty Pressure Mean data: %d: %s', ii, ConditionTag)
        continue
    end
    
    % Average plot
    plot(TimeScale, Pressure_Mean,'color',Color{1},'LineWidth', LineWidth)
    fill(TimeScaleClosedLoop,DataInBetween,Color{1},'facealpha',0.2,'edgecolor', 'none', 'edgealpha', 0.001,'HandleVisibility','off') 
    xlim([0,15])
    xlabel('Time [ms]')
    ylabel('Pressure [bar]')
    %title(ConditionTag)
    
    % Individual plot
    Pressure_Indie = S_StatisticalAnalysis(ii).P_Indie;
    RunNumbers = S_StatisticalAnalysis(ii).RunNumbers;
    Pressure_GObject = gobjects(length(Pressure_Indie),1);
    for jj = 1:length(Pressure_Indie)
        Pressure_GObject(jj) = plot(TimeScale, Pressure_Indie{jj},'LineStyle','-','LineWidth',1,'Color',Color{2});
        Pressure_GObject(jj).ButtonDownFcn = {@ButtonDownCallBack, Color{2}};
        Pressure_GObject(jj).XDataSource = sprintf('%03d', RunNumbers(jj));
    end

    print(f_Pressure,fullfile(FigureOutput,['Pressure-',ConditionTag]),'-dpng')
    fprintf('Image saved: %s\n', ConditionTag)

    % close all
end

%% -- Section: Heat release rate
close all

for ii = 1:length(S_StatisticalAnalysis)
    f_HRR = figure; grid on; hold on; box on; set(f_HRR,'color','white')
    ConditionTag = Class_TableManagement.ConditionTagGenerator(S_StatisticalAnalysis, obj, ii);
    % HRR
    HRR_Mean = S_StatisticalAnalysis(ii).HRR_Mean;
    HRR_Std = S_StatisticalAnalysis(ii).HRR_Std;
    DataInBetween = [HRR_Mean - HRR_Std,...
        fliplr(HRR_Mean + HRR_Std)];

    if isempty(HRR_Mean)
        warning('Empty HRR Mean data: %d: %s', ii, ConditionTag)
        close all
        continue
    end

    % Individual plot
    HRR_Indie = S_StatisticalAnalysis(ii).HRR_Indie;
    RunNumbers = S_StatisticalAnalysis(ii).RunNumbers;
    HRR_GObject = gobjects(length(HRR_Indie),1);
    for jj = 1:length(HRR_Indie)
        HRR_GObject(jj) = plot(TimeScale, HRR_Indie{jj},'LineStyle','-','LineWidth',1,'Color',Color{2});
        HRR_GObject(jj).ButtonDownFcn = {@ButtonDownCallBack, Color{2}};
        HRR_GObject(jj).XDataSource = sprintf('%03d', RunNumbers(jj));
    end
    % Average plot
    plot(TimeScale, HRR_Mean,'color',Color{1},'LineWidth', LineWidth)
    fill(TimeScaleClosedLoop,DataInBetween,Color{1},'facealpha',0.2,'edgecolor', 'none', 'edgealpha', 0.001,'HandleVisibility','off') 
    
    % Label
    xlim([0,20])
    ylim([-50,1200])
    xlabel('Time [ms]')
    ylabel('HRR [J/ms]')
    %title(ConditionTag)
    
    % Format
    objFig = CLASS_FormatFigure;
        objFig.PositionAxis = round([60 60 270 270] * 0.8);
        objFig.PositionFigure = round([50 50 340 360] * 0.8);
    objFig.SingleYAxis(f_HRR);
    
    print(f_HRR,fullfile(FigureOutput,['HRR-',ConditionTag,'.png']),'-dpng','-r500')
    fprintf('Image saved: %s\n', ConditionTag)
    keyboard
    close all
end

%% -- Section: HRR Ignition Delay
close all
f_IgnitionDelay = figure; grid on; hold on; box on; set(f_IgnitionDelay,'color','white')
IgnDelay_GObject = gobjects(length(S_StatisticalAnalysis),50); % add if runs at one condition > 50
for ii = 1 :length(S_StatisticalAnalysis) % ii is condition index
    
    % Filename
    ConditionTag = Class_TableManagement.ConditionTagGenerator(S_StatisticalAnalysis, obj, ii);

    % Ignition delay
    IgnitionDelay_Mean = S_StatisticalAnalysis(ii).HRRIgnitionDelay_Mean;
    IgnitionDelay_Std = S_StatisticalAnalysis(ii).HRRIgnitionDelay_Std;

    if isempty(IgnitionDelay_Mean)
        warning('Empty IgnitionDelay Mean data: %d: %s', ii, ConditionTag)
        continue
    end

    % Average plot
    plot(ii, IgnitionDelay_Mean,'color',Color{1},'Marker', 'o','MarkerSize',10)
    errorbar(ii, IgnitionDelay_Mean, IgnitionDelay_Std,'k')

    % Individual plot
    IgnitionDelay_Indie = S_StatisticalAnalysis(ii).HRRIgnitionDelay_Indie;
    RunNumbers = S_StatisticalAnalysis(ii).RunNumbers;
    for jj = 1:length(IgnitionDelay_Indie)  % jj is run index at a certain condition
        if ~isempty(IgnitionDelay_Indie{jj})
            IgnDelay_GObject(ii, jj) = plot(ii, IgnitionDelay_Indie{jj},'Marker','x','Color',Color{2});
            IgnDelay_GObject(ii, jj).ButtonDownFcn = {@ButtonDownCallBack,Color{2}};
            IgnDelay_GObject(ii, jj).XDataSource = sprintf('%03d', RunNumbers(jj));
        end
    end
end

xlabel('Case')
ylabel('Ignition delay [ms]')
xlim([0,length(S_StatisticalAnalysis)+1])
ylim([0,20])
thisAxis = gca;
set(thisAxis,'XTick',1:length(S_StatisticalAnalysis))
set(thisAxis,'XTicklabel',1:length(S_StatisticalAnalysis))
set(f_IgnitionDelay,'position',[50,50,1000,600])

keyboard
print(f_IgnitionDelay,fullfile(FigureOutput,['IgnitionDelay-',ConditionTag]),'-dpng')
fprintf('Image saved: %s\n', 'Ignition delay')

%% -- Section: Photodiode
close all
for ii = [11,10,12]
    f_Photodiode = figure; grid on; hold on; box on; set(f_Photodiode,'color','white')
    % Filename
    % Filename % Modifiled for DH2 test
    ConditionTag = Class_TableManagement.ConditionTagGenerator(S_StatisticalAnalysis, obj, ii);

    % Photodiode
    Photodiode_Mean = S_StatisticalAnalysis(ii).Photodiode_Mean;
    Photodiode_Std = S_StatisticalAnalysis(ii).Photodiode_Std;
    DataInBetween = [Photodiode_Mean - Photodiode_Std,...
        fliplr(Photodiode_Mean + Photodiode_Std)];

    if isempty(Photodiode_Mean)
        warning('Empty Photodiode Mean data: %d: %s', ii, ConditionTag)
        close all
        continue
    end

    % Average plot
    plot(TimeScale, Photodiode_Mean,'color',Color{1},'LineWidth', LineWidth)
    fill(TimeScaleClosedLoop,DataInBetween,Color{1},'facealpha',0.2,'edgecolor', 'none', 'edgealpha', 0.001,'HandleVisibility','off') 
    xlim([0,15])
    %ylim([0,10])
    xlabel('Time [ms]')
    ylabel('Photodiode [V]')
    title(ConditionTag)
    
    % Individual plot
    Photodiode_Indie = S_StatisticalAnalysis(ii).Photodiode_Indie;
    RunNumbers = S_StatisticalAnalysis(ii).RunNumbers;
    Photodiode_GObject = gobjects(length(Photodiode_Indie),1);
    for jj = 1:length(Photodiode_Indie)
        thisSet = Photodiode_Indie{jj};
        if ~isempty(thisSet)
            Photodiode_GObject(jj) = plot(TimeScale, thisSet,'LineStyle','-','LineWidth',1,'Color',Color{2});
            Photodiode_GObject(jj).ButtonDownFcn = {@ButtonDownCallBack, Color{2}};
            Photodiode_GObject(jj).XDataSource = sprintf('%03d', RunNumbers(jj));
        end
    end
    
    keyboard
    print(f_Photodiode,fullfile(FigureOutput,['Photodiode-',ConditionTag,'.png']),'-dpng')
    fprintf('Image saved: %s\n', ConditionTag)

    close all
end

%% -- Section: By Conditions - Ignition Delay
VariationGroups = {[3,2,1], [3,11,10], [3,7,6] };
SaveNameVariations = {'journalFig_Autoignition_H2CH4_IgnitionDelay_T_Var',...
    'journalFig_Autoignition_H2CH4_IgnitionDelay_O_Var',...
    'journalFig_Autoignition_H2CH4_IgnitionDelay_P_Var'};
XLabel = {'Ambient gas temperature [K]', 'Ambient oxygen concentration [vol.%]', 'Injection pressure [MPa]'};
XTickValue = {[1200, 1140, 1060], [21, 15, 10],[20, 15, 10]};
close all
for mm = 1:numel(VariationGroups)
    f_IgnitionDelay = figure; grid on; hold on; box on; set(f_IgnitionDelay,'color','white')
    % SaveNameCell = cell(length(S_StatisticalAnalysis), 1);
    % IgnDelay_GObject = gobjects(length(S_StatisticalAnalysis),50); % add if runs at one condition > 50
    dataPosition = 1;
    for ii = VariationGroups{mm} % ii is condition index
        
        ConditionTag = SaveNameVariations{mm};
        % Ignition delay
        IgnitionDelay_Mean = S_StatisticalAnalysis(ii).HRRIgnitionDelay_Mean;
        IgnitionDelay_Std = S_StatisticalAnalysis(ii).HRRIgnitionDelay_Std;

%         % Individual plot at +.5 location
%         IgnitionDelay_Indie = S_StatisticalAnalysis(ii).IgnitionDelay_Indie;
%         RunNumbers = S_StatisticalAnalysis(ii).RunNumber;
%         for jj = 1:length(IgnitionDelay_Indie)  % jj is run index at a certain condition
%             IgnDelay_GObject(ii, jj) = plot(dataPosition + 0.5, IgnitionDelay_Indie{jj},'color',Color{2},'Marker','x','LineStyle','none','LineWidth',1.0,...
%                 'HandleVisibility', 'off');
%             IgnDelay_GObject(ii, jj).ButtonDownFcn = {@ButtonDownCallBack,Color{2}};
%             IgnDelay_GObject(ii, jj).XDataSource = sprintf('%03d', RunNumbers(jj));
%         end
        
        % Average plot
        plot(dataPosition, IgnitionDelay_Mean,'color','#BD0000','Marker','o','LineStyle','none','LineWidth',1.5)
        errorbar(dataPosition, IgnitionDelay_Mean, IgnitionDelay_Std,'vertical','Color','k','LineWidth',1.1)
%         plot(10000, 10000,'color',Color{2},'Marker','x','LineStyle','none','LineWidth',1.0)
        xlabel(XLabel{mm})
        ylabel('Ignition delay [ms]')
        
        LL = legend('Mean value','Uncertainty','location','Northwest');
%         LL = legend('Mean value','Uncertainty','Individual cases','location','Northwest');
        dataPosition = dataPosition + 1;
    end
    
    xlim([0,4])
    ylim([0,9])
    thisAxis = gca;
    set(thisAxis,'XTick',1:length(S_StatisticalAnalysis))
    set(thisAxis,'XTick',[1,2,3])
    set(thisAxis,'XTicklabel',XTickValue{mm})
        objFig = CLASS_FormatFigure;
            objFig.PositionAxis = [60 60 270 270];
            objFig.PositionFigure = [50 50 340 340];
        objFig.SingleYAxis(f_IgnitionDelay);
        CLASS_Utilis.SaveFigureToPDF(ConditionTag, f_IgnitionDelay)
        fprintf('Image saved: case: %d\n', mm)
end

%% -- Section: By Conditions - Shifted HRR DH2 (1030k)
clc; close all

VariationGroups = {[1,2,3,7]};
SaveNameVariations = {'journalFig_DH2_IntensityIgnitionDelay_SHRR_1030k_2ndASOI'};
LineColor = {'#A71B1B','#1A51C6','k','#196F3D'}; LineColor = fliplr(LineColor);
FaceColor = {'r','b','k','g'}; FaceColor = fliplr(FaceColor);
LegendLabel = {{'1 ms - 1 ms - 1 ms', '1 ms - 2 ms - 1 ms', '1 ms - 3 ms - 1 ms', '2 ms - 1 ms - 1 ms'}}
FacePosition =  {[4, 4,-100,-35], [4, 4,-100,-35], [4, 4,-100,-35],[4, 4,-100,-35]};
ErrorbarHeight = [-50,-75,-100,-125];

close all
for mm = 1:numel(VariationGroups)
    f_SHRR = figure; grid minor; hold on; box on; set(f_SHRR,'color','white')
    ConditionTag = SaveNameVariations{mm};
    %SaveName = [SaveName,'wide'];
    SHRR_GObject = gobjects(length(S_StatisticalAnalysis),50); % add if runs at one condition > 50
    idx = 1;

    for ii = VariationGroups{mm} % ii is condition index
        this_S_StatisticalAnalysis = S_StatisticalAnalysis(ii);


        % SHRR and Ignition delay
        ShiftedHRR_Mean = this_S_StatisticalAnalysis.ShiftedHRR_IntensityID_Mean;
        ShiftedHRR_Std = this_S_StatisticalAnalysis.ShiftedHRR_IntensityID_Std;
        IgnitionDelay_Mean = this_S_StatisticalAnalysis.IntensityIgnitionDelay_Mean;
        IgnitionDelay_Std = this_S_StatisticalAnalysis.IntensityIgnitionDelay_Std;

        plot(TimeScale, ShiftedHRR_Mean,'color',LineColor{idx},'LineWidth',1.3)
        DataInBetween = [ShiftedHRR_Mean - ShiftedHRR_Std,...
                fliplr(ShiftedHRR_Mean + ShiftedHRR_Std)];
        fill(TimeScaleClosedLoop,DataInBetween,FaceColor{idx},'facealpha',0.2,'edgecolor', 'none', 'edgealpha', 0.001,'HandleVisibility','off') 
    
         idx_ID = find(TimeScale>IgnitionDelay_Mean,1,'first');
         point_HRR = ShiftedHRR_Mean(idx_ID);
         point_time = TimeScale(idx_ID);

         plot([IgnitionDelay_Mean,IgnitionDelay_Mean], [ErrorbarHeight(idx), -200], 'color', LineColor{idx},'LineWidth',1, 'linestyle','-',...
             'HandleVisibility','off')
%          plot(point_time,point_HRR, 'color', LineColor{idx},'marker','o','LineWidth',1.3,'linestyle','none',...
%              'HandleVisibility','off');
%          plot([IgnitionDelay_Mean,IgnitionDelay_Mean], [ErrorbarHeight(idx), point_HRR], 'color', LineColor{idx},'LineWidth',1, 'linestyle','-',...
%              'HandleVisibility','off')
         errorbar(IgnitionDelay_Mean, ErrorbarHeight(idx),IgnitionDelay_Std,'horizontal','color', LineColor{idx},'LineWidth',1.3,...
             'HandleVisibility','off');
         % EOI
        %  plot([7.86,7.86], [-200, 1000], 'color', 'k','LineWidth',1, 'linestyle','-.',...
        %      'HandleVisibility','off')
        % annotation(f_SHRR,'arrow',[0.7 0.6],...
        %     [0.252 0.252]);
        % text(10, -100,'EOI','fontsize',9,'FontWeight','bold','FontName','helvetica')

         idx = idx + 1;
    end
    xlim([0,12])
    ylim([-150,500])
    xlabel('Time [ms]')
    ylabel('AHRR [J/ms]')
    LG = legend(LegendLabel{mm},'Location','northeast');
    set(LG,'Units','pixels')
    
    objFig = CLASS_FormatFigure;
        objFig.PositionAxis = [60 60 270 270];
        objFig.PositionFigure = [50 50 340 340];
%                 objFig.PositionAxis = [60 60 400 170];
%         objFig.PositionFigure = [50 50 500 240];
    objFig.SingleYAxis(f_SHRR);
    CreateRectangle_4Elements(f_SHRR,LG.Position,FacePosition{mm});
    CLASS_Utilis.SaveFigureToPDF(ConditionTag, f_SHRR)
    fprintf('Image saved: case: %d\n', mm)
   
end

%% -- Section: By Conditions - Shifted HRR DH2_Autoignition 
clc; close all

VariationGroups = {[30,31,32,33]};
SaveNameVariations = {'journalFig_DH2_IntensityIgnitionDelay_SHRR_1030k_'};
LineColor = {'#A71B1B','#1A51C6','k','#196F3D'}; LineColor = fliplr(LineColor);
FaceColor = {'r','b','k','g'}; FaceColor = fliplr(FaceColor);
LegendLabel = {{'1 ms 1030 K', '7 ms 1030 K', '1 ms 1000 K', '7 ms 1000 K'}}
FacePosition =  {[4, 4,-100,-35], [4, 4,-100,-35], [4, 4,-100,-35],[4, 4,-100,-35]};
ErrorbarHeight = [-50,-75,-100,-125];

close all
for mm = 1:numel(VariationGroups)
    f_SHRR = figure; grid minor; hold on; box on; set(f_SHRR,'color','white')
    ConditionTag = SaveNameVariations{mm};
    %SaveName = [SaveName,'wide'];
    SHRR_GObject = gobjects(length(S_StatisticalAnalysis),50); % add if runs at one condition > 50
    idx = 1;

    for ii = VariationGroups{mm} % ii is condition index
        this_S_StatisticalAnalysis = S_StatisticalAnalysis(ii);


        % SHRR and Ignition delay
        ShiftedHRR_Mean = this_S_StatisticalAnalysis.ShiftedHRR_IntensityID_Mean;
        ShiftedHRR_Std = this_S_StatisticalAnalysis.ShiftedHRR_IntensityID_Std;
        IgnitionDelay_Mean = this_S_StatisticalAnalysis.IntensityIgnitionDelay_Mean;
        IgnitionDelay_Std = this_S_StatisticalAnalysis.IntensityIgnitionDelay_Std;

        plot(TimeScale, ShiftedHRR_Mean,'color',LineColor{idx},'LineWidth',1.3)
        DataInBetween = [ShiftedHRR_Mean - ShiftedHRR_Std,...
                fliplr(ShiftedHRR_Mean + ShiftedHRR_Std)];
        fill(TimeScaleClosedLoop,DataInBetween,FaceColor{idx},'facealpha',0.2,'edgecolor', 'none', 'edgealpha', 0.001,'HandleVisibility','off') 
    
         idx_ID = find(TimeScale>IgnitionDelay_Mean,1,'first');
         point_HRR = ShiftedHRR_Mean(idx_ID);
         point_time = TimeScale(idx_ID);

         plot([IgnitionDelay_Mean,IgnitionDelay_Mean], [ErrorbarHeight(idx), -200], 'color', LineColor{idx},'LineWidth',1, 'linestyle','-',...
             'HandleVisibility','off')
%          plot(point_time,point_HRR, 'color', LineColor{idx},'marker','o','LineWidth',1.3,'linestyle','none',...
%              'HandleVisibility','off');
%          plot([IgnitionDelay_Mean,IgnitionDelay_Mean], [ErrorbarHeight(idx), point_HRR], 'color', LineColor{idx},'LineWidth',1, 'linestyle','-',...
%              'HandleVisibility','off')
         errorbar(IgnitionDelay_Mean, ErrorbarHeight(idx),IgnitionDelay_Std,'horizontal','color', LineColor{idx},'LineWidth',1.3,...
             'HandleVisibility','off');
         % EOI
        %  plot([7.86,7.86], [-200, 1000], 'color', 'k','LineWidth',1, 'linestyle','-.',...
        %      'HandleVisibility','off')
        % annotation(f_SHRR,'arrow',[0.7 0.6],...
        %     [0.252 0.252]);
        % text(10, -100,'EOI','fontsize',9,'FontWeight','bold','FontName','helvetica')

         idx = idx + 1;
    end
    xlim([0,18])
    ylim([-150,500])
    xlabel('Time [ms]')
    ylabel('AHRR [J/ms]')
    LG = legend(LegendLabel{mm},'Location','northeast');
    set(LG,'Units','pixels')
    
    objFig = CLASS_FormatFigure;
        objFig.PositionAxis = [60 60 270 270];
        objFig.PositionFigure = [50 50 340 340];
%                 objFig.PositionAxis = [60 60 400 170];
%         objFig.PositionFigure = [50 50 500 240];
    objFig.SingleYAxis(f_SHRR);
    CreateRectangle_4Elements(f_SHRR,LG.Position,FacePosition{mm});
    CLASS_Utilis.SaveFigureToPDF(ConditionTag, f_SHRR)
    fprintf('Image saved: case: %d\n', mm)
   
end
%% -- Section: By Conditions - Shifted HRR DH2 (1000k)
clc; close all

VariationGroups = {[4,5,6],[8,9,10],[42,43,44]};
SaveNameVariations = {'journalFig_DH2_IntensityIgnitionDelay_SHRR_1000k_1ms_2ndASOI',...
    'journalFig_DH2_IntensityIgnitionDelay_SHRR_1000k_2ms_2ndASOI',...
    'journalFig_DH2_IntensityIgnitionDelay_SHRR_1000k_3ms_2ndASOI'};
LineColor = {'#A71B1B','#1A51C6','k'}; LineColor = fliplr(LineColor);
FaceColor = {'r','b','k'}; FaceColor = fliplr(FaceColor);
LegendLabel = {{'1 ms - 1 ms - 1 ms', '1 ms - 2 ms - 1 ms', '1 ms - 3 ms - 1 ms'},...
    {'2 ms - 1 ms - 1 ms', '2 ms - 2 ms - 1 ms', '2 ms - 3 ms - 1 ms'},...
    {'3 ms - 1 ms - 1 ms', '3 ms - 2 ms - 1 ms', '3 ms - 3 ms - 1 ms'}};
FacePosition =  {[4, 4,-100,-35], [4, 4,-100,-35], [4, 4,-100,-35],[4, 4,-100,-35]};
ErrorbarHeight = [-50,-75,-100,-125];

close all
for mm = 1:numel(VariationGroups)
    f_SHRR = figure; grid minor; hold on; box on; set(f_SHRR,'color','white')
    ConditionTag = SaveNameVariations{mm};
    %SaveName = [SaveName,'wide'];
    SHRR_GObject = gobjects(length(S_StatisticalAnalysis),50); % add if runs at one condition > 50
    idx = 1;

    for ii = VariationGroups{mm} % ii is condition index
        this_S_StatisticalAnalysis = S_StatisticalAnalysis(ii);


        % SHRR and Ignition delay
        ShiftedHRR_Mean = this_S_StatisticalAnalysis.ShiftedHRR_IntensityID_Mean;
        ShiftedHRR_Std = this_S_StatisticalAnalysis.ShiftedHRR_IntensityID_Std;
        IgnitionDelay_Mean = this_S_StatisticalAnalysis.IntensityIgnitionDelay_Mean;
        IgnitionDelay_Std = this_S_StatisticalAnalysis.IntensityIgnitionDelay_Std;

        plot(TimeScale, ShiftedHRR_Mean,'color',LineColor{idx},'LineWidth',1.3)
        DataInBetween = [ShiftedHRR_Mean - ShiftedHRR_Std,...
                fliplr(ShiftedHRR_Mean + ShiftedHRR_Std)];
        fill(TimeScaleClosedLoop,DataInBetween,FaceColor{idx},'facealpha',0.2,'edgecolor', 'none', 'edgealpha', 0.001,'HandleVisibility','off') 
    
         idx_ID = find(TimeScale>IgnitionDelay_Mean,1,'first');
         point_HRR = ShiftedHRR_Mean(idx_ID);
         point_time = TimeScale(idx_ID);

         plot([IgnitionDelay_Mean,IgnitionDelay_Mean], [ErrorbarHeight(idx), -200], 'color', LineColor{idx},'LineWidth',1, 'linestyle','-',...
             'HandleVisibility','off')
%          plot(point_time,point_HRR, 'color', LineColor{idx},'marker','o','LineWidth',1.3,'linestyle','none',...
%              'HandleVisibility','off');
%          plot([IgnitionDelay_Mean,IgnitionDelay_Mean], [ErrorbarHeight(idx), point_HRR], 'color', LineColor{idx},'LineWidth',1, 'linestyle','-',...
%              'HandleVisibility','off')
         errorbar(IgnitionDelay_Mean, ErrorbarHeight(idx),IgnitionDelay_Std,'horizontal','color', LineColor{idx},'LineWidth',1.3,...
             'HandleVisibility','off');
         % EOI
        %  plot([7.86,7.86], [-200, 1000], 'color', 'k','LineWidth',1, 'linestyle','-.',...
        %      'HandleVisibility','off')
        % annotation(f_SHRR,'arrow',[0.7 0.6],...
        %     [0.252 0.252]);
        % text(10, -100,'EOI','fontsize',9,'FontWeight','bold','FontName','helvetica')

         idx = idx + 1;
    end
    xlim([0,12])
    ylim([-150,500])
    xlabel('Time [ms]')
    ylabel('AHRR [J/ms]')
    LG = legend(LegendLabel{mm},'Location','northeast');
    set(LG,'Units','pixels')
    
    objFig = CLASS_FormatFigure;
        objFig.PositionAxis = [60 60 270 270];
        objFig.PositionFigure = [50 50 340 340];
%                 objFig.PositionAxis = [60 60 400 170];
%         objFig.PositionFigure = [50 50 500 240];
    objFig.SingleYAxis(f_SHRR);
    CreateRectangle_3Element(f_SHRR,LG.Position,FacePosition{mm});
    CLASS_Utilis.SaveFigureToPDF(ConditionTag, f_SHRR)
    fprintf('Image saved: case: %d\n', mm)
   
end
%% -- Section: By Conditions - Shifted HRR (All-in-one figure)
clc; close all

VariationGroups = {[1,2,3]};
ConditionTag = 'journalFig_DoubleInjection_H2_SHRR_AllInOne';
LineColor = {'#A71B1B','#1A51C6','k'}; LineColor = fliplr(LineColor);
FaceColor = {'r','b','k'}; FaceColor = fliplr(FaceColor);
LegendLabel = {{'1200 K', '1140 K', '1060 K'}, {'21%', '15%', '10%'}, {'20 MPa', '15 MPa', '10 MPa'}};
FacePosition =  {[4, 4,-45,-35], [4, 4,-33,-35], [4, 4,-48,-35]};
ErrorbarHeight = [-50,-100,-150] -25;
close all

% Build axes matrix
AxisResizeRatio = 1;
AxesHandles = CLASS_AxesHandleStore;
AxesHandles.RowNumber = 3;     % make this flexibe 
AxesHandles.ColumnNumber = 1;  % make this flexibe 
AxesHandles.MarginRight = 10;
AxesHandles.MarginTop = 10;
AxesHandles.GapRow = 3;                 % gap between rows
AxesHandles.GapColumn = 20;            % gap between columns
AxesHandles.MarginLeft = 50;
AxesHandles.MarginBottom = 40;
AxesHandles.AxesWidth = floor((340-60) * AxisResizeRatio);
AxesHandles.AxesHeight = floor(170 * AxisResizeRatio);

AxesHandles = AxesHandles.ConstuctAxes();
AxesMatrix = AxesHandles.GetAxesHandleMatrix();

%
for mm = 1:numel(VariationGroups)
    axes(AxesMatrix(mm));
    grid minor; hold on; box on; set(gcf,'color','white')
    %SaveName = [SaveName,'wide'];
    SHRR_GObject = gobjects(length(S_StatisticalAnalysis),50); % add if runs at one condition > 50
    idx = 1;

    for ii = VariationGroups{mm} % ii is condition index
        this_S_StatisticalAnalysis = S_StatisticalAnalysis(ii);
        
        % SHRR and Ignition delay
        ShiftedHRR_Mean = this_S_StatisticalAnalysis.ShiftedHRR_PressureID_Mean;
        ShiftedHRR_Std = this_S_StatisticalAnalysis.ShiftedHRR_PressureID_Std;
        IgnitionDelay_Mean = this_S_StatisticalAnalysis.PressureIgnitionDelay_Mean;
        IgnitionDelay_Std = this_S_StatisticalAnalysis.PressureIgnitionDelay_Std;

        plot(TimeScale, ShiftedHRR_Mean,'color',LineColor{idx},'LineWidth',1.3)
        DataInBetween = [ShiftedHRR_Mean - ShiftedHRR_Std,...
                fliplr(ShiftedHRR_Mean + ShiftedHRR_Std)];
        fill(TimeScaleClosedLoop,DataInBetween,FaceColor{idx},'facealpha',0.2,'edgecolor', 'none', 'edgealpha', 0.001,'HandleVisibility','off') 
    
         idx_ID = find(TimeScale>IgnitionDelay_Mean,1,'first');
         point_HRR = ShiftedHRR_Mean(idx_ID);
         point_time = TimeScale(idx_ID);

         plot([IgnitionDelay_Mean,IgnitionDelay_Mean], [ErrorbarHeight(idx), -250], 'color', LineColor{idx},'LineWidth',1, 'linestyle','-',...
             'HandleVisibility','off')
%          plot(point_time,point_HRR, 'color', LineColor{idx},'marker','o','LineWidth',1.3,'linestyle','none',...
%              'HandleVisibility','off');
%          plot([IgnitionDelay_Mean,IgnitionDelay_Mean], [ErrorbarHeight(idx), point_HRR], 'color', LineColor{idx},'LineWidth',1, 'linestyle','-',...
%              'HandleVisibility','off')
         errorbar(IgnitionDelay_Mean, ErrorbarHeight(idx),IgnitionDelay_Std,'horizontal','color', LineColor{idx},'LineWidth',1.3,...
             'HandleVisibility','off');
         % EOI
         plot([7.86,7.86], [-250, 1000], 'color', 'k','LineWidth',1, 'linestyle','-.',...
             'HandleVisibility','off')
         idx = idx + 1;
    end
    xlim([0,15])
    ylim([-250,700])
    
    if mm == 2
        ylabel('AHRR [J/ms]')
    elseif mm == 3
        xlabel('Time [ms]')
        text(10, -150,'EOI','fontsize',9,'FontWeight','bold','FontName','helvetica')
        annotation(gcf,'arrow',[0.68 0.58], [0.1 0.1]);
    end
    if mm ~= 3
        AxesMatrix(mm).XTickLabel = [];
    end
    AxesMatrix(mm).YTickLabel{1} = '';
    LG = legend(LegendLabel{mm},'Location','northeast');
    set(LG,'Units','pixels')
    
    CreateRectangle_4Elements(gcf,LG.Position, FacePosition{mm});
end
CLASS_Utilis.SaveFigureToPDF(ConditionTag, gcf)
fprintf('Image saved: case: %d\n', mm)
%% -- Section: By Conditions - Photodiode
clc; close all
VariationGroups = {[11, 10 ,12]};
SaveNameVariations = {'journalFig_DF_CH4_PD_T_Var',...
    'journalFig_DF_CH4_PD_O_Var',...
    'journalFig_DF_CH4_PD_P_Var'};
LineColor = {'#A71B1B','#1A51C6','k'}; LineColor = fliplr(LineColor);
FaceColor = {'r','b','k'}; FaceColor = fliplr(FaceColor);
LegendLabel = {{'4.4 MPa', '5.2 MPa', '6.3 MPa'}, {'21%', '15%', '10%'}, {'20 MPa', '15 MPa', '10 MPa'}};
FacePosition =  {[4, 4,-50,-35], [4, 4,-33,-35], [4, 4,-48,-35]};
close all

for mm = 1:numel(VariationGroups)  % only plot P inj var
    f_PD = figure; grid on; hold on; box on; set(f_PD,'color','white')
    ConditionTag = SaveNameVariations{mm};
    PD_GObject = gobjects(length(S_StatisticalAnalysis),50); % add if runs at one condition > 50
    idx = 1;

    for ii = VariationGroups{mm} % ii is condition index
        this_S_StatisticalAnalysis = S_StatisticalAnalysis(ii);
        
        % SHRR and Ignition delay
        Photodiode_Mean = this_S_StatisticalAnalysis.Photodiode_Mean;
        Photodiode_Std = this_S_StatisticalAnalysis.Photodiode_Std;
        plot(TimeScale, Photodiode_Mean,'color',LineColor{idx},'LineWidth',1.3)
        DataInBetween = [Photodiode_Mean - Photodiode_Std,...
                fliplr(Photodiode_Mean + Photodiode_Std)];
        fill(TimeScaleClosedLoop,DataInBetween,FaceColor{idx},'facealpha',0.2,'edgecolor', 'none', 'edgealpha', 0.001,'HandleVisibility','off') 

        % % EOI
        % plot([7.86,7.86], [-200, 1000], 'color', 'k','LineWidth',1, 'linestyle','-.',...
        %      'HandleVisibility','off')
        % annotation(f_PD,'arrow',[0.7 0.6],...
        %     [0.252 0.252]);
        % text(10, -0.3,'EOI','fontsize',9,'FontWeight','bold','FontName','helvetica')

         idx = idx + 1;
    end
    xlim([0,10])
    ylim([-1,3])
    xlabel('Time [ms]')
    ylabel('Intensity [V]')
    LG = legend(LegendLabel{mm},'Location','northeast');
    set(LG,'Units','pixels')
    
    objFig = CLASS_FormatFigure;
        objFig.PositionAxis = [60 60 270 270] * 0.8;
        objFig.PositionFigure = [50 50 340 340] * 0.8;
    objFig.SingleYAxis(f_PD);
    CreateRectangle_3Element(f_PD,LG.Position, FacePosition{mm});
    CLASS_Utilis.SaveFigureToPDF(fullfile(Configs.FigureConditionOutput, ConditionTag), f_PD)
    fprintf('Image saved: case: %d\n', mm)
end
%
%% -- Section: By Conditions - Shifted PD
clc; close all
VariationGroups = {[3,2,1], [3,11,10], [3,7,6] };
SaveNameVariations = {'journalFig_Autoignition_H2CH4_SPD_T_Var',...
    'journalFig_Autoignition_H2CH4_SPD_O_Var',...
    'journalFig_Autoignition_H2CH4_SPD_P_Var'};
LineColor = {'#A71B1B','#1A51C6','k'}; LineColor = fliplr(LineColor);
FaceColor = {'r','b','k'}; FaceColor = fliplr(FaceColor);
LegendLabel = {{'1200 K', '1140 K', '1060 K'}, {'21%', '15%', '10%'}, {'20 MPa', '15 MPa', '10 MPa'}};
FacePosition =  {[4, 4,-45,-35], [4, 4,-33,-35], [4, 4,-48,-35]};
ErrorbarHeight = [-0.2,-0.5,-0.8];
close all
for mm = 1:3
    f_SPD = figure; grid minor; hold on; box on; set(f_SPD,'color','white')
    ConditionTag = SaveNameVariations{mm};
    %SaveName = [SaveName,'wide'];
    SPD_GObject = gobjects(length(S_StatisticalAnalysis),50); % add if runs at one condition > 50
    idx = 1;

    for ii = VariationGroups{mm} % ii is condition index
        this_S_StatisticalAnalysis = S_StatisticalAnalysis(ii);
        
        % SPD and Ignition delay
        ShiftedPD_Mean = this_S_StatisticalAnalysis.ShiftedPD_Mean;
        ShiftedPD_Std = this_S_StatisticalAnalysis.ShiftedPD_Std;
        IgnitionDelay_Mean = this_S_StatisticalAnalysis.IgnitionDelay_Mean;
        IgnitionDelay_Std = this_S_StatisticalAnalysis.IgnitionDelay_Std;

        plot(TimeScale, ShiftedPD_Mean,'color',LineColor{idx},'LineWidth',1.3)
        DataInBetween = [ShiftedPD_Mean - ShiftedPD_Std,...
                fliplr(ShiftedPD_Mean + ShiftedPD_Std)];
        fill(TimeScaleClosedLoop,DataInBetween,FaceColor{idx},'facealpha',0.2,'edgecolor', 'none', 'edgealpha', 0.001,'HandleVisibility','off') 
    
         idx_ID = find(TimeScale>IgnitionDelay_Mean,1,'first');
         point_PD = ShiftedPD_Mean(idx_ID);
         point_time = TimeScale(idx_ID);

         plot([IgnitionDelay_Mean,IgnitionDelay_Mean], [ErrorbarHeight(idx), -1], 'color', LineColor{idx},'LineWidth',1, 'linestyle','-',...
             'HandleVisibility','off')
%          plot(point_time,point_PD, 'color', LineColor{idx},'marker','o','LineWidth',1.3,'linestyle','none',...
%              'HandleVisibility','off');
%          plot([IgnitionDelay_Mean,IgnitionDelay_Mean], [ErrorbarHeight(idx), point_PD], 'color', LineColor{idx},'LineWidth',1, 'linestyle','-',...
%              'HandleVisibility','off')
         errorbar(IgnitionDelay_Mean, ErrorbarHeight(idx),IgnitionDelay_Std,'horizontal','color', LineColor{idx},'LineWidth',1.3,...
             'HandleVisibility','off');
         % EOI
         plot([7.86,7.86], [-200, 1000], 'color', 'k','LineWidth',1, 'linestyle','-.',...
             'HandleVisibility','off')
        annotation(f_SPD,'arrow',[0.7 0.6],...
            [0.252 0.252]);
        text(10, -0.3,'EOI','fontsize',9,'FontWeight','bold','FontName','helvetica')

         idx = idx + 1;
    end
    xlim([0,15])
    ylim([-1,6])
    xlabel('Time [ms]')
    ylabel('Intensity [V]')
    LG = legend(LegendLabel{mm},'Location','northeast');
    set(LG,'Units','pixels')
    
    objFig = CLASS_FormatFigure;
        objFig.PositionAxis = [60 60 270 270];
        objFig.PositionFigure = [50 50 340 340];
%                 objFig.PositionAxis = [60 60 400 170];
%         objFig.PositionFigure = [50 50 500 240];
    objFig.SingleYAxis(f_SPD);
    CreateRectangle_4Elements(f_SPD,LG.Position, FacePosition{mm});
    CLASS_Utilis.SaveFigureToPDF(ConditionTag, f_SPD)
    fprintf('Image saved: case: %d\n', mm)
end

%
fprintf('Script >> %s << is finished\n', mfilename());
%% Internal functions
%% -- Click callback
function ButtonDownCallBack(hObj,~,ColorNonActive) % src is ~

    thisFig = gcf;
    SelectType = thisFig.SelectionType;

    if strcmp(SelectType,'normal')
        fprintf('Clicked data is #%s\n',hObj.XDataSource)
        hObj.Color = [1,0,0];
    elseif strcmp(SelectType,'alt')
        hObj.Color = ColorNonActive;
    end
end

%% -- Calculate mean and std from a cell
function [MeanData,StdData] = StatisticFromCell(C_Input)
    MeanData = [];
    StdData = [];
    % Check if the input cell is empty
    Flags = ~cellfun(@isempty, C_Input);

    if sum(Flags) == 0
        fprintf("- - - Input cell is empty\n")
        return
    end
    
    % Find the data size
    idx = find(Flags == 1,1,'first');
    SumData = zeros(numel(C_Input), length(C_Input{idx}));
    for ii = 1:numel(C_Input)
        if ~isempty(C_Input{ii})
            SumData(ii,:) = C_Input{ii};
        else
            SumData(ii,:) = NaN;
            % fprintf('Empty...\n')
        end
    end
    MeanData = mean(SumData,1,'omitnan');
    StdData = std(SumData,1,'omitnan');

%     f = figure; grid on; hold on; box on; set(f,'color','white')
%     for ii = 1:numel(C_Input)
%         plot(C_Input{ii});
%     end
%     plot(MeanData,'color','k','LineWidth',3,'LineStyle','--','Marker','o','MarkerSize',1);
end

%% -- Shift HRR
function newHRR_Indie = ShiftHRR_PressureIgnitionDelay(DataPassIn)
    SampingRate = 200000/1000; % sampling per ms
    IgnitionDelay_Mean = DataPassIn.PressureIgnitionDelay_Mean;
    IgnitionDelay_Indie = DataPassIn.PressureIgnitionDelay_Indie;
    HRR_Indie = DataPassIn.HRR_Indie;
    newHRR_Indie = cell(size(HRR_Indie));
    RunAmount = length(IgnitionDelay_Indie);
    % Circshift HRR
    for jj = 1:RunAmount
            thisIgnitionDelay = IgnitionDelay_Indie{jj};
            if isempty(thisIgnitionDelay)
                warning('Error in shifting HRR, Pressure ignition delay is empty')
                return
            end                
            thisHRR_Indie = HRR_Indie{jj};
            ID_difference = floor((IgnitionDelay_Mean - thisIgnitionDelay) * SampingRate);
            newHRR = circshift(thisHRR_Indie, ID_difference);   
            newHRR_Indie{jj} = newHRR;
    end
end

function newHRR_Indie = ShiftHRR_IntensityIgnitionDelay(DataPassIn)
    SampingRate = 200000/1000; % sampling per ms
    IgnitionDelay_Mean = DataPassIn.IntensityIgnitionDelay_Mean;
    IgnitionDelay_Indie = DataPassIn.IntensityIgnitionDelay_Indie;
    HRR_Indie = DataPassIn.HRR_Indie;
    newHRR_Indie = cell(size(HRR_Indie));
    RunAmount = length(IgnitionDelay_Indie);
    % Circshift HRR
    for jj = 1:RunAmount
            thisIgnitionDelay = IgnitionDelay_Indie{jj};
            if  isempty(thisIgnitionDelay)
                    warning('Error in shifting HRR, Intensity ignition delay is empty')
                    return
            end
            thisHRR_Indie = HRR_Indie{jj};
            ID_difference = floor((IgnitionDelay_Mean - thisIgnitionDelay) * SampingRate);
            newHRR = circshift(thisHRR_Indie, ID_difference);   
            newHRR_Indie{jj} = newHRR;
    end
end

%% -- Shift PD
function newPD_Indie = ShiftPD(DataPassIn)
    SampingRate = 200000/1000; % sampling per ms
    IgnitionDelay_Mean = DataPassIn.PressureIgnitionDelay_Mean;
    IgnitionDelay_Indie = DataPassIn.PressureIgnitionDelay_Indic;
    PD_Indie = DataPassIn.Photodiode_Indie;
    newPD_Indie = cell(size(PD_Indie));
    RunAmount = length(IgnitionDelay_Indie);
    % Circshift PD
    for jj = 1:RunAmount
            thisIgnitionDelay = IgnitionDelay_Indie{jj};
            if ~isempty(thisIgnitionDelay)
                thisPD_Indie = PD_Indie{jj};
                ID_difference = floor((IgnitionDelay_Mean - thisIgnitionDelay) * SampingRate);
                newPD = circshift(thisPD_Indie, ID_difference);   
                newPD_Indie{jj} = newPD;
            else
                fprintf('Error in shifting PD, Ignition delay is empty')
            end
    end
end

%% -- CreateRectangle For 4 elements
function CreateRectangle_4Elements(HandleFigure, L_position, offset)
    %CREATERECTANGLE(figure1)
    %  FIGURE1:  annotation figure

    %  Auto-generated by MATLAB on 08-Jul-2021 16:15:43

    % Create rectangle
    color_shade = {'r','b','k','g'};
    Step_height = 12.5;

    %color_shade = fliplr(color_shade);
    L_position = L_position + offset;

    L_base = L_position(2);
    stepY = 14;
    for ii = 1:length(color_shade)

        L_position(2) = L_base + stepY*(ii-1);
        L_position(4) = Step_height;
        
        annotation(HandleFigure,'rectangle',...
            'unit','pixel','position',...
            L_position,...
            'LineStyle','none',...
            'FaceColor', color_shade{ii},...
            'FaceAlpha',0.2);
    end
end
%% -- CreateRectangle For 3 elements
function CreateRectangle_3Element(HandleFigure, L_position, offset)
    %CREATERECTANGLE(figure1)
    %  FIGURE1:  annotation figure

    %  Auto-generated by MATLAB on 08-Jul-2021 16:15:43

    % Create rectangle
    color_shade = {'r','b','k'};

    %color_shade = fliplr(color_shade);
    L_position = L_position + offset;

    L_base = L_position(2);
    stepY = 14;
    for ii = 1:length(color_shade)

        L_position(2) = L_base + stepY*(ii-1);
        
        annotation(HandleFigure,'rectangle',...
            'unit','pixel','position',...
            L_position,...
            'LineStyle','none',...
            'FaceColor', color_shade{ii},...
            'FaceAlpha',0.2);
    end
end