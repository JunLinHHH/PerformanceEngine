classdef CLASS_PressureProcess
    % Read/Process raw pressure data (.tdms file)
    % Read run conditions from .tdms log files
    % Calculate:
    %   Corrected pressure
    %   HRR
    %   Ignition delay

    properties
        FileName

        SamplingRate = 200000;
        InjectionDelay = 1000 * 11 / 37500; % unit: ms !!! hard-coded (11 frames at 37500 FPS)
        VoltageToBarRatio = 15
        Specific_Heat_Ratio = 1.35;
        Volume_CVCC = 1.4e-3;  %m^3

        TimeLowerLimit = -10;
        TimeUpperLimit = 20;
        
        RawPressureSignal
        FilteredPressureSignal
        RawMainTriggerSignal
        FilteredMainTriggerSignal
        reZeroTimeScale

        bPlotIgnitionDelay = 0
        bPlotPressureCorrected = 0
    end

    properties (Access = private)
        PressureInfo
    end

    methods
        %% Constructor
        function this = CLASS_PressureProcess()
            % User input
            S_PressureInfo.InjectionDelay = this.InjectionDelay;
            % From log file 
            S_PressureInfo.TriggerPressureSet = [];
            S_PressureInfo.O2Percentage = [];
            S_PressureInfo.PressureOffset = [];
            % Processed data
            S_PressureInfo.TimeScale = [];
            S_PressureInfo.Pressure = [];
            S_PressureInfo.HRR = [];
            S_PressureInfo.P_Corrected = [];
            S_PressureInfo.IgnitionDelay = [];
            S_PressureInfo.TriggerPressureTrue = [];

            this.PressureInfo = S_PressureInfo;
        end
        %% Getter
        function S_PressureInfo = GetPressureInfo(this)
             % access the PressureInfo
             S_PressureInfo = this.PressureInfo;
        end
        
        %% Non-Static methods
        % Input must be a File Directory or a one row Struct includes field named
        % 'File_Pressure'
        function this = PressureProcess(this, dir_tdms)
            S_PressureInfo = this.PressureInfo;
            % Process .tdms pressure file
            % Read files


            % if isempty(dir_file)   % Open browsing
            %     [file,path] = uigetfile('*.tdms');
            %     FilePressure_Tdms = fullfile(path, file);
            % else
            %     if istable(dir_file)         % input is one row table
            %         isTableCol = @(t, thisCol) ismember(thisCol, t.Properties.VariableNames);
            %         if isTableCol(dir_file,'File_Pressure')
            %             FilePressure_Tdms = dir_file.('File_Pressure'){1};
            %         else
            %             warning("Input table has no column named File_Pressure ")
            %         end
            %     elseif isstruct(dir_file) % input is one row struct
            %         if isfield(dir_file,'File_Pressure')
            %             FilePressure_Tdms = dir_file.File_Pressure;
            %         else
            %             warning("Input struct has no field named File_Pressure ")
            %             return
            %         end
            %     elseif isfile(dir_file)     % input is full file name
            %         FilePressure_Tdms = dir_file;
            %     end
            % end
        
            if isempty(dir_tdms) || ~isfile(dir_tdms)
                warning('Input is not valid pressure file (.tdms)')
                return
            else
                this.FileName = dir_tdms;    % Store FileName to Object
            end
            
            % Get information from .log file
            % LogFile = strrep(FilePressure_Tdms,'.tdms','.log');
            % 
            % try
            %     S_LogInfo = this.FUNCTION_ReadPressureLogFile(LogFile);
            %     PressureOffset = S_LogInfo.PressureOffset;
            %     TriggerPressure = S_LogInfo.TriggerPressure;
            %     O2Percentage = S_LogInfo.O2Percentage;
            % 
            %     % Outputs store from .log file
            %     S_PressureInfo.TriggerPressureSet = TriggerPressure;
            %     S_PressureInfo.O2Percentage = O2Percentage;
            % 
            % catch
            %     warning('Error for reading log file ( .log):\n %s', LogFile);
            %     return
            % end

            [PressureOffset, TriggerPressure, O2Percentage] = this.ReadPressureLogFile();
            % Outputs store from .log file
            S_PressureInfo.PressureOffset = PressureOffset;
            S_PressureInfo.TriggerPressureSet = TriggerPressure;
            S_PressureInfo.O2Percentage = O2Percentage;
            
            % Read tdms file
            try
                this = this.ReadPressureFile_TDMS();
            catch
                warning('Error for reading pressure file (.tdms):\n %s', dir_tdms);
                return
            end

            % Process pressure data -> Get corrected pressure & HRR
            try
                [TimeScale_SOI, Pressure_SOI, HRR_SOI, Pressure_SOI_Corrected, Pressure_TrigTrue,TriggerSignal_AtSOI] =...
                    this.ReadPressureSignal(PressureOffset);
                
                % Outputs store from .tdms file
                S_PressureInfo.TimeScale = TimeScale_SOI;
                S_PressureInfo.Pressure = Pressure_SOI;
                S_PressureInfo.HRR = HRR_SOI;
                S_PressureInfo.P_Corrected = Pressure_SOI_Corrected;
                S_PressureInfo.TriggerPressureTrue = Pressure_TrigTrue;
                S_PressureInfo.Trigger_Signal = TriggerSignal_AtSOI;
            catch
                warning('Error for pressure processing: %s', dir_tdms);
                return
            end
            
            % Calculate ignition delay
            try
                IgnitionDelay = this.IgnitionDelayFromCorrectedPressure(Pressure_SOI_Corrected);
                S_PressureInfo.IgnitionDelay = IgnitionDelay;
            catch
                warning('Error for ignition delay processing: %s', dir_tdms);
                return
            end
            
            this.PressureInfo = S_PressureInfo;
        end
        %-------------------------------------------------------------------------------------%
        % Read tdms file
        function this = ReadPressureFile_TDMS(this)
                DataRaw = CLASS_Utilis.ReadTdmsData(this.FileName); % read data from .tdms file (needs tdms read package)
                
                % Pre-proccesing from raw data
                    % Covert voltage to bar (15 bar/V)
                SignalInjection= DataRaw.Untitled.Injection_trigger.data;             % filter and find trigger signal
                SignalPressure = double(DataRaw.Untitled.Pressure__bar_.data*15);     % convert voltage to bar
                
                this.RawPressureSignal = SignalPressure;
                this.RawMainTriggerSignal = SignalInjection;  % Keep raw noisy digital data
                
        end
        
        % This somehow gives Nan at the end of pressure trace
        function this = ReadPressureFile_TDMS_NEW(this)
            dir_tdms = this.FileName;
            data_tdms = tdmsread(dir_tdms); % read data from .tdms file (needs tdms read package)
            data_tdms_table = data_tdms{1};
            % Pre-proccesing from raw data
                % Covert voltage to bar (15 bar/V)
            SignalInjection = data_tdms_table.("Injection trigger");
            SignalPressure = data_tdms_table.("Pressure [bar]");                            % filter and find trigger signal
            SignalPressure = double(SignalPressure* this.VoltageToBarRatio);     % convert voltage to bar
            
            this.RawPressureSignal = SignalPressure;
            this.RawMainTriggerSignal = SignalInjection;  % Keep raw noisy digital data
        end

        %-------------------------------------------------------------------------------------%
        % Process HRR
        function [TimeScale_AtSOI, Pressure_AtSOI, HRR_AtSOI, Pressure_AtSOI_Corrected, Pressure_TrigTrue, TriggerSignal_AtSOI] =...
            ReadPressureSignal(this, PressureOffset)
                % PressureOffset: bar

                SamplingRate_Local = this.SamplingRate;   % sample/sec
                InjectionDelay_Local = this.InjectionDelay; % ms
                
                if ~isempty(PressureOffset)
                    SignalPressure = this.RawPressureSignal + PressureOffset;
                else
                    SignalPressure = this.RawPressureSignal;
                    fprintf('Pressure offset from log file is empty. This can cause incorrect pressure trace.\n')
                end
                

                % Filter Main trigger signal
                this.FilteredMainTriggerSignal = medfilt1(double(this.RawMainTriggerSignal),10);        % Filter trigger signal % store to property

                SignalMainTrigger = this.FilteredMainTriggerSignal;
                SignalTimeScale = (0:length(SignalPressure)-1)/SamplingRate_Local;                         % time step with sampling rate (unit:s)
                
                % Find electrical injection trigger time
                % Reset time zero to time aSOI
                MaxSignalInjection = max(SignalMainTrigger); 
                MinSignalInjection = min(SignalMainTrigger);
                ThresholdTrigger = MinSignalInjection + (MaxSignalInjection - MinSignalInjection)/2;

                SignalTimeScale_reZero = SignalTimeScale - SignalTimeScale(find(SignalMainTrigger > ThresholdTrigger,1,'first'));      
                IndexTrigger = find(SignalTimeScale_reZero == 0);
                IndexSOI = IndexTrigger + floor(SamplingRate_Local * InjectionDelay_Local / 1000);
                SignalTimeScale_reZero = SignalTimeScale_reZero - SignalTimeScale_reZero(IndexSOI); % reset time 0 to SOI 
                this.reZeroTimeScale = SignalTimeScale_reZero;   % store to property
                
                % Filtering pressure signal
                SignalPressure_Filtered = bandstop(SignalPressure,[9000 15000],SamplingRate_Local);
                SignalPressure_Filtered = bandstop(SignalPressure_Filtered,[35000 60000],SamplingRate_Local);
                SignalPressure_Filtered = sgolayfilt(SignalPressure_Filtered,5,201);
                SignalPressure_Filtered = bandstop(SignalPressure_Filtered,[3500 4000],SamplingRate_Local);
                this.FilteredPressureSignal = SignalPressure_Filtered; % store to property
        
                % Calculate & filter HRR
                HeatReleaseRate = 1/(specific_heat_ratio-1) * volume_CVCC * diff(SignalPressure_Filtered) * 1e5 * SamplingRate_Local / 1000;
                HeatReleaseRate = sgolayfilt(HeatReleaseRate,3,201);
                HeatReleaseRate = [0, HeatReleaseRate];
                
                % Calculate corrected pressure
                SignalPressure_Corrected = this.PressureCorrection();
                
                % Crop data to -10 to 20 ms centred at SOI (zero) 
                timeLowerLimit = this.TimeLowerLimit * SamplingRate_Local/1000;
                timeUpperLimit = this.TimeUpperLimit * SamplingRate_Local/1000;
        
                % Return
                if (IndexSOI+timeLowerLimit) < 0
                    warning('Trggering index is incorrect, check raw trigger signal')
                end
                TimeScale_AtSOI = SignalTimeScale_reZero(IndexSOI+timeLowerLimit : IndexSOI+timeUpperLimit);
                Pressure_AtSOI = SignalPressure_Filtered(IndexSOI+timeLowerLimit : IndexSOI+timeUpperLimit);
                HRR_AtSOI = HeatReleaseRate(IndexSOI+timeLowerLimit : IndexSOI+timeUpperLimit);
                Pressure_AtSOI_Corrected = SignalPressure_Corrected(IndexSOI+timeLowerLimit : IndexSOI+timeUpperLimit);
                Pressure_TrigTrue = SignalPressure(IndexSOI); % Pressure at SOI
                TriggerSignal_AtSOI = SignalMainTrigger(IndexSOI+timeLowerLimit : IndexSOI+timeUpperLimit);
        end
        
        %%
        %-------------------------------------------------------------------------------------%
        % Isolate fuel injection/combustion event, rule out heat loss to wall 
        function Pressure_Corrected = PressureCorrection(this)
                % Input is filtered pressure, entire time range (re zero)
                SignalTimeScale_reZero = this.reZeroTimeScale;  % in sec
                SignalPressure_Filtered = this.FilteredPressureSignal;

                SamplingRate_Local = this.SamplingRate;
                FitRatio = 0.1; % Decide how many data points are used to fit
                        
                % Find the index of pressure rising due to fuel combustion 
                Index_Start = this.FindPressureRisingIndex();
                
                FitIndexRange = Index_Start-SamplingRate_Local*FitRatio : Index_Start;
                FitData_Time = SignalTimeScale_reZero(FitIndexRange);
                FitData_Pressure = SignalPressure_Filtered(FitIndexRange);
                
                % Be aware of the matrix shape
                [FitResult, ~] = fit(FitData_Time', FitData_Pressure', 'exp1');
                Pressure_Corrected = SignalPressure_Filtered - FitResult(SignalTimeScale_reZero)';
                        
                % Plot (Optional)
                if this.bPlotPressureCorrected
                    figure; hold on; grid on; box on; set(gcf,'color','white'); set(gca,'FontSize',15);
                    f = gcf; % f.Name = file(ii).name;
                    set(f,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);
            
                    yyaxis left
                    plot(SignalTimeScale_reZero, SignalPressure_Filtered,'k-')  % Entire time range
                    plot(SignalTimeScale_reZero(FitIndexRange), SignalPressure_Filtered(FitIndexRange),'bx')
                    plot(SignalTimeScale_reZero(Index_Start), SignalPressure_Filtered(Index_Start),'ro')
                    plot(SignalTimeScale_reZero, FitResult(SignalTimeScale_reZero),'r-.')
                    plot(SignalTimeScale_reZero, this.FilteredMainTriggerSignal* 100, 'Color','#555EA5','LineStyle','-','Marker','none') % Entire time range
                    
                    ylabel('Original pressure [bar]');
                    
                    yyaxis right
                    plot(SignalTimeScale_reZero, Pressure_Corrected,'Color','#809A80','LineStyle','-')
                    ylabel('Corrected pressure [bar]');
                    %xlim([-300,100])
                    legend('Pressure (filtered)','Region for fitting','Pressure rising start','Exp-fit','Main trigger signal', 'Corrected pressure', 'location','northeast')
                    xlabel('Time [ms]')
                    set(f,'Position',[50,50,800,800])

                    [pks, locs] = findpeaks(-SignalPressure_Filtered,SignalTimeScale_reZero,'MinPeakProminence',0.4,'Annotate','extents','MinPeakDistance',1);
                    plot(SignalTimeScale_reZero, -SignalPressure_Filtered,'g-')
                    plot(locs, pks, 'rv'); 
                    legend('Pressure (filtered)','Region for fitting','Pressure rising start','Exp-fit','Main trigger signal', 'Corrected pressure', ...
                        'Inverted pressure',...
                        'Peaks (must be 3)', ...
                        'location','northeast')
                    titleName =  strrep(this.FileName(strfind(this.FileName, 'run')-9 :strfind(this.FileName, 'run')+6),'_','-');
                    title(titleName)
                    drawnow
                end
        end
        
        %%
        % Find the index of pressure rising due to fuel combustion 
        function PressureRisingIndex = FindPressureRisingIndex(this)
            % Input is filtered pressure, entire time range (re zero)
            SignalTimeScale_reZero = this.reZeroTimeScale;
            SignalPressure_Filtered = this.FilteredPressureSignal;

            % Filp the pressure and find peaks
            % 0.5 is around pressure rising magnitude
            % 1 is for remove peaks on noisy regions
            [~,locs,~,~] =  findpeaks(-SignalPressure_Filtered,SignalTimeScale_reZero,'MinPeakProminence',0.4,'Annotate','extents','MinPeakDistance',1);
            % Run
            % findpeaks(-Pressure_Input,Time_Input,'MinPeakProminence',0.5,'Annotate','extents','MinPeakDistance',1)
            % without; to show plot
            % prevent collpse
            if isempty(locs)
                    PressureRisingIndex = 1; 
                    return
            end
            
            % Criterion 1: peaks at start and end of pressure history due to Noise and the Triggering
            if length(locs) == 3
                    Time_Trigger = locs(2);
                    PressureRisingIndex = find(SignalTimeScale_reZero > Time_Trigger,1,'first');
            % Criterion 2: find the one close to time 0 (time zero should at triggering)
            else
                    locs(locs<0) = 10; % remove negative number, 
                    [~,idx_locs] = min(locs);
                    Time_Trigger = locs(idx_locs);
                    PressureRisingIndex = find(SignalTimeScale_reZero > Time_Trigger,1,'first');
            end
        end
        
        %%
        % Find Ignition delay from Corrected Pressure data
        function IgnitionDelayinMS = IgnitionDelayFromCorrectedPressure(this, Pressure_Corrected)
            timeLowerLimit = this.TimeLowerLimit;
            SamplingRate_Local = this.SamplingRate; 
            IndexForTrigger = SamplingRate_Local/1000 * abs(timeLowerLimit);  % pressure data contains -10 ms from trigger 
            
            NoiseBeforeTriggering = mean(Pressure_Corrected(1:IndexForTrigger));
            STDNoiseBeforeTriggering = std(Pressure_Corrected(1:IndexForTrigger));
            Threshold2STD = 2*STDNoiseBeforeTriggering + NoiseBeforeTriggering + 0.05; % 0.05 is original thershold offset 
            
            IgnitionDelayinMS = (find(Pressure_Corrected>Threshold2STD,1,'first')-IndexForTrigger) / (SamplingRate_Local / 1000);
        
            if this.bPlotIgnitionDelay
                figure; hold on; grid on; box on; set(gcf,'color','white'); set(gca,'FontSize',9);
                plot( (1:length(Pressure_Corrected))/(SamplingRate_Local/1000) + timeLowerLimit, Pressure_Corrected);
                xline(IgnitionDelayinMS);
                xlabel('Time [ms]')
                ylabel('Pressure [bar]')
                xlim([0,15])
                ylim([-1,4])
                legend('Pressure corrected','Ignition delay')
            end
        end
    
        %% Plot
        function ShowPlot(this)
                S_PressureInfo = this.PressureInfo;
%                 S_PressureInfo.TriggerPressureSet = [];
%                 S_PressureInfo.O2Percentage = [];
%                 S_PressureInfo.InjectionDelay = [];
%                 
%                 S_PressureInfo.TimeScale = [];
%                 S_PressureInfo.Pressure = [];
%                 S_PressureInfo.HRR = [];
%                 S_PressureInfo.P_Corrected = [];
%                 S_PressureInfo.IgnitionDelay = [];
%                 S_PressureInfo.TriggerPressureTrue = [];
    
                f1 = figure; box on; set(f1,'color','white')
                f1.Name = sprintf('Trigger Pressure: %02d, O2: %02d', S_PressureInfo.TriggerPressureSet, S_PressureInfo.O2Percentage);
                hold on
                Ax_signalPlot = gca;
                yyaxis left
                Ax_signalPlot.YColor = '#1A51C6';
                plot(S_PressureInfo.TimeScale* 1000, S_PressureInfo.HRR,'color','#1A51C6','LineWidth',1.2);
                plot([S_PressureInfo.IgnitionDelay,S_PressureInfo.IgnitionDelay],[-1000,1000],'k-');
                ylim([-100,900])
                ylabel('HRR [J/ms]')

                yyaxis right
                Ax_signalPlot.YColor = '#A71B1B';
                plot(S_PressureInfo.TimeScale * 1000, S_PressureInfo.P_Corrected,'color','#A71B1B','LineWidth',1.2);
                ylim([-1,9])
                xlim([-10,20])
                ylabel('Pressure [bar]')
                xlabel('TIme [ms]')
                legend('HRR','Ignition delay','Pressure (corrected)')
                objFig = CLASS_FormatFigure;
                    objFig.PositionAxis = [60 60 270 270];
                    objFig.PositionFigure = [50 50 380 340];
                objFig.SingleYAxis(f1);
        end
        
        % Read log file
        function  [PressureOffset, TriggerPressure, O2Percentage] = ReadPressureLogFile(this)
            % Get information from .log file
            this_pressure_file = this.FileName;
            if isempty(this_pressure_file)
                warning('No pressure file')
                return
            end
            dir_log_file = strrep(this_pressure_file,'.tdms','.log');
            
            S_LogInfo = this.FUNCTION_ReadPressureLogFile(dir_log_file);
            PressureOffset = S_LogInfo.PressureOffset;
            TriggerPressure = S_LogInfo.TriggerPressure;
            O2Percentage = S_LogInfo.O2Percentage;
        end
    end

    %% Static Methods
    methods (Static)
        function S_LogInfo = FUNCTION_ReadPressureLogFile(FileName)
            S_LogInfo.PressureOffset = [];
            S_LogInfo.O2Percentage = [];
            S_LogInfo.TriggerPressure = [];
        
            S_LogInfo.P_C2H2 = [];
            S_LogInfo.P_H2 = [];
            S_LogInfo.P_O2 = [];
            S_LogInfo.P_N2 = [];
        
            % Return Empty 
            f = fopen(FileName);
            if f == -1
                fprintf('File: %s, is not valid\n', FileName)
                return
            end
            
            % Get O2 condition from target N2 percentage
            LogText = textscan(f,'%s','Delimiter','\t');
            N2Percentage = str2double(LogText{1}{30});
            switch(N2Percentage)
                case 49.2
                    O2Percentage = 21;
                case 53.97
                    O2Percentage = 18;
                case 58.72
                    O2Percentage = 15;
                case 66.68
                    O2Percentage = 10;
                case 82.55
                    O2Percentage = 0;
                otherwise
                    O2Percentage = [];
            end
            %%
            % Outpus
            S_LogInfo.PressureOffset = str2double(LogText{1}{15});
            S_LogInfo.O2Percentage = O2Percentage;
            S_LogInfo.TriggerPressure = str2double(LogText{1}{41});
            % Get actual pressure for different injected gas
            S_LogInfo.P_C2H2 = str2double(LogText{1}{20});
            S_LogInfo.P_H2 = str2double(LogText{1}{24});
            S_LogInfo.P_O2 = str2double(LogText{1}{28});
            S_LogInfo.P_N2 = str2double(LogText{1}{32});
        end
    end
end