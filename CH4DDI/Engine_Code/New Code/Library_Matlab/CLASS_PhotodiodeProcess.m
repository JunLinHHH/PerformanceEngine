classdef CLASS_PhotodiodeProcess
    properties
        FileName
        SamplingRate = 200000;
        InjectionDelay = 1000 * 11 / 37500; % unit: ms !!! hard-coded (11 frames at 37500 FPS)
        Gain_PD = 0;

        RawInjectionSignal
        RawPhotodiodeSignal
    end

    properties (Access = private)
        PhotodiodeInfo
    end

    methods
        %% Constructor
        function this = CLASS_PhotodiodeProcess()
            S_PhotodiodeInfo.TimeScale = [];
            S_PhotodiodeInfo.Photodiode = [];

            this.PhotodiodeInfo = S_PhotodiodeInfo;
        end

        %% Getter
        function S_PhotodiodeInfo = GetPhotodiodeInfo(this)
             % access the PhotodiodeInfo
             S_PhotodiodeInfo = this.PhotodiodeInfo;
        end
        
        %% Non-Static methods
        % Input must be a File Directory or a one row Struct includes field named
        % 'File_Photodiode'
        function this = PhotodiodeProcess(this, dir_tdms)
            S_PhotodiodeInfo = this.PhotodiodeInfo;
            Gain_PD_Local = this.Gain_PD;
            % Process .tdms Photodiode file
            % Read files
   
            % if istable(DataInput) % input is one row table
            %     isTableCol = @(t, thisCol) ismember(thisCol, t.Properties.VariableNames);
            %     if isTableCol(DataInput,'File_Photodiode')
            %         FilePhotodiode_Tdms = DataInput.('File_Photodiode'){1};
            %     else
            %         warning("Input table has no column named File_Photodiode ")
            %     end
            % elseif isstruct(DataInput) % input is one row struct
            %     if isfield(DataInput,'File_Photodiode')
            %         FilePhotodiode_Tdms = DataInput.File_Photodiode;
            %     else
            %         warning("Input struct has no field named File_Photodiode ")
            %         return
            %     end
            % elseif isfile(DataInput) % input is full file name
            %     FilePhotodiode_Tdms = DataInput;
            % end
        
            if isempty(dir_tdms) || ~isfile(dir_tdms)
                warning('Input is not valid photodiode file (.tdms)')
                return
            else
                this.FileName = dir_tdms;    % Store FileName to Object
            end
            
            % Read tdms file
            try
                this = this.ReadTdmsFile();
            catch
                warning('Error for reading .tdms: %s', dir_tdms);
            end

            % Process tdms
            try
                [TimeScale_AtSOI, Photodiode_AtSOI] = ReadPhotodiodeSignal(this, Gain_PD_Local);
                S_PhotodiodeInfo.TimeScale = TimeScale_AtSOI;
                S_PhotodiodeInfo.Photodiode = Photodiode_AtSOI;
            catch
                warning('Error for photodiode file processing: %s', dir_tdms);
            end

            % Update
            this.PhotodiodeInfo = S_PhotodiodeInfo;
        end

        %-------------------------------------------------------------------------------------%
        % Read tdms file
        function this = ReadTdmsFile(this)
                DataRaw = CLASS_Utilis.ReadTdmsData(this.FileName); % read data from .tdms file (needs tdms read package)
                
                % Pre-proccesing from raw data
                % Old channel names
                % SignalInjection= medfilt1(double(DataRaw.Digital_channels.Untitled_7.data),10);             % filter and find trigger signal
                % SignalPhotodiode = double(DataRaw.Analog_channels.Photodiode.data);                        % convert voltage to bar

                SignalInjection= medfilt1(double(DataRaw.Digital_Channels.Dev0_Di7.data),10);             % filter and find trigger signal
                SignalPhotodiode = double(DataRaw.Analog_Channels.Dev0_Ai0.data);                      

                this.RawPhotodiodeSignal = SignalPhotodiode;
                this.RawInjectionSignal = SignalInjection;
        end
        
        %-------------------------------------------------------------------------------------%
        % Process PD signal
        function [TimeScale_AtSOI, Photodiode_AtSOI] = ReadPhotodiodeSignal(this, Gain_PD)
            SamplingRate_Local = this.SamplingRate;   % sample/sec
            InjectionDelay_Local = this.InjectionDelay; % ms
            SignalPhotodiode = this.RawPhotodiodeSignal;
            SignalInjection = this.RawInjectionSignal;

            SignalTimeScale = (0:length(SignalPhotodiode)-1)/SamplingRate_Local;                         % time step with sampling rate (unit:s)

            % >>>>> Reset time zero to time aSOI
            MaxSignalInjection = max(SignalInjection); 
            MinSignalInjection = min(SignalInjection);
            ThresholdTrigger = MinSignalInjection + (MaxSignalInjection - MinSignalInjection)/2;

            SignalTimeScale = SignalTimeScale - SignalTimeScale(find(SignalInjection > ThresholdTrigger,1,'first'));      
            IndexTrigger = find(SignalTimeScale == 0);
            IndexSOI = IndexTrigger + floor(SamplingRate_Local * InjectionDelay_Local / 1000);
            SignalTimeScale = SignalTimeScale - SignalTimeScale(IndexSOI);                                  % reset time 0 to SOI 
            
            % >>>>> Filtering
            % % Unfilter
            SignalPhotodiode_Filtered = SignalPhotodiode;
           
            % >>>>> output
            % Crop data to -10 to 20 ms centred at SOI (zero) 
            TimeLowerLimit = -10 * SamplingRate_Local/1000;
            TimeUpperLimit = 20 * SamplingRate_Local/1000;
            TimeScale_AtSOI = SignalTimeScale(IndexSOI+TimeLowerLimit : IndexSOI+TimeUpperLimit);
            Photodiode_AtSOI = SignalPhotodiode_Filtered(IndexSOI+TimeLowerLimit : IndexSOI+TimeUpperLimit);
            
            % >>>>> Normalise gain
            Photodiode_AtSOI = Photodiode_AtSOI/(sqrt(10)^(Gain_PD/10));
        end
    end
end