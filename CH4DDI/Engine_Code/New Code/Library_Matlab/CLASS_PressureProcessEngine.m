classdef CLASS_PressureProcessEngine
    % CLASS_PRESSUREPROCESSENGINE - Process ONE Set (multiple Takes)
    % Can discard individual Takes
    
    properties
        SetNumber                            % Set number
        TakeFiles                            % Cell array of TDMS file paths
        DiscardedTakes                       % Logical array: true = discarded
        
        % Channel names in TDMS
        PressureChannel = 'Ai7';              % Pressure signal
        PilotFuelTriggerChannel = 'Di3';      % Diesel injection signal
        MainFuelTriggerChannel = 'Di4';       % H2 injection signal
        EncoderDegreeChannel = 'Di6';         % Encoder degree
        EncoderZPulseChannel = 'Di7';         % Encoder Z-pulse
        
        % Engine parameters (fixed for this engine)
        DV = 0.005890/6;    % Displacement volume [m^3]
        CR = 17.4;          % Compression ratio
        crlen = 0.2221;     % Connecting rod length [m]
        crad = 0.0625;      % Crank radius [m] (half of stroke)
        B = 0.100;          % Bore [m]
        TDC_shift = -67;    % TDC position shift
        
        % Derived parameters
        Lambda              % crad / crlen
        CA                  % Crank angle array: (0:0.1:719.9)-360
        Vol                 % Volume vs crank angle
    end
    
    properties (Access = private)
        PressureInfo        % Structure array: one element per Take
    end
    
    methods
        %% Constructor
        function this = CLASS_PressureProcessEngine(tdmsFiles, setNumber, varargin)
            % Constructor - Create processor for ONE Set
            %
            % Usage:
            %   processor = CLASS_PressureProcessEngine(tdmsFiles)
            %   processor = CLASS_PressureProcessEngine(tdmsFiles, setNumber)
            %   processor = CLASS_PressureProcessEngine(tdmsFiles, setNumber, 'PropertyName', value, ...)
            %
            % Inputs:
            %   tdmsFiles - Cell array of TDMS file paths for this Set
            %   setNumber - (optional) Set number for tracking
            %
            % Optional Parameters:
            %   'PressureChannel' - Default: 'Ai7'
            %   'PilotFuelTriggerChannel' - Default: 'Di3'
            %   'MainFuelTriggerChannel' - Default: 'Di4'
            %   'EncoderDegreeChannel' - Default: 'Di6'
            %   'EncoderZPulseChannel' - Default: 'Di7'
            %   'DV' - Displacement volume, Default: 0.005890/6
            %   'CR' - Compression ratio, Default: 17.4
            %   'crlen' - Connecting rod length, Default: 0.2221
            %   'crad' - Crank radius, Default: 0.0625
            %   'B' - Bore, Default: 0.100
            %   'TDC_shift' - TDC shift, Default: -67
            %
            % Examples:
            %   % Default settings
            %   processor = CLASS_PressureProcessEngine(tdmsFiles, 1);
            %
            %   % Custom channels
            %   processor = CLASS_PressureProcessEngine(tdmsFiles, 1, ...
            %       'PressureChannel', 'Ai0', ...
            %       'PilotFuelTriggerChannel', 'Di0');
            %
            %   % Custom engine parameters
            %   processor = CLASS_PressureProcessEngine(tdmsFiles, 1, ...
            %       'DV', 0.00098, ...
            %       'CR', 17.0);
            
            % Validate input
            if nargin < 1 || isempty(tdmsFiles)
                error('Must provide TDMS files (cell array)');
            end
            
            if ~iscell(tdmsFiles)
                tdmsFiles = {tdmsFiles};
            end
            
            % Store files
            this.TakeFiles = tdmsFiles;
            this.SetNumber = [];
            if nargin >= 2 && ~isempty(setNumber)
                this.SetNumber = setNumber;
            end
            
            % Parse optional parameters
            if ~isempty(varargin)
                for i = 1:2:length(varargin)
                    paramName = varargin{i};
                    paramValue = varargin{i+1};
                    
                    % Check if property exists
                    if isprop(this, paramName)
                        this.(paramName) = paramValue;
                    else
                        warning('Unknown property: %s', paramName);
                    end
                end
            end
            
            % Initialize discard array (all valid initially)
            this.DiscardedTakes = false(length(tdmsFiles), 1);
            
            % Calculate derived engine parameters
            this.Lambda = this.crad / this.crlen;
            this.CA = (0:0.1:719.9) - 360;
            this.Vol = this.DV * (1/(this.CR-1) + 0.5*(1 + 1/this.Lambda * ...
                (1 - sqrt(1 - (this.Lambda*sind(this.CA)).^2)) - cosd(this.CA)));
            
            % Initialize PressureInfo as structure array
            this.PressureInfo = struct('DataRaw', {}, 'FileName', {}, 'TakeNumber', {});
        end
        
        %% Read Pressure Files from TDMS
        function this = ReadPressureFile_TDMS(this)
            % ReadPressureFile_TDMS - Read all Takes for this Set
            
            nTakes = length(this.TakeFiles);
            fprintf('Reading %d Takes for Set %s...\n', nTakes, mat2str(this.SetNumber));
            
            % Read each Take
            for j = 1:nTakes
                fprintf('  Take %02d: ', j-1);
                
                try
                    % Read TDMS file
                    DataRaw = CLASS_Utilis.ReadTdmsData(this.TakeFiles{j});
                    
                    % Store in PressureInfo
                    this.PressureInfo(j).DataRaw = DataRaw;
%                     this.PressureInfo(j).FileName = this.TakeFiles{j};
%                     this.PressureInfo(j).TakeNumber = j - 1;  % 0-indexed
                    
                    fprintf('✓\n');
                    
                catch ME
                    fprintf('✗ Failed: %s\n', ME.message);
                    
                    % Store empty
                    this.PressureInfo(j).DataRaw = [];
%                     this.PressureInfo(j).FileName = this.TakeFiles{j};
%                     this.PressureInfo(j).TakeNumber = j - 1;
                end
            end
            
            fprintf('Done\n\n');
        end
        
        %% Discard Take
        function this = DiscardTake(this, takeNumber)
            % DiscardTake - Mark a Take as discarded
            %
            % Usage:
            %   processor = processor.DiscardTake(0)  % Discard Take 0
            %   processor = processor.DiscardTake(1)  % Discard Take 1
            
            % Convert 0-indexed to 1-indexed
            idx = takeNumber + 1;
            
            if idx < 1 || idx > length(this.TakeFiles)
                warning('Take %d does not exist', takeNumber);
                return;
            end
            
            this.DiscardedTakes(idx) = true;
            fprintf('Take %02d marked as DISCARDED\n', takeNumber);
        end
        
        %% Recover Take
        function this = RecoverTake(this, takeNumber)
            % RecoverTake - Unmark a Take (make it valid again)
            %
            % Usage:
            %   processor = processor.RecoverTake(0)  % Recover Take 0
            
            % Convert 0-indexed to 1-indexed
            idx = takeNumber + 1;
            
            if idx < 1 || idx > length(this.TakeFiles)
                warning('Take %d does not exist', takeNumber);
                return;
            end
            
            this.DiscardedTakes(idx) = false;
            fprintf('Take %02d marked as VALID\n', takeNumber);
        end
        
        %% Get Valid Takes
        function validTakes = GetValidTakes(this)
            % GetValidTakes - Get PressureInfo for non-discarded Takes only
            %
            % Returns:
            %   Structure array with only valid Takes
            
            validIdx = ~this.DiscardedTakes;
            validTakes = this.PressureInfo(validIdx);
        end
        
        %% Get Discarded Takes
        function discardedTakes = GetDiscardedTakes(this)
            % GetDiscardedTakes - Get PressureInfo for discarded Takes
            
            discardedIdx = this.DiscardedTakes;
            discardedTakes = this.PressureInfo(discardedIdx);
        end
        
        %% Display Status
        function DisplayStatus(this)
            % DisplayStatus - Show Take status
            
            nTotal = length(this.TakeFiles);
            nDiscarded = sum(this.DiscardedTakes);
            nValid = nTotal - nDiscarded;
            
            fprintf('\n');
            fprintf(repmat('=', 1, 60));
            fprintf('\nSet %s - Take Status\n', mat2str(this.SetNumber));
            fprintf(repmat('=', 1, 60));
            fprintf('\n');
            fprintf('Total Takes:     %d\n', nTotal);
            fprintf('Valid Takes:     %d\n', nValid);
            fprintf('Discarded Takes: %d\n', nDiscarded);
            fprintf(repmat('=', 1, 60));
            fprintf('\n\n');
        end
        
        %% Getter - All Takes (including discarded)
        function S_PressureInfo = GetPressureInfo(this)
            % GetPressureInfo - Access ALL PressureInfo (including discarded)
            % Use GetValidTakes() to get only valid Takes
            
            S_PressureInfo = this.PressureInfo;
        end
    end
end
