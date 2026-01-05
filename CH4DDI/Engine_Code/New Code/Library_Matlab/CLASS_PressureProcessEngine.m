classdef CLASS_PressureProcessEngine < handle
    % CLASS_PRESSUREPROCESSENGINE - Engine pressure data processing for multi-take datasets
    %
    % Purpose:
    %   Process engine combustion pressure data from TDMS files for ONE Set
    %   containing multiple Takes. Supports per-take discard/recovery and
    %   comprehensive cycle-resolved thermodynamic analysis.
    %
    % Main Responsibilities:
    %   - Read TDMS files containing pressure, encoder, and injection signals
    %   - Detect TDC from Z-pulse and encoder signals
    %   - Extract individual cycles aligned to crank angle
    %   - Calculate IMEP, CoV, peak pressure, pressure rise rate
    %   - Calculate apparent HRR, cumulative HRR, and combustion phasing (CA10/50/90)
    %   - Support per-cycle metrics and ensemble averaging
    %
    % Data Flow:
    %   Constructor → ReadPressureFile_TDMS → ProcessAll:
    %     1. ProcessPressureTraces (filter)
    %     2. FindTDC (encoder-based detection)
    %     3. ExtractCycles (CA-aligned cycles)
    %     4. CalculateIMEP (thermodynamic work)
    %     5. CalculateHRR (combustion phasing)
    %
    % Output Storage:
    %   All results stored in PressureInfo struct array (one per Take)
    %   Access via GetPressureInfo() or GetValidTakes()
    %
    % See also: CLASS_PressureProcess (CVCC variant)
    
    properties
        % File and dataset identification
        SetNumber                            % Set number for this dataset
        TakeFiles                            % Cell array of TDMS file paths
        DiscardedTakes                       % Logical array marking discarded takes
        
        % Channel configuration
        PressureChannel = 'Ai5';             % Pressure analog channel name
        PilotFuelTriggerChannel = 'Di3';     % Pilot (diesel) injection digital signal
        MainFuelTriggerChannel = 'Di4';      % Main (H2) injection digital signal
        EncoderDegreeChannel = 'Di6';        % 0.2° encoder pulse channel
        EncoderZPulseChannel = 'Di7';        % Z-pulse (once-per-rev) channel
        
        % Engine geometry (fixed for this engine)
        DV = 0.005890/6;    % Displacement volume [m^3]
        CR = 17.4;          % Compression ratio
        crlen = 0.2221;     % Connecting rod length [m]
        crad = 0.0625;      % Crank radius [m] (half of stroke)
        B = 0.100;          % Bore [m]
        TDC_shift = -68;    % TDC position shift (hardware-specific) [CAD]
        
        % Derived engine parameters (computed in constructor)
        Lambda              % Crank ratio: crad / crlen
        CA                  % Crank angle array: -360 to +360 [CAD]
        Vol                 % Cylinder volume vs crank angle [m^3]

        % Signal processing parameters
        Filter_CutoffFreq = 10000;   % Pressure filter cutoff frequency [Hz]
        aHRRFilterOrder = 9;         % Savitzky-Golay filter order for aHRR
        aHRRFilterWindow = 201;      % Savitzky-Golay filter window size [samples]
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
            %   'TDC_shift' - TDC shift, Default: -68
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

            % enforce column vectors for consistent math
            this.CA = this.CA(:);
            this.Vol = this.Vol(:);
            
            % Initialize PressureInfo as structure array
            this.PressureInfo = struct('DataRaw', {}, 'FileName', {}, 'TakeNumber', {});
        end
        
        %%  File I/O Methods
        function this = ReadPressureFile_TDMS(this)
            % ReadPressureFile_TDMS - Read all TDMS files for this Set
            % Loads raw data from TDMS files and stores in PressureInfo
            
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
                    fprintf('✓\n');
                    
                catch ME
                    fprintf('✗ Failed: %s\n', ME.message);
                    
                    % Store empty
                    this.PressureInfo(j).DataRaw = [];
                end
            end
            
            fprintf('Done\n\n');
        end
        
        %% Signal Processing Methods
        function P_Corrected = PressureFilter(this, P_Raw, Filter_CutoffFreq, rpm)
            % PressureFilter - Apply 2nd-order Butterworth low-pass filter
            % Legacy implementation matching TDMS_Processing.m lines 190-193
            %
            % Inputs:
            %   P_Raw             - Raw pressure trace [bar]
            %   Filter_CutoffFreq - Cutoff frequency [Hz]
            %   rpm               - Engine speed for normalization
            %
            % Output:
            %   P_Corrected - Filtered pressure [bar]

            if nargin < 4 || isempty(rpm) || rpm <= 0
                rpm = 1400;
            end

            fc_norm = Filter_CutoffFreq /(7200*(rpm/60/2));  % normalized cutoff

            [b, a] = butter(2, fc_norm);
            P_Corrected = filtfilt(b, a, P_Raw);
        end

        function [raw_traces, filt_traces, time_traces] = ProcessPressureTraces(this, tableRow)
            % ProcessPressureTraces - Extract and filter pressure for all Takes
            % Handles channel auto-discovery and applies Butterworth filtering
            %
            % Input:
            %   tableRow - Current row from T_DataAll (contains RPM, etc.)
            %
            % Outputs:
            %   raw_traces  - Cell array of raw pressure traces
            %   filt_traces - Cell array of filtered pressure traces
            %   time_traces - Cell array of time vectors

            rpm = tableRow.RPM;
            if isempty(rpm) || ~isnumeric(rpm) || isnan(rpm)
                rpm = 1400;
                warning('Invalid RPM in table; defaulting to 1400');
            end

            nTakes = numel(this.PressureInfo);
            raw_traces = cell(nTakes, 1);
            filt_traces = cell(nTakes, 1);
            time_traces = cell(nTakes, 1);

            for jj = 1:nTakes
                dataRaw = this.PressureInfo(jj).DataRaw;
                requestedChannel = this.PressureChannel;
                channelName = requestedChannel;

                % If caller passed Dev0_Ai5, also try base name Ai5
                baseChannel = channelName;
                if startsWith(channelName, 'Dev0_') && numel(channelName) > 5
                    baseChannel = channelName(6:end);
                end

                % Auto-discover channel if the requested one is missing
                analogFields = [];
                if isfield(dataRaw, 'Analog_channels')
                    analogFields = fieldnames(dataRaw.Analog_channels);
                end

                raw = [];
                fs = [];

                % Try exact name
                if isfield(dataRaw.Analog_channels, channelName)
                    raw = double(dataRaw.Analog_channels.(channelName).data) * 20;
                    try
                        fs = 1 ./ dataRaw.Analog_channels.(channelName).Props.wf_increment;
                    catch
                    end
                end

                % Try base name if we were given Dev0_* but exact failed
                if isempty(raw) && ~strcmp(baseChannel, channelName) && isfield(dataRaw.Analog_channels, baseChannel)
                    raw = double(dataRaw.Analog_channels.(baseChannel).data) * 20;
                    try
                        fs = 1 ./ dataRaw.Analog_channels.(baseChannel).Props.wf_increment;
                    catch
                    end
                    channelName = baseChannel;
                end

                % Try Dev0_ prefix if not found
                if isempty(raw)
                    prefixed = ['Dev0_' baseChannel];
                    if isfield(dataRaw.Analog_channels, prefixed)
                        raw = double(dataRaw.Analog_channels.(prefixed).data) * 20;
                        try
                            fs = 1 ./ dataRaw.Analog_channels.(prefixed).Props.wf_increment;
                        catch
                        end
                        channelName = prefixed;  % remember which one we used
                    end
                end

                % Fallback: pick first analog channel if nothing yet
                if isempty(raw) && ~isempty(analogFields)
                    fallbackChan = analogFields{1};
                    raw = double(dataRaw.Analog_channels.(fallbackChan).data) * 20;
                    try
                        fs = 1 ./ dataRaw.Analog_channels.(fallbackChan).Props.wf_increment;
                    catch
                    end
                    warning('Set %s Take %02d: pressure channel %s not found; using %s instead', ...
                        mat2str(this.SetNumber), jj-1, requestedChannel, fallbackChan);
                    channelName = fallbackChan;
                end

                % Final guard: if still empty, skip take but set Fs for downstream code
                if isempty(raw)
                    warning('Set %s Take %02d: pressure channel %s not found', mat2str(this.SetNumber), jj-1, requestedChannel);
                    this.PressureInfo(jj).Fs = 200000; % fallback to keep TDC running
                    continue
                end

                if isempty(fs)
                    fs = 200000; % fallback
                end

                t = (0:numel(raw)-1).' ./ fs;
                filt = this.PressureFilter(raw, this.Filter_CutoffFreq, rpm);

                raw_traces{jj} = raw;
                filt_traces{jj} = filt;
                time_traces{jj} = t;

                % Persist into PressureInfo for downstream use
                this.PressureInfo(jj).Raw = raw;
                this.PressureInfo(jj).Filtered = filt;
                this.PressureInfo(jj).Time = t;
                this.PressureInfo(jj).RPM = rpm;
                this.PressureInfo(jj).Fs = fs;
            end
        end

        %% Main Processing Pipeline
        function this = ProcessAll(this, tableRow)
            % ProcessAll - Main processing orchestrator
            % Executes complete analysis pipeline in proper sequence
            %
            % Input:
            %   tableRow - Current row from T_DataAll with metadata (RPM, etc.)
            %
            % Processing Steps:
            %   1. ProcessPressureTraces - Filter pressure signals
            %   2. FindTDC               - Detect TDC from encoder
            %   3. ExtractCycles         - Extract CA-aligned cycles
            %   4. CalculateIMEP         - Compute thermodynamic work
            %   5. CalculateHRR          - Compute combustion phasing
            
            % Step 1: Process pressure traces (filter, extract raw/time)
            this.ProcessPressureTraces(tableRow);
            
            % Step 2: Find TDC and compute crank angle
            this.FindTDC(tableRow);
            
            % Step 3: Extract individual cycles and interpolate to CA basis
            this.ExtractCycles(tableRow);
            
            % Step 4: Calculate IMEP and statistics
            this.CalculateIMEP(tableRow);
            
            % Step 5: Calculate HRR and combustion phasing
            this.CalculateHRR(tableRow);
        end
        
        %% Cycle Analysis Methods
        function this = FindTDC(this, tableRow)
            % FindTDC - Detect TDC from Z-pulse and encoder signals
            % Uses falling edge of Z-pulse as TDC marker with hardware offset
            %
            % Input:
            %   tableRow - Current row from T_DataAll (contains RPM)
            %
            % Algorithm:
            %   1. Detect Z-pulse falling edges (TDC markers)
            %   2. Apply hardware-specific phase shift (TDC_shift property)
            %   3. Compute CA from 0.2° encoder pulses
            %   4. Smooth and align CA to TDC positions
            
            rpm = tableRow.RPM;
            if isempty(rpm) || ~isnumeric(rpm) || isnan(rpm)
                rpm = 1400;
                warning('Invalid RPM in table; defaulting to 1400 for TDC calculation');
            end
            
            nTakes = numel(this.PressureInfo);
            
            for jj = 1:nTakes
                dataRaw = this.PressureInfo(jj).DataRaw;
                
                if isempty(dataRaw)
                    continue;
                end
                
                try
                    % Z-pulse rising edge (0->1) for raw TDCs
                    zpulseField = this.EncoderZPulseChannel;
                    Zpulse = [];
                    try
                        Zpulse = double(dataRaw.Digital_channels.(zpulseField).data);
                    catch
                        Zpulse = double(dataRaw.Digital_channels.(['Untitled_' zpulseField(end)]).data);
                    end
                    Zpulse = Zpulse(:);

                    Zprev = [1; Zpulse(1:end-1)];
                    isFallZ = (Zpulse == 0) & (Zprev == 1); % falling edge = TDC
                    isRiseZ = (Zpulse == 1) & (Zprev == 0); % rising edge = cycle boundary
                    TDC_raw = find(isFallZ);
                    Z_rise_raw = find(isRiseZ);
                    if numel(TDC_raw) < 2 || numel(Z_rise_raw) < 2
                        warning('Set %s Take %02d: Not enough Z-pulses detected', mat2str(this.SetNumber), jj-1);
                        continue;
                    end

                    % Apply phase shift only for alignment (not for cycle completeness)
                    TDC_diff = mean(diff(TDC_raw));
                    shift_samples = round((this.TDC_shift) / 720 * TDC_diff);
                    TDC_shifted = TDC_raw - shift_samples;
                    Z_rise_shifted = Z_rise_raw - shift_samples;

                    % Keep paired raw/shifted indices that remain in bounds after shift
                    in_bounds_fall = TDC_shifted >= 1 & TDC_shifted <= length(Zpulse);
                    TDC_raw = TDC_raw(in_bounds_fall);
                    TDC_shifted = TDC_shifted(in_bounds_fall);
                    in_bounds_rise = Z_rise_shifted >= 1 & Z_rise_shifted <= length(Zpulse);
                    Z_rise_raw = Z_rise_raw(in_bounds_rise);
                    Z_rise_shifted = Z_rise_shifted(in_bounds_rise);
                    if numel(TDC_raw) < 2 || numel(Z_rise_raw) < 2
                        warning('Set %s Take %02d: Not enough in-bound Z-pulses after shift', mat2str(this.SetNumber), jj-1);
                        continue;
                    end

                    % Encoder pulses rising edges -> CA
                    encField = this.EncoderDegreeChannel;
                    PulsTrace = [];
                    try
                        PulsTrace = double(dataRaw.Digital_channels.(encField).data);
                    catch
                        PulsTrace = double(dataRaw.Digital_channels.(['Untitled_' encField(end)]).data);
                    end
                    PulsTrace = PulsTrace(:);
                    CA_Enc = cumsum(double(PulsTrace > 0 & [1; PulsTrace(1:end-1)] == 0)) * 0.2;

                    % Smooth CA
                    fs_local = this.PressureInfo(jj).Fs;
                    if isempty(fs_local) || fs_local <= 0
                        fs_local = 200000;  % fallback sampling rate
                    end
                    wf_increment = 1 / fs_local;
                    bLen = round((60/rpm/1800) / wf_increment * 5);
                    bLen = max(bLen,1);
                    b = ones(bLen,1)/bLen;
                    CA_Enc = filtfilt(b, 1, CA_Enc);

                    % Normalize CA so that TDC positions are at CA=0 (not at wrap points)
                    % First unwrap to continuous angle
                    CA_Enc_unwrapped = unwrap(CA_Enc * pi/180) * 180/pi;
                    
                    % Calculate offset so TDC_curr positions are at multiples of 720
                    CA_offset = mean(CA_Enc_unwrapped(TDC_shifted) - (1:length(TDC_shifted))' * 720);
                    CA_Enc_aligned = CA_Enc_unwrapped - CA_offset;
                    
                    % Wrap back to -360 to +360 range
                    CA_Enc = mod(CA_Enc_aligned + 360, 720) - 360;

                    % Store both raw and shifted; shifted are primary cycle boundaries and CA=0
                    this.PressureInfo(jj).TDC_indices = TDC_shifted;
                    this.PressureInfo(jj).TDC_indices_shifted = TDC_shifted;
                    this.PressureInfo(jj).TDC_indices_raw = TDC_raw;
                    this.PressureInfo(jj).Z_rise_indices = Z_rise_raw;
                    this.PressureInfo(jj).Z_rise_indices_shifted = Z_rise_shifted;
                    this.PressureInfo(jj).CA_Enc = CA_Enc;
                    this.PressureInfo(jj).TDC_from_Z = TDC_shifted;
                    this.PressureInfo(jj).TDC_from_Z_raw = TDC_raw;

                catch ME
                    warning('Set %s Take %02d: TDC calculation failed - %s', ...
                        mat2str(this.SetNumber), jj-1, ME.message);
                end
            end
        end
        
        function this = ExtractCycles(this, tableRow)
            % ExtractCycles - Extract individual cycles aligned to crank angle
            % Centers each cycle on TDC (CA=0), extracts ±360° window
            %
            % Input:
            %   tableRow - Current row from T_DataAll
            %
            % Processing:
            %   1. Determine samples_per_degree from TDC spacing
            %   2. Extract ±360° pressure window around each TDC
            %   3. Build CA grid (encoder-based or uniform fallback)
            %   4. Interpolate pressure to standard CA grid
            %   5. Align cycles for ensemble averaging
            
            nTakes = numel(this.PressureInfo);
            
            for jj = 1:nTakes
                if ~isfield(this.PressureInfo(jj), 'Filtered') || isempty(this.PressureInfo(jj).Filtered) || ...
                   ~isfield(this.PressureInfo(jj), 'TDC_indices_shifted') || isempty(this.PressureInfo(jj).TDC_indices_shifted)
                    continue;
                end
                
                try
                    P_filt = this.PressureInfo(jj).Filtered;
                    TDC_idx_shifted = this.PressureInfo(jj).TDC_indices_shifted;
                    Z_rise_shifted = [];
                    if isfield(this.PressureInfo(jj), 'Z_rise_indices_shifted')
                        Z_rise_shifted = this.PressureInfo(jj).Z_rise_indices_shifted;
                    end
                    
                    % Initialize CA warnings counter for this take
                    n_CA_warnings = 0;
                    
                    % Estimate samples_per_degree from boundary spacing
                    % CRITICAL: Use TDC_idx_shifted (FALLING edge boundaries) for cycle extraction
                    % Z_rise_shifted is more stable but NOT used for extraction windows
                    if ~isempty(Z_rise_shifted) && numel(Z_rise_shifted) > 1
                        samples_per_cycle_nom = mean(diff(Z_rise_shifted));
                        samples_per_cycle_std = std(diff(Z_rise_shifted));
                    else
                        samples_per_cycle_nom = mean(diff(TDC_idx_shifted));
                        samples_per_cycle_std = std(diff(TDC_idx_shifted));
                    end
                    samples_per_degree = samples_per_cycle_nom / 720;
                    
                    % === CRITICAL FIX: Check TDC boundary variance (NOT Z-rise variance) ===
                    % TDC boundaries define cycle extraction windows, so their stability matters
                    tdc_spacing_std = std(diff(TDC_idx_shifted));
                    tdc_spacing_mean = mean(diff(TDC_idx_shifted));
                    tdc_boundary_variance = tdc_spacing_std / tdc_spacing_mean;
                    
                    % Store TDC boundary variance for CA encoder fallback logic
                    this.PressureInfo(jj).BoundarySpacingVariance = tdc_boundary_variance;
                    
                    if tdc_boundary_variance > 0.05  % >5% variance = severe problem
                        warning('Set %s Take %02d: High TDC boundary spacing variance (%.1f%%), cycle extraction may be unreliable', ...
                            mat2str(this.SetNumber), jj-1, tdc_boundary_variance * 100);
                    end
                    
                    % Extract valid cycles around each TDC (skip edges lacking full ±360)
                    nCycles = numel(TDC_idx_shifted) - 2;
                    if nCycles < 1
                        continue;
                    end
                    
                    % Storage for individual cycles interpolated to standard CA grid
                    P_cycles = zeros(length(this.CA), nCycles);
                    CA_start = nan(nCycles,1);
                    CA_end = nan(nCycles,1);
                    idx_zero_list = nan(nCycles,1);
                    boundary_start = nan(nCycles,1);
                    boundary_end = nan(nCycles,1);


                    for cc = 2:(numel(TDC_idx_shifted)-1)
                        idx_tdc = TDC_idx_shifted(cc);
                        CA_cycle = [];
                        
                        % Calculate indices for -360 to +360 around this TDC
                        idx_start = round(idx_tdc - 360 * samples_per_degree);
                        idx_end = round(idx_tdc + 360 * samples_per_degree);
                        
                        % Require indices to be in bounds
                        if idx_start < 1 || idx_end > length(P_filt)
                            continue;
                        end
                        
                        % Extract pressure for this cycle
                        P_cycle = P_filt(idx_start:idx_end);
                        
                        % Build CA grid from encoder if available; fall back to uniform
                        use_encoder_ca = isfield(this.PressureInfo(jj), 'CA_Enc') && ...
                                         ~isempty(this.PressureInfo(jj).CA_Enc) && ...
                                         numel(this.PressureInfo(jj).CA_Enc) >= idx_end;
                        
                        % === FIX: Detect bad encoder CA and fallback to uniform ===
                        % If boundary spacing variance > 0.5%, encoder CA is unreliable
                        % (Takes with high TDC variance show broken CA mapping)
                        force_uniform_ca = false;
                        if isfield(this.PressureInfo(jj), 'BoundarySpacingVariance')
                            if this.PressureInfo(jj).BoundarySpacingVariance > 0.005
                                force_uniform_ca = true;
                            end
                        end
                        
                        if use_encoder_ca && ~force_uniform_ca
                            CA_cycle_raw = this.PressureInfo(jj).CA_Enc(idx_start:idx_end);
                            CA_cycle = CA_cycle_raw - this.PressureInfo(jj).CA_Enc(idx_tdc); % zero at TDC

                            % === FIX CA WRAP-AROUND WITHIN CYCLE (Task B.4) ===
                            % Detect and unwrap CA discontinuities within a single cycle
                            % Expected: CA progresses smoothly from ~-360 → 0 (TDC) → +360
                            % Problem: Encoder wraps at ±360, causing jumps like +358 → -359
                            
                            dCA = diff(CA_cycle);
                            big_neg_jumps = find(dCA < -300);  % e.g., +358 → -359 (jump ≈ -717)
                            big_pos_jumps = find(dCA > 300);   % e.g., -358 → +359 (jump ≈ +717)
                            
                            % Unwrap: add 720° after negative jumps, subtract 720° after positive jumps
                            CA_cycle_unwrapped = CA_cycle;
                            for kk = 1:numel(big_neg_jumps)
                                idx = big_neg_jumps(kk);
                                CA_cycle_unwrapped((idx+1):end) = CA_cycle_unwrapped((idx+1):end) + 720;
                            end
                            for kk = 1:numel(big_pos_jumps)
                                idx = big_pos_jumps(kk);
                                CA_cycle_unwrapped((idx+1):end) = CA_cycle_unwrapped((idx+1):end) - 720;
                            end
                            
                            CA_cycle = CA_cycle_unwrapped;
                            % === END CA UNWRAP ===

                            % === CA INTEGRITY DEBUG & CHECK (Task A) ===
                            n_original = numel(CA_cycle);
                            [CA_cycle_unique, uniq_idx] = unique(CA_cycle, 'stable');
                            P_cycle_unique = P_cycle(uniq_idx);
                            n_unique = numel(CA_cycle_unique);
                            n_duplicates = n_original - n_unique;
                            
                            % Detailed diagnostic for first 3 cycles (after unwrap)
                            if cc <= 4 && jj == 1
                                fprintf('\n=== DEBUG Cycle %d ===\n', cc-1);
                                fprintf('  CA_raw range: [%.2f, %.2f]\n', min(CA_cycle_raw), max(CA_cycle_raw));
                                fprintf('  CA_zeroed+unwrapped range: [%.2f, %.2f]\n', min(CA_cycle), max(CA_cycle));
                                fprintf('  n_samples: %d, n_unique: %d, duplicates: %d\n', n_original, n_unique, n_duplicates);
                                fprintf('  Negative jumps corrected: %d, Positive jumps corrected: %d\n', ...
                                    numel(big_neg_jumps), numel(big_pos_jumps));
                                
                                % Check monotonicity after unique + unwrap
                                dCA_unique = diff(CA_cycle_unique);
                                neg_steps = find(dCA_unique <= 0);
                                if ~isempty(neg_steps)
                                    fprintf('  ✗ Still non-monotonic after unwrap: %d negative steps\n', numel(neg_steps));
                                else
                                    fprintf('  ✓ Monotonic after unwrap\n');
                                end
                            end
                            
                            if n_duplicates > 0
                                dup_pct = 100 * n_duplicates / n_original;
                                if dup_pct > 1.0
                                    warning('Set %s Take %02d Cycle %d: %d CA duplicates (%.1f%%), encoder issue', ...
                                        mat2str(this.SetNumber), jj-1, cc-1, n_duplicates, dup_pct);
                                end
                            end
                            
                            % Use the unique, unwrapped CA/Pressure vectors
                            CA_cycle = CA_cycle_unique;
                            P_cycle = P_cycle_unique;
                            
                            % Check monotonicity
                            dCA_check = diff(CA_cycle);
                            if any(dCA_check <= 0)
                                % Cycle is broken - reject it
                                warning('Set %s Take %02d Cycle %d: Non-monotonic CA after unwrap - REJECTING CYCLE', ...
                                    mat2str(this.SetNumber), jj-1, cc-1);
                                n_CA_warnings = n_CA_warnings + 1;  % Increment counter
                                continue;  % Skip this cycle
                            end
                            
                            % === NEW: Detect unrealistic linear pressure profiles ===
                            % Physical principle: Real cycles show curved compression/expansion
                            % Linear profile indicates wrong boundary or CA assignment
                            if numel(P_cycle) > 10
                                % Fit line to pressure vs CA and check R²
                                CA_fit = CA_cycle(:);
                                P_fit = P_cycle(:);
                                % Simple linear fit
                                p = polyfit(CA_fit, P_fit, 1);
                                P_fitted = polyval(p, CA_fit);
                                ss_res = sum((P_fit - P_fitted).^2);
                                ss_tot = sum((P_fit - mean(P_fit)).^2);
                                R_squared = 1 - (ss_res / ss_tot);
                                
                                % If R² > 0.95, cycle is very linear (likely wrong extraction)
                                if R_squared > 0.95
                                    warning('Set %s Take %02d Cycle %d: Unrealistic linear pressure profile (R²=%.3f) - CHECK BOUNDARY DETECTION', ...
                                        mat2str(this.SetNumber), jj-1, cc-1, R_squared);
                                end
                            end
                            % === END LINEAR PROFILE CHECK ===
                        end

                        % Fallback: if CA_cycle wasn't created or has wrong length, use uniform grid
                        if ~exist('CA_cycle','var') || isempty(CA_cycle) || numel(CA_cycle) ~= numel(P_cycle)
                            n_samples = idx_end - idx_start + 1;
                            CA_cycle = linspace(-360, 360, n_samples)';
                        end

                        % Interpolate to standard CA grid (this.CA: -360 to +360)
                        P_interp = interp1(CA_cycle, P_cycle, this.CA, 'linear', 'extrap');
                        if numel(P_interp) ~= numel(this.CA)
                            warning('Set %s Take %02d Cycle %d: interp length %d != CA %d, skipping cycle', ...
                                mat2str(this.SetNumber), jj-1, cc-1, numel(P_interp), numel(this.CA));
                            continue;
                        end
                        
                        P_cycles(:, cc-1) = P_interp;
                        CA_start(cc-1) = CA_cycle(1);
                        CA_end(cc-1) = CA_cycle(end);
                        idx_zero_list(cc-1) = idx_tdc;
                        boundary_start(cc-1) = idx_start;
                        boundary_end(cc-1) = idx_end;
                    end
                    
                    % Align cycles to TDC before averaging
                    % Each cycle has a different TDC position - we need to shift them so TDC aligns
                    idx_tdc_global = find(abs(this.CA) < 0.05, 1);  % TDC should be around index 3611
                    if isempty(idx_tdc_global)
                        idx_tdc_global = find(this.CA > 0, 1);
                    end
                    
                    P_cycles_aligned = P_cycles;  % Start with original
                    for cc = 1:size(P_cycles, 2)
                        idx_tdc_cycle = idx_zero_list(cc);
                        shift_amount = idx_tdc_global - idx_tdc_cycle;
                        
                        if shift_amount ~= 0 && abs(shift_amount) < 500  % Only shift if reasonable
                            P_shifted = circshift(P_cycles(:, cc), shift_amount);
                            P_cycles_aligned(:, cc) = P_shifted;
                        end
                    end
                    
                    % Store results
                    this.PressureInfo(jj).P_cycles = P_cycles;  % Individual cycles (original)
                    this.PressureInfo(jj).P_cycles_aligned = P_cycles_aligned;  % Aligned cycles
                    this.PressureInfo(jj).P_mean = mean(P_cycles_aligned, 2);  % Mean of aligned cycles
                    this.PressureInfo(jj).Cycle_CA_start = CA_start;
                    this.PressureInfo(jj).Cycle_CA_end = CA_end;
                    this.PressureInfo(jj).Cycle_zero_sample = idx_zero_list;
                    this.PressureInfo(jj).Cycle_bound_start = boundary_start;
                    this.PressureInfo(jj).Cycle_bound_end = boundary_end;
                    this.PressureInfo(jj).CA_warnings = n_CA_warnings;  % Store warning count
                    
                catch ME
                    warning('Set %s Take %02d: Cycle extraction failed - %s', ...
                        mat2str(this.SetNumber), jj-1, ME.message);
                end
            end
        end
        
        %% Thermodynamic Analysis Methods
        function this = CalculateIMEP(this, tableRow)
            % CalculateIMEP - Calculate indicated mean effective pressure and statistics
            %
            % Input:
            %   tableRow - Current row from T_DataAll
            %
            % Outputs (stored in PressureInfo):
            %   IMEP_cycles - IMEP for each individual cycle [bar]
            %   IMEP        - Mean IMEP across all cycles [bar]
            %   CoV_IMEP    - Coefficient of variation (std/mean)
            %   P_max       - Peak pressure [bar]
            %   PRR_max     - Maximum pressure rise rate [bar/CAD]
            
            nTakes = numel(this.PressureInfo);
            
            for jj = 1:nTakes
                % default empty outputs
                this.PressureInfo(jj).IMEP_cycles = [];
                this.PressureInfo(jj).IMEP = [];
                this.PressureInfo(jj).CoV_IMEP = [];
                this.PressureInfo(jj).P_max = [];
                this.PressureInfo(jj).PRR_max = [];

                if ~isfield(this.PressureInfo(jj), 'P_cycles') || isempty(this.PressureInfo(jj).P_cycles)
                    continue;
                end
                
                try
                    P_cycles = this.PressureInfo(jj).P_cycles;
                    if isfield(this.PressureInfo(jj), 'P_mean_corrected') && ~isempty(this.PressureInfo(jj).P_mean_corrected)
                        P_mean = this.PressureInfo(jj).P_mean_corrected;
                    else
                        P_mean = this.PressureInfo(jj).P_mean;
                    end

                    if size(P_cycles,1) ~= numel(this.CA) || numel(P_mean) ~= numel(this.Vol)
                        warning('Set %s Take %02d: P length mismatch (P_cycles rows %d, CA %d, Vol %d). Skipping IMEP.', ...
                            mat2str(this.SetNumber), jj-1, size(P_cycles,1), numel(this.CA), numel(this.Vol));
                        continue;
                    end

                    % Apply intake offset so IMEP uses absolute-like pressure (~0.9 bar intake)
                    CA_range_intake = find(this.CA > -210, 1) : find(this.CA > -180, 1);
                    if isempty(CA_range_intake)
                        warning('Set %s Take %02d: intake CA range missing, skipping IMEP correction.', mat2str(this.SetNumber), jj-1);
                        continue;
                    end
                    P_offset = mean(P_mean(CA_range_intake)); % TDMS_Processing.m:253
                    P_mean_corr = P_mean - P_offset + 0.9;
                    P_cycles_corr = P_cycles - P_offset + 0.9;
                    
                    % Calculate IMEP for each individual cycle
                    nCycles = size(P_cycles, 2);
                    if nCycles == 0
                        continue;
                    end
                    IMEP_cycles = zeros(nCycles, 1);
                    
                    dV = diff(this.Vol);  % column vector (because Vol is column)
                    for cc = 1:nCycles
                        % IMEP = ∫P·dV / DV (TDMS_Processing.m:233)
                        % Subtract 1.0 bar offset: The pressure transducer is vented to atmosphere
                        % (reads gauge pressure), so raw readings have a negative baseline at intake
                        % (~-14.7 bar). When integrated over the cycle, this creates a systematic
                        % offset that must be removed to get the true indicated work.
                        % The 1.0 bar = 0.9 bar (target absolute intake) + 0.1 bar (polytropic offset)
                        IMEP_cycles(cc) = sum(P_cycles(1:end-1, cc) .* dV) / this.DV - 1.0;
                    end
                    
                    % Mean IMEP
                    IMEP_mean = mean(IMEP_cycles, 'omitnan');
                    
                    % Coefficient of Variation (CoV)
                    CoV_IMEP = std(IMEP_cycles, 'omitnan') / IMEP_mean;
                    
                    % Peak pressure (relative to minimum + motoring pressure ~0.8 bar)
                    P_max = max(P_mean_corr) - min(P_mean_corr) + 0.8;
                    
                    % Peak Pressure Rise Rate (PRR)
                    dPdCA = diff(P_mean_corr);
                    PRR_max = max(dPdCA) * 10;  % Factor of 10 for units conversion
                    
                    % Store results
                    this.PressureInfo(jj).IMEP_cycles = IMEP_cycles;
                    this.PressureInfo(jj).IMEP = IMEP_mean;
                    this.PressureInfo(jj).CoV_IMEP = CoV_IMEP;
                    this.PressureInfo(jj).P_max = P_max;
                    this.PressureInfo(jj).PRR_max = PRR_max;
                    this.PressureInfo(jj).P_mean_corrected = P_mean_corr;
                    
                catch ME
                    warning('Set %s Take %02d: IMEP calculation failed - %s', ...
                        mat2str(this.SetNumber), jj-1, ME.message);
                end
            end
        end
        
        function this = CalculateHRR(this, tableRow)
            % CalculateHRR - Calculate heat release rate and combustion phasing
            % Based on Process_thermodynamics function from TDMS_Processing.m
            %
            % Input:
            %   tableRow - Current row from T_DataAll
            %
            % Outputs (stored in PressureInfo):
            %   aHRR     - Apparent heat release rate [J/CAD]
            %   cumHRR   - Cumulative heat release [J]
            %   CA10     - 10% burn location [CAD]
            %   CA50     - 50% burn location (combustion phasing) [CAD]
            %   CA90     - 90% burn location [CAD]
            %   BurnDur  - Burn duration (CA90-CA10) [CAD]
            %   IgnDelay - Ignition delay (CA10-SOI) [CAD]
            %
            % Algorithm:
            %   1. Fit motoring pressure to determine heat loss
            %   2. Apply pressure correction (intake offset)
            %   3. Calculate aHRR from dP/dCA and dV/dCA
            %   4. Integrate to get cumulative HRR
            %   5. Determine CA10/50/90 from cumulative HRR percentiles
            %   6. Extract SOI from injection signal for ignition delay
            
            nTakes = numel(this.PressureInfo);
            
            for jj = 1:nTakes
                % default empty outputs
                this.PressureInfo(jj).P_corrected_motoring = [];
                this.PressureInfo(jj).TDC_P = [];
                this.PressureInfo(jj).aHRR = [];
                this.PressureInfo(jj).cumHRR = [];
                this.PressureInfo(jj).SumHRR = [];
                this.PressureInfo(jj).CA10 = [];
                this.PressureInfo(jj).CA50 = [];
                this.PressureInfo(jj).CA90 = [];
                this.PressureInfo(jj).BurnDur = [];
                this.PressureInfo(jj).CA10_50 = [];
                this.PressureInfo(jj).CA50_90 = [];
                this.PressureInfo(jj).CA10_90 = [];
                this.PressureInfo(jj).IgnDelay = [];
                this.PressureInfo(jj).SOI_CA = [];

                if ~isfield(this.PressureInfo(jj), 'P_mean') || isempty(this.PressureInfo(jj).P_mean)
                    continue;
                end
                
                try
                    P_mean = this.PressureInfo(jj).P_mean;

                    lenP = numel(P_mean);
                    lenCA = numel(this.CA);
                    lenVol = numel(this.Vol);

                    if lenP ~= lenCA || lenP ~= lenVol
                        warning('Set %s Take %02d: P_mean length mismatch (P_mean %d, CA %d, Vol %d). Skipping HRR.', ...
                            mat2str(this.SetNumber), jj-1, lenP, lenCA, lenVol);
                        continue;
                    end
                    
                    % Constants
                    CV = this.DV / (this.CR - 1);  % clearance volume [m^3]
                    Kappa = 1.40;  % specific heat ratio
                    
                    % Step 1: Fit motoring pressure (P0) in compression range
                    % Range: -143 to -110 CAD
                    CA_range_fit = find(this.CA > -143, 1) : find(this.CA > -110, 1);
                    if isempty(CA_range_fit)
                        error('CA_range_fit empty');
                    end
                    
                    % Objective: fit P = P0(1) + P0(2) * (V_TDC/V)^Kappa
                    fun = @(P0) sum((P_mean(CA_range_fit) - P0(1) - P0(2) * ((CV + this.DV) ./ this.Vol(CA_range_fit)).^Kappa).^2);
                    P0_guess = [-1, 1];
                    P0 = fminsearch(fun, P0_guess);
                    
                    % Determine offset from intake pressure (-210 to -180 CAD) (TDMS_Processing.m:252)
                    CA_range_intake = find(this.CA > -210, 1) : find(this.CA > -180, 1);
                    if isempty(CA_range_intake)
                        error('CA_range_intake empty');
                    end
                    % Legacy formula: P0(1) = mean(P_intake) - 0.1
                    % Then: P_corrected = P - P0(1) + 0.9 = P - mean + 1.0
                    P_offset = mean(P_mean(CA_range_intake)) - 0.1;  % TDMS_Processing.m:252
                    P_corrected = P_mean - P_offset + 0.9;
                    
                    % Apply Savitzky-Golay filter to Pmean (TDMS_Processing.m:326)
                    P_corrected = sgolayfilt(P_corrected, 7, 91);
                    
                    P_cycles_corr = [];
                    if isfield(this.PressureInfo(jj), 'P_cycles') && ~isempty(this.PressureInfo(jj).P_cycles)
                        % Apply offset correction to all per-cycle pressure traces
                        P_cycles_corr = this.PressureInfo(jj).P_cycles - P_offset + 0.9;
                        
                        % Apply Savitzky-Golay filter to each cycle (TDMS_Processing.m:326)
                        nCycles_temp = size(P_cycles_corr, 2);
                        for cc_filt = 1:nCycles_temp
                            P_cycles_corr(:, cc_filt) = sgolayfilt(P_cycles_corr(:, cc_filt), 7, 91);
                        end
                    end
                    
                    % TDC pressure
                    [~, idx_TDC] = min(abs(this.CA));
                    TDC_P = P_corrected(idx_TDC);

                    % Step 2: Calculate apparent HRR with Kappa = 1.35
                    Kappa = 1.35;
                    dV = diff(this.Vol);
                    dCA = this.CA(2) - this.CA(1);
                    
                    % Zero cumHRR reference point
                    idx_zero = find(this.CA > -50, 1);
                    if isempty(idx_zero)
                        error('Cannot find -50 CA in this.CA array');
                    end
                    
                    % === STEP 3: PER-CYCLE CA METRICS (Validation Task 4) ===
                    % Physical principle: Cycle-to-cycle variability requires per-cycle analysis
                    % Compute aHRR, cumHRR, CA10/CA50/CA90 per cycle FIRST, then average
                    
                    if ~isempty(P_cycles_corr)
                        nCycles = size(P_cycles_corr, 2);
                        CA10_cycles = nan(nCycles, 1);
                        CA50_cycles = nan(nCycles, 1);
                        CA90_cycles = nan(nCycles, 1);
                        BurnDur_cycles = nan(nCycles, 1);
                        SumHRR_cycles = nan(nCycles, 1);
                        aHRR_cycles = zeros(lenP-1, nCycles);
                        cumHRR_cycles = zeros(lenP-1, nCycles);
                        
                        for cc = 1:nCycles
                            % Compute aHRR for this individual cycle
                            dP_c = diff(P_cycles_corr(:, cc));
                            aHRR_c = (Kappa/(Kappa-1)) * P_cycles_corr(1:end-1, cc) .* dV + ...
                                     (1/(Kappa-1)) * this.Vol(1:end-1) .* dP_c;
                            aHRR_c = aHRR_c * 1e5 / dCA;
                            
                            % Cumulative HRR for this cycle
                            cumHRR_c = cumsum(aHRR_c) * dCA;
                            cumHRR_c = cumHRR_c - cumHRR_c(idx_zero);  % Zero at -50 CAD
                            
                            % Compute SumHRR for this cycle
                            CA_range_min = find(this.CA > -50, 1) : find(this.CA > -30, 1);
                            CA_range_max = find(this.CA > 0, 1) : find(this.CA > 180, 1);
                            minCum_c = min(cumHRR_c(CA_range_min));
                            maxCum_c = max(cumHRR_c(CA_range_max));
                            SumHRR_c = maxCum_c - minCum_c;
                            SumHRR_cycles(cc) = SumHRR_c;
                            
                            % === COMBUSTION DETECTION (Validation Task B.5) ===
                            % Physical principle: Firing cycles have high peak pressure (combustion release)
                            % Motoring/misfire cycles have low peak pressure (compression only)
                            % Threshold: P_max > 45 bar for firing (empirical for this engine)
                            
                            P_max_cycle = max(P_cycles_corr(:, cc));
                            combustion_threshold_P = 45;  % bar
                            
                            if P_max_cycle > combustion_threshold_P
                                % Firing cycle: compute CA10/CA50/CA90
                                CA_range = find(this.CA > -50, 1) : find(this.CA > 150, 1);
                                
                                idx_10 = find(cumHRR_c(CA_range) > minCum_c + 0.1 * SumHRR_c, 1, 'first');
                                idx_50 = find(cumHRR_c(CA_range) > minCum_c + 0.5 * SumHRR_c, 1, 'first');
                                idx_90 = find(cumHRR_c(CA_range) > minCum_c + 0.9 * SumHRR_c, 1, 'first');
                                
                                if ~isempty(idx_10)
                                    CA10_cycles(cc) = this.CA(find(this.CA > -50, 1) + idx_10 - 1);
                                end
                                if ~isempty(idx_50)
                                    CA50_cycles(cc) = this.CA(find(this.CA > -50, 1) + idx_50 - 1);
                                end
                                if ~isempty(idx_90)
                                    CA90_cycles(cc) = this.CA(find(this.CA > -50, 1) + idx_90 - 1);
                                    BurnDur_cycles(cc) = CA90_cycles(cc) - CA10_cycles(cc);
                                end
                            else
                                % Non-firing cycle (motoring or misfire): set metrics to NaN
                                % Do not compute CA10/50/90 from near-zero HRR
                                CA10_cycles(cc) = NaN;
                                CA50_cycles(cc) = NaN;
                                CA90_cycles(cc) = NaN;
                                BurnDur_cycles(cc) = NaN;
                            end
                            
                            % Store per-cycle aHRR and cumHRR
                            aHRR_cycles(:, cc) = aHRR_c;
                            cumHRR_cycles(:, cc) = cumHRR_c;
                        end
                        
                        % === CREATE FIRING INDICATOR ARRAY ===
                        % 1 = firing cycle (P_max > 45 bar), 0 = misfired/motoring (P_max <= 45 bar)
                        firing_indicator = ~isnan(CA10_cycles);
                        this.PressureInfo(jj).Firing_indicator = firing_indicator;
                        
                        % Count firing vs non-firing cycles
                        n_firing = sum(firing_indicator);
                        n_total = nCycles;
                        firing_rate = 100 * n_firing / n_total;
                        
                        % Report combustion statistics
                        if firing_rate < 80
                            warning('Set %s Take %02d: Low firing rate %.1f%% (%d/%d cycles)', ...
                                mat2str(this.SetNumber), jj-1, firing_rate, n_firing, n_total);
                        end
                        
                        % Average across all cycles
                        CA10 = mean(CA10_cycles, 'omitnan');
                        CA50 = mean(CA50_cycles, 'omitnan');
                        CA90 = mean(CA90_cycles, 'omitnan');
                        BurnDur = mean(BurnDur_cycles, 'omitnan');
                        SumHRR = mean(SumHRR_cycles, 'omitnan');
                        
                        % Standard deviations for uncertainty reporting
                        CA10_std = std(CA10_cycles, 'omitnan');
                        CA50_std = std(CA50_cycles, 'omitnan');
                        CA90_std = std(CA90_cycles, 'omitnan');
                        
                        % Mean aHRR/cumHRR for plotting (ensemble average)
                        aHRR = mean(aHRR_cycles, 2);
                        cumHRR = mean(cumHRR_cycles, 2);
                        
                    else
                        % Fallback: no per-cycle data, use averaged pressure (legacy behavior)
                        warning('Set %s Take %02d: No P_cycles available, using averaged metrics', ...
                            mat2str(this.SetNumber), jj-1);
                        
                        % Compute from averaged pressure
                        dP = diff(P_corrected(:));
                        aHRR = (Kappa/(Kappa-1)) * P_corrected(1:end-1)' .* dV + ...
                               (1/(Kappa-1)) * this.Vol(1:end-1) .* dP;
                        aHRR = aHRR * 1e5 / dCA;
                        
                        cumHRR = cumsum(aHRR) * dCA;
                        cumHRR = cumHRR - cumHRR(idx_zero);
                        
                        % Compute CA metrics from averaged signals
                        CA_range_min = find(this.CA > -50, 1) : find(this.CA > -30, 1);
                        CA_range_max = find(this.CA > 0, 1) : find(this.CA > 180, 1);
                        minCum = min(cumHRR(CA_range_min));
                        maxCum = max(cumHRR(CA_range_max));
                        SumHRR = maxCum - minCum;
                        
                        % Check if combustion occurred (fallback path - use P_max as indicator)
                        P_max_Take = max(P_corrected);
                        combustion_threshold_P = 45;  % bar
                        if P_max_Take < combustion_threshold_P
                            warning('Set %s Take %02d: No combustion detected (P_max=%.1f bar, threshold=%.1f bar)', ...
                                mat2str(this.SetNumber), jj-1, P_max_Take, combustion_threshold_P);
                            CA10 = NaN; CA50 = NaN; CA90 = NaN; BurnDur = NaN;
                            CA10_std = NaN; CA50_std = NaN; CA90_std = NaN;
                        else
                            CA_range = find(this.CA > -50, 1) : find(this.CA > 150, 1);
                            idx_10 = find(cumHRR(CA_range) > minCum + 0.1 * SumHRR, 1, 'first');
                            idx_50 = find(cumHRR(CA_range) > minCum + 0.5 * SumHRR, 1, 'first');
                            idx_90 = find(cumHRR(CA_range) > minCum + 0.9 * SumHRR, 1, 'first');
                            
                            CA10 = this.CA(find(this.CA > -50, 1) + idx_10 - 1); 
                            CA50 = this.CA(find(this.CA > -50, 1) + idx_50 - 1);
                            CA90 = this.CA(find(this.CA > -50, 1) + idx_90 - 1);
                            BurnDur = CA90 - CA10;
                            
                            CA10_std = NaN; CA50_std = NaN; CA90_std = NaN;  % No per-cycle data
                        end
                        
                        aHRR_cycles = [];
                        cumHRR_cycles = [];
                        CA10_cycles = [];
                        CA50_cycles = [];
                        CA90_cycles = [];
                    end
                    % === END PER-CYCLE METRICS ===
                    
                    CA10_50 = CA50 - CA10;
                    CA50_90 = CA90 - CA50;
                    
                    % Ignition delay: CA10 - SOI (Start Of Injection)
                    % Extract SOI from injection signal (matches TDMS_Processing)
                    SOI_CA = NaN;  % Default if injection not found
                    try
                        dataRaw = this.PressureInfo(jj).DataRaw;
                        if ~isempty(dataRaw) && isfield(dataRaw, 'Digital_channels')
                            % Priority: Di3 (pilot injection) is earliest physical signal
                            InjectionSignal = [];
                            channelName = '';
                            
                            % Try channels in priority order: Di3 > Untitled_3 > Di5 > Untitled_5
                            if isfield(dataRaw.Digital_channels, this.PilotFuelTriggerChannel)
                                InjectionSignal = double(dataRaw.Digital_channels.(this.PilotFuelTriggerChannel).data);
                                channelName = this.PilotFuelTriggerChannel;
                            elseif isfield(dataRaw.Digital_channels, 'Untitled_3')
                                InjectionSignal = double(dataRaw.Digital_channels.Untitled_3.data);
                                channelName = 'Untitled_3';
                            else
                                % List available channels for debugging
                                availChans = fieldnames(dataRaw.Digital_channels);
                                warning('Set %d Take %02d: No injection channel found. Available: %s', ...
                                    this.SetNumber, jj-1, strjoin(availChans, ', '));
                            end
                            
                            if ~isempty(InjectionSignal)
                                % Find rising edges (injection start)
                                % Ensure column vector for concatenation
                                InjectionSignal = InjectionSignal(:);
                                diff_inj = [0; diff(InjectionSignal)];
                                inj_edges = find(diff_inj > 0.5);
                                
                                if ~isempty(inj_edges)
                                    % Use first injection edge in the measurement
                                    inj_sample_idx = inj_edges(1);
                                    
                                    % Map sample index to CA using CA_Enc
                                    if isfield(this.PressureInfo(jj), 'CA_Enc') && ...
                                       ~isempty(this.PressureInfo(jj).CA_Enc) && ...
                                       inj_sample_idx <= numel(this.PressureInfo(jj).CA_Enc)
                                        SOI_CA = this.PressureInfo(jj).CA_Enc(inj_sample_idx);
                                        % fprintf('  Set %d Take %02d: SOI detected at %.1f CAD using %s\n', ...
                                        %     this.SetNumber, jj-1, SOI_CA, channelName);
                                    end
                                end
                            end
                        end
                    catch ME
                        warning('Set %d Take %02d: Injection signal extraction failed: %s', ...
                            this.SetNumber, jj-1, ME.message);
                    end
                    
                    % Compute ignition delay as CA10 relative to SOI
                    if ~isnan(SOI_CA)
                        IgnDelay = CA10 - SOI_CA;
                    else
                        IgnDelay = CA10;  % Fallback: absolute CA10
                    end

                    % Store results (including per-cycle metrics)
                    this.PressureInfo(jj).P_corrected_motoring = P_corrected;
                    this.PressureInfo(jj).TDC_P = TDC_P;
                    this.PressureInfo(jj).aHRR = aHRR;
                    this.PressureInfo(jj).cumHRR = cumHRR;
                    this.PressureInfo(jj).aHRR_cycles = aHRR_cycles;
                    this.PressureInfo(jj).cumHRR_cycles = cumHRR_cycles;
                    this.PressureInfo(jj).SumHRR = SumHRR;
                    this.PressureInfo(jj).CA10 = CA10;
                    this.PressureInfo(jj).CA50 = CA50;
                    this.PressureInfo(jj).CA90 = CA90;
                    this.PressureInfo(jj).BurnDur = BurnDur;
                    this.PressureInfo(jj).CA10_50 = CA10_50;
                    this.PressureInfo(jj).CA50_90 = CA50_90;
                    this.PressureInfo(jj).CA10_90 = BurnDur;
                    this.PressureInfo(jj).IgnDelay = IgnDelay;
                    this.PressureInfo(jj).SOI_CA = SOI_CA;
                    
                    % Store per-cycle metrics (Validation Task 4)
                    if exist('CA10_cycles', 'var')
                        this.PressureInfo(jj).CA10_cycles = CA10_cycles;
                        this.PressureInfo(jj).CA50_cycles = CA50_cycles;
                        this.PressureInfo(jj).CA90_cycles = CA90_cycles;
                        this.PressureInfo(jj).SumHRR_cycles = SumHRR_cycles;
                        this.PressureInfo(jj).CA10_std = CA10_std;
                        this.PressureInfo(jj).CA50_std = CA50_std;
                        this.PressureInfo(jj).CA90_std = CA90_std;
                    end
                    
                catch ME
                    warning('Set %s Take %02d: HRR calculation failed - %s', ...
                        mat2str(this.SetNumber), jj-1, ME.message);
                end
            end
        end
        
        %% Take Management Methods
        function this = DiscardTake(this, takeNumber)
            % DiscardTake - Mark a Take as discarded (exclude from analysis)
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
        
        function this = RecoverTake(this, takeNumber)
            % RecoverTake - Unmark a Take (restore to valid status)
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
        
        %% Getter Methods
        function validTakes = GetValidTakes(this)
            % GetValidTakes - Get PressureInfo for non-discarded Takes only
            %
            % Returns:
            %   Structure array with only valid Takes
            
            validIdx = ~this.DiscardedTakes;
            validTakes = this.PressureInfo(validIdx);
        end
        
        function discardedTakes = GetDiscardedTakes(this)
            % GetDiscardedTakes - Get PressureInfo for discarded Takes only
            
            discardedIdx = this.DiscardedTakes;
            discardedTakes = this.PressureInfo(discardedIdx);
        end
        
        %% Display and Reporting Methods
        function DisplayStatus(this)
            % DisplayStatus - Print summary of Take validity status
            
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
        
        function S_PressureInfo = GetPressureInfo(this)
            % GetPressureInfo - Access complete PressureInfo struct array
            % Returns ALL Takes (including discarded)
            %
            % Note: Use GetValidTakes() to retrieve only non-discarded Takes
            
            S_PressureInfo = this.PressureInfo;
        end
    end
end
