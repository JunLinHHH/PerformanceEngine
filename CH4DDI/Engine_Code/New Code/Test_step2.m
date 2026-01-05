%% Diagnostic - Check Table Structure
clearvars; close all; clc;

Configs = Configuration();
Loaded = load(Configs.DataObjReadDir);
obj = Loaded.obj;
T_DataAll = obj.DataMatrix;

fprintf('=== TABLE DIAGNOSTICS ===\n\n');

% Check Set column type
fprintf('Set column class: %s\n', class(T_DataAll.Set));
fprintf('Set column size: %s\n', mat2str(size(T_DataAll.Set)));
fprintf('Is cell array: %s\n', mat2str(iscell(T_DataAll.Set)));
fprintf('Is numeric: %s\n', mat2str(isnumeric(T_DataAll.Set)));

% Show first few values
fprintf('\nFirst 3 Set values:\n');
for i = 1:min(3, height(T_DataAll))
    if iscell(T_DataAll.Set)
        fprintf('  Row %d: Set = %s (cell)\n', i, mat2str(T_DataAll.Set{i}));
    else
        fprintf('  Row %d: Set = %s\n', i, mat2str(T_DataAll.Set(i)));
    end
end

% Check File_TDMS_AllTakes column
fprintf('\nFile_TDMS_AllTakes column class: %s\n', class(T_DataAll.File_TDMS_AllTakes));

fprintf('\n=== END DIAGNOSTICS ===\n');

%% 
filePath = 'C:\Users\jun-y\OneDrive\桌面\20250909\20250909\_Set08_1400_30_90\_Set08_1400_30_90_Take02.tdms';
DataRaw  = TDMS_getStruct(filePath);

PressureRaw = double(DataRaw.Analog_channels.Dev0_Ai5.data)*20;
SamplingRate_Local = 200000;
SignalTimeScale = (0:length(PressureRaw)-1)/SamplingRate_Local; 

%% Properties
TDC_shift = -68;

%% Build pressure signal & time vector
SignalPressure = PressureRaw;

%% --- Method 1: Bandstop + Savitzky-Golay ---
SignalPressure_Method1 = bandstop(SignalPressure,[8800 11200],SamplingRate_Local);
SignalPressure_Method1 = bandstop(SignalPressure_Method1,[36000 42000],SamplingRate_Local);
SignalPressure_Method1 = sgolayfilt(SignalPressure_Method1, 3, 101);
SignalPressure_Method1 = bandstop(SignalPressure_Method1,[3400 3800],SamplingRate_Local);

%% --- Method 2: Butterworth only ---
f_cutoff = 10000;  
RPM = 1400;
fc_norm  = f_cutoff / (7200*(RPM/60/2));        
[b_bw, a_bw] = butter(2, fc_norm);           
SignalPressure_Method2 = filtfilt(b_bw, a_bw, SignalPressure);

%% --- Method 3: Combined (Bandstop + Savitzky-Golay + Butterworth) ---
SignalPressure_Method3 = bandstop(SignalPressure,[8800 11200],SamplingRate_Local);
SignalPressure_Method3 = bandstop(SignalPressure_Method3,[36000 42000],SamplingRate_Local);
SignalPressure_Method3 = bandstop(SignalPressure_Method3,[3400 3800],SamplingRate_Local);
SignalPressure_Method3 = sgolayfilt(SignalPressure_Method3, 3, 51);   % shorter SG window to reduce ripple
SignalPressure_Method3 = filtfilt(b_bw, a_bw, SignalPressure_Method3);

% For backward compatibility
SignalPressure_Filtered = SignalPressure_Method3;


%% Plot time-domain: raw vs filtered (full)
t = SignalTimeScale;

figure;
subplot(2,1,1);
plot(t, SignalPressure, 'k', 'LineWidth', 1.5); hold on;
plot(t, SignalPressure_Method1, 'r', 'LineWidth', 1);
plot(t, SignalPressure_Method2, 'b', 'LineWidth', 1);
plot(t, SignalPressure_Method3, 'm', 'LineWidth', 1);
xlabel('Time [s]');
ylabel('Pressure [bar]');
title('Pressure Signal Comparison (Full)');
legend('Raw','Method1 (Bandstop+SG)','Method2 (Butterworth)','Method3 (Combined)','Location','best');
grid on;

% 👀 Zoom region
t_min = 0.005;
t_max = 0.015;
idx = t >= t_min & t <= t_max;

subplot(2,1,2);
plot(t(idx), SignalPressure(idx), 'k', 'LineWidth', 1.5); hold on;
plot(t(idx), SignalPressure_Method1(idx), 'r', 'LineWidth', 1);
plot(t(idx), SignalPressure_Method2(idx), 'b', 'LineWidth', 1);
plot(t(idx), SignalPressure_Method3(idx), 'm', 'LineWidth', 1);
xlabel('Time [s]');
ylabel('Pressure [bar]');
title(sprintf('Zoomed Region (%.4f–%.4f s)', t_min, t_max));
legend('Raw','Method1','Method2','Method3','Location','best');
grid on;


%% Frequency-domain: PSD before vs after
[PSD_raw, f_psd]  = pwelch(SignalPressure,[],[],[],SamplingRate_Local);
[PSD_m1, ~]       = pwelch(SignalPressure_Method1,[],[],[],SamplingRate_Local);
[PSD_m2, ~]       = pwelch(SignalPressure_Method2,[],[],[],SamplingRate_Local);
[PSD_m3, ~]       = pwelch(SignalPressure_Method3,[],[],[],SamplingRate_Local);

figure;
loglog(f_psd, PSD_raw, 'k', 'LineWidth', 1.5); hold on;
loglog(f_psd, PSD_m1, 'r', 'LineWidth', 1);
loglog(f_psd, PSD_m2, 'b', 'LineWidth', 1);
loglog(f_psd, PSD_m3, 'm', 'LineWidth', 1);
xlabel('Frequency [Hz]');
ylabel('PSD');
title('Power Spectral Density Comparison');
legend('Raw','Method1 (Bandstop+SG)','Method2 (Butterworth)','Method3 (Combined)','Location','best');
grid on;
xlim([1 1e5]);


%% Noise estimate for each method
noise_est_m1 = SignalPressure - SignalPressure_Method1;
noise_est_m2 = SignalPressure - SignalPressure_Method2;
noise_est_m3 = SignalPressure - SignalPressure_Method3;

[PSD_noise_m1, f_noise] = pwelch(noise_est_m1,[],[],[],SamplingRate_Local);
[PSD_noise_m2, ~] = pwelch(noise_est_m2,[],[],[],SamplingRate_Local);
[PSD_noise_m3, ~] = pwelch(noise_est_m3,[],[],[],SamplingRate_Local);

figure;
loglog(f_psd, PSD_raw, 'k', 'LineWidth', 1.5); hold on;
loglog(f_noise, PSD_noise_m1, 'r', 'LineWidth', 1);
loglog(f_noise, PSD_noise_m2, 'b', 'LineWidth', 1);
loglog(f_noise, PSD_noise_m3, 'm', 'LineWidth', 1);
xlabel('Frequency [Hz]');
ylabel('PSD');
title('Noise Spectrum (Raw - Filtered)');
legend('Raw Signal','Noise M1','Noise M2','Noise M3','Location','best');
grid on;


%% Lag check for each method
sig_raw = SignalPressure - mean(SignalPressure);
sig_m1 = SignalPressure_Method1 - mean(SignalPressure_Method1);
sig_m2 = SignalPressure_Method2 - mean(SignalPressure_Method2);
sig_m3 = SignalPressure_Method3 - mean(SignalPressure_Method3);

[acor_m1, lag] = xcorr(sig_m1, sig_raw, 'coeff');
[~, I1] = max(acor_m1);
lag_samples_m1 = lag(I1);
lag_time_m1 = lag_samples_m1 / SamplingRate_Local;

[acor_m2, ~] = xcorr(sig_m2, sig_raw, 'coeff');
[~, I2] = max(acor_m2);
lag_samples_m2 = lag(I2);
lag_time_m2 = lag_samples_m2 / SamplingRate_Local;

[acor_m3, ~] = xcorr(sig_m3, sig_raw, 'coeff');
[~, I3] = max(acor_m3);
lag_samples_m3 = lag(I3);
lag_time_m3 = lag_samples_m3 / SamplingRate_Local;

fprintf('\n=== Filter Lag Analysis ===\n');
fprintf('Method 1 lag: %d samples (%.6f s)\n', lag_samples_m1, lag_time_m1);
fprintf('Method 2 lag: %d samples (%.6f s)\n', lag_samples_m2, lag_time_m2);
fprintf('Method 3 lag: %d samples (%.6f s)\n', lag_samples_m3, lag_time_m3);


%% dP/dt (unchanged)
dPdt = gradient(SignalPressure_Filtered, SignalTimeScale);
figure;
plot(dPdt);
grid on;

%% Diesel signal
DieselSignal = DataRaw.Digital_channels.Untitled_3.data*20;
MianFuelSignal = DataRaw.Digital_channels.Untitled_4.data*20;

%% Encoder processing test

% Pick one dataset (if you already have DataRaw, we can work directly on it)

% Digital channels from DataRaw
zPulse    = double(DataRaw.Digital_channels.Untitled_7.data);  % Z pulse
PulsTrace = double(DataRaw.Digital_channels.Untitled_6.data);  % encoder pulses

% Make sure they are column vectors
zPulse    = zPulse(:);
PulsTrace = PulsTrace(:);

% 1. Z-pulse: find raw TDC markers

% Rising edge = 0 → 1 = TDC marker
Zprev   = [0; zPulse(1:end-1)];                 % dummy 0 at start
isRiseZ = (zPulse == 0) & (Zprev == 1);
TDC_raw = find(isRiseZ);

if numel(TDC_raw) < 2
    error('Not enough Z-pulses to estimate samplesPerCycle.');
end

% Samples per 720° cycle
samplesPerCycle = mean(diff(TDC_raw));

% Define TDC_shift if you haven't already (in crank degrees)
% e.g. TDC_shift = 0;
% TDC_shift should exist before this block
sampleShift = round((TDC_shift / 720) * samplesPerCycle);
TDC_shifted = TDC_raw - sampleShift;

% Remove invalid indices
TDC_shifted(TDC_shifted < 1 | TDC_shifted > numel(zPulse)) = [];


% 2. Encoder pulses → crank angle

% Rising edges of encoder steps
Pprev        = [0; PulsTrace(1:end-1)];
isRiseEncoder = (PulsTrace > 0) & (Pprev == 0);

% Convert pulses to degrees
pulsesPerRev = 1800;                 % (0.2° per pulse)
degPerPulse  = 360 / pulsesPerRev;   % = 0.2
CA_raw       = cumsum(double(isRiseEncoder)) * degPerPulse;


% 3. Smooth crank angle

wf_increment = DataRaw.Analog_channels.Dev0_Ai5.Props.wf_increment;  % [s/sample]

% seconds per encoder pulse = (60/rpm)/pulsesPerRev
% samples per pulse = time_per_pulse / wf_increment
samplesPerPulse = (60 / RPM / pulsesPerRev) / wf_increment;
winLen          = max(1, round(samplesPerPulse * 5));  % ~5-pulse moving average

b = ones(winLen,1) / winLen;
CA_smooth = filtfilt(b, 1, CA_raw);


%4. Align crank angle so TDC = 0, 720, 1440…

nCycles = numel(TDC_shifted);
if nCycles < 1
    error('No valid TDC_shifted indices after bounds check.');
end

offset = mean( CA_smooth(TDC_shifted) - (1:nCycles) * 720 );
CA_aligned = mod(CA_smooth - offset, 720) - 360;


%5. Final TDC from angle wrap
ca       = CA_aligned(:);
ca_next  = [ca(2:end); ca(end)];
TDC_final = find(ca > ca_next);   % local maxima of CA (wrap points)
validIdx = TDC_final(TDC_final >= 1 & TDC_final <= length(SignalPressure));


figure;
subplot(4,1,1)
plot(zPulse); hold on;
plot(TDC_raw,    zPulse(TDC_raw),    'ro','MarkerSize',6);
plot(TDC_shifted,zPulse(TDC_shifted),'gx','MarkerSize',8);
plot(SignalPressure, 'k');
title('Z-pulse: raw TDC (red) and shifted TDC (green)');
xlabel('Samples'); ylabel('Digital Level');
grid on;

subplot(4,1,2)
plot(CA_raw); hold on;
plot(CA_smooth);
legend('Raw CA', 'Smoothed CA');
ylabel('Crank Angle (deg)');
title('Encoder angle reconstruction');
grid on;

subplot(4,1,3)
plot(CA_aligned); hold on;
plot(TDC_final, CA_aligned(TDC_final), 'ro');
ylabel('Aligned CA'); xlabel('Samples');
title('Final TDC positions from CA wrap');
grid on;

subplot(4,1,4)
plot(t, SignalPressure, 'k'); hold on;
plot(t, SignalPressure_Filtered, 'r');
plot(t(validIdx), SignalPressure(validIdx), 'gx', 'MarkerSize', 10, 'LineWidth', 2);
ylabel('Pressure'); xlabel('Time');
title('TDC positions on pressure trace');
grid on;

% validation
% Wrap your DataRaw in a cell so it matches the function’s expected input
handles = {DataRaw};

% Pick an rpm and TDC_shift (deg). For now TDC_shift = 0; later you calibrate it.
[TDCs, CATimeScale, ~, ~, ~] = Find_TDC_with_TimeScale(handles, RPM, TDC_shift);

TDC_idx = TDCs{1};           % indices for this dataset

% Time & pressure
P  = SignalPressure;         % raw or filtered, as you like

figure;
plot(P, 'k'); hold on;
plot(TDC_idx, zPulse(TDC_idx), 'gx', 'MarkerSize', 10, 'LineWidth', 2);
plot(zPulse)
plot(DieselSignal)
plot(MianFuelSignal)
ylabel('Pressure [bar]');
grid on;

%% Plot check
keyboard

dig   = DataRaw.Digital_channels;
fn    = fieldnames(dig);

figure;
plotRow = 1;

for i = 1:numel(fn)
    f = fn{i};

    % Only process fields that actually have a `.data` field
    if isstruct(dig.(f)) && isfield(dig.(f), 'data')
        subplot( numel(fn)-2, 1, plotRow );   % -2 to ignore 'name' and 'Props'
        plot(dig.(f).data);
        
        title(f);
        plotRow = plotRow + 1;
    end
end


%% Testing Function 
function [TDCs, Pressure_all, Injection_all] = Find_TDC(handles,rpm,TDC_shift)
    
    TDCs = cell(length(handles),1);
    Pressure_all = cell(length(handles),1);
    Injection_all = cell(length(handles),1);

    for i = 1:length(handles)
        Zpulse = double(handles{i,1}.Digital_channels.Untitled_7.data);

        % Find the TDC positions
        Zfilt = Zpulse;
        TDC_curr = find(Zfilt==1 & ([1 Zfilt(1:end-1)])==0);
        TDC_diff = mean(diff(TDC_curr));
        TDC_curr = TDC_curr - round((TDC_shift) / 720 * TDC_diff);

        %get rid of the chipped cycles
        if(TDC_curr(1)<1)
            TDC_curr = TDC_curr(2:end);
        end
        if (TDC_curr(end)>length(Zfilt))
            TDC_curr=TDC_curr(1:end-1);
        end
        TDCs{i}=TDC_curr;

        %extract the CA degree from the encoder trace
        PulsTrace = handles{i,1}.Digital_channels.Untitled_6.data;
        CA_Enc_curr{i} = cumsum(double(PulsTrace>0 & ([1 PulsTrace(1:end-1)])==0))*0.2;
        a=1; 
        b=round((60/rpm/1800)/handles{i,1}.Analog_channels.Dev0_Ai5.Props.wf_increment*5);
        b=ones(b,1)/b;
        CA_Enc_curr{i} =filtfilt(b,a,CA_Enc_curr{i});
        CA_Enc_curr{i} = mod((CA_Enc_curr{i}-mean(CA_Enc_curr{i}(TDC_curr) - (1:length(TDC_curr))*720)),720)-360;
        TDCs{i} = find(CA_Enc_curr{i}>[CA_Enc_curr{i}(2:end) CA_Enc_curr{i}(end)]);
    end
end


function [TDCs, CA_Enc, TimeScale, Pressure_all, Injection_all] = Find_TDC_with_TimeScale(handles, rpm, TDC_shift)
% Find_TDC_with_TimeScale - Finds TDC positions and generates CA timescale
%
% Inputs:
%   handles    - Cell array of data structures with digital/analog channels
%   rpm        - Engine RPM
%   TDC_shift  - TDC shift in crank angle degrees (0-720)
%
% Outputs:
%   TDCs          - Cell array of TDC sample indices
%   CA_Enc        - Cell array of crank angle traces (degrees, -360 to 360)
%   TimeScale     - Cell array of time vectors (seconds)
%   Pressure_all  - Cell array for pressure data (placeholder)
%   Injection_all - Cell array for injection data (placeholder)

    % Initialize output cell arrays
    TDCs = cell(length(handles), 1);
    CA_Enc = cell(length(handles), 1);
    TimeScale = cell(length(handles), 1);
    Pressure_all = cell(length(handles), 1);
    Injection_all = cell(length(handles), 1);
    
    for i = 1:length(handles)
        % Extract Z-pulse data (TDC marker signal)
        Zpulse = double(handles{i,1}.Digital_channels.Untitled_7.data);
        
        % Find the TDC positions (falling edge detection)
        Zfilt = Zpulse;
        TDC_curr = find(Zfilt == 0 & ([1 Zfilt(1:end-1)]) == 1);
        
        % Calculate mean cycle length
        TDC_diff = mean(diff(TDC_curr));
        
        % Apply TDC shift correction
        TDC_curr = TDC_curr - round((TDC_shift) / 720 * TDC_diff);
        
        % Remove invalid cycles at boundaries
        if (TDC_curr(1) < 1)
            TDC_curr = TDC_curr(2:end);
        end
        if (TDC_curr(end) > length(Zfilt))
            TDC_curr = TDC_curr(1:end-1);
        end
        
        TDCs{i} = TDC_curr;
        
        % Extract the CA degree from the encoder trace
        PulsTrace = handles{i,1}.Digital_channels.Untitled_6.data;
        
        % Cumulative angle calculation (0.2 deg per pulse)
        CA_Enc_curr = cumsum(double(PulsTrace > 0 & ([1 PulsTrace(1:end-1)]) == 0)) * 0.2;
        
        % Design low-pass filter for CA smoothing
        a = 1;
        % Filter width based on RPM and sampling rate
        dt = handles{i,1}.Analog_channels.Dev0_Ai5.Props.wf_increment;
        filter_width = round((60/rpm/1800) / dt * 5);
        b = ones(filter_width, 1) / filter_width;
        CA_Enc_curr = filtfilt(b, a, CA_Enc_curr);
        
        % Normalize CA to -360 to 360 degree range centered on TDC
        CA_offset = mean(CA_Enc_curr(TDC_curr) - (1:length(TDC_curr)) * 720);
        CA_Enc_curr = mod((CA_Enc_curr - CA_offset), 720) - 360;
        
        % Find updated TDC positions based on CA peaks
        TDCs{i} = find(CA_Enc_curr > [CA_Enc_curr(2:end) CA_Enc_curr(end)]);
        
        % Store CA trace
        CA_Enc{i} = CA_Enc_curr;
        
        % Generate TimeScale
        N_samples = length(CA_Enc_curr);
        TimeScale{i} = (0:N_samples-1)' * dt;
    end
end
