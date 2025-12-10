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
filePath = 'C:\Users\jun-y\OneDrive\桌面\20250909\20250909\_Set17_1400_30_90\_Set17_1400_30_90_Take00.tdms';
DataRaw  = TDMS_getStruct(filePath);
PressureRaw = double(DataRaw.Analog_channels.Dev0_Ai5.data)*20;
SamplingRate_Local = 200000;
SignalTimeScale = (0:length(PressureRaw)-1)/SamplingRate_Local; 

%% Filter settings
% f_cutoff = 10000;                % 10 kHz cutoff
% fc_norm = f_cutoff / (SamplingRate_Local/2);   % normalised cutoff (0–1)
% 
% [b,a] = butter(2, fc_norm);      % 2nd-order Butterworth filter
% Use PressureRaw as the input signal
%% Build pressure signal & time vector
SignalPressure = PressureRaw;

%% --- Your filter chain ---
SignalPressure_Filtered = bandstop(SignalPressure,[8800 11200],SamplingRate_Local);
SignalPressure_Filtered = bandstop(SignalPressure_Filtered,[36000 42000],SamplingRate_Local);
SignalPressure_Filtered = sgolayfilt(SignalPressure_Filtered, 3, 101);   % order 5, frame 201
SignalPressure_Filtered = bandstop(SignalPressure_Filtered,[3400 3800],SamplingRate_Local);

% If you're inside a class method:
% this.FilteredPressureSignal = SignalPressure_Filtered;

%% Plot time-domain: raw vs filtered
t = SignalTimeScale;

figure;
subplot(2,1,1);
plot(t, SignalPressure, 'k'); hold on;
plot(t, SignalPressure_Filtered, 'r');
xlabel('Time [s]');
ylabel('Pressure [bar]');
title('Pressure Signal: Raw vs Filtered (Full)');
legend('Raw','Filtered');
grid on;

% 👀 Zoom into a region of interest (edit t_min/t_max)
t_min = 0.005;   % change to where your main event is
t_max = 0.015;

idx = t >= t_min & t <= t_max;

subplot(2,1,2);
plot(t(idx), SignalPressure(idx), 'k'); hold on;
plot(t(idx), SignalPressure_Filtered(idx), 'r');
xlabel('Time [s]');
ylabel('Pressure [bar]');
title(sprintf('Zoomed-in Region (%.4f–%.4f s)', t_min, t_max));
legend('Raw','Filtered');
grid on;

%% Frequency-domain: PSD before vs after
[PSD_raw, f_psd]  = pwelch(SignalPressure,[],[],[],SamplingRate_Local);
[PSD_filt, ~]     = pwelch(SignalPressure_Filtered,[],[],[],SamplingRate_Local);

figure;
loglog(f_psd, PSD_raw); hold on;
loglog(f_psd, PSD_filt);
xlabel('Frequency [Hz]');
ylabel('PSD');
title('Power Spectral Density: Raw vs Filtered');
legend('Raw','Filtered');
grid on;
xlim([10 1e5]);   % adjust if needed

noise_est = SignalPressure - SignalPressure_Filtered;

[PSD_noise, f_noise] = pwelch(noise_est,[],[],[],SamplingRate_Local);

figure;
loglog(f_psd, PSD_raw, 'k'); hold on;
loglog(f_noise, PSD_noise, 'r');
legend('Raw','Raw - Filtered (noise estimate)');

sig_raw = SignalPressure - mean(SignalPressure);
sig_flt = SignalPressure_Filtered - mean(SignalPressure_Filtered);

[acor, lag] = xcorr(sig_flt, sig_raw, 'coeff');
[~, I] = max(acor);
lag_samples = lag(I);
lag_time    = lag_samples / SamplingRate_Local;   % seconds

disp(lag_time);
dPdt = gradient(SignalPressure_Filtered, SignalTimeScale);   % same size as PressureRaw
figure
plot(dPdt)
grid on
