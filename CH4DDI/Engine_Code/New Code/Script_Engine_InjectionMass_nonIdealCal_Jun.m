% H2 Injection Mass Calculation: Non-ideal (SRK) vs Ideal Gas Law

clear; clc; close all

% Add Library_Matlab to path for CLASS_SRK_EOS
addpath('Library_Matlab');

% -----------------------------------------
%saveName = 'C:\Users\z5058464\OneDrive - UNSW\CH4DDI\Engine\MassflowRate\Trial';
% -----------------------------------------
% Load data from Excel
fileName = 'C:\Users\jun-y\OneDrive - UNSW\CH4DDI\Engine\MassflowRate\20251215_OldChamber_FlowRate.xlsx'; % Replace with the name of your Excel file
sheetnames = 'H2_350_ST02_NonIdeal_B0_T2';
data = readtable(fileName,'Sheet',sheetnames);
% -------------------------------------------------------------------------

% Given constants
% === Gas Properties ===
R = 8.314; % Universal gas constant, J/(mol·K)
R_bar = 0.08314; % Universal gas constant, L·bar/(mol·K)
M_H2 = 2.016; % Molar mass of hydrogen, g/mol
M_N2 = 28.014; % Molar mass of nitrogen, g/mol

% === Critical Properties ===
% Hydrogen
Pc_H2 = 1.29e6; % Critical pressure, Pa (1.29 MPa)
Tc_H2 = 33.2; % Critical temperature, K
omega_H2 = -0.22; % Acentric factor

% Nitrogen
Pc_N2 = 3.39e6; % Critical pressure, Pa (3.39 MPa)
Tc_N2 = 126.2; % Critical temperature, K
omega_N2 = 0.04; % Acentric factor

% === SRK EOS Coefficients ===
SRK_a_coeff = 0.42748; % SRK attraction parameter coefficient
SRK_b_coeff = 0.08664; % SRK repulsion parameter coefficient
SRK_m_c0 = 0.48; % SRK m parameter constant
SRK_m_c1 = 1.574; % SRK m parameter linear coefficient
SRK_m_c2 = 0.176; % SRK m parameter quadratic coefficient

% === Chamber and Pressure ===
V = 0.0014; % Chamber volume, m^3
P_atm_bar = 1.01325; % Atmospheric pressure to convert gauge -> absolute, bar

% === Conversion Factors ===
bar_to_Pa = 1e5; % Conversion: bar to Pa
bar_to_MPa = 0.1; % Conversion: bar to MPa
MPa_to_Pa = 1e6; % Conversion: MPa to Pa
m3_to_L = 1e3; % Conversion: m^3 to L
mol_to_mg = 1000; % Conversion: mol to mg (when multiplied by molar mass)

% === Numerical Tolerances ===
iter_tol = 1e-5; % Tolerance for iterative convergence
iter_max = 100; % Maximum iterations for convergence
imag_tol = 1e-6; % Tolerance for imaginary part in root finding

% -------------------------------------------------------------------------

% Initialize results
mass_H2_injected = zeros(height(data), 1); % Mass of hydrogen injected (mg)
flowRate_H2 = zeros(height(data), 1); % Flow rate (mg/ms)
% Ideal gas counterparts
mass_H2_injected_ideal = zeros(height(data), 1); % Ideal-gas injected mass (mg)
flowRate_H2_ideal = zeros(height(data), 1); % Ideal-gas flow rate (mg/ms)
% (diagnostics removed)

for i = 1:height(data)
    % Extract data for the current row
    T = data.InitialTemperature_K_T_1(i); % Temperature (K)
    P_diff_bar = data.PressureRiseDueToH2Injection_bar_P_diff(i); % Pressure rise (gauge, bar)
    Pb_gauge_bar = data.BackPressure_bar_P_b(i); % back pressure (gauge, bar)
    P_total_gauge_bar = P_diff_bar + Pb_gauge_bar; % gauge pressure after injection (bar)

    % Convert to absolute pressures for all thermodynamic calculations
    Pb_abs_bar = Pb_gauge_bar + P_atm_bar;
    P_total_abs_bar = P_total_gauge_bar + P_atm_bar;
    % store for cross-check later (init vectors on first iteration)
    if i == 1
        P_total_abs_bar_vec = zeros(height(data),1);
        Pb_abs_bar_vec = zeros(height(data),1);
        T_vec = zeros(height(data),1);
        inj_qty_vec = zeros(height(data),1);
        n_total_nonideal_vec = zeros(height(data),1);
        n_N2_nonideal_vec = zeros(height(data),1);
    end
    P_total_abs_bar_vec(i) = P_total_abs_bar;
    Pb_abs_bar_vec(i) = Pb_abs_bar;
    T_vec(i) = T;
    injection_duration = data.InjectorEnergisingTime_ms_(i); % Injection duration (ms)
    injection_count = data.InjectionQuantity(i);
    inj_qty_vec(i) = injection_count;
    
    % === NON-IDEAL: Use SRK mixture Z calculation ===
    Pc_MPa_vec = [Pc_H2, Pc_N2] / 1e6; % convert Pa to MPa for SRK
    Tc_vec = [Tc_H2, Tc_N2];
    omega_vec = [omega_H2, omega_N2];
    [Z_mix, x_H2] = CLASS_SRK_EOS.CalculateBinaryMixtureZ(T, P_total_abs_bar, V, Pb_abs_bar, ...
        Pc_MPa_vec, Tc_vec, omega_vec);
    
    % Total moles after injection
    n_total = (P_total_abs_bar * bar_to_Pa * V) / (Z_mix * R * T);
    
    % Initial N2 moles (pure N2 at back pressure) using SRK
    Z_N2 = CLASS_SRK_EOS.CalculateSingleComponentZ(Pb_abs_bar * bar_to_Pa, T, Pc_N2, Tc_N2, omega_N2);
    n_N2 = (Pb_abs_bar * bar_to_Pa * V) / (Z_N2 * R * T);
    n_total_nonideal_vec(i) = n_total;
    n_N2_nonideal_vec(i) = n_N2;
    
    % H2 moles injected (guard against numerical negatives)
    n_H2 = max(n_total - n_N2, 0);
    
    % Calculate injected hydrogen mass
    mass_H2_injected(i) = (n_H2 * M_H2 * mol_to_mg) / injection_count;
    
    % Calculate flow rate
    flowRate_H2(i) = mass_H2_injected(i) / injection_duration;

    % === IDEAL: Use ideal gas law (absolute pressures) ===
    % Initial N2 moles (pure N2 at back pressure)
    n_N2_ideal = (Pb_abs_bar * bar_to_Pa * V) / (R * T);
    
    % Total moles after injection (mixture ideal)
    n_total_ideal = (P_total_abs_bar * bar_to_Pa * V) / (R * T);
    
    % Mole difference gives injected H2 directly under ideal assumption
    n_H2_ideal = max(n_total_ideal - n_N2_ideal, 0);
    % store for cross-check
    if i == 1
        n_H2_ideal_vec = zeros(height(data),1);
    end
    n_H2_ideal_vec(i) = n_H2_ideal;
    
    % Molar-difference path (source of truth)
    mass_H2_injected_ideal(i) = (n_H2_ideal * M_H2 * mol_to_mg) / injection_count;
    flowRate_H2_ideal(i) = mass_H2_injected_ideal(i) / injection_duration;
end

% Append results to the table and save
data.Mass_H2_Injected_mg = mass_H2_injected;
data.FlowRate_H2_mg_per_ms = flowRate_H2;
% Append ideal-gas results for comparison
data.Mass_H2_Injected_mg_Ideal = mass_H2_injected_ideal;
data.FlowRate_H2_Ideal_mg_per_ms = flowRate_H2_ideal;
% (diagnostics removed)

% Compute cumulative injection counts per sequence (e.g., raw 20,20,20 -> 20,40,60)
durations = data.InjectorEnergisingTime_ms_;
uDur = unique(durations);
data.InjectionQuantity_Cum = NaN(height(data),1);
for ui = 1:numel(uDur)
    rows = find(durations == uDur(ui));
    qty_seq = inj_qty_vec(rows);
    cum_seq = cumsum(qty_seq);
    data.InjectionQuantity_Cum(rows) = cum_seq;
    % Validate cumulative is strictly increasing
    if any(diff(cum_seq) <= 0)
        warning('Sequence %g: Cumulative injection count not strictly increasing. Raw: %s', uDur(ui), mat2str(qty_seq(:)'));
    end
end

% Compute per-injection mass using from-initial cumulative (Method C) for non-ideal (SRK)
data.mH2_NonIdeal_PerInj_Blockwise = mass_H2_injected;
data.mH2_NonIdeal_PerInj_FromInitial = NaN(height(data),1);
for ui = 1:numel(uDur)
    rows = find(durations == uDur(ui));
    % Order by cumulative injection count (handles repeated raw entries)
    [~, ord] = sort(data.InjectionQuantity_Cum(rows));
    rows = rows(ord);
    cum_seq = data.InjectionQuantity_Cum(rows);
    r0 = rows(1);
    nN2_base = n_N2_nonideal_vec(r0);
    N_cum_nonideal = n_total_nonideal_vec(rows) - nN2_base;
    m_perinj_fromInit = NaN(size(rows));
    idx_pos = cum_seq > 0;
    m_perinj_fromInit(idx_pos) = (N_cum_nonideal(idx_pos) * M_H2 * mol_to_mg) ./ cum_seq(idx_pos);
    data.mH2_NonIdeal_PerInj_FromInitial(rows) = m_perinj_fromInit;
end

% Compute per-injection mass using from-initial cumulative (Method C) for ideal
data.mH2_Ideal_PerInj_Blockwise = mass_H2_injected_ideal;
data.mH2_Ideal_PerInj_FromInitial = NaN(height(data),1);
for ui = 1:numel(uDur)
    rows = find(durations == uDur(ui));
    % Order by cumulative injection count (handles repeated raw entries)
    [~, ord] = sort(data.InjectionQuantity_Cum(rows));
    rows = rows(ord);
    cum_seq = data.InjectionQuantity_Cum(rows);
    r0 = rows(1);
    P0_abs = Pb_abs_bar_vec(r0);
    T0 = T_vec(r0);
    N_cum = (P_total_abs_bar_vec(rows) * bar_to_Pa * V) ./ (R * T_vec(rows)) - (P0_abs * bar_to_Pa * V) / (R * T0);
    m_perinj_fromInit = NaN(size(rows));
    idx_pos = cum_seq > 0;
    m_perinj_fromInit(idx_pos) = (N_cum(idx_pos) * M_H2 * mol_to_mg) ./ cum_seq(idx_pos);
    data.mH2_Ideal_PerInj_FromInitial(rows) = m_perinj_fromInit;
end

% Calculate statistics for both non-ideal and ideal results
unique_durations = unique(data.InjectorEnergisingTime_ms_);
n_durations = length(unique_durations);

% Pre-allocate arrays
mean_flow_rate = zeros(n_durations, 1);
std_injected_mass = zeros(n_durations, 1);
mean_injected_mass = zeros(n_durations, 1);
mean_flow_rate_ideal = zeros(n_durations, 1);
std_injected_mass_ideal = zeros(n_durations, 1);
mean_injected_mass_ideal = zeros(n_durations, 1);
mean_injected_mass_nonideal_fromInit = zeros(n_durations, 1);
std_injected_mass_nonideal_fromInit = zeros(n_durations, 1);

% Calculate statistics in single loop
for i = 1:n_durations
    rows = data.InjectorEnergisingTime_ms_ == unique_durations(i);
    
    % Non-ideal statistics
    mean_flow_rate(i) = mean(flowRate_H2(rows));
    std_injected_mass(i) = std(mass_H2_injected(rows));
    mean_injected_mass(i) = mean(mass_H2_injected(rows));
    valid_nonideal_fromInit = rows & (data.InjectionQuantity_Cum > 0) & ~isnan(data.mH2_NonIdeal_PerInj_FromInitial);
    mean_injected_mass_nonideal_fromInit(i) = mean(data.mH2_NonIdeal_PerInj_FromInitial(valid_nonideal_fromInit));
    std_injected_mass_nonideal_fromInit(i) = std(data.mH2_NonIdeal_PerInj_FromInitial(valid_nonideal_fromInit));
    
    % Ideal statistics
    mean_flow_rate_ideal(i) = mean(flowRate_H2_ideal(rows));
    std_injected_mass_ideal(i) = std(mass_H2_injected_ideal(rows));
    mean_injected_mass_ideal(i) = mean(mass_H2_injected_ideal(rows));
end

% Create summary table for mean and std results
summary_table = table(unique_durations, mean_flow_rate, std_injected_mass, mean_injected_mass, ...
    'VariableNames', {'Injection_Duration_ms', 'Mean_FlowRate_mg_per_ms', 'Std_Injected_Mass_mg', 'Mean_Injected_Mass_mg'});

% Ideal-gas summary table
summary_table_ideal = table(unique_durations, mean_flow_rate_ideal, std_injected_mass_ideal, mean_injected_mass_ideal, ...
    'VariableNames', {'Injection_Duration_ms', 'Mean_FlowRate_Ideal_mg_per_ms', 'Std_Injected_Mass_Ideal_mg', 'Mean_Injected_Mass_Ideal_mg'});

% Non-ideal from-initial (Method C) summary table
summary_table_nonideal_fromInit = table(unique_durations, mean_injected_mass_nonideal_fromInit, std_injected_mass_nonideal_fromInit, ...
    'VariableNames', {'Injection_Duration_ms', 'Mean_Injected_Mass_NonIdeal_FromInit_mg', 'Std_Injected_Mass_NonIdeal_FromInit_mg'});

% Create individual mass table for each injection
injected_mass_table = table(data.InjectorEnergisingTime_ms_, mass_H2_injected, ...
    'VariableNames', {'Injection_Duration_ms', 'Injected_Mass_mg'});

% Individual per-injection from-initial tables
injected_mass_table_nonideal_fromInit = table(data.InjectorEnergisingTime_ms_, data.mH2_NonIdeal_PerInj_FromInitial, ...
    'VariableNames', {'Injection_Duration_ms', 'Injected_Mass_FromInit_mg'});
injected_mass_table_ideal_fromInit = table(data.InjectorEnergisingTime_ms_, data.mH2_Ideal_PerInj_FromInitial, ...
    'VariableNames', {'Injection_Duration_ms', 'Injected_Mass_FromInit_mg'});

% Ideal-gas individual mass table
injected_mass_table_ideal = table(data.InjectorEnergisingTime_ms_, mass_H2_injected_ideal, ...
    'VariableNames', {'Injection_Duration_ms', 'Injected_Mass_mg'});

% Linear regression (filter duration >= 0.6 ms)
filtered_idx = summary_table.Injection_Duration_ms >= 0.6;

% Non-ideal fit
x_filtered1 = summary_table.Injection_Duration_ms(filtered_idx);
y_filtered1 = summary_table.Mean_Injected_Mass_mg(filtered_idx);
std_filtered1 = summary_table.Std_Injected_Mass_mg(filtered_idx);
coefficients1 = polyfit(x_filtered1, y_filtered1, 1);
y_fit1 = polyval(coefficients1, x_filtered1);
SS_res1 = sum((y_filtered1 - y_fit1).^2);
SS_tot1 = sum((y_filtered1 - mean(y_filtered1)).^2);
R_squared1 = 1 - (SS_res1 / SS_tot1);
y_upper1 = y_fit1 + std_filtered1;
y_lower1 = y_fit1 - std_filtered1;
fitEqnNonIdeal = sprintf('Non-ideal: y = %.3fx + %.3f\nR^2 = %.4f', coefficients1(1), coefficients1(2), R_squared1)

% Fit for Non-ideal (from-initial average per injection, Method C)
x_filtered1c = summary_table_nonideal_fromInit.Injection_Duration_ms(filtered_idx);
y_filtered1c = summary_table_nonideal_fromInit.Mean_Injected_Mass_NonIdeal_FromInit_mg(filtered_idx);
std_filtered1c = summary_table_nonideal_fromInit.Std_Injected_Mass_NonIdeal_FromInit_mg(filtered_idx);
coefficients1c = polyfit(x_filtered1c, y_filtered1c, 1);
y_fit1c = polyval(coefficients1c, x_filtered1c);
SS_res1c = sum((y_filtered1c - y_fit1c).^2);
SS_tot1c = sum((y_filtered1c - mean(y_filtered1c)).^2);
R_squared1c = 1 - (SS_res1c / SS_tot1c);
y_upper1c = y_fit1c + std_filtered1c;
y_lower1c = y_fit1c - std_filtered1c;
fitEqnNonIdeal_FromInitAvg = sprintf('Non-ideal (from-initial avg): y = %.3fx + %.3f\nR^2 = %.4f', coefficients1c(1), coefficients1c(2), R_squared1c)

% Ideal fit (same filter)
x_filtered2 = summary_table_ideal.Injection_Duration_ms(filtered_idx);
y_filtered2 = summary_table_ideal.Mean_Injected_Mass_Ideal_mg(filtered_idx);
std_filtered2 = summary_table_ideal.Std_Injected_Mass_Ideal_mg(filtered_idx);
coefficients2 = polyfit(x_filtered2, y_filtered2, 1);
y_fit2 = polyval(coefficients2, x_filtered2);
SS_res2 = sum((y_filtered2 - y_fit2).^2);
SS_tot2 = sum((y_filtered2 - mean(y_filtered2)).^2);
R_squared2 = 1 - (SS_res2 / SS_tot2);
y_upper2 = y_fit2 + std_filtered2;
y_lower2 = y_fit2 - std_filtered2;
fitEqnIdeal = sprintf('Ideal: y = %.3fx + %.3f\nR^2 = %.4f', coefficients2(1), coefficients2(2), R_squared2)

% Fit for Ideal (from-initial average per injection, Method C)
mean_injected_mass_ideal_fromInit = zeros(n_durations, 1);
std_injected_mass_ideal_fromInit = zeros(n_durations, 1);
for i = 1:n_durations
    rows = data.InjectorEnergisingTime_ms_ == unique_durations(i);
    valid = rows & (data.InjectionQuantity_Cum > 0) & ~isnan(data.mH2_Ideal_PerInj_FromInitial);
    mean_injected_mass_ideal_fromInit(i) = mean(data.mH2_Ideal_PerInj_FromInitial(valid));
    std_injected_mass_ideal_fromInit(i) = std(data.mH2_Ideal_PerInj_FromInitial(valid));
end
summary_table_ideal_fromInit = table(unique_durations, mean_injected_mass_ideal_fromInit, std_injected_mass_ideal_fromInit, ...
    'VariableNames', {'Injection_Duration_ms', 'Mean_Injected_Mass_Ideal_FromInit_mg', 'Std_Injected_Mass_Ideal_FromInit_mg'});

% Linear regression (same duration filter)
x_filtered3 = summary_table_ideal_fromInit.Injection_Duration_ms(filtered_idx);
y_filtered3 = summary_table_ideal_fromInit.Mean_Injected_Mass_Ideal_FromInit_mg(filtered_idx);
std_filtered3 = summary_table_ideal_fromInit.Std_Injected_Mass_Ideal_FromInit_mg(filtered_idx);
coefficients3 = polyfit(x_filtered3, y_filtered3, 1);
y_fit3 = polyval(coefficients3, x_filtered3);
SS_res3 = sum((y_filtered3 - y_fit3).^2);
SS_tot3 = sum((y_filtered3 - mean(y_filtered3)).^2);
R_squared3 = 1 - (SS_res3 / SS_tot3);
y_upper3 = y_fit3 + std_filtered3;
y_lower3 = y_fit3 - std_filtered3;
fitEqnIdeal_FromInitAvg = sprintf('Ideal (from-initial avg): y = %.3fx + %.3f\nR^2 = %.4f', coefficients3(1), coefficients3(2), R_squared3)

%% Cross-check: block-wise vs cumulative-from-initial (ideal gas)
% We compare per-block moles computed two ways:
%   Method A (block-wise): n_H2_ideal_vec (current) = (P_after - P_before) V / (R T)
%   Method B (from-initial): N_k = (P_k/T_k - P_0/T_0) V/R; per-block = N_k - N_{k-1}
% Groups are split by energizing time to avoid mixing separate sequences.

% Pre-allocate outputs
data.nH2_Ideal_Blockwise = n_H2_ideal_vec;
data.nH2_Ideal_FromInitial = NaN(height(data),1);
data.nH2_Ideal_Diff = NaN(height(data),1);
data.nH2_Ideal_DiffRel = NaN(height(data),1);

durations = data.InjectorEnergisingTime_ms_;
uDur = unique(durations);
for ui = 1:numel(uDur)
    rows = find(durations == uDur(ui));
    % sort rows by cumulative InjectionQuantity (e.g., 0,20,40,...)
    [~, ord] = sort(data.InjectionQuantity_Cum(rows));
    rows = rows(ord);
    % baseline (assumed first in sequence)
    r0 = rows(1);
    P0_abs = Pb_abs_bar_vec(r0); % initial absolute back pressure
    T0 = T_vec(r0);
    % cumulative from initial for each row
    S = (P_total_abs_bar_vec(rows) * bar_to_Pa * V) ./ (R * T_vec(rows));
    S0 = (P0_abs * bar_to_Pa * V) / (R * T0);
    N_cum = S - S0;
    % per-block from cumulative: diff with previous row in group
    n_from_initial = NaN(size(rows));
    for k = 2:numel(rows)
        n_from_initial(k) = N_cum(k) - N_cum(k-1);
    end
    % write back
    data.nH2_Ideal_FromInitial(rows) = n_from_initial;
    % differences where both are defined
    idx_valid = rows(2:end);
    diff_loc = n_from_initial(2:end) - n_H2_ideal_vec(idx_valid);
    data.nH2_Ideal_Diff(idx_valid) = diff_loc;
    data.nH2_Ideal_DiffRel(idx_valid) = diff_loc ./ max(1e-12, n_H2_ideal_vec(idx_valid));
end

% Report cross-check statistics
valid = ~isnan(data.nH2_Ideal_Diff);
if any(valid)
    absDiff = abs(data.nH2_Ideal_Diff(valid));
    maxAbs = max(absDiff);
    meanAbs = mean(absDiff);
    absDiffSorted = sort(absDiff);
    p95Idx = max(1, round(0.95 * numel(absDiffSorted)));
    p95Abs = absDiffSorted(p95Idx);
    maxRel = max(abs(data.nH2_Ideal_DiffRel(valid)));
    fprintf('Cross-check (ideal): maxAbs=%.3e mol, meanAbs=%.3e mol, p95Abs=%.3e mol, maxRel=%.3e (N=%d)\n',...
        maxAbs, meanAbs, p95Abs, maxRel, numel(absDiff));
end

%% Cross-check 2a: per-injection from-initial avg vs block-wise (non-ideal SRK)
% Method C (from-initial avg): m_per_inj(k) = [ n_total(k) - n_N2(base) ] * M_H2 / InjQty_cum_k
% Method A (block-wise): existing mass_H2_injected (per-injection within block)

valid2_non = ~isnan(data.mH2_NonIdeal_PerInj_FromInitial) & ~isnan(data.mH2_NonIdeal_PerInj_Blockwise);
if any(valid2_non)
    d2n = data.mH2_NonIdeal_PerInj_FromInitial(valid2_non) - data.mH2_NonIdeal_PerInj_Blockwise(valid2_non);
    maxAbs2n = max(abs(d2n));
    meanAbs2n = mean(abs(d2n));
    d2ns = sort(abs(d2n));
    p95Abs2n = d2ns(max(1, round(0.95*numel(d2ns))));
    maxRel2n = max(abs(d2n) ./ max(1e-12, data.mH2_NonIdeal_PerInj_Blockwise(valid2_non)));
    fprintf('Cross-check (non-ideal per-inj avg vs block): maxAbs=%.3e mg, meanAbs=%.3e mg, p95Abs=%.3e mg, maxRel=%.3e (N=%d)\n',...
        maxAbs2n, meanAbs2n, p95Abs2n, maxRel2n, numel(d2n));
end

%% Cross-check 2: per-injection from-initial avg vs block-wise (ideal gas)
% Method C (from-initial avg): m_per_inj(k) = [ (P_k/T_k - P_0/T_0) V/R ] * M_H2 / InjQty_k
% Method A (block-wise): existing mass_H2_injected_ideal (per-injection within block)

% Differences vs block-wise per-injection
valid2 = ~isnan(data.mH2_Ideal_PerInj_FromInitial) & ~isnan(data.mH2_Ideal_PerInj_Blockwise);
if any(valid2)
    d2 = data.mH2_Ideal_PerInj_FromInitial(valid2) - data.mH2_Ideal_PerInj_Blockwise(valid2);
    maxAbs2 = max(abs(d2));
    meanAbs2 = mean(abs(d2));
    d2s = sort(abs(d2));
    p95Abs2 = d2s(max(1, round(0.95*numel(d2s))));
    maxRel2 = max(abs(d2) ./ max(1e-12, data.mH2_Ideal_PerInj_Blockwise(valid2)));
    fprintf('Cross-check (per-inj avg vs block): maxAbs=%.3e mg, meanAbs=%.3e mg, p95Abs=%.3e mg, maxRel=%.3e (N=%d)\n',...
        maxAbs2, meanAbs2, p95Abs2, maxRel2, numel(d2));
end
%% Plot
% Create a new figure using custom axes class
obj_ax = CLASS_AxesHandleStore;
obj_ax.MarginLeft = 45;
obj_ax.MarginRight = 30;
obj_ax.MarginBottom = 50;
obj_ax.MarginTop = 30;
obj_ax.AxesWidth = 300;
obj_ax.AxesHeight = 300;
obj_ax.RowNumber = 2;
obj_ax.ColumnNumber = 1;
GapRow = ones(1, obj_ax.RowNumber) * 15;
GapCol = ones(1, obj_ax.ColumnNumber) * 5;
obj_ax.GapRow = GapRow;
obj_ax.GapColumn = GapCol;
obj_ax = obj_ax.ConstuctAxesUnevenGap();
ax_all = obj_ax.GetAxesHandleMatrix();
ax_non = ax_all(1); % non-ideal plot
ax_ideal = ax_all(2); % ideal plot

% Non-ideal plot (top)
hold(ax_non, "on");
fill(ax_non, [x_filtered1; flipud(x_filtered1)], [y_upper1; flipud(y_lower1)], [0.82 0.93 0.98], 'EdgeColor', 'none');
plot(ax_non, injected_mass_table.Injection_Duration_ms, injected_mass_table.Injected_Mass_mg, 'o', 'MarkerFaceColor', [0.53 0.74 0.87], 'MarkerEdgeColor', 'none');
plot(ax_non, summary_table.Injection_Duration_ms, summary_table.Mean_Injected_Mass_mg, 'o', 'MarkerFaceColor', [0.05 0.33 0.65], 'MarkerEdgeColor', 'none');
plot(ax_non, x_filtered1, y_fit1, 'Color', [0.05 0.33 0.65], 'LineWidth', 1.5);
xticks(ax_non, [0.6 0.8 1 2 3 4 5 6 8 10]);
xtickangle(ax_non, 45);
xticklabels(ax_non, {});
yticks(ax_non, [0 2 4 6 8 10 12 14]);
xlim(ax_non, [0.45, 6.25]);
ylim(ax_non, [0, 12]);
CLASS_Utilis.InsertFigureText(ax_non, 0.3, 0.78, fitEqnNonIdeal);
title(ax_non, sheetnames, 'Interpreter', 'none');

% Ideal plot (bottom)
hold(ax_ideal, "on");
fill(ax_ideal, [x_filtered2; flipud(x_filtered2)], [y_upper2; flipud(y_lower2)], [0.97 0.89 0.82], 'EdgeColor', 'none');
plot(ax_ideal, injected_mass_table_ideal.Injection_Duration_ms, injected_mass_table_ideal.Injected_Mass_mg, 'o', 'MarkerFaceColor', [0.96 0.69 0.25], 'MarkerEdgeColor', 'none');
plot(ax_ideal, summary_table_ideal.Injection_Duration_ms, summary_table_ideal.Mean_Injected_Mass_Ideal_mg, 'o', 'MarkerFaceColor', [0.83 0.33 0], 'MarkerEdgeColor', 'none');
plot(ax_ideal, x_filtered2, y_fit2, 'Color', [0.83 0.33 0], 'LineWidth', 1.5);
xticks(ax_ideal, [0.6 0.8 1 2 3 4 5 6 8 10]);
xtickangle(ax_ideal, 45);
yticks(ax_ideal, [0 2 4 6 8 10 12 14]);
xlim(ax_ideal, [0.45, 6.25]);
ylim(ax_ideal, [0, 12]);
CLASS_Utilis.InsertFigureText(ax_ideal, 0.3, 0.78, fitEqnIdeal);
xlabel(ax_ideal, 'Energizing time [ms]');
ylabel(ax_ideal, 'Injected mass [mg]');

% From-initial average plot (duplicate layout)
obj_ax2 = CLASS_AxesHandleStore;
obj_ax2.MarginLeft = 45;
obj_ax2.MarginRight = 30;
obj_ax2.MarginBottom = 50;
obj_ax2.MarginTop = 30;
obj_ax2.AxesWidth = 300;
obj_ax2.AxesHeight = 300;
obj_ax2.RowNumber = 2;
obj_ax2.ColumnNumber = 1;
GapRow2 = ones(1, obj_ax2.RowNumber) * 15;
GapCol2 = ones(1, obj_ax2.ColumnNumber) * 5;
obj_ax2.GapRow = GapRow2;
obj_ax2.GapColumn = GapCol2;
obj_ax2 = obj_ax2.ConstuctAxesUnevenGap();
ax_all2 = obj_ax2.GetAxesHandleMatrix();
ax_non_c = ax_all2(1); % non-ideal from-initial plot
ax_ideal_c = ax_all2(2); % ideal from-initial plot

% Non-ideal from-initial plot (top)
hold(ax_non_c, "on");
fill(ax_non_c, [x_filtered1c; flipud(x_filtered1c)], [y_upper1c; flipud(y_lower1c)], [0.82 0.93 0.98], 'EdgeColor', 'none');
plot(ax_non_c, injected_mass_table_nonideal_fromInit.Injection_Duration_ms, injected_mass_table_nonideal_fromInit.Injected_Mass_FromInit_mg, 'o', 'MarkerFaceColor', [0.53 0.74 0.87], 'MarkerEdgeColor', 'none');
plot(ax_non_c, summary_table_nonideal_fromInit.Injection_Duration_ms, summary_table_nonideal_fromInit.Mean_Injected_Mass_NonIdeal_FromInit_mg, 'o', 'MarkerFaceColor', [0.05 0.33 0.65], 'MarkerEdgeColor', 'none');
plot(ax_non_c, x_filtered1c, y_fit1c, 'Color', [0.05 0.33 0.65], 'LineWidth', 1.5);
xticks(ax_non_c, [0.6 0.8 1 2 3 4 5 6 8 10]);
xtickangle(ax_non_c, 45);
xticklabels(ax_non_c, {});
yticks(ax_non_c, [0 2 4 6 8 10 12 14]);
xlim(ax_non_c, [0.45, 6.25]);
ylim(ax_non_c, [0, 12]);
CLASS_Utilis.InsertFigureText(ax_non_c, 0.05, 0.85, fitEqnNonIdeal_FromInitAvg);
title(ax_non_c, [sheetnames ' (From-initial avg)'], 'Interpreter', 'none');

% Ideal from-initial plot (bottom)
hold(ax_ideal_c, "on");
fill(ax_ideal_c, [x_filtered3; flipud(x_filtered3)], [y_upper3; flipud(y_lower3)], [0.97 0.89 0.82], 'EdgeColor', 'none');
plot(ax_ideal_c, injected_mass_table_ideal_fromInit.Injection_Duration_ms, injected_mass_table_ideal_fromInit.Injected_Mass_FromInit_mg, 'o', 'MarkerFaceColor', [0.96 0.69 0.25], 'MarkerEdgeColor', 'none');
plot(ax_ideal_c, summary_table_ideal_fromInit.Injection_Duration_ms, summary_table_ideal_fromInit.Mean_Injected_Mass_Ideal_FromInit_mg, 'o', 'MarkerFaceColor', [0.83 0.33 0], 'MarkerEdgeColor', 'none');
plot(ax_ideal_c, x_filtered3, y_fit3, 'Color', [0.83 0.33 0], 'LineWidth', 1.5);
xticks(ax_ideal_c, [0.6 0.8 1 2 3 4 5 6 8 10]);
xtickangle(ax_ideal_c, 45);
yticks(ax_ideal_c, [0 2 4 6 8 10 12 14]);
xlim(ax_ideal_c, [0.45, 6.25]);
ylim(ax_ideal_c, [0, 12]);
CLASS_Utilis.InsertFigureText(ax_ideal_c, 0.05, 0.85, fitEqnIdeal_FromInitAvg);
xlabel(ax_ideal_c, 'Energizing time [ms]');
ylabel(ax_ideal_c, 'Injected mass (from-initial avg) [mg]');

%% SRK calculations now handled by CLASS_SRK_EOS
% No inline functions needed - see Library_Matlab/CLASS_SRK_EOS.m

