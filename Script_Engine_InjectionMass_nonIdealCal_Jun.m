% H2 Injection Mass Calculation: Non-ideal (SRK) vs Ideal Gas Law

clear; clc; close all

% Add Library_Matlab to path for CLASS_SRK_EOS
addpath('Library_Matlab');

% -----------------------------------------
saveName = 'C:\Users\z5058464\OneDrive - UNSW\CH4DDI\Engine\MassflowRate\Trial';
% -----------------------------------------
% Load data from Excel
fileName = 'C:\Users\jun-y\OneDrive - UNSW\CH4DDI\Engine\MassflowRate\20251215_OldChamber_FlowRate.xlsx'; % Replace with the name of your Excel file
sheetnames = 'H2_350_ST02_NonIdeal_B25_T2';
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
    injection_duration = data.InjectorEnergisingTime_ms_(i); % Injection duration (ms)
    injection_count = data.InjectionQuantity(i);
    
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

% Calculate statistics in single loop
for i = 1:n_durations
    rows = data.InjectorEnergisingTime_ms_ == unique_durations(i);
    
    % Non-ideal statistics
    mean_flow_rate(i) = mean(flowRate_H2(rows));
    std_injected_mass(i) = std(mass_H2_injected(rows));
    mean_injected_mass(i) = mean(mass_H2_injected(rows));
    
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

% Create individual mass table for each injection
injected_mass_table = table(data.InjectorEnergisingTime_ms_, mass_H2_injected, ...
    'VariableNames', {'Injection_Duration_ms', 'Injected_Mass_mg'});

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

%% SRK calculations now handled by CLASS_SRK_EOS
% No inline functions needed - see Library_Matlab/CLASS_SRK_EOS.m

