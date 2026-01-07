% CH4 Injection Mass Calculation: Non-ideal (PR) vs Ideal Gas Law
clear; clc; close all;

% Add Library_Matlab to path for CLASS_SRK_EOS
addpath('Library_Matlab');

% -----------------------------------------
saveName = 'C:\Users\z5058464\OneDrive - UNSW\CH4DDI\Engine\MassflowRate\Trial'; %#ok<NASGU>
% -----------------------------------------

% Load data (CH4 sheet)
fileName   = 'C:\Users\jun-y\OneDrive - UNSW\CH4DDI\Engine\MassflowRate\20251215_OldChamber_FlowRate.xlsx';
sheetnames = 'CH4_350_ST02_NonIdeal_B50_T2';
data = readtable(fileName, 'Sheet', sheetnames);

% Given constants (calculation only)
R = 8.314;           % J/(mol*K)
R_bar = 0.08314;     % L·bar/(mol*K)
V = 0.0014;          % m^3
M_CH4 = 16.04;       % g/mol
M_N2  = 28.014;      % g/mol
P_atm_bar = 1.01325; % bar

% Critical properties (CH4, N2)
Pc_CH4 = 4.599e6; Tc_CH4 = 190.56; omega_CH4 = 0.01142; % Pa, K, -
Pc_N2  = 3.39e6;  Tc_N2  = 126.2;  omega_N2  = 0.04;   % Pa, K, -

% Conversions
bar_to_Pa = 1e5; mol_to_mg = 1000; bar_to_MPa = 0.1; MPa_to_Pa = 1e6; m3_to_L = 1e3;

% Numerical tolerances
iter_tol = 1e-5; iter_max = 100; imag_tol = 1e-6;

% Debug flag
DEBUG_SRK = false;

% Initialize results
mass_CH4_injected = zeros(height(data), 1);
flowRate_CH4 = zeros(height(data), 1);
mass_CH4_injected_ideal = zeros(height(data), 1);
flowRate_CH4_ideal = zeros(height(data), 1);

for i = 1:height(data)
    T = data.InitialTemperature_K_T_1(i);               % K
    % Allow either CH4 or legacy H2 column name
    if ismember('PressureRiseDueToCH4Injection_bar_P_diff', data.Properties.VariableNames)
        P_diff_bar = data.PressureRiseDueToCH4Injection_bar_P_diff(i);
    else
        P_diff_bar = data.PressureRiseDueToH2Injection_bar_P_diff(i);
    end
    Pb_gauge_bar = data.BackPressure_bar_P_b(i);        % bar (gauge)
    P_total_gauge_bar = P_diff_bar + Pb_gauge_bar;      % bar (gauge after injection)
    injection_duration = data.InjectorEnergisingTime_ms_(i); % ms
    injection_count = data.InjectionQuantity(i);

    Pb_abs_bar = Pb_gauge_bar + P_atm_bar;
    P_total_abs_bar = P_total_gauge_bar + P_atm_bar;

    % Non-ideal (SRK mixture, CH4-N2) using CLASS_SRK_EOS
    Pc_MPa_vec = [Pc_CH4, Pc_N2] / 1e6; % convert Pa to MPa
    Tc_vec = [Tc_CH4, Tc_N2];
    omega_vec = [omega_CH4, omega_N2];
    [Z_mix, x_CH4] = CLASS_SRK_EOS.CalculateBinaryMixtureZ(T, P_total_abs_bar, V, Pb_abs_bar, ...
        Pc_MPa_vec, Tc_vec, omega_vec);
    n_total = (P_total_abs_bar * bar_to_Pa * V) / (Z_mix * R * T);

    Z_N2 = CLASS_SRK_EOS.CalculateSingleComponentZ(Pb_abs_bar * bar_to_Pa, T, Pc_N2, Tc_N2, omega_N2);
    n_N2 = (Pb_abs_bar * bar_to_Pa * V) / (Z_N2 * R * T);

    n_CH4 = max(n_total - n_N2, 0);
    mass_CH4_injected(i) = (n_CH4 * M_CH4 * mol_to_mg) / injection_count;
    flowRate_CH4(i) = mass_CH4_injected(i) / injection_duration;

    % Ideal comparison
    n_N2_ideal = (Pb_abs_bar * bar_to_Pa * V) / (R * T);
    n_total_ideal = (P_total_abs_bar * bar_to_Pa * V) / (R * T);
    n_CH4_ideal = max(n_total_ideal - n_N2_ideal, 0);
    mass_CH4_injected_ideal(i) = (n_CH4_ideal * M_CH4 * mol_to_mg) / injection_count;
    flowRate_CH4_ideal(i) = mass_CH4_injected_ideal(i) / injection_duration;

    % Debug output
    if DEBUG_SRK && (i <= 3 || mod(i,10)==0)
        fprintf('Row %d: T=%.1f K, P_total=%.2f bar, P_back=%.2f bar\n', i, T, P_total_abs_bar, Pb_abs_bar);
        fprintf('  Z_mix=%.4f, Z_N2=%.4f\n', Z_mix, Z_N2);
        fprintf('  n_total(real)=%.6f, n_total(ideal)=%.6f\n', n_total, n_total_ideal);
        fprintf('  n_CH4(real)=%.6f, n_CH4(ideal)=%.6f\n', n_CH4, n_CH4_ideal);
        fprintf('  mass(real)=%.3f mg, mass(ideal)=%.3f mg, ratio=%.4f\n\n', ...
                mass_CH4_injected(i), mass_CH4_injected_ideal(i), mass_CH4_injected(i)/mass_CH4_injected_ideal(i));
    end
end

data.Mass_CH4_Injected_mg = mass_CH4_injected;
data.FlowRate_CH4_mg_per_ms = flowRate_CH4;
data.Mass_CH4_Injected_mg_Ideal = mass_CH4_injected_ideal;
data.FlowRate_CH4_Ideal_mg_per_ms = flowRate_CH4_ideal;

% Statistics
unique_durations = unique(data.InjectorEnergisingTime_ms_);
n_durations = numel(unique_durations);
mean_flow_rate = zeros(n_durations, 1);
std_injected_mass = zeros(n_durations, 1);
mean_injected_mass = zeros(n_durations, 1);
mean_flow_rate_ideal = zeros(n_durations, 1);
std_injected_mass_ideal = zeros(n_durations, 1);
mean_injected_mass_ideal = zeros(n_durations, 1);

for i = 1:n_durations
    rows = data.InjectorEnergisingTime_ms_ == unique_durations(i);
    mean_flow_rate(i) = mean(flowRate_CH4(rows));
    std_injected_mass(i) = std(mass_CH4_injected(rows));
    mean_injected_mass(i) = mean(mass_CH4_injected(rows));
    mean_flow_rate_ideal(i) = mean(flowRate_CH4_ideal(rows));
    std_injected_mass_ideal(i) = std(mass_CH4_injected_ideal(rows));
    mean_injected_mass_ideal(i) = mean(mass_CH4_injected_ideal(rows));
end

summary_table = table(unique_durations, mean_flow_rate, std_injected_mass, mean_injected_mass, ...
    'VariableNames', {'Injection_Duration_ms', 'Mean_FlowRate_mg_per_ms', 'Std_Injected_Mass_mg', 'Mean_Injected_Mass_mg'});
summary_table_ideal = table(unique_durations, mean_flow_rate_ideal, std_injected_mass_ideal, mean_injected_mass_ideal, ...
    'VariableNames', {'Injection_Duration_ms', 'Mean_FlowRate_Ideal_mg_per_ms', 'Std_Injected_Mass_Ideal_mg', 'Mean_Injected_Mass_Ideal_mg'});

injected_mass_table = table(data.InjectorEnergisingTime_ms_, mass_CH4_injected, ...
    'VariableNames', {'Injection_Duration_ms', 'Injected_Mass_mg'});
injected_mass_table_ideal = table(data.InjectorEnergisingTime_ms_, mass_CH4_injected_ideal, ...
    'VariableNames', {'Injection_Duration_ms', 'Injected_Mass_mg'});

% Linear regression (duration >= 0.6 ms)
filtered_idx = summary_table.Injection_Duration_ms >= 0.8;

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
ax_non = ax_all(1);
ax_ideal = ax_all(2);

hold(ax_non, "on");
fill(ax_non, [x_filtered1; flipud(x_filtered1)], [y_upper1; flipud(y_lower1)], [0.82 0.93 0.98], 'EdgeColor', 'none');
plot(ax_non, injected_mass_table.Injection_Duration_ms, injected_mass_table.Injected_Mass_mg, 'o', 'MarkerFaceColor', [0.53 0.74 0.87], 'MarkerEdgeColor', 'none');
plot(ax_non, summary_table.Injection_Duration_ms, summary_table.Mean_Injected_Mass_mg, 'o', 'MarkerFaceColor', [0.05 0.33 0.65], 'MarkerEdgeColor', 'none');
plot(ax_non, x_filtered1, y_fit1, 'Color', [0.05 0.33 0.65], 'LineWidth', 1.5);
xticks(ax_non, [0.6 0.8 1 2 3 4 5 6 8 10]);
xtickangle(ax_non, 45);
xticklabels(ax_non, {});
yticks(ax_non, [0 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30]);
xlim(ax_non, [0.45, 6.25]);
ylim(ax_non, [0, 30.5]);
CLASS_Utilis.InsertFigureText(ax_non, 0.3, 0.85, fitEqnNonIdeal);
title(ax_non, sheetnames, 'Interpreter', 'none');

hold(ax_ideal, "on");
fill(ax_ideal, [x_filtered2; flipud(x_filtered2)], [y_upper2; flipud(y_lower2)], [0.97 0.89 0.82], 'EdgeColor', 'none');
plot(ax_ideal, injected_mass_table_ideal.Injection_Duration_ms, injected_mass_table_ideal.Injected_Mass_mg, 'o', 'MarkerFaceColor', [0.96 0.69 0.25], 'MarkerEdgeColor', 'none');
plot(ax_ideal, summary_table_ideal.Injection_Duration_ms, summary_table_ideal.Mean_Injected_Mass_Ideal_mg, 'o', 'MarkerFaceColor', [0.83 0.33 0], 'MarkerEdgeColor', 'none');
plot(ax_ideal, x_filtered2, y_fit2, 'Color', [0.83 0.33 0], 'LineWidth', 1.5);
xticks(ax_ideal, [0.6 0.8 1 2 3 4 5 6 8 10]);
xtickangle(ax_ideal, 45);
yticks(ax_ideal, [0 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30]);
xlim(ax_ideal, [0.45, 6.25]);
ylim(ax_ideal, [0, 30.5]);
CLASS_Utilis.InsertFigureText(ax_ideal, 0.3, 0.85, fitEqnIdeal);
xlabel(ax_ideal, 'Energizing time [ms]');
ylabel(ax_ideal, 'Injected mass [mg]');

%% SRK calculations now handled by CLASS_SRK_EOS
% No inline functions needed - see Library_Matlab/CLASS_SRK_EOS.m


