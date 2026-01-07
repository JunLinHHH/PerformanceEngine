% Peng-Robinson EOS parameters
% EOS parameters
clear; clc; close all

% -----------------------------------------
saveName = 'C:\Users\z5058464\OneDrive - UNSW\CH4DDI\Engine\MassflowRate\Trial';
% -----------------------------------------


% Given constants
R = 8.314; % Universal gas constant, J/(mol·K)
V = 0.0014; % Chamber volume, m^3
M_CH4 = 16.04; % Molar mass of methane, g/mol

% Load data from Excel
fileName = 'C:\Users\z5058464\OneDrive - UNSW\CH4DDI\Engine\MassflowRate\20251215_OldChamber_FlowRate.xlsx'; % Replace with the name of your Excel file
sheetnames = 'CH4_350_ST02_NonIdeal_B0_T2';
data = readtable(fileName,'Sheet',sheetnames);

% Initialize results
mass_CH4_injected = zeros(height(data), 1); % Mass of methane injected (mg)
flowRate_CH4 = zeros(height(data), 1); % Flow rate (mg/ms)

for i = 1:height(data)
    % Extract data for the current row
    T = data.InitialTemperature_K_T_1(i); % Temperature (K)
    P_diff = data.PressureRiseDueToH2Injection_bar_P_diff(i)*1e5; % Pressure rise due to CH4 injection
    Pb = data.BackPressure_bar_P_b(i) * 1e5; % back pressure (Pa)
    P_total = P_diff + Pb; % pressure after inejection
    injection_duration = data.InjectorEnergisingTime_ms_(i); % Injection duration (ms)
    injection_count = data.InjectionQuantity(i);
    
    % Critical properties for CH4
    P_c_CH4 = 4.599e6; % Critical pressure (Pa)
    T_c_CH4 = 190.56; % Critical temperature (K)
    omega_CH4 = 0.01142; % Acentric factor
    
    % Final moles in the chamber
    Z = calculateZCH4(P_total, T, V, P_c_CH4, T_c_CH4, omega_CH4); % Compressibility factor at P2
    n = (P_diff * V) / (Z * R * T); % Final moles
    
    % Calculate injected methane mass
    mass_CH4_injected(i) = (n * M_CH4 * 1000)/injection_count; % Convert to mg
    
    % Calculate flow rate
    flowRate_CH4(i) = mass_CH4_injected(i) / injection_duration; % Flow rate in mg/ms
end

% Append results to the table and save
data.Mass_CH4_Injected_mg = mass_CH4_injected;
data.FlowRate_CH4_mg_per_ms = flowRate_CH4;

% Save the updated data table
writetable(data, [saveName '_processed.xlsx']);

% Calculate mean and standard deviation of flow rates for each injection duration
unique_durations = unique(data.InjectorEnergisingTime_ms_); % Unique injection durations
mean_flow_rate = zeros(length(unique_durations), 1); % Mean flow rates
std_flow_rate = zeros(length(unique_durations), 1); % Standard deviations

for i = 1:length(unique_durations)
    % Get rows with the current injection duration
    current_duration = unique_durations(i);
    rows = data.InjectorEnergisingTime_ms_ == current_duration;
    
    % Calculate mean and standard deviation for the current group
    mean_flow_rate(i) = mean(flowRate_CH4(rows));
    std_flow_rate(i) = std(flowRate_CH4(rows));
    std_injected_mass(i,1) = std(mass_CH4_injected(rows));
    mean_injected_mass(i,1) = mean(mass_CH4_injected(rows)); % Calculate mean injected mass
end

% Create summary table for mean and std results
summary_table = table(unique_durations, mean_flow_rate, std_flow_rate, std_injected_mass, mean_injected_mass, ...
    'VariableNames', {'Injection_Duration_ms', 'Mean_FlowRate_mg_per_ms', 'Std_FlowRate_mg_per_ms', 'Std_Injected_Mass_mg', 'Mean_Injected_Mass_mg'});

% Create individual mass table for each injection
injected_mass_table = table(data.InjectorEnergisingTime_ms_, mass_CH4_injected, ...
    'VariableNames', {'Injection_Duration_ms', 'Injected_Mass_mg'});

% Line of best fit
% 1. When injection duration >= 1
x1 = summary_table.Injection_Duration_ms;
y1 = summary_table.Mean_Injected_Mass_mg;
std1 = summary_table.Std_Injected_Mass_mg;

filtered_idx1 = x1 >= 0.8;
x_filtered1 = x1(filtered_idx1);
y_filtered1 = y1(filtered_idx1);
std_filtered1 = std1(filtered_idx1);
coefficients1 = polyfit(x_filtered1, y_filtered1, 1); % Linear fit (1st order polynomial)
y_fit1 = polyval(coefficients1, x_filtered1); % Evaluate the fit

SS_res1 = sum((y_filtered1 - y_fit1).^2); % Residual sum of squares
SS_tot1 = sum((y_filtered1 - mean(y_filtered1)).^2); % Total sum of squares
R_squared1 = 1 - (SS_res1 / SS_tot1);

y_upper1 = y_fit1 + std_filtered1;
y_lower1 = y_fit1 - std_filtered1;

fitEqn = sprintf('y = %.3fx + %.3f\nR^2 = %.4f', coefficients1(1), coefficients1(2), R_squared1);
%% Plot
% Create a new figure using custom axes class
obj_ax = CLASS_AxesHandleStore;
obj_ax.MarginLeft = 45;
obj_ax.MarginRight = 30;
obj_ax.MarginBottom = 50;
obj_ax.MarginTop = 30;
obj_ax.AxesWidth = 500; % Example value, adjust as needed
obj_ax.AxesHeight = 500; % Example value, adjust as needed
obj_ax.RowNumber = 1;
obj_ax.ColumnNumber = 1;
GapRow = ones(1, obj_ax.RowNumber) * 5;
GapCol = ones(1,obj_ax.ColumnNumber) * 5;
obj_ax.GapRow = GapRow;
obj_ax.GapColumn = GapCol;
obj_ax = obj_ax.ConstuctAxesUnevenGap();
handle_ax = obj_ax.GetAxesHandleMatrix();

hold(handle_ax,"on");
fill(handle_ax,[x_filtered1; flipud(x_filtered1)], [y_upper1; flipud(y_lower1)], ...
    [0.9, 0.9, 0.9], 'EdgeColor', 'none','FaceColor', '#d2edf9');
plot(handle_ax,injected_mass_table.Injection_Duration_ms,injected_mass_table.Injected_Mass_mg,'o','MarkerFaceColor','#87BDDD','MarkerEdgeColor','none');
plot(handle_ax,summary_table.Injection_Duration_ms,summary_table.Mean_Injected_Mass_mg,'o','MarkerFaceColor','#0C54A5','MarkerEdgeColor','none');
plot(x_filtered1,y_fit1,'Color','#0C54A5','LineWidth',1.5);
CLASS_Utilis.InsertFigureText(handle_ax,0.3,0.78,fitEqn);
xticks([0.6 0.8 1 2 3 4 5 6 8 10]);
yticks([0 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40]);
xlim([0.45,6.25])
title(handle_ax, sheetnames,'Interpreter','none')
xlabel('Energizing time [ms]')
ylabel('Injected mass [mg]')

%% Function
% Input
% Pressure (Pa), Temp (K), volume (m^3), Critical pressure (Pa), Critical Temp
% (k), Acentric factor
function Z = calculateZCH4(P, T, ~, P_c, T_c, omega)
    % Peng-Robinson parameters
    R = 8.314; % Universal gas constant, J/(mol·K)
    Tr = T / T_c; % Reduced temperature
    Pr = P / P_c; % Reduced pressure
    a = 0.45724 * (R^2 * T_c^2) / P_c; % Attraction parameter
    b = 0.07780 * (R * T_c) / P_c; % Repulsion parameter
    alpha = (1 + (0.480 + 1.574 * omega - 0.176 * omega^2) * (1 - sqrt(Tr)))^2; % Temperature-dependent factor
    A = (a * alpha * P) / (R^2 * T^2); % Dimensionless A
    B = (b * P) / (R * T); % Dimensionless B
    
    % Cubic equation coefficients
     coeffs = [1, -(1 - B), A - 2 * B - 3 * B^2, -(A * B - B^2 - B^3)]; % PR EOS
%     coeffs = [1, -(1 - B), A - B - B^2, -A * B]; % SRK EOS

    % Solve cubic equation for Z
    roots_Z = roots(coeffs);
    Z = max(real(roots_Z)); % Select the largest real root
end
