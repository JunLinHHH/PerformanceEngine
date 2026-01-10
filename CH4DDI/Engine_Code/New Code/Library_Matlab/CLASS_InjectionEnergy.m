% filepath: g:\Fork\PerformanceEngine\PerformanceEngine\Library_Matlab\CLASS_EnergyLinearMap.m
classdef CLASS_InjectionEnergy
    % Maps fuel injection rate to cumulative energy input over injection duration.
    %
    % Purpose:
    %   Given total energy required, fuel lower heating value (LHV), and a linear
    %   injection rate model (mass_rate = a*t + b), compute a new linear function
    %   that relates injection time to cumulative energy as a percentage of total.
    %
    % Usage Example:
    %   % Define injection rate: mass_rate(t) = 5*t + 10 [mg/ms]
    %   mapper = CLASS_EnergyLinearMap(100, 43.5, 5, 10);
    %   % Returns linear fit: energy_pct(t) = slope*t + intercept
    %   [slope, intercept, t_100pct] = mapper.GetEnergyLinearFunction();
    %
    % Properties (constant, no external modification needed):
    %   E_total_J       - Total energy required [J]
    %   LHV_fuel_MJ_kg  - Fuel lower heating value [MJ/kg]
    %   a_rate          - Slope of injection rate [mg/ms²]
    %   b_rate          - Intercept of injection rate [mg/ms]
    %
    % Methods:
    %   GetEnergyLinearFunction - Returns linear fit slope, intercept, and time to 100%
    %   GetCumulativeEnergy     - Query cumulative energy at any time t
    %   GetEnergyPercent        - Query energy percentage (0–100%) at time t
    %
    % Notes:
    %   - Injection rate units: mg/ms → convert to kg/s for LHV energy calculation
    %   - Cumulative energy is quadratic in time; linear fit is best-fit approximation
    %   - Time to reach 100% energy is solved implicitly from total mass requirement
    
    properties (Constant)
        % Physical constants (no modification)
    end
    
    properties
        E_total_J       % Total energy required [J]
        LHV_fuel_MJ_kg  % Fuel lower heating value [MJ/kg] (e.g., 43.5 for diesel)
        a_rate          % Slope: dm/dt = a*t + b [mg/ms²]
        b_rate          % Intercept [mg/ms]
    end
    
    methods
        function obj = CLASS_EnergyLinearMap(E_total_J, LHV_fuel_MJ_kg, a_rate, b_rate)
            % Constructor: Initialize energy mapper with injection parameters.
            %
            % Inputs:
            %   E_total_J       - Total energy required [J]
            %   LHV_fuel_MJ_kg  - Fuel lower heating value [MJ/kg]
            %   a_rate          - Slope of injection rate dm/dt = a*t + b [mg/ms²]
            %   b_rate          - Intercept of injection rate [mg/ms]
            
            obj.E_total_J = E_total_J;
            obj.LHV_fuel_MJ_kg = LHV_fuel_MJ_kg;
            obj.a_rate = a_rate;
            obj.b_rate = b_rate;
        end
        
        function [slope, intercept, t_100pct] = GetEnergyLinearFunction(obj)
            % Returns best-fit linear function: energy_pct(t) = slope*t + intercept
            %
            % Physical model:
            %   - Injection rate: m_dot(t) = a*t + b [mg/ms]
            %   - Cumulative mass: M(t) = 0.5*a*t² + b*t [mg]
            %   - Cumulative energy: E(t) = M(t) * LHV = (0.5*a*t² + b*t) * LHV * 1e-6 [J]
            %     (where 1e-6 converts mg→kg and MJ→J: 1e-3 * 1e6)
            %   - Energy %: E_pct(t) = 100 * E(t) / E_total
            %
            % For linear fit, sample energy % at t=0, t=0.5*t_100%, t=t_100%
            % and fit: E_pct = slope*t + intercept
            
            % Step 1: Compute time to inject all required fuel
            % Total mass needed: M_total = E_total / (LHV * 1e6) [kg] = E_total / (LHV * 1e6) * 1e6 [mg]
            M_total_mg = obj.E_total_J / (obj.LHV_fuel_MJ_kg * 1e6) * 1e6;
            
            % Solve: 0.5*a*t² + b*t = M_total
            % Quadratic: 0.5*a*t² + b*t - M_total = 0
            if abs(obj.a_rate) < 1e-12
                % Linear case: b*t = M_total
                if abs(obj.b_rate) < 1e-12
                    error('CLASS_EnergyLinearMap: injection rate is zero (a≈0, b≈0)');
                end
                t_100pct = M_total_mg / obj.b_rate;
            else
                % Quadratic case
                A_coeff = 0.5 * obj.a_rate;
                B_coeff = obj.b_rate;
                C_coeff = -M_total_mg;
                discriminant = B_coeff^2 - 4*A_coeff*C_coeff;
                if discriminant < 0
                    error('CLASS_EnergyLinearMap: negative discriminant (injection rate insufficient)');
                end
                t_sol1 = (-B_coeff + sqrt(discriminant)) / (2*A_coeff);
                t_sol2 = (-B_coeff - sqrt(discriminant)) / (2*A_coeff);
                t_100pct = max(t_sol1, t_sol2); % Physical solution (positive time)
            end
            
            % Step 2: Sample energy % at three points for linear fit
            t_samples = [0, 0.5*t_100pct, t_100pct];
            E_pct_samples = zeros(1, 3);
            
            for i = 1:3
                t = t_samples(i);
                E_pct_samples(i) = obj.GetEnergyPercent(t);
            end
            
            % Step 3: Least-squares linear fit E_pct = slope*t + intercept
            % Use all three points; MATLAB's polyfit (degree 1) gives best fit
            p = polyfit(t_samples, E_pct_samples, 1);
            slope = p(1);      % Energy % per ms
            intercept = p(2);  % Energy % at t=0
        end
        
        function E_cum_J = GetCumulativeEnergy(obj, t)
            % Returns cumulative energy injected at time t.
            %
            % Input:
            %   t - Time [ms]
            %
            % Output:
            %   E_cum_J - Cumulative energy [J]
            
            % Cumulative mass: M(t) = 0.5*a*t² + b*t [mg]
            M_cum_mg = 0.5 * obj.a_rate * t^2 + obj.b_rate * t;
            
            % Convert to kg and apply LHV
            M_cum_kg = M_cum_mg * 1e-6;
            E_cum_J = M_cum_kg * obj.LHV_fuel_MJ_kg * 1e6; % [J]
        end
        
        function E_pct = GetEnergyPercent(obj, t)
            % Returns cumulative energy as a percentage of total (0–100%).
            %
            % Input:
            %   t - Time [ms]
            %
            % Output:
            %   E_pct - Energy percentage [0, 100]
            
            E_cum = obj.GetCumulativeEnergy(t);
            E_pct = 100 * E_cum / obj.E_total_J;
        end
    end
end