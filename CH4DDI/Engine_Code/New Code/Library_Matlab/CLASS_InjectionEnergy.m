classdef CLASS_InjectionEnergy
    % CLASS_InjectionEnergy (static methods only)
    % Utilities to map a linear injection-rate model (cumulative mass: M(t) = 0.5*a*t^2 + b*t + c)
    % to cumulative energy and percent-of-total energy, with crank angle conversion.
    %
    % Conventions / units expected:
    %  - a: acceleration coefficient of cumulative mass [mg/ms^2]
    %  - b: velocity coefficient of cumulative mass [mg/ms]
    %  - c: initial cumulative mass offset [mg] (optional, default 0)
    %  - t: time in [us] (microseconds) for public APIs; internally converted to [ms]
    %  - LHV_fuel_MJ_kg: lower heating value in [MJ/kg]
    %  - E_total_J: total target energy in [J]
    %  - rpm: engine speed in [RPM] (default 1400)
    %
    % Public static methods:
    %  - GetEnergyLinearFunction(E_total_J, LHV_MJ_kg, a, b, c)
    %      -> [slope_pct_per_us, intercept_pct, t_100pct_us]
    %      Returns linear approximation of energy percent vs time (us)
    %
    %  - GetEnergyLinearCoefficients(E_total_J, LHV_MJ_kg, a, b, c)
    %      -> [slope_J_per_us, intercept_J, t_100pct_us]
    %      Returns linear coefficients for Energy [J] vs time [us]
    %
    %  - GetCumulativeEnergy(t_us, LHV_MJ_kg, a, b, c)
    %      -> E_cum_J (vectorized)
    %      Computes cumulative energy at specified time(s)
    %
    %  - GetEnergyPercent(t_us, E_total_J, LHV_MJ_kg, a, b, c)
    %      -> percent [0..100] (vectorized)
    %      Computes percent of total energy reached at specified time(s)
    %
    %  - GetCrankAngle(t_us, rpm)
    %      -> CA_deg (vectorized)
    %      Converts time [us] to crank angle [degrees] at given RPM
    %
    %  - GetPercentTable(E_total_J, LHV_MJ_kg, a, b, c, step_percent, rpm)
    %      -> table with columns: Pct, Energy_J, Time_us, Time_ms, CA_deg
    %      Generates table of energy milestones at regular percent intervals

    methods (Static)
        function [slope, intercept, t_100pct_us] = GetEnergyLinearFunction(E_total_J, LHV_fuel_MJ_kg, a, b, c)
            % Compute a simple linear approximation of cumulative energy percent vs time.
            %
            % Inputs (scalars):
            %  E_total_J        - total energy target [J], must be > 0
            %  LHV_fuel_MJ_kg   - fuel LHV [MJ/kg], must be > 0
            %  a                - slope of mass rate [mg/ms^2]
            %  b                - intercept of mass rate [mg/ms]
            %
            % Outputs:
            %  slope            - percent-per-us (linear fit on microsecond axis)
            %  intercept        - percent at t=0
            %  t_100pct_us      - time [us] at which cumulative energy reaches 100%

            % Provide default total energy if not supplied
            if nargin < 1 || isempty(E_total_J)
                E_total_J = 1400; % default total energy [J]
            end

            % Input validation
            validateattributes(E_total_J, {'numeric'}, {'scalar','real','positive'}, mfilename, 'E_total_J');
            validateattributes(LHV_fuel_MJ_kg, {'numeric'}, {'scalar','real','positive'}, mfilename, 'LHV_fuel_MJ_kg');
            validateattributes(a, {'numeric'}, {'scalar','real','finite','nonnan'}, mfilename, 'a');
            validateattributes(b, {'numeric'}, {'scalar','real','finite','nonnan'}, mfilename, 'b');
            if nargin < 5 || isempty(c)
                c = 0; % initial cumulative mass at t=0 [mg]
            end
            validateattributes(c, {'numeric'}, {'scalar','real','finite','nonnan'}, mfilename, 'c');

            % Convert total energy to total mass required, expressed in mg.
            % Derivation: mass_kg = E_total_J / (LHV_J_per_kg), LHV_J_per_kg = LHV_MJ_kg * 1e6
            % mass_mg = mass_kg * 1e6 -> mass_mg = E_total_J / LHV_MJ_kg
            M_total_mg = E_total_J / LHV_fuel_MJ_kg; % [mg]

            % Solve 0.5*a*t^2 + b*t + c = M_total_mg for t >= 0
            if abs(a) < 1e-15
                % effectively linear mass rate: b*t = M_total_mg
                if abs(b) < 1e-15
                    error('CLASS_InjectionEnergy:ZeroRate','Both a and b are zero; injection rate is zero.');
                end
                t_100pct = (M_total_mg - c) / b;
            else
                A = 0.5 * a;
                B = b;
                C = c - M_total_mg;
                disc = B.^2 - 4*A.*C;
                if disc < 0
                    error('CLASS_InjectionEnergy:InsufficientRate','Negative discriminant: injection parameters cannot reach required mass.');
                end
                t1 = (-B + sqrt(disc)) / (2*A);
                t2 = (-B - sqrt(disc)) / (2*A);
                t_100pct = max(t1, t2);
                if t_100pct < 0
                    error('CLASS_InjectionEnergy:NegativeTime','Computed positive root is negative; check parameters.');
                end
            end

            % Return time in microseconds; internal calculation used ms.
            t_100pct_us = t_100pct * 1e3; % ms -> us

            % Sample percent at three points (times in us) and fit linear model
            t_samples = [0, 0.5*t_100pct_us, t_100pct_us];
            E_pct_samples = CLASS_InjectionEnergy.GetEnergyPercent(t_samples, E_total_J, LHV_fuel_MJ_kg, a, b, c);

            p = polyfit(t_samples(:), E_pct_samples(:), 1);
            slope = p(1);    % percent per us
            intercept = p(2);

            % Numerical guard: clamp very small negative intercepts to zero
            if intercept < 0 && intercept > -1e-9
                intercept = 0;
            end
        end

        function [slope_J_per_us, intercept_J, t_100pct_us] = GetEnergyLinearCoefficients(E_total_J, LHV_fuel_MJ_kg, a, b, c)
            % Return linear coefficients for Energy (J) vs time (us) using
            % the internal percent-based linearization.
            if nargin < 5 || isempty(c)
                c = 0;
            end
            [slope_pct_per_us, intercept_pct, t_100pct_us] = CLASS_InjectionEnergy.GetEnergyLinearFunction(E_total_J, LHV_fuel_MJ_kg, a, b, c);
            % Convert percent-per-us to J per us
            slope_J_per_us = (slope_pct_per_us ./ 100) .* E_total_J;
            intercept_J = (intercept_pct ./ 100) .* E_total_J;
        end

        function E_cum_J = GetCumulativeEnergy(t_us, LHV_fuel_MJ_kg, a, b, c)
            % Vectorized cumulative energy in Joules for times t (microseconds).
            %
            % t_us can be scalar or vector. Returns same-shaped E_cum_J.

            validateattributes(t_us, {'numeric'}, {'real','finite'}, mfilename, 't_us');
            validateattributes(LHV_fuel_MJ_kg, {'numeric'}, {'scalar','real','positive'}, mfilename, 'LHV_fuel_MJ_kg');
            validateattributes(a, {'numeric'}, {'scalar','real','finite','nonnan'}, mfilename, 'a');
            validateattributes(b, {'numeric'}, {'scalar','real','finite','nonnan'}, mfilename, 'b');

            % convert input microseconds to ms because a and b are in mg/ms units
            t_ms = double(t_us) * 1e-3;
            % cumulative mass in mg: M(t_ms) = 0.5*a*t_ms.^2 + b.*t_ms + c
            if nargin < 5 || isempty(c)
                c = 0;
            end
            M_cum_mg = 0.5 * a .* (t_ms.^2) + b .* t_ms + c;
            % guard small negative rounding errors
            M_cum_mg = max(M_cum_mg, 0);

            % Convert mg -> kg and apply LHV (MJ/kg -> J/kg)
            M_cum_kg = M_cum_mg * 1e-6;
            LHV_J_per_kg = LHV_fuel_MJ_kg * 1e6;
            E_cum_J = M_cum_kg .* LHV_J_per_kg;
        end

        function E_pct = GetEnergyPercent(t_us, E_total_J, LHV_fuel_MJ_kg, a, b, c)
            % Vectorized percent of total energy reached at times t (microseconds).
            % Usage: E_pct = CLASS_InjectionEnergy.GetEnergyPercent(t_us, E_total_J, LHV, a, b)

            % default total energy when not provided
            if nargin < 2 || isempty(E_total_J)
                E_total_J = 1400; % default total energy [J]
            end
            validateattributes(E_total_J, {'numeric'}, {'scalar','real','positive'}, mfilename, 'E_total_J');
            if nargin < 6 || isempty(c)
                c = 0;
            end
            E_cum = CLASS_InjectionEnergy.GetCumulativeEnergy(t_us, LHV_fuel_MJ_kg, a, b, c);
            E_pct = 100 .* (E_cum ./ E_total_J);
            % clamp to [0,100] with tolerance
            E_pct = min(max(E_pct, 0), 100);
        end

        function CA_deg = GetCrankAngle(t_us, rpm)
            % GetCrankAngle  Convert time (microseconds) to crank angle (degrees) for given rpm.
            %
            % Usage:
            %   CA_deg = CLASS_InjectionEnergy.GetCrankAngle(t_us)
            %   CA_deg = CLASS_InjectionEnergy.GetCrankAngle(t_us, rpm)
            %
            % Inputs:
            %   t_us   - time in microseconds (scalar or vector)
            %   rpm - engine speed in RPM (optional, default 1400)
            %
            % Output:
            %   CA_deg - crank angle in degrees (same shape as t_us)

            if nargin < 2 || isempty(rpm)
                rpm = 1400; % default RPM
            end
            validateattributes(t_us, {'numeric'}, {'real','finite'}, mfilename, 't_us');
            validateattributes(rpm, {'numeric'}, {'scalar','real','finite','positive'}, mfilename, 'rpm');
            % CA [deg] = 360 * revolutions = 360 * (rpm/60) * t_seconds
            % t provided in microseconds -> t_seconds = t_us * 1e-6
            % Therefore CA = 360*(rpm/60)*t_us*1e-6 = 6e-6 * rpm * t_us
            CA_deg = 6e-6 .* double(rpm) .* double(t_us);
        end

        function tbl = GetPercentTable(E_total_J, LHV_fuel_MJ_kg, a, b, c, step_percent, rpm)
            % GetPercentTable  Build a table of percent steps vs energy/time/CA.
            %
            % Returns a table with columns: Pct, Energy_J, Time_us, Time_ms, CA_deg
            if nargin < 1 || isempty(E_total_J)
                E_total_J = 1400;
            end
            if nargin < 5 || isempty(c)
                c = 0;
            end
            if nargin < 6 || isempty(step_percent)
                step_percent = 10;
            end
            if nargin < 7 || isempty(rpm)
                rpm = 1400;
            end

            validateattributes(E_total_J, {'numeric'}, {'scalar','real','positive'}, mfilename, 'E_total_J');
            validateattributes(LHV_fuel_MJ_kg, {'numeric'}, {'scalar','real','positive'}, mfilename, 'LHV_fuel_MJ_kg');
            validateattributes(a, {'numeric'}, {'scalar','real','finite','nonnan'}, mfilename, 'a');
            validateattributes(b, {'numeric'}, {'scalar','real','finite','nonnan'}, mfilename, 'b');
            validateattributes(step_percent, {'numeric'}, {'scalar','integer','>=',1,'<=',100}, mfilename, 'step_percent');
            validateattributes(rpm, {'numeric'}, {'scalar','real','finite','positive'}, mfilename, 'rpm');
            validateattributes(c, {'numeric'}, {'scalar','real','finite','nonnan'}, mfilename, 'c');

            pcts = 0:step_percent:100;
            n = numel(pcts);
            Energy_J = (pcts/100) .* E_total_J;
            Time_us = nan(size(Energy_J));
            Time_ms = nan(size(Energy_J));
            CA_deg = nan(size(Energy_J));

            for ii = 1:n
                E_i = Energy_J(ii);
                M_i_mg = E_i / LHV_fuel_MJ_kg; % mg (see class convention)

                if abs(a) < 1e-15
                    if abs(b) < 1e-15
                        Time_ms(ii) = NaN;
                    else
                        Time_ms(ii) = (M_i_mg - c) / b;
                    end
                else
                    A = 0.5 * a;
                    B = b;
                    C = c - M_i_mg;
                    disc = B.^2 - 4*A.*C;
                    if disc < 0
                        Time_ms(ii) = NaN;
                    else
                        t1 = (-B + sqrt(disc)) / (2*A);
                        t2 = (-B - sqrt(disc)) / (2*A);
                        tpos = max(t1,t2);
                        if tpos < 0
                            Time_ms(ii) = NaN;
                        else
                            Time_ms(ii) = tpos;
                        end
                    end
                end

                Time_us(ii) = Time_ms(ii) * 1e3;
                CA_deg(ii) = CLASS_InjectionEnergy.GetCrankAngle(Time_us(ii), rpm);
            end

            tbl = table(pcts(:), Energy_J(:), Time_us(:), Time_ms(:), CA_deg(:), ...
                'VariableNames', {'Pct','Energy_J','Time_us','Time_ms','CA_deg'});
        end
    end
end