classdef CLASS_SRK_EOS
    % CLASS_SRK_EOS: Soave-Redlich-Kwong equation of state utility
    % Provides static methods for single-component and binary mixture Z calculations
    %
    % Example usage:
    %   % Single component Z calculation
    %   Z_N2 = CLASS_SRK_EOS.CalculateSingleComponentZ(P_Pa, T, Pc_N2, Tc_N2, omega_N2);
    %
    %   % Binary mixture with default parameters
    %   Pc_MPa = [Pc_H2, Pc_N2] / 1e6;
    %   Tc = [Tc_H2, Tc_N2];
    %   omega = [omega_H2, omega_N2];
    %   [Z_mix, x_H2] = CLASS_SRK_EOS.CalculateBinaryMixtureZ(T, P_total_bar, V, P_back_bar, ...
    %                                                          Pc_MPa, Tc, omega);
    %
    %   % Binary mixture with custom tolerance
    %   [Z_mix, x_H2] = CLASS_SRK_EOS.CalculateBinaryMixtureZ(T, P_total_bar, V, P_back_bar, ...
    %                                                          Pc_MPa, Tc, omega, 1e-6);
    %
    %   % Binary mixture with custom options
    %   opts.max_iter = 200;
    %   opts.R_bar = 0.08314;
    %   [Z_mix, x_H2] = CLASS_SRK_EOS.CalculateBinaryMixtureZ(T, P_total_bar, V, P_back_bar, ...
    %                                                          Pc_MPa, Tc, omega, [], opts);
    
    properties (Constant)
        % SRK EOS universal coefficients
        SRK_A_COEFF = 0.42748;  % Attraction parameter coefficient
        SRK_B_COEFF = 0.08664;  % Repulsion parameter coefficient
        SRK_M_C0 = 0.48;        % m parameter constant
        SRK_M_C1 = 1.574;       % m parameter linear coefficient
        SRK_M_C2 = 0.176;       % m parameter quadratic coefficient
        
        % Default physical constants and conversion factors
        DEFAULT_R_BAR = 0.08314;        % L·bar/(mol·K)
        DEFAULT_BAR_TO_MPA = 0.1;       % Conversion: bar to MPa
        DEFAULT_MPA_TO_PA = 1e6;        % Conversion: MPa to Pa
        DEFAULT_M3_TO_L = 1e3;          % Conversion: m^3 to L
        
        % Default numerical parameters
        DEFAULT_TOL = 1e-5;             % Iteration tolerance
        DEFAULT_MAX_ITER = 100;         % Maximum iterations
        DEFAULT_IMAG_TOL = 1e-6;        % Imaginary tolerance for roots
    end
    
    methods (Static)
        function Z = CalculateSingleComponentZ(P, T, P_c, T_c, omega, imag_tol)
            % Calculate compressibility factor Z for a single component using SRK EOS
            % Inputs:
            %   P         - Pressure [Pa]
            %   T         - Temperature [K]
            %   P_c       - Critical pressure [Pa]
            %   T_c       - Critical temperature [K]
            %   omega     - Acentric factor [-]
            %   imag_tol  - Optional: Imaginary tolerance for root selection (default: 1e-6)
            % Output:
            %   Z         - Compressibility factor (largest real root)
            
            % Set default for optional parameter
            if nargin < 6
                imag_tol = CLASS_SRK_EOS.DEFAULT_IMAG_TOL;
            end
            
            R = 8.314; % Universal gas constant, J/(mol·K)
            Tr = T / T_c; % Reduced temperature
            
            % SRK temperature-dependent parameter
            m = CLASS_SRK_EOS.SRK_M_C0 + CLASS_SRK_EOS.SRK_M_C1 * omega - CLASS_SRK_EOS.SRK_M_C2 * omega^2;
            alpha = (1 + m * (1 - sqrt(Tr)))^2;
            
            % SRK attraction and repulsion parameters
            a = CLASS_SRK_EOS.SRK_A_COEFF * (R^2 * T_c^2) / P_c * alpha;
            b = CLASS_SRK_EOS.SRK_B_COEFF * (R * T_c) / P_c;
            
            % Dimensionless parameters
            A = (a * P) / (R^2 * T^2);
            B = (b * P) / (R * T);
            
            % SRK cubic equation: Z^3 - Z^2 + (A - B - B^2)Z - AB = 0
            coeffs = [1, -1, A - B - B^2, -A * B];
            roots_Z = roots(coeffs);
            
            % Select largest real root (vapor phase)
            Z = max(real(roots_Z(abs(imag(roots_Z)) < imag_tol)));
        end
        
        function [Z_final, x_inj_final] = CalculateBinaryMixtureZ(T, P_total_bar, V_m3, P_back_bar, Pc_MPa, Tc, omega, tol, options)
            % Iteratively solve Z and injected gas mole fraction for binary mixture using SRK EOS
            % Inputs:
            %   T            - Temperature [K]
            %   P_total_bar  - Final absolute pressure after injection [bar]
            %   V_m3         - Chamber volume [m^3]
            %   P_back_bar   - Initial absolute back pressure from background gas [bar]
            %   Pc_MPa       - Critical pressures [MPa] (vector: [injected, background])
            %   Tc           - Critical temperatures [K] (vector: [injected, background])
            %   omega        - Acentric factors (vector: [injected, background])
            %   tol          - Optional: Iteration tolerance (default: 1e-5)
            %   options      - Optional: Parameters struct with fields:
            %                  .R_bar (default: 0.08314 L·bar/mol·K)
            %                  .bar_to_MPa (default: 0.1)
            %                  .MPa_to_Pa (default: 1e6)
            %                  .m3_to_L (default: 1e3)
            %                  .max_iter (default: 100)
            %                  .imag_tol (default: 1e-6)
            % Outputs:
            %   Z_final      - Final compressibility factor Z
            %   x_inj_final  - Final mole fraction of injected gas in mixture
            
            % Set default for tol
            if nargin < 8 || isempty(tol)
                tol = CLASS_SRK_EOS.DEFAULT_TOL;
            end
            
            % Set default values for optional parameters
            if nargin < 9
                options = struct();
            end
            if ~isfield(options, 'R_bar'), options.R_bar = CLASS_SRK_EOS.DEFAULT_R_BAR; end
            if ~isfield(options, 'bar_to_MPa'), options.bar_to_MPa = CLASS_SRK_EOS.DEFAULT_BAR_TO_MPA; end
            if ~isfield(options, 'MPa_to_Pa'), options.MPa_to_Pa = CLASS_SRK_EOS.DEFAULT_MPA_TO_PA; end
            if ~isfield(options, 'm3_to_L'), options.m3_to_L = CLASS_SRK_EOS.DEFAULT_M3_TO_L; end
            if ~isfield(options, 'max_iter'), options.max_iter = CLASS_SRK_EOS.DEFAULT_MAX_ITER; end
            if ~isfield(options, 'imag_tol'), options.imag_tol = CLASS_SRK_EOS.DEFAULT_IMAG_TOL; end
            
            R_bar = options.R_bar;
            bar_to_MPa = options.bar_to_MPa;
            MPa_to_Pa = options.MPa_to_Pa;
            m3_to_L = options.m3_to_L;
            max_iter = options.max_iter;
            imag_tol = options.imag_tol;
            
            % Convert pressures to MPa
            P_total = P_total_bar * bar_to_MPa;
            P_back = P_back_bar * bar_to_MPa;
            
            % Calculate initial background gas moles (pure background at back pressure) using SRK
            i_bg = 2; % Background gas index
            Z_bg = CLASS_SRK_EOS.CalculateSingleComponentZ(P_back * MPa_to_Pa, T, Pc_MPa(i_bg) * MPa_to_Pa, Tc(i_bg), omega(i_bg));
            n_bg_initial = (P_back_bar * V_m3 * m3_to_L) / (Z_bg * R_bar * T); % Fixed background moles
            
            % Initial guess for injected gas mole fraction
            x_inj = (P_total_bar - P_back_bar) / P_total_bar;
            
            % Iterative solution for mixture Z and composition
            for iter = 1:max_iter
                y = [x_inj, 1 - x_inj]; % [y_injected, y_background]
                a_i = zeros(1,2);
                b_i = zeros(1,2);
                
                % Calculate component-specific SRK parameters
                for i = 1:2
                    Tr = T / Tc(i);
                    m = CLASS_SRK_EOS.SRK_M_C0 + CLASS_SRK_EOS.SRK_M_C1 * omega(i) - CLASS_SRK_EOS.SRK_M_C2 * omega(i)^2;
                    alpha = (1 + m * (1 - sqrt(Tr)))^2;
                    a_i(i) = CLASS_SRK_EOS.SRK_A_COEFF * R_bar^2 * Tc(i)^2 / Pc_MPa(i) * alpha;
                    b_i(i) = CLASS_SRK_EOS.SRK_B_COEFF * R_bar * Tc(i) / Pc_MPa(i);
                end
                
                % Mixing rules (van der Waals, no binary interaction parameter)
                a_mix = 0;
                for i = 1:2
                    for j = 1:2
                        a_mix = a_mix + y(i) * y(j) * sqrt(a_i(i) * a_i(j));
                    end
                end
                b_mix = sum(y .* b_i);
                
                % Dimensionless mixture parameters
                A = a_mix * P_total / (R_bar^2 * T^2);
                B = b_mix * P_total / (R_bar * T);
                
                % SRK cubic equation: Z^3 - Z^2 + (A - B - B^2)Z - AB = 0
                coeffs = [1, -1, A - B - B^2, -A * B];
                roots_Z = roots(coeffs);
                Z_real = real(roots_Z(abs(imag(roots_Z)) < imag_tol));
                Z = max(Z_real);
                
                % Calculate total moles and update injected gas mole fraction
                n_total = (P_total_bar * V_m3 * m3_to_L) / (Z * R_bar * T);
                x_inj_new = (n_total - n_bg_initial) / n_total;
                
                % Check convergence
                if abs(x_inj_new - x_inj) < tol
                    Z_final = Z;
                    x_inj_final = x_inj_new;
                    return;
                end
                x_inj = x_inj_new;
            end
            
            error('CLASS_SRK_EOS:NoConvergence', 'SRK mixture iteration did not converge after %d iterations.', max_iter);
        end
    end
end
