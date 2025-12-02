classdef LV_model
    % Generalized Lotka–Volterra (GLV) for 2 species with individual K
    % where A_LV encodes self-limitation and cross-competition terms.

    % This also has the function to embed the LV parameters into the 3x3
    % embed matrix and then fit the Replicator equation to the LV dynamics.

    %% =================== Properties ===================
    properties
        G          % (2x1) Growth rates
        K          % (2x1) Carrying capacities
        B          % (2x1) Tunable Interaction Parameters
        A_LV       % (2x2) Interaction/self-limiting matrix
        tspan      % (1x2 or vector) time span or evaluation grid
        N0         % (scalar) initial total density
        fractions  % (nx1) initial x1 fractions to simulate (x2=1-x1)

        % Simulation outputs
        N_data     % cell{n} each: [T x 2] populations
        X_data     % cell{n} each: [T x 2] fractions
        T_data     % cell{n} each: [T x 1] time grid
        
    end

    %% =================== Methods ===================
    methods
        function obj = LV_model(G,K,B_row,N0,varargin)
            % LV_model(G,K,B_row,N0, [tspan], [fractions])
            % G,K: 2x1 column vectors
            % B_row: [b12  b21], dimensionless cross-competition “pressure”
            % N0: initial total density (cells/mL)
            %
            % A_LV is constructed as: This is the interaction parameter
            %   a11 = -G1/K1 (self), a22 = -G2/K2 (self)
            %   a12 = -G1*b12/K1 (effect of 2 on 1)
            %   a21 = -G2*b21/K2 (effect of 1 on 2)

            obj.G = G(:);
            obj.K = K(:);
            obj.B = B_row.';% Since B_row is a 1x2 row vector

            a11 = -obj.G(1)/obj.K(1);
            a22 = -obj.G(2)/obj.K(2);
            a12 = -obj.G(1)*B_row(1,1)/obj.K(1); % effect of 2 on 1
            a21 = -obj.G(2)*B_row(1,2)/obj.K(2); % effect of 1 on 2
            obj.A_LV = [a11, a12; a21, a22];

            obj.N0 = N0;

            % Optional tspan & fractions
            if ~isempty(varargin)
                obj.tspan = varargin{1};
            else
                % default: evaluate each hour from 0..24
                obj.tspan = 0:1:24;
            end
            if numel(varargin) >= 2 && ~isempty(varargin{2})
                obj.fractions = varargin{2}(:);
            else
                obj.fractions = linspace(0.1,0.9,9).';
            end

            % Preallocate outputs
            nF = numel(obj.fractions);
            obj.N_data = cell(nF,1);
            obj.X_data = cell(nF,1);
            obj.T_data = cell(nF,1);
        end


        function obj = run_lotka_hours(obj, T_end, dt, do_plots)
        % RUN_LOTKA_HOURS  simulate GLV for any horizon and step, with upfront analytical guards.
        % Usage: obj = obj.run_lotka_hours(24, 1, true);
        %        obj = obj.run_lotka_hours(200, 1, false);
        
            if nargin < 2 || isempty(T_end),  T_end = 24; end
            if nargin < 3 || isempty(dt),     dt    = 1;  end
            if nargin < 4,                    do_plots = true; end
        
            % ----------------------------------------------------------
            % 1) ANALYTICAL CHECKS 
            % ----------------------------------------------------------
            K1 = obj.K(1); K2 = obj.K(2);
            G1 = obj.G(1); G2 = obj.G(2);
            a11 = obj.A_LV(1,1); a12 = obj.A_LV(1,2);
            a21 = obj.A_LV(2,1); a22 = obj.A_LV(2,2);
        
            % Infer b12, b21 from your mapping (only for diagnostics/thresholds)
            b12 = -a12 * K1 / G1;
            b21 = -a21 * K2 / G2;
        
            % Invasion thresholds
            th12 = K1 / K2;   % vertical line  b12 = K1/K2
            th21 = K2 / K1;   % horizontal line b21 = K2/K1
            can1 = (b12 < th12);   % species 1 invades (0,K2)
            can2 = (b21 < th21);   % species 2 invades (K1,0)
        
            % Interior candidate
            D  = 1 - b12*b21;
            y1_int = (K1 - b12*K2) / D;
            y2_int = (K2 - b21*K1) / D;
        
            % Hard guard 1: blow-up / non-finite interior
            if D == 0
                error('LV:BadParams', ...
                     'Analytical check failed: D = 1 - b12*b21 = %.4g <= 0 (b12=%.4g, b21=%.4g).', D, b12, b21);
            end
        
            % Hard guard 2: if coexistence is predicted, interior must be positive
            if can1 && can2
                if ~(isfinite(y1_int) && isfinite(y2_int)) || y1_int <= 0 || y2_int <= 0
                    error('LV:BadParams', ...
                         ['Analytical check failed: coexistence predicted (b12<th12 & b21<th21) ', ...
                          'but interior populations are non-positive: [y1*, y2*] = [%.4g, %.4g].'], ...
                          y1_int, y2_int);
                end
            end
        
            % ----------------------------------------------------------
            % 2) SIMULATION 
            % ----------------------------------------------------------
            obj.tspan = 0:dt:T_end;  % evaluation grid
        
            nF = numel(obj.fractions);
            for i = 1:nF
                p1 = obj.fractions(i);
                y0 = [p1; 1-p1] * obj.N0;
        
                % Solve
                [t, Y] = obj.solve_lotka(y0, obj.tspan);
        
                % Fractions
                Nsum = Y(:,1) + Y(:,2);
                x = zeros(size(Y));
                nz = Nsum > 0;
                x(nz,1) = Y(nz,1) ./ Nsum(nz);
                x(nz,2) = Y(nz,2) ./ Nsum(nz);
        
                % Store per-inoculum trajectories
                obj.N_data{i} = Y;
                obj.X_data{i} = x;
                obj.T_data{i} = t;
            end
        
            if do_plots
                obj.plot_populations();
                obj.plot_fractions();
            end
        end


        

        function [t, N] = simulate_logistic(obj, r, K, N0, t_end, dt,col)
        % SIMULATE_LOGISTIC
        % Simulate a 1-species logistic differential equation:
        %   dN/dt = r * N * (1 - N/K)
        %
        % Inputs:
        %   r     – growth rate (1/hour)
        %   K     – carrying capacity (cells/mL)
        %   N0    – initial density (cells/mL)
        %   t_end – end time (hours)
        %   dt    – time step (hours)
        %
        % Outputs:
        %   t     – time vector
        %   N     – population trajectory
        %
        % Also makes a log10(N) vs time figure.
        
            if nargin < 6 || isempty(dt), dt = 0.1; end
            if nargin < 5 || isempty(t_end), t_end = 24; end
        
            % Time grid
            t = 0:dt:t_end;
        
            % Logistic ODE RHS
            rhs = @(t, N) r * N * (1 - N / K);
        
            % Solve ODE
            opts = odeset('RelTol',1e-9,'AbsTol',1e-12,'NonNegative',1);
            [t_out, N_out] = ode45(rhs, t, N0, opts);
        
            % Ensure non-negative
            N = max(N_out, eps);
            t = t_out;
        
            % ----- Plot log10 populations -----
            figure('Color','w');
            plot(t, log10(N), '-', 'LineWidth', 2.5, 'Color', col);
        
            %plot(t, log10(N), col, 'LineWidth', 2);
            xlabel('Time (hours)');
            ylabel('log_{10} N (cells/mL)');
            title(sprintf('1-Species Logistic Growth: r = %.2f, K = %.2g', r, K));
            ylim([5 9.5]);
            grid on;
        
        end



%% Plotting Functions
        function plot_populations(obj)
            % LOG10(populations) vs time for all inocula
            nF = numel(obj.fractions);
            figure('Color','w','Name','LV: log_{10} populations vs time (all inocula)');
            tiledlayout(ceil(nF/3), 3, 'Padding','compact', 'TileSpacing','compact');

            for i = 1:nF
                t = obj.T_data{i};
                Y = obj.N_data{i};
                Ysafe = max(Y, eps);

                nexttile; hold on; box on;
                plot(t, log10(Ysafe(:,1)), '-',  'LineWidth', 1.5, 'DisplayName','Species 1');
                plot(t, log10(Ysafe(:,2)), '--', 'LineWidth', 1.5, 'DisplayName','Species 2');
                xlabel('Time (hours)');
                ylabel('log_{10} N (cells/mL)');
                title(sprintf('Init frac = %.2f / %.2f', obj.fractions(i), 1-obj.fractions(i)));
                grid on;
                if i == 1, legend('Location','best'); end
            end
        end

         function plot_fractions(obj)
            % Fractions vs time for all inocula
            nF = numel(obj.fractions);
            figure('Color','w','Name','LV: Fractions over time');
            tiledlayout(ceil(nF/3), 3, 'Padding','compact', 'TileSpacing','compact');

            for i = 1:nF
                t = obj.T_data{i};
                x = obj.X_data{i};

                nexttile; hold on; box on;
                plot(t, x(:,1), '-',  'LineWidth', 1.7, 'DisplayName','Species 1');
                plot(t, x(:,2), '--', 'LineWidth', 1.7, 'DisplayName','Species 2');
                ylim([0 1]); xlim([t(1) t(end)]);
                title(sprintf('Init frac = %.2f', x(1,1)));
                xlabel('Time (h)'); ylabel('Population fraction');
                if i == 1, legend('Location','best'); end
            end
        end
    end

    %% ======= Helper / private methods =======
    methods (Access = private)
        function [t, Y] = solve_lotka(obj, y0, tspan)
            % Solve GLV using ode45; supports vector tspan (evaluation grid) or [t0 tf]
            if numel(tspan) == 2
                t0 = tspan(1); tf = tspan(2);
                opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
                ode = @(t,y) obj.glv_rhs(t,y);
                [t, Y] = ode45(ode, [t0 tf], y0, opts);
            else
                % fixed evaluation grid
                t_eval = tspan(:);
                t0 = t_eval(1); tf = t_eval(end);
                opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
                ode = @(t,y) obj.glv_rhs(t,y);
                [t_raw, Y_raw] = ode45(ode, [t0 tf], y0, opts);
                % interpolate onto t_eval
                Y = interp1(t_raw, Y_raw, t_eval, 'pchip');
                t = t_eval;
            end

            % Enforce non-negativity (numerical safety)
            Y = max(Y, 0);
        end

        function dydt = glv_rhs(obj, ~, y)
            % GLV RHS: dy/dt = diag(G).*y + [a11 a12; a21 a22] * (y .* y)
            % Derivation for logistic-competition form:
            %   dy1/dt = G1*y1 + a11*y1^2 + a12*y1*y2, with a11=-G1/K1, a12=-G1*b12/K1, etc.
            y1 = y(1); y2 = y(2);
            G1 = obj.G(1); G2 = obj.G(2);
            a11 = obj.A_LV(1,1); a12 = obj.A_LV(1,2);
            a21 = obj.A_LV(2,1); a22 = obj.A_LV(2,2);
            dydt = [
                G1*y1 + a11*y1^2 + a12*y1*y2;
                G2*y2 + a21*y1*y2 + a22*y2^2
            ];
        end
    end
end
