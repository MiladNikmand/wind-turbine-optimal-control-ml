function p = params(overrides)
%DDP.PARAMS  Every constant the DDP-HJB controller uses, in one place.
%
%   p = ddp.params()
%   p = ddp.params(struct('max_iter', 20, 'Cp_opt', 0.42))
%
% In Main_Control_Script these were scattered through the segment loop and
% redefined on every pass. Nothing here depends on the segment, so they are
% hoisted out. Values are EXACTLY as the thesis code had them.
%
% All of these are adjustable. Change them here, not in the solver.

%% ---------------- turbine plant ----------------
p.R    = 4.4;      % rotor radius (m)
p.Kt   = 52;       % drivetrain damping
p.Jt   = 16;       % total inertia (single mass)
p.Rs   = 0.01;     % generator stator resistance
p.Ls   = 0.01;     % generator stator inductance
p.Rho  = 1.25;     % air density (kg/m^3)

%% ---------------- optimal operating point ----------------
% NOTE: the Cp model used by the predictors and by project 2,
%   Cp = 0.22*(116/lambda - 5)*exp(-12.5/lambda),
% peaks at lambda = 8.12 with Cp = 0.4382. The values below are the ones
% the thesis code shipped with. Both are adjustable; see README.
p.lambda_opt = 8.11;
p.Cp_opt     = 0.435;

%% ---------------- DDP iteration ----------------
p.max_iter       = 50;
p.RMSE_threshold = 0.001;   % early-exit on tracking RMSE
p.alpha          = 0.1;     % control update step (change limiter)
p.alpha_filter   = 0.50;    % first-order filter on U between steps
p.coeff_init     = 0.95;    % first-segment state initialisation scale
p.U_log_num      = 10;      % sentinel that marks "never written this step"

%% ---------------- cost weights ----------------
p.n = 4;   % state dimension [Wr; lambda; Cp; Tg]
p.m = 1;   % control dimension

% PRESERVED QUIRK: Main_Control_Script assigned Q and Qf twice. The first
% pair was immediately overwritten by the second, so the weights actually
% used are the Q1..Q4 set below. The dead assignment is kept as a comment
% so the history is visible rather than lost.
%
%   Q  = diag([5e3, 1e3, 1e3, 1e0]);   % <- overwritten, never used
%   Qf = diag([5e3, 1e3, 1e3, 1e0]);   % <- overwritten, never used
%
p.Q1 = 1*5e9;
p.Q2 = 1*1e7;
p.Q3 = 1*1e7;
p.Q4 = 0.05*5000;
p.Q  = diag([p.Q1, p.Q2, p.Q3, p.Q4]);
p.Qf = diag([p.Q1, p.Q2, p.Q3, p.Q4]);

% PRESERVED FLAW: the backward pass computes the control-cost derivatives as
%       l_u  = R * (U_log(k) - 0);
%       l_uu = R;
% where R is the ROTOR RADIUS (4.4), not Rmat. Rmat below is defined and
% used only in the reported true_cost, so the cost the gains are derived
% from disagrees with the cost that gets logged. Kept verbatim -- project 1
% is a faithful refactor. See README "Known issues".
p.Rmat = 1e-1;

%% ---------------- numerical caps ----------------
p.cap_P  = 1e2;
p.cap_s  = 1e2;
p.cap_s0 = 1e2;
p.reg    = 1e-6;    % Q_uu regularisation

%% ---------------- actuator limits ----------------
p.U_rate  = 15;         % max change per step, relative to previous U
p.U_abs   = [-1e3 1e3]; % applied first
p.U_sat   = [0 1e2];    % applied second, so this is the effective limit

%% ---------------- disturbance ----------------
p.tau_d_amp  = 0.1;   % tau_d = amp * sin(2*pi*k*dt)
p.high_level_sec = 0.01;   % hierarchical Tg_des update period

%% ---------------- selection ----------------
p.selection_criterion = 'RMSE';   % 'RMSE' | 'COST'

%% ---------------- apply overrides ----------------
if nargin > 0 && ~isempty(overrides)
    f = fieldnames(overrides);
    for i = 1:numel(f)
        if ~isfield(p, f{i})
            warning('ddp:params:unknownField', ...
                'Override "%s" is not a known parameter -- adding it anyway.', f{i});
        end
        p.(f{i}) = overrides.(f{i});
    end
    % keep Q consistent if the individual weights were overridden
    if any(ismember({'Q1','Q2','Q3','Q4'}, f)) && ~ismember('Q', f)
        p.Q  = diag([p.Q1, p.Q2, p.Q3, p.Q4]);
        p.Qf = diag([p.Q1, p.Q2, p.Q3, p.Q4]);
    end
end
end
