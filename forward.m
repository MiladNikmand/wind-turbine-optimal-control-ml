function S = forward(S, ref, p, dt_controller, iter)
%DDP.FORWARD  One forward pass over a segment.
%
%   S = ddp.forward(S, ref, p, dt_controller, iter)
%
% Verbatim port of the forward-pass block from Main_Control_Script's
% segment loop. Rolls the plant forward under the current policy
% (k_ff + K*delta_X on iterations 2+, an open-loop guess on iteration 1),
% applying the hierarchical Tg_des update, the rate limiter, the filter
% and both saturations in the original order.
%
% S carries the mutable per-iteration state:
%   S.X, S.U_log, S.U_filtered, S.Wr, S.Tg, S.Tg_des, S.TAU_d,
%   S.Kmat, S.k_ff, S.carryover_U

    L  = ref.L;
    Rs = p.Rs;
    Ls = p.Ls;

    % hierarchical (slow) loop period, in steps
    high_level_steps = floor(p.high_level_sec / dt_controller);
    if high_level_steps < 1, high_level_steps = 1; end

    S.Tg_des = ref.Tg_opt;   % reset the desired-Tg trajectory each iteration

    for k = 1:L-1
        % ---- high-level control (slow loop) ----
        if mod(k, high_level_steps) == 0
            Wr_error = S.Wr(k) - ref.Wr_opt(k);
            S.Tg_des(k:min(k+high_level_steps-1, L)) = ref.Tg_opt(k) + 0.5*Wr_error;
        end

        wind_k = max(S.wind_seg(k), 1e-3);
        S.TAU_d(k) = p.tau_d_amp * sin(2*pi*k*dt_controller);

        % ---- control ----
        if iter == 1
            if isempty(S.carryover_U)
                % PRESERVED: dt (the wind step size) appears here, not
                % dt_controller. Original line was
                %   U_k = (Rs * Tg_opt(k)) / (1 + dt*Rs/Ls);
                U_k = (Rs * ref.Tg_opt(k)) / (1 + S.dt_wind*Rs/Ls);
            else
                U_k = S.carryover_U;
            end
            S.U_log(k) = U_k;
        else
            X_ref   = [ref.Wr_opt(k); ref.lambda_opt(k); ref.Cp_opt(k); S.Tg_des(k)];
            delta_X = S.X(:,k) - X_ref;
            delta_u = S.k_ff(:,k) + S.Kmat(:,:,k) * delta_X;
            U_k     = S.U_log(k) - p.alpha * delta_u;

            if k == 1
                S.U_filtered(k) = U_k;
            else
                U_k = p.alpha_filter * S.U_filtered(k-1) + (1 - p.alpha_filter) * U_k;
                S.U_filtered(k) = U_k;
            end
        end

        % ---- clamping, in the original order ----
        U_k = max(min(U_k, S.U_log(k) + p.U_rate), S.U_log(k) - p.U_rate);
        U_k = max(min(U_k, p.U_abs(2)), p.U_abs(1));
        U_k = max(min(U_k, p.U_sat(2)), p.U_sat(1));
        S.U_log(k) = U_k;

        % ---- plant step ----
        S.X(:,k+1) = ddp.rk4(S.X(:,k), U_k, S.Tg_des(k), wind_k, k, dt_controller);
        S.Wr(k+1)  = S.X(1,k+1);
        S.Tg(k+1)  = S.X(4,k+1);
    end
end
