function S = backward(S, ref, p, dt_controller, iter)
%DDP.BACKWARD  One backward pass over a segment.
%
%   S = ddp.backward(S, ref, p, dt_controller, iter)
%
% Verbatim port of the backward-pass block from Main_Control_Script,
% including the terminal conditions, the pre-pass smoothing of X, the
% linearised A/B, the second-order tensor contractions, the caps on P and
% s, and the true_cost accumulation.
%
% PRESERVED FLAWS (see README "Known issues") -- do not silently repair:
%
%   1. l_u = R*(U_log(k)-0) and l_uu = R, where R is the ROTOR RADIUS
%      (4.4), not the control weight Rmat. The gains are therefore derived
%      from a different cost than the one reported in true_cost.
%
%   2. A32 = A22*dCp_dlambda + d2Cp_dlambda2 * R * Wr / wind_seg
%      uses mrdivide on two ROW VECTORS, producing a least-squares scalar
%      rather than the elementwise Wr(k)/wind_k the commented-out line
%      above it intended.
%
%   3. dCp/dlambda here uses exp(0.4375/(beta^3+1) - 12.5/D), whose Cp
%      peaks near lambda = 11.9, while the forward pass integrates
%      Cp = 0.22*(116/lambda - 5)*exp(-12.5/lambda), peaking at 8.12.
%      The gradient does not match the plant it differentiates.

    L  = ref.L;
    n  = p.n;
    Rs = p.Rs;  Ls = p.Ls;  R = p.R;
    Kt = p.Kt;  Jt = p.Jt;  Rho = p.Rho;

    % ---- terminal conditions ----
    S.Pmat(:,:,L) = p.Qf;

    kL = L - 1;   % the original used the loop-exit value of k here
    X_ref   = [ref.Wr_opt(kL); ref.lambda_opt(kL); ref.Cp_opt(kL); S.Tg_des(kL)];
    X_ref_N = [ref.Wr_opt(L);  ref.lambda_opt(L);  ref.Cp_opt(L);  S.Tg_des(L)];

    S.s(:,L)  = p.Qf * (S.X(:,L) - X_ref);
    S.s0(L)   = (S.X(:,L) - X_ref_N)' * p.Qf * (S.X(:,L) - X_ref_N);

    % ---- smooth the trajectory before linearising ----
    % Window is 2*dt/dt = 2 for Wr and 1 for the rest, i.e. essentially a
    % no-op for lambda/Cp/Tg. Kept exactly as written.
    S.X(1,:) = smoothdata(S.X(1,:), 'movmean', 2*dt_controller/dt_controller);
    S.X(2,:) = smoothdata(S.X(2,:), 'movmean', dt_controller/dt_controller);
    S.X(3,:) = smoothdata(S.X(3,:), 'movmean', dt_controller/dt_controller);
    S.X(4,:) = smoothdata(S.X(4,:), 'movmean', dt_controller/dt_controller);
    S.X(1,:) = max(S.X(1,:), 1e-2);
    S.X(2,:) = max(S.X(2,:), 1e-2);
    S.X(3,:) = max(S.X(3,:), 1e-2);
    S.X(4,:) = max(S.X(4,:), 1e-2);

    true_cost = 0;

    for k = L-1:-1:1
        wind_k     = max(S.wind_seg(k), 1e-3);
        Wr_k       = S.Wr(k);
        lambda_val = S.X(2,k);
        Cp_k       = S.X(3,k);
        Beta       = 0;

        % ---- dCp/dlambda (see PRESERVED FLAW 3) ----
        D = 0.08*Beta + lambda_val;
        D = max(D, 0.01);
        expPart = exp((0.4375 / (Beta^3 + 1)) - (12.5 / D));
        term1 = -25.52 / D^2;
        term2 = (25.52 / D - 0.088 * Beta - 0.1) * (12.5 / D^2);
        dCp_dlambda = expPart * (term1 + term2);
        if isnan(dCp_dlambda)
            error('ddp:backward:nanGradient', 'NaN in dCp/dlambda at k = %d', k);
        end

        % ---- d2Cp/dlambda2 ----
        piece01 = 25 - (156.25 / D);
        piece02 = 0.088*Beta + 1.1 + (0.8932 / (Beta^3 + 1)) - (25.52 / D);
        piece03 = 51.04 - (638 / D);
        piece04 = exp((0.4375 / (Beta^3 + 1)) - (12.5 / D));
        piece05 = D^3;
        d2Cp_dlambda2 = ((piece01 * piece02) + piece03) * piece04 / piece05;
        d2Cp_dlambda2 = min(max(d2Cp_dlambda2, -50), 50);

        % ---- A matrix ----
        A11 = -Kt / Jt;
        A12 = -(1.5 * Rho * pi * R^3 * wind_k^2 * Cp_k) / (Jt * lambda_val^4);
        A13 = (0.5 * Rho * pi * R^3 * wind_k^2) / (Jt * lambda_val^3);
        A14 = -1 / Jt;

        A21 = A11 * R / wind_k;
        A22 = A12 * R / wind_k;
        A23 = A13 * R / wind_k;
        A24 = A14 * R / wind_k;

        A31 = A21 * dCp_dlambda;
        % PRESERVED FLAW 2: mrdivide of two row vectors. The intended form
        % was  (dCp_dlambda*R/wind_k)*A12 + d2Cp_dlambda2*R*Wr_k/wind_k.
        A32 = A22 * dCp_dlambda + d2Cp_dlambda2 * R * S.Wr / S.wind_seg;
        A33 = A23 * dCp_dlambda;
        A34 = A24 * dCp_dlambda;

        A41 = 0; A42 = 0; A43 = 0; A44 = -Rs / Ls;

        A = [A11 A12 A13 A14 ;
             A21 A22 A23 A24 ;
             A31 A32 A33 A34 ;
             A41 A42 A43 A44 ];

        % ---- B matrix ----
        B1 = -1 / (Jt * Ls);
        B2 = -R / (Jt * Ls * wind_k);
        B3 = B2 * dCp_dlambda;
        B4 = 1 / Ls;
        B  = [B1; B2; B3; B4];

        S.Amat(:,:,k,iter) = A;
        S.Bmat(:,k,iter)   = B;

        X_ref   = [ref.Wr_opt(k); ref.lambda_opt(k); ref.Cp_opt(k); S.Tg_des(k)];
        delta_X = S.X(:,k) - X_ref;

        % ---- second-order terms ----
        [f_xx, f_ux, f_uu] = ddp.second_order(S.X(:,k), S.wind_seg(k), S.U_log(k));

         % if getappdata(0,'use_v2')
         %     [f_xx, f_ux, f_uu] = ddp.second_order_v2(S.X(:,k), S.wind_seg(k), S.U_log(k), p);
         % else
         %     [f_xx, f_ux, f_uu] = ddp.second_order(S.X(:,k), S.wind_seg(k), S.U_log(k));
         % end

        % ---- cost derivatives (see PRESERVED FLAW 1: R is the radius) ----
        l_x  = p.Q * delta_X;
        l_u  = R * (S.U_log(k) - 0);
        l_xx = p.Q;
        l_uu = R;

        P_next = S.Pmat(:,:,k+1);
        s_next = S.s(:,k+1);

        Q_x = l_x + A' * s_next;
        Q_u = l_u + B' * s_next;

        Q_xx = l_xx + A'*P_next*A;
        for i = 1:4
            H_i = squeeze(f_xx(i,:,:,:));
            for j = 1:4
                Q_xx = Q_xx + s_next(i) * squeeze(H_i(j,:,:));
            end
        end

        Q_ux = B'*P_next*A;
        for i = 1:4
            Q_ux = Q_ux + s_next(i)*squeeze(f_ux(i,:,:))';
        end

        Q_uu = l_uu + B'*P_next*B;
        for i = 1:4
            Q_uu = Q_uu + s_next(i)*squeeze(f_uu(i,:,:));
        end

        Q_uu_reg = Q_uu + p.reg*eye(size(Q_uu));

        K_k   = -Q_uu_reg \ Q_ux;
        kff_k = +Q_uu_reg \ Q_u;

        S.Kmat(:,:,k) = K_k';
        S.k_ff(:,k)   = kff_k;

        % PRESERVED ORDERING: the cap is applied to Pmat(:,:,k) BEFORE the
        % Riccati update writes it, so it acts on the previous iteration's
        % value (or zeros on the first). Kept as-is.
        if norm(S.Pmat(:,:,k)) > p.cap_P
            S.Pmat(:,:,k) = S.Pmat(:,:,k) * p.cap_P / norm(S.Pmat(:,:,k));
        end

        S.Pmat(:,:,k) = Q_xx - Q_ux' * (Q_uu_reg \ Q_ux);
        S.Pmat(:,:,k) = (S.Pmat(:,:,k) + S.Pmat(:,:,k)') / 2;
        min_eig = min(real(eig(S.Pmat(:,:,k))));
        if min_eig < 1e-6
            S.Pmat(:,:,k) = S.Pmat(:,:,k) + (1e-6 - min_eig) * eye(size(S.Pmat(:,:,k)));
        end

        S.s(:,k) = Q_x - K_k'*Q_uu*kff_k + Q_ux'*kff_k + K_k'*Q_u;
        if norm(S.s(:,k)) > p.cap_s
            S.s(:,k) = S.s(:,k) * p.cap_s / norm(S.s(:,k));
        end

        S.s0(k) = S.s0(k+1) + 0.5*kff_k'*Q_uu*kff_k + kff_k'*Q_u;
        if norm(S.s0(k)) > p.cap_s0
            S.s0(k) = S.s0(k) * p.cap_s0 / norm(S.s0(k));
        end

        S.trace_P(k)  = trace(S.Pmat(:,:,k));
        S.cond_Quu(k) = cond(Q_uu);

        dx = S.X(:,k) - X_ref;
        true_cost = true_cost + (0.5 * dx' * p.Q * dx + 0.5 * p.Rmat * S.U_log(k)^2);

        S.Amat_cl(:,:,k,iter)   = A + B*K_k;
        S.eig_Amat_cl(:,k,iter) = eig(A + B*K_k);
    end

    S.true_cost = true_cost;
end
