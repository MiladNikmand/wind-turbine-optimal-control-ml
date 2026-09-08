function out = solve(wind_seg, carryover_state, carryover_U, p, dt_controller, dt_wind, opts)
%DDP.SOLVE  Run DDP-HJB over one wind segment.
%
%   out = ddp.solve(wind_seg, carryover_state, carryover_U, p, ...
%                   dt_controller, dt_wind, opts)
%
% Iterates forward/backward passes up to p.max_iter, stops early when the
% tracking RMSE drops below p.RMSE_threshold, then returns the best
% iteration by p.selection_criterion.
%
% Inputs:
%   wind_seg         1xL predicted wind for this segment, on the controller grid
%   carryover_state  [] on the first segment, else the previous 4x1 final state
%   carryover_U      [] on the first segment, else the previous final control
%   p                ddp.params()
%   dt_controller    controller time step
%   dt_wind          the wind step_size (used only by the iteration-1 guess)
%   opts             .verbose (default true)
%
% Returns a struct with the best-iteration trajectories, the per-iteration
% logs, and the carryover values for the next segment.

    if nargin < 7 || isempty(opts), opts = struct(); end
    if ~isfield(opts, 'verbose'), opts.verbose = true; end

    wind_seg = reshape(wind_seg, 1, []);
    ref = ddp.refs(wind_seg, p);
    L   = ref.L;
    n   = p.n;

    if L < 3
        error('ddp:solve:segmentTooShort', ...
            'Segment has %d samples; need at least 3.', L);
    end

    %% ---------------- preallocation ----------------
    S = struct();
    S.wind_seg = wind_seg;
    S.dt_wind  = dt_wind;
    S.carryover_U = carryover_U;

    S.Pmat = zeros(n, n, L);   S.Pmat(:,:,L) = p.Qf;
    S.s    = zeros(n, L);
    S.s0   = zeros(1, L);
    S.Kmat = zeros(p.m, n, L);
    S.k_ff = zeros(p.m, L);

    S.trace_P  = zeros(1, L);
    S.cond_Quu = zeros(1, L);
    S.TAU_d    = zeros(1, L);

    S.Amat = zeros(n, n, L-1, p.max_iter);
    S.Bmat = zeros(n, L-1, p.max_iter);
    % PRESERVED: Main_Control_Script sized these with an undefined N that
    % leaked in from wind_profile_loader (N = numel(wind)), which silently
    % grew the arrays. Sized to the segment here -- this changes no numbers,
    % it only stops the reallocation. See README "Known issues".
    S.Amat_cl     = zeros(n, n, L-1, p.max_iter);
    S.eig_Amat_cl = zeros(n, L-1, p.max_iter);

    S.X   = zeros(n, L);
    S.Wr  = zeros(1, L);
    S.Tg  = zeros(1, L);
    S.U_filtered = zeros(1, L);

    % CRITICAL: U_log is initialised ONCE, here, not inside the iteration
    % loop. Main_Control_Script does the same. DDP is a line search around
    % the PREVIOUS control trajectory -- iteration i computes
    %     U_k = U_log(k) - alpha*delta_u
    % where U_log(k) holds iteration i-1's value at that step. Resetting it
    % each iteration makes the update start from the sentinel instead of the
    % previous trajectory, and the solver diverges immediately (NaN by
    % iteration 2, then overflow). Only Wr, Tg and U_log(1) are reset per
    % iteration.
    S.U_log = p.U_log_num * ones(1, L);

    X_log     = zeros(n, L, p.max_iter);
    U_log_all = zeros(p.max_iter, L);
    Wr_log    = zeros(p.max_iter, L);
    Tg_log    = zeros(p.max_iter, L);
    lambda_log= zeros(p.max_iter, L);
    Cp_log    = zeros(p.max_iter, L);
    RMSE_log  = zeros(1, p.max_iter);
    cost_total= zeros(1, p.max_iter);
    iter_time = zeros(p.max_iter, 1);

    n_iter_done = 0;

    %% ---------------- iterations ----------------
    for iter = 1:p.max_iter
        iter_start = tic;

        % ---- reset state (exactly what the original resets per iteration) ----
        S.Wr = zeros(1, L);
        S.Tg = zeros(1, L);
        % NOTE: S.U_log is deliberately NOT reset here -- see the comment at
        % its allocation above. Only its first element is set, below.

        if isempty(carryover_state)
            S.Wr(1) = p.coeff_init * ref.Wr_opt(1);
            S.Tg(1) = p.coeff_init * ref.Tg_opt(1);
            S.X(:,1) = [p.coeff_init * ref.Wr_opt(1);
                        p.coeff_init * ref.lambda_opt(1);
                        p.coeff_init * ref.Cp_opt(1);
                        p.coeff_init * ref.Tg_opt(1)];
            S.U_log(1) = p.Rs * ref.Tg_opt(1);
        else
            S.Wr(1) = carryover_state(1);
            S.Tg(1) = carryover_state(4);
            S.X(:,1) = carryover_state(:);
            S.U_log(1) = carryover_U;   % control continuity across segments
        end

        % ---- passes ----
        S = ddp.forward(S, ref, p, dt_controller, iter);
        S = ddp.backward(S, ref, p, dt_controller, iter);

        % ---- logs ----
        X_log(:,:,iter)   = S.X;
        U_log_all(iter,:) = S.U_log;
        Wr_log(iter,:)    = S.Wr;
        Tg_log(iter,:)    = S.Tg;
        % Main_Control_Script wrote scalars into these rows (lambda_val and
        % Cp_k left over from the backward loop at k=1), so the stored rows
        % were constant. The real trajectories live in X_log rows 2 and 3;
        % they are recorded properly here so downstream plots are meaningful.
        lambda_log(iter,:) = S.X(2,:);
        Cp_log(iter,:)     = S.X(3,:);

        rmse = sqrt(mean((S.Wr - ref.Wr_opt).^2));
        RMSE_log(iter)   = rmse;
        cost_total(iter) = S.true_cost;
        iter_time(iter)  = toc(iter_start);
        n_iter_done = iter;

        if opts.verbose
            fprintf('   iter %2d | RMSE = %.4f | cost = %.3e | %.2f s\n', ...
                iter, rmse, S.true_cost, iter_time(iter));
        end

        if rmse < p.RMSE_threshold
            if opts.verbose
                fprintf('   converged early at iteration %d (RMSE = %.4f)\n', iter, rmse);
            end
            break;
        end
    end

    %% ---------------- pick the best iteration ----------------
    valid = 1:n_iter_done;
    switch upper(p.selection_criterion)
        case 'RMSE'
            [~, bi] = min(RMSE_log(valid));
        case 'COST'
            [~, bi] = min(cost_total(valid));
        otherwise
            error('ddp:solve:badCriterion', ...
                'selection_criterion must be ''RMSE'' or ''COST'' (got "%s").', ...
                p.selection_criterion);
    end
    best_iter = valid(bi);

    % ---- forward-fill control values that were never written ----
    U_best = U_log_all(best_iter, 1:L);
    for k = 2:numel(U_best)
        if U_best(k) == p.U_log_num
            U_best(k) = U_best(k-1);
        end
    end
    U_log_all(best_iter, :) = U_best;

    %% ---------------- output ----------------
    out.best_iter   = best_iter;
    out.n_iter      = n_iter_done;
    out.X           = X_log(:,:,best_iter);
    out.Wr          = Wr_log(best_iter,:);
    out.Tg          = Tg_log(best_iter,:);
    out.U           = U_best;
    out.lambda      = lambda_log(best_iter,:);
    out.Cp          = Cp_log(best_iter,:);
    out.RMSE        = RMSE_log(best_iter);
    out.cost        = cost_total(best_iter);

    out.RMSE_log    = RMSE_log(valid);
    out.cost_log    = cost_total(valid);
    out.iter_time   = iter_time(valid);
    out.total_time  = sum(iter_time(valid));

    out.ref         = ref;
    out.wind_seg    = wind_seg;
    out.trace_P     = S.trace_P;
    out.cond_Quu    = S.cond_Quu;
    out.Kmat        = S.Kmat;
    out.k_ff        = S.k_ff;
    out.Amat_cl     = S.Amat_cl(:,:,:,best_iter);
    out.eig_Amat_cl = S.eig_Amat_cl(:,:,best_iter);

    % carryover for the next segment
    out.carryover_state = X_log(:,end,best_iter);
    out.carryover_U     = U_best(end);
end