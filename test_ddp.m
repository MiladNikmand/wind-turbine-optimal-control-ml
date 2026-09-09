%% ========================================================================
%  TEST_DDP  --  smoke test for the +ddp package
%
%      proj1/
%        +windlib/  +predictors/  +ddp/
%        test_ddp.m
%
%      >> cd C:\path\to\proj1
%      >> test_ddp
%
%  Runs the solver on real predicted wind and checks the pieces behave.
%  Uses a small max_iter so it finishes quickly.
% ========================================================================
clc;
clearvars -except keep   % scripts share a workspace; avoid stale leftovers

fprintf('========================================================\n');
fprintf(' +ddp package smoke test\n');
fprintf('========================================================\n\n');

pass = 0; fail = 0;

try
    ddp.params(); windlib.registry(); predictors.registry();
catch ME
    fprintf(2, 'FAIL: cannot reach the packages (%s)\n', ME.message);
    fprintf(2, '      Folder: %s\n', pwd);
    return;
end

step_size = 0.01; history_sec = 15; predict_sec = 1.5; stride_sec = 1.5;
dt_controller = 0.01;

%% ---- 1. params ----
fprintf('--- 1. params ---\n');
p = ddp.params();
[pass, fail] = chk(abs(p.lambda_opt - 8.11) < 1e-12, sprintf('lambda_opt = %g', p.lambda_opt), pass, fail);
[pass, fail] = chk(abs(p.Cp_opt - 0.435) < 1e-12, sprintf('Cp_opt = %g', p.Cp_opt), pass, fail);
[pass, fail] = chk(isequal(diag(p.Q)', [5e9 1e7 1e7 250]), ...
    sprintf('Q = diag([%g %g %g %g]) -- the second assignment, not the dead one', ...
    p.Q(1,1), p.Q(2,2), p.Q(3,3), p.Q(4,4)), pass, fail);

p2 = ddp.params(struct('max_iter', 5, 'Cp_opt', 0.42));
[pass, fail] = chk(p2.max_iter == 5 && abs(p2.Cp_opt - 0.42) < 1e-12, ...
    'overrides applied', pass, fail);

p3 = ddp.params(struct('Q1', 1e3));
[pass, fail] = chk(p3.Q(1,1) == 1e3, 'Q rebuilt when Q1 overridden', pass, fail);
fprintf('\n');

%% ---- 2. refs ----
fprintf('--- 2. refs ---\n');
wind_seg = 10 + 0.5*sin(linspace(0, 2*pi, 150));
ref = ddp.refs(wind_seg, p);
[pass, fail] = chk(numel(ref.Wr_opt) == 150, 'Wr_opt has one entry per step', pass, fail);
[pass, fail] = chk(all(ref.Wr_opt > 0), 'Wr_opt positive throughout', pass, fail);

expect_Wr = p.lambda_opt * wind_seg(1) / p.R;
[pass, fail] = chk(abs(ref.Wr_opt(1) - expect_Wr) < 1e-12, ...
    sprintf('Wr_opt(1) = lambda*v/R = %.4f', ref.Wr_opt(1)), pass, fail);

K = 0.5*p.Rho*pi*p.R^2*(p.R^3*p.Cp_opt/p.lambda_opt^3);
[pass, fail] = chk(abs(ref.Tg_opt(1) - K*ref.Wr_opt(1)^2) < 1e-9, ...
    sprintf('Tg_opt = K_opt*Wr^2 = %.2f', ref.Tg_opt(1)), pass, fail);
fprintf('\n');

%% ---- 3. rk4 ----
fprintf('--- 3. rk4 single step ---\n');
x0 = [ref.Wr_opt(1); p.lambda_opt; p.Cp_opt; ref.Tg_opt(1)];
x1 = ddp.rk4(x0, 50, ref.Tg_opt(1), wind_seg(1), 1, dt_controller);
[pass, fail] = chk(numel(x1) == 4 && all(isfinite(x1)), ...
    sprintf('x_next = [%.3f %.3f %.4f %.1f]', x1(1), x1(2), x1(3), x1(4)), pass, fail);
[pass, fail] = chk(abs(x1(2) - p.R*x1(1)/wind_seg(1)) < 1e-9, ...
    'lambda is recomputed consistently as R*Wr/v', pass, fail);
fprintf('\n');

%% ---- 4. second_order ----
fprintf('--- 4. second_order tensors ---\n');
[f_xx, f_ux, f_uu] = ddp.second_order(x0, wind_seg(1), 50);
[pass, fail] = chk(isequal(size(f_xx), [4 4 4 4]), sprintf('f_xx is %s', mat2str(size(f_xx))), pass, fail);
[pass, fail] = chk(size(f_ux,1) == 4, sprintf('f_ux is %s', mat2str(size(f_ux))), pass, fail);
[pass, fail] = chk(size(f_uu,1) == 4, sprintf('f_uu is %s', mat2str(size(f_uu))), pass, fail);
fprintf('\n');

%% ---- 5. resample_segment (the one kept fix) ----
fprintf('--- 5. resample_segment ---\n');
v = linspace(8, 9, 150)';
for dtc = [0.01 0.005 0.002 0.001]
    y = ddp.resample_segment(v, step_size, dtc);
    need = round(150 * step_size / dtc);
    [pass, fail] = chk(numel(y) == need, ...
        sprintf('dt=%-6g -> %5d samples (stride-correct, need %d)', dtc, numel(y), need), pass, fail);
end
y = ddp.resample_segment(v, 0.01, 0.01);
[pass, fail] = chk(numel(y) == 150 && max(abs(y - v)) < 1e-12, ...
    'identity when dt_controller == step_size', pass, fail);
fprintf('\n');

%% ---- 6. solve on real predicted wind ----
fprintf('--- 6. solve on a real predicted segment ---\n');
w = windlib.load('wind03', step_size); w = w(:);
[future_pred, ~, ~] = predictors.rbf(w, step_size, 5, history_sec, predict_sec, stride_sec, false, false);
seg = ddp.resample_segment(future_pred, step_size, dt_controller)';
fprintf('    segment: %d samples, wind %.2f .. %.2f m/s\n', numel(seg), min(seg), max(seg));

pq = ddp.params(struct('max_iter', 8));
out = ddp.solve(seg, [], [], pq, dt_controller, step_size, struct('verbose', true));

[pass, fail] = chk(numel(out.Wr) == numel(seg), 'Wr spans the segment', pass, fail);
[pass, fail] = chk(all(isfinite(out.Wr)), 'Wr is finite throughout', pass, fail);
[pass, fail] = chk(all(isfinite(out.U)), 'U is finite throughout', pass, fail);

% --- divergence checks: the BEST iteration being finite is not enough, the
% --- whole iteration sequence has to stay sane. A reset bug in the control
% --- update showed up here as NaN at iteration 2 and 1e63 by iteration 3.
[pass, fail] = chk(all(isfinite(out.RMSE_log)), ...
    sprintf('every iteration RMSE is finite  [%s]', ...
    strjoin(arrayfun(@(v) sprintf('%.3g', v), out.RMSE_log, 'UniformOutput', false), ' ')), pass, fail);
[pass, fail] = chk(all(isfinite(out.cost_log)), 'every iteration cost is finite', pass, fail);
[pass, fail] = chk(max(out.RMSE_log) < 1e4, ...
    sprintf('no iteration blows up (max RMSE = %.4g)', max(out.RMSE_log)), pass, fail);
[pass, fail] = chk(out.RMSE_log(end) <= out.RMSE_log(1) * 10, ...
    sprintf('RMSE does not run away: %.4f -> %.4f', out.RMSE_log(1), out.RMSE_log(end)), pass, fail);

[pass, fail] = chk(~any(out.U == pq.U_log_num), 'no unwritten sentinel values left in U', pass, fail);
[pass, fail] = chk(all(out.U >= pq.U_sat(1) - 1e-9 & out.U <= pq.U_sat(2) + 1e-9), ...
    sprintf('U respects saturation [%g %g]', pq.U_sat(1), pq.U_sat(2)), pass, fail);
[pass, fail] = chk(out.best_iter >= 1 && out.best_iter <= out.n_iter, ...
    sprintf('best_iter = %d of %d', out.best_iter, out.n_iter), pass, fail);
[pass, fail] = chk(out.RMSE == min(out.RMSE_log), ...
    sprintf('best iteration really is the min-RMSE one (%.4f)', out.RMSE), pass, fail);
[pass, fail] = chk(numel(out.carryover_state) == 4, 'carryover_state is 4x1', pass, fail);
[pass, fail] = chk(isscalar(out.carryover_U), 'carryover_U is scalar', pass, fail);
[pass, fail] = chk(~all(out.lambda == out.lambda(1)), ...
    'lambda log is a real trajectory, not a constant row', pass, fail);
fprintf('\n');

%% ---- 7. carryover continuity ----
fprintf('--- 7. two segments in sequence ---\n');
[fp2, ~, ~] = predictors.rbf(w, step_size, 6, history_sec, predict_sec, stride_sec, false, false);
seg2 = ddp.resample_segment(fp2, step_size, dt_controller)';
out2 = ddp.solve(seg2, out.carryover_state, out.carryover_U, pq, dt_controller, step_size, ...
    struct('verbose', false));

[pass, fail] = chk(abs(out2.X(1,1) - out.carryover_state(1)) < 1e-12, ...
    'segment 2 starts exactly where segment 1 ended', pass, fail);
jump = abs(out2.Wr(1) - out.Wr(end));
fprintf('    seam jump in Wr: %.6f rad/s\n', jump);
[pass, fail] = chk(isfinite(jump), 'seam jump is finite', pass, fail);
fprintf('\n');

%% ---- 8. timing vs the real-time budget ----
fprintf('--- 8. real-time budget ---\n');
fprintf('    DDP time for %d iterations: %.2f s\n', out.n_iter, out.total_time);
fprintf('    per iteration             : %.3f s\n', mean(out.iter_time));
fprintf('    horizon (predict_sec)     : %.2f s\n', predict_sec);
full_p = ddp.params();
projected = mean(out.iter_time) * full_p.max_iter;
fprintf('    projected at max_iter=%d   : %.2f s\n', full_p.max_iter, projected);
if projected > predict_sec
    fprintf(2, '    => %.1fx over the horizon budget (prediction not yet included)\n', projected/predict_sec);
else
    fprintf('    => fits inside the horizon budget\n');
end
fprintf('\n');

fprintf('========================================================\n');
if fail == 0
    fprintf(' ALL %d CHECKS PASSED\n', pass);
else
    fprintf(2, ' %d passed, %d FAILED\n', pass, fail);
end
fprintf('========================================================\n');

%% ========================================================================
function [pass, fail] = chk(cond, label, pass, fail)
    if cond
        pass = pass + 1;
        fprintf('  PASS  %s\n', label);
    else
        fail = fail + 1;
        fprintf(2, '  FAIL  %s\n', label);
    end
end
