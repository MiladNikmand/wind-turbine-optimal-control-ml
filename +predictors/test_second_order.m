%% ========================================================================
%  TEST_SECOND_ORDER  --  audit the DDP second-order tensors
%
%      >> cd C:\path\to\proj1
%      >> test_second_order
%
%  Diagnoses the suspected state-ordering mismatch in ddp.second_order.
%  Changes nothing. Run it, read the verdict, then decide.
%
%  The project's state vector is  x = [Wr; lambda; Cp; Tg]  -- see
%  ddp.rk4 lines 31-35, ddp.forward (S.Tg <- X(4,:)), ddp.backward
%  (X_ref = [Wr_opt; lambda_opt; Cp_opt; Tg_des]) and ddp.params (p.n note).
%  ddp.second_order reads x(3) as Tg and x(4) as Cp, which is the reverse.
%
%  TEST 1 is the one that matters and it needs no reference model:
%
%      Tg_dot = (1/Ls)*(-Rs*Tg + U)
%
%  is exactly affine in x and in u. Its second derivatives are therefore
%  identically zero, whatever Cp model you use, whatever the turbine
%  constants are. So row 4 of f_xx and f_ux MUST be all zeros. If it isn't,
%  something nonlinear has been written into the generator-torque row --
%  which is exactly what a 3<->4 swap would do.
%
%  Tests 2-4 add supporting evidence. Test 5 is a finite-difference check
%  and is the only one that depends on choosing a reference model, so its
%  result is reported as advisory rather than as a verdict.
% ========================================================================
clc;

fprintf('=====================================================\n');
fprintf(' SECOND-ORDER TENSOR AUDIT\n');
fprintf('=====================================================\n\n');

%% ---- operating point ----
% A physically sensible point: lambda near optimal, Cp near peak.
p = ddp.params();
wind_test = 10;                        % m/s
lambda0   = p.lambda_opt;              % 8.11
Wr0       = lambda0 * wind_test / p.R; % consistent with ddp.refs
Cp0       = p.Cp_opt;                  % 0.435
Tg0       = 50;                        % mid-range, inside p.U_sat
U0        = 50;

x0 = [Wr0; lambda0; Cp0; Tg0];

fprintf('Operating point (project ordering [Wr; lambda; Cp; Tg]):\n');
fprintf('  Wr     = %8.4f rad/s\n', x0(1));
fprintf('  lambda = %8.4f\n',        x0(2));
fprintf('  Cp     = %8.4f\n',        x0(3));
fprintf('  Tg     = %8.4f\n',        x0(4));
fprintf('  wind   = %8.4f m/s,  U = %.4f\n\n', wind_test, U0);

[f_xx, f_ux, f_uu] = ddp.second_order(x0, wind_test, U0);

fprintf('Returned shapes: f_xx %s, f_ux %s, f_uu %s\n\n', ...
    mat2str(size(f_xx)), mat2str(size(f_ux)), mat2str(size(f_uu)));

tol = 1e-12;
verdict = struct('swap_evidence', 0, 'tests_run', 0);

%% =====================================================================
%  TEST 1  --  the generator-torque row must be identically zero
%  =====================================================================
fprintf('--- 1. Tg_dot row (row 4) must be all zeros ---\n');
fprintf('    Tg_dot = (1/Ls)*(-Rs*Tg + U) is affine, so d2(Tg_dot) = 0.\n\n');

row4_xx = f_xx(4,:,:,:);
row4_ux = f_ux(4,:,:);
n4 = nnz(abs(row4_xx) > tol) + nnz(abs(row4_ux) > tol);
verdict.tests_run = verdict.tests_run + 1;

if n4 == 0
    fprintf('    PASS  row 4 is zero.\n\n');
else
    fprintf(2, '    FAIL  %d nonzero entries in the Tg_dot row.\n', n4);
    verdict.swap_evidence = verdict.swap_evidence + 1;
    report_nonzeros(f_xx, f_ux, 4, tol);
end

%% =====================================================================
%  TEST 2  --  the Cp row should NOT be zero
%  =====================================================================
fprintf('--- 2. Cp_dot row (row 3) should be nonzero ---\n');
fprintf('    Cp_dot = (dCp/dlambda)*lambda_dot carries d2Cp/dlambda2,\n');
fprintf('    so it is the one genuinely nonlinear row.\n\n');

row3_xx = f_xx(3,:,:,:);
row3_ux = f_ux(3,:,:);
n3 = nnz(abs(row3_xx) > tol) + nnz(abs(row3_ux) > tol);
verdict.tests_run = verdict.tests_run + 1;

if n3 == 0
    fprintf(2, '    FAIL  the Cp_dot row is entirely zero.\n');
    fprintf(2, '          The curvature that DDP exists to exploit is missing\n');
    fprintf(2, '          from this row -- consistent with it living in row 4.\n\n');
    verdict.swap_evidence = verdict.swap_evidence + 1;
else
    fprintf('    PASS  %d nonzero entries.\n\n', n3);
    report_nonzeros(f_xx, f_ux, 3, tol);
end

%% =====================================================================
%  TEST 3  --  does swapping rows 3 and 4 restore the expected pattern?
%  =====================================================================
fprintf('--- 3. Row-swap probe ---\n');
fprintf('    If rows 3 and 4 are transposed, swapping them should make\n');
fprintf('    row 4 vanish and row 3 become nonlinear.\n\n');

g_xx = f_xx;  g_ux = f_ux;
g_xx([3 4],:,:,:) = f_xx([4 3],:,:,:);
g_ux([3 4],:,:)   = f_ux([4 3],:,:);

s4 = nnz(abs(g_xx(4,:,:,:)) > tol) + nnz(abs(g_ux(4,:,:)) > tol);
s3 = nnz(abs(g_xx(3,:,:,:)) > tol) + nnz(abs(g_ux(3,:,:)) > tol);
verdict.tests_run = verdict.tests_run + 1;

fprintf('    as shipped : row3 = %d nonzeros, row4 = %d nonzeros\n', n3, n4);
fprintf('    if swapped : row3 = %d nonzeros, row4 = %d nonzeros\n\n', s3, s4);

if s4 == 0 && s3 > 0 && ~(n4 == 0 && n3 > 0)
    fprintf(2, '    Swapping rows 3 and 4 produces the mathematically\n');
    fprintf(2, '    expected pattern; the shipped version does not.\n\n');
    verdict.swap_evidence = verdict.swap_evidence + 1;
else
    fprintf('    Swapping does not improve the pattern.\n\n');
end

%% =====================================================================
%  TEST 4  --  hardcoded constants vs ddp.params
%  =====================================================================
fprintf('--- 4. Constant drift between second_order.m and ddp.params ---\n\n');
verdict.tests_run = verdict.tests_run + 1;

src = fileread(fullfile('+ddp','second_order.m'));
names = {'R','Kt','Jt','Rs','Rho','Ls'};
drift = 0;
for i = 1:numel(names)
    tok = regexp(src, ['^\s*' names{i} '\s*=\s*([0-9.eE+-]+)\s*;'], ...
                 'tokens','lineanchors','once');
    if isempty(tok)
        fprintf('    %-4s not found as a literal\n', names{i});
        continue;
    end
    lit = str2double(tok{1});
    ref = p.(names{i});
    if abs(lit - ref) > 1e-12
        fprintf(2, '    %-4s literal %-10g  ddp.params %-10g   DRIFTED\n', ...
            names{i}, lit, ref);
        drift = drift + 1;
    else
        fprintf('    %-4s literal %-10g  matches ddp.params\n', names{i}, lit);
    end
end
if drift == 0
    fprintf('\n    In sync today, but the values are duplicated rather than\n');
    fprintf('    shared -- editing ddp.params alone will silently desync them.\n\n');
else
    fprintf(2, '\n    %d constant(s) already out of sync.\n\n', drift);
end

%% =====================================================================
%  TEST 5  --  finite-difference cross-check   (ADVISORY)
%  =====================================================================
fprintf('--- 5. Finite-difference cross-check (advisory) ---\n');
fprintf('    Compares against a reference model built from the RK4 stage-1\n');
fprintf('    equations. ddp.backward''s analytic A matrix assumes a different\n');
fprintf('    aerodynamic torque than ddp.rk4 evaluates, so disagreement here\n');
fprintf('    may reflect that known split rather than the ordering. Read the\n');
fprintf('    per-row pattern, not the absolute numbers.\n\n');

h = 1e-5;
fd = zeros(4,4,4);                       % fd(i,j,k) = d2 f_i / dx_j dx_k
for j = 1:4
    for k = 1:4
        ep = zeros(4,1); eq = zeros(4,1);
        ep(j) = h; eq(k) = h;
        fpp = ref_dynamics(x0+ep+eq, U0, wind_test, p);
        fpm = ref_dynamics(x0+ep-eq, U0, wind_test, p);
        fmp = ref_dynamics(x0-ep+eq, U0, wind_test, p);
        fmm = ref_dynamics(x0-ep-eq, U0, wind_test, p);
        fd(:,j,k) = (fpp - fpm - fmp + fmm) / (4*h*h);
    end
end

% collapse the shipped 4x4x4x4 the way ddp.backward consumes it: sum over j
shipped = squeeze(sum(f_xx, 2));         % 4x4x4

fprintf('    row   max|FD|      max|shipped|   note\n');
labels = {'Wr_dot','lambda_dot','Cp_dot','Tg_dot'};
for i = 1:4
    a = max(abs(reshape(fd(i,:,:),1,[])));
    b = max(abs(reshape(shipped(i,:,:),1,[])));
    note = '';
    if a < 1e-6 && b > 1e-6, note = '<- FD says linear, shipped says curved'; end
    if a > 1e-6 && b < 1e-6, note = '<- FD says curved, shipped says linear'; end
    fprintf('    %-4d %-12.4g %-14.4g %s   (%s)\n', i, a, b, note, labels{i});
end
fprintf('\n');

%% =====================================================================
fprintf('=====================================================\n');
fprintf(' VERDICT\n');
fprintf('=====================================================\n');
if verdict.swap_evidence >= 2
    fprintf(2, ' %d of 3 structural tests point to rows 3 and 4 being\n', verdict.swap_evidence);
    fprintf(2, ' transposed relative to [Wr; lambda; Cp; Tg].\n\n');
    fprintf(' Next: swap line 122 of +ddp/backward.m to call\n');
    fprintf('   ddp.second_order_v2(S.X(:,k), S.wind_seg(k), S.U_log(k), p)\n');
    fprintf(' and re-run one profile to see whether it changes the result.\n');
elseif verdict.swap_evidence == 1
    fprintf(' Mixed signals -- 1 of 3 structural tests flagged. Read the\n');
    fprintf(' detail above before changing anything.\n');
else
    fprintf(' No ordering problem detected. The tensors are consistent\n');
    fprintf(' with x = [Wr; lambda; Cp; Tg].\n');
end
fprintf('\n');

%% ========================================================================
function report_nonzeros(f_xx, f_ux, row, tol)
%REPORT_NONZEROS  Print every nonzero entry in one row, with its meaning.
    fprintf('          entry                value        reads as\n');
    nm = {'Wr','lambda','Cp','Tg'};
    for j = 1:4
        for k = 1:4
            for l = 1:4
                v = f_xx(row,j,k,l);
                if abs(v) > tol
                    fprintf('          f_xx(%d,%d,%d,%d) %14.6g   d2 f_%s / d%s d%s\n', ...
                        row, j, k, l, v, nm{row}, nm{k}, nm{l});
                end
            end
        end
    end
    for j = 1:4
        v = f_ux(row,1,j);
        if abs(v) > tol
            fprintf('          f_ux(%d,1,%d)   %14.6g   d2 f_%s / du d%s\n', ...
                row, j, v, nm{row}, nm{j});
        end
    end
    fprintf('\n');
end

% ------------------------------------------------------------------------
function dx = ref_dynamics(x, U, wind, p)
%REF_DYNAMICS  Continuous-time reference, x = [Wr; lambda; Cp; Tg].
%   Mirrors the RK4 stage-1 equations in ddp.rk4, with Cp treated as a
%   free state (which is what a 4-state Jacobian implies) rather than
%   recomputed from lambda.
    Wr = x(1);  lambda = max(real(x(2)), 0.5);  Cp = x(3);  Tg = x(4);

    Ta = 0.5 * p.Rho * pi * p.R^3 * wind^2 * (Cp / lambda);

    Wr_dot     = (1/p.Jt) * (Ta - p.Kt*Wr + (p.Rs/p.Ls)*Tg - (1/p.Ls)*U);
    lambda_dot = (p.R / wind) * Wr_dot;

    Beta = 0;
    D    = max(0.08*Beta + lambda, 1e-2);
    expP = exp((0.4375/(Beta^3+1)) - (12.5/D));
    dCp_dlambda = expP * (-25.52/D^2 + (25.52/D - 0.088*Beta - 0.1)*(12.5/D^2));

    Cp_dot = dCp_dlambda * lambda_dot;
    Tg_dot = (1/p.Ls) * (-p.Rs*Tg + U);

    dx = [Wr_dot; lambda_dot; Cp_dot; Tg_dot];
end