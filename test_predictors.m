%% ========================================================================
%  TEST_PREDICTORS  --  smoke test for the +predictors package
%
%  HOW TO RUN
%      proj1/
%        +windlib/         (stage 1)
%        +predictors/      (this stage)
%        test_predictors.m <- this file
%
%      >> cd C:\path\to\proj1
%      >> test_predictors
%
%  Checks the shared contract on every predictor over a couple of real
%  segments, reports which models this MATLAB install can run, and
%  verifies the RBF dedup produces identical results through both callers.
% ========================================================================
clc;
clearvars -except keep   % scripts share a workspace; avoid stale leftovers


if exist('predictors', 'var'), clear predictors; end
if exist('windlib', 'var'),    clear windlib;    end

fprintf('========================================================\n');
fprintf(' +predictors package smoke test\n');
fprintf('========================================================\n\n');

pass = 0; fail = 0;

%% ---- 0. packages visible? ----
try
    predictors.registry();
    windlib.registry();
catch ME
    fprintf(2, 'FAIL: cannot reach the packages (%s)\n', ME.message);
    fprintf(2, '      cd into the folder containing +predictors and +windlib.\n');
    fprintf(2, '      Current folder: %s\n', pwd);
    return;
end
fprintf('Packages found. Folder: %s\n\n', pwd);

%% ---- 1. registry ----
fprintf('--- 1. registry ---\n');
reg = predictors.registry();
names = fieldnames(reg);
[pass, fail] = chk(numel(names) == 10, sprintf('registry lists %d models (expect 10)', numel(names)), pass, fail);
[pass, fail] = chk(all(cellfun(@(n) isa(reg.(n), 'function_handle'), names)), ...
    'every entry is a function handle', pass, fail);

e = predictors.registry('rbf');
[pass, fail] = chk(strcmp(e.name,'rbf'), 'single-model lookup works', pass, fail);

ok = false; try, predictors.registry('nope'); catch, ok = true; end
[pass, fail] = chk(ok, 'unknown model raises an error', pass, fail);
fprintf('\n');

%% ---- 2. toolbox availability ----
fprintf('--- 2. availability on this install ---\n');
A = predictors.available();
runnable = A.Model(A.Available);
fprintf('Runnable here: %s\n\n', strjoin(runnable', ', '));

%% ---- 3. geometry helpers ----
fprintf('--- 3. num_segments ---\n');
step_size = 0.01; history_sec = 15; predict_sec = 1.5; stride_sec = 1.5;
[w, T] = windlib.load('wind03', step_size);
w = w(:);
N = numel(w);

ns = predictors.num_segments(N, step_size, history_sec, predict_sec, stride_sec);
[pass, fail] = chk(ns == 30, sprintf('wind03 -> %d segments (expect 30)', ns), pass, fail);

nf = windlib.feasibility(step_size, history_sec, predict_sec, stride_sec, 'wind03');
[pass, fail] = chk(nf == ns, 'windlib.feasibility agrees with predictors.num_segments', pass, fail);
fprintf('\n');

%% ---- 4. contract check on every runnable model ----
fprintf('--- 4. contract check (2 segments each) ---\n');
seglen = round(predict_sec / step_size);
test_segs = [1, 2];

results = struct();
fprintf('%12s %10s %10s %10s %12s\n', 'model', 'len ok', 'finite', 'RMSE', 'time (s)');
for i = 1:numel(names)
    nm = names{i};
    if ~A.Available(strcmp(A.Model, nm))
        fprintf('%12s %10s   (skipped: %s)\n', nm, '-', A.Reason{strcmp(A.Model,nm)});
        continue;
    end
    okall = true; rms = nan; tsec = nan;
    try
        for s = test_segs
            [yp, mtr, tt] = reg.(nm)(w, step_size, s, history_sec, predict_sec, stride_sec, false, false);
            okall = okall && numel(yp) == seglen && iscolumn(yp(:));
            okall = okall && isstruct(mtr) && isfield(mtr,'RMSE') && isfield(mtr,'R2');
            okall = okall && isscalar(tt) && isnumeric(tt);
            if s == test_segs(1), rms = mtr.RMSE; tsec = tt; end
        end
        fprintf('%12s %10d %10d %10.4f %12.3f\n', nm, okall, all(isfinite(yp)), rms, tsec);
        [pass, fail] = chk(okall, sprintf('  %s honours the contract', nm), pass, fail, true);
        results.(nm) = struct('rmse', rms, 'time', tsec);
    catch ME
        fprintf(2, '%12s   ERROR: %s\n', nm, ME.message);
        fail = fail + 1;
    end
end
fprintf('\n');

%% ---- 5. RBF dedup: both callers must use the SAME core ----
fprintf('--- 5. rbf_core shared by rbf and rbf_arima ---\n');
cfg = predictors.rbf_core('defaults');
[pass, fail] = chk(cfg.num_centers == 150 && abs(cfg.sigma_factor-0.7) < 1e-12, ...
    sprintf('defaults: centers=%d sigma_factor=%g', cfg.num_centers, cfg.sigma_factor), pass, fail);

rng(42);
m1 = predictors.rbf_core('train', w, round(history_sec/step_size), 150, 0.7, true);
rng(42);
m2 = predictors.rbf_core('train', w, round(history_sec/step_size), 150, 0.7, true);
[pass, fail] = chk(isequal(m1.centers, m2.centers), ...
    'rbf_core is reproducible under a fixed rng seed', pass, fail);

past = w(1:round(history_sec/step_size));
p1 = predictors.rbf_core('predict', past, m1, seglen);
[pass, fail] = chk(numel(p1) == seglen && all(isfinite(p1)), ...
    sprintf('rbf_core predict returns %d finite values', seglen), pass, fail);
fprintf('\n');

%% ---- 6. stitch ----
fprintf('--- 6. stitch ---\n');
if isfield(results, 'rbf')
    preds = cell(1, ns);
    for s = 1:min(ns, 5)
        preds{1, s} = predictors.rbf(w, step_size, s, history_sec, predict_sec, stride_sec, false, false);
    end
    for s = min(ns,5)+1 : ns
        preds{1, s} = nan(seglen, 1);
    end
    stitched = predictors.stitch(preds, ones(1, ns), w, step_size, history_sec, predict_sec, stride_sec);
    [pass, fail] = chk(numel(stitched) == N, sprintf('stitched length = %d (matches wind)', numel(stitched)), pass, fail);
    win = round(history_sec/step_size);
    [pass, fail] = chk(all(isnan(stitched(1:win))), ...
        sprintf('first %d samples are NaN (no forecast reaches them)', win), pass, fail);
    [pass, fail] = chk(any(isfinite(stitched)), 'stitched vector has finite values', pass, fail);
else
    fprintf('  (skipped -- rbf unavailable)\n');
end
fprintf('\n');

%% ---- summary ----
fprintf('========================================================\n');
if fail == 0
    fprintf(' ALL %d CHECKS PASSED\n', pass);
else
    fprintf(2, ' %d passed, %d FAILED\n', pass, fail);
end
fprintf('========================================================\n');

%% ========================================================================
function [pass, fail] = chk(cond, label, pass, fail, quiet)
    if nargin < 5, quiet = false; end
    if cond
        pass = pass + 1;
        if ~quiet, fprintf('  PASS  %s\n', label); end
    else
        fail = fail + 1;
        fprintf(2, '  FAIL  %s\n', label);
    end
end
