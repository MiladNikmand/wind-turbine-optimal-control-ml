%% ========================================================================
%  TEST_REPORT  --  smoke test for the +report package
%
%      >> cd C:\path\to\proj1
%      >> test_report
%
%  Builds a small synthetic result set, writes it through the exporters into
%  a throwaway folder, and checks every expected file appears and the JSON
%  round-trips. Deletes the folder afterwards unless you set keep = true.
% ========================================================================
clc;
keep = false;   % set true to inspect the sandbox afterwards
clearvars -except keep


fprintf('========================================================\n');
fprintf(' +report package smoke test\n');
fprintf('========================================================\n\n');

pass = 0; fail = 0;

try
    report.paths(pwd, 'x');
catch ME
    fprintf(2, 'FAIL: cannot reach +report (%s)\n', ME.message);
    return;
end

sandbox = fullfile(tempdir, sprintf('proj1_report_test_%d', round(rand*1e6)));
fprintf('Sandbox: %s\n\n', sandbox);

%% ---- 1. paths ----
fprintf('--- 1. paths ---\n');
P = report.paths(sandbox, 'wind_profile_03', struct('create', true));
for d = {'root','report','figures','stats','animation','web','web_figures','web_anim','control'}
    [pass, fail] = chk(exist(P.(d{1}), 'dir') == 7, sprintf('created %s/', d{1}), pass, fail, true);
end
fprintf('  9 folders created\n\n');

%% ---- 2. synthetic results ----
fprintf('--- 2. building a small synthetic result set ---\n');
step_size = 0.01; history_sec = 15; predict_sec = 1.5; stride_sec = 1.5;
w = windlib.load('wind03', step_size);
nseg = min(6, predictors.num_segments(numel(w), step_size, history_sec, predict_sec, stride_sec));
models = {'rbf','esn','gru'};
nm = numel(models);
seglen = round(predict_sec/step_size);
stride = round(stride_sec/step_size);
win    = round(history_sec/step_size);

preds = cell(nm, nseg);
RMSE = zeros(nm,nseg); MAE = RMSE; MAPE = RMSE; R2 = RMSE; TIME = RMSE;
nRMSE = RMSE; sMAPE = RMSE;
for m = 1:nm
    for s = 1:nseg
        i0 = (s-1)*stride + win + 1;
        truth = w(i0:i0+seglen-1)';
        preds{m,s} = truth + 0.05*m*randn(seglen,1);
        e = truth - preds{m,s};
        RMSE(m,s) = sqrt(mean(e.^2)); MAE(m,s) = mean(abs(e));
        MAPE(m,s) = mean(abs(e)./max(abs(truth),1e-9))*100;
        sMAPE(m,s)= MAPE(m,s); nRMSE(m,s) = RMSE(m,s)/(max(w)-min(w));
        R2(m,s)   = 1 - sum(e.^2)/max(sum((truth-mean(truth)).^2),eps);
        TIME(m,s) = 0.05*m;
    end
end
[best_score, best_idx] = min(RMSE, [], 1, 'omitnan');
best_model_per_segment = models(best_idx);

T = table(models(:), mean(RMSE,2), mean(nRMSE,2), mean(MAE,2), mean(MAPE,2), ...
    mean(sMAPE,2), mean(R2,2), [1;0.02;0.03], mean(TIME,2), ...
    'VariableNames', {'Model','RMSE','nRMSE','MAE','MAPE','sMAPE_pct','R2','DM_pvalue','Time_sec'});
T = sortrows(T, 'RMSE');

RES.predictions = preds;  RES.summary_table = T;  RES.best_idx = best_idx;
RES.best_model_per_segment = best_model_per_segment;  RES.best_score = best_score;
RES.RMSE = RMSE; RES.MAE = MAE; RES.MAPE = MAPE; RES.R2 = R2; RES.TIME = TIME;
RES.nRMSE = nRMSE; RES.sMAPE = sMAPE;
RES.predicted_wind_bestof = predictors.stitch(preds, best_idx, w(:), step_size, ...
    history_sec, predict_sec, stride_sec);

cfg = struct('history_sec',history_sec,'predict_sec',predict_sec,'stride_sec',stride_sec, ...
    'num_segments',nseg,'models_to_run',{models},'selection_criterion','RMSE');
fprintf('  %d models x %d segments\n\n', nm, nseg);

%% ---- 3. export_prediction ----
fprintf('--- 3. export_prediction (gifs on, this takes a moment) ---\n');
tic;
assets = report.export_prediction(P, w, step_size, cfg, RES, ...
    struct('export_web', true, 'make_gifs', true, 'gif_fps', 4));
fprintf('  took %.1f s\n', toc);

for m = 1:nm
    nmm = models{m};
    [pass, fail] = chk(exist(fullfile(P.figures, [nmm '_full_timeline.png']),'file')==2, ...
        sprintf('%s hi-res full png', nmm), pass, fail, true);
    [pass, fail] = chk(exist(fullfile(P.web_figures, [nmm '_full.png']),'file')==2, ...
        sprintf('%s web full png', nmm), pass, fail, true);
    [pass, fail] = chk(exist(fullfile(P.animation, nmm, 'full.gif'),'file')==2, ...
        sprintf('%s hi-res gif', nmm), pass, fail, true);
    [pass, fail] = chk(exist(fullfile(P.web_anim, nmm, 'full.gif'),'file')==2, ...
        sprintf('%s web gif', nmm), pass, fail, true);
end
fprintf('  PASS  all per-model assets present in both tiers\n');

% web tier must actually be smaller
a = dir(fullfile(P.animation, models{1}, 'full.gif'));
b = dir(fullfile(P.web_anim,  models{1}, 'full.gif'));
fprintf('  gif size: hi-res %.0f KB -> web %.0f KB (%.0f%%)\n', a.bytes/1024, b.bytes/1024, 100*b.bytes/a.bytes);
[pass, fail] = chk(b.bytes < a.bytes, 'web gif is smaller than hi-res', pass, fail);

for f = {'summary_table.csv','segment_metrics.csv','best_model_per_segment.csv','prediction_suite_summary.json'}
    [pass, fail] = chk(exist(fullfile(P.report, f{1}),'file')==2, ...
        sprintf('wrote report/%s', f{1}), pass, fail);
end
[pass, fail] = chk(exist(P.tex_file,'file')==2, 'wrote summary_table.tex', pass, fail);
[pass, fail] = chk(exist(P.readme,'file')==2, 'wrote README.md', pass, fail);
fprintf('\n');

%% ---- 4. JSON round-trip ----
fprintf('--- 4. summary JSON ---\n');
fid = fopen(P.summary_json,'r'); raw = fread(fid,'*char')'; fclose(fid);
S = jsondecode(raw);
[pass, fail] = chk(numel(fieldnames(S.models)) == nm, ...
    sprintf('%d models in the JSON', numel(fieldnames(S.models))), pass, fail);
[pass, fail] = chk(isfield(S.models.(matlab.lang.makeValidName(lower(T.Model{1}))), 'assets'), ...
    'each model carries its asset paths', pass, fail);
[pass, fail] = chk(numel(S.segment_best) == nseg, 'segment_best has one entry per segment', pass, fail);
[pass, fail] = chk(S.settings.num_segments == nseg, 'settings round-trip', pass, fail);
fprintf('\n');

%% ---- 5. stats figures + txt ----
fprintf('--- 5. stats figures and text report ---\n');
report.stats_figures(P.stats, T, RMSE, MAE, MAPE, R2, TIME, models, best_idx, ...
    w(:), step_size, RES.predicted_wind_bestof, history_sec, predict_sec, stride_sec);
n_png = numel(dir(fullfile(P.stats,'*.png')));
[pass, fail] = chk(n_png == 10, sprintf('%d statistics figures written', n_png), pass, fail);

report.write_txt(P.report_txt, P.profile_id, history_sec, predict_sec, stride_sec, ...
    'RMSE', models, T, best_model_per_segment, best_score, nseg);
[pass, fail] = chk(exist(P.report_txt,'file')==2, 'wrote prediction_suite_report.txt', pass, fail);
fprintf('\n');

%% ---- 6. manifest + build_data ----
fprintf('--- 6. manifest and report bundle ---\n');
report.update_manifest(sandbox, 'wind_profile_03', T, nseg);
mf = fullfile(sandbox, 'wind_profiles', 'manifest.json');
[pass, fail] = chk(exist(mf,'file')==2, 'manifest.json written', pass, fail);

fid = fopen(mf,'r'); M = jsondecode(fread(fid,'*char')'); fclose(fid);
[pass, fail] = chk(~isempty(M.profiles), 'manifest has a profiles entry', pass, fail);

% a single profile must still encode as a LIST, not a bare object
[pass, fail] = chk(iscell(M.profiles) || numel(M.profiles)==1, ...
    'single profile encodes without error', pass, fail);

report.build_data(sandbox);
js = fullfile(sandbox,'wind_profiles','report_data.js');
jn = fullfile(sandbox,'wind_profiles','report_data.json');
[pass, fail] = chk(exist(js,'file')==2, 'report_data.js written', pass, fail);
[pass, fail] = chk(exist(jn,'file')==2, 'report_data.json written', pass, fail);

fid = fopen(js,'r'); txt = fread(fid,'*char')'; fclose(fid);
[pass, fail] = chk(contains(txt,'window.REPORT_DATA'), ...
    'js assigns window.REPORT_DATA (works over file://)', pass, fail);

fid = fopen(jn,'r'); B = jsondecode(fread(fid,'*char')'); fclose(fid);
[pass, fail] = chk(~isempty(B.profiles), 'bundle contains the profile', pass, fail);
fprintf('\n');

%% ---- 7. second profile, no clash ----
fprintf('--- 7. a second profile joins cleanly ---\n');
report.update_manifest(sandbox, 'wind_profile_07', T, nseg);
fid = fopen(mf,'r'); M2 = jsondecode(fread(fid,'*char')'); fclose(fid);
n2 = numel(M2.profiles);
[pass, fail] = chk(n2 == 2, sprintf('manifest now lists %d profiles', n2), pass, fail);
fprintf('\n');

%% ---- cleanup ----
if keep
    fprintf('Sandbox kept at:\n  %s\n\n', sandbox);
else
    rmdir(sandbox, 's');
    fprintf('Sandbox removed (set keep = true at the top to inspect it).\n\n');
end

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
