%% ========================================================================
%  RUN_PREDICTION  --  wind prediction benchmark for one profile
%
%  Runs every model in models_to_run over every sliding-window segment of a
%  wind profile, ranks them, and writes the full asset set report.html reads.
%
%  HOW TO RUN
%      >> cd C:\path\to\proj1
%      >> run_prediction
%
%  Everything is configured in the CONFIG block below -- there are no
%  prompts, so this can be driven from a loop (see run_all.m).
%
%  REUSING AN EARLIER RUN
%      Set reuse_existing = true to skip prediction entirely and rebuild the
%      figures, JSON, CSV and report from an existing
%      wind_profiles/<id>/prediction_suite_results.mat. Useful when you only
%      want to regenerate outputs, or when you already have results from the
%      old flat-layout scripts. The .mat must contain: predictions,
%      models_to_run, best_idx, RMSE, MAE, MAPE, R2, TIME.
% ========================================================================
clc; close all;

%% ======================= CONFIG =======================
if ~exist('profile_id','var'), profile_id = 'wind03'; end        % windlib id, or a number 1..7
if ~exist('project_root','var'), project_root = pwd; end

if ~exist('step_size','var'), step_size = 0.01; end            % resampling grid for the wind profile
if ~exist('history_sec','var'), history_sec = 15.0; end            % training window
if ~exist('predict_sec','var'), predict_sec = 1.5; end             % forecast horizon
if ~exist('stride_sec','var'), stride_sec = 1.5; end             % window step

if ~exist('selection_criterion','var'), selection_criterion = 'RMSE'; end   % 'RMSE' | 'MAE'

if ~exist('models_to_run','var'), models_to_run = {'rbf','rbf_arima','arima','svr','gp','esn','gru','lstm','bilstm','tcn'}; end
% models_to_run = {'rbf','esn','gru'};      % quick subset

if ~exist('use_parallel','var'), use_parallel = false; end           % parfor across segments
                                % NOTE: each worker keeps its own RBF cache,
                                % so W workers means W training passes.

if ~exist('reuse_existing','var'), reuse_existing = false; end         % skip prediction, rebuild outputs from .mat

% ---- export options ----
if ~exist('export_web_assets','var'), export_web_assets = true; end  % downscaled tier for report.html
if ~exist('make_gifs','var'), make_gifs = true; end
if ~exist('export_segment_figures','var'), export_segment_figures = false; end % ~4 PNGs per (model, segment): 7,600 files
                                % for wind01. Off by default; see .gitignore.
if ~exist('export_stats_figures','var'), export_stats_figures = true; end  % the ten cross-model comparison figures
if ~exist('zoom_types','var'), zoom_types = {'actual_vs_pred','residual','parity'}; end
if ~exist('gif_fps','var'), gif_fps = 2; end

%% ======================= SETUP =======================
t_run = tic;
fprintf('=====================================================\n');
fprintf(' WIND PREDICTION  --  %s\n', profile_id);
fprintf('=====================================================\n');

[wind, T_dur, dt, winfo] = windlib.load(profile_id, step_size);
wind = wind(:);
N = numel(wind);
folder_id = winfo.folder;

num_segments = predictors.num_segments(N, step_size, history_sec, predict_sec, stride_sec);
if num_segments == 0
    windlib.feasibility(step_size, history_sec, predict_sec, stride_sec);
    error('run_prediction:noSegments', ...
        ['%s is %g s long and yields no segments at history=%g predict=%g. ' ...
         'Shorten history_sec (see the table above).'], profile_id, T_dur, history_sec, predict_sec);
end

window_size    = round(history_sec / step_size);
segment_length = round(predict_sec / step_size);
segment_stride = round(stride_sec  / step_size);

fprintf(' profile   : %s (%s)\n', winfo.id, winfo.name);
fprintf(' samples   : %d over %g s at dt = %g\n', N, T_dur, step_size);
fprintf(' geometry  : history %d | horizon %d | stride %d samples\n', ...
    window_size, segment_length, segment_stride);
fprintf(' segments  : %d\n', num_segments);

P = report.paths(project_root, folder_id, struct('create', true));

%% ======================= AVAILABILITY =======================
if ~reuse_existing
    avail = predictors.available(models_to_run);
    if ~all(avail)
        dropped = models_to_run(~avail);
        fprintf('\n');
        predictors.available();
        fprintf(2, 'Excluding %d unavailable model(s): %s\n\n', ...
            numel(dropped), strjoin(dropped, ', '));
        models_to_run = models_to_run(avail);
    end
    if isempty(models_to_run)
        error('run_prediction:noModels', 'No requested model can run on this MATLAB install.');
    end
end

%% ======================= PREDICT (or reuse) =======================
if reuse_existing
    fprintf('\nReusing %s\n', P.mat_file);
    if exist(P.mat_file, 'file') ~= 2
        error('run_prediction:noMat', ...
            'reuse_existing is set but %s does not exist.', P.mat_file);
    end
    Sm = load(P.mat_file);
    need = {'predictions','models_to_run','RMSE','MAE','MAPE','R2','TIME'};
    for r = need
        if ~isfield(Sm, r{1})
            error('run_prediction:badMat', '%s is missing "%s".', P.mat_file, r{1});
        end
    end
    predictions   = Sm.predictions;
    models_to_run = Sm.models_to_run(:)';
    RMSE = Sm.RMSE;  MAE = Sm.MAE;  MAPE = Sm.MAPE;  R2 = Sm.R2;  TIME = Sm.TIME;
    if isfield(Sm,'nRMSE'), nRMSE = Sm.nRMSE; else, nRMSE = RMSE / max(max(wind)-min(wind), eps); end
    if isfield(Sm,'sMAPE'), sMAPE = Sm.sMAPE; else, sMAPE = MAPE; end
    num_segments = size(predictions, 2);
    fprintf('  %d models x %d segments loaded\n', numel(models_to_run), num_segments);

else
    num_models = numel(models_to_run);
    predictions = cell(num_models, num_segments);
    RMSE  = nan(num_models, num_segments);
    MAE   = nan(num_models, num_segments);
    MAPE  = nan(num_models, num_segments);
    sMAPE = nan(num_models, num_segments);
    nRMSE = nan(num_models, num_segments);
    R2    = nan(num_models, num_segments);
    TIME  = nan(num_models, num_segments);
    WALL  = nan(num_models, num_segments);   % externally measured, comparable

    wind_range = max(wind) - min(wind);
    reg = predictors.registry();
    models_ok = false(1, num_models);

    if use_parallel
        pool = gcp('nocreate');
        if isempty(pool)
            parpool('Threads', min(4, maxNumCompThreads));
        end
    end

    fprintf('\n');
    for m = 1:num_models
        name = models_to_run{m};
        fn   = reg.(name);
        fprintf('--- %-10s ', name);
        t_model = tic;

        pred_seg = cell(1, num_segments);
        rm = nan(1,num_segments); ma = rm; mp = rm; sm = rm; nr = rm; r2 = rm; tm = rm; wl = rm;

        if use_parallel
            parfor s = 1:num_segments
                [pred_seg{s}, rm(s), ma(s), mp(s), sm(s), nr(s), r2(s), tm(s), wl(s)] = ...
                    local_one(fn, wind, step_size, s, history_sec, predict_sec, stride_sec, ...
                              window_size, segment_length, segment_stride, wind_range);
            end
        else
            for s = 1:num_segments
                [pred_seg{s}, rm(s), ma(s), mp(s), sm(s), nr(s), r2(s), tm(s), wl(s)] = ...
                    local_one(fn, wind, step_size, s, history_sec, predict_sec, stride_sec, ...
                              window_size, segment_length, segment_stride, wind_range);
            end
        end

        predictions(m,:) = pred_seg;
        RMSE(m,:)=rm; MAE(m,:)=ma; MAPE(m,:)=mp; sMAPE(m,:)=sm;
        nRMSE(m,:)=nr; R2(m,:)=r2; TIME(m,:)=tm; WALL(m,:)=wl;

        models_ok(m) = any(~isnan(rm));
        if models_ok(m)
            fprintf('done  %6.1f s total, %.3f s/segment\n', toc(t_model), mean(wl,'omitnan'));
        else
            fprintf('FAILED on every segment (%.1f s)\n', toc(t_model));
        end
    end

    % drop models that produced nothing usable
    if ~all(models_ok)
        fprintf(2, '\nExcluding %s -- no usable output.\n', strjoin(models_to_run(~models_ok), ', '));
        models_to_run = models_to_run(models_ok);
        predictions = predictions(models_ok,:);
        RMSE=RMSE(models_ok,:); MAE=MAE(models_ok,:); MAPE=MAPE(models_ok,:);
        sMAPE=sMAPE(models_ok,:); nRMSE=nRMSE(models_ok,:); R2=R2(models_ok,:);
        TIME=TIME(models_ok,:); WALL=WALL(models_ok,:);
    end
    if isempty(models_to_run)
        error('run_prediction:allFailed', 'Every model failed. Nothing to report.');
    end
end

num_models = numel(models_to_run);
if ~exist('WALL','var'), WALL = TIME; end

%% ======================= AGGREGATE =======================
avg_RMSE  = mean(RMSE, 2,'omitnan');
avg_nRMSE = mean(nRMSE,2,'omitnan');
avg_MAE   = mean(MAE,  2,'omitnan');
avg_MAPE  = mean(MAPE, 2,'omitnan');
avg_sMAPE = mean(sMAPE,2,'omitnan');
avg_R2    = mean(R2,   2,'omitnan');
avg_TIME  = mean(WALL, 2,'omitnan');   % wall clock, comparable across models

% Diebold-Mariano against the best model
[~, best_m] = min(avg_RMSE);
DM = nan(num_models,1);
e_best = RMSE(best_m,:).^2;
for m = 1:num_models
    if m == best_m, DM(m) = 1; continue; end
    d = (RMSE(m,:).^2 - e_best)';
    d = d(~isnan(d));
    if ~isempty(d) && var(d) > 0
        DM(m) = 2*(1 - normcdf(abs(mean(d)/sqrt(var(d)/numel(d)))));
    end
end

summary_table = table(models_to_run(:), avg_RMSE, avg_nRMSE, avg_MAE, avg_MAPE, ...
    avg_sMAPE, avg_R2, DM, avg_TIME, ...
    'VariableNames', {'Model','RMSE','nRMSE','MAE','MAPE','sMAPE_pct','R2','DM_pvalue','Time_sec'});
summary_table = sortrows(summary_table, 'RMSE');

fprintf('\n=================== SUMMARY (by RMSE) ===================\n');
disp(summary_table);

%% ======================= BEST PER SEGMENT =======================
switch upper(selection_criterion)
    case 'RMSE', score = RMSE;
    case 'MAE',  score = MAE;
    otherwise, error('run_prediction:badCriterion', 'Use ''RMSE'' or ''MAE''.');
end
[best_score, best_idx] = min(score, [], 1, 'omitnan');
dead = isnan(best_score);
if any(dead)
    fprintf(2, '%d segment(s) had no usable prediction -- falling back to the overall best model.\n', nnz(dead));
    best_idx(dead) = best_m;
end
best_model_per_segment = models_to_run(best_idx);

predicted_wind_bestof = predictors.stitch(predictions, best_idx, wind, step_size, ...
    history_sec, predict_sec, stride_sec);

%% ======================= REAL-TIME BUDGET =======================
fprintf('\n=================== REAL-TIME BUDGET ===================\n');
fprintf(' horizon = %.2f s per segment\n\n', predict_sec);
fprintf(' %-12s %14s %10s\n', 'model', 'wall s/segment', 'fits?');
for i = 1:height(summary_table)
    ts = summary_table.Time_sec(i);
    fprintf(' %-12s %14.3f %10s\n', summary_table.Model{i}, ts, ternary_str(ts <= predict_sec));
end
fprintf('\n (prediction only -- DDP time is added in run_control)\n');

%% ======================= EXPORT =======================
fprintf('\n=================== EXPORT ===================\n');
cfg = struct('history_sec',history_sec, 'predict_sec',predict_sec, ...
    'stride_sec',stride_sec, 'num_segments',num_segments, ...
    'models_to_run',{models_to_run}, 'selection_criterion',selection_criterion);

RES = struct('predictions',{predictions}, 'summary_table',summary_table, ...
    'best_idx',best_idx, 'best_model_per_segment',{best_model_per_segment}, ...
    'best_score',best_score, 'predicted_wind_bestof',predicted_wind_bestof, ...
    'RMSE',RMSE, 'MAE',MAE, 'MAPE',MAPE, 'R2',R2, 'TIME',WALL, ...
    'nRMSE',nRMSE, 'sMAPE',sMAPE);

report.export_prediction(P, wind, step_size, cfg, RES, ...
    struct('export_web',export_web_assets, 'make_gifs',make_gifs, 'gif_fps',gif_fps));

if export_stats_figures
    fprintf('  - statistics figures\n');
    report.stats_figures(P.stats, summary_table, RMSE, MAE, MAPE, R2, WALL, ...
        models_to_run, best_idx, wind, step_size, predicted_wind_bestof, ...
        history_sec, predict_sec, stride_sec);
end

if export_segment_figures
    fprintf('  - per-segment figures (%d files)...\n', num_models*num_segments*(1+numel(zoom_types)));
    for m = 1:num_models
        for s = 1:num_segments
            report.segment_figures(P.figures, models_to_run{m}, s, wind, step_size, ...
                window_size, segment_length, segment_stride, predictions{m,s}, zoom_types);
        end
    end
end

if make_gifs
    fprintf('  - all-model comparison gifs\n');
    report.animation(fullfile(P.animation,'all_models_full_view.gif'), wind, step_size, ...
        window_size, segment_length, segment_stride, predictions, models_to_run, 'full', best_idx, gif_fps);
    report.animation(fullfile(P.animation,'all_models_zoomed.gif'), wind, step_size, ...
        window_size, segment_length, segment_stride, predictions, models_to_run, 'zoomed', best_idx, gif_fps);
end

report.write_txt(P.report_txt, folder_id, history_sec, predict_sec, stride_sec, ...
    selection_criterion, models_to_run, summary_table, best_model_per_segment, ...
    best_score, num_segments);

save(P.mat_file, 'summary_table','RMSE','MAE','MAPE','sMAPE','nRMSE','R2','TIME','WALL', ...
    'predictions','best_idx','best_model_per_segment','best_score', ...
    'predicted_wind_bestof','models_to_run','history_sec','predict_sec','stride_sec', ...
    'step_size','selection_criterion','num_segments','profile_id');

report.update_manifest(project_root, folder_id, summary_table, num_segments);
report.build_data(project_root);

fprintf('\nDone in %.1f s.\n', toc(t_run));
fprintf('  %s\n', P.root);
fprintf('  open report.html to view (enable GitHub Pages for the online version)\n');

%% ========================================================================
function [yp, rmse, mae, mape, smape, nrmse, r2, tsec, wall] = local_one( ...
    fn, wind, step_size, s, history_sec, predict_sec, stride_sec, ...
    window_size, segment_length, segment_stride, wind_range)
% One (model, segment) call. Wall clock is measured HERE, outside the
% predictor, because the self-reported segment_time_sec is not comparable
% across models: rbf_arima starts its timer after RBF training, and lstm
% hardcodes it to zero.
    t0 = tic;
    try
        [yp, mtr, tsec] = fn(wind, step_size, s, history_sec, predict_sec, stride_sec, false, false);
        wall = toc(t0);
        yp = yp(:);

        i0 = (s-1)*segment_stride + 1;
        idx = (i0+window_size) : (i0+window_size+segment_length-1);
        y = wind(idx);

        denom = (abs(y) + abs(yp))/2;
        smape = mean(abs(yp - y) ./ max(denom,1e-6)) * 100;
        nrmse = mtr.RMSE / max(wind_range, eps);
        rmse = mtr.RMSE; mae = mtr.MAE; mape = mtr.MAPE; r2 = mtr.R2;
    catch
        wall = toc(t0);
        yp = nan(segment_length,1);
        rmse=nan; mae=nan; mape=nan; smape=nan; nrmse=nan; r2=nan; tsec=nan;
    end
end

function s = ternary_str(c)
    if c, s = 'yes'; else, s = 'NO'; end
end
