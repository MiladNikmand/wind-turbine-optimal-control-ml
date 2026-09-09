%% ========================================================================
%  RUN_CONTROL  --  DDP-HJB control driven by cached wind forecasts
%
%  Reads predictions{m,seg} from prediction_suite_results.mat, solves the
%  DDP-HJB problem over each predicted segment, carries state across segment
%  boundaries, and exports everything report.html needs.
%
%  HOW TO RUN
%      >> cd C:\path\to\proj1
%      >> run_prediction          % once per profile, first
%      >> run_control
%
%  COMPARING FORECAST SOURCES
%      Run it three times, changing forecast_mode each time:
%          'truth'   perfect foresight -- the upper bound
%          'bestof'  the per-segment winner
%          'model'   one named predictor, e.g. 'rbf'
%      Each becomes a selectable run in section 5 of the report, so the cost
%      of prediction error is the gap between them.
% ========================================================================
clc; close all;

%% ======================= CONFIG =======================
if ~exist('profile_id','var'), profile_id = 'wind03'; end      % windlib id, or a number 1..7
if ~exist('project_root','var'), project_root = pwd; end

% must match the run_prediction settings for this profile, or the forecast
% source will refuse to start (segment indices would not line up)
if ~exist('step_size','var'), step_size = 0.01; end
if ~exist('history_sec','var'), history_sec = 15.0; end
if ~exist('predict_sec','var'), predict_sec = 1.5; end
if ~exist('stride_sec','var'), stride_sec = 1.5; end

if ~exist('forecast_mode','var'), forecast_mode = 'bestof'; end    % 'bestof' | 'model' | 'truth'
if ~exist('forecast_model','var'), forecast_model = 'rbf'; end       % used only when forecast_mode == 'model'
if ~exist('on_missing','var'), on_missing = 'persistence'; end   % 'persistence' | 'bestof' | 'error'

if ~exist('dt_controller','var'), dt_controller = 0.01; end        % controller step. Equal to step_size means the
                              % segment resampling is an exact identity.

if ~exist('max_iter','var'), max_iter = 50; end          % DDP iterations per segment
if ~exist('selection_criterion','var'), selection_criterion = 'RMSE'; end % 'RMSE' | 'COST'

if ~exist('max_segments','var'), max_segments = Inf; end         % cap for a quick trial run, e.g. 5
if ~exist('verbose_iter','var'), verbose_iter = false; end       % per-iteration prints (noisy over many segments)

% ---- export ----
if ~exist('export_web_assets','var'), export_web_assets = true; end
if ~exist('make_gifs','var'), make_gifs = true; end
if ~exist('gif_frames','var'), gif_frames = 40; end

%% ======================= SETUP =======================
t_run = tic;
fprintf('=====================================================\n');
fprintf(' DDP-HJB CONTROL  --  %s  [%s]\n', profile_id, forecast_mode);
fprintf('=====================================================\n');

[wind, T_dur, ~, winfo] = windlib.load(profile_id, step_size);
wind = wind(:);
folder_id = winfo.folder;

P = report.paths(project_root, folder_id, struct('create', true));

fsrc = predictors.cached(P.root, wind, step_size, history_sec, predict_sec, stride_sec, ...
    struct('mode', forecast_mode, 'model', forecast_model, ...
           'on_missing', on_missing, 'verbose', true));

num_segments = min(fsrc.num_segments, max_segments);
p = ddp.params(struct('max_iter', max_iter, 'selection_criterion', selection_criterion));

fprintf('\n segments   : %d\n', num_segments);
fprintf(' dt_control : %g  (wind dt = %g)\n', dt_controller, step_size);
fprintf(' max_iter   : %d\n\n', p.max_iter);

%% ======================= SEGMENT LOOP =======================
best_time_all = [];
best_Wr_all   = [];
best_Tg_all   = [];
best_U_all    = [];
best_lam_all  = [];
best_Cp_all   = [];
Wr_ref_all    = [];
Tg_ref_all    = [];
wind_pred_all = [];

segment_times = zeros(1, num_segments);
segment_span  = zeros(2, num_segments);
RMSE_all      = zeros(1, num_segments);
cost_all      = zeros(1, num_segments);
ddp_time      = zeros(1, num_segments);
n_iter_all    = zeros(1, num_segments);
fmodel_all    = cell(1, num_segments);

carryover_state = [];
carryover_U     = [];

t_offset = fsrc.t_offset;

for seg = 1:num_segments
    span = fsrc.seg_idx(seg);
    segment_span(:, seg) = span(:);
    segment_times(seg)   = (span(1) - 1) * step_size;

    % ---- forecast for this segment ----
    [future_pred, finfo] = fsrc.get(seg);
    fmodel_all{seg} = finfo.model;

    % ---- onto the controller grid ----
    wind_seg = ddp.resample_segment(future_pred, step_size, dt_controller)';

    % ---- solve ----
    t_seg = tic;
    out = ddp.solve(wind_seg, carryover_state, carryover_U, p, ...
        dt_controller, step_size, struct('verbose', verbose_iter));
    ddp_time(seg) = toc(t_seg);

    % ---- accumulate on a continuous timeline ----
    L = numel(out.Wr);
    t_seg_vec = segment_times(seg) + (0:L-1) * dt_controller;

    best_time_all = [best_time_all, t_seg_vec];      %#ok<AGROW>
    best_Wr_all   = [best_Wr_all,   out.Wr];         %#ok<AGROW>
    best_Tg_all   = [best_Tg_all,   out.Tg];         %#ok<AGROW>
    best_U_all    = [best_U_all,    out.U];          %#ok<AGROW>
    best_lam_all  = [best_lam_all,  out.lambda];     %#ok<AGROW>
    best_Cp_all   = [best_Cp_all,   out.Cp];         %#ok<AGROW>
    Wr_ref_all    = [Wr_ref_all,    out.ref.Wr_opt]; %#ok<AGROW>
    Tg_ref_all    = [Tg_ref_all,    out.ref.Tg_opt]; %#ok<AGROW>
    wind_pred_all = [wind_pred_all, wind_seg];       %#ok<AGROW>

    RMSE_all(seg)  = out.RMSE;
    cost_all(seg)  = out.cost;
    n_iter_all(seg)= out.n_iter;

    carryover_state = out.carryover_state;
    carryover_U     = out.carryover_U;

    if isempty(finfo.fallback)
        tag = finfo.model;
    else
        tag = sprintf('%s->%s', finfo.model, finfo.fallback);
    end
    fprintf(' seg %3d/%d | %-14s | RMSE %8.4f | %2d iter | %5.2f s\n', ...
        seg, num_segments, tag, out.RMSE, out.n_iter, ddp_time(seg));
end

fprintf('\n');

%% ======================= SUMMARY =======================
err = best_Wr_all - Wr_ref_all;
fprintf('=================== CONTROL SUMMARY ===================\n');
fprintf(' forecast source : %s\n', fsrc.label);
fprintf(' segments        : %d\n', num_segments);
fprintf(' timeline        : %.2f .. %.2f s (starts at history_sec)\n', ...
    best_time_all(1), best_time_all(end));
fprintf(' RMSE  (Wr)      : %.4f rad/s\n', sqrt(mean(err.^2)));
fprintf(' max |error|     : %.4f rad/s\n', max(abs(err)));
fprintf(' mean per-segment RMSE : %.4f\n', mean(RMSE_all));
fprintf(' energy          : %.2f kJ\n', sum(abs(best_Tg_all .* best_Wr_all))*dt_controller*1e-3);

% seam discontinuity -- the number that quantifies the stitching artefacts
seam = [];
for i = 2:num_segments
    k = find(best_time_all >= segment_times(i), 1, 'first');
    if ~isempty(k) && k > 1, seam(end+1) = abs(best_Wr_all(k) - best_Wr_all(k-1)); end %#ok<SAGROW>
end
if ~isempty(seam)
    fprintf(' seam jump       : mean %.4f, max %.4f rad/s\n', mean(seam), max(seam));
end

%% ======================= REAL-TIME BUDGET =======================
fprintf('\n=================== REAL-TIME BUDGET ===================\n');
mean_ddp = mean(ddp_time);
fprintf(' horizon              : %.2f s per segment\n', predict_sec);
fprintf(' DDP per segment      : %.3f s (mean over %d segments, %d iter)\n', ...
    mean_ddp, num_segments, round(mean(n_iter_all)));

pred_cost = NaN;
if exist(P.mat_file, 'file') == 2
    Sm = load(P.mat_file, 'summary_table');
    if isfield(Sm,'summary_table')
        if strcmpi(forecast_mode,'model')
            r = find(strcmpi(Sm.summary_table.Model, forecast_model), 1);
        else
            r = 1;   % the best model, for a like-for-like figure
        end
        if ~isempty(r), pred_cost = Sm.summary_table.Time_sec(r); end
    end
end
if ~isnan(pred_cost)
    fprintf(' prediction per segment: %.3f s\n', pred_cost);
    fprintf(' total per segment     : %.3f s\n', mean_ddp + pred_cost);
    if mean_ddp + pred_cost <= predict_sec
        fprintf(' => FITS inside the horizon (%.0f%% of budget)\n', ...
            100*(mean_ddp+pred_cost)/predict_sec);
    else
        fprintf(2, ' => %.2fx OVER the horizon budget\n', (mean_ddp+pred_cost)/predict_sec);
    end
else
    fprintf(' (prediction cost unavailable -- run run_prediction to record it)\n');
end

%% ======================= EXPORT =======================
fprintf('\n=================== EXPORT ===================\n');

t_wind = (0:numel(wind)-1) * step_size;

C = struct();
C.time          = best_time_all;
C.Wr            = best_Wr_all;
C.Wr_ref        = Wr_ref_all;
C.Tg            = best_Tg_all;
C.Tg_ref        = Tg_ref_all;
C.U             = best_U_all;
C.lambda        = best_lam_all;
C.Cp            = best_Cp_all;
C.wind_time     = t_wind;
C.wind_true     = wind(:)';
C.wind_pred     = wind_pred_all;
C.segment_times = segment_times;
C.RMSE_all      = RMSE_all;
C.cost_all      = cost_all;
C.forecast_mode = forecast_mode;
C.forecast_model= forecast_model;
C.dt_controller = dt_controller;
C.settings      = struct( ...
    'history_sec', history_sec, 'predict_sec', predict_sec, ...
    'stride_sec', stride_sec, 'step_size', step_size, ...
    'max_iter', p.max_iter, 'selection_criterion', selection_criterion, ...
    't_control_offset', t_offset, ...
    'ddp_sec_per_segment', mean_ddp, ...
    'predict_sec_per_segment', pred_cost);

if strcmpi(forecast_mode, 'model')
    run_label = forecast_model;
else
    run_label = forecast_mode;
end

report.export_control(P.root, folder_id, run_label, C, ...
    struct('export_web', export_web_assets, 'make_gifs', make_gifs, ...
           'gif_frames', gif_frames));

report.build_data(project_root);

%% ======================= SAVE =======================
ctrl_mat = fullfile(P.control, run_label, 'control_results.mat');
save(ctrl_mat, 'C', 'segment_span', 'ddp_time', 'n_iter_all', 'fmodel_all', ...
     'RMSE_all', 'cost_all', 'p', 'forecast_mode', 'forecast_model', ...
     'dt_controller', 'num_segments');

fprintf('\nDone in %.1f s.\n', toc(t_run));
fprintf('  %s\n', fullfile(P.control, run_label));
fprintf('  open report.html -- section 5 now has the "%s" run\n', run_label);
