%% ========================================================================
%  RUN_ALL  --  drive the whole pipeline unattended
%
%      >> cd C:\path\to\proj1
%      >> run_all
%
%  Loops over profiles and forecast modes with no prompts, skipping work
%  that is already done. Safe to re-run: it picks up where it left off.
%
%  A full sweep is long. Start with one profile and one mode, confirm the
%  report looks right, then widen it.
% ========================================================================
clc; close all;
clearvars

%% ======================= CONFIG =======================
project_root = pwd;

% profiles = {'wind04'};
% profiles = {'wind02','wind03','wind05','wind07'};        % the 60-100 s ones
profiles = {'wind01','wind02','wind03','wind04','wind05','wind07'};  % all feasible

step_size   = 0.05;
history_sec = 15.0;
predict_sec = 1.5;
stride_sec  = 1.5;

models_to_run = {'rbf','rbf_arima','arima','svr','gp','esn','gru','lstm','bilstm','tcn'};
models_to_run = {'rbf','rbf_arima','arima'};

% Forecast sources for the control stage. 'truth' is the perfect-foresight
% upper bound; the gap between it and the others is the cost of forecast error.
control_runs = { ...
    struct('mode','truth',  'model',''), ...
    struct('mode','bestof', 'model',''), ...
    struct('mode','model',  'model','rbf') };

dt_controller = 0.01;
max_iter      = 50;
max_segments  = 5;      % set small (e.g. 5) for a smoke run

do_prediction = true;
do_control    = true;
skip_done     = true;     % skip a stage whose output already exists

%% ======================= FEASIBILITY =======================
fprintf('=====================================================\n');
fprintf(' RUN ALL\n');
fprintf('=====================================================\n\n');
windlib.feasibility(step_size, history_sec, predict_sec, stride_sec);

keep = true(1, numel(profiles));
for i = 1:numel(profiles)
    if windlib.feasibility(step_size, history_sec, predict_sec, stride_sec, profiles{i}) < 1
        fprintf(2, 'Skipping %s -- no segments at this configuration.\n', profiles{i});
        keep(i) = false;
    end
end
profiles = profiles(keep);
if isempty(profiles)
    error('run_all:noProfiles', 'No profile yields any segment. Shorten history_sec.');
end

fprintf('Profiles : %s\n', strjoin(profiles, ', '));
fprintf('Models   : %d\n', numel(models_to_run));
fprintf('Control  : %d run(s) per profile\n\n', numel(control_runs));

t_all = tic;
log = {};

%% ======================= PREDICTION =======================
if do_prediction
    for i = 1:numel(profiles)
        pid = profiles{i};
        info = windlib.registry(pid);
        P = report.paths(project_root, info.folder);

        if skip_done && exist(P.mat_file, 'file') == 2
            fprintf('[skip] prediction for %s (results already exist)\n', pid);
            log{end+1} = sprintf('%-8s prediction  SKIPPED', pid); %#ok<SAGROW>
            continue;
        end

        fprintf('\n>>> PREDICTION: %s\n', pid);
        t0 = tic;
        try
            run_prediction_for(project_root, pid, step_size, history_sec, ...
                predict_sec, stride_sec, models_to_run);
            log{end+1} = sprintf('%-8s prediction  %.0f s', pid, toc(t0)); %#ok<SAGROW>
        catch ME
            fprintf(2, '  FAILED: %s\n', ME.message);
            log{end+1} = sprintf('%-8s prediction  FAILED (%s)', pid, ME.message); %#ok<SAGROW>
        end
    end
end

%% ======================= CONTROL =======================
if do_control
    for i = 1:numel(profiles)
        pid = profiles{i};
        info = windlib.registry(pid);
        P = report.paths(project_root, info.folder);

        if exist(P.mat_file, 'file') ~= 2
            fprintf(2, '[skip] control for %s -- no prediction results yet\n', pid);
            continue;
        end

        for r = 1:numel(control_runs)
            cr = control_runs{r};
            if strcmpi(cr.mode, 'model')
                label = cr.model;
            else
                label = cr.mode;
            end

            done_marker = fullfile(P.control, label, 'control_results.mat');
            if skip_done && exist(done_marker, 'file') == 2
                fprintf('[skip] control %s/%s (already done)\n', pid, label);
                log{end+1} = sprintf('%-8s control %-8s SKIPPED', pid, label); %#ok<SAGROW>
                continue;
            end

            fprintf('\n>>> CONTROL: %s [%s]\n', pid, label);
            t0 = tic;
            try
                run_control_for(project_root, pid, step_size, history_sec, predict_sec, ...
                    stride_sec, cr.mode, cr.model, dt_controller, max_iter, max_segments);
                log{end+1} = sprintf('%-8s control %-8s %.0f s', pid, label, toc(t0)); %#ok<SAGROW>
            catch ME
                fprintf(2, '  FAILED: %s\n', ME.message);
                log{end+1} = sprintf('%-8s control %-8s FAILED (%s)', pid, label, ME.message); %#ok<SAGROW>
            end
        end
    end
end

%% ======================= BUNDLE + SUMMARY =======================
report.build_data(project_root);

fprintf('\n=====================================================\n');
fprintf(' RUN LOG\n');
fprintf('=====================================================\n');
for i = 1:numel(log)
    fprintf('  %s\n', log{i});
end
fprintf('\nTotal: %.1f min\n', toc(t_all)/60);
fprintf('Open report.html to view everything.\n');

%% ========================================================================
%  Thin wrappers. Each sets the variables run_prediction / run_control read,
%  then runs the script inside this function's own workspace. The config
%  blocks in those scripts are guarded with  if ~exist(...)  so whatever is
%  set here wins, and each call starts from a clean scope.
% ========================================================================
function run_prediction_for(project_root, profile_id, step_size, history_sec, ...
        predict_sec, stride_sec, models_to_run)
    reuse_existing = false;
    use_parallel = false;
    selection_criterion = 'RMSE';
    export_web_assets = true;
    make_gifs = true;
    export_segment_figures = false;
    export_stats_figures = true;
    zoom_types = {'actual_vs_pred','residual','parity'};
    gif_fps = 2;
    run('run_prediction.m');
end

function run_control_for(project_root, profile_id, step_size, history_sec, predict_sec, ...
        stride_sec, forecast_mode, forecast_model, dt_controller, max_iter, max_segments)
    on_missing = 'persistence';
    selection_criterion = 'RMSE';
    verbose_iter = false;
    export_web_assets = true;
    make_gifs = true;
    gif_frames = 40;
    run('run_control.m');
end
