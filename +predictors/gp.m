function [future_pred, metrics, segment_time_sec] = gp( ...
    wind_data, step_size, seg, history_sec, predict_sec, stride_sec, ...
    show_plots, show_text)
%PREDICTORS.GP  Segment predictor -- ported verbatim from predict_with_gp_segment
%
% Contract (identical for every predictor in this package):
%   [future_pred, metrics, segment_time_sec] = predictors.gp( ...
%       wind_data, step_size, seg, history_sec, predict_sec, stride_sec, ...
%       show_plots, show_text)
%
% Body is unchanged from the thesis code. Only the function name and, where
% noted, the RBF helper calls have moved.
% =========================================================================
% Gaussian Process Wind Predictor — per-segment, RBF-compatible I/O
% - Autoregressive GP (lags as features), one-step recursive rollout
% - Z-score normalization on a longer fit window
% - Matern 5/2 kernel, compact optimization budget
% - Full + zoomed plots, robust fallbacks
% Requires: Statistics and Machine Learning Toolbox (fitrgp).
% =========================================================================

% ---------- indices ----------
dt               = step_size;
window_size      = round(history_sec  / dt);
segment_length   = round(predict_sec  / dt);
segment_stride   = round(stride_sec   / dt);

t = (0:length(wind_data)-1)*dt;
idx_start  = (seg - 1)*segment_stride + 1;
idx_end    = idx_start + window_size - 1;
future_end = idx_end   + segment_length;

if future_end > length(wind_data)
    future_pred = nan(segment_length,1);
    metrics = struct('RMSE',nan,'MAE',nan,'MAPE',nan,'R2',nan,'time',0);
    segment_time_sec = 0; return;
end

past_segment  = wind_data(idx_start:idx_end);
future_actual = wind_data(idx_end+1:future_end);

% ---------- build a longer fit window ----------
fit_window = max(4*window_size, 6*segment_length);
fit_start  = max(1, idx_end - fit_window + 1);
y_fit_full = wind_data(fit_start:idx_end);
y_fit_full = y_fit_full(:);

% scrub
if any(~isfinite(y_fit_full))
    good = isfinite(y_fit_full);
    y_fit_full = y_fit_full(good);
end

Tfit = numel(y_fit_full);
if Tfit < (window_size + segment_length + 10)
    % not enough clean data -> persistence
    future_pred = repmat(wind_data(idx_end), segment_length, 1);
    segment_time_sec = 0;
    metrics = local_metrics(future_actual, future_pred, 0);
    if show_text
        fprintf('[GP] seg %d: too little data → persistence.\n', seg);
    end
    local_plots(); return;
end

% ---------- normalization (z-score on fit window) ----------
mu  = mean(y_fit_full, 'omitnan');
sig = std(y_fit_full, 'omitnan') + 1e-9;
y_norm = (y_fit_full - mu) / sig;

% ---------- autoregressive design (lags as features) ----------
% choose lag order proportional to window, but capped
p = min(20, max(5, round(0.15*window_size)));
Ntr = Tfit - p;                % number of training rows
if Ntr < 30
    % fallback
    future_pred = repmat(wind_data(idx_end), segment_length, 1);
    segment_time_sec = 0;
    metrics = local_metrics(future_actual, future_pred, 0);
    if show_text
        fprintf('[GP] seg %d: not enough AR rows (N=%d) → persistence.\n', seg, Ntr);
    end
    local_plots(); return;
end

Xtr = zeros(Ntr, p);
Ytr = zeros(Ntr, 1);
for i = 1:Ntr
    Xtr(i,:) = y_norm(i : i+p-1).';
    Ytr(i)   = y_norm(i+p);
end

% final scrub
good = all(isfinite(Xtr),2) & isfinite(Ytr);
Xtr = Xtr(good,:); Ytr = Ytr(good);
if size(Xtr,1) < 25
    future_pred = repmat(wind_data(idx_end), segment_length, 1);
    segment_time_sec = 0;
    metrics = local_metrics(future_actual, future_pred, 0);
    if show_text
        fprintf('[GP] seg %d: too few clean GP rows → persistence.\n', seg);
    end
    local_plots(); return;
end

% ---------- fit GP (Matern 5/2, modest HPO) ----------
kernel = 'matern52';
optsHPO = struct( ...
    'ShowPlots', false, ...
    'Verbose', 0, ...
    'MaxObjectiveEvaluations', min(30, 5 + 2*p));   % keep it snappy

try
    tic;
    gprMdl = fitrgp(Xtr, Ytr, ...
        'KernelFunction', kernel, ...
        'BasisFunction', 'constant', ...
        'SigmaLowerBound', 1e-6, ...
        'Standardize', false, ...       % we normalized ourselves
        'OptimizeHyperparameters', {'KernelScale','Sigma','BasisFunction'}, ...
        'HyperparameterOptimizationOptions', optsHPO);
    segment_time_sec = toc;
catch ME
    % fallback without HPO
    warning('[GP] seg %d: HPO failed (%s). Using default hyperparams.', seg, ME.message);
    tic;
    gprMdl = fitrgp(Xtr, Ytr, ...
        'KernelFunction', kernel, ...
        'BasisFunction', 'constant', ...
        'Sigma', 0.1, ...
        'KernelScale', 'auto', ...
        'Standardize', false);
    segment_time_sec = toc;
end

% ---------- recursive rollout (normalized space) ----------
% start from the last p normalized samples from the *past segment*
past_norm_all = (past_segment(:) - mu) / sig;
if numel(past_norm_all) < p
    % pad with the earliest value if window shorter than p
    start_vec = [repmat(past_norm_all(1), p-numel(past_norm_all), 1); past_norm_all];
else
    start_vec = past_norm_all(end-p+1:end);
end

xf = start_vec(:).';
y_future_norm = zeros(segment_length,1);

for k = 1:segment_length
    yhat = predict(gprMdl, xf);
    y_future_norm(k) = yhat;
    % shift-in the new prediction
    xf = [xf(2:end) yhat];
end

% ---------- denormalize ----------
future_pred = y_future_norm * sig + mu;

% ---------- metrics + plots ----------
metrics = local_metrics(future_actual, future_pred, segment_time_sec);
if show_text
    fprintf('[GP] seg %d  RMSE=%.4f  R2=%.4f  time=%.2fs  (p=%d, rows=%d)\n', ...
        seg, metrics.RMSE, metrics.R2, metrics.time, p, size(Xtr,1));
end
local_plots();

% ================= helpers =================
    function m = local_metrics(y_true, y_hat, tsec)
        y_true = y_true(:); y_hat = y_hat(:);
        err    = y_true - y_hat;
        m.RMSE = sqrt(mean(err.^2));
        m.MAE  = mean(abs(err));
        m.MAPE = mean(abs(err) ./ max(abs(y_true),1e-9))*100;
        m.R2   = 1 - sum(err.^2)/max(sum((y_true-mean(y_true)).^2), eps);
        m.time = tsec;
    end

    function local_plots()
        if ~show_plots, return; end
        fig = figure(282); clf(fig); set(fig,'Name',sprintf('GP Segment %d', seg));

        % full view
        subplot(2,1,1);
        plot(t, wind_data, 'k-', 'LineWidth', 1); hold on;
        xline(t(idx_start), 'b--'); xline(t(idx_end), 'b--');
        xline(t(idx_end+1), 'r--'); xline(t(future_end), 'r--');
        plot(t(idx_end+1:future_end), future_pred, 'm.-', 'LineWidth', 1.2);
        title(sprintf('Full Wind View — Segment %d (GP-AR)', seg));
        legend('Wind Data','Train Start','Train End','Pred Start','Pred End','GP Prediction');
        xlim([t(max(1,idx_start-5)), t(min(length(t),future_end+5))]);
        grid on; grid minor;

        % zoomed
        subplot(2,1,2);
        plot(t(idx_end+1:future_end), future_actual, 'k-', 'LineWidth', 1.5); hold on;
        plot(t(idx_end+1:future_end), future_pred,  'm--', 'LineWidth', 1.5);
        title('Zoomed-In Segment View');
        xlabel('Time (s)'); ylabel('Wind Speed (m/s)');
        legend('True Wind','GP Prediction');
        grid on; grid minor;
    end
end
