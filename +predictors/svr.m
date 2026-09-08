function [future_pred, metrics, segment_time_sec] = svr( ...
    wind_data, step_size, seg, history_sec, predict_sec, stride_sec, ...
    show_plots, show_text)
%PREDICTORS.SVR  Segment predictor -- ported verbatim from predict_with_svr_segment
%
% Contract (identical for every predictor in this package):
%   [future_pred, metrics, segment_time_sec] = predictors.svr( ...
%       wind_data, step_size, seg, history_sec, predict_sec, stride_sec, ...
%       show_plots, show_text)
%
% Body is unchanged from the thesis code. Only the function name and, where
% noted, the RBF helper calls have moved.
% =========================================================================
% SVR Wind Predictor (per-segment, AR with lags)
% Requires: Statistics and Machine Learning Toolbox (fitrsvm).
% =========================================================================

% ---------- indices ----------
dt               = step_size;
window_size      = round(history_sec  / dt);
segment_length   = round(predict_sec  / dt);
segment_stride   = round(stride_sec   / dt);

t = (0:length(wind_data)-1) * dt;
idx_start  = (seg - 1) * segment_stride + 1;
idx_end    = idx_start + window_size - 1;
future_end = idx_end + segment_length;

if future_end > length(wind_data)
    future_pred = nan(segment_length,1);
    metrics = struct('RMSE',nan,'MAE',nan,'MAPE',nan,'R2',nan,'time',0);
    segment_time_sec = 0;
    return;
end

past_segment  = wind_data(idx_start:idx_end);
future_actual = wind_data(idx_end+1:future_end);

% ---------- fit window ----------
fit_window = max(4*window_size, 6*segment_length);
fit_start  = max(1, idx_end - fit_window + 1);
y_fit_full = wind_data(fit_start:idx_end);
y_fit_full = y_fit_full(:);

% scrub non-finite
if any(~isfinite(y_fit_full))
    y_fit_full = y_fit_full(isfinite(y_fit_full));
end

Tfit = numel(y_fit_full);
if Tfit < (window_size + segment_length + 10)
    future_pred = repmat(wind_data(idx_end), segment_length, 1);
    segment_time_sec = 0;
    metrics = local_metrics(future_actual, future_pred, 0);
    if show_text
        fprintf('[SVR] seg %d: too little data -> persistence.\n', seg);
    end
    local_plots();
    return;
end

% ---------- normalization (z-score on fit window) ----------
mu  = mean(y_fit_full, 'omitnan');
sig = std(y_fit_full, 'omitnan') + 1e-9;
y_norm = (y_fit_full - mu) / sig;

% ---------- autoregressive design (lags as features) ----------
p = min(24, max(6, round(0.2*window_size)));   % lag order
Ntr = Tfit - p;
if Ntr < 30
    future_pred = repmat(wind_data(idx_end), segment_length, 1);
    segment_time_sec = 0;
    metrics = local_metrics(future_actual, future_pred, 0);
    if show_text
        fprintf('[SVR] seg %d: not enough AR rows (N=%d) -> persistence.\n', seg, Ntr);
    end
    local_plots();
    return;
end

Xtr = zeros(Ntr, p);
Ytr = zeros(Ntr, 1);
for i = 1:Ntr
    Xtr(i,:) = y_norm(i:i+p-1).';
    Ytr(i)   = y_norm(i+p);
end

good = all(isfinite(Xtr),2) & isfinite(Ytr);
Xtr = Xtr(good,:); 
Ytr = Ytr(good);

if size(Xtr,1) < 25
    future_pred = repmat(wind_data(idx_end), segment_length, 1);
    segment_time_sec = 0;
    metrics = local_metrics(future_actual, future_pred, 0);
    if show_text
        fprintf('[SVR] seg %d: too few clean rows -> persistence.\n', seg);
    end
    local_plots();
    return;
end

% ---------- train SVR (RBF kernel) ----------
% Initial heuristics for BoxConstraint and Epsilon
boxC   = 10;          % try 5..50 if underfitting
epsVal = 0.05;        % try 0.02..0.10 depending on noise
kerScale = 'auto';    % or a numeric value

try
    tic;
    svrMdl = fitrsvm(Xtr, Ytr, ...
        'KernelFunction','gaussian', ...
        'KernelScale', kerScale, ...
        'BoxConstraint', boxC, ...
        'Epsilon', epsVal, ...
        'Standardize', false, ...
        'Solver','SMO', ...
        'Verbose', 0);
    segment_time_sec = toc;
catch ME
    % fallback with simplified options
    if show_text
        warning('[SVR] seg %d: fitrsvm failed (%s). Using coarse defaults.', seg, ME.message);
    end
    tic;
    svrMdl = fitrsvm(Xtr, Ytr, 'KernelFunction','gaussian', 'Standardize', false);
    segment_time_sec = toc;
end

% ---------- recursive rollout (normalized space) ----------
past_norm = (past_segment(:) - mu) / sig;
if numel(past_norm) < p
    start_vec = [repmat(past_norm(1), p-numel(past_norm), 1); past_norm];
else
    start_vec = past_norm(end-p+1:end);
end

xf = start_vec(:).';
y_future_norm = zeros(segment_length,1);

for k = 1:segment_length
    yhat = predict(svrMdl, xf);
    y_future_norm(k) = yhat;
    xf = [xf(2:end) yhat];
end

% ---------- denormalize ----------
future_pred = y_future_norm * sig + mu;

% ---------- metrics + plots ----------
metrics = local_metrics(future_actual, future_pred, segment_time_sec);
if show_text
    fprintf('[SVR] seg %d  RMSE=%.4f  R2=%.4f  time=%.2fs  (p=%d, rows=%d)\n', ...
        seg, metrics.RMSE, metrics.R2, metrics.time, p, size(Xtr,1));
end
local_plots();

% ================= helpers =================
    function m = local_metrics(y_true, y_hat, tsec)
        y_true = y_true(:); 
        y_hat  = y_hat(:);
        err    = y_true - y_hat;
        m.RMSE = sqrt(mean(err.^2));
        m.MAE  = mean(abs(err));
        m.MAPE = mean(abs(err) ./ max(abs(y_true),1e-9)) * 100;
        m.R2   = 1 - sum(err.^2) / max(sum((y_true - mean(y_true)).^2), eps);
        m.time = tsec;
    end

    function local_plots()
        if ~show_plots, return; end
        fig = figure(292); clf(fig); set(fig,'Name',sprintf('SVR Segment %d', seg));

        % full view
        subplot(2,1,1);
        plot(t, wind_data, 'k-', 'LineWidth', 1); hold on;
        xline(t(idx_start), 'b--'); 
        xline(t(idx_end), 'b--');
        xline(t(idx_end+1), 'r--'); 
        xline(t(future_end), 'r--');
        plot(t(idx_end+1:future_end), future_pred, 'm.-', 'LineWidth', 1.2);
        title(sprintf('Full Wind View - Segment %d (SVR-AR)', seg));
        legend('Wind Data','Train Start','Train End','Pred Start','Pred End','SVR Prediction');
        xlim([t(max(1,idx_start-5)), t(min(length(t),future_end+5))]);
        grid on; grid minor;

        % zoomed
        subplot(2,1,2);
        plot(t(idx_end+1:future_end), future_actual, 'k-', 'LineWidth', 1.5); hold on;
        plot(t(idx_end+1:future_end), future_pred,  'm--', 'LineWidth', 1.5);
        title('Zoomed-In Segment View');
        xlabel('Time (s)'); ylabel('Wind Speed (m/s)');
        legend('True Wind','SVR Prediction');
        grid on; grid minor;
    end
end
