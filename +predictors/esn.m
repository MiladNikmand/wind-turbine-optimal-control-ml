function [future_pred, metrics, segment_time_sec] = esn( ...
    wind_data, step_size, seg, history_sec, predict_sec, stride_sec, ...
    show_plots, show_text)
%PREDICTORS.ESN  Segment predictor -- ported verbatim from predict_with_esn_segment
%
% Contract (identical for every predictor in this package):
%   [future_pred, metrics, segment_time_sec] = predictors.esn( ...
%       wind_data, step_size, seg, history_sec, predict_sec, stride_sec, ...
%       show_plots, show_text)
%
% Body is unchanged from the thesis code. Only the function name and, where
% noted, the RBF helper calls have moved.
% =========================================================================
% ESN Wind Predictor — One segment at a time (RBF-compatible I/O)
% - Reservoir fixed (persistent), readout trained per segment (ridge)
% - Teacher-forced next-step training, AR rollout for horizon
% - Full + zoomed plots like your RBF
% =========================================================================

% -------------------- indexing (same as RBF) --------------------
dt               = step_size;
window_size      = round(history_sec  / dt);
segment_length   = round(predict_sec  / dt);
segment_duration = round(stride_sec   / dt);

t = (0:length(wind_data)-1)*dt;
idx_start  = (seg - 1)*segment_duration + 1;
idx_end    = idx_start + window_size - 1;
future_end = idx_end   + segment_length;

if future_end > length(wind_data)
    future_pred = nan(segment_length,1);
    metrics = struct('RMSE',nan,'MAE',nan,'MAPE',nan,'R2',nan,'time',0);
    segment_time_sec = 0; return;
end

past_segment  = wind_data(idx_start:idx_end);
future_actual = wind_data(idx_end+1:future_end);

% -------------------- normalization (min-max) -------------------
min_val = min(wind_data);
max_val = max(wind_data);
scale   = max(max_val - min_val, eps);

x_seq   = (past_segment(:)  - min_val) / scale;    % [T×1], T=window_size
y_true  = (future_actual(:) - min_val) / scale;    % for metrics only

% -------------------- ESN params (tweak later) ------------------
reservoir_size  = 600;     % try 600–1200
spectral_radius = 0.9;     % 0.8–0.98 (lower = smoother)
input_scale     = 0.8;     % 0.5–1.0
leak            = 0.2;     % 0.1–0.3 (lower = smoother/slower)
ridge_lambda    = 1e-3;    % 1e-4–1e-2
use_residuals   = true;    % predict deltas; usually better
clip_to_unit    = true;    % clamp AR output in [0,1] to avoid drift

% -------------------- persistent reservoir ----------------------
persistent Win W bias init_done
if isempty(init_done)
    rng(42);  % reproducible
    Win  = (rand(reservoir_size,1)-0.5)*2 * input_scale;
    Wraw = sprandn(reservoir_size, reservoir_size, 0.03); % sparse helps speed
    % scale to desired spectral radius
    opts.disp = 0;
    e = eigs(Wraw,1,'lm',opts);
    W = (Wraw / max(abs(e),1e-6)) * spectral_radius;
    bias = randn(reservoir_size,1) * 0.01;
    init_done = true;
end

% -------------------- collect states over window ----------------
tic;
x = zeros(reservoir_size,1);
H = zeros(window_size-1, reservoir_size); % state at t for predicting x(t+1)
for k = 1:window_size
    u  = x_seq(k);
    x  = (1 - leak)*x + leak * tanh(Win*u + W*x + bias);
    if k < window_size
        H(k,:) = x.';
    end
end

% targets
if use_residuals
    Y = diff(x_seq);           % x(t+1) - x(t), length T-1
else
    Y = x_seq(2:end);          % absolute next
end

% -------------------- train readout (ridge) ---------------------
Phi = [H, ones(size(H,1),1)];          % add bias column
I   = eye(size(Phi,2)); I(end,end)=0;  % don't regularize bias
theta = (Phi.'*Phi + ridge_lambda*I) \ (Phi.'*Y);
Wout = theta(1:reservoir_size).';
bout = theta(end);

% -------------------- AR rollout for horizon --------------------
future_norm = zeros(segment_length,1);
x_last = x_seq(end);  % last known normalized sample; AR seed

for k = 1:segment_length
    % predict next delta/absolute from current state
    yhat = Wout*x + bout;
    if use_residuals
        x_next = x_last + yhat;
    else
        x_next = yhat;
    end
    if clip_to_unit
        x_next = min(max(x_next,0),1);
    end
    future_norm(k) = x_next;

    % advance reservoir with predicted input
    x = (1 - leak)*x + leak * tanh(Win*x_next + W*x + bias);
    x_last = x_next;
end
segment_time_sec = toc;

% -------------------- denorm + metrics ------------------------------------
future_pred = future_norm*scale + min_val;
err = future_actual(:) - future_pred(:);
metrics.RMSE = sqrt(mean(err.^2));
metrics.MAE  = mean(abs(err));
metrics.MAPE = mean(abs(err) ./ max(abs(future_actual(:)),1e-9))*100;
metrics.R2   = 1 - sum(err.^2) / max(sum((future_actual(:)-mean(future_actual(:))).^2),eps);
metrics.time = segment_time_sec;

if show_text
    fprintf('[ESN] seg %d  RMSE=%.4f  R2=%.4f  time=%.3fs\n', seg, metrics.RMSE, metrics.R2, metrics.time);
end

% -------------------- plots: FULL + ZOOMED --------------------------------
if show_plots
    fig = figure(212); clf(fig); set(fig,'Name',sprintf('ESN Segment %d', seg));

    % full view
    subplot(2,1,1);
    plot(t, wind_data, 'k-', 'LineWidth', 1); hold on;
    xline(t(idx_start), 'b--'); xline(t(idx_end), 'b--');
    xline(t(idx_end+1), 'r--'); xline(t(future_end), 'r--');
    plot(t(idx_end+1:future_end), future_pred, 'm.-', 'LineWidth', 1.2);
    title(sprintf('Full Wind View — Segment %d', seg));
    legend('Wind Data','Train Start','Train End','Pred Start','Pred End','ESN Prediction');
    xlim([t(max(1,idx_start-5)) , t(min(length(t),future_end+5))]);
    grid on; grid minor;

    % zoomed
    subplot(2,1,2);
    plot(t(idx_end+1:future_end), future_actual, 'k-', 'LineWidth', 1.5); hold on;
    plot(t(idx_end+1:future_end), future_pred,  'm--', 'LineWidth', 1.5);
    title('Zoomed-In Segment View');
    xlabel('Time (s)'); ylabel('Wind Speed (m/s)');
    legend('True Wind','ESN Prediction');
    grid on; grid minor;
end
end
