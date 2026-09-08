function [future_pred, metrics, segment_time_sec] = gru( ...
    wind_data, step_size, seg, history_sec, predict_sec, stride_sec, ...
    show_plots, show_text)
%PREDICTORS.GRU  Segment predictor -- ported verbatim from predict_with_gru_segment
%
% Contract (identical for every predictor in this package):
%   [future_pred, metrics, segment_time_sec] = predictors.gru( ...
%       wind_data, step_size, seg, history_sec, predict_sec, stride_sec, ...
%       show_plots, show_text)
%
% Body is unchanged from the thesis code. Only the function name and, where
% noted, the RBF helper calls have moved.
% =========================================================================
% GRU Reservoir + Linear Readout (AR decoding, residual option)
% RBF-compatible I/O: returns future_pred (column), metrics, segment_time_sec
% =========================================================================

% ------------------------ Setup & indexing (as RBF) ------------------------
dt               = step_size;
window_size      = round(history_sec  / dt);
segment_length   = round(predict_sec  / dt);
segment_duration = round(stride_sec   / dt);

t = (0:length(wind_data)-1) * dt;
idx_start  = (seg - 1)*segment_duration + 1;
idx_end    = idx_start + window_size - 1;
future_end = idx_end   + segment_length;

if future_end > length(wind_data)
    future_pred = nan(segment_length,1);
    metrics = struct('RMSE', nan, 'MAE', nan, 'MAPE', nan, 'R2', nan, 'time', 0);
    segment_time_sec = 0; return;
end

past_segment   = wind_data(idx_start:idx_end);
future_actual  = wind_data(idx_end+1:future_end);

% ------------------------ Normalization (min-max on entire series) --------
min_val = min(wind_data);
max_val = max(wind_data);
scale   = max(max_val - min_val, eps);

x_seq   = (past_segment(:)  - min_val) / scale;          % T×1 in [0,1]
y_next  = (future_actual(:) - min_val) / scale;          % H×1 for metrics

% ------------------------ Options ------------------------------------------
use_residuals = true;  % predict deltas in normalized space, then integrate
hidden_size   = 128;
leak          = 0.98;
ridge_lambda  = 1e-2;
clip_to_unit  = true;

% ------------------------ GRU "reservoir" (persistent init) ----------------
persistent Wz Wr Wh Wy by init_done
if isempty(init_done)
    rng(42);
    Wz = 0.05 * randn(hidden_size, 1 + hidden_size);
    Wr = 0.05 * randn(hidden_size, 1 + hidden_size);
    Wh = 0.05 * randn(hidden_size, 1 + hidden_size);
    % tame recurrent part a bit
    scale_rec = 0.9;
    Wz(:,2:end) = Wz(:,2:end) * (scale_rec / max(norm(Wz(:,2:end)),1e-6));
    Wr(:,2:end) = Wr(:,2:end) * (scale_rec / max(norm(Wr(:,2:end)),1e-6));
    Wh(:,2:end) = Wh(:,2:end) * (scale_rec / max(norm(Wh(:,2:end)),1e-6));
    % readout (trained per segment, but keep shape persistent)
    Wy = zeros(1, hidden_size);
    by = 0;
    init_done = true;
end
sigmoid = @(x) 1./(1+exp(-x));

% ------------------------ Encode past window -> states ---------------------
tic;
h = zeros(hidden_size,1);
H = zeros(window_size-1, hidden_size);  % state(t) for predicting next step
for tstep = 1:window_size
    xt     = x_seq(tstep);
    concat = [xt; h];
    z      = sigmoid(Wz*concat);
    r      = sigmoid(Wr*concat);
    h_til  = tanh(Wh * [xt; r.*h]);
    h      = (1 - z).*h + z.*h_til;
    h      = leak * h;
    if tstep < window_size
        H(tstep,:) = h.';
    end
end

% targets for readout
if use_residuals
    % target is delta: x(t+1) - x(t), length = window_size-1
    Y = diff(x_seq);                            % [T-1×1]
else
    % target is absolute next: x(t+1)
    Y = x_seq(2:end);                           % [T-1×1]
end

% ------------------------ Train readout (ridge) ----------------------------
Phi = [H, ones(size(H,1),1)];
I   = eye(size(Phi,2)); I(end,end)=0;          % don't regularize bias
theta = (Phi.'*Phi + ridge_lambda*I) \ (Phi.'*Y);
Wy = theta(1:hidden_size).';
by = theta(end);

% ------------------------ Autoregressive rollout ---------------------------
future_pred_norm = zeros(segment_length,1);

% start from last known sample in normalized space:
x_last = x_seq(end);  % last past value (normalized)
% h is already the hidden state after consuming x_seq(end)

for k = 1:segment_length
    yhat = Wy*h + by;       % either delta or absolute next in [~0..1]
    if use_residuals
        x_next = x_last + yhat;       % integrate residual
    else
        x_next = yhat;                 % absolute
    end
    if clip_to_unit
        x_next = min(max(x_next, 0), 1);
    end
    future_pred_norm(k) = x_next;

    % feed back x_next and update hidden
    xt     = x_next;
    concat = [xt; h];
    z      = sigmoid(Wz*concat);
    r      = sigmoid(Wr*concat);
    h_til  = tanh(Wh * [xt; r.*h]);
    h      = (1 - z).*h + z.*h_til;
    h      = leak * h;

    x_last = x_next;
end

% ------------------------ Denormalize & metrics ----------------------------
future_pred = future_pred_norm*scale + min_val;
segment_time_sec = toc;

err = future_actual(:) - future_pred(:);
metrics.RMSE = sqrt(mean(err.^2));
metrics.MAE  = mean(abs(err));
metrics.MAPE = mean(abs(err) ./ max(abs(future_actual(:)), 1e-9)) * 100;
metrics.R2   = 1 - sum(err.^2) / max(sum((future_actual(:) - mean(future_actual(:))).^2), eps);
metrics.time = segment_time_sec;

if show_text
    fprintf('[GRU-AR%s] seg %d  RMSE=%.4f  R2=%.4f  time=%.3fs\n', ...
        ternary(use_residuals,'-res',''), seg, metrics.RMSE, metrics.R2, metrics.time);
end

% ------------------------ Plots: FULL + ZOOMED -----------------------------
if show_plots
    fig = figure(202); clf(fig); set(fig,'Name',sprintf('GRU Segment %d', seg));

    % FULL VIEW (like RBF)
    subplot(2,1,1);
    plot(t, wind_data, 'k-', 'LineWidth', 1); hold on;
    xline(t(idx_start), 'b--'); xline(t(idx_end), 'b--');
    xline(t(idx_end+1), 'r--'); xline(t(future_end), 'r--');
    plot(t(idx_end+1:future_end), future_pred, 'm.-', 'LineWidth', 1.2);
    title(sprintf('Full Wind View — Segment %d', seg));
    legend('Wind Data','Train Start','Train End','Pred Start','Pred End','GRU Prediction');
    xlim([t(max(1,idx_start-5)) , t(min(length(t),future_end+5))]);
    grid on; grid minor;

    % ZOOMED VIEW
    subplot(2,1,2);
    plot(t(idx_end+1:future_end), future_actual, 'k-', 'LineWidth', 1.5); hold on;
    plot(t(idx_end+1:future_end), future_pred,  'm--', 'LineWidth', 1.5);
    title('Zoomed-In Segment View');
    xlabel('Time (s)'); ylabel('Wind Speed (m/s)');
    legend('True Wind','GRU Prediction');
    grid on; grid minor;
end
end

% tiny helper
function out = ternary(cond, a, b), if cond, out=a; else, out=b; end, end
