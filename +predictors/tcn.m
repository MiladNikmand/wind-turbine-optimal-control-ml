function [future_pred, metrics, segment_time_sec] = tcn( ...
    wind_data, step_size, seg, history_sec, predict_sec, stride_sec, ...
    show_plots, show_text)
%PREDICTORS.TCN  Segment predictor -- ported verbatim from predict_with_tcn_segment
%
% Contract (identical for every predictor in this package):
%   [future_pred, metrics, segment_time_sec] = predictors.tcn( ...
%       wind_data, step_size, seg, history_sec, predict_sec, stride_sec, ...
%       show_plots, show_text)
%
% Body is unchanged from the thesis code. Only the function name and, where
% noted, the RBF helper calls have moved.
% =========================================================================
% TCN Wind Predictor — One segment at a time (RBF-compatible I/O)
% - Causal dilated conv, kernel_size=2
% - Readout trained per segment (ridge), AR rollout for horizon
% - Full + zoomed plots
% =========================================================================

% ------------------------ indexing (as RBF) --------------------------------
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

% ------------------------ normalization -----------------------------------
min_val = min(wind_data);
max_val = max(wind_data);
scale   = max(max_val - min_val, eps);

x_seq   = (past_segment(:)  - min_val) / scale;     % [T×1], T=window_size
y_true  = (future_actual(:) - min_val) / scale;     % for metrics only

% ------------------------ TCN hyperparams ---------------------------------
kernel_size   = 2;
num_filters   = 64;      % 64–128 are good starts
num_layers    = 4;       % dilations: 1,2,4,8
ridge_lambda  = 1e-3;    % readout regularization
use_residuals = true;    % predict deltas, then integrate
clip_to_unit  = true;    % clamp AR outputs to [0,1]

% ------------------------ persistent weights -------------------------------
persistent W b init_done
if isempty(init_done)
    rng(42);
    W = cell(num_layers,1);
    b = cell(num_layers,1);
    for L = 1:num_layers
        % For kernel_size=2: weights(:,1) for x(t), weights(:,2) for x(t-d)
        W{L} = 0.05 * randn(num_filters, kernel_size);
        b{L} = zeros(num_filters,1);
    end
    init_done = true;
end
relu  = @(x) max(0,x);

% ------------------------ helper: forward pass over a sequence -------------
% Returns H_T (T×F) final-layer features for each timestep (causal).
    function H = tcn_forward(x)
        T = numel(x);
        H_l = cell(num_layers,1);
        H_in = x(:).';             % row vector length T
        for L = 1:num_layers
            d = 2^(L-1);           % dilation
            Y = zeros(num_filters, T);
            x_now = H_in;
            x_del = [zeros(1,min(d,T)), H_in(1:max(0,T-d))];  % causal pad
            WX = W{L}(:,1) * x_now + W{L}(:,2) * x_del + b{L};
            Y = relu(WX);
            H_l{L} = Y;            % F×T
            H_in = mean(Y,1);      % simple channel pooling to feed next layer (1×T)
        end
        H = H_l{end}.';            % T×F (final layer)
    end

% ------------------------ Train readout (teacher forcing) ------------------
tic;
H_full = tcn_forward(x_seq);         % T×F
H = H_full(1:end-1,:);               % features at t predict t+1
if use_residuals
    Y = diff(x_seq);                 % x(t+1)-x(t)
else
    Y = x_seq(2:end);
end

Phi = [H, ones(size(H,1),1)];        % add bias
I   = eye(size(Phi,2)); I(end,end)=0;
theta = (Phi.'*Phi + ridge_lambda*I) \ (Phi.'*Y);
Wout = theta(1:end-1);
bout = theta(end);

% ------------------------ AR rollout for horizon ---------------------------
future_norm = zeros(segment_length,1);
x_last = x_seq(end);
x_work = x_seq;                      % will append predictions

for k = 1:segment_length
    Hk = tcn_forward(x_work);        % recompute features (small horizons OK)
    hT = Hk(end,:).';                % last timestep feature (F×1)
    yhat = Wout.'*hT + bout;         % delta or absolute
    if use_residuals
        x_next = x_last + yhat;
    else
        x_next = yhat;
    end
    if clip_to_unit
        x_next = min(max(x_next,0),1);
    end
    future_norm(k) = x_next;
    x_work = [x_work; x_next];       % append for next step
    x_last = x_next;
end
segment_time_sec = toc;

% ------------------------ denorm + metrics ---------------------------------
future_pred = future_norm*scale + min_val;
err = future_actual(:) - future_pred(:);
metrics.RMSE = sqrt(mean(err.^2));
metrics.MAE  = mean(abs(err));
metrics.MAPE = mean(abs(err) ./ max(abs(future_actual(:)),1e-9))*100;
metrics.R2   = 1 - sum(err.^2) / max(sum((future_actual(:)-mean(future_actual(:))).^2),eps);
metrics.time = segment_time_sec;

if show_text
    fprintf('[TCN] seg %d  RMSE=%.4f  R2=%.4f  time=%.3fs\n', seg, metrics.RMSE, metrics.R2, metrics.time);
end

% ------------------------ Plots: FULL + ZOOMED -----------------------------
if show_plots
    fig = figure(222); clf(fig); set(fig,'Name',sprintf('TCN Segment %d', seg));

    % Full view
    subplot(2,1,1);
    plot(t, wind_data, 'k-', 'LineWidth', 1); hold on;
    xline(t(idx_start), 'b--'); xline(t(idx_end), 'b--');
    xline(t(idx_end+1), 'r--'); xline(t(future_end), 'r--');
    plot(t(idx_end+1:future_end), future_pred, 'm.-', 'LineWidth', 1.2);
    title(sprintf('Full Wind View — Segment %d', seg));
    legend('Wind Data','Train Start','Train End','Pred Start','Pred End','TCN Prediction');
    xlim([t(max(1,idx_start-5)), t(min(length(t), future_end+5))]);
    grid on; grid minor;

    % Zoomed view
    subplot(2,1,2);
    plot(t(idx_end+1:future_end), future_actual, 'k-', 'LineWidth', 1.5); hold on;
    plot(t(idx_end+1:future_end), future_pred,  'm--', 'LineWidth', 1.5);
    title('Zoomed-In Segment View');
    xlabel('Time (s)'); ylabel('Wind Speed (m/s)');
    legend('True Wind','TCN Prediction');
    grid on; grid minor;
end
end
