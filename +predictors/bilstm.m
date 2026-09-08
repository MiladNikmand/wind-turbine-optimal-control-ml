function [future_pred, metrics, segment_time_sec] = bilstm( ...
    wind_data, step_size, seg, history_sec, predict_sec, stride_sec, ...
    show_plots, show_text)
%PREDICTORS.BILSTM  Segment predictor -- ported verbatim from predict_with_bilstm_segment
%
% Contract (identical for every predictor in this package):
%   [future_pred, metrics, segment_time_sec] = predictors.bilstm( ...
%       wind_data, step_size, seg, history_sec, predict_sec, stride_sec, ...
%       show_plots, show_text)
%
% Body is unchanged from the thesis code. Only the function name and, where
% noted, the RBF helper calls have moved.
% =========================================================================
% Bi-LSTM Wind Predictor — per-segment, RBF-compatible I/O
% - Builds a small supervised set from a longer rolling fit window
% - Trains a compact BiLSTM encoder -> FC head to predict the whole horizon
% - Scrubs NaNs/Inf from inputs & responses; uses numeric matrix for Y
% - Robust fallbacks when data is insufficient or training fails
% Requires: Deep Learning Toolbox.
% =========================================================================

% ----- indices (like RBF) -----
dt               = step_size;
window_size      = round(history_sec  / dt);   % encoder input length
segment_length   = round(predict_sec  / dt);   % forecast horizon
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

future_actual = wind_data(idx_end+1:future_end);

% ----- longer fit window to assemble training samples -----
fit_window = max(8*window_size, 8*segment_length);
fit_start  = max(1, idx_end - fit_window + 1);
y_fit_full = wind_data(fit_start:idx_end);
Tfit       = numel(y_fit_full);

% Need at least window_size + segment_length
if Tfit < window_size + segment_length + 5
    % fallback to naive persistence if not enough data
    future_pred = repmat(wind_data(idx_end), segment_length, 1);
    segment_time_sec = 0;
    metrics = local_metrics(future_actual, future_pred, 0);
    if show_text
        fprintf('[BiLSTM] seg %d: too little data, returning persistence.\n', seg);
    end
    local_plots(); return;
end

% ----- normalization (z-score on fit window) -----
mu  = mean(y_fit_full, 'omitnan');
sig = std(y_fit_full, 'omitnan') + 1e-9;
y_norm = (y_fit_full - mu) / sig;

% sanitize any weirdness from upstream
y_norm(~isfinite(y_norm)) = 0;

% ----- build supervised samples from fit window -----
% inputs: length window_size,  targets: row vector length segment_length
Xseq = {};        % cell: each is 1×window_size sequence
Ymat = [];        % numeric matrix: numObs × segment_length

last_start = Tfit - (window_size + segment_length) + 1;
for i = 1:last_start
    xin  = y_norm(i : i+window_size-1);
    yout = y_norm(i+window_size : i+window_size+segment_length-1);

    if any(~isfinite(xin)) || any(~isfinite(yout))
        continue; % skip dirty sample
    end
    Xseq{end+1} = xin(:)';          %#ok<AGROW>
    Ymat = [Ymat; yout(:)'];        %#ok<AGROW>
end

numObs = numel(Xseq);
if numObs < 8
    % fallback if too few clean samples
    future_pred = repmat(wind_data(idx_end), segment_length, 1);
    segment_time_sec = 0;
    metrics = local_metrics(future_actual, future_pred, 0);
    if show_text
        fprintf('[BiLSTM] seg %d: too few clean samples (%d). Using persistence.\n', seg, numObs);
    end
    local_plots(); return;
end

% final scrub (belt & suspenders)
good = all(isfinite(Ymat), 2);
Xseq = Xseq(good);
Ymat = Ymat(good, :);
numObs = numel(Xseq);

% ----- network definition -----
numHidden   = 256;   % power knob: 96–192
drop        = 0.1;
miniBatch   = min(64, max(8, floor(numObs/4)));
maxEpochs   = 25;
learnRate   = 2e-3;

if isempty(ver('nnet')) && isempty(ver('Deep Learning Toolbox'))
    error('BiLSTM requires Deep Learning Toolbox. If unavailable, consider GRU/TCN/ESN/RBF.');
end

layers = [ ...
    sequenceInputLayer(1, "Name","input")
    bilstmLayer(numHidden,"OutputMode","last","Name","bilstm")
    dropoutLayer(drop,"Name","drop")
    fullyConnectedLayer(segment_length,"Name","fc")
    regressionLayer("Name","reg")];

opts = trainingOptions('adam', ...
    'MaxEpochs', maxEpochs, ...
    'MiniBatchSize', miniBatch, ...
    'InitialLearnRate', learnRate, ...
    'GradientThreshold', 1.0, ...
    'Shuffle','every-epoch', ...
    'Verbose', false);

% ----- train with robust fallback -----
try
    tic;
    net = trainNetwork(Xseq, Ymat, layers, opts);
    segment_time_sec = toc;
catch ME
    warning('[BiLSTM] seg %d: trainNetwork failed (%s). Falling back to persistence.', seg, ME.message);
    future_pred = repmat(wind_data(idx_end), segment_length, 1);
    segment_time_sec = 0;
    metrics = local_metrics(future_actual, future_pred, 0);
    local_plots(); return;
end

% ----- predict for current segment -----
past_segment = wind_data(idx_start:idx_end);
xin = (past_segment - mu) / sig;
xin(~isfinite(xin)) = 0;

yhat_norm = predict(net, {xin(:)'});
yhat_norm = yhat_norm(:);
future_pred = yhat_norm * sig + mu;

% ----- metrics + plots -----
metrics = local_metrics(future_actual, future_pred, segment_time_sec);
if show_text
    fprintf('[BiLSTM] seg %d  RMSE=%.4f  R2=%.4f  time=%.2fs  (train %d samples)\n', ...
        seg, metrics.RMSE, metrics.R2, metrics.time, numObs);
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
        fig = figure(272); clf(fig); set(fig,'Name',sprintf('BiLSTM Segment %d', seg));
        % full view
        subplot(2,1,1);
        plot(t, wind_data, 'k-', 'LineWidth', 1); hold on;
        xline(t(idx_start), 'b--'); xline(t(idx_end), 'b--');
        xline(t(idx_end+1), 'r--'); xline(t(future_end), 'r--');
        plot(t(idx_end+1:future_end), future_pred, 'm.-', 'LineWidth', 1.2);
        title(sprintf('Full Wind View — Segment %d (BiLSTM)', seg));
        legend('Wind Data','Train Start','Train End','Pred Start','Pred End','BiLSTM Prediction','Location', 'west');
        xlim([t(max(1,idx_start-5)), t(min(length(t),future_end+5))]);
        grid on; grid minor;

        % zoomed
        subplot(2,1,2);
        plot(t(idx_end+1:future_end), future_actual, 'k-', 'LineWidth', 1.5); hold on;
        plot(t(idx_end+1:future_end), future_pred,  'm--', 'LineWidth', 1.5);
        title('Zoomed-In Segment View');
        xlabel('Time (s)'); ylabel('Wind Speed (m/s)');
        legend('True Wind','BiLSTM Prediction','Location', 'west');
        grid on; grid minor;
    end
end
