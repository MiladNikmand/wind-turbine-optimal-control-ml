function [future_pred, metrics, segment_time_sec] = lstm( ...
    wind_data, step_size, seg, history_sec, predict_sec, stride_sec, ...
    show_plots, show_text)
%PREDICTORS.LSTM  Segment predictor -- ported verbatim from predict_with_lstm_segment
%
% Contract (identical for every predictor in this package):
%   [future_pred, metrics, segment_time_sec] = predictors.lstm( ...
%       wind_data, step_size, seg, history_sec, predict_sec, stride_sec, ...
%       show_plots, show_text)
%
% Body is unchanged from the thesis code. Only the function name and, where
% noted, the RBF helper calls have moved.
% =========================================================================
% LSTM (fixed random) + ridge readout, improved
% - z-score on longer fit window
% - multi-tap readout features: [h_t,c_t,h_{t-1},c_{t-1},h_{t-2},c_{t-2}]
% - residualization vs local moving-average baseline
% - ridge lambda chosen by tiny validation on tail of window
% - AR rollout with clipping and optional light smoothing
% =========================================================================

% ---------- indices ----------
dt               = step_size;
T_win            = round(history_sec  / dt);
H                 = round(predict_sec  / dt);
stride            = round(stride_sec   / dt);

t = (0:length(wind_data)-1)*dt;
idx_start  = (seg - 1)*stride + 1;
idx_end    = idx_start + T_win - 1;
future_end = idx_end   + H;

if future_end > length(wind_data)
    future_pred = nan(H,1);
    metrics = struct('RMSE',nan,'MAE',nan,'MAPE',nan,'R2',nan,'time',0);
    segment_time_sec = 0; return;
end

past_segment  = wind_data(idx_start:idx_end);
future_actual = wind_data(idx_end+1:future_end);

% ---------- longer fit window for normalization ----------
fit_window = max(3*T_win, 5*H);
fit_start  = max(1, idx_end - fit_window + 1);
y_fit_full = wind_data(fit_start:idx_end);

% z-score using fit window
mu  = mean(y_fit_full);
sig = std(y_fit_full) + 1e-9;

x_seq   = (past_segment(:)  - mu) / sig;       % T×1
y_true  = (future_actual(:) - mu) / sig;       % normalized only for metrics calc later

% ---------- local baseline (moving average) ----------
ma_len = max(3, round(0.1*T_win));            % 5% of window, at least 3
base_seq = movmean(x_seq, ma_len);

% the model predicts RESIDUALS wrt baseline: r_t = x_t - base_t
r_seq = x_seq - base_seq;

% ---------- LSTM hyperparams ----------
hidden_size  = 512;            % more capacity
clip_to_unit = false;          % we’re in z-score space; no [0,1] clamp
smooth_out   = true;           % final light smoothing

% ---------- persistent random LSTM (better init) ----------
persistent Wi Wf Wo Wc Ui Uf Uo Uc bi bf bo bc init_done
if isempty(init_done)
    rng(42);
    s = 0.6; r = 1.0;                      % input / recurrent scales
    Wi = s*randn(hidden_size,1);  Ui = r*orth(randn(hidden_size));
    Wf = s*randn(hidden_size,1);  Uf = r*orth(randn(hidden_size));
    Wo = s*randn(hidden_size,1);  Uo = r*orth(randn(hidden_size));
    Wc = s*randn(hidden_size,1);  Uc = r*orth(randn(hidden_size));
    bi = zeros(hidden_size,1);
    bf = ones(hidden_size,1)*2.0;          % strong memory
    bo = zeros(hidden_size,1);
    bc = zeros(hidden_size,1);
    init_done = true;
end
sigmoid = @(x) 1./(1+exp(-x));

% ---------- forward over window; collect multi-tap features ----------
T = numel(r_seq);
h = zeros(hidden_size,1);
c = zeros(hidden_size,1);

% store h,c history to build taps
Hstore = zeros(hidden_size, T);
Cstore = zeros(hidden_size, T);

for k = 1:T
    xk = r_seq(k);                 % residual input
    i = sigmoid(Wi*xk + Ui*h + bi);
    f = sigmoid(Wf*xk + Uf*h + bf);
    o = sigmoid(Wo*xk + Uo*h + bo);
    g = tanh(   Wc*xk + Uc*h + bc);
    c = f.*c + i.*g;
    h = o.*tanh(c);
    Hstore(:,k) = h;  Cstore(:,k) = c;
end

% build multi-tap features at t to predict residual at t+1
% taps: t, t-1, t-2  (use zeros for missing at start)
make_feat = @(k) [ ...
    Hstore(:,k); Cstore(:,k); ...
    (k-1>=1)*Hstore(:,max(k-1,1)); (k-1>=1)*Cstore(:,max(k-1,1)); ...
    (k-2>=1)*Hstore(:,max(k-2,1)); (k-2>=1)*Cstore(:,max(k-2,1)) ];

Phi = zeros(T-1, 6*hidden_size + 1);     % +1 bias
Y    = zeros(T-1, 1);                    % target: residual delta (r_{t+1})
for k = 1:(T-1)
    feat = make_feat(k);
    Phi(k,1:end-1) = feat.';
    Phi(k,end) = 1;                      % bias
    Y(k) = r_seq(k+1);                   % next residual
end

% ---------- choose ridge λ via tiny validation tail ----------
val_len = max( min( round(0.2*(T-1)),  max(5, round(0.5*H)) ), 5 ); % ~20% or ~0.5*H
tr_len  = (T-1) - val_len;
Phi_tr = Phi(1:tr_len,:);  Y_tr = Y(1:tr_len);
Phi_v  = Phi(tr_len+1:end,:); Y_v = Y(tr_len+1:end);

lambdas = [1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1];
bestL = lambdas(1); bestLoss = inf; bestTheta = [];
for L = lambdas
    I = eye(size(Phi_tr,2)); I(end,end)=0;            % don’t penalize bias
    theta = (Phi_tr.'*Phi_tr + L*I) \ (Phi_tr.'*Y_tr);
    pred_v = Phi_v*theta;
    loss = mean((pred_v - Y_v).^2);
    if loss < bestLoss
        bestLoss = loss; bestL = L; bestTheta = theta;
    end
end
theta = bestTheta;
segment_time_sec = 0; % training is tiny

% ---------- AR rollout on residuals over horizon ----------
% continue LSTM from last state, and roll r_{T+1..T+H}
r_last = r_seq(end);
h = Hstore(:,end);
c = Cstore(:,end);

r_future = zeros(H,1);
for k = 1:H
    % features at current step (t = end)
    feat_t = [h; c; Hstore(:,max(T-1,1)); Cstore(:,max(T-1,1)); Hstore(:,max(T-2,1)); Cstore(:,max(T-2,1))];
    % Note: reusing last taps; for even better fidelity, you can maintain a tap buffer of generated states
    phi_t  = [feat_t.' 1];
    r_next = phi_t*theta;
    if clip_to_unit
        r_next = max(min(r_next, 3), -3); % in z-score space; keep sane range
    end
    r_future(k) = r_next;

    % feed the predicted residual back into LSTM
    xk = r_next;
    i = sigmoid(Wi*xk + Ui*h + bi);
    f = sigmoid(Wf*xk + Uf*h + bf);
    o = sigmoid(Wo*xk + Uo*h + bo);
    g = tanh(   Wc*xk + Uc*h + bc);
    c = f.*c + i.*g;
    h = o.*tanh(c);

    % update tap memory (slide)
    Hstore = [Hstore(:,2:end) h];
    Cstore = [Cstore(:,2:end) c];
end

% reconstruct x = baseline + residual, then denormalize
% baseline for future: extend with last baseline value (persistence)
base_last = base_seq(end);
base_future = repmat(base_last, H, 1);

x_future_norm = base_future + r_future;
future_pred = x_future_norm * sig + mu;

% optional light smoothing
if smooth_out
    future_pred = smoothdata(future_pred, 'movmean', 3);
end

% ---------- metrics ----------
err = future_actual(:) - future_pred(:);
metrics.RMSE = sqrt(mean(err.^2));
metrics.MAE  = mean(abs(err));
metrics.MAPE = mean(abs(err) ./ max(abs(future_actual(:)),1e-9))*100;
metrics.R2   = 1 - sum(err.^2) / max(sum((future_actual(:)-mean(future_actual(:))).^2), eps);
metrics.time = segment_time_sec;

if show_text
    fprintf('[LSTM-imp] seg %d  RMSE=%.4f  R2=%.4f  (λ=%.0e, taps=3, H=%d)\n', ...
        seg, metrics.RMSE, metrics.R2, bestL, H);
end

% ---------- plots ----------
if show_plots
    fig = figure(262); clf(fig); set(fig,'Name',sprintf('LSTM Segment %d', seg));

    % Full view
    subplot(2,1,1);
    plot(t, wind_data, 'k-', 'LineWidth', 1); hold on;
    xline(t(idx_start), 'b--'); xline(t(idx_end), 'b--');
    xline(t(idx_end+1), 'r--'); xline(t(future_end), 'r--');
    plot(t(idx_end+1:future_end), future_pred, 'm.-', 'LineWidth', 1.2);
    title(sprintf('Full Wind View — Segment %d (LSTM-improved)', seg));
    legend('Wind Data','Train Start','Train End','Pred Start','Pred End','LSTM Prediction','Location', 'west');
    xlim([t(max(1,idx_start-5)) , t(min(length(t),future_end+5))]);
    grid on; grid minor;

    % Zoomed view
    subplot(2,1,2);
    plot(t(idx_end+1:future_end), future_actual, 'k-', 'LineWidth', 1.5); hold on;
    plot(t(idx_end+1:future_end), future_pred,  'm--', 'LineWidth', 1.5);
    title('Zoomed-In Segment View');
    xlabel('Time (s)'); ylabel('Wind Speed (m/s)');
    legend('True Wind','LSTM Prediction','Location', 'west');
    grid on; grid minor;
end
end
