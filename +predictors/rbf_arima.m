function [future_pred, metrics, segment_time_sec] = rbf_arima( ...
    wind_data, step_size, seg, history_sec, predict_sec, stride_sec, ...
    show_plots, show_text)
%PREDICTORS.RBF_ARIMA  Segment predictor -- ported verbatim from predict_with_rbf_arima_segment
%
% Contract (identical for every predictor in this package):
%   [future_pred, metrics, segment_time_sec] = predictors.rbf_arima( ...
%       wind_data, step_size, seg, history_sec, predict_sec, stride_sec, ...
%       show_plots, show_text)
%
% Body is unchanged from the thesis code. Only the function name and, where
% noted, the RBF helper calls have moved.
% =========================================================================
% RBF + ARIMA(residual) — Segment predictor (RBF-compatible I/O)
% - Train a global RBF (fast) to get base forecasts
% - On a longer rolling fit window, compute 1-step residuals of RBF
% - Fit a small ARIMA to residuals, forecast residuals over horizon
% - Final = RBF_horizon + ARIMA(residual)_horizon
% - Full + zoomed plots; safe fallbacks
% Requires Econometrics Toolbox for ARIMA.
% =========================================================================

% ------------------------ indexing (as RBF) --------------------------------
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
    metrics = struct('RMSE',nan,'MAE',nan,'MAPE',nan,'R2',nan,'time',0);
    segment_time_sec = 0; return;
end

past_segment   = wind_data(idx_start:idx_end);
future_actual  = wind_data(idx_end+1:future_end);

% ------------------------ RBF settings (same flavor as yours) --------------
num_centers    = 150;
sigma_factor   = 0.7;
use_lsqminnorm = true;

% ------------------------ Train one RBF model (global) ---------------------
rbf_model = predictors.rbf_core('train', wind_data, window_size, num_centers, sigma_factor, use_lsqminnorm);

% ------------------------ RBF horizon forecast (base) ----------------------
base_pred = predictors.rbf_core('predict', past_segment, rbf_model, segment_length);
base_pred = base_pred(:);

% ------------------------ Build residuals on a longer fit window -----------
fit_window = max(3*window_size, 5*segment_length);           % more data for ARIMA
fit_start  = max(1, idx_end - fit_window + 1);
y_fit      = wind_data(fit_start:idx_end);
Tfit       = numel(y_fit);

% 1-step-ahead base preds over the fit window (teacher forcing)
% We need past windows of length window_size to predict the next point.
i0 = window_size;                      % last index of first past window within y_fit
if Tfit <= i0
    % too short -> skip ARIMA and return base_pred only
    future_pred = base_pred;
    segment_time_sec = 0;
    [metrics,~] = compute_metrics(future_actual, future_pred);
    if show_text
        fprintf('[RBF+ARIMA] seg %d: short fit window, using pure RBF. RMSE=%.4f R2=%.4f\n', ...
            seg, metrics.RMSE, metrics.R2);
    end
    do_plots(); return;
end

base_fit_pred = zeros(Tfit - i0, 1);
for k = i0:(Tfit-1)
    past = y_fit(k - window_size + 1 : k);           % length window_size
    yhat = predictors.rbf_core('predict', past, rbf_model, 1);  % 1-step
    base_fit_pred(k - i0 + 1) = yhat;
end
y_fit_targets = y_fit(i0+1:end);                     % align to predicted next
resid_fit     = y_fit_targets(:) - base_fit_pred(:); % residuals for ARIMA

% ------------------------ Fit small ARIMA on residuals ---------------------
segment_tic = tic;
bestMdl = [];
bestAIC = inf;
cands = [];
for p=0:2, for d=0:1, for q=0:2
    if p==0 && q==0, continue; end
    cands = [cands; p d q]; %#ok<AGROW>
end, end, end

for i = 1:size(cands,1)
    p = cands(i,1); d = cands(i,2); q = cands(i,3);
    try
        mdl = arima('Constant',NaN,'ARLags',1:p,'D',d,'MALags',1:q);
        est = estimate(mdl, resid_fit, 'Display','off', 'TolCon',1e-6, 'TolFun',1e-6);
        S   = summarize(est);
        logL = S.LogLikelihood;
        k    = numel(est.AR) + numel(est.MA) + ~isempty(est.Constant) + 1; % +1 variance
        AIC  = 2*k - 2*logL;
        if AIC < bestAIC
            bestAIC = AIC; bestMdl = est;
        end
    catch
        % skip unstable fits
    end
end

% Forecast residuals over the horizon
if isempty(bestMdl)
    resid_fore = zeros(segment_length,1);     % fallback: no residual correction
else
    try
        resid_fore = forecast(bestMdl, segment_length, resid_fit);
    catch
        resid_fore = zeros(segment_length,1);
    end
end

% ------------------------ Final forecast = base + residual -----------------
future_pred = base_pred + resid_fore(:);
segment_time_sec = toc(segment_tic);

% ------------------------ Metrics + optional plots -------------------------
[metrics,err] = compute_metrics(future_actual, future_pred);
if show_text
    fprintf('[RBF+ARIMA] seg %d  RMSE=%.4f  R2=%.4f  time=%.3fs\n', ...
        seg, metrics.RMSE, metrics.R2, segment_time_sec);
end
do_plots();

% ====================== helpers ======================
    function [m,err_] = compute_metrics(y_true, y_hat)
        y_true = y_true(:); y_hat = y_hat(:);
        err_   = y_true - y_hat;
        m.RMSE = sqrt(mean(err_.^2));
        m.MAE  = mean(abs(err_));
        m.MAPE = mean(abs(err_) ./ max(abs(y_true),1e-9)) * 100;
        m.R2   = 1 - sum(err_.^2) / max(sum((y_true - mean(y_true)).^2), eps);
        m.time = segment_time_sec;
    end

    function do_plots()
        if ~show_plots, return; end
        fig = figure(241); clf(fig); set(fig,'Name',sprintf('RBF+ARIMA Segment %d', seg));

        % Full view
        subplot(2,1,1);
        plot(t, wind_data, 'k-', 'LineWidth', 1); hold on;
        xline(t(idx_start), 'b--'); xline(t(idx_end), 'b--');
        xline(t(idx_end+1), 'r--'); xline(t(future_end), 'r--');
        plot(t(idx_end+1:future_end), future_pred, 'm.-', 'LineWidth', 1.2);
        title(sprintf('Full Wind View — Segment %d (RBF+ARIMA)', seg));
        legend('Wind Data','Train Start','Train End','Pred Start','Pred End','Hybrid Prediction');
        xlim([t(max(1,idx_start-5)) , t(min(length(t),future_end+5))]);
        grid on; grid minor;

        % Zoomed view
        subplot(2,1,2);
        plot(t(idx_end+1:future_end), future_actual, 'k-', 'LineWidth', 1.5); hold on;
        plot(t(idx_end+1:future_end), future_pred,  'm--', 'LineWidth', 1.5);
        title('Zoomed-In Segment View');
        xlabel('Time (s)'); ylabel('Wind Speed (m/s)');
        legend('True Wind','Hybrid Prediction');
        grid on; grid minor;
    end
end
