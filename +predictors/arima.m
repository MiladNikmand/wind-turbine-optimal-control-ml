function [future_pred, metrics, segment_time_sec] = arima( ...
    wind_data, step_size, seg, history_sec, predict_sec, stride_sec, ...
    show_plots, show_text)
%PREDICTORS.ARIMA  Segment predictor -- ported verbatim from predict_with_arima_segment
%
% Contract (identical for every predictor in this package):
%   [future_pred, metrics, segment_time_sec] = predictors.arima( ...
%       wind_data, step_size, seg, history_sec, predict_sec, stride_sec, ...
%       show_plots, show_text)
%
% Body is unchanged from the thesis code. Only the function name and, where
% noted, the RBF helper calls have moved.
% =========================================================================
% ARIMA Wind Predictor — One segment at a time (RBF-compatible I/O)
% Improvements:
%   - Fit on a longer rolling window (more data -> stabler params)
%   - Box–Cox variance stabilization (lambda sweep)
%   - Broader (p,d,q) grid, AIC model selection
%   - Full + zoomed plots, safe fallbacks
% Requires Econometrics Toolbox.
% =========================================================================

% ---- indexing (same as RBF) ----
dt               = step_size;
window_size      = round(history_sec  / dt);
segment_length   = round(predict_sec  / dt);
segment_duration = round(stride_sec   / dt);

t = (0:length(wind_data)-1)*dt;
idx_start  = (seg - 1)*segment_duration + 1;
idx_end    = idx_start + window_size - 1;
future_end = idx_end   + segment_length;

% guard
if future_end > length(wind_data)
    future_pred = nan(segment_length,1);
    metrics = struct('RMSE',nan,'MAE',nan,'MAPE',nan,'R2',nan,'time',0);
    segment_time_sec = 0; return;
end

future_actual = wind_data(idx_end+1:future_end);

% ---- choose a longer fit window than the control history ----
fit_window = max(3*window_size, 5*segment_length);           % tweak if needed
fit_start  = max(1, idx_end - fit_window + 1);
y_fit_full = wind_data(fit_start:idx_end);

% ---- Box–Cox transform (lambda sweep) ----
% 0 => log transform. Pick lambda minimizing var(diff(y_bc)).
y_pos = max(y_fit_full, 1e-6);                                % safety for log
lambda_grid = [-0.5 0 0.25 0.5 1];
bestL = 1; bestVar = inf;
for L = lambda_grid
    if L == 0
        y_bc = log(y_pos);
    else
        y_bc = (y_pos.^L - 1)/L;
    end
    v = var(diff(y_bc));
    if isfinite(v) && v < bestVar
        bestVar = v; bestL = L;
    end
end
L = bestL;
if L == 0
    y_fit = log(y_pos);
else
    y_fit = (y_pos.^L - 1)/L;
end
y_fit = y_fit(:);

% ---- candidate (p,d,q) grid; small but decent ----
cands = [];
for p = 0:3
    for d = 0:2
        for q = 0:3
            if p==0 && q==0, continue; end
            cands = [cands; p d q]; %#ok<AGROW>
        end
    end
end

% ---- estimate models, pick lowest AIC ----
bestAIC = inf; bestMdl = [];
segment_tic = tic;
for i = 1:size(cands,1)
    p = cands(i,1); d = cands(i,2); q = cands(i,3);
    try
        mdl = arima('Constant',NaN,'ARLags',1:p,'D',d,'MALags',1:q);
        est = estimate(mdl, y_fit, 'Display','off', 'TolCon',1e-6, 'TolFun',1e-6);
        S = summarize(est);
        logL = S.LogLikelihood;
        % crude parameter count: AR + MA + Constant + variance + (maybe) d
        k = numel(est.AR) + numel(est.MA) + ~isempty(est.Constant) + 1; % +1 for variance
        AIC = 2*k - 2*logL;
        if AIC < bestAIC
            bestAIC = AIC; bestMdl = est;
        end
    catch
        % skip unstable/failed fits
    end
end

% fallback if all candidates failed
if isempty(bestMdl)
    try
        bestMdl = estimate(arima(1,1,1), y_fit, 'Display','off');
    catch
        % ultimate fallback: naive hold-last on original scale
        future_pred = repmat(wind_data(idx_end), segment_length, 1);
        segment_time_sec = toc(segment_tic);
        do_metrics_and_plots(); return;
    end
end

% ---- forecast in transformed space ----
try
    yhat_f = forecast(bestMdl, segment_length, y_fit);
catch
    % fallback to naive in transformed space
    yhat_f = repmat(y_fit(end), segment_length, 1);
end

% ---- invert Box–Cox ----
if L == 0
    future_pred = exp(yhat_f);
else
    future_pred = max(L*yhat_f + 1, 0).^(1/L);
end

segment_time_sec = toc(segment_tic);

% ---- metrics + plots ----
do_metrics_and_plots();

% ================= nested helper =================
    function do_metrics_and_plots()
        pred = future_pred(:);
        err  = future_actual(:) - pred;
        metrics.RMSE = sqrt(mean(err.^2));
        metrics.MAE  = mean(abs(err));
        metrics.MAPE = mean(abs(err) ./ max(abs(future_actual(:)),1e-9)) * 100;
        metrics.R2   = 1 - sum(err.^2) / max(sum((future_actual(:)-mean(future_actual(:))).^2), eps);
        metrics.time = segment_time_sec;

        if show_text
            fprintf('[ARIMA+] seg %d  RMSE=%.4f  R2=%.4f  time=%.3fs  (L=%.2f, fit_window=%d)\n', ...
                seg, metrics.RMSE, metrics.R2, metrics.time, L, numel(y_fit_full));
        end

        if show_plots
            fig = figure(231); clf(fig); set(fig,'Name',sprintf('ARIMA+ Segment %d', seg));

            % Full view
            subplot(2,1,1);
            plot(t, wind_data, 'k-', 'LineWidth', 1); hold on;
            xline(t(idx_start), 'b--'); xline(t(idx_end), 'b--');
            xline(t(idx_end+1), 'r--'); xline(t(future_end), 'r--');
            plot(t(idx_end+1:future_end), pred, 'm.-', 'LineWidth', 1.2);
            title(sprintf('Full Wind View — Segment %d (ARIMA+)', seg));
            legend('Wind Data','Train Start','Train End','Pred Start','Pred End','ARIMA+ Prediction');
            xlim([t(max(1,idx_start-5)) , t(min(length(t),future_end+5))]);
            grid on; grid minor;

            % Zoomed view
            subplot(2,1,2);
            plot(t(idx_end+1:future_end), future_actual, 'k-', 'LineWidth', 1.5); hold on;
            plot(t(idx_end+1:future_end), pred,  'm--', 'LineWidth', 1.5);
            title('Zoomed-In Segment View');
            xlabel('Time (s)'); ylabel('Wind Speed (m/s)');
            legend('True Wind','ARIMA+ Prediction');
            grid on; grid minor;
        end
    end
end
