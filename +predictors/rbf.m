function [future_pred, metrics, segment_time_sec] = rbf(wind_data, step_size, seg, ...
                                history_sec, predict_sec, stride_sec, ...
                                show_plots, show_text)
%PREDICTORS.RBF  Segment predictor -- ported verbatim from predict_with_rbf_segment
%
% Contract (identical for every predictor in this package):
%   [future_pred, metrics, segment_time_sec] = predictors.rbf( ...
%       wind_data, step_size, seg, history_sec, predict_sec, stride_sec, ...
%       show_plots, show_text)
%
% Body is unchanged from the thesis code. Only the function name and, where
% noted, the RBF helper calls have moved.
% =========================================================================
% RBF Wind Predictor — Predict One Segment at a Time (Self-contained)
% =========================================================================

% ------------------------ Hyperparameters ----------------------------------
dt               = step_size;
window_size      = round(history_sec / dt);
segment_length   = round(predict_sec / dt);
segment_duration = round(stride_sec / dt);
num_centers      = 150;
sigma_factor     = 0.7;
use_lsqminnorm   = true;

t = (0:length(wind_data)-1) * dt;

% -------------------- Segment Indexing -------------------------------------
idx_start   = (seg - 1)*segment_duration + 1;
idx_end     = idx_start + window_size - 1;
future_end  = idx_end + segment_length;

% Out-of-range guard
% if future_end > length(wind_data)
%     future_pred = nan(segment_length,1);
%     metrics = struct('RMSE', nan, 'MAE', nan, 'MAPE', nan, 'R2', nan);
%     return;
% end

if future_end > length(wind_data)
    future_pred = nan(segment_length,1);
    metrics = struct('RMSE', nan, 'MAE', nan, 'MAPE', nan, 'R2', nan);
    segment_time_sec = 0;   % <-- make sure this is set
    return;
end

past_segment   = wind_data(idx_start:idx_end);
future_actual  = wind_data(idx_end+1:future_end);




% -----------------------------------------------------------------------
segment_start_time = tic;
% ---------------------- Train RBF Model ------------------------------------
rbf_model = predictors.rbf_core('train', wind_data, window_size, num_centers, sigma_factor, use_lsqminnorm);
% ---------------------- Predict This Segment -------------------------------
future_pred = predictors.rbf_core('predict', past_segment, rbf_model, segment_length);
% -----------------------------------------------------------------------
segment_time_sec = toc(segment_start_time);
% -----------------------------------------------------------------------



% ---------------------- Metrics --------------------------------------------
err = future_actual - future_pred;
metrics.RMSE = sqrt(mean(err.^2));
metrics.MAE  = mean(abs(err));
metrics.MAPE = mean(abs(err) ./ max(abs(future_actual), 1e-9)) * 100;
metrics.R2   = 1 - sum(err.^2) / max(sum((future_actual - mean(future_actual)).^2), eps);

% ---------------------- Visualization --------------------------------------
if show_text
    fprintf('--------------------------------------------------------------- \n')
    fprintf('Segment #%d | Train: %.2f–%.2f s | Predict: %.2f–%.2f s\n', ...
        seg, t(idx_start), t(idx_end), t(idx_end+1), t(future_end));
    fprintf(' \n')
    fprintf('   RMSE = %.4f | MAE = %.4f | MAPE = %.2f%% | R² = %.4f\n', ...
        metrics.RMSE, metrics.MAE, metrics.MAPE, metrics.R2);
    fprintf('\n')
    fprintf("Segment %d took %.3f seconds to process.\n", seg, segment_time_sec);
    fprintf('--------------------------------------------------------------- \n')
end

if show_plots
    fig = figure(99); clf(fig); set(fig,'Name','RBF Segment Prediction');

    subplot(2,1,1);
    plot(t, wind_data, 'k-', 'LineWidth', 1); hold on;
    xline(t(idx_start), 'b--'); xline(t(idx_end), 'b--');
    xline(t(idx_end+1), 'r--'); xline(t(future_end), 'r--');
    plot(t(idx_end+1:future_end), future_pred, 'm.-', 'LineWidth', 1);
    title(sprintf('Full Wind View — Segment %d', seg));
    legend('Wind Data', 'Train Start', 'Train End', 'Pred Start', 'Pred End', 'Prediction');
    xlim([t(idx_start)-5, t(future_end)+5]);
    grid on; grid minor;

    subplot(2,1,2);
    plot(t(idx_end+1:future_end), future_actual, 'k-', 'LineWidth', 1.5); hold on;
    plot(t(idx_end+1:future_end), future_pred, 'm--', 'LineWidth', 1.5);
    title('Zoomed-In Segment View');
    xlabel('Time (s)'); ylabel('Wind Speed (m/s)');
    legend('True Wind','Predicted Wind');
    grid on; grid minor;
end

end
