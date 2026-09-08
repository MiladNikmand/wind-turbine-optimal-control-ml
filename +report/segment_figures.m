function segment_figures(figures_root, model_name, seg, wind_data, step_size, ...
    window_size, segment_length, segment_stride, future_pred, zoom_types)
% =========================================================================
%REPORT.SEGMENT_FIGURES  (was save_segment_figures)
% Every (model, segment) prediction gets its own figures on disk:
%   - one general/full-series view
%   - a handful of zoomed-in views for detail (actual-vs-predicted,
%     residual/error, parity plot)
% Each is saved as an individual PNG (no subplots) so every figure gets
% its own file to be viewed or shared on its own. Figures render off-
% screen (Visible='off') and are closed immediately after saving, since
% a full run can generate hundreds of these.
% =========================================================================

    if isempty(future_pred) || any(~isfinite(future_pred))
        return;   % nothing meaningful to plot for a failed/skipped segment
    end

    model_dir = fullfile(figures_root, model_name);
    if ~exist(model_dir, 'dir')
        mkdir(model_dir);
    end

    N = length(wind_data);
    t = (0:N-1) * step_size;

    idx_start    = (seg-1)*segment_stride + 1;
    idx_end      = idx_start + window_size - 1;
    future_start = idx_end + 1;
    future_end   = idx_end + segment_length;

    if future_end > N
        return;
    end

    future_actual = wind_data(future_start:future_end);
    t_future      = t(future_start:future_end);
    pred          = future_pred(:);
    actual        = future_actual(:);

    tag = sprintf('seg_%03d', seg);

    % ---------------- General view: full series with windows highlighted ----------------
    fig = figure('Visible','off','Position',[100 100 900 400]);
    plot(t, wind_data, 'k-', 'LineWidth', 1); hold on;
    xline(t(idx_start), 'b--'); xline(t(idx_end), 'b--');
    xline(t(future_start), 'r--'); xline(t(future_end), 'r--');
    plot(t_future, pred, 'm.-', 'LineWidth', 1.3, 'MarkerSize', 10);
    title(sprintf('%s -- Segment %d -- Full View', upper(model_name), seg), 'Interpreter','none');
    xlabel('Time (s)'); ylabel('Wind speed (m/s)');
    legend('Wind data','Train start','Train end','Pred start','Pred end','Prediction', ...
        'Location','best');
    xlim([t(max(1,idx_start-5)) t(min(N, future_end+5))]);
    grid on; grid minor;
    saveas(fig, fullfile(model_dir, [tag '_general.png']));
    close(fig);

    % ---------------- Zoom: actual vs predicted over the horizon ----------------
    if ismember('actual_vs_pred', zoom_types)
        fig = figure('Visible','off','Position',[100 100 700 400]);
        plot(t_future, actual, 'k-', 'LineWidth', 1.6); hold on;
        plot(t_future, pred, 'm--o', 'LineWidth', 1.4, 'MarkerSize', 4);
        title(sprintf('%s -- Segment %d -- Zoom: Actual vs Predicted', upper(model_name), seg), ...
            'Interpreter','none');
        xlabel('Time (s)'); ylabel('Wind speed (m/s)');
        legend('True wind','Prediction','Location','best');
        grid on; grid minor;
        saveas(fig, fullfile(model_dir, [tag '_zoom_actual_vs_pred.png']));
        close(fig);
    end

    % ---------------- Zoom: residual (error) over the horizon ----------------
    if ismember('residual', zoom_types)
        err = actual - pred;
        fig = figure('Visible','off','Position',[100 100 700 350]);
        stem(t_future, err, 'filled', 'Color',[0.85 0.1 0.1]);
        yline(0, 'k-');
        title(sprintf('%s -- Segment %d -- Zoom: Residual (Actual - Predicted)', upper(model_name), seg), ...
            'Interpreter','none');
        xlabel('Time (s)'); ylabel('Error (m/s)');
        grid on; grid minor;
        saveas(fig, fullfile(model_dir, [tag '_zoom_residual.png']));
        close(fig);
    end

    % ---------------- Zoom: parity plot (predicted vs actual) ----------------
    if ismember('parity', zoom_types)
        fig = figure('Visible','off','Position',[100 100 450 450]);
        scatter(actual, pred, 40, 'filled'); hold on;
        lims = [min([actual; pred]), max([actual; pred])];
        if diff(lims) == 0
            lims = lims + [-1 1];
        end
        plot(lims, lims, 'k--', 'LineWidth', 1);
        axis equal; xlim(lims); ylim(lims);
        title(sprintf('%s -- Segment %d -- Zoom: Parity Plot', upper(model_name), seg), ...
            'Interpreter','none');
        xlabel('Actual (m/s)'); ylabel('Predicted (m/s)');
        grid on; grid minor;
        saveas(fig, fullfile(model_dir, [tag '_zoom_parity.png']));
        close(fig);
    end
end
