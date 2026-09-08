function animation(gif_path, wind_data, step_size, window_size, ...
    segment_length, segment_stride, predictions, models_to_run, view_mode, best_idx, fps)
% =========================================================================
%REPORT.ANIMATION  (was save_prediction_animation)
% One frame per segment. Each frame shows the true wind trace, the
% current training window + forecast horizon markers, and every model's
% prediction for that segment overlaid in its own color -- so you can
% watch, frame by frame, how each ANN's forecast compares to the ground
% truth and to every other model as the window slides forward. The
% segment's winning model (per best_idx, if supplied) is drawn thicker
% so you can also watch the "winner" change over time.
%
%   view_mode : 'full'   -- fixed x-limits spanning the whole series;
%                           watch the forecast segment travel across a
%                           wide view of the whole wind trace.
%               'zoomed' -- x-limits scroll to frame the current window
%                           + forecast horizon each frame; better for
%                           actually seeing model-to-model differences.
%   best_idx  : optional 1xnum_segments vector (index into models_to_run)
%               of the winning model per segment. Pass [] to skip the
%               highlight.
%   fps       : optional frames per second for the GIF (default 2).
% =========================================================================

    if nargin < 10
        best_idx = [];
    end
    if nargin < 11 || isempty(fps)
        fps = 2;
    end

    num_models   = numel(models_to_run);
    num_segments = size(predictions, 2);
    N = length(wind_data);
    t = (0:N-1) * step_size;

    colors = lines(num_models);
    delay  = 1 / fps;
    wrote_first_frame = false;

    for seg = 1:num_segments
        idx_start    = (seg-1)*segment_stride + 1;
        idx_end      = idx_start + window_size - 1;
        future_start = idx_end + 1;
        future_end   = idx_end + segment_length;

        if future_end > N
            continue;
        end

        fig = figure('Visible','off','Position',[100 100 950 450]);
        plot(t, wind_data, 'k-', 'LineWidth', 1.1); hold on;
        xline(t(idx_start),    'b--', 'HandleVisibility','off');
        xline(t(idx_end),      'b--', 'HandleVisibility','off');
        xline(t(future_start), 'r--', 'HandleVisibility','off');
        xline(t(future_end),   'r--', 'HandleVisibility','off');

        t_future = t(future_start:future_end);
        legend_entries = {'Wind data'};

        for m = 1:num_models
            pred = predictions{m, seg};
            if isempty(pred) || any(~isfinite(pred))
                continue;
            end
            lw = 1.4;
            if ~isempty(best_idx) && best_idx(seg) == m
                lw = 3.2;   % highlight this segment's winning model
            end
            plot(t_future, pred(:), '-', 'Color', colors(m,:), 'LineWidth', lw);
            legend_entries{end+1} = models_to_run{m}; %#ok<AGROW>
        end

        legend(legend_entries, 'Location','eastoutside', 'Interpreter','none');
        xlabel('Time (s)'); ylabel('Wind speed (m/s)');
        title(sprintf('Segment %d / %d -- all models', seg, num_segments), 'Interpreter','none');
        grid on; grid minor;

        switch view_mode
            case 'zoomed'
                xlim([t(max(1,idx_start-5)) t(min(N, future_end+5))]);
            otherwise
                xlim([t(1) t(end)]);
        end

        frame = getframe(fig);
        [imind, cm] = rgb2ind(frame2im(frame), 256);
        if ~wrote_first_frame
            imwrite(imind, cm, gif_path, 'gif', 'LoopCount', Inf, 'DelayTime', delay);
            wrote_first_frame = true;
        else
            imwrite(imind, cm, gif_path, 'gif', 'WriteMode', 'append', 'DelayTime', delay);
        end
        close(fig);
    end
end
