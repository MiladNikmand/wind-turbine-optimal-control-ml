function stats_figures(stats_dir, summary_table, RMSE, MAE, MAPE, R2, TIME, ...
    models_to_run, best_idx, wind_data, step_size, predicted_wind_bestof, ...
    history_sec, predict_sec, stride_sec)
% =========================================================================
%REPORT.STATS_FIGURES  (was save_statistical_figures)
% A broader statistical comparison than a single bar + boxplot: ten
% figures, each its own PNG, each with its own "limelight" -- no
% subplots crammed together.
%
%   01  mean RMSE per model (sorted, bar)
%   02  mean MAE per model (bar)
%   03  mean R^2 per model (bar)
%   04  mean compute time per model (bar, log scale)
%   05  RMSE spread per model (boxplot)
%   06  MAE spread per model (boxplot)
%   07  RMSE heatmap, model x segment
%   08  "segments won" count per model (bar)
%   09  RMSE across segments, one line per model
%   10  best-of stitched prediction vs. true wind (kept, as requested)
% =========================================================================

    cats = categorical(models_to_run, models_to_run);   % fixes x-axis order

    % ---- 01: mean RMSE, sorted ascending ----
    fig = figure('Visible','off','Position',[100 100 700 400]);
    bar(categorical(summary_table.Model, summary_table.Model), summary_table.RMSE);
    ylabel('Mean RMSE'); title('Mean RMSE per model (sorted, lower is better)');
    grid on; grid minor;
    saveas(fig, fullfile(stats_dir, 'stats_01_mean_rmse_bar.png')); close(fig);

    % ---- 02: mean MAE ----
    fig = figure('Visible','off','Position',[100 100 700 400]);
    bar(cats, mean(MAE, 2, 'omitnan'));
    ylabel('Mean MAE'); title('Mean MAE per model');
    grid on; grid minor;
    saveas(fig, fullfile(stats_dir, 'stats_02_mean_mae_bar.png')); close(fig);

    % ---- 03: mean R^2 ----
    fig = figure('Visible','off','Position',[100 100 700 400]);
    bar(cats, mean(R2, 2, 'omitnan'));
    ylabel('Mean R^2'); title('Mean R^2 per model (closer to 1 is better)');
    grid on; grid minor;
    saveas(fig, fullfile(stats_dir, 'stats_03_mean_r2_bar.png')); close(fig);

    % ---- 04: mean compute time (log scale -- some models are orders of magnitude slower) ----
    fig = figure('Visible','off','Position',[100 100 700 400]);
    bar(cats, mean(TIME, 2, 'omitnan'));
    set(gca, 'YScale', 'log');
    ylabel('Mean time per segment (s, log scale)'); title('Mean computation time per model');
    grid on; grid minor;
    saveas(fig, fullfile(stats_dir, 'stats_04_mean_time_bar.png')); close(fig);

    % ---- 05: RMSE boxplot ----
    fig = figure('Visible','off','Position',[100 100 800 400]);
    boxplot(RMSE', 'Labels', models_to_run);
    ylabel('RMSE per segment'); title('RMSE spread across segments');
    grid on; grid minor;
    saveas(fig, fullfile(stats_dir, 'stats_05_rmse_boxplot.png')); close(fig);

    % ---- 06: MAE boxplot ----
    fig = figure('Visible','off','Position',[100 100 800 400]);
    boxplot(MAE', 'Labels', models_to_run);
    ylabel('MAE per segment'); title('MAE spread across segments');
    grid on; grid minor;
    saveas(fig, fullfile(stats_dir, 'stats_06_mae_boxplot.png')); close(fig);

    % ---- 07: RMSE heatmap, model x segment ----
    fig = figure('Visible','off','Position',[100 100 900 400]);
    imagesc(RMSE); colorbar;
    set(gca, 'YTick', 1:numel(models_to_run), 'YTickLabel', models_to_run);
    xlabel('Segment'); title('RMSE heatmap (model x segment)');
    saveas(fig, fullfile(stats_dir, 'stats_07_rmse_heatmap.png')); close(fig);

    % ---- 08: segments won per model ----
    win_counts = zeros(numel(models_to_run), 1);
    for m = 1:numel(models_to_run)
        win_counts(m) = sum(best_idx == m);
    end
    fig = figure('Visible','off','Position',[100 100 700 400]);
    bar(cats, win_counts);
    ylabel('# segments won'); title('Segment "wins" per model (best score per segment)');
    grid on; grid minor;
    saveas(fig, fullfile(stats_dir, 'stats_08_win_count_bar.png')); close(fig);

    % ---- 09: RMSE across segments, one line per model ----
    fig = figure('Visible','off','Position',[100 100 900 450]);
    hold on;
    for m = 1:numel(models_to_run)
        plot(1:size(RMSE,2), RMSE(m,:), '-o', 'DisplayName', models_to_run{m});
    end
    xlabel('Segment'); ylabel('RMSE');
    title('RMSE across segments, per model');
    legend('Location','bestoutside');
    grid on; grid minor;
    saveas(fig, fullfile(stats_dir, 'stats_09_rmse_over_segments.png')); close(fig);

    % ---- 10: best-of stitched prediction vs. true wind (kept) ----
    t = (0:length(wind_data)-1) * step_size;
    fig = figure('Visible','off','Position',[100 100 1000 400]);
    plot(t, wind_data, 'k-', 'LineWidth', 1); hold on;
    plot(t, predicted_wind_bestof, 'm--', 'LineWidth', 1.3);
    legend('True wind', 'Best-of predicted wind (per-segment winner)', 'Location', 'best');
    xlabel('Time (s)'); ylabel('Wind speed (m/s)');
    title(sprintf('Best-of stitched prediction | history=%.2fs predict=%.2fs stride=%.2fs', ...
        history_sec, predict_sec, stride_sec));
    grid on; grid minor;
    saveas(fig, fullfile(stats_dir, 'stats_10_bestof_stitched.png')); close(fig);
end
