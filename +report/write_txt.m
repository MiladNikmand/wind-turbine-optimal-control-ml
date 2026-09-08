function write_txt(report_txt_path, run_name, history_sec, predict_sec, ...
    stride_sec, selection_criterion, models_to_run, summary_table, best_model_per_segment, ...
    best_score, num_segments)
% =========================================================================
%REPORT.WRITE_TXT  (was write_prediction_report)
% Everything that used to only print to the console gets its own durable
% copy here, in plain text: settings, the ranked summary table, and the
% best-model-per-segment breakdown, plus a map of what's in the rest of
% the run folder.
% =========================================================================

    fid = fopen(report_txt_path, 'w');
    if fid == -1
        warning('write_prediction_report:cannotOpen', 'Could not open %s for writing.', report_txt_path);
        return;
    end

    fprintf(fid, '=====================================================\n');
    fprintf(fid, ' WIND PREDICTION SUITE -- RUN REPORT\n');
    fprintf(fid, ' Run: %s\n', run_name);
    fprintf(fid, '=====================================================\n\n');

    fprintf(fid, 'Settings\n');
    fprintf(fid, '--------\n');
    fprintf(fid, ' history_sec          : %.3f\n', history_sec);
    fprintf(fid, ' predict_sec          : %.3f\n', predict_sec);
    fprintf(fid, ' stride_sec           : %.3f\n', stride_sec);
    fprintf(fid, ' selection_criterion  : %s\n', selection_criterion);
    fprintf(fid, ' num_segments         : %d\n', num_segments);
    fprintf(fid, ' models_to_run        : %s\n', strjoin(models_to_run, ', '));
    fprintf(fid, '\n');

    fprintf(fid, 'Summary (sorted by RMSE)\n');
    fprintf(fid, '------------------------\n');
    fprintf(fid, '%-12s %10s %10s %10s %10s %12s\n', 'Model','RMSE','MAE','MAPE','R2','Time_sec');
    for i = 1:height(summary_table)
        fprintf(fid, '%-12s %10.5g %10.5g %10.5g %10.4g %12.5g\n', ...
            summary_table.Model{i}, summary_table.RMSE(i), summary_table.MAE(i), ...
            summary_table.MAPE(i), summary_table.R2(i), summary_table.Time_sec(i));
    end
    fprintf(fid, '\n');

    fprintf(fid, 'Best model per segment (criterion: %s)\n', selection_criterion);
    fprintf(fid, '----------------------------------------\n');
    for seg = 1:num_segments
        fprintf(fid, ' segment %3d -> %-10s (%s = %.5g)\n', ...
            seg, best_model_per_segment{seg}, selection_criterion, best_score(seg));
    end
    fprintf(fid, '\n');

    fprintf(fid, 'Output layout (relative to this run folder)\n');
    fprintf(fid, '--------------------------------------------\n');
    fprintf(fid, ' figures/<model>/seg_NNN_general.png             full-series view, per segment\n');
    fprintf(fid, ' figures/<model>/seg_NNN_zoom_actual_vs_pred.png  zoomed actual-vs-predicted\n');
    fprintf(fid, ' figures/<model>/seg_NNN_zoom_residual.png        zoomed residual (error) plot\n');
    fprintf(fid, ' figures/<model>/seg_NNN_zoom_parity.png          zoomed predicted-vs-actual scatter\n');
    fprintf(fid, ' statistics/stats_01..10_*.png                    cross-model comparison figures\n');
    fprintf(fid, ' animation/all_models_full_view.gif               all models, fixed full-series view\n');
    fprintf(fid, ' animation/all_models_zoomed.gif                  all models, camera scrolls with window\n');
    fprintf(fid, ' report/prediction_suite_report.txt               this file\n');
    fprintf(fid, ' prediction_suite_results.mat                     all raw numbers behind this report\n');

    fclose(fid);
end
