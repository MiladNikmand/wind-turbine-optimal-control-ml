function assets = export_prediction(P, wind, step_size, cfg, R, opts)
%REPORT.EXPORT_PREDICTION  Write the per-model prediction assets and summaries.
%
%   assets = report.export_prediction(P, wind, step_size, cfg, R, opts)
%
%   P     report.paths(...) with create=true
%   wind  1xN wind vector
%   cfg   .history_sec .predict_sec .stride_sec .num_segments
%         .models_to_run .selection_criterion
%   R     results: .predictions .summary_table .best_idx .best_model_per_segment
%         .best_score .predicted_wind_bestof .RMSE .MAE .MAPE .R2 .TIME
%         .nRMSE .sMAPE
%   opts  .export_web (true) .make_gifs (true) .gif_fps (2)
%         .dpi_hires (150) .dpi_web (96) .web_width (640) .web_colors (128)
%
% Produces, per model: a full-timeline PNG and a zoomed PNG in both tiers,
% plus full and zoomed GIFs in both tiers. The hi-res tier is thesis-grade;
% the web tier is what report.html loads.

    if nargin < 6 || isempty(opts), opts = struct(); end
    opts = setdef(opts, 'export_web', true);
    opts = setdef(opts, 'make_gifs',  true);
    opts = setdef(opts, 'gif_fps',    2);
    opts = setdef(opts, 'dpi_hires',  150);
    opts = setdef(opts, 'dpi_web',    96);
    opts = setdef(opts, 'web_width',  640);
    opts = setdef(opts, 'web_colors', 128);
    opts = setdef(opts, 'gif_render_px', [1280 640]);

    wind = reshape(wind, 1, []);
    N = numel(wind);
    t = (0:N-1) * step_size;

    window_size    = round(cfg.history_sec / step_size);
    segment_length = round(cfg.predict_sec / step_size);
    segment_stride = round(cfg.stride_sec  / step_size);

    models = cfg.models_to_run;
    nm = numel(models);

    y_min = min(wind) - 0.1*std(wind);
    y_max = max(wind) + 0.1*std(wind);

    mid_seg = max(1, round(cfg.num_segments / 2));
    t_start = (mid_seg - 1) * cfg.stride_sec;
    t_end   = t_start + cfg.history_sec + cfg.predict_sec*3;

    assets = struct();

    %% ---------------- static figures ----------------
    fig = figure('Visible','off','Position',[100 100 1200 500],'GraphicsSmoothing','off');
    ax  = axes('Parent', fig);
    plot(ax, t, wind, 'k-', 'LineWidth', 1.2, 'DisplayName', 'Original Wind');
    hold(ax, 'on');
    hp = plot(ax, nan, nan, 'r--', 'LineWidth', 1.2);
    grid(ax, 'on'); xlabel(ax, 'Time (s)'); ylabel(ax, 'Wind Speed (m/s)');
    legend(ax, 'Location', 'best');

    for m = 1:nm
        name = models{m};
        stitched = predictors.stitch(R.predictions, repmat(m,1,cfg.num_segments), ...
            wind(:), step_size, cfg.history_sec, cfg.predict_sec, cfg.stride_sec);

        set(hp, 'XData', t, 'YData', stitched(:)', 'DisplayName', sprintf('%s forecast', name));
        ylim(ax, [y_min y_max]);

        xlim(ax, [t(1) t(end)]);
        title(ax, sprintf('%s -- full view', upper(name)), 'Interpreter','none');
        exportgraphics(fig, fullfile(P.figures, sprintf('%s_full_timeline.png', name)), ...
            'Resolution', opts.dpi_hires);
        if opts.export_web
            exportgraphics(fig, fullfile(P.web_figures, sprintf('%s_full.png', name)), ...
                'Resolution', opts.dpi_web);
        end

        xlim(ax, [t_start t_end]);
        title(ax, sprintf('%s -- zoomed view', upper(name)), 'Interpreter','none');
        exportgraphics(fig, fullfile(P.figures, sprintf('%s_zoomed_timeline.png', name)), ...
            'Resolution', opts.dpi_hires);
        if opts.export_web
            exportgraphics(fig, fullfile(P.web_figures, sprintf('%s_zoomed.png', name)), ...
                'Resolution', opts.dpi_web);
        end

        assets.(matlab.lang.makeValidName(name)) = ...
            asset_paths(P.profile_id, name, opts.export_web, opts.make_gifs);
        fprintf('  - %s static\n', name);
    end
    close(fig);

    %% ---------------- per-model animations ----------------
    if opts.make_gifs
        for m = 1:nm
            name = models{m};
            adir = fullfile(P.animation, name);
            if ~exist(adir,'dir'), mkdir(adir); end
            g_full = fullfile(adir, 'full.gif');
            g_zoom = fullfile(adir, 'zoomed.gif');
            if opts.export_web
                wdir = fullfile(P.web_anim, name);
                if ~exist(wdir,'dir'), mkdir(wdir); end
                gw_full = fullfile(wdir, 'full.gif');
                gw_zoom = fullfile(wdir, 'zoomed.gif');
            end

            gf = figure('Visible','off','Color','w', ...
                'Position',[100 100 opts.gif_render_px(1) opts.gif_render_px(2)]);
            ga = axes('Parent', gf);
            plot(ga, t, wind, 'Color',[0.2 0.2 0.2], 'LineWidth',1, 'DisplayName','Original Wind');
            hold(ga,'on');
            hh = plot(ga, nan, nan, 'b-', 'LineWidth',1.5, 'DisplayName','History window');
            hq = plot(ga, nan, nan, 'r-', 'LineWidth',2.0, 'DisplayName',sprintf('%s prediction', name));
            grid(ga,'on'); ylim(ga,[y_min y_max]); legend(ga,'Location','northeast');

            ref_h = []; ref_w = [];
            for seg = 1:cfg.num_segments
                i0 = (seg-1)*segment_stride + 1;
                hist_idx = i0 : i0+window_size-1;
                pred_idx = (i0+window_size) : (i0+window_size+segment_length-1);
                if pred_idx(end) > N, break; end

                yp = R.predictions{m, seg};
                if isempty(yp), yp = nan(segment_length,1); end

                set(hh, 'XData', t(hist_idx), 'YData', wind(hist_idx));
                set(hq, 'XData', t(pred_idx), 'YData', yp(:)');

                xlim(ga, [0 t(end)]);
                title(ga, sprintf('%s | seg %d/%d (full)', upper(name), seg, cfg.num_segments), ...
                    'Interpreter','none');
                cd1 = getframe(gf); cd1 = cd1.cdata;
                if isempty(ref_h), ref_h = size(cd1,1); ref_w = size(cd1,2); end
                cd1 = lock_frame(cd1, ref_h, ref_w);
                [im,cm] = rgb2ind(cd1, 256, 'nodither');
                write_gif(g_full, im, cm, seg==1, opts.gif_fps);
                if opts.export_web
                    sm = shrink(cd1, opts.web_width);
                    [iw,cw] = rgb2ind(sm, opts.web_colors, 'nodither');
                    write_gif(gw_full, iw, cw, seg==1, opts.gif_fps);
                end

                lo = t(i0);
                hi = t(min(pred_idx(end) + segment_stride*2, N));
                xlim(ga, [lo hi]);
                title(ga, sprintf('%s | seg %d/%d (zoomed)', upper(name), seg, cfg.num_segments), ...
                    'Interpreter','none');
                cd2 = lock_frame(getframe(gf).cdata, ref_h, ref_w);
                [im,cm] = rgb2ind(cd2, 256, 'nodither');
                write_gif(g_zoom, im, cm, seg==1, opts.gif_fps);
                if opts.export_web
                    sm = shrink(cd2, opts.web_width);
                    [iw,cw] = rgb2ind(sm, opts.web_colors, 'nodither');
                    write_gif(gw_zoom, iw, cw, seg==1, opts.gif_fps);
                end
            end
            close(gf);
            fprintf('  - %s gifs\n', name);
        end
    end

    %% ---------------- JSON summary ----------------
    write_summary_json(P, cfg, R, assets, wind, step_size);

    %% ---------------- CSV ----------------
    writetable(R.summary_table, fullfile(P.report, 'summary_table.csv'));

    [Mg, Sg] = ndgrid(1:nm, 1:cfg.num_segments);
    seg_tbl = table(reshape(models(Mg(:)),[],1), Sg(:), R.RMSE(:), R.nRMSE(:), ...
        R.MAE(:), R.MAPE(:), R.sMAPE(:), R.R2(:), R.TIME(:), ...
        'VariableNames', {'Model','Segment','RMSE','nRMSE','MAE','MAPE','sMAPE_pct','R2','Time_sec'});
    writetable(seg_tbl, fullfile(P.report, 'segment_metrics.csv'));

    writetable(table((1:cfg.num_segments)', R.best_model_per_segment(:), R.best_score(:), ...
        'VariableNames', {'Segment','Model', cfg.selection_criterion}), ...
        fullfile(P.report, 'best_model_per_segment.csv'));

    %% ---------------- LaTeX ----------------
    write_tex(P.tex_file, R.summary_table);

    %% ---------------- Markdown ----------------
    write_markdown(P, cfg, R, opts.export_web, wind);
end

%% ========================= helpers =========================
function o = setdef(o, f, v)
    if ~isfield(o, f) || isempty(o.(f)), o.(f) = v; end
end

function a = asset_paths(profile_id, name, web, gifs)
    base = ['wind_profiles/' profile_id '/prediction'];
    a.static_full_hires = sprintf('%s/figures/%s_full_timeline.png',   base, name);
    a.static_zoom_hires = sprintf('%s/figures/%s_zoomed_timeline.png', base, name);
    if web
        a.static_full   = sprintf('%s/web/figures/%s_full.png',   base, name);
        a.static_zoomed = sprintf('%s/web/figures/%s_zoomed.png', base, name);
    else
        a.static_full   = a.static_full_hires;
        a.static_zoomed = a.static_zoom_hires;
    end
    if gifs
        a.gif_full_hires   = sprintf('%s/animation/%s/full.gif',   base, name);
        a.gif_zoomed_hires = sprintf('%s/animation/%s/zoomed.gif', base, name);
        if web
            a.gif_full   = sprintf('%s/web/animation/%s/full.gif',   base, name);
            a.gif_zoomed = sprintf('%s/web/animation/%s/zoomed.gif', base, name);
        else
            a.gif_full   = a.gif_full_hires;
            a.gif_zoomed = a.gif_zoomed_hires;
        end
    end
end

function write_gif(path, im, cm, first, fps)
    if first
        imwrite(im, cm, path, 'gif', 'Loopcount', inf, 'DelayTime', 1/fps);
    else
        imwrite(im, cm, path, 'gif', 'WriteMode', 'append', 'DelayTime', 1/fps);
    end
end

function out = lock_frame(cd, ref_h, ref_w)
% getframe can return off-by-one geometry between calls on some renderers,
% which makes imwrite('append') fail mid-GIF. Crop or edge-pad instead.
    [h,w,~] = size(cd);
    if h==ref_h && w==ref_w, out = cd; return; end
    out = cd(1:min(h,ref_h), 1:min(w,ref_w), :);
    if size(out,1) < ref_h, out = cat(1, out, repmat(out(end,:,:), ref_h-size(out,1),1,1)); end
    if size(out,2) < ref_w, out = cat(2, out, repmat(out(:,end,:), 1, ref_w-size(out,2),1)); end
end

function out = shrink(cd, target_w)
    w = size(cd,2);
    if w <= target_w, out = cd; return; end
    if exist('imresize','file')==2 || exist('imresize','builtin')==5
        out = imresize(cd, target_w/w, 'bilinear');
    else
        k = max(1, round(w/target_w));
        out = cd(1:k:end, 1:k:end, :);
    end
end

function write_summary_json(P, cfg, R, assets, wind, step_size)
    S.run    = P.profile_id;
    S.schema = 'wind-prediction-suite/1.1';
    S.generated_utc = char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss''Z'''));
    S.settings = struct( ...
        'history_sec', cfg.history_sec, 'predict_sec', cfg.predict_sec, ...
        'stride_sec', cfg.stride_sec, 'step_size', step_size, ...
        'selection_criterion', cfg.selection_criterion, ...
        'num_segments', cfg.num_segments, 'signal_length', numel(wind), ...
        'duration_sec', (numel(wind)-1)*step_size, ...
        'wind_min', min(wind), 'wind_max', max(wind), ...
        'wind_mean', mean(wind), 'wind_std', std(wind), ...
        'models_to_run', {cfg.models_to_run});

    S.model_order = R.summary_table.Model(:)';

    mm = struct();
    for i = 1:height(R.summary_table)
        nm  = R.summary_table.Model{i};
        key = matlab.lang.makeValidName(lower(nm));
        e = struct('name', nm, 'rank', i, ...
            'rmse', R.summary_table.RMSE(i), 'nrmse', R.summary_table.nRMSE(i), ...
            'mae', R.summary_table.MAE(i), 'mape', R.summary_table.MAPE(i), ...
            'smape', R.summary_table.sMAPE_pct(i), 'r2', R.summary_table.R2(i), ...
            'dm_pvalue', R.summary_table.DM_pvalue(i), 'time_sec', R.summary_table.Time_sec(i));
        if isfield(assets, key), e.assets = assets.(key); end
        mm.(key) = e;
    end
    S.models = mm;

    sb = struct('segment',{},'model',{},'rmse',{});
    for s = 1:cfg.num_segments
        sb(s).segment = s;
        sb(s).model   = R.best_model_per_segment{s};
        sb(s).rmse    = R.best_score(s);
    end
    S.segment_best = sb;

    fid = fopen(P.summary_json, 'w');
    fwrite(fid, jsonencode(S, 'PrettyPrint', true), 'char');
    fclose(fid);
end

function write_tex(f, T)
    fid = fopen(f, 'w');
    if fid == -1, return; end
    fprintf(fid, '\\begin{table}[htbp]\n\\centering\n');
    fprintf(fid, '\\caption{Wind Prediction Suite Benchmarking Metrics}\n');
    fprintf(fid, '\\label{tab:wind_prediction_summary}\n');
    fprintf(fid, '\\begin{tabular}{lrrrrrrrr}\n\\toprule\n');
    fprintf(fid, 'Model & RMSE & nRMSE & MAE & MAPE (\\%%) & sMAPE (\\%%) & $R^2$ & DM $p$ & Time (s) \\\\\n');
    fprintf(fid, '\\midrule\n');
    for i = 1:height(T)
        fprintf(fid, '%s & %.4f & %.4f & %.4f & %.2f & %.2f & %.4f & %.4f & %.2f \\\\\n', ...
            upper(T.Model{i}), T.RMSE(i), T.nRMSE(i), T.MAE(i), T.MAPE(i), ...
            T.sMAPE_pct(i), T.R2(i), T.DM_pvalue(i), T.Time_sec(i));
    end
    fprintf(fid, '\\bottomrule\n\\end{tabular}\n\\end{table}\n');
    fclose(fid);
end

function write_markdown(P, cfg, R, web, wind)
    if web
        fig_dir = 'prediction/web/figures'; sfx = '_full.png';
    else
        fig_dir = 'prediction/figures';     sfx = '_full_timeline.png';
    end
    fid = fopen(P.readme, 'w');
    if fid == -1, return; end
    T = R.summary_table;

    fprintf(fid, '# %s\n\n', strrep(P.profile_id, '_', ' '));
    fprintf(fid, 'Wind prediction benchmark. Generated automatically -- do not edit by hand.\n\n');

    fprintf(fid, '## Run settings\n\n| Parameter | Value |\n|---|---|\n');
    fprintf(fid, '| History window | %.2f s |\n', cfg.history_sec);
    fprintf(fid, '| Forecast horizon | %.2f s |\n', cfg.predict_sec);
    fprintf(fid, '| Stride | %.2f s |\n', cfg.stride_sec);
    fprintf(fid, '| Segments | %d |\n', cfg.num_segments);
    fprintf(fid, '| Wind range | %.2f to %.2f m/s |\n', min(wind), max(wind));
    fprintf(fid, '| Selection criterion | %s |\n\n', cfg.selection_criterion);

    fprintf(fid, '## Model ranking\n\n');
    fprintf(fid, '| Rank | Model | RMSE | nRMSE | MAE | MAPE %% | sMAPE %% | R2 | DM p | Time (s) |\n');
    fprintf(fid, '|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|\n');
    for i = 1:height(T)
        fprintf(fid, '| %d | `%s` | %.4f | %.4f | %.4f | %.2f | %.2f | %.4f | %.4f | %.2f |\n', ...
            i, T.Model{i}, T.RMSE(i), T.nRMSE(i), T.MAE(i), T.MAPE(i), ...
            T.sMAPE_pct(i), T.R2(i), T.DM_pvalue(i), T.Time_sec(i));
    end

    fprintf(fid, '\n## Per-model forecasts\n\n');
    for i = 1:height(T)
        nm = T.Model{i};
        fprintf(fid, '### %s (rank %d, RMSE %.4f)\n\n', upper(nm), i, T.RMSE(i));
        fprintf(fid, '![%s](%s/%s%s)\n\n', nm, fig_dir, nm, sfx);
    end

    fprintf(fid, '## Artifacts\n\n');
    fprintf(fid, '- `report/prediction_suite_report.txt` -- full written report\n');
    fprintf(fid, '- `report/prediction_suite_summary.json` -- machine-readable summary\n');
    fprintf(fid, '- `report/*.csv` -- metrics as CSV (GitHub renders these as tables)\n');
    fprintf(fid, '- `summary_table.tex` -- LaTeX table for the thesis\n');
    fprintf(fid, '- `prediction/figures/`, `prediction/animation/` -- full-resolution assets\n');
    if web
        fprintf(fid, '- `prediction/web/` -- downscaled assets used by report.html\n');
    end
    fprintf(fid, '- `prediction_suite_results.mat` -- raw numbers behind all of the above\n');
    fclose(fid);
end
