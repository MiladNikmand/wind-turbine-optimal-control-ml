function export_control(profile_root, profile_id, run_label, C, opts)
%REPORT.EXPORT_CONTROL  Write DDP-HJB control results for report.html.
%   (was export_control_results)
%
%   Mirrors the asset layout of run_prediction_suite: a full-resolution tier
%   for the thesis and a downscaled tier for the web report, plus a JSON
%   summary the HTML reads.
%
%   report.export_control(profile_root, profile_id, run_label, C, opts)
%
%   run_label  short id for this run, e.g. 'bestof', 'rbf', 'truth'.
%              Multiple runs coexist; each becomes a selectable chip in the
%              report so forecast sources can be compared side by side.
%
%   C  struct of control results (all vectors same length unless noted):
%     .time           1xM  control timeline (s)
%     .Wr             1xM  achieved rotor speed
%     .Wr_ref         1xM  reference rotor speed, time-aligned
%     .Tg             1xM  generator torque
%     .Tg_ref         1xM  reference torque, time-aligned
%     .U              1xM  control input
%     .lambda         1xM  tip-speed ratio      (optional)
%     .Cp             1xM  power coefficient    (optional)
%     .wind_time      1xP  wind timeline (s)
%     .wind_true      1xP  true wind
%     .wind_pred      1xM  forecast actually fed to the controller (optional)
%     .segment_times  1xS  start time of each segment (for seam markers)
%     .RMSE_all       1xS  per-segment RMSE     (optional)
%     .cost_all       1xS  per-segment cost     (optional)
%     .forecast_mode  char 'bestof' | 'model' | 'truth'
%     .forecast_model char predictor name when mode == 'model'
%     .dt_controller  scalar
%     .settings       struct of anything else worth recording (optional)
%
%   opts (optional):
%     .export_web    logical, build the downscaled tier   (default true)
%     .make_gifs     logical, render animations           (default true)
%     .gif_fps       frames per second                    (default 4)
%     .gif_frames    max frames per animation             (default 40)
%     .dpi_hires     static PNG resolution                (default 150)
%     .dpi_web       web PNG resolution                   (default 96)
%     .web_width     web GIF width in px                  (default 640)
%     .web_colors    web GIF palette size                 (default 128)

if nargin < 5 || isempty(opts), opts = struct(); end
opts = set_default(opts, 'export_web', true);
opts = set_default(opts, 'make_gifs',  true);
opts = set_default(opts, 'gif_fps',    4);
opts = set_default(opts, 'gif_frames', 40);
opts = set_default(opts, 'dpi_hires',  150);
opts = set_default(opts, 'dpi_web',    96);
opts = set_default(opts, 'web_width',  640);
opts = set_default(opts, 'web_colors', 128);

run_key = matlab.lang.makeValidName(run_label);

%% ---------------------- folders ----------------------
P.root      = profile_root;
P.control   = fullfile(profile_root, 'control', run_label);
P.figures   = fullfile(P.control, 'figures');
P.animation = fullfile(P.control, 'animation');
P.web       = fullfile(P.control, 'web');
P.web_fig   = fullfile(P.web, 'figures');
P.web_anim  = fullfile(P.web, 'animation');
P.report    = fullfile(profile_root, 'report');

dirs = {P.control, P.figures, P.animation, P.report};
if opts.export_web, dirs = [dirs, {P.web, P.web_fig, P.web_anim}]; end
for d = dirs
    if ~exist(d{1}, 'dir'), mkdir(d{1}); end
end

fprintf('[control-export] %s / %s\n', profile_id, run_label);

%% ---------------------- normalise inputs ----------------------
C = ensure_row(C, {'time','Wr','Wr_ref','Tg','Tg_ref','U','lambda','Cp', ...
                   'wind_time','wind_true','wind_pred','segment_times', ...
                   'RMSE_all','cost_all'});

M = numel(C.time);
if M < 2
    error('report:export_control:emptyRun', 'C.time has %d samples -- nothing to plot.', M);
end

% Trim every co-indexed signal to the shortest common length. Segment
% concatenation can leave these off by a sample or two and an uncaught
% mismatch would kill the plot loop halfway through.
co = {'time','Wr','Wr_ref','Tg','Tg_ref','U','lambda','Cp','wind_pred'};
lens = [];
for f = co
    if isfield(C, f{1}) && ~isempty(C.(f{1})), lens(end+1) = numel(C.(f{1})); end %#ok<AGROW>
end
Mcommon = min(lens);
if Mcommon < M
    fprintf('[control-export] trimming signals to common length %d (was up to %d)\n', Mcommon, max(lens));
end
for f = co
    if isfield(C, f{1}) && ~isempty(C.(f{1})), C.(f{1}) = C.(f{1})(1:Mcommon); end
end
M = Mcommon;

%% ---------------------- panel definitions ----------------------
% Each panel becomes one row in the report: static | animated | zoomed.
panels = {};
panels{end+1} = mk('wr',     'Rotor speed tracking', '\omega_r (rad/s)', ...
    {C.Wr_ref, C.Wr}, {'Reference','Achieved'}, {'k--','b-'});
panels{end+1} = mk('tg',     'Generator torque',     'T_g (N\cdotm)', ...
    {C.Tg_ref, C.Tg}, {'Reference','Achieved'}, {'k--','r-'});
panels{end+1} = mk('u',      'Control input',        'u', ...
    {C.U}, {'u(t)'}, {'m-'});

if isfield(C,'wind_pred') && ~isempty(C.wind_pred)
    wt = interp1(C.wind_time, C.wind_true, C.time, 'linear', 'extrap');
    panels{end+1} = mk('wind', 'Wind: forecast vs truth', 'v (m/s)', ...
        {wt, C.wind_pred}, {'True wind','Forecast fed to controller'}, {'k-','c-'});
end
if isfield(C,'Cp') && ~isempty(C.Cp)
    panels{end+1} = mk('cp',   'Power coefficient',    'C_p', ...
        {C.Cp}, {'C_p'}, {'g-'});
end
if isfield(C,'lambda') && ~isempty(C.lambda)
    panels{end+1} = mk('lambda','Tip-speed ratio',     '\lambda', ...
        {C.lambda}, {'\lambda'}, {'b-'});
end

%% ---------------------- tracking error metrics ----------------------
err = C.Wr - C.Wr_ref;
Mt.rmse_wr  = sqrt(mean(err.^2, 'omitnan'));
Mt.mae_wr   = mean(abs(err), 'omitnan');
Mt.max_wr   = max(abs(err));
denom       = max(abs(C.Wr_ref), 1e-9);
Mt.mape_wr  = mean(abs(err) ./ denom, 'omitnan') * 100;
rng_ref     = max(C.Wr_ref) - min(C.Wr_ref);
Mt.nrmse_wr = Mt.rmse_wr / max(rng_ref, eps);

terr = C.Tg - C.Tg_ref;
Mt.rmse_tg  = sqrt(mean(terr.^2, 'omitnan'));

dtc = C.dt_controller;
Mt.energy_kj    = sum(abs(C.Tg .* C.Wr)) * dtc * 1e-3;
Mt.control_effort = sum(C.U.^2) * dtc;
Mt.du_max       = max(abs(diff(C.U))) / dtc;
if isfield(C,'Cp') && ~isempty(C.Cp), Mt.avg_cp = mean(C.Cp, 'omitnan'); else, Mt.avg_cp = NaN; end

% Seam discontinuity: the jump in Wr at each segment boundary. This is the
% number that quantifies the spikiness, so it belongs in the report.
seam_jump = [];
if isfield(C,'segment_times') && numel(C.segment_times) > 1
    for i = 2:numel(C.segment_times)
        k = find(C.time >= C.segment_times(i), 1, 'first');
        if ~isempty(k) && k > 1 && k <= M
            seam_jump(end+1) = abs(C.Wr(k) - C.Wr(k-1)); %#ok<AGROW>
        end
    end
end
if isempty(seam_jump)
    Mt.seam_jump_mean = NaN; Mt.seam_jump_max = NaN;
else
    Mt.seam_jump_mean = mean(seam_jump);
    Mt.seam_jump_max  = max(seam_jump);
end

fprintf('[control-export] RMSE(Wr) = %.4f | energy = %.2f kJ | seam jump max = %.4f\n', ...
    Mt.rmse_wr, Mt.energy_kj, Mt.seam_jump_max);

%% ---------------------- static figures ----------------------
assets = struct();
t0 = C.time(1); t1 = C.time(end);

% zoom window: centred on the worst tracking error, three segments wide
[~, kworst] = max(abs(err));
if isfield(C,'segment_times') && numel(C.segment_times) > 1
    seg_dur = median(diff(C.segment_times));
else
    seg_dur = (t1 - t0) / 10;
end
zc = C.time(min(max(kworst,1), M));
zoom_lo = max(t0, zc - 1.5*seg_dur);
zoom_hi = min(t1, zc + 1.5*seg_dur);
if zoom_hi - zoom_lo < seg_dur, zoom_hi = min(t1, zoom_lo + 2*seg_dur); end

fig = figure('Visible','off','Color','w','Position',[100 100 1200 500]);
for i = 1:numel(panels)
    p = panels{i};
    clf(fig);
    ax = axes('Parent', fig); hold(ax,'on');
    for j = 1:numel(p.series)
        plot(ax, C.time, p.series{j}, p.styles{j}, 'LineWidth', 1.3, 'DisplayName', p.names{j});
    end
    draw_seams(ax, C);
    grid(ax,'on'); xlabel(ax,'Time (s)'); ylabel(ax, p.ylabel);
    legend(ax, 'Location','best', 'Interpreter','tex');

    ylo = min(cellfun(@(v) min(v(isfinite(v))), p.series));
    yhi = max(cellfun(@(v) max(v(isfinite(v))), p.series));
    pad = 0.08 * max(yhi - ylo, eps);
    ylim(ax, [ylo-pad, yhi+pad]);

    % full view
    xlim(ax, [t0 t1]);
    title(ax, sprintf('%s — %s (full)', p.title, run_label), 'Interpreter','tex');
    exportgraphics(fig, fullfile(P.figures, [p.key '_full.png']), 'Resolution', opts.dpi_hires);
    if opts.export_web
        exportgraphics(fig, fullfile(P.web_fig, [p.key '_full.png']), 'Resolution', opts.dpi_web);
    end

    % zoomed view -> doubles as the poster frame for the zoomed GIF
    xlim(ax, [zoom_lo zoom_hi]);
    title(ax, sprintf('%s — %s (zoom)', p.title, run_label), 'Interpreter','tex');
    exportgraphics(fig, fullfile(P.figures, [p.key '_zoomed.png']), 'Resolution', opts.dpi_hires);
    if opts.export_web
        exportgraphics(fig, fullfile(P.web_fig, [p.key '_zoomed.png']), 'Resolution', opts.dpi_web);
    end

    assets.(p.key) = build_asset_paths(profile_id, run_label, p.key, opts.export_web, opts.make_gifs);
    fprintf('  - %s static\n', p.key);
end
close(fig);

%% ---------------------- animations ----------------------
if opts.make_gifs
    nframes = min(opts.gif_frames, max(2, M-1));
    kf = unique(round(linspace(2, M, nframes)));

    for i = 1:numel(panels)
        p = panels{i};
        gif_full   = fullfile(P.animation, [p.key '_full.gif']);
        gif_zoom   = fullfile(P.animation, [p.key '_zoomed.gif']);
        if opts.export_web
            gif_wfull = fullfile(P.web_anim, [p.key '_full.gif']);
            gif_wzoom = fullfile(P.web_anim, [p.key '_zoomed.gif']);
        end

        gfig = figure('Visible','off','Color','w','Position',[100 100 1280 640]);
        gax  = axes('Parent', gfig); hold(gax,'on');
        h = gobjects(1, numel(p.series));
        for j = 1:numel(p.series)
            h(j) = plot(gax, nan, nan, p.styles{j}, 'LineWidth', 1.6, 'DisplayName', p.names{j});
        end
        draw_seams(gax, C);
        grid(gax,'on'); xlabel(gax,'Time (s)'); ylabel(gax, p.ylabel);
        legend(gax, 'Location','northeast', 'Interpreter','tex');
        ylo = min(cellfun(@(v) min(v(isfinite(v))), p.series));
        yhi = max(cellfun(@(v) max(v(isfinite(v))), p.series));
        pad = 0.08 * max(yhi - ylo, eps);
        ylim(gax, [ylo-pad, yhi+pad]);

        ref_h = []; ref_w = [];
        for fi = 1:numel(kf)
            k = kf(fi);
            for j = 1:numel(p.series)
                set(h(j), 'XData', C.time(1:k), 'YData', p.series{j}(1:k));
            end

            % full
            xlim(gax, [t0 t1]);
            title(gax, sprintf('%s — %s  (t = %.2f s)', p.title, run_label, C.time(k)), 'Interpreter','tex');
            cd1 = getframe(gfig); cd1 = cd1.cdata;
            if isempty(ref_h), ref_h = size(cd1,1); ref_w = size(cd1,2); end
            cd1 = lock_frame(cd1, ref_h, ref_w);
            [im,cm] = rgb2ind(cd1, 256, 'nodither');
            write_gif(gif_full, im, cm, fi==1, opts.gif_fps);
            if opts.export_web
                sm = shrink(cd1, opts.web_width);
                [iw,cw] = rgb2ind(sm, opts.web_colors, 'nodither');
                write_gif(gif_wfull, iw, cw, fi==1, opts.gif_fps);
            end

            % zoomed: window slides to follow the leading edge
            half = 1.5*seg_dur;
            lo = max(t0, C.time(k) - 2*half);
            hi = min(t1, max(lo + 2*half, C.time(k) + 0.2*half));
            if hi <= lo, hi = lo + max(seg_dur, eps); end
            xlim(gax, [lo hi]);
            title(gax, sprintf('%s — %s  (zoom, t = %.2f s)', p.title, run_label, C.time(k)), 'Interpreter','tex');
            cd2 = lock_frame(getframe(gfig).cdata, ref_h, ref_w);
            [im,cm] = rgb2ind(cd2, 256, 'nodither');
            write_gif(gif_zoom, im, cm, fi==1, opts.gif_fps);
            if opts.export_web
                sm = shrink(cd2, opts.web_width);
                [iw,cw] = rgb2ind(sm, opts.web_colors, 'nodither');
                write_gif(gif_wzoom, iw, cw, fi==1, opts.gif_fps);
            end
        end
        close(gfig);
        fprintf('  - %s gifs (%d frames)\n', p.key, numel(kf));
    end
end

%% ---------------------- JSON summary ----------------------
json_path = fullfile(P.report, 'control_summary.json');

S = struct();
if exist(json_path, 'file')
    try
        fid = fopen(json_path,'r'); raw = fread(fid,'*char')'; fclose(fid);
        S = jsondecode(raw);
    catch ME
        warning('control_summary.json unreadable (%s) -- rebuilding.', ME.message);
        S = struct();
    end
end
S.schema = 'wind-control-suite/1.0';
S.profile = profile_id;
S.generated_utc = char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss''Z'''));

% existing runs -> cell array so field sets can differ without erroring
runs = {};
if isfield(S,'runs') && ~isempty(S.runs)
    if iscell(S.runs), runs = S.runs(:)'; else, runs = num2cell(S.runs(:)'); end
end

entry = struct();
entry.id    = run_label;
entry.key   = run_key;
entry.label = pretty_label(C);
entry.forecast_mode  = getdef(C,'forecast_mode','');
entry.forecast_model = getdef(C,'forecast_model','');
entry.dt_controller  = dtc;
entry.num_segments   = numel(getdef(C,'segment_times',[]));
entry.duration_sec   = t1 - t0;
entry.t_start        = t0;
entry.metrics        = Mt;
entry.panels         = panel_meta(panels);
entry.assets         = assets;
if isfield(C,'RMSE_all') && ~isempty(C.RMSE_all), entry.rmse_per_segment = C.RMSE_all; end
if isfield(C,'cost_all') && ~isempty(C.cost_all), entry.cost_per_segment = C.cost_all; end
if isfield(C,'settings'), entry.settings = C.settings; end

hit = [];
for k = 1:numel(runs)
    if isfield(runs{k},'id') && strcmp(runs{k}.id, run_label), hit = k; break; end
end
if isempty(hit), runs{end+1} = entry; else, runs{hit} = entry; end
S.runs = runs;

fid = fopen(json_path, 'w');
fwrite(fid, jsonencode(S, 'PrettyPrint', true), 'char');
fclose(fid);
fprintf('[control-export] summary -> %s (%d run(s))\n', json_path, numel(runs));

end

%% ======================= local helpers =======================
function o = set_default(o, f, v)
    if ~isfield(o, f) || isempty(o.(f)), o.(f) = v; end
end

function v = getdef(S, f, d)
    if isfield(S, f) && ~isempty(S.(f)), v = S.(f); else, v = d; end
end

function C = ensure_row(C, fields)
    for f = fields
        if isfield(C, f{1}) && ~isempty(C.(f{1}))
            C.(f{1}) = reshape(C.(f{1}), 1, []);
        end
    end
end

function p = mk(key, ttl, ylab, series, names, styles)
    p.key = key; p.title = ttl; p.ylabel = ylab;
    p.series = series; p.names = names; p.styles = styles;
end

function m = panel_meta(panels)
    m = struct('key',{},'title',{},'ylabel',{});
    for i = 1:numel(panels)
        m(i).key = panels{i}.key;
        m(i).title = panels{i}.title;
        m(i).ylabel = strrep(strrep(panels{i}.ylabel,'\',''),'cdot','.');
    end
end

function draw_seams(ax, C)
% Segment boundaries are where the stitching artefacts live, so mark them.
    if ~isfield(C,'segment_times') || numel(C.segment_times) < 2, return; end
    st = C.segment_times;
    if numel(st) > 60, st = st(round(linspace(1, numel(st), 60))); end
    for i = 2:numel(st)
        xline(ax, st(i), 'Color', [0.85 0.85 0.85], 'HandleVisibility','off');
    end
end

function a = build_asset_paths(profile_id, run_label, key, web, gifs)
    base = sprintf('wind_profiles/%s/control/%s', profile_id, run_label);
    a = struct();
    a.static_full_hires = sprintf('%s/figures/%s_full.png',   base, key);
    a.static_zoom_hires = sprintf('%s/figures/%s_zoomed.png', base, key);
    if web
        a.static_full   = sprintf('%s/web/figures/%s_full.png',   base, key);
        a.static_zoomed = sprintf('%s/web/figures/%s_zoomed.png', base, key);
    else
        a.static_full   = a.static_full_hires;
        a.static_zoomed = a.static_zoom_hires;
    end
    if gifs
        a.gif_full_hires   = sprintf('%s/animation/%s_full.gif',   base, key);
        a.gif_zoomed_hires = sprintf('%s/animation/%s_zoomed.gif', base, key);
        if web
            a.gif_full   = sprintf('%s/web/animation/%s_full.gif',   base, key);
            a.gif_zoomed = sprintf('%s/web/animation/%s_zoomed.gif', base, key);
        else
            a.gif_full   = a.gif_full_hires;
            a.gif_zoomed = a.gif_zoomed_hires;
        end
    end
end

function s = pretty_label(C)
    mode = '';
    if isfield(C,'forecast_mode'), mode = C.forecast_mode; end
    switch lower(mode)
        case 'truth',  s = 'Perfect foresight';
        case 'bestof', s = 'Best-of-segment forecast';
        case 'model'
            if isfield(C,'forecast_model') && ~isempty(C.forecast_model)
                s = sprintf('Forecast: %s', upper(C.forecast_model));
            else
                s = 'Single-model forecast';
            end
        otherwise, s = 'Control run';
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
    [h,w,~] = size(cd);
    if h == ref_h && w == ref_w, out = cd; return; end
    out = cd(1:min(h,ref_h), 1:min(w,ref_w), :);
    if size(out,1) < ref_h, out = cat(1, out, repmat(out(end,:,:), ref_h-size(out,1), 1, 1)); end
    if size(out,2) < ref_w, out = cat(2, out, repmat(out(:,end,:), 1, ref_w-size(out,2), 1)); end
end

function out = shrink(cd, target_w)
    w = size(cd,2);
    if w <= target_w, out = cd; return; end
    if exist('imresize','file') == 2 || exist('imresize','builtin') == 5
        out = imresize(cd, target_w/w, 'bilinear');
    else
        k = max(1, round(w/target_w));
        out = cd(1:k:end, 1:k:end, :);
    end
end
