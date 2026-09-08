function P = paths(project_root, profile_id, opts)
%REPORT.PATHS  The output folder layout for one wind profile.
%
%   P = report.paths(project_root, profile_id)
%   P = report.paths(project_root, profile_id, struct('create', true))
%
% Layout:
%   wind_profiles/<profile_id>/
%     README.md
%     summary_table.tex
%     prediction_suite_results.mat
%     report/     prediction_suite_report.txt, *.json, *.csv
%     prediction/
%       figures/      full-resolution PNGs (and per-segment, if enabled)
%       statistics/   the ten cross-model comparison figures
%       animation/    full-resolution GIFs
%       web/          downscaled mirror consumed by report.html
%     control/<run>/  same two-tier split, per forecast mode

    if nargin < 3 || isempty(opts), opts = struct(); end
    if ~isfield(opts, 'create'), opts.create = false; end

    P.project_root = project_root;
    P.profile_id   = profile_id;

    P.root      = fullfile(project_root, 'wind_profiles', profile_id);
    P.report    = fullfile(P.root, 'report');
    P.figures   = fullfile(P.root, 'prediction', 'figures');
    P.stats     = fullfile(P.root, 'prediction', 'statistics');
    P.animation = fullfile(P.root, 'prediction', 'animation');
    P.web         = fullfile(P.root, 'prediction', 'web');
    P.web_figures = fullfile(P.web, 'figures');
    P.web_anim    = fullfile(P.web, 'animation');
    P.control     = fullfile(P.root, 'control');

    P.mat_file     = fullfile(P.root, 'prediction_suite_results.mat');
    P.tex_file     = fullfile(P.root, 'summary_table.tex');
    P.readme       = fullfile(P.root, 'README.md');
    P.report_txt   = fullfile(P.report, 'prediction_suite_report.txt');
    P.summary_json = fullfile(P.report, 'prediction_suite_summary.json');
    P.control_json = fullfile(P.report, 'control_summary.json');

    if opts.create
        d = {P.root, P.report, P.figures, P.stats, P.animation, ...
             P.web, P.web_figures, P.web_anim, P.control};
        for i = 1:numel(d)
            if ~exist(d{i}, 'dir'), mkdir(d{i}); end
        end
    end
end
