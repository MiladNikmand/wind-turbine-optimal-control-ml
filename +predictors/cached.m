function src = cached(profile_root, wind, step_size, history_sec, predict_sec, stride_sec, opts)
%PREDICTORS.CACHED  Serve saved per-segment forecasts to the control loop.
%
%   src = predictors.cached(profile_root, wind, step_size, ...
%                           history_sec, predict_sec, stride_sec, opts)
%
% Reads <profile_root>/prediction_suite_results.mat, written by
% run_prediction, and hands the control loop one forecast column per
% segment. No predictor is re-run, so a full controller sweep across all ten
% models costs nothing beyond the original prediction pass.
%
% opts:
%   .mode        'bestof' (default) | 'model' | 'truth'
%                  bestof -> the per-segment winner from best_idx
%                  model  -> always opts.model
%                  truth  -> perfect foresight, the upper-bound baseline
%   .model       predictor name, required when mode == 'model'
%   .on_missing  'persistence' (default) | 'bestof' | 'error'
%   .verbose     logical, default true
%
% Returns:
%   src.get(seg)       -> [future_pred, info]
%   src.seg_idx(seg)   -> [i0 i1] wind indices the forecast covers
%   src.num_segments, .window_size, .segment_length, .segment_stride
%   src.t_offset       time of the first forecast sample (= history_sec)
%   src.models, .mode, .label, .describe()
%
% There is deliberately no forecast before t = history_sec: segment 1 needs
% a full history window before it can predict anything, so the controlled
% timeline starts at src.t_offset, not at zero.

if nargin < 7 || isempty(opts), opts = struct(); end
if ~isfield(opts,'mode'),       opts.mode = 'bestof';         end
if ~isfield(opts,'on_missing'), opts.on_missing = 'persistence'; end
if ~isfield(opts,'verbose'),    opts.verbose = true;          end
if ~isfield(opts,'model'),      opts.model = '';              end

wind = wind(:);

D.wind       = wind;
D.mode       = lower(opts.mode);
D.on_missing = lower(opts.on_missing);
D.verbose    = logical(opts.verbose);

% round(), matching every predict_with_* and predictors.num_segments
D.window_size    = round(history_sec / step_size);
D.segment_length = round(predict_sec / step_size);
D.segment_stride = round(stride_sec  / step_size);
D.step_size      = step_size;

if ~strcmp(D.mode, 'truth')
    mat_file = fullfile(profile_root, 'prediction_suite_results.mat');
    if exist(mat_file,'file') ~= 2
        error('predictors:cached:noResults', ...
            ['No prediction_suite_results.mat in:\n  %s\n' ...
             'Run run_prediction for this profile first, or use opts.mode = ''truth''.'], ...
            profile_root);
    end

    S = load(mat_file, 'predictions','models_to_run','best_idx','num_segments', ...
                       'history_sec','predict_sec','stride_sec','step_size');

    for r = {'predictions','models_to_run','best_idx','num_segments'}
        if ~isfield(S, r{1})
            error('predictors:cached:badResults', ...
                '%s is missing field "%s".', mat_file, r{1});
        end
    end

    % A silent timing mismatch is the worst failure here: every array would
    % be the right SIZE but every segment index would point at a different
    % slice of wind than the suite predicted.
    check_timing('history_sec', history_sec, S);
    check_timing('predict_sec', predict_sec, S);
    check_timing('stride_sec',  stride_sec,  S);
    check_timing('step_size',   step_size,   S);

    D.predictions  = S.predictions;
    D.models       = S.models_to_run(:)';
    D.best_idx     = S.best_idx(:)';
    D.num_segments = double(S.num_segments);

    [nm, ns] = size(D.predictions);
    if ns < D.num_segments
        warning('predictors:cached:segmentCount', ...
            'predictions has %d columns but num_segments = %d. Using %d.', ns, D.num_segments, ns);
        D.num_segments = ns;
    end

    switch D.mode
        case 'model'
            if isempty(opts.model)
                error('predictors:cached:noModel', ...
                    'mode ''model'' needs opts.model. Available: %s', strjoin(D.models, ', '));
            end
            hit = find(strcmpi(D.models, opts.model), 1);
            if isempty(hit)
                error('predictors:cached:unknownModel', ...
                    'Model "%s" is not in this results file. Available: %s', ...
                    opts.model, strjoin(D.models, ', '));
            end
            D.model_idx = hit;
            D.label = D.models{hit};
        case 'bestof'
            D.model_idx = [];
            D.label = 'best-of-segment';
        otherwise
            error('predictors:cached:badMode', ...
                'mode must be ''model'', ''bestof'' or ''truth'' (got "%s").', D.mode);
    end

    if numel(D.best_idx) < D.num_segments
        error('predictors:cached:shortBestIdx', ...
            'best_idx has %d entries but there are %d segments.', numel(D.best_idx), D.num_segments);
    end
    if any(D.best_idx < 1 | D.best_idx > nm)
        error('predictors:cached:badBestIdx', 'best_idx has entries outside 1..%d.', nm);
    end
else
    D.predictions  = {};
    D.models       = {'truth'};
    D.best_idx     = [];
    D.model_idx    = [];
    D.label        = 'truth (perfect foresight)';
    D.num_segments = floor((numel(wind) - D.window_size - D.segment_length) / D.segment_stride) + 1;
end

max_seg = floor((numel(wind) - D.window_size - D.segment_length) / D.segment_stride) + 1;
if max_seg < 1
    error('predictors:cached:tooShort', ...
        'Wind profile is too short: %d samples, need %d (history %d + horizon %d).', ...
        numel(wind), D.window_size + D.segment_length, D.window_size, D.segment_length);
end
if D.num_segments > max_seg
    if D.verbose
        fprintf('[forecast] trimming %d -> %d segments\n', D.num_segments, max_seg);
    end
    D.num_segments = max_seg;
end

src.get            = @(seg) local_forecast(seg, D);
src.seg_idx        = @(seg) local_seg_idx(seg, D);
src.describe       = @() local_describe(D);
src.num_segments   = D.num_segments;
src.window_size    = D.window_size;
src.segment_length = D.segment_length;
src.segment_stride = D.segment_stride;
src.step_size      = D.step_size;
src.t_offset       = D.window_size * D.step_size;
src.models         = D.models;
src.mode           = D.mode;
src.label          = D.label;

if D.verbose
    fprintf('[forecast] %s\n', local_describe(D));
    fprintf('[forecast] %d segments | history %d | horizon %d | stride %d samples\n', ...
        D.num_segments, D.window_size, D.segment_length, D.segment_stride);
    fprintf('[forecast] control timeline starts at t = %.3f s\n', src.t_offset);
end
end

% =========================================================================
function check_timing(name, controller_value, S)
    if ~isfield(S, name) || isempty(S.(name)), return; end
    stored = S.(name);
    if abs(stored - controller_value) > 1e-9 * max(1, abs(stored))
        error('predictors:cached:timingMismatch', ...
            ['%s mismatch: the control script uses %.6g but the cached\n' ...
             'predictions were generated with %.6g. Segment indices would not line up.\n' ...
             'Re-run run_prediction with matching settings, or change the controller to match.'], ...
            name, controller_value, stored);
    end
end

% =========================================================================
function idx = local_seg_idx(seg, D)
    i0 = (seg-1)*D.segment_stride + 1;
    i1 = i0 + D.window_size - 1;
    idx = [i1 + 1, i1 + D.segment_length];
end

% =========================================================================
function [future_pred, info] = local_forecast(seg, D)
    if seg < 1 || seg > D.num_segments
        error('predictors:cached:segOutOfRange', ...
            'Segment %d requested but only 1..%d exist.', seg, D.num_segments);
    end

    span = local_seg_idx(seg, D);
    info = struct('segment',seg, 'idx_start',span(1), 'idx_end',span(2), ...
                  'model','', 'fallback','');

    if strcmp(D.mode, 'truth')
        future_pred = D.wind(span(1):span(2));
        info.model = 'truth';
        return;
    end

    if isempty(D.model_idx)
        m = D.best_idx(seg);
    else
        m = D.model_idx;
    end
    info.model = D.models{m};

    future_pred = D.predictions{m, seg};
    future_pred = future_pred(:);

    bad = isempty(future_pred) || numel(future_pred) ~= D.segment_length || ...
          ~all(isfinite(future_pred));

    if bad
        switch D.on_missing
            case 'error'
                error('predictors:cached:missingSegment', ...
                    'Model "%s" produced no usable forecast for segment %d.', info.model, seg);
            case 'bestof'
                mb  = D.best_idx(seg);
                alt = D.predictions{mb, seg};
                alt = alt(:);
                if ~isempty(alt) && numel(alt) == D.segment_length && all(isfinite(alt))
                    future_pred = alt;
                    info.fallback = sprintf('bestof:%s', D.models{mb});
                else
                    future_pred = repmat(D.wind(span(1)-1), D.segment_length, 1);
                    info.fallback = 'persistence';
                end
            otherwise
                future_pred = repmat(D.wind(span(1)-1), D.segment_length, 1);
                info.fallback = 'persistence';
        end
        if D.verbose
            fprintf('[forecast] seg %d: "%s" unusable -> %s\n', seg, info.model, info.fallback);
        end
    end
end

% =========================================================================
function s = local_describe(D)
    switch D.mode
        case 'truth'
            s = 'forecast source: TRUTH (perfect foresight upper bound)';
        case 'bestof'
            s = sprintf('forecast source: cached best-of-segment across %d models', numel(D.models));
        otherwise
            s = sprintf('forecast source: cached model "%s"', D.label);
    end
end
