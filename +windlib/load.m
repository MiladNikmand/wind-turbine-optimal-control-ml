function [w, T, dt, info] = load(profile_id, step_size)
%WINDLIB.LOAD  Load a wind profile, resampled onto a uniform step_size grid.
%
%   [w, T, dt] = windlib.load('wind03', 0.01)
%   [w, T, dt] = windlib.load(3, 0.01)              % by menu number
%   [w, T, dt, info] = windlib.load(...)            % also return metadata
%
% Replaces wind_profile_loader.m. Same numbers out, but:
%   * no input() prompts        (was: step_size prompt + 7-way menu)
%   * no `close all; clear; clc` -- the old script wiped the caller's
%     workspace, which made it unusable from any driver
%   * no plotting                (use windlib.preview if you want the figure)
%   * data lives in +windlib/data/*.csv, not 950 lines of pasted coordinates
%
% Returns:
%   w     1xN row vector of wind speed (m/s)
%   T     duration in seconds
%   dt    = step_size, echoed for convenience
%   info  the windlib.registry entry, plus N and the realised grid spacing

    if nargin < 2 || isempty(step_size)
        error('windlib:load:noStepSize', ...
            'step_size is required, e.g. windlib.load(''wind03'', 0.01).');
    end
    if ~isscalar(step_size) || ~isfinite(step_size) || step_size <= 0
        error('windlib:load:badStepSize', 'step_size must be a positive scalar.');
    end

    R = windlib.registry(profile_id);
    T = R.T;

    switch R.source
        case 'csv'
            raw = local_read_csv(R.csv);
            % Every profile takes column 2. In the original loader, case 1/4
            % did W_data(:,2) at the call site while data_bringer02..05 did
            % it internally -- same result either way.
            wind_speed = raw(:, 2);
            [~, w] = windlib.resample(wind_speed, step_size, T);

        case 'synthetic'
            % ---- verbatim port of wind_profile_loader case 5 ----
            % Note this uses floor(T/step_size)+1, NOT the T/step_size+1 rule
            % that Intp_data uses. Both give the same count for integer
            % T/step_size, so the two paths agree in practice.
            N = floor(T / step_size) + 1;
            t = linspace(0, T, N);
            w = 12*ones(1,N) + sin(0.7*t);
            % 7-second moving-mean window. On a 100 s signal this flattens
            % most of the sin(0.7t) ripple -- intentional, kept as-is.
            w = smoothdata(w, 'movmean', 7/step_size);

        otherwise
            error('windlib:load:badSource', 'Unknown source "%s".', R.source);
    end

    % ---- optional min-max rescale (profile 4) ----
    if ~isempty(R.rescale)
        L = R.rescale(1);
        U = R.rescale(2);
        % PORTED VERBATIM, including the original's quirk: the loader wrote
        %     x_min = min(wind01);  x_max = max(wind);
        % i.e. min of one variable and max of another. At that point in the
        % script they held the same data, so the result is correct -- but the
        % asymmetry is preserved here deliberately rather than "fixed",
        % because project 1 is a faithful refactor.
        x_min = min(w);
        x_max = max(w);
        if x_max == x_min
            w = ones(size(w)) * (L + U) / 2;
        else
            w = L + (w - x_min) * (U - L) / (x_max - x_min);
        end
    end

    w  = reshape(w, 1, []);
    dt = step_size;

    if nargout > 3
        info = R;
        info.N = numel(w);
        info.dt_realised = T / (numel(w) - 1);
        info.step_size = step_size;
    end
end

% =========================================================================
function M = local_read_csv(name)
    here = fileparts(mfilename('fullpath'));
    f = fullfile(here, 'data', [name '.csv']);
    if exist(f, 'file') ~= 2
        error('windlib:load:missingData', ...
            'Wind data file not found:\n  %s\nIs +windlib/data/ intact?', f);
    end
    if exist('readmatrix', 'file')
        M = readmatrix(f);
    else
        M = csvread(f, 1, 0); %#ok<CSVRD>  % pre-R2019a fallback
    end
    M = M(all(isfinite(M), 2), :);
    if size(M, 2) < 2
        error('windlib:load:badData', '%s should have 2 columns, found %d.', f, size(M,2));
    end
end