function out = feasibility(step_size, history_sec, predict_sec, stride_sec, profile_id)
%WINDLIB.FEASIBILITY  How many segments each profile supports at a given config.
%
%   windlib.feasibility(0.01, 15, 1.5, 1.5)              print a table
%   T = windlib.feasibility(0.01, 15, 1.5, 1.5)          return it
%   n = windlib.feasibility(0.01, 15, 1.5, 1.5, 'wind06')  one profile
%
% Exists because a too-short profile otherwise fails deep inside the suite.
% wind06 is only 10 s long, so at the default 15/1.5/1.5 config it yields
% zero segments and the run dies. This tells you up front, and reports the
% longest history window that would work.
%
% Uses exactly the geometry every predict_with_*_segment assumes:
%   window = round(history/dt), horizon = round(predict/dt),
%   stride = round(stride/dt), segments = floor((N-window-horizon)/stride)+1

    if nargin < 5, profile_id = []; end

    if ~isempty(profile_id)
        R = windlib.registry(profile_id);
    else
        R = windlib.registry();
    end

    dt     = step_size;
    window = round(history_sec / dt);
    horizon= round(predict_sec / dt);
    stride = round(stride_sec  / dt);

    if stride <= 0
        error('windlib:feasibility:badStride', 'stride_sec/step_size must round to a positive integer.');
    end

    n = numel(R);
    id = cell(n,1); T = zeros(n,1); N = zeros(n,1);
    segs = zeros(n,1); max_hist = zeros(n,1); ok = false(n,1);

    for i = 1:n
        id{i} = R(i).id;
        T(i)  = R(i).T;
        N(i)  = round(R(i).T / dt) + 1;
        segs(i) = max(0, floor((N(i) - window - horizon) / stride) + 1);
        ok(i)   = segs(i) > 0;
        % longest history_sec that still yields at least one segment
        max_hist(i) = max(0, (N(i) - horizon - 1)) * dt;
    end

    out = table(id, T, N, segs, max_hist, ok, ...
        'VariableNames', {'Profile','T_sec','N_samples','Segments','MaxHistory_sec','Feasible'});

    if nargout == 0
        fprintf('\nFeasibility at dt=%g, history=%g s, predict=%g s, stride=%g s\n', ...
            dt, history_sec, predict_sec, stride_sec);
        fprintf('(window=%d, horizon=%d, stride=%d samples)\n\n', window, horizon, stride);
        disp(out);
        bad = out.Profile(~out.Feasible);
        if ~isempty(bad)
            fprintf(2, 'Infeasible: %s\n', strjoin(bad', ', '));
            for i = find(~ok)'
                fprintf(2, '  %s is %g s long -- needs history_sec <= %.2f at this horizon.\n', ...
                    id{i}, T(i), max_hist(i));
            end
            fprintf('\n');
        end
        clear out;
    elseif ~isempty(profile_id)
        out = segs(1);
    end
end
