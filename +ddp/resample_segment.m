function y = resample_segment(v, original_dt, new_dt)
%DDP.RESAMPLE_SEGMENT  Re-grid one predicted segment onto the controller step.
%
%   y = ddp.resample_segment(v, original_dt, new_dt)
%
% THE ONE NUMERICAL FIX CARRIED INTO PROJECT 1.
%
% The original resample_vector_to_dt treated a segment as spanning
% (N-1)*original_dt. But when segments are laid end to end, each one
% occupies N*original_dt of wall time, not (N-1). Resampling onto the
% shorter span dropped samples and the control timeline drifted against the
% wind timeline whenever dt_controller ~= step_size:
%
%   dt_controller   needed   old      lost/seg   drift over 30 segments
%   0.01            150      150      0          0.000 s
%   0.005           300      299      1          0.150 s
%   0.002           750      746      4          0.240 s
%   0.001           1500     1491     9          0.270 s
%
% When dt_controller == step_size the two agree exactly, so runs at the
% default settings are unaffected.

    v = v(:);
    N_orig = numel(v);

    if N_orig < 2
        y = repmat(v, max(1, round(N_orig * original_dt / new_dt)), 1);
        return;
    end

    t_old = (0:N_orig-1).' * original_dt;
    n_new = max(1, round(N_orig * original_dt / new_dt));
    t_new = (0:n_new-1).' * new_dt;

    % The final target sample can land up to one original interval past the
    % last source sample; linear extrapolation over that span is well posed.
    y = interp1(t_old, v, t_new, 'linear', 'extrap');

    if any(~isfinite(y))
        y = fillmissing(y, 'nearest');
    end
    y = y(:);
end
