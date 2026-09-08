function predicted_wind = stitch(predictions, chosen_idx, wind_data, ...
    step_size, history_sec, predict_sec, stride_sec)
% =========================================================================
%PREDICTORS.STITCH  (was stitch_best_predictions)
% Rebuilds one continuous predicted-wind vector (same length as wind_data)
% out of the per-segment forecasts, using whichever model index is given
% in chosen_idx(seg) for each segment (this is normally the winner picked
% by RMSE/MAE, but you can pass a fixed model index for every segment too).
%
% If stride_sec < predict_sec, consecutive prediction windows overlap in
% time. Overlapping timestamps are simply averaged across the segments
% that cover them, which tends to smooth out single-segment glitches --
% consistent with the "moving window as memory" idea: newer windows see
% more recent data, but older windows still contribute where they overlap.
%
% Outputs:
%   predicted_wind  Nx1, NaN wherever no segment's forecast reaches.
% =========================================================================

    dt              = step_size;
    window_size     = round(history_sec / dt);
    segment_length  = round(predict_sec / dt);
    segment_stride  = round(stride_sec  / dt);

    N            = length(wind_data);
    num_segments = size(predictions, 2);

    accum  = zeros(N,1);
    counts = zeros(N,1);

    for seg = 1:num_segments
        m = chosen_idx(seg);
        if m < 1 || isnan(m)
            continue;
        end
        pred = predictions{m, seg};
        if isempty(pred) || any(~isfinite(pred))
            continue;
        end

        idx_start    = (seg-1)*segment_stride + 1;
        idx_end      = idx_start + window_size - 1;
        future_start = idx_end + 1;
        future_end   = idx_end + segment_length;

        if future_end > N
            continue;
        end

        accum(future_start:future_end)  = accum(future_start:future_end)  + pred(:);
        counts(future_start:future_end) = counts(future_start:future_end) + 1;
    end

    predicted_wind = nan(N,1);
    has_pred = counts > 0;
    predicted_wind(has_pred) = accum(has_pred) ./ counts(has_pred);
end
