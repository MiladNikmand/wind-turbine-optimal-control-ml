function n_seg = num_segments(N, step_size, history_sec, predict_sec, stride_sec)
% =========================================================================
%PREDICTORS.NUM_SEGMENTS  (was segment_initiator)
% Computes how many sliding-window segments fit inside a wind series of
% length N, using the *same* indexing convention every predict_with_*_segment
% function already assumes:
%
%   window_size      = round(history_sec / dt)
%   segment_length   = round(predict_sec / dt)
%   segment_stride   = round(stride_sec  / dt)
%
%   idx_start   = (seg-1)*segment_stride + 1
%   idx_end     = idx_start + window_size - 1
%   future_end  = idx_end   + segment_length
%
% A segment is only valid if future_end <= N. This function is model-
% agnostic: it doesn't know or care which predictor will consume the
% segments, it only knows the windowing geometry.
% =========================================================================

    dt              = step_size;
    window_size     = round(history_sec / dt);
    segment_length  = round(predict_sec / dt);
    segment_stride  = round(stride_sec  / dt);

    if segment_stride <= 0
        error('predictors:num_segments:badStride', 'stride_sec/step_size must round to a positive integer.');
    end

    max_seg = floor((N - window_size - segment_length) / segment_stride) + 1;
    n_seg = max(max_seg, 0);

    if n_seg == 0
        warning('predictors:num_segments:noSegments', ...
            ['Wind series too short for the requested history_sec/predict_sec/stride_sec. ' ...
             'Need at least %d samples (window %d + horizon %d), have %d.'], ...
             window_size + segment_length, window_size, segment_length, N);
    end
end
