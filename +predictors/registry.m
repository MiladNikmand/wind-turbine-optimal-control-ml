function reg = registry(name)
%PREDICTORS.REGISTRY  Single source of truth: model name -> function handle.
%
%   reg = predictors.registry()        struct of all handles (as before)
%   e   = predictors.registry('rbf')   one entry, with metadata
%
% Every entry satisfies the shared contract:
%
%   [future_pred, metrics, segment_time_sec] = fn( ...
%       wind_data, step_size, seg, history_sec, predict_sec, stride_sec, ...
%       show_plots, show_text)
%
% To add a model: drop the .m file in +predictors and add one line to the
% table below. Nothing else in the suite needs to change.
%
% The 'requires' field lets predictors.available() tell you up front which
% models this MATLAB install can actually run, instead of discovering it
% halfway through a batch.

    T = { ...
    %   name          handle                        requires              note
        'rbf',        @predictors.rbf,              {'stats'},            'baseline; kmeans/pdist2 need Statistics'
        'rbf_arima',  @predictors.rbf_arima,        {'stats','econ'},     'RBF base + ARIMA on residuals'
        'arima',      @predictors.arima,            {'econ'},             'Box-Cox + (p,d,q) grid, AIC selection'
        'svr',        @predictors.svr,              {'stats'},            'fitrsvm, gaussian kernel'
        'gp',         @predictors.gp,               {'stats'},            'fitrgp, matern52'
        'esn',        @predictors.esn,              {},                   'echo state net, hand-rolled'
        'gru',        @predictors.gru,              {},                   'GRU reservoir + ridge readout'
        'lstm',       @predictors.lstm,             {},                   'fixed LSTM + multi-tap ridge readout'
        'bilstm',     @predictors.bilstm,           {'deep_legacy'},      'trainNetwork -- REMOVED in R2024a+'
        'tcn',        @predictors.tcn,              {},                   'causal dilated conv + ridge readout'
    };

    if nargin == 0 || isempty(name)
        % Back-compatible shape: a plain struct of handles, exactly like
        % the original predictor_registry().
        reg = struct();
        for i = 1:size(T,1)
            reg.(T{i,1}) = T{i,2};
        end
        return;
    end

    hit = find(strcmpi(T(:,1), name), 1);
    if isempty(hit)
        error('predictors:registry:unknownModel', ...
            'Unknown model "%s". Available: %s', name, strjoin(T(:,1)', ', '));
    end
    reg = struct('name', T{hit,1}, 'fn', T{hit,2}, ...
                 'requires', {T{hit,3}}, 'note', T{hit,4});
end
