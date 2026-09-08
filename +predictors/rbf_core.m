function varargout = rbf_core(action, varargin)
%PREDICTORS.RBF_CORE  The one and only RBF implementation, with model caching.
%
%   model = predictors.rbf_core('train', wind_data, window_size, ...
%                               num_centers, sigma_factor, use_lsqminnorm)
%   yhat  = predictors.rbf_core('predict', past_segment, model, segment_length)
%   cfg   = predictors.rbf_core('defaults')
%
%   predictors.rbf_core('cache', 'off')     disable caching (exact old behaviour)
%   predictors.rbf_core('cache', 'on')      re-enable (this is the default)
%   predictors.rbf_core('cache', 'clear')   drop the cached model
%   predictors.rbf_core('cache', 'reset')   drop it and zero the counters
%   s = predictors.rbf_core('stats')        hits / misses / trains
%
% ------------------------------------------------------------------------
% WHY THIS FILE EXISTS
%   train_rbf_model_fast and predict_wind_segment_fast used to be local
%   functions duplicated byte-for-byte at the bottom of BOTH
%   predict_with_rbf_segment.m and predict_with_rbf_arima_segment.m --
%   47 identical lines in two places. MATLAB scopes local functions to
%   their file, so rbf_arima physically could not call rbf's copy. Change
%   num_centers in one and the hybrid would silently keep the old value,
%   so you would be benchmarking RBF against a *different* RBF hidden
%   inside RBF+ARIMA. One definition now, both callers share it.
%
% ------------------------------------------------------------------------
% WHY CACHING
%   Measured on wind03 at dt=0.01: 11.0 s per segment, of which kmeans is
%   98%. The design matrix is 4501 x 1500 and kmeans runs 150 clusters x 5
%   replicates over it.
%
%   But look at what 'train' receives: wind_data, window_size, num_centers,
%   sigma_factor, use_lsqminnorm. NONE of them depend on the segment index.
%   Only past_segment -- used in 'predict' -- changes. So the suite was
%   rebuilding an identical global model 30 times for wind03 and 190 times
%   for wind01. (It is segment-independent precisely because 'train' uses
%   the whole wind_data vector, the look-ahead behaviour preserved from the
%   thesis code.)
%
%   Caching keys on the training inputs, so a genuine change to any of them
%   still forces a retrain.
%
%       wind03  30 segs:  330 s -> 15 s   (23x)
%       wind07  56 segs:  617 s -> 18 s   (35x)
%       wind01 190 segs: 2092 s -> 34 s   (62x)
%
%   ONE BEHAVIOURAL NOTE. kmeans starts from a random initialisation, so
%   the old code drew slightly different centres on every retrain. With the
%   cache, one set of centres serves the whole run. Predictions shift a
%   little as a result -- regenerate results rather than mixing cached and
%   uncached numbers. Under the suite's parfor the results were already
%   non-deterministic (each worker trained independently), and the cache is
%   per-worker for exactly the same reason.
%
%   Use predictors.rbf_core('cache','off') to restore the old
%   retrain-every-segment behaviour exactly.
% ------------------------------------------------------------------------

    persistent CACHE_ON CACHE_KEY CACHE_MODEL N_HIT N_MISS N_TRAIN
    if isempty(CACHE_ON)
        CACHE_ON = true;
        N_HIT = 0; N_MISS = 0; N_TRAIN = 0;
    end

    switch lower(action)
        case 'defaults'
            varargout{1} = struct( ...
                'num_centers',    150, ...
                'sigma_factor',   0.7, ...
                'use_lsqminnorm', true);

        case 'cache'
            mode = lower(varargin{1});
            switch mode
                case 'on',    CACHE_ON = true;
                case 'off',   CACHE_ON = false; CACHE_KEY = []; CACHE_MODEL = [];
                case 'clear', CACHE_KEY = []; CACHE_MODEL = [];
                case 'reset', CACHE_KEY = []; CACHE_MODEL = []; N_HIT = 0; N_MISS = 0; N_TRAIN = 0;
                otherwise
                    error('predictors:rbf_core:badCacheMode', ...
                        'Use ''on'', ''off'', ''clear'' or ''reset''.');
            end
            if nargout > 0, varargout{1} = CACHE_ON; end

        case 'stats'
            varargout{1} = struct('enabled', CACHE_ON, 'hits', N_HIT, ...
                'misses', N_MISS, 'trains', N_TRAIN, 'cached', ~isempty(CACHE_MODEL));

        case 'train'
            wind_data      = varargin{1};
            window_size    = varargin{2};
            num_centers    = varargin{3};
            sigma_factor   = varargin{4};
            use_lsqminnorm = varargin{5};

            if ~CACHE_ON
                N_TRAIN = N_TRAIN + 1;
                varargout{1} = local_train(wind_data, window_size, num_centers, ...
                    sigma_factor, use_lsqminnorm);
                return;
            end

            key = local_key(wind_data, window_size, num_centers, sigma_factor, use_lsqminnorm);
            if ~isempty(CACHE_KEY) && isequal(key, CACHE_KEY)
                N_HIT = N_HIT + 1;
                varargout{1} = CACHE_MODEL;
                return;
            end

            N_MISS  = N_MISS + 1;
            N_TRAIN = N_TRAIN + 1;
            CACHE_MODEL = local_train(wind_data, window_size, num_centers, ...
                sigma_factor, use_lsqminnorm);
            CACHE_KEY = key;
            varargout{1} = CACHE_MODEL;

        case 'predict'
            varargout{1} = local_predict(varargin{:});

        otherwise
            error('predictors:rbf_core:badAction', ...
                ['Unknown action "%s". Use ''train'', ''predict'', ''defaults'', ' ...
                 '''cache'' or ''stats''.'], action);
    end
end

% =========================================================================
function k = local_key(wind_data, window_size, num_centers, sigma_factor, use_lsqminnorm)
% Cheap O(N) fingerprint of the training inputs. A full element-wise
% comparison would cost more than it saves on a 6001-sample vector queried
% once per segment; these moments separate any two wind profiles in this
% project comfortably.
    w = wind_data(:);
    k = [ numel(w), sum(w), sum(w.^2), w(1), w(end), min(w), max(w), ...
          window_size, num_centers, sigma_factor, double(use_lsqminnorm) ];
end

% =========================================================================
function rbf_model = local_train(wind_data, window_size, num_centers, sigma_factor, use_lsqminnorm)
% Verbatim port of train_rbf_model_fast.
    min_val   = min(wind_data);
    max_val   = max(wind_data);
    norm_data = (wind_data - min_val) / max(max_val - min_val, eps);

    n_samples = length(norm_data) - window_size;
    X = zeros(n_samples, window_size);
    y = zeros(n_samples, 1);

    for i = 1:n_samples
        X(i,:) = norm_data(i:i+window_size-1)';
        y(i)   = norm_data(i+window_size);
    end

    [~, centers] = kmeans(X, num_centers, 'MaxIter', 300, 'Replicates', 5);
    Dcc   = pdist2(centers, centers);
    medcc = median(Dcc(Dcc>0));
    sigma = medcc * sigma_factor;

    DX  = pdist2(X, centers);
    phi = exp(-(DX.^2) / (2 * sigma^2));

    if use_lsqminnorm
        weights = lsqminnorm(phi, y);
    else
        weights = phi \ y;
    end

    rbf_model.centers     = centers;
    rbf_model.weights     = weights;
    rbf_model.sigma       = sigma;
    rbf_model.min_val     = min_val;
    rbf_model.max_val     = max_val;
    rbf_model.window_size = window_size;
end

% =========================================================================
function predicted = local_predict(past_segment, rbf_model, segment_length)
% Verbatim port of predict_wind_segment_fast.
    current_input = (past_segment(:) - rbf_model.min_val) / max(rbf_model.max_val - rbf_model.min_val, eps);
    predicted = zeros(segment_length,1);

    for k = 1:segment_length
        DX  = sum((rbf_model.centers - current_input').^2, 2);
        phi = exp(-DX / (2 * rbf_model.sigma^2));
        norm_pred = phi' * rbf_model.weights;
        predicted(k) = norm_pred * (rbf_model.max_val - rbf_model.min_val) + rbf_model.min_val;
        current_input = [current_input(2:end); norm_pred];
    end
end